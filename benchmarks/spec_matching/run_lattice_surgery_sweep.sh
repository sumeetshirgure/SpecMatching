#!/bin/bash
# The lattice-surgery sweep of docs/lattice_surgery_profiler.md §2.3, as cut for the paper:
# p in {5e-4, 1e-3} (no 3e-3), d in {17, 25}, plus N = 1 as the q reference.
#
#   benchmarks/spec_matching/run_lattice_surgery_sweep.sh <results_dir> [shots] [jobs]
#
# One profiler process per (N, d, p, orientation) corpus, each running its k x W cells into
# <results_dir>/N{N}_d{d}_p{p}_{orientation}/. Up to `jobs` (default 5, the M5 Pro's performance-core
# count) run at once, largest first. Resumable: a corpus whose run.log ends in status=ok is skipped,
# so re-running the same command after an interruption only redoes what did not finish. When every
# corpus is done, the tables and figures go into <results_dir> itself.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
RESULTS="$1"
SHOTS="${2:-50000}"
JOBS="${3:-5}"
BIN="$ROOT/obj/lattice_surgery_profiler"
PY="$HOME/.virtualenvs/qq/bin/python"
mkdir -p "$RESULTS"
RESULTS="$(cd "$RESULTS" && pwd)"

(cd "$ROOT/obj" && cmake .. > /dev/null && make lattice_surgery_profiler -j8 | tail -1)

# Largest first, so the long corpora are not the last ones started.
CORPORA=()
for d in 25 17; do for N in 8 4 2 1; do for p in 0.001 0.0005; do for o in long short; do
    CORPORA+=("$N $d $p $o")
done; done; done; done

run_corpus() {
    read -r N d p o <<< "$1"
    dir="$RESULTS/N${N}_d${d}_p${p}_${o}"
    if [ -f "$dir/run.log" ] && grep -q "^status=ok" "$dir/run.log"; then
        echo "skip  $dir (done)"
        return 0
    fi
    rm -rf "$dir"
    echo "start $dir $(date +%T)"
    if "$BIN" --N "$N" --k 4,7 --d "$d" --p "$p" --T 2 --W 2,4 --orientation "$o" \
        --shots "$SHOTS" --warmup 1000 --verify --jobs-log sample --out "$dir" > "$dir.stdout" 2>&1; then
        echo "done  $dir $(date +%T)"
    else
        echo "FAIL  $dir $(date +%T) (see $dir.stdout and $dir/run.log)"
    fi
}
export -f run_corpus
export RESULTS SHOTS BIN

{
    echo "lattice_surgery sweep: shots=$SHOTS jobs=$JOBS started=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "git=$(cd "$ROOT" && git rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "cpu=$(sysctl -n machdep.cpu.brand_string) perf_cores=$(sysctl -n hw.perflevel0.physicalcpu)"
    echo "grid: N in {1,2,4,8}, k in {4,7}, d in {17,25}, p in {5e-4,1e-3}, T=2, W in {2,4}, orientation in {long,short}"
} >> "$RESULTS/sweep.log"

# caffeinate -i keeps the machine from idle-sleeping for as long as the sweep runs.
printf '%s\n' "${CORPORA[@]}" | caffeinate -i xargs -P "$JOBS" -I{} bash -c 'run_corpus "$@"' _ {} | tee -a "$RESULTS/sweep.log"

if grep -q "^FAIL" "$RESULTS/sweep.log"; then
    echo "some corpora failed; fix and re-run this command to resume" | tee -a "$RESULTS/sweep.log"
fi
"$PY" "$HERE/analyze_lattice_surgery.py" "$RESULTS" --out "$RESULTS"
"$PY" "$HERE/plot_lattice_surgery.py" "$RESULTS" --out "$RESULTS/figures" --formats png,pdf
echo "finished=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$RESULTS/sweep.log"
