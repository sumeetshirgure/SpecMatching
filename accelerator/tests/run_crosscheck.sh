#!/usr/bin/env bash
# End-to-end cross-check of the h_construct RTL against the software ball-graph build.
#
#   accelerator/tests/run_crosscheck.sh [hc_dump options...]
#
# Steps: build the table compiler, compile the tables and dump the vectors and the reference H,
# verilate the wrapper against the unmodified RTL, simulate every shot, then diff the two H's and
# the RTL's per-shot cycle counts against hc_dump's cycle model.
#
# Any arguments are passed to hc_dump (e.g. `--d 7 --rounds 21 --p 0.003 --shots 50`, or
# `--no-prune` / `--no-dedup`); the default is the checked-in vector set: the first 20 shots of a
# d = 5, 5-round, p = 5e-3, 60-shot corpus.
#
# Environment overrides:
#   HC_BUILD_DIR   cmake / verilator build dir        (default ${TMPDIR:-/tmp}/hc_dump_build)
#   HC_VEC_DIR     vector / output dir                (default accelerator/tests/vectors)
#   PYTHON         interpreter for compare.py         (default python3; stdlib only)
#
# Exits non-zero if any shot's H or cycle count differs.

set -euo pipefail

# A stale VERILATOR_ROOT in the environment makes verilator refuse to run at all
# ("VERILATOR_ROOT is set to inconsistent path"), so it is cleared for this script. The installed
# verilator finds its own root.
unset VERILATOR_ROOT

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"

build_dir="${HC_BUILD_DIR:-${TMPDIR:-/tmp}/hc_dump_build}"
vec_dir="${HC_VEC_DIR:-$root/accelerator/tests/vectors}"
python="${PYTHON:-python3}"

if [ "$#" -eq 0 ]; then
    set -- --d 5 --rounds 5 --p 0.005 --shots 60 --keep 20 --seed 7
fi

for tool in cmake verilator "$python"; do
    command -v "$tool" >/dev/null 2>&1 || { echo "missing required tool: $tool" >&2; exit 2; }
done

echo "== 1/6 build the table compiler"
# Reuse the main build's dependency checkouts when they exist, so configuring needs no network.
fetch_flags=()
for dep in stim googletest; do
    if [ -d "$root/obj/_deps/$dep-src" ]; then
        upper="$(echo "$dep" | tr '[:lower:]' '[:upper:]')"
        fetch_flags+=("-DFETCHCONTENT_SOURCE_DIR_$upper=$root/obj/_deps/$dep-src")
    fi
done
cmake -S "$root/accelerator/tablegen" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release -Wno-dev ${fetch_flags[@]+"${fetch_flags[@]}"} >/dev/null
cmake --build "$build_dir" --target hc_dump -j "$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" >/dev/null

echo "== 2/6 compile tables, dump vectors and reference H"
mkdir -p "$vec_dir"
rm -f "$vec_dir"/dut_*.txt "$vec_dir"/syn_*.txt "$vec_dir"/ref_*.txt
"$build_dir/hc_dump" "$vec_dir" "$@"

# HC_NUM_SHOTS is written by hc_dump, so the shot count always matches the vectors on disk.
num_shots="$(sed -n 's/^`define HC_NUM_SHOTS \([0-9]*\).*/\1/p' "$vec_dir/params.vh")"
if [ -z "$num_shots" ] || [ "$num_shots" -lt 1 ]; then
    echo "could not read HC_NUM_SHOTS from $vec_dir/params.vh" >&2
    exit 2
fi

echo "== 3/6 verilate the wrapper + RTL (-Wall)"
# Every source is named explicitly rather than globbed, so adding or removing a file under
# Verilog/tb never silently changes what is simulated.
#
# The two -Wno flags cover pre-existing warnings in the RTL, which no longer carries the inline
# `lint_off` waivers it once did: UNUSEDSIGNAL for the unread high bits of hc_addr_gen's `ch_full`
# and `size_m1`, and UNSIGNED for h_construct's `LANE_NO <= batch_hi` lane bound, whose comparison
# is constant at the two end lanes. Both are intended by the design. They are passed here instead
# of edited into the .v files; the wrapper itself is clean under plain -Wall.
#
# -CFLAGS -std=gnu++17 is not optional here: this verilator install reports an empty
# CFG_CXXFLAGS_STD, so it adds no -std of its own, and the host compiler defaults to C++98
# (__cplusplus 199711L), under which verilator's own verilatedos.h does not compile.
#
# The model is rebuilt per run because its port widths come from the generated params.vh.
obj_dir="$build_dir/obj_dir"
rm -rf "$obj_dir"
verilator --cc --exe --build -Wall -Wno-UNUSEDSIGNAL -Wno-UNSIGNED \
    -CFLAGS -std=gnu++17 \
    "+incdir+$vec_dir" \
    --top-module hc_top \
    --Mdir "$obj_dir" \
    -o hc_sim \
    "$root/accelerator/Verilog/tb/hc_top.v" \
    "$root/accelerator/Verilog/h_construct.v" \
    "$root/accelerator/Verilog/hc_addr_gen.v" \
    "$root/accelerator/Verilog/hc_bnd_check.v" \
    "$root/accelerator/Verilog/hc_collect.v" \
    "$root/accelerator/Verilog/hc_fifo.v" \
    "$root/accelerator/Verilog/hc_lane.v" \
    "$root/accelerator/Verilog/ball_sram.v" \
    "$here/sim_main.cc" >/dev/null

echo "== 4/6 simulate $num_shots shots"
# The simulator runs from the vector dir: ball_sram's INIT_FILE reaches $readmemh as a bare file
# name and is read when the model is constructed, and the shot files and outputs live beside it.
(cd "$vec_dir" && "$obj_dir/hc_sim" . "$num_shots") | tail -3

echo "== 5/6 compare H"
"$python" "$here/compare.py" "$vec_dir" "$num_shots"

echo "== 6/6 compare cycles against the model"
if diff "$vec_dir/model_cycles.txt" "$vec_dir/dut_cycles.txt" >/dev/null; then
    echo "CYCLES PASS $num_shots/$num_shots shots"
else
    echo "CYCLES FAIL: RTL cycle counts differ from hc_dump's model (model < > rtl):" >&2
    diff "$vec_dir/model_cycles.txt" "$vec_dir/dut_cycles.txt" | head -20 >&2
    exit 1
fi
