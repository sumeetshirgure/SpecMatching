#!/usr/bin/env bash
# End-to-end cross-check of the h_construct RTL against the software ball-graph build.
#
#   accelerator/tests/run_crosscheck.sh
#
# Steps: build the vector generator, dump the tables and the reference H, verilate the wrapper
# against the unmodified RTL, simulate every shot, then diff the two H's.
#
# Environment overrides:
#   HC_BUILD_DIR   cmake / verilator build dir        (default ${TMPDIR:-/tmp}/hc_dump_build)
#   HC_VEC_DIR     vector / output dir                (default accelerator/tests/vectors)
#   PYTHON         interpreter for compare.py         (default python3; stdlib only)
#
# Exits non-zero if any shot's H differs.

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

for tool in cmake verilator "$python"; do
    command -v "$tool" >/dev/null 2>&1 || { echo "missing required tool: $tool" >&2; exit 2; }
done

echo "== 1/5 build the vector generator"
cmake -S "$root/accelerator/tools" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$build_dir" -j "$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" >/dev/null

echo "== 2/5 dump tables and reference H"
mkdir -p "$vec_dir"
rm -f "$vec_dir"/dut_*.txt
"$build_dir/hc_dump" "$vec_dir"

# HC_NUM_SHOTS is written by hc_dump, so the shot count always matches the vectors on disk.
num_shots="$(sed -n 's/^`define HC_NUM_SHOTS \([0-9]*\).*/\1/p' "$vec_dir/params.vh")"
if [ -z "$num_shots" ] || [ "$num_shots" -lt 1 ]; then
    echo "could not read HC_NUM_SHOTS from $vec_dir/params.vh" >&2
    exit 2
fi

echo "== 3/5 verilate the wrapper + RTL (-Wall)"
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
obj_dir="$build_dir/obj_dir"
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

echo "== 4/5 simulate $num_shots shots"
# The simulator runs from the vector dir: ball_sram's INIT_FILE reaches $readmemh as a bare file
# name and is read when the model is constructed, and the shot files and outputs live beside it.
(cd "$vec_dir" && "$obj_dir/hc_sim" . "$num_shots")

echo "== 5/5 compare"
"$python" "$here/compare.py" "$vec_dir" "$num_shots"
