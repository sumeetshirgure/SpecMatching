# `h_construct` accelerator

RTL front end that builds the ball graph `H` for one shot: given the shot's syndrome and the
ball and boundary tables, it emits `H`'s defect-defect edges and boundary edges.

## Tables: time-shift encoding

The tables are compiled offline by `accelerator/tablegen/hc_dump.cc` from the software ball tables
(`compile_ball_tables`), and are plain SRAM contents, so the accelerator stays reprogrammable for
any DEM.

- **Ball table**: one fixed-width word per *unique row*. An entry is `{w, delta}` with
  `delta = target - source`, so a row does not name its source detector. When the DEM's weights do
  not change over time, the balls of one stabilizer in different bulk rounds are shifts of each other
  and hold the *same* row. Rows are shared only between detectors with the **same space coordinate**
  (every stim detector coordinate except the time one, `--time-coord`). Detectors at different space
  positions never share a row, even if their `{w, delta}` lists coincide. Lists are compared
  exactly, so nothing is assumed about the DEM: a DEM that is not time invariant just yields more
  rows, and the result stays correct.
- **Pruning**: entries with `target <= source` are dropped offline. `H` keeps each pair once, from
  its lower endpoint, and that test depends only on the row, so the lanes do not need it. Pruning
  halves the rows and the chunks streamed. The lanes still pass only `delta > 0` (a sign-bit test),
  so the same RTL runs unpruned tables (`--no-prune`). Padding entries are written as
  `delta = -1` so they never pass.
- **Boundary table**: one word per detector, `{row, hb, bw, size}`: the detector's ball-table row,
  its boundary edge (`hb` set when `bcost <= T`), and its row length. It is the one table indexed
  by detector id, so it is O(detectors), at about 40 bits each.

At d = 31, p = 1e-3, the 29,760 detectors share 7,888 rows (8.2 per space position), and the tables
shrink from 236.9 Mbit to 46.6 Mbit (5.1x). The effect grows with d, because the number of rows per
space position is set by the radius, not by the round count. See `sweep_tshift/table_sizes.md`.

## Modules

`accelerator/Verilog/`:

| file | purpose |
|---|---|
| `h_construct.v` | Top level: slices the boundary word once (`{row, hb, bw, size}`), instantiates the address generator, boundary check, FIFO, `NUM_LANES` search lanes and the collector, and carries the `{valid, u}` shadow pipeline that keeps each batch's sideband aligned with its lanes. |
| `hc_addr_gen.v` | Walks the syndrome: issues one boundary-table read per defect, latches the defect's row pointer, and issues one ball-table read (at that row) per non-empty ball. It then issues that ball's chunks one per unstalled cycle, with the lane bound `(size - 1) mod NUM_LANES` on the last chunk. |
| `hc_bnd_check.v` | One registered boundary decision per ball, straight off the held boundary word: emits `{u, bcost_w}` when the word's `hb` bit is set. |
| `hc_collect.v` | Serialises one batch's lane hits into the edge stream, one record per cycle, holding the batch and back-pressuring the lanes until its hits have drained. |
| `hc_fifo.v` | Control FSM plus the batch FIFO: cuts each chunk out of the held ball word, buffers batches with their sideband (including the source detector), generates `oGo`/`oStall`, and raises `oDone` once every stage is quiet. |
| `hc_lane.v` | One search lane: forms `target = iDet + sext(delta)`, filters on `delta > 0` and the lane bound, then looks the target up in the sorted syndrome with a pipelined binary search, one probe level per stage. |
| `ball_sram.v` | Single-port synchronous read SRAM model, `$readmemh`-initialised from `INIT_FILE`, holding its output while `iEn` is low. |
| `tb/hc_top.v` | Verilator top wrapper: `h_construct` plus the two table SRAMs, exposing only the build interface so the C++ harness never sees the SRAM read ports. |

New parameters: `LOG_NUM_ROWS` (ball-table address width) and `DELTA_W` (signed offset width,
`<= LOG_NUM_DET`). `LANE_W = WEIGHT_W + DELTA_W`. `WEIGHT_W` is sized by the compiler from 2T, with
16 bits as the floor (17 at p = 5e-4).

## What the test compares against

The cross-check diffs the `H` the RTL emits (defect-defect edges and boundary edges) against the
software reference built by `build_ball_graph(..., BallGraphBuildMode::SCAN)` in
`src/specmatching/spec_matching/manifold/ball_graph.cc`, on the **uncompressed** ball table and the
same syndrome. Both sides are dumped as `E u v w` / `B u w` lines, and `accelerator/tests/compare.py`
canonicalises each into a sorted record list before diffing, so emission order does not matter but
multiplicity does.

It also diffs the RTL's per-shot cycle count against `hc_dump`'s cycle model
(`model_cycles.txt` vs `dut_cycles.txt`). The model is a register-level replay of the control path,
and it is what `sweep_tshift/` reports at sizes where full RTL simulation is slow.

Configurations checked (H and cycles, all exact):

| corpus | tables | lanes |
|---|---|---|
| d = 5, 5 rounds, p = 5e-3, 20 shots (checked-in vectors) | default | 64 |
| d = 5, 15 rounds, p = 5e-3, 60 shots | default, `--no-prune`, `--no-dedup --no-prune` | 64 |
| d = 7, 21 rounds, p = 3e-3, 60 shots | default | 16, 128 |
| d = 9, 27 rounds, p = 2e-3 / 5e-4 (17-bit weights), 60 shots | default | 32, 64 |
| d = 17, 17 rounds, p = 1e-3 / 3e-3, 100 / 40 shots | default | 64 |

## How to run

```sh
accelerator/tests/run_crosscheck.sh                                   # checked-in vectors
accelerator/tests/run_crosscheck.sh --d 7 --rounds 21 --p 0.003 --shots 60 --lanes 16
python3 accelerator/tablegen/sweep.py <build-dir>/hc_dump accelerator/sweep_tshift
```

The runner builds `hc_dump` against the repository's own library (`accelerator/tablegen/CMakeLists.txt`
pulls in the root project `EXCLUDE_FROM_ALL`, so only `libspecmatching` and `libstim` are compiled,
using `obj/_deps` checkouts when present). It then compiles the tables, dumps the vectors and the
reference `H`, verilates `tb/hc_top.v` against the unmodified RTL, simulates every shot, and diffs
both `H` and cycles. Any arguments go to `hc_dump`. It exits non-zero if any shot differs.

The table compiler lives in `tablegen/` rather than `tools/` because the root `.gitignore` ignores
every `tools` path. That is why the previous generator was never committed.

Requires `cmake`, `verilator` (5.054 here), a C++17 compiler and `python3` (`compare.py` and
`sweep.py` are standard library only). The runner clears a stale `VERILATOR_ROOT`, which would
otherwise make verilator refuse to start, and passes `-CFLAGS -std=gnu++17`, because this verilator
install adds no `-std` of its own and the host compiler defaults to C++98.

Environment overrides:

| variable | meaning | default |
|---|---|---|
| `HC_BUILD_DIR` | cmake / verilator build dir | `${TMPDIR:-/tmp}/hc_dump_build` |
| `HC_VEC_DIR` | vector and output dir | `accelerator/tests/vectors` |
| `PYTHON` | interpreter for `compare.py` | `python3` |
