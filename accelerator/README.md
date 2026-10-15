# `h_construct` accelerator

RTL front end that builds the ball graph `H` for one shot: given the shot's syndrome and the
per-detector ball and boundary tables, it emits `H`'s defect-defect edges and boundary edges.

## Modules

`accelerator/Verilog/`:

| file | purpose |
|---|---|
| `h_construct.v` | Top level: slices the boundary word once, instantiates the address generator, boundary check, FIFO, `NUM_LANES` search lanes and the collector, and carries the `{valid, u}` shadow pipeline that keeps each batch's sideband aligned with its lanes. |
| `hc_addr_gen.v` | Walks the syndrome: issues one boundary-table read per defect and one ball-table read per non-empty ball, then issues that ball's chunks one per unstalled cycle, with the lane bound `(size - 1) mod NUM_LANES` on the last chunk. |
| `hc_bnd_check.v` | One registered boundary decision per ball, straight off the held boundary word: emits `{u, bcost_w}` when the word's `hb` bit is set. |
| `hc_collect.v` | Serialises one batch's lane hits into the edge stream, one record per cycle, holding the batch and back-pressuring the lanes until its hits have drained. |
| `hc_fifo.v` | Control FSM plus the batch FIFO: cuts each chunk out of the held ball word, buffers batches with their sideband, generates `oGo`/`oStall`, and raises `oDone` once every stage is quiet. |
| `hc_lane.v` | One search lane: filters an entry on `e_tgt > iDet` and the lane bound, then looks the target up in the sorted syndrome with a pipelined binary search, one probe level per stage. |
| `ball_sram.v` | Single-port synchronous read SRAM model, `$readmemh`-initialised from `INIT_FILE`, holding its output while `iEn` is low. |
| `tb/hc_top.v` | Verilator top wrapper: `h_construct` plus the two table SRAMs, exposing only the build interface so the C++ harness never sees the SRAM read ports. |

## What the test compares against

The cross-check diffs the `H` the RTL emits — defect-defect edges and boundary edges — against the
software reference built by `build_ball_graph(..., BallGraphBuildMode::SCAN)` in
`src/specmatching/spec_matching/manifold/ball_graph.cc`, on the same ball table and the same
syndrome. Both sides are dumped as `E u v w` / `B u w` lines, and `accelerator/tests/compare.py`
canonicalises each into a sorted record list before diffing, so emission order does not matter but
multiplicity does. Current result: `CROSSCHECK PASS 20/20 shots`.

## How to run

```sh
accelerator/tests/run_crosscheck.sh
```

It builds the vector generator, dumps the tables and the reference `H`, verilates `tb/hc_top.v`
against the unmodified RTL, simulates every shot, and diffs. It exits non-zero if any shot's `H`
differs.

Requires `cmake`, `verilator` (5.052 here), a C++17 compiler and `python3` (`compare.py` is standard
library only). The runner clears a stale `VERILATOR_ROOT`, which would otherwise make verilator
refuse to start, and passes `-CFLAGS -std=gnu++17`, because this verilator install adds no `-std` of
its own and the host compiler defaults to C++98.

Environment overrides:

| variable | meaning | default |
|---|---|---|
| `HC_BUILD_DIR` | cmake / verilator build dir | `${TMPDIR:-/tmp}/hc_dump_build` |
| `HC_VEC_DIR` | vector and output dir | `accelerator/tests/vectors` |
| `PYTHON` | interpreter for `compare.py` | `python3` |
