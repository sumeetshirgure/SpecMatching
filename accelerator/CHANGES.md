# Changes from `main`: time-shift ball tables

Branch `hopt-tshift`. Every change is inside `accelerator/`.

## Summary

The ball table used to hold one full ball per detector, with absolute target ids. Its size grew as
O(d³) (rounds × stabilizers), and it was ~99.8% of the accelerator's area. The table now stores
target **offsets** (`delta = target − source`). Balls that are time shifts of one another become
identical and are stored once. The boundary table gains a per-detector row pointer.

Headline (`sweep_tshift/table_sizes.md`, 10,000 shots per point, L = 64):

| d | detectors | unique rows | tables, `main` → this branch | reduction | mean cycles at p = 1e-3 |
|--:|--:|--:|--:|--:|--:|
| 17 | 4,896 | 2,414 | 36.5 → 13.2 Mbit | 2.8× | 302 → 218 |
| 25 | 15,600 | 5,158 | 120.2 → 29.3 Mbit | 4.1× | 1,018 → 710 |
| 31 | 29,760 | 7,888 | 236.9 → 46.6 Mbit | 5.1× | 1,998 → 1,376 |

The number of rows per space position stays at ~8.2–8.9 at every d. It is set by the ball radius,
not the round count, so the tables now grow as O(T·d²).

## Table format

| | `main` | this branch |
|---|---|---|
| Ball-table word | one per detector, addressed by detector id | one per **unique row**, addressed by a row pointer |
| Ball entry | `{w, target}`, `WEIGHT_W + LOG_NUM_DET` bits | `{w, delta}`, `WEIGHT_W + DELTA_W` bits, `delta` two's complement |
| Entries kept | the whole ball (both directions of every pair) | only `target > source` (**pruned** offline); `--no-prune` keeps all |
| Padding entry | zero | `delta = −1`, which never passes the lane filter |
| Boundary word | `{hb, bw, size}` | `{row, hb, bw, size}`, with `row` in the top `LOG_NUM_ROWS` bits |
| `WEIGHT_W` | fixed at 16 | sized from 2T, minimum 16 (17 at p = 5e-4, where 16 bits overflow for both layouts) |

**Deduplication rule.** Two detectors share a row only if both of these hold:

1. They have the same **space coordinate**: every stim detector coordinate except the time one
   (`--time-coord`, default last). A detector with no coordinates never shares a row.
2. Their `{w, delta}` lists are **exactly equal**.

Detectors at different space positions never share a row, even when their lists coincide. This
keeps the hardware reprogrammable for DEMs whose space-like weights vary within a layer. Nothing
about the DEM is assumed: a DEM that is not time invariant just yields more rows, and the result is
still correct. Rows with no entries are not stored, since no ball read is issued for them.

## RTL (`Verilog/`)

| file | change |
|---|---|
| `h_construct.v` | New parameters `LOG_NUM_ROWS` and `DELTA_W` (`LANE_W = WEIGHT_W + DELTA_W`). `iBndTabRdata` is widened by `LOG_NUM_ROWS` and the `row` field is sliced out of it. `oBallAddr` is now `LOG_NUM_ROWS` bits wide. Both parameters are passed down to the address generator and lanes. |
| `hc_addr_gen.v` | New `iBndRow` input. On each `load` it latches the defect's row pointer into `row`, next to `ccnt`/`bhi`. `oBallAddr` is `row` instead of `det_i`. |
| `hc_lane.v` | An entry is `{w, delta}`. The lane computes `target = iDet + sext(delta)`, which replaces the `e_tgt > iDet` comparator with an adder of the same width. The filter becomes `!delta[MSB]` (delta is never 0). The binary search is unchanged. |
| `tb/hc_top.v` | Passes the new parameters through. The two SRAMs get separate depths (`BND_DEPTH`, `BALL_DEPTH`) in place of the shared `SRAM_DEPTH`. |

`hc_fifo.v`, `hc_collect.v`, `hc_bnd_check.v` and `ball_sram.v` are unchanged. The FIFO sideband
already carries the source detector id, which is what the lanes add `delta` to. Syndrome ingestion
and the edge emission order are unchanged:

- **Defect-defect edges:** `(u, w, v)` order, with `u < v`.
- **Boundary edges:** ascending `u`, on their own port.

## Table compiler (`tablegen/`, new)

The vector generator that `run_crosscheck.sh` called on `main` (`accelerator/tools/hc_dump.cc`) was
never committed, because the root `.gitignore` ignores every `tools` path. It is rewritten here, in
`tablegen/`.

- **`hc_dump.cc`** generates a rotated-surface-code memory-X corpus, runs `compile_ball_tables`
  (T = T_max = 2 lattice edge weights, R = 2T), and re-encodes the result.
  - **Outputs:** `ball.hex`, `bnd.hex`, `params.vh`, `syn_<k>.txt`, `ref_<k>.txt`, `model_cycles.txt`
    and `stats.txt`. The reference H comes from `build_ball_graph(..., SCAN)` on the
    **uncompressed** tables.
  - **Options:** `--d --rounds --p --shots --keep --seed --lanes --time-coord --no-prune --no-dedup
    --stats-only`.
  - **Cycle model:** a register-level replay of the control path. It covers the scan slots, the
    SRAM stage, the FIFO and its `cnt >= 3` stall, the frozen lane pipeline, the collector drain
    and the done FSM's quiet test. A closed-form sum `7 + clog2(cap+1) + fetch + drain` was not
    enough: it is one cycle off whenever trailing scan slots hide under a drain, for example a final
    defect with an empty forward row, which pruning makes common.
  - **Legacy layout:** the same model is also evaluated on the old table layout over the same
    shots. Its numbers reproduce `sweep_dist/` to within 0.3%.
- **`CMakeLists.txt`** builds `hc_dump` against the repository's own library. The root project is
  added `EXCLUDE_FROM_ALL`, so only `libspecmatching` and `libstim` are compiled.
- **`sweep.py`** runs `hc_dump --stats-only` over d ∈ {17, 25, 31} × p ∈ {5e-4, 1e-3, 3e-3} and
  writes `sweep_tshift/stats.csv` and `sweep_tshift/table_sizes.md`. It uses the standard library
  only.

## Tests (`tests/`)

- **`run_crosscheck.sh`**
  - Builds `tablegen/` instead of the missing `tools/`. When `obj/_deps` checkouts exist it reuses
    them, so configuring needs no network access.
  - Passes any arguments through to `hc_dump`. The default is the checked-in vector set: the first
    20 shots of a d = 5, 5-round, p = 5e-3, 60-shot corpus, seed 7.
  - Re-verilates on every run, because the port widths come from `params.vh`.
  - Adds a sixth step that diffs the RTL's per-shot cycle counts against the model.
- **`sim_main.cc`**
  - Writes `dut_cycles.txt`.
  - `iM` is assigned at its Verilated width; it used to be cast to `uint8_t`, which truncates
    syndromes of more than 255 defects.
  - The watchdog is raised to 200,000 cycles, for d = 17, p = 3e-3 shots.
- **`compare.py`** is unchanged.
- **`vectors/`** is regenerated in the new format.
  - Shots 0–17 have the same syndromes and reference H as on `main`.
  - Shots 18–19 differ. The old generator skipped one shot of the 60-shot draw for a reason it did
    not record, and this one does not.
  - New files: `model_cycles.txt`, `dut_cycles.txt`, `stats.txt`.

## Verification

Every configuration below matches on H and on cycles for every shot:

| corpus | tables | lanes |
|---|---|---|
| d = 5, 5 rounds, p = 5e-3 (checked-in vectors) | default | 64 |
| d = 5, 15 rounds, p = 5e-3, 60 shots | default, `--no-prune`, `--no-dedup --no-prune` | 64 |
| d = 7, 21 rounds, p = 3e-3, 60 shots | default | 16, 128 |
| d = 9, 27 rounds, p = 2e-3, 60 shots | default | 32 |
| d = 9, 27 rounds, p = 5e-4, 60 shots (17-bit weights) | default | 64 |
| d = 17, 17 rounds, p = 1e-3 / 3e-3, 100 / 40 shots | default | 64 |

The largest case, d = 17 at p = 3e-3, has up to 271 defects per shot.

## Results and caveats

- **Where the latency gain comes from:** pruning, not deduplication. Pruning halves the chunks
  fetched, but it packs the same edges into fewer batches, so the collector's drain grows. The
  speedup is 1.4–1.6× at L = 64 for p ≤ 1e-3, but only 1.14× at p = 3e-3, where the
  one-edge-per-cycle collector bounds the latency.
- **Wasted word width:** the ball word is still sized by the longest row (145–158 entries, i.e. 3
  chunks at L = 64, down from 4). Shorter rows leave the rest of the word unused.
- **Boundary table:** it is the one table indexed by detector id, so it stays O(d³). It is about
  40 bits per detector (1.1 Mbit at d = 31), roughly 2% of the tables.
- **Not yet re-done:** area and power have not been re-synthesized or re-run through CACTI for the
  new tables.
