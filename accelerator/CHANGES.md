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
| 17 | 4,896 | 2,414 | 36.5 → 13.2 Mbit | 2.8× | 328 → 235 |
| 25 | 15,600 | 5,158 | 120.2 → 29.3 Mbit | 4.1× | 1,105 → 767 |
| 31 | 29,760 | 7,888 | 236.9 → 46.6 Mbit | 5.1× | 2,171 → 1,486 |

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
| `hc_lane.v` | An entry is `{w, delta}`. The lane computes `target = iDet + sext(delta)`, which replaces the `e_tgt > iDet` comparator with an adder of the same width. The filter becomes `!delta[MSB]` (delta is never 0). The binary search keeps its one-probe-per-stage shape, but each level now registers the candidate *pair* its successor will probe (`p_pair`), so a probe is a 2:1 pick (`csel`) in front of the compare instead of a `2**gs : 1` mux off the path. Levels past `MUX_T` narrow their pair across several precompute stages (`g_look`), and `p_path_b`/`p_path_c`, `p_go_b`/`p_go_c`, `p_hot_b`/`p_hot_c` hold both polarities of those selects. Pipeline depth is unchanged. |
| `hc_fifo.v` | Pipelined. The chunk cut out of `iBallRdata` has its select replicated over `SEL_COPIES = 8` slices and is registered (`chunk_q`, `rd_sb2`, `rd_v2`), so two push stages sit between the SRAM read and the FIFO write. The FIFO reads out through a first-word-fall-through output register (`o_v`/`o_d`), and `mem_empty`/`stall_q` are registered off `cnt_next` instead of being decoded from `cnt`. `FIFO_HIGH` is `FIFO_DEPTH - 2`, so `oStall` rises two slots early, since two pushes are in flight when it does. The `quiet` test counts both push stages and a batch still held in the output register. |
| `tb/hc_top.v` | Passes the new parameters through. The two SRAMs get separate depths (`BND_DEPTH`, `BALL_DEPTH`) in place of the shared `SRAM_DEPTH`. |

`hc_collect.v`, `hc_bnd_check.v` and `ball_sram.v` are unchanged. The FIFO sideband already carries
the source detector id, which is what the lanes add `delta` to. Syndrome ingestion and the edge
emission order are unchanged:

- **Defect-defect edges:** `(u, w, v)` order, with `u < v`.
- **Boundary edges:** ascending `u`, on their own port.

The pipelining in `hc_fifo.v` and `hc_lane.v` is there for 1 GHz timing closure in ASAP7;
`hc_max_frequency.md` reports the post-route slack and the f_max it implies, with `hc_lane` on the
critical path. The commit that introduced it, 72ffb71, dropped the `reg rd_v;` declaration together
with the comment block around it, leaving `rd_v` undeclared where `hc_fifo.v` still drives and reads
it; that one line is restored.

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
    --stats-only --dist-csv`.
  - **`--dist-csv <path>`** writes one row per shot for the *legacy* layout at per-shot capacity:
    `shot,n_defects,edges`, then `cycles_L<L>`, `fetch_L<L>` and `drain_L<L>` for each lane count.
    It is the per-shot form of the aggregated `legacy_*` fields and changes no other output.
  - **Cycle model:** a register-level replay of the control path. It covers the scan slots, the two
    cycles between the SRAM read and the FIFO write (the sideband stage and the registered chunk
    cut), the FIFO with its `cnt >= 2` stall and its registered output slot, the frozen lane
    pipeline, the collector drain and the done FSM's quiet test. A closed-form sum
    `9 + clog2(cap+1) + fetch + drain` is not enough: it is one cycle off whenever trailing scan
    slots hide under a drain, for example a final defect with an empty forward row, which pruning
    makes common.
  - **Legacy layout:** the same model is also evaluated on the old table layout over the same
    shots, so both sides of every latency comparison are the same pipelined RTL.
- **`CMakeLists.txt`** builds `hc_dump` against the repository's own library. The root project is
  added `EXCLUDE_FROM_ALL`, so only `libspecmatching` and `libstim` are compiled.
- **`sweep.py`** runs `hc_dump --stats-only` over d ∈ {17, 25, 31} × p ∈ {5e-4, 1e-3, 3e-3} and
  writes `sweep_tshift/stats.csv` and `sweep_tshift/table_sizes.md`. It uses the standard library
  only.
- **`sweep_dist.py`** runs `hc_dump --stats-only --dist-csv` over the same nine points and writes
  `sweep_dist/`: a per-shot `dist_d<d>_p<p>.csv` per point, then `dist_stats.csv` (mean, median,
  p95, p99, p99.9 and max per (d, p, lane count)) and `latency_table_dist.md`, both derived from
  those CSVs and from nothing else. Standard library only.
- **`plot_dist.py`** draws the CCDF grid from the same per-shot CSVs, one panel per (d, p) and one
  curve per lane count, to `sweep_dist/latency_dist.png` and `.pdf`. It needs matplotlib.

## Latency distribution (`sweep_dist/`)

`sweep_dist/` is a fresh run: nine (d, p) points × 10,000 shots, seed 20260903, the **legacy** ball
layout under the pipelined RTL, at per-shot lane-pipeline capacity `clog2(M + 1)`, from the
register-level replay that `tests/run_crosscheck.sh` checks cycle-exact. It is produced by
`tablegen/sweep_dist.py` on `hc_dump --dist-csv`, with the CCDFs from `tablegen/plot_dist.py`; every
cell of `latency_table_dist.md` is an aggregate of the per-shot CSVs in the same directory.

It supersedes the October 2026 run that previously filled this directory and shares no numbers with
it: that one used the closed form `7 + clog2(max(M,1)+1) + t_fetch + (edges - H)` against the
pre-pipeline RTL, from a generator (`latency_dist_bsearch.py`) that is in no tree of this repository
and whose shot draw `hc_dump` does not reproduce. It is the same model on the same grid as the
`legacy_*` columns of `sweep_tshift/stats.csv`, and agrees with them to `hc_dump`'s printed
precision on all 36 (d, p, lane count) cells: the maxima are identical and the means match within
the six significant digits `hc_dump` prints. `table_sizes.md` renders those means to one decimal, so
a cell can disagree in its last digit — d = 25, p = 1e-3, L = 32 reads 1786.5 there, against an
exact mean of 1786.4467. What this run adds is the distribution the headline table cannot carry:
median, p95, p99, p99.9 and max per cell.

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
  speedup is 1.4–1.6× at L = 64 for p ≤ 1e-3, but only 1.19× at p = 3e-3, where the
  one-edge-per-cycle collector bounds the latency.
- **What the pipelining costs in cycles:** +2 per shot wherever the design never stalls — every
  checked-in d = 5 vector moves by exactly +2 — and up to ~10% where the FIFO stall bounds the run
  (the small-L, p = 3e-3 points), because `FIFO_HIGH` tightened by one slot. Both sides of every
  legacy comparison are the pipelined RTL, so the table-size and speedup story is unaffected.
- **Wasted word width:** the ball word is still sized by the longest row (145–158 entries, i.e. 3
  chunks at L = 64, down from 4). Shorter rows leave the rest of the word unused.
- **Boundary table:** it is the one table indexed by detector id, so it stays O(d³). It is about
  40 bits per detector (1.1 Mbit at d = 31), roughly 2% of the tables.
- **Not yet re-done:** area and power have not been re-synthesized or re-run through CACTI for the
  new tables.
