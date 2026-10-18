# H-Build Latency Table from the 10000-Shot Distribution Run (bsearch lane)

Built from `sweep_dist/` ONLY: every cell below is an aggregate of
`sweep_dist/dist_stats.csv`, i.e. the 10000-shot distribution run of
`latency_dist_bsearch.py` (seed 20260903) for the bsearch lane at per-shot
capacity. No number is taken from the 1000-shot sweeps (`latency_table.md`,
`sweep_bsearch/`) or their frozen TSVs; the CAM twin of this run lives in
`sweep_dist_cam/` and is deliberately not repeated here.

Model per shot (see `bsearch/latency_table.md` for derivation and the RTL
anchors): cycles = 7 + clog2(max(M,1)+1) + t_fetch + (edges - H). Both
validation gates passed on all 9 points (`sweep_dist/run.log`); the p = 3e-3
bundles are the `bundles_sweep_dist/` set (provenance in `latency_dist.md`).

Caveats carried over: the capacity caveat (a deployed fixed-capacity instance
pays a constant clog2(cap+1) - 1 every shot instead of the per-shot term),
the batch-size caveat (this 10000-shot batch is its own seeded draw), and the
grid here is p = 5e-4/1e-3/3e-3, so the two RTL-confirmed p = 1e-4 cells of
the 1000-shot table do not appear; no cell below is RTL-measured. No
zero-defect shots occurred at any point.

Full per-shot distributions: `sweep_dist/dist_<point>.csv`; median, p95 and
max per cell: `sweep_dist/dist_stats.csv`; CCDFs: `sweep_dist/latency_dist.png`.

---

## Latency (cycles), 10000 shots per cell, model-derived

### p = 0.0005

Mean +- std (sample, ddof = 1):

| d \ LANES | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| 17 | 438.974 +- 118.74 | 239.403 +- 64.74 | 143.194 +- 39.73 | 98.049 +- 29.28 |
| 25 | 1528.009 +- 234.44 | 815.883 +- 127.02 | 471.467 +- 77.03 | 302.672 +- 54.80 |
| 31 | 3044.820 +- 334.79 | 1617.545 +- 181.11 | 925.276 +- 109.21 | 580.345 +- 76.71 |

Tails, p99 / p99.9:

| d \ LANES | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| 17 | 746.0 / 879.0 | 410.0 / 481 | 251 / 303.0 | 181 / 221 |
| 25 | 2111 / 2317.0 | 1135.0 / 1241.0 | 666 / 734.0 | 443 / 499.0 |
| 31 | 3865 / 4175.0 | 2061 / 2241.0 | 1197.0 / 1293.0 | 774.0 / 843.0 |

### p = 0.001

Mean +- std (sample, ddof = 1):

| d \ LANES | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| 17 | 881.962 +- 172.11 | 488.420 +- 99.39 | 301.404 +- 68.04 | 221.119 +- 57.39 |
| 25 | 3078.384 +- 334.66 | 1682.147 +- 192.13 | 1016.240 +- 130.44 | 718.586 +- 108.61 |
| 31 | 6122.521 +- 486.81 | 3333.062 +- 278.85 | 2000.848 +- 188.62 | 1394.486 +- 155.60 |

Tails, p99 / p99.9:

| d \ LANES | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| 17 | 1310 / 1508 | 742 / 875.0 | 481.0 / 576.0 | 380 / 463.0 |
| 25 | 3882.1 / 4189.0 | 2150.0 / 2327.0 | 1346.0 / 1468 | 1000.0 / 1106 |
| 31 | 7282.0 / 7707.0 | 4009.0 / 4262.0 | 2470.0 / 2639.0 | 1788.0 / 1951.0 |

### p = 0.003

Mean +- std (sample, ddof = 1):

| d \ LANES | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| 17 | 2751.448 +- 330.85 | 1687.980 +- 239.40 | 1258.018 +- 215.93 | 1117.067 +- 213.85 |
| 25 | 9624.896 +- 651.96 | 5879.937 +- 472.72 | 4370.541 +- 428.03 | 3862.970 +- 425.64 |
| 31 | 19091.768 +- 931.98 | 11649.327 +- 676.85 | 8653.130 +- 614.86 | 7635.033 +- 612.88 |

Tails, p99 / p99.9:

| d \ LANES | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| 17 | 3550.0 / 3951.0 | 2284.0 / 2599.0 | 1816.0 / 2094.0 | 1676.0 / 1946.0 |
| 25 | 11144 / 11748.0 | 7004.0 / 7459.1 | 5408.0 / 5823.0 | 4909.0 / 5332.0 |
| 31 | 21296.0 / 22001.0 | 13261.0 / 13768.0 | 10147.0 / 10597.0 | 9138.0 / 9568.0 |
