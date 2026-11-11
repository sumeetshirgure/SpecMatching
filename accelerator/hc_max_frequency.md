# Maximum frequency (based on hc_lane)

Synthesis was run with Cadence Genus and place-and-route with Innovus, in ASAP7 at a 1 ns target clock (1 GHz). Only the RTL was changed; the constraints are the flow's stock template SDC, which defines the clock on `iClk` only. Every slack below comes from the final post-route optimization report (`postRouteOpt`).

The maximum frequency is based on **hc_lane**. Its critical path starts at the enable input `iEn`, which fans out to every pipeline flop.

## Post-route slack at 1 GHz

| Block   | Depth  | WNS (ps) | Setup vio | max_cap | max_tran |
|---------|--------|---------:|----------:|--------:|---------:|
| hc_lane | d = 17 | +72      | 0         | 0       | 0        |
| hc_lane | d = 25 | -39      | 6         | 0       | 0        |

## Maximum frequency (hc_lane)

T_crit = T_clk - WNS, and f_max = 1 / T_crit.

| Depth   | NUM_SYND | WNS (ps) | T_crit (ps) | f_max (GHz) |
|---------|---------:|---------:|------------:|------------:|
| d = 17  | 123      | +72      | 928         | 1.08        |
| d = 25  | 340      | -39      | 1039        | 0.96        |
| d = 31* | 590      | –        | 1099–1167   | 0.86–0.91   |

\* Extrapolated from d = 17 and d = 25, not placed and routed. The two bounds assume the critical path grows logarithmically (a buffer tree) or linearly (wire length) with NUM_SYND.

## Notes

- f_max = 1 / (T_clk - WNS) is a first-order estimate; a run at the relaxed clock would be optimized differently.
- The flow's worst-case and best-case timing views both use the ASAP7 RVT fast (FF) library, so these frequencies are optimistic for slow-corner silicon.
- The d = 17 lane result is from run `rtl12_hc_lane_d17`, the block version with an input register rank. The d = 25 result is from `best_hc_lane_d25`, built from `best_rtl/hc_lane.v`.
