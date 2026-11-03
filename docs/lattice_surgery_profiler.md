# `lattice_surgery_profiler` — design

Status: **implemented** in `benchmarks/spec_matching/lattice_surgery_profiler.cc` (generator:
`lattice_surgery_circuit.h`; analysis: `analyze_lattice_surgery.py`, `plot_lattice_surgery.py`). No
campaign has been run. Where the implementation departs from this design — most importantly, slab
defects are phase-1 defects rather than obstacles (§5.1 as written cannot hold `Y(a) <= d(a, b) - T`
for a cube defect inside a neighbour's ball), and `H`-only defects are re-solved in the speculative
`H` instance, which makes §10 items 3 and 5 unnecessary — is listed as `deviation_*` in the binary's
header and in every `run.log`. This document records the
design discussion behind the lattice-surgery experiment for the evaluation section, and specifies
the profiler binary that measures it. It replaces the cascaded-fallback sketch of the paper's
§4.6 (`arxiv/main.tex`, "Lattice surgery"), whose problems are listed in §1.

Contents

- [0. Summary](#s0)
- [1. Why §4.6 of the paper is replaced](#s1)
- [2. The experiment](#s2)
- [3. The modelled system](#s3)
- [4. The no-fallback path: per-cube solve and depth-1 seam fusion](#s4)
- [5. The fallback path: fusing a G state with an H state](#s5)
- [6. Exactness](#s6)
- [7. Implementation without enumerating G or H](#s7)
- [8. The binary: execution model, CLI, outputs](#s8)
- [9. Verification and offline checks](#s9)
- [10. Required code changes](#s10)
- [11. Alternatives considered and rejected](#s11)
- [12. Open questions and risks](#s12)
- [Appendix A. Notation](#sA)

---

## <a name="s0"></a>0. Summary

A `1×N` strip of distance-`d` surface-code patches is merged by lattice surgery and decoded for `d`
rounds, which gives a detector graph `G` of size `O(d × Nd × d)`. The strip is cut into `N`
spacetime **cubes** `G_1 … G_N`, one per patch, each `d × d × d`, meeting at `N − 1` spacelike
**seams**. Each cube has its own SpecMatching instance: an accelerator, a manager core, `k` solver
cores and a fallback core. A **global manager** core schedules the per-cube managers and the seam
fusions.

The two questions the experiment answers:

1. **No fallback anywhere.** Can the `N` speculative results be fused exactly across the seams,
   and at what latency? Answer by design: yes, with the existing fusion machinery, and in
   **constant depth** — every seam fuses at once, because blossom on `H` splits exactly into
   connected components, and two seam fusions conflict only when one component of `H` spans a
   whole cube (§4).
2. **Some cube falls back.** Can its fallback's `G` state be fused with the neighbours' `H` (or
   `G`) states **locally**, without escalating the whole strip and without enumerating `G` or `H`?
   Answer by design: yes. The fallback is given a seam condition (**obstacles**) that makes its
   state compatible with *any* neighbour `H` state. A hybrid resume then runs on a **full-resolution
   territory `Ω`** that starts as the cube plus a slab of width `W` and grows only on demand (§5).

The second point is the one that matters at scale. A real computation has billions of lattice
surgery operations. If one truncated component forced the whole computation to escalate, the
system would fail. The design keeps every `H → G` transition **rare, at the per-cube fallback rate
`q`, and local** to the cube that fell back plus its seam slabs, apart from extensions of `Ω` that
the experiment measures.

---

## <a name="s1"></a>1. Why §4.6 of the paper is replaced

The review of §4.6 as it stood found these problems. They motivate the design below.

1. **The fallback cannot be fused.** The fallback is stock pymatching
   (`spec_matching/escalation/escalate.cc:48`, `pm::decode_detection_events`), which returns
   observables and a weight, but no primal-dual state. §2.2 of the paper itself says sparse
   blossom's clock-based state cannot be handed over as an intermediate state. Fusing fallback
   results needs a fallback that exports its state (§10).
2. **`G` and `H` are conflated.** "Edges of `G` that cross a seam" have detector endpoints, which
   are usually not defects, but the Step 1 dummy mechanism lives on `H`, whose nodes are defects.
   The correct speculative seam cut is the set of **`H`-edges** across the seam: defect pairs with
   `d_G ≤ 2T`. Computing them needs syndrome from both sides of the seam.
3. **Ball tables must be compiled on the full `G`.** Per-cube tables would give `d_{G_i} ≥ d_G`
   when a shortest path passes through a neighbour cube, and the §4.5 certificate would break.
4. **The cascade is contradictory.** The stopping rule "until the union of the duals is feasible"
   is the wrong criterion, because fusion keeps feasibility by construction; what fails is the
   `Y ≤ T` cap. As written, the signal climbs every ancestor up to the root, which contradicts
   "only the seams adjacent to a truncated block pay".
5. **"Internal nodes are seams" is wrong for cyclic block layouts.** One fusion can close several
   seams at once, as with the quotient multiplicities `c_PQ` of Step 2.
6. **"Regions far from it are untouched" is empirical, not guaranteed.** The same is true in
   Fusion Blossom.
7. **"Rare, therefore localized" is a non sequitur.** Table 1's rates are for `d`-round memory
   experiments, not for merged blocks.

Cross-reference: §5.3 of the paper ("fuses adjacent blocks exactly, with work local to the seam")
depends on §4.6, and must be updated together with it once this experiment has numbers.

---

## <a name="s2"></a>2. The experiment

### 2.1 Circuit

- **Merged patch.** A rotated surface code on a `d × L` rectangle, with `L = N·d` (`--seam-gap 0`)
  or `L = N·d + (N − 1)` (`--seam-gap 1`, one routing column between neighbouring patches), run for
  `d` rounds.
- **Merge initialisation** (`--merge-init`, default on): the seam stabilizers, the ones that exist
  only in the merged patch, get **no detector in their first round**, because their outcomes are
  random. The product of their first-round outcomes is the joint logical measurement, recorded as
  an observable. With `--merge-init off`, the circuit is a plain rectangular memory experiment,
  kept as a control.
- **Orientation** (`--orientation long|short`). It decides which boundaries the decoded detector
  type sees.
  - `long`: that type's boundaries run along the long sides, so `b_G(u) ≲ d/2`.
  - `short`: its boundaries are the two short ends, so `b_G(u)` can reach `≈ Nd/2`.

  The locality claims of §5 are much stronger for `long` (§12), so both orientations must be run.
- **Generator.** `stim::generate_surface_code_circuit` makes square patches only. A custom
  generator `benchmarks/spec_matching/lattice_surgery_circuit.h` is needed, built with the same
  noise channels as `profiler_util.h`'s `generate` (after-Clifford depolarization,
  after-reset flip and before-measure flip, all at `p`), so that numbers are comparable with the
  memory-experiment profilers. `run.log` records the exact generator call.

### 2.2 Cubes and seams

- **Cube assignment.** Detector `v` belongs to cube `i = cube_of(v)`, determined by its spatial
  coordinate along the strip. With `--seam-gap 1`, detectors in a routing column are assigned to
  the cube on their left. `run.log` records the rule and the detector count of each cube.
- **Seams.** Seam `s_i` lies between cubes `i` and `i+1`. The **face** `F_i` of cube `i` toward
  seam `s_i` is the set of cube-`i` detectors with a `G`-edge into cube `i+1`.
- **Timelike seams** (successive windows of rounds) are **out of scope**: there is a single window
  of `d` rounds. §12 notes what changes for streaming.

### 2.3 Sweep

| parameter | values for the paper | flag |
|---|---|---|
| patches `N` | 2, 4, 8 | `--N` |
| solver cores per patch `k` | 4, 7 | `--k` |
| distance `d` | 17, 25 | `--d` |
| physical error rate `p` | `5e-4`, `1e-3`, `3e-3` | `--p` |
| cutoff `T` | 2 | `--T` |
| slab width `W` | `T`, `2T` | `--W` |
| orientation | long, short | `--orientation` |

`T` and `W` are in **multiples of `w_p = ln((1−p)/p)`**, as `--T` is in every other binary in this
directory. `run.log` records `T`, `T_weight_units`, `T_int`, and the same three fields for `W`.

### 2.4 What is reported

1. **Seam-fusion depth.** The fraction of shots, and of seams, where a component of `H` spans a
   whole cube and forces sequential seam fusion (§4.3).
2. **Latency** of the no-fallback path, and of the full system with local fallback resolution,
   against the baselines in §8.6. Reported per shot and per round (divided by `d`).
3. **Fallback statistics.** The per-cube truncation rate against `q` from Table 1, and the
   strip-level rate of "any cube truncated", which is roughly `N·q`.
4. **Locality of resolution.** Size of the final `Ω` (in detectors, and in cubes touched), number
   of `(t1)`/`(t2)` extensions, number of reset `H`-components and reset defects, and how often `Ω`
   leaves the preloaded slab (`handoff`, §7.3).
5. **Exactness.** Under `--verify`, every shot's matching weight equals stock pymatching on the
   full strip `G`.

---

## <a name="s3"></a>3. The modelled system

```
                             ┌───────────────────────┐
                             │    global manager     │  seam job graph, seam-fusion placement,
                             │      (1 core)         │  phase-2 orchestration, combine + certify
                             └──────────┬────────────┘
           ┌────────────────────────────┼────────────────────────────┐
   ┌───────┴────────┐           ┌───────┴────────┐           ┌───────┴────────┐
   │ cube 1         │  seam s_1 │ cube 2         │  seam s_2 │ cube N         │
   │ accelerator    │◄─────────►│ accelerator    │◄── … ────►│ accelerator    │
   │ manager core   │           │ manager core   │           │ manager core   │
   │ k solver cores │           │ k solver cores │           │ k solver cores │
   │ fallback core  │           │ fallback core  │           │ fallback core  │
   └────────────────┘           └────────────────┘           └────────────────┘
```

Total cores: `N·(k + 2) + 1`.

**Accelerator (assumed, not designed).** Each cube's accelerator does what it does now: it emits
the `H`-edges between defects of its own cube. It can also be extended, easily, to emit:

- the **cross-seam `H`-edges**: defect pairs `(a, b)`, with `a` in cube `i` and `b` in cube `i±1`,
  and `d_G(a, b) ≤ 2T`. Both accelerators of a seam emit the same edge; by convention the cube on
  the left owns it, and the global manager deduplicates;
- the neighbour's **syndrome bits inside the preloaded slab** (§7.3), forwarded to the fallback
  core for obstacle setup (§5.1). This is a few bits per detector.

The profiler stands in for the accelerators with ball tables compiled **on the full strip `G`**
(§1.3). It splits the edge list into per-cube and per-seam lists, outside every timed region. The
accelerator's latency is modelled as a constant per cube (`--accel-ns`, default 0; set it from the
preprocessor numbers of Table 2 to include it).

**Per-cube manager.** Runs Steps 1–5 of §4.3 of the paper on its own cube's `H_i`. It treats the
cross-seam `H`-edges exactly as Step 1 crossing edges: masked, with a dummy of `floor(w_min/2)`
rounded down to even on each endpoint, which is `fusion::even_dummy_weight`.

**Fallback core.** Holds `G_i` plus a **preloaded slab** of each neighbour (§7.3). It **never holds
ball tables**: the profiler must enforce that phase 1 and phase 2 code paths take only `G`. From
the start of every shot it runs phase 1 (§5.1), on every cube, whether or not the speculative path
will truncate. This is today's fallback timing.

**Global manager.** Keeps the seam job graph, runs the depth-1 seam conflict check (§4.3), places
seam fusions, orchestrates phase 2 for each truncated cube (§5.2), and combines the observables and
weights of every root into the strip's correction.

---

## <a name="s4"></a>4. The no-fallback path: per-cube solve and depth-1 seam fusion

### 4.1 Per-cube solve

This is unchanged from `fusion_lb_profiler`. Bounded union-find cut, fusion tree, leaves, fusions
and extraction run within each cube, on that cube's `k` solver cores. The cross-seam edges stay
masked throughout.

### 4.2 The `H`-components across the strip

Blossom on `H` splits **exactly** into independent connected components of the full strip's `H`,
counting all cross-seam edges. A seam fusion reopens only the components that have crossing edges
at that seam. Two consequences:

- Fusions at different seams **commute exactly** when no component touches both seams.
- PWD and modular decoding need two layers and buffers of width `~d` because they have no dual.
  Here the dummy inequality `b(u) + b(v) ≤ w` already makes the joined state feasible, so no buffer
  is needed and the result stays exact.

### 4.3 Depth-1 schedule and the conflict check

- Each per-cube manager labels every local component `L` (touches seam `s_{i−1}`) and/or `R`
  (touches seam `s_i`). This is O(1) extra work per crossing edge, inside the union-find pass.
- The global manager joins the labels across seams and finds the **spanning chains**: the
  components that touch two or more seams. Only those need sequential seam fusions, in seam order.
  Every other seam component is fused **concurrently with all others, in one layer**.
- A seam fusion is an ordinary §4.4 fuse job. It unmasks that component's crossing edges, raises
  or removes the dummies, releases the regions they held, and resumes `T`-truncated blossom. It is
  placed on the solver core of the side whose leaf finished last, which is the existing rule,
  extended across cubes.
- Spanning a whole cube needs a chain of defects spaced `≤ 2T` apart across distance `d`. That is a
  percolation event, exponentially rare below threshold. The experiment measures it (§2.4.1).

The user's two-layer odd/even idea (fuse 1–2, 3–4, 5–6, then 2–3, 4–5) is valid, but strictly
weaker. With component labels, all seams go in one layer, and only spanning chains are sequenced.

### 4.4 State placement and communication

Each cube's solver state lives in its own arena. A seam fusion reads the seam components' state
from both cubes. The profiler models the transfer cost as `--link-ns` per transferred defect
(default 0) and logs the number of defects moved, so that a link model can be applied afterwards.

### 4.5 Truncation is per component

Every job of the per-cube solve, and every seam fusion, works on a single component of `H`.
So `TRUNCATED` is a **per-component** event. The other components of the cube still complete, and
the rest of the strip still fuses. Two consequences for the implementation:

- The truncated component's regions are **evicted** from the shared per-cube instance, so that later
  jobs on that instance never process its events. Today's `fusion_lb_profiler` abandons the shot
  at the first truncation (its deviation 5), and that has to change here.
- Its `H` state is discarded. Phase 2 replaces it with the fallback's `G` state.

---

## <a name="s5"></a>5. The fallback path: fusing a G state with an H state

Cube `i` is **truncated** if any job touching one of its defects reports `TRUNCATED`. Its result
then comes from its fallback core. The task is to fuse that `G` state with the states of cubes
`i ± 1`, which may be `H` states or, if those cubes also fell back, `G` states, **cheaply and
locally**.

### 5.1 Phase 1: fallback with obstacles (concurrent, from the start of the shot)

**Graph.** The fallback core of cube `i` solves on `G_i ∪ B_i`, where `B_i` is the slab of the
neighbour cubes within graph distance `W` of the faces of cube `i`.

**Obstacles.** The fallback does not know the neighbour's `H` solution; it only knows that *every*
`H` state has `Y(b) ≤ T`. So it treats each neighbour bubble as solid and at its maximum size.
Every neighbour defect `b` whose `T`-ball reaches the slab, meaning those within `W + T` of the
face, becomes a frozen ball of radius `T`:

- Run one Dijkstra from `b` on `G`, out to radius `T` plus one edge.
- For each **shell node** `x` (one with `d(x, b) ≥ T` and a neighbour inside the ball), add a
  virtual-boundary edge with weight `d(x, b) − T`, rounded down to even (§7.4).

**Outer edge.** The outer layer of `B_i` is a virtual boundary, as in Fusion Blossom.

**Guarantee.** Stock pymatching on this graph gives, for every cube-`i` defect `a` and every
neighbour defect `b`:

```
Y(a) ≤ d(a, b) − T
```

So the result is compatible with whatever the neighbour's `H` solve produced.

**Cost.** One ball of radius `T` plus one edge per obstacle, a few thousand node visits in total,
**off the critical path**. The fallback still starts at the beginning of the shot, as today. The
modification to pymatching is per-shot boundary edges plus state export (§10).

**The `G`–`G` case comes for free.** If both cubes `i` and `i+1` fell back, each side satisfies
`Y ≤ d − T` against the other, so `Y(a) + Y(b) ≤ d(a, b) − 2T`. Their union is feasible with `T`
to spare, and the fusion is Fusion Blossom with ball obstacles in place of virtual-boundary
vertices. The fallback never needs to know in advance whether its neighbour will fall back.

### 5.2 Phase 2: hybrid resume (after phase 1 of cube `i` finishes)

**Territory `Ω`.** Intuitively, `Ω` is the region where the decoder has **switched to full
resolution**. Inside `Ω` it works on `G`: explicit detectors and edges, uncapped regions, and
stock pymatching geometry. Outside `Ω` everything stays in `H`: defects only, each region a bubble
of radius at most `T`, and distances from the ball tables (via the accelerators' edges).
Initially `Ω = G_i ∪ B_i`. It **only grows**, and only where the solve pushes against its edge.

**Slab width `W`.** `W` is the thickness of the initial slab, and the step size when `Ω` grows. It
sets how far into the neighbour the fallback can see in `G`.

- `W ≈ T` is the natural minimum, because the obstacle shells already need about `T` of the
  neighbour's `G`.
- A larger `W` makes `(t1)` fire less often, but phase 1 has more obstacles to set up. The sweep
  (§2.3) measures the trade-off.

It is Fusion Blossom's partition boundary, except the boundary can move outward.

**Hybrid graph `X`.**

- `G[Ω]`: detector nodes carrying the fallback's state.
- **`H`-only defects**: the neighbour defects outside `Ω`, with their `H` state, their `H`-edges
  among themselves, and the dummies of their own other seams.
- The cross-seam `H`-edges from the accelerators. These are redundant with the slab geometry but
  exact, so they are harmless.

All nodes are addressed by **detector id** (§7.2).

**Release.**

- **`G` side.** Unmatch the regions matched to obstacle shells or to the outer virtual boundary,
  and remove those edges. The regions keep their radii and become roots.
- **`H` side.** Unmatch the regions matched to the dummies of seam `s_i` (or `s_{i−1}`), and remove
  those dummy edges. The regions keep their radii and become roots.
- **Reset.** Every `H`-component that intersects `Ω`, or comes within `T` of it, is **reset**: its
  duals are set to zero and all its matches are removed. Its defects inside `Ω` become ordinary
  `G` nodes, and those outside stay `H`-only with `Y = 0`. Lowering duals is always feasible, so
  this is a valid intermediate state, and it avoids converting nested `H`-blossoms into `G`
  shells. These components sit at the seam and would have been reopened by the fusion anyway.
  This also covers any seam fusion already done speculatively between cube `i` and its neighbour:
  those components contain cube-`i` defects, so they intersect `Ω` and are reset.

**Closure.** Establish the invariant before resuming:

> **Invariant.** Every `H`-only defect is at distance **at least `T`** from `Ω`.

Find defects within `T` of `Ω` with a multi-source Dijkstra outward from `Ω`'s outer layer, out to
radius `T`, using the syndrome bitmap. Extend `Ω` around each one (by `W`), reset its component,
and repeat until none is left. Defects are sparse, so this terminates after a few local steps.

**Resume.** Run blossom on `X`. `G`-territory defects are uncapped, and `H`-only defects keep the
cap `T`. Two triggers extend `Ω`; both are O(1) flag checks in the flooder:

- **(t1)** A region claims a node on `Ω`'s outer layer. Extend `Ω` by `W` from the preloaded `G`
  around that node, then run the closure.
- **(t2)** An `H`-only defect hits its cap `T`. Extend `Ω` around it, then run the closure.

**Why the margin is `T`, with no extra edge.** Only outer-layer nodes have edges that leave `Ω`.
`(t1)` fires the moment such a node is claimed, before any part of an edge leaving `Ω` is covered,
so no region ever reaches past `Ω`. (An earlier draft used `T + w_max`; that was over-conservative.
`T` and `w_max` are in the same units, `w_max ≈ 1.39 w_p` for stim's depolarizing model at
`p = 1e-3`, so the margin was never `T + 1` either.)

**Several cubes.** If several nearby cubes truncated, run one phase 2 over the union of their `Ω`s.
A cube with two neighbours gets both slabs in one `Ω` (one resume, not two).

**Where it runs.** Phase 2 runs on the fallback core of cube `i`, inside its own pymatching
instance, which already holds the phase 1 state (§7.1). When several cubes share one phase 2, it
runs on the fallback core of the lowest-index truncated cube.

### 5.3 Pause instead of abort

The per-cube spec solve keeps running after a component truncates (§4.5), so the rest of the strip
completes and fuses at the `H` level while phase 1 is still running. Phase 2 starts when **both**
of these are true:

1. phase 1 of cube `i` has finished, and
2. every `H` job whose components phase 2 might reset has finished. In practice these are the seam
   components of seams `s_{i−1}` and `s_i`.

The global manager enforces this ordering with an ownership claim on each component. A component
claimed by phase 2 is not handed to a seam fusion, and phase 2 waits for any fusion already
running on it.

---

## <a name="s6"></a>6. Exactness

**Certificate.** The final dual is a feasible dual of the metric-closure LP on all defects of the
strip, and every match is tight. So the matching is a minimum-weight perfect matching on the full
strip's `G`, as in §4.5 of the paper. The three kinds of pair:

| pair | why feasible |
|---|---|
| both in `G` territory | exact blossom on `G[Ω]`, given geodesic convexity (§9.2) |
| both `H`-only | `Y ≤ T` on both sides, and pairs missing from `H` have `d > 2T` (§4.5 of the paper) |
| `a` in `Ω`, `b` `H`-only | let `x` be the outer-layer node where the `a–b` geodesic exits `Ω`. Then `Y(a) ≤ d(a, x)` (from `(t1)`) and `d(x, b) ≥ T ≥ Y(b)` (from the invariant), so `Y(a) + Y(b) ≤ d(a, b)` |

**Initial feasibility of phase 2.** Phase 1 gave `Y(a) ≤ d(a, b) − T ≤ d(a, b) − Y(b)`. Dummy
releases keep radii, and resets only lower duals.

**Termination.** `Ω` only grows. In the worst case it becomes all of `G`, and then the run is
stock blossom on `G`. No separate whole-strip safety net is needed for correctness. The fallback
pymatching on the full strip is run only as the `--verify` reference and the latency baseline.

**What is not guaranteed.** No exact decoder can guarantee locality in the worst case, because one
augmenting path can cross the whole computation. Here locality is probabilistic: an extension past
the neighbour cubes needs a chain of length about `d`, which is suppressed like a logical failure.
The experiment measures this distribution (§2.4.4). If a hard bound is ever required, cap `Ω` at
about `d` around the truncated cube and accept a PWD-style approximate result beyond it, at a rate
of about the logical error rate. The profiler counts how often such a cap would have been hit
(`omega_cubes > 3`) but never applies it.

---

## <a name="s7"></a>7. Implementation without enumerating G or H

The constraint: no per-shot work proportional to `|G_i|` or to `|E(G_i)|`. Sparse blossom already
satisfies this, because its work scales with the volume its regions explore. The pieces below keep
it that way.

### 7.1 Phase 2 runs in the fallback's own instance

The fallback core's pymatching instance is built **once**, over the preloaded graph (§7.3). Phase 1
leaves its state there. Phase 2 **edits that instance in place**:

- It removes the obstacle and outer-boundary edges (it zeroes their slots, below).
- It appends the `H`-only defects and their `H`-edges (a small, per-shot addition).
- It unmasks the slab nodes that join `Ω`.

Nothing is rebuilt, and the per-shot reset touches only nodes the flood touched, as stock
pymatching already does.

### 7.2 One index space

Every node in the fallback instance and in `X` is addressed by **detector id**. The `H` nodes of
the speculative path are positions in the sorted syndrome list, which map to detector ids through
the list itself.

- An `H`-only defect **inside** the preloaded slab uses its own detector node. Its `G`-edges stay
  masked, and its `H`-edges go in a per-shot overflow adjacency.
- An `H`-only defect **beyond** the preloaded slab uses a slot from a fixed pool of spare nodes.
  `run.log` records the pool size, and running out of slots is an error counted per shot.

### 7.3 The preloaded slab

The fallback core preloads `G_i` plus a slab `P_i ⊇ B_i` of each neighbour, with
`--preload-width`, default `max(2W, W + 2T)`. The slab nodes outside `Ω` are **masked** with
`GraphFlooder`'s existing `edge_mask`.

- **Extending `Ω`** is unmasking, plus recomputing the outer-layer flag for the newly unmasked
  nodes and their neighbours. Both are local.
- **Extending `Ω` past `P_i`** would be a hand-off to the neighbour's fallback core in hardware.
  The profiler keeps going on the full `G`, logs the event as a `handoff`, and charges
  `--link-ns` per node moved.

### 7.4 Obstacles as per-shot boundary edges

Every node of `P_i` is preallocated **one boundary-edge slot**, with infinite weight by default.
Phase 1 sets the slot weight on each obstacle shell node, and on the outer layer of `B_i`, and
phase 2 sets them back to infinite. A shell node that lies on several obstacle shells keeps the
**minimum**, because a smaller boundary weight only strengthens `Y(a) ≤ d(a, b) − T`.

Weights obey the solver's even-integer rule (`fusion_solve.h:58-71`): `d(x, b) − T` is rounded
**down** to even, which again only strengthens the inequality. `T_int` is even by construction
(`truncation/horizon.h`).

### 7.5 Triggers are flags

- `(t1)` is a per-node `outer_layer` bit, checked where the flooder claims an empty node.
- `(t2)` is the existing dual-cap check (`graph_flooder.cc:501`), applied to `H`-only nodes only.
- The closure Dijkstra runs only when `Ω` changes, and only from the new outer layer.

No coordinate grid is used anywhere (§11.2).

---

## <a name="s8"></a>8. The binary: execution model, CLI, outputs

### 8.1 Location and build

`benchmarks/spec_matching/lattice_surgery_profiler.cc`, with a `CMakeLists.txt` target alongside
`fusion_lb_profiler`, using the same three provenance defines (`SPECMATCHING_GIT_HASH`,
`SPECMATCHING_BUILD_FLAGS`, `SPECMATCHING_STIM_VERSION`). A `lattice_surgery_profiler_debug` twin
is built with `-UNDEBUG -g`, so the feasibility asserts of §9.1 are live. It is not in the python
wheel.

### 8.2 Execution model

This follows `fusion_lb_profiler`'s `fusion_k_core_model`:

- **One thread; no thread is ever spawned.** Every unit of work runs on the profiling thread and is
  timed with `hires_timer.h`'s thread-CPU clock.
- **A simulated list scheduler** places those measurements on `N·(k + 2) + 1` cores: per cube, one
  manager, `k` solvers and one fallback; plus the global manager. The per-shot latency is the
  makespan of that schedule.
- **Dependencies** are the job graph of §8.3. Each job's start time is the maximum of its
  dependencies' end times and its core's free time, plus `--link-ns` per transferred defect when
  the job reads another cube's state.
- **Fallbacks run every shot**, from `t = 0` on their own cores, as in the real design. Their
  measured time is always logged, even when the result is unused.
- **The same caveats as `fusion_lb_profiler`**: transfer cost is zero unless `--link-ns` is set;
  the warm-cache effect of running jobs back to back is not modelled; nothing outside a timed
  region is on the critical path. That covers edge-list construction (the accelerator stand-in),
  instance resets, the `--verify` reference and output writing.

### 8.3 Jobs

Every job is one row of `jobs_*.csv` (§8.5).

| `kind` | core | depends on | what is timed |
|---|---|---|---|
| `accel` | cube accelerator | — | constant `--accel-ns` (not measured) |
| `mgr_uf` | cube manager | `accel` | Step 1 (union-find, cut, `L`/`R` seam labels) |
| `mgr_tree` | cube manager | `mgr_uf` | Step 2 (fusion tree) |
| `mgr_scatter` | cube manager | `mgr_tree` | Step 3 (scatter) |
| `leaf` | cube solver | `mgr_scatter` | build and solve of one piece |
| `fuse` | cube solver | its two children | intra-cube fusion |
| `g_seam_check` | global manager | every cube's `mgr_uf` | joining labels, finding spanning chains |
| `seam_fuse` | solver of the last-finishing side | `g_seam_check` and the roots of both sides' components | one seam component's fusion |
| `extract` | the root's core | the root | reading off matches |
| `fb_phase1` | cube fallback | — (starts at `t = 0`) | obstacle setup and stock solve on `G_i ∪ B_i` |
| `fb_obstacles` | cube fallback | — | obstacle setup alone (also part of `fb_phase1`; logged separately) |
| `p2_release` | fallback of the truncated cube | `fb_phase1` and the §5.3 claims | release and reset |
| `p2_closure` | same | `p2_release`, or an extension | closure Dijkstra |
| `p2_resume` | same | `p2_closure` | blossom on `X` until done or a trigger fires |
| `p2_extend` | same | a trigger | unmasking, outer-layer flag update |
| `g_combine` | global manager | every `extract` and every phase 2 | XOR of observables, sum of weights |

One phase 2 alternates `p2_resume → p2_extend → p2_closure → p2_resume`. Each segment is its own
row, so the extension count and its cost can be read from the job log.

### 8.4 CLI

```
lattice_surgery_profiler --N 2,4,8 --k 4,7 --d 17,25 --p 1e-3 --T 2
                         [--W 2] [--orientation long,short] [--seam-gap 0|1]
                         [--merge-init on|off] [--shots S] [--warmup 1000] [--seed SEED]
                         [--preload-width W'] [--accel-ns 0] [--link-ns 0]
                         [--sort] [--verify] [--check-determinism] [--jobs-log all|sample|none]
                         [--out DIR]
```

- `--N`, `--k`, `--d`, `--p`, `--T`, `--W` and `--orientation` take comma-separated lists. The
  binary runs their Cartesian product, one **cell** per combination.
- `--shots` is the number of shots per cell (required, positive). `--warmup` shots are decoded and
  discarded first.
- `--sort` has the same meaning as in `fusion_lb_profiler` (sort each edge list by weight).
- `--jobs-log sample` writes job rows only for every shot that truncates, plus the first 1000
  shots, to keep `jobs_*.csv` bounded at `10^6` shots. The default is `all`.
- Validation follows `fusion_lb_profiler`: `d` is odd and at least 3, `p` is in `(0, 1)`, `T`,
  `W`, `N` and `k` are positive.

### 8.5 Outputs (all to `--out`)

**`run.log`**: one per invocation.

- Provenance: git hash, build flags, stim version, the full command line.
- The circuit generator call (§2.1), the cube assignment rule, and the detector count of each cube.
- `w_p`, `w_min`, `w_max` and the number of distinct weights of the DEM, and `T`, `W`, and the
  preload width in float and integer units.
- The modelled core count `N·(k + 2) + 1`, `--accel-ns`, `--link-ns`, and the spare-pool size
  (§7.2).
- The model choices: one thread, simulated schedule, zero transfer cost unless set, evicting
  truncated components.
- The result of the offline convexity check (§9.2) for each cube.

**`shots_N{N}_k{k}_d{d}_p{p}_T{T}_W{W}_{orientation}.csv`**: one row per shot.

| column | meaning |
|---|---|
| `shot` | index |
| `n_def`, `n_def_cube_{min,max}` | defects in the strip, and the min/max over cubes |
| `n_edges`, `n_cross_edges` | `H`-edges in total, and across seams |
| `n_seam_components` | components with at least one crossing edge |
| `n_spanning`, `max_span` | components touching two or more seams, and the most seams any one component touches |
| `seam_depth` | layers of seam fusion used: 1 unless a spanning chain exists |
| `n_cubes_truncated` | cubes with at least one truncated job |
| `n_trunc_components` | truncated components |
| `trunc_where` | `leaf`, `fuse` or `seam` for the first truncation (empty if none) |
| `omega_nodes`, `omega_cubes` | size of the final `Ω` in detectors, and the number of cubes it touches |
| `n_t1`, `n_t2`, `n_closure_steps` | extension counts |
| `n_reset_components`, `n_reset_defects` | phase 2 resets |
| `n_obstacles` | obstacles set up across all cubes |
| `n_handoff`, `handoff_nodes` | extensions past the preloaded slab |
| `n_released_g`, `n_released_h` | regions released on each side |
| `n_moved_defects` | defects transferred between cubes (for the link model) |
| `t_spec_ns` | makespan of the speculative path alone (no truncation case) |
| `t_fb_max_ns` | the slowest `fb_phase1` in the strip |
| `t_p2_ns` | from the start of the first `p2_release` to the end of the last phase-2 job |
| `t_sys_ns` | makespan to `g_combine` end: the system latency |
| `t_strip_pm_ns` | baseline B1 (§8.6) |
| `t_global_escalate_ns` | baseline B2 |
| `weight`, `ref_weight` | the matching weight, and the reference weight under `--verify` |
| `obs_xor` | the observable mask produced |

**`jobs_…csv`**: one row per scheduled job, with columns `shot, job_id, kind, cube, seam,
component, core, deps, start_ns, end_ns, measured_ns, n_items, status`. `n_items` is the job's size:
defects for a leaf, crossing edges for a fusion, nodes unmasked for `p2_extend`. `status` is one of
`COMPLETE`, `TRUNCATED`, `EVICTED` and `TRIGGER_T1` / `TRIGGER_T2`.

**`phase2_…csv`**: one row per phase-2 run (per group of truncated cubes per shot), with columns
`shot, cubes, start_ns, end_ns, omega_nodes_initial, omega_nodes_final, omega_cubes, n_t1, n_t2,
n_closure_steps, n_reset_components, n_reset_defects, n_handoff, max_Y_over_T, resume_events`.

**`agg.json`**: per cell, the means and 99.9th percentiles of every `t_*` column (taken over
non-truncated shots, truncated shots and all shots separately); the strip-level truncation rate and
the rate per cube, against the `q` the existing `sparse_graph_stats` binary measures on a `d × d × d`
memory experiment; histograms of `seam_depth`, `omega_cubes`, `n_t1 + n_t2` and `n_reset_defects`;
the counts of `handoff` events and of `omega_cubes > 3` (§6); and the `--verify` mismatch count.

### 8.6 Baselines

| name | what it is |
|---|---|
| **B1**, `strip_pm` | stock pymatching on the full strip's `G`, on one core. This is the exact reference. |
| **B2**, `global_escalate` | the paper's §4.5 semantics applied to the strip: if any cube truncates, the strip waits for B1. Otherwise it uses `t_spec`. |
| **B3**, `cube_pm_max` | the maximum of stock pymatching per cube, with no seam handling. It is not exact; it is a reference point for the parallel decoding of `N` independent cubes, and is reported only as a latency number. |

The headline comparison is `t_sys` against B1 and B2, as a function of `N`. The B2 gap is what local
resolution buys, and should grow with `N`, since B2's escalation rate grows like `N·q`.

---

## <a name="s9"></a>9. Verification and offline checks

### 9.1 Per shot (`--verify`)

- **Weight identity.** `weight == ref_weight`, where `ref_weight` is B1's weight. This is the
  exactness test, and it counts mismatches rather than aborting, as in `fusion_lb_profiler`. The
  observable mask is also compared with B1's, and a mismatch is counted separately, because
  degenerate matchings can share a weight.
- **Debug twin.** In `lattice_surgery_profiler_debug`, asserts check:
  - dual feasibility of the phase 2 initial state against every `(a, b)` pair within `2T + W`;
  - the invariant of §5.2 after every closure;
  - that no `H`-only defect exceeds `T`;
  - that every match is tight.
- **`--check-determinism`.** It decodes every shot twice and compares rows, ignoring timing columns.

### 9.2 Offline, once per DEM

- **Geodesic convexity.** For pairs of detectors in the face slab of each cube, compare the
  distance over `G_i ∪ P_i` with the distance over `G`, up to radius `2T + W`. This is expected to
  hold in the surface-code bulk. Cover the places where it could fail: the merge seam (first-round
  seam stabilizers have no detectors), the time boundaries, and hook edges that run diagonally to
  the face. Results go in `run.log`, and the binary refuses to run a cell that fails unless
  `--allow-nonconvex` is given.
- **Weights.** `w_min`, `w_max` and the distinct weights, logged in `run.log`. The obstacle Dijkstra
  radius depends on them (§5.1); the margin does not (§5.2).

---

## <a name="s10"></a>10. Required code changes

1. **Fallback that exports its state.** It is stock pymatching (`escalate.cc:48`) today. Phase 2
   needs the fallback's `pm::Mwpm` instance alive after the solve: its regions, matches and the
   flooder's queue. This is a new entry point in the spec driver that runs
   `process_timeline_until_completion` without the final `reset`, which is what
   `process_timeline_until_completion_or_report` already allows (`docs/upstream_edits.md`).
2. **Per-shot boundary-edge slots on preloaded slab nodes** (§7.4). This is a vendored-graph
   change, so it needs a row in `docs/upstream_edits.md`.
3. **An overflow adjacency for `H`-edges on detector nodes, and a spare node pool** (§7.2). This is
   also vendored, and needs a row as well.
4. **Outer-layer flag and `(t1)` hook** where the flooder claims an empty node (§7.5). Also
   vendored, also a row.
5. **Per-node dual cap applied only to `H`-only nodes** in the phase 2 instance. The cap is
   currently per instance (`dual_cap` in `graph_flooder.h`).
6. **Evicting a truncated component** from the shared per-cube instance (§4.5), replacing
   `fusion_lb_profiler`'s abandon-at-first-truncation (its deviation 5).
7. **Component reset**: set duals to zero and unmatch, for a set of defects in the speculative
   instance, and transfer their `H` adjacency into the phase 2 instance.
8. **The cross-seam edge stand-in**: ball tables on the full strip's `G`, and splitting the edge
   lists by cube and seam (§3).
9. **The `1×N` circuit generator** (§2.1).
10. **The scheduler** over `N·(k + 2) + 1` cores with the job graph of §8.3, extended from
    `fusion_lb_profiler`'s list scheduler.

Items 2–5 are the only edits to vendored files.

---

## <a name="s11"></a>11. Alternatives considered and rejected

1. **Grafting `G` of the truncated cube onto the neighbours' `H`, with bridge edges** from cube
   detectors to neighbour defects taken from the ball tables. It is correct with an adaptive graft,
   but it needs ball-table rows for non-defect detectors at run time, and the ball table does not
   fit in the fallback core's memory. Replaced by obstacles (§5.1) plus territory growth (§5.2).
2. **A coordinate-grid proximity guard.** Using `d_G(v, b) ≥ w_min·⌈cheb(v, b)/s⌉`, with
   `w_min ≈ 0.78 w_p`, a coordinate step `s = 4` and margin `T + w_max`, the guard covers about 5
   edges, which is roughly `11³ ≈ 1300` cells to scan for every visited node. That is too slow.
   (Stamping cells from each defect would make it O(1), but it is unnecessary given the
   outer-layer flag.)
3. **Promoting only the truncated component, by resetting it and re-solving it on `G` without
   waiting for the fallback.** It is cheaper on paper, but it duplicates work the fallback core is
   already doing, and it needs its own proximity guard. Kept as a possible later variant, for
   comparison against §5.
4. **Promoting all of `H` in the cube to `G` from the truncated `H` state.** Correct, but it costs
   about as much as pymatching on the cube, and needs `H`-blossoms converted into `G` shells.
5. **Fallback with Fusion-Blossom virtual-boundary vertices at the face.** It is not compatible with
   an arbitrary `H` state on the other side. A neighbour defect `b` close to the face can have
   `Y(b)` larger than `d(b, face)`, which breaks feasibility. Obstacles fix this by using the
   worst-case bubble.
6. **The odd/even two-layer seam schedule** (PWD style). It is valid, but dominated by the
   component-labelled depth-1 schedule (§4.3).

---

## <a name="s12"></a>12. Open questions and risks

1. **Orientation `short`.** The decoded type's nearest boundary can be `≈ Nd/2` away, so a
   truncated region's reach is not bounded by one cube. Expect larger `Ω` there; the experiment has
   to show how much larger.
2. **`(t2)` has no worst-case bound.** Expect it to be small, because truncated defects usually need
   only a little more than `T`.
3. **Phase 2 latency sits behind phase 1.** The fallback of a `d × d × d` cube is the same cost as
   today's fallback, but in the strip it is on the critical path for every truncated cube. Compare
   with alternative 3 of §11 if it dominates.
4. **Pool exhaustion and hand-offs** (§7.2–7.3) are counted, not handled in hardware terms.
5. **Streaming (timelike seams).** A seam to a window that has not yet arrived can't be fused. The
   rule would be to commit a component once it no longer touches a future seam. That needs a
   left-deep schedule in time, and is a separate experiment.
6. **Real hardware link latency** between cubes is unknown. `--link-ns` and `n_moved_defects` let
   it be applied afterwards.

---

## <a name="sA"></a>Appendix A. Notation

| symbol | meaning |
|---|---|
| `G`, `G_i` | detector graph of the full strip, and of cube `i` |
| `H`, `H_i` | sparse graph on defects at cutoff `T`, for the strip and for cube `i` |
| `s_i`, `F_i` | seam between cubes `i` and `i+1`, and cube `i`'s face toward it |
| `T` | dual cutoff, in units of `w_p = ln((1−p)/p)`; `T = 2` |
| `Y(u)` | nested dual sum of defect `u` (the radius explored from `u`) |
| `b(u)` | dummy boundary weight, `floor(w_min(u)/2)` rounded down to even |
| `B_i` | initial slab of the neighbours, within `W` of the faces of cube `i` |
| `P_i` | preloaded slab, `P_i ⊇ B_i`, width `--preload-width` |
| `W` | slab width and extension step, in units of `w_p` |
| `Ω` | full-resolution (`G`) territory of a phase 2 run; starts as `G_i ∪ B_i`, only grows |
| obstacle | a neighbour defect treated, in phase 1, as a frozen ball of radius `T` |
| `(t1)`, `(t2)` | extension triggers: an outer-layer node is claimed; an `H`-only defect hits `T` |
| `q` | per-cube fallback (truncation) rate, from Table 1 of the paper |
