// Copyright 2026 SpecMatching contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/// `lattice_surgery_profiler` — speculative decoding of a `1×N` lattice-surgery strip, cut into `N`
/// spacetime cubes, with depth-1 seam fusion and local fallback resolution.
///
/// The design is `docs/lattice_surgery_profiler.md`; section numbers below refer to it. In brief:
/// each cube runs `fusion_lb_profiler`'s divide-and-conquer solve on its own `H_i`, with the
/// cross-seam `H`-edges masked behind dummies (§4.1); a global manager fuses every seam component in
/// one layer and sequences only the spanning chains (§4.3); every cube's fallback core runs phase 1
/// (§5.1) from `t = 0`; and a cube whose solve truncates is resolved by a phase-2 hybrid resume
/// (§5.2) that starts from its fallback's `G` state and grows a full-resolution territory `Ω` only
/// where the solve pushes against its edge. The output is exact MWPM on the full strip's `G` on
/// every shot, which `--verify` checks against stock pymatching on the strip (baseline B1).
///
/// ## Execution model
///
/// **One thread; no thread is ever spawned.** A shot is executed in a fixed dependency-respecting
/// order on the profiling thread — per-cube managers, the global seam check, every leaf, every
/// intra-cube fusion, every seam fusion, every cube's phase 1, each phase-2 group, the extractions,
/// the combine — and every unit of work is timed with `hires_timer.h`'s thread-CPU clock. A list
/// scheduler then places those measurements on `N (k + 2) + 1` simulated cores (per cube: a manager,
/// `k` solvers and a fallback; plus the global manager) and the per-shot latency is the makespan of
/// that schedule. Every job's computation is independent of where and when it is scheduled, which
/// is what makes executing first and scheduling afterwards equivalent to scheduling online; the
/// scheduler's own dispatch work (the ready-heap pop and the core pick) is timed *as it runs*, and
/// charged to the dispatching manager's core, as in `fusion_lb_profiler`.
///
/// The caveats of `fusion_lb_profiler` all travel with these numbers: transfer cost is zero unless
/// `--link-ns` is set; the warm-cache effect of running jobs back to back is not modelled; nothing
/// outside a timed region is on the critical path (the accelerator stand-in, the instance resets,
/// the baselines' own bookkeeping and output writing).
///
/// ## Where this implementation departs from the design, and why
///
/// All of these are written to `run.log` as `deviation_*` lines.
///
///  1. **Slab defects are phase-1 defects, not obstacles (§5.1).** As written, every neighbour
///     defect within `W + T` of the face is a frozen ball. A cube-`i` defect `a` within `T` of such a
///     `b` sits *inside* `b`'s ball, where `Y(a) <= d(a, b) - T` cannot hold, and phase 2 would later
///     have to inject `b` onto a node `a`'s region already covers. That is the commonest seam event
///     there is — one error on a seam-crossing edge. So the defects *inside* the initial territory
///     `Ω_0 = G_i ∪ B_i` are solved in phase 1 like cube-`i` defects, the obstacles are the defects
///     *outside* `Ω_0` whose `T`-ball reaches it, and an obstacle whose ball contains a phase-1
///     defect is **promoted**: `Ω_0` grows by `W` around it and it becomes a phase-1 defect too.
///  2. **The outer layer of `Ω_0` is a weight-0 virtual boundary.** "As in Fusion Blossom" is read
///     as: a region that reaches an outer-layer node may match there at no further cost. That gives
///     `Y(a) <= d(a, x)` for every outer node `x`, which is the same inequality `(t1)` maintains in
///     phase 2, so the phase-1 state is a valid phase-2 starting point without further argument.
///  3. **`H`-only defects are re-solved in the speculative `H` instance, not the fallback's (§7.2).**
///     An `H`-only defect shares no edge with `G` territory — the cross `H`-edges are redundant by
///     §5.2's own argument, and leaving them out is what lets a reset component be torn down without
///     touching `G` — so the two solves interact only through `(t1)`/`(t2)`, and the `H` side can run
///     where its adjacency already lives. §10 items 3 and 5 (overflow adjacency, spare node pool,
///     per-node dual cap) are therefore not needed and not implemented; the time is still charged to
///     the phase-2 job on the fallback core, and `n_moved_defects` counts the transfer.
///  4. **Claims.** Phase 2 claims (resets) every component with a defect in `Ω` or within `T` of it
///     (§5.2), **and** every component the speculative path left unfinished — one with a truncated or
///     skipped job — which the design leaves implicit. A seam component that contains a truncated
///     cube's defects is not handed to a seam fusion (§5.3's claim), as if the truncation were known
///     before any of its seam fusions started; the work a real system might have wasted on it is not
///     charged.
///  5. **Groups of truncated cubes (§5.2, "several cubes").** Adjacent truncated cubes, and cubes
///     joined by an unfinished component, form one group, whose phase 2 runs in the lowest cube's
///     fallback instance from that cube's phase-1 state. The other cubes' phase-1 states are not
///     merged in — their territories overlap, and a merge of two `G` states is not defined here — so
///     their unfinished components enter through the general mechanism: claimed, re-solved in `H` at
///     `Y = 0`, and pulled into `G` by `(t2)` where they need more than `T`. Their own `fb_phase1`
///     rows are still timed and logged. Two groups whose territories or claims meet are a
///     **conflict**: counted (`p2_conflict`), and the shot takes B1's answer and latency.
///  6. **Injected defects at an odd clock.** A defect injected into `G` territory mid-resume starts
///     with radius 1 when the shared clock is odd, so that every growing region keeps an even
///     `y`-intercept (the parity rule of `fusion::even_dummy_weight`). It is injected only on a node
///     that no region covers, with slack of at least one edge, so the extra unit is always feasible.
///  7. **`T`, `W` and the preload width are in lattice edge weights** — the median discretised edge
///     of `G`, `profiler_util.h`'s `edge_weight_units` — as `--T` is in every other binary in this
///     directory. The design's `w_p = ln((1-p)/p)` is logged beside it.
///  8. **The schedule is a time-ordered list schedule**: ready jobs are taken in `(ready time,
///     larger estimated cost, id)` order, one heap per dispatcher, rather than `fusion_lb_profiler`'s
///     pure priority order, because cubes are independent and a priority across them is meaningless.
///  9. **`p2_*` segments.** A phase-2 run is `p2_release`, then alternating `p2_extend`,
///     `p2_closure` and `p2_resume` rows. A resume runs the `G` side until it settles or `(t1)` fires,
///     then each component that needs an `H` solve until it completes or `(t2)` fires; the first
///     trigger ends the segment. `resume_events` counts the `G` side's matcher events.
/// 10. **H-side release is subsumed by reset.** Every component that carries a dummy toward a
///     phase-2 territory is unfinished or within `T` of it, so it is claimed and reset; no `H` region
///     is released in place, and `n_released_h` is always 0.
///
/// ## Usage
///
///   lattice_surgery_profiler --N 2,4,8 --k 4,7 --d 17,25 --p 1e-3 --T 2 --shots S
///                            [--W 2] [--orientation long,short] [--seam-gap 0|1] [--merge-init on|off]
///                            [--warmup 1000] [--seed SEED] [--preload-width W'] [--accel-ns 0]
///                            [--link-ns 0] [--sort] [--verify] [--check-determinism]
///                            [--jobs-log all|sample|none] [--allow-nonconvex]
///                            [--convexity-sources 512] [--out DIR]

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "benchmarks/spec_matching/hires_timer.h"
#include "benchmarks/spec_matching/lattice_surgery_circuit.h"
#include "benchmarks/spec_matching/profiler_util.h"
#include "specmatching/spec_matching/fusion/fusion_solve.h"
#include "specmatching/spec_matching/manifold/ball_decoding.h"

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

#ifndef SPECMATCHING_GIT_HASH
#define SPECMATCHING_GIT_HASH "unknown"
#endif
#ifndef SPECMATCHING_BUILD_FLAGS
#define SPECMATCHING_BUILD_FLAGS "unknown"
#endif
#ifndef SPECMATCHING_STIM_VERSION
#define SPECMATCHING_STIM_VERSION "unknown"
#endif

using namespace pm::spec_matching;
using namespace pm::spec_matching::profiler;
namespace fus = pm::spec_matching::fusion;

namespace {

constexpr size_t SHOT_BATCH = 4096;
constexpr uint32_t NONE = UINT32_MAX;
constexpr int64_t DIST_INF = INT64_MAX / 4;
constexpr pm::weight_int WEIGHT_INF = fus::NO_BOUNDARY_WEIGHT;
constexpr size_t JOBS_SAMPLE_FIRST = 1000;
constexpr size_t INITIAL_EDGE_CAPACITY = 4096;

enum class JobsLog { ALL, SAMPLE, NONE };

// ---------------------------------------------------------------------------------------------
// §8.4. Options.
// ---------------------------------------------------------------------------------------------

struct Options {
    std::vector<size_t> patches = {2};
    std::vector<size_t> core_counts = {4};
    std::vector<size_t> distances = {5};
    std::vector<double> error_rates = {1e-3};
    std::vector<double> horizons = {2.0};
    std::vector<double> slab_widths = {2.0};
    /// 1 = long, 0 = short.
    std::vector<int> orientations = {1};
    size_t seam_gap = 0;
    bool merge_init = true;
    size_t shots = 0;
    size_t warmup = 1000;
    uint64_t seed = 20260910;
    /// In lattice edge weights; negative means the design's default, `max(2W, W + 2T)`, per cell.
    double preload_width = -1;
    double accel_ns = 0;
    double link_ns = 0;
    bool sort_edges = false;
    bool verify = false;
    bool check_determinism = false;
    JobsLog jobs_log = JobsLog::ALL;
    bool allow_nonconvex = false;
    /// Sources per cube for the §9.2 convexity check; 0 means every face-slab node.
    size_t convexity_sources = 512;
    std::string out_dir = "benchmarks/spec_matching/results/lattice_surgery_profiler";
    std::string command_line;
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 0; i < argc; i++) {
        if (i)
            options.command_line += " ";
        options.command_line += argv[i];
    }
    for (int i = 1; i < argc; i++) {
        std::string flag = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument("missing value for " + flag);
            return argv[++i];
        };
        if (flag == "--N") {
            options.patches = parse_list<size_t>(next());
        } else if (flag == "--k") {
            options.core_counts = parse_list<size_t>(next());
        } else if (flag == "--d") {
            options.distances = parse_list<size_t>(next());
        } else if (flag == "--p") {
            options.error_rates = parse_list<double>(next());
        } else if (flag == "--T") {
            options.horizons = parse_list<double>(next());
        } else if (flag == "--W") {
            options.slab_widths = parse_list<double>(next());
        } else if (flag == "--orientation") {
            options.orientations.clear();
            for (const std::string& item : parse_list<std::string>(next())) {
                if (item == "long")
                    options.orientations.push_back(1);
                else if (item == "short")
                    options.orientations.push_back(0);
                else
                    throw std::invalid_argument("--orientation takes long and/or short, got " + item);
            }
        } else if (flag == "--seam-gap") {
            options.seam_gap = std::stoul(next());
        } else if (flag == "--merge-init") {
            std::string v = next();
            if (v == "on")
                options.merge_init = true;
            else if (v == "off")
                options.merge_init = false;
            else
                throw std::invalid_argument("--merge-init takes on or off");
        } else if (flag == "--shots") {
            options.shots = std::stoul(next());
        } else if (flag == "--warmup") {
            options.warmup = std::stoul(next());
        } else if (flag == "--seed") {
            options.seed = std::stoull(next());
        } else if (flag == "--preload-width") {
            options.preload_width = std::stod(next());
        } else if (flag == "--accel-ns") {
            options.accel_ns = std::stod(next());
        } else if (flag == "--link-ns") {
            options.link_ns = std::stod(next());
        } else if (flag == "--sort") {
            options.sort_edges = true;
        } else if (flag == "--verify") {
            options.verify = true;
        } else if (flag == "--check-determinism") {
            options.check_determinism = true;
        } else if (flag == "--jobs-log") {
            std::string v = next();
            if (v == "all")
                options.jobs_log = JobsLog::ALL;
            else if (v == "sample")
                options.jobs_log = JobsLog::SAMPLE;
            else if (v == "none")
                options.jobs_log = JobsLog::NONE;
            else
                throw std::invalid_argument("--jobs-log takes all, sample or none");
        } else if (flag == "--allow-nonconvex") {
            options.allow_nonconvex = true;
        } else if (flag == "--convexity-sources") {
            options.convexity_sources = std::stoul(next());
        } else if (flag == "--out" || flag == "--out-dir") {
            options.out_dir = next();
        } else {
            throw std::invalid_argument("unrecognised flag " + flag);
        }
    }
    if (options.patches.empty() || options.core_counts.empty() || options.distances.empty() ||
        options.error_rates.empty() || options.horizons.empty() || options.slab_widths.empty() ||
        options.orientations.empty())
        throw std::invalid_argument("--N, --k, --d, --p, --T, --W and --orientation each need at least one value");
    if (options.shots == 0)
        throw std::invalid_argument("--shots is required and must be positive");
    for (size_t d : options.distances)
        if (d < 3 || d % 2 == 0)
            throw std::invalid_argument("--d values must be odd and at least 3");
    for (double p : options.error_rates)
        if (!(p > 0) || !(p < 1))
            throw std::invalid_argument("--p values must lie in (0, 1)");
    for (double t : options.horizons)
        if (!(t > 0))
            throw std::invalid_argument("--T values must be positive");
    for (double w : options.slab_widths)
        if (!(w > 0))
            throw std::invalid_argument("--W values must be positive");
    for (size_t n : options.patches)
        if (n == 0)
            throw std::invalid_argument("--N values must be positive");
    for (size_t k : options.core_counts)
        if (k == 0)
            throw std::invalid_argument("--k values must be positive");
    if (options.seam_gap > 1)
        throw std::invalid_argument("--seam-gap must be 0 or 1");
    if (options.accel_ns < 0 || options.link_ns < 0)
        throw std::invalid_argument("--accel-ns and --link-ns must be non-negative");
    return options;
}

std::string fmt_g(double value) {
    char buffer[40];
    std::snprintf(buffer, sizeof(buffer), "%.10g", value);
    return std::string(buffer);
}

std::string stamp() {
    std::time_t now = std::time(nullptr);
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&now));
    return std::string(buffer);
}

struct Digest {
    uint64_t value{1469598103934665603ull};
    template <typename T>
    void add(const T& item) {
        const uint8_t* at = (const uint8_t*)&item;
        for (size_t i = 0; i < sizeof(T); i++) {
            value ^= at[i];
            value *= 1099511628211ull;
        }
    }
};

bool stderr_is_terminal() {
#if defined(__unix__) || defined(__APPLE__)
    if (isatty(fileno(stderr)) != 1)
        return false;
    const char* term = std::getenv("TERM");
    return term != nullptr && std::string(term) != "dumb";
#else
    return false;
#endif
}

void draw_progress(bool interactive, const std::string& label, size_t cell, size_t cells, size_t done, size_t total) {
    double fraction = total ? (double)done / (double)total : 1.0;
    if (!interactive) {
        if (done == 0 || done == total || done % (SHOT_BATCH * 16) == 0) {
            std::fprintf(stderr, "[%zu/%zu] %s  %5.1f%%  %zu/%zu shots\n", cell, cells, label.c_str(), 100 * fraction, done,
                         total);
            std::fflush(stderr);
        }
        return;
    }
    char bar[25];
    size_t filled = (size_t)(fraction * 24);
    for (size_t i = 0; i < 24; i++)
        bar[i] = i < filled ? '#' : '.';
    bar[24] = '\0';
    std::fprintf(stderr, "\r[%zu/%zu] %s  [%s] %5.1f%%  %zu/%zu shots\x1b[K", cell, cells, label.c_str(), bar,
                 100 * fraction, done, total);
    std::fflush(stderr);
}

// ---------------------------------------------------------------------------------------------
// The strip's detector graph `G`, as plain CSR, for every distance computation this binary makes.
// ---------------------------------------------------------------------------------------------

/// `G` with its boundary half-edges split out. The ball tables, the fallback instances and the
/// Dijkstras below all read the same integer weights, so every distance is in one metric.
struct StripGraph {
    size_t n{0};
    size_t num_cubes{0};
    size_t num_observables{0};
    double normalising_constant{0};
    std::vector<uint32_t> off;
    std::vector<uint32_t> nbr;
    /// `rev[e]` is the position of the reverse entry within `nbr[e]`'s own adjacency.
    std::vector<uint32_t> rev;
    std::vector<pm::weight_int> w;
    std::vector<pm::obs_int> obs;
    std::vector<uint8_t> has_bnd;
    std::vector<pm::weight_int> bnd_w;
    std::vector<pm::obs_int> bnd_obs;
    std::vector<uint32_t> cube;
    std::vector<uint32_t> cube_off;
    pm::weight_int w_min{0};
    pm::weight_int w_max{0};
    size_t n_distinct{0};

    inline uint32_t degree(uint32_t u) const {
        return off[u + 1] - off[u];
    }

    void build(const pm::MatchingGraph& g, const LatticeSurgeryCircuit& lc) {
        n = g.nodes.size();
        if (n != lc.num_detectors)
            throw std::runtime_error("the DEM's detector count does not match the circuit's");
        num_cubes = lc.params.patches;
        num_observables = g.num_observables;
        normalising_constant = g.normalising_constant;
        off.assign(n + 1, 0);
        has_bnd.assign(n, 0);
        bnd_w.assign(n, 0);
        bnd_obs.assign(n, 0);
        std::vector<pm::weight_int> all;
        for (size_t u = 0; u < n; u++) {
            const pm::DetectorNode& node = g.nodes[u];
            for (size_t k = 0; k < node.neighbors.size(); k++) {
                all.push_back(node.neighbor_weights[k]);
                if (node.neighbors[k] == nullptr) {
                    has_bnd[u] = 1;
                    bnd_w[u] = node.neighbor_weights[k];
                    bnd_obs[u] = node.neighbor_observables[k];
                } else {
                    off[u + 1]++;
                }
            }
        }
        for (size_t u = 0; u < n; u++)
            off[u + 1] += off[u];
        nbr.resize(off[n]);
        w.resize(off[n]);
        obs.resize(off[n]);
        rev.resize(off[n]);
        for (size_t u = 0; u < n; u++) {
            const pm::DetectorNode& node = g.nodes[u];
            uint32_t at = off[u];
            for (size_t k = 0; k < node.neighbors.size(); k++) {
                if (node.neighbors[k] == nullptr)
                    continue;
                nbr[at] = (uint32_t)(node.neighbors[k] - g.nodes.data());
                w[at] = node.neighbor_weights[k];
                obs[at] = node.neighbor_observables[k];
                at++;
            }
        }
        for (uint32_t u = 0; u < n; u++) {
            for (uint32_t e = off[u]; e < off[u + 1]; e++) {
                uint32_t v = nbr[e];
                uint32_t found = NONE;
                for (uint32_t f = off[v]; f < off[v + 1]; f++) {
                    if (nbr[f] == u) {
                        found = f - off[v];
                        break;
                    }
                }
                if (found == NONE)
                    throw std::runtime_error("G adjacency is not symmetric");
                rev[e] = found;
            }
        }
        cube = lc.detector_cube;
        cube_off = lc.cube_offsets;
        std::sort(all.begin(), all.end());
        w_min = all.empty() ? 0 : all.front();
        w_max = all.empty() ? 0 : all.back();
        n_distinct = (size_t)(std::unique(all.begin(), all.end()) - all.begin());
    }
};

/// A bounded Dijkstra over `G`, reset through its own touched list so that a call costs what it
/// visits and not `|G|`.
struct Dijkstra {
    std::vector<int64_t> dist;
    std::vector<uint32_t> touched;
    std::vector<std::pair<int64_t, uint32_t>> heap;

    void configure(size_t n) {
        dist.assign(n, DIST_INF);
        touched.reserve(1024);
        heap.reserve(1024);
    }
    void clear() {
        for (uint32_t t : touched)
            dist[t] = DIST_INF;
        touched.clear();
        heap.clear();
    }

    /// Settles every node within `limit` (inclusive) of the sources, through nodes `pass` accepts
    /// (sources always count), calling `visit(v, d)` once per settled node.
    template <typename Pass, typename Visit>
    void run(const StripGraph& g, const uint32_t* sources, size_t count, int64_t limit, const Pass& pass,
             const Visit& visit) {
        clear();
        auto cmp = std::greater<std::pair<int64_t, uint32_t>>();
        for (size_t i = 0; i < count; i++) {
            uint32_t s = sources[i];
            if (dist[s] == 0)
                continue;
            if (dist[s] == DIST_INF)
                touched.push_back(s);
            dist[s] = 0;
            heap.push_back({0, s});
        }
        std::make_heap(heap.begin(), heap.end(), cmp);
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end(), cmp);
            auto [d, u] = heap.back();
            heap.pop_back();
            if (d > dist[u])
                continue;
            visit(u, d);
            for (uint32_t e = g.off[u]; e < g.off[u + 1]; e++) {
                uint32_t v = g.nbr[e];
                int64_t nd = d + (int64_t)g.w[e];
                if (nd > limit || nd >= dist[v] || !pass(v))
                    continue;
                if (dist[v] == DIST_INF)
                    touched.push_back(v);
                dist[v] = nd;
                heap.push_back({nd, v});
                std::push_heap(heap.begin(), heap.end(), cmp);
            }
        }
    }
};

// ---------------------------------------------------------------------------------------------
// §2.2, §5.1, §7.3. The static regions of each cube, per `(T, W, preload)`.
// ---------------------------------------------------------------------------------------------

struct CubeRegion {
    /// `S_c = G_c ∪ B_c`, the initial territory, and its outer layer.
    std::vector<uint32_t> s_nodes;
    std::vector<uint32_t> outer_nodes;
    std::vector<uint8_t> in_s;
    /// `P_c`, the preloaded region: everything within the preload width of `G_c`.
    std::vector<uint8_t> in_p;
    /// Outside `S_c` and strictly within `T` of it: where an obstacle's ball can reach `S_c`.
    std::vector<uint8_t> in_band;
    /// The face slab of §9.2: cube nodes within `W` of a face, plus `B_c`.
    std::vector<uint32_t> face_slab;
    std::vector<uint8_t> in_face_slab;
    /// The cubes whose nodes `P_c` reaches, which is where phase 1 looks for neighbour defects.
    uint32_t reach_lo{0};
    uint32_t reach_hi{0};
    size_t b_nodes{0};
    size_t p_nodes{0};
    size_t band_nodes{0};
    std::vector<uint8_t> s_cubes;
    /// §9.2.
    uint64_t convexity_pairs{0};
    uint64_t convexity_mismatches{0};
    size_t convexity_sources_used{0};
};

void compute_cube_region(const StripGraph& g, uint32_t c, int64_t W_int, int64_t T_int, int64_t P_int,
                         Dijkstra& dij, CubeRegion& r) {
    size_t n = g.n;
    r.in_s.assign(n, 0);
    r.in_p.assign(n, 0);
    r.in_band.assign(n, 0);
    r.in_face_slab.assign(n, 0);
    r.s_nodes.clear();
    r.outer_nodes.clear();
    r.face_slab.clear();
    r.s_cubes.assign(g.num_cubes, 0);
    std::vector<uint32_t> cube_nodes;
    for (uint32_t v = g.cube_off[c]; v < g.cube_off[c + 1]; v++)
        cube_nodes.push_back(v);
    r.b_nodes = 0;
    r.p_nodes = 0;
    r.reach_lo = c;
    r.reach_hi = c;
    dij.run(g, cube_nodes.data(), cube_nodes.size(), P_int, [](uint32_t) { return true; },
            [&](uint32_t v, int64_t d) {
                r.in_p[v] = 1;
                r.p_nodes++;
                r.reach_lo = std::min(r.reach_lo, g.cube[v]);
                r.reach_hi = std::max(r.reach_hi, g.cube[v]);
                if (g.cube[v] == c || d <= W_int) {
                    r.in_s[v] = 1;
                    r.s_nodes.push_back(v);
                    r.s_cubes[g.cube[v]] = 1;
                    if (g.cube[v] != c)
                        r.b_nodes++;
                }
            });
    std::sort(r.s_nodes.begin(), r.s_nodes.end());
    for (uint32_t v : r.s_nodes) {
        for (uint32_t e = g.off[v]; e < g.off[v + 1]; e++) {
            if (!r.in_s[g.nbr[e]]) {
                r.outer_nodes.push_back(v);
                break;
            }
        }
    }
    r.band_nodes = 0;
    dij.run(g, r.s_nodes.data(), r.s_nodes.size(), T_int - 1, [](uint32_t) { return true; },
            [&](uint32_t v, int64_t) {
                if (!r.in_s[v]) {
                    r.in_band[v] = 1;
                    r.band_nodes++;
                }
            });
    // Face slab: cube nodes within W of a cube node with an edge out of the cube, plus B_c.
    std::vector<uint32_t> faces;
    for (uint32_t v : cube_nodes) {
        for (uint32_t e = g.off[v]; e < g.off[v + 1]; e++) {
            if (g.cube[g.nbr[e]] != c) {
                faces.push_back(v);
                break;
            }
        }
    }
    dij.run(g, faces.data(), faces.size(), W_int, [&](uint32_t v) { return g.cube[v] == c; },
            [&](uint32_t v, int64_t) {
                if (!r.in_face_slab[v]) {
                    r.in_face_slab[v] = 1;
                    r.face_slab.push_back(v);
                }
            });
    for (uint32_t v : r.s_nodes) {
        if (g.cube[v] != c && !r.in_face_slab[v]) {
            r.in_face_slab[v] = 1;
            r.face_slab.push_back(v);
        }
    }
    std::sort(r.face_slab.begin(), r.face_slab.end());
}

/// §9.2's geodesic convexity, once per cube: for face-slab pairs within `2T + W`, the distance over
/// `P_c` (which contains `G_c`) against the distance over `G`.
void check_convexity(const StripGraph& g, int64_t limit, size_t max_sources, Dijkstra& full, Dijkstra& restricted,
                     CubeRegion& r) {
    r.convexity_pairs = 0;
    r.convexity_mismatches = 0;
    r.convexity_sources_used = 0;
    size_t total = r.face_slab.size();
    if (total == 0)
        return;
    size_t step = (max_sources == 0 || total <= max_sources) ? 1 : (total + max_sources - 1) / max_sources;
    for (size_t at = 0; at < total; at += step) {
        uint32_t x = r.face_slab[at];
        r.convexity_sources_used++;
        restricted.run(g, &x, 1, limit, [&](uint32_t v) { return r.in_p[v] != 0; }, [](uint32_t, int64_t) {});
        full.run(g, &x, 1, limit, [](uint32_t) { return true; }, [&](uint32_t v, int64_t d) {
            if (!r.in_face_slab[v] || v == x)
                return;
            r.convexity_pairs++;
            if (restricted.dist[v] != d)
                r.convexity_mismatches++;
        });
    }
}

// ---------------------------------------------------------------------------------------------
// §7.1, §7.3, §7.4. A cube's fallback instance: `G` with a boundary slot on every node, masked.
// ---------------------------------------------------------------------------------------------

enum SlotKind : uint8_t { SLOT_NONE = 0, SLOT_REAL = 1, SLOT_OBST = 2, SLOT_OUTER = 3 };

/// The whole strip's `G`, built **once** per corpus, with a boundary half-edge at index 0 of every
/// node (§7.4: "one boundary-edge slot ... infinite weight by default", made a masked slot here) and
/// an edge mask that hides everything outside `Ω`. The baseline is `Ω = S_c`, the outer layer a
/// weight-0 virtual boundary, and the real code boundary inside `S_c`; a shot only ever edits the
/// nodes it touches, and `undo` puts exactly those back.
struct FallbackInstance {
    pm::Mwpm mwpm;
    const StripGraph* g{nullptr};
    uint32_t cube{0};
    std::vector<uint32_t> moff;
    std::vector<uint8_t> mask;
    std::vector<uint8_t> base_mask;
    std::vector<uint8_t> in_omega;
    std::vector<uint8_t> base_in_omega;
    std::vector<uint8_t> outer;
    std::vector<uint8_t> base_outer;
    std::vector<uint8_t> slot_kind;
    std::vector<uint8_t> base_slot_kind;
    std::vector<pm::weight_int> base_slot_w;
    std::vector<pm::obs_int> base_slot_obs;
    std::vector<uint8_t> touched_flag;
    std::vector<uint32_t> touched;
    /// Nodes whose slot holds an obstacle or outer weight this shot, for the phase-2 release.
    std::vector<uint32_t> virtual_slots;
    /// Old-`Ω` nodes adjacent to nodes just added, whose events have to be rescheduled.
    std::vector<uint32_t> resched;
    std::vector<uint32_t> outer_dirty;
    size_t omega_size{0};
    size_t base_omega_size{0};

    void build(const StripGraph& graph, uint32_t c) {
        g = &graph;
        cube = c;
        size_t n = graph.n;
        pm::MatchingGraph mg(n, graph.num_observables, graph.normalising_constant);
        moff.assign(n + 1, 0);
        for (size_t u = 0; u < n; u++)
            moff[u + 1] = moff[u] + 1 + graph.degree((uint32_t)u);
        for (size_t u = 0; u < n; u++) {
            pm::DetectorNode& node = mg.nodes[u];
            node.neighbors.push_back(nullptr);
            node.neighbor_weights.push_back(graph.has_bnd[u] ? graph.bnd_w[u] : 0);
            node.neighbor_observables.push_back(graph.has_bnd[u] ? graph.bnd_obs[u] : 0);
            for (uint32_t e = graph.off[u]; e < graph.off[u + 1]; e++) {
                node.neighbors.push_back(&mg.nodes[graph.nbr[e]]);
                node.neighbor_weights.push_back(graph.w[e]);
                node.neighbor_observables.push_back(graph.obs[e]);
            }
            node.neighbor_implied_weights.resize(node.neighbors.size());
        }
        mwpm = pm::Mwpm(pm::GraphFlooder(std::move(mg)));
        mask.assign(moff[n], 1);
        in_omega.assign(n, 0);
        outer.assign(n, 0);
        slot_kind.assign(n, SLOT_NONE);
        touched_flag.assign(n, 0);
        touched.reserve(4096);
        mwpm.flooder.edge_mask = mask.data();
        mwpm.flooder.edge_mask_offsets = moff.data();
        mwpm.flooder.horizon = pm::NO_HORIZON;
        mwpm.flooder.dual_cap = pm::NO_HORIZON;
        mwpm.flooder.outer_layer = nullptr;
    }

    inline pm::DetectorNode& node(uint32_t u) {
        return mwpm.flooder.graph.nodes[u];
    }

    /// Installs the baseline for territory `r`. Setup-time work, `O(|G|)`, once per `(T, W)`.
    void set_baseline(const CubeRegion& r) {
        size_t n = g->n;
        std::fill(mask.begin(), mask.end(), 1);
        std::fill(in_omega.begin(), in_omega.end(), 0);
        std::fill(outer.begin(), outer.end(), 0);
        for (uint32_t v : r.s_nodes)
            in_omega[v] = 1;
        for (uint32_t v : r.outer_nodes)
            outer[v] = 1;
        for (size_t u = 0; u < n; u++) {
            pm::DetectorNode& nd = node((uint32_t)u);
            nd.neighbor_weights[0] = g->has_bnd[u] ? g->bnd_w[u] : 0;
            nd.neighbor_observables[0] = g->has_bnd[u] ? g->bnd_obs[u] : 0;
            slot_kind[u] = SLOT_NONE;
            if (!in_omega[u])
                continue;
            for (uint32_t e = g->off[u]; e < g->off[u + 1]; e++)
                mask[moff[u] + 1 + (e - g->off[u])] = in_omega[g->nbr[e]] ? 0 : 1;
            if (outer[u]) {
                // The weight-0 virtual boundary of deviation 2.
                nd.neighbor_weights[0] = 0;
                nd.neighbor_observables[0] = 0;
                slot_kind[u] = SLOT_OUTER;
                mask[moff[u]] = 0;
            } else if (g->has_bnd[u]) {
                slot_kind[u] = SLOT_REAL;
                mask[moff[u]] = 0;
            }
        }
        base_mask = mask;
        base_in_omega = in_omega;
        base_outer = outer;
        base_slot_kind = slot_kind;
        base_slot_w.assign(n, 0);
        base_slot_obs.assign(n, 0);
        for (size_t u = 0; u < n; u++) {
            base_slot_w[u] = node((uint32_t)u).neighbor_weights[0];
            base_slot_obs[u] = node((uint32_t)u).neighbor_observables[0];
        }
        base_omega_size = r.s_nodes.size();
        omega_size = base_omega_size;
        touched.clear();
        std::fill(touched_flag.begin(), touched_flag.end(), 0);
        virtual_slots.clear();
    }

    inline void touch(uint32_t u) {
        if (!touched_flag[u]) {
            touched_flag[u] = 1;
            touched.push_back(u);
        }
    }

    inline void set_slot(uint32_t u, pm::weight_int weight, pm::obs_int obs_mask, uint8_t kind) {
        touch(u);
        pm::DetectorNode& nd = node(u);
        nd.neighbor_weights[0] = weight;
        nd.neighbor_observables[0] = obs_mask;
        mask[moff[u]] = 0;
        slot_kind[u] = kind;
    }

    /// The slot as the final problem has it: the real code boundary, or none.
    inline void restore_real_slot(uint32_t u) {
        touch(u);
        pm::DetectorNode& nd = node(u);
        if (g->has_bnd[u]) {
            nd.neighbor_weights[0] = g->bnd_w[u];
            nd.neighbor_observables[0] = g->bnd_obs[u];
            mask[moff[u]] = 0;
            slot_kind[u] = SLOT_REAL;
        } else {
            mask[moff[u]] = 1;
            slot_kind[u] = SLOT_NONE;
        }
    }

    /// Adds `x` to `Ω`: unmasks every edge between `x` and `Ω`, gives `x` its real boundary slot,
    /// and records the outer-layer and reschedule bookkeeping. `new_outer` receives `x` if it is
    /// on the outer layer once added.
    void add_node(uint32_t x, std::vector<uint32_t>& new_outer) {
        if (in_omega[x])
            return;
        touch(x);
        in_omega[x] = 1;
        omega_size++;
        restore_real_slot(x);
        bool is_outer = false;
        for (uint32_t e = g->off[x]; e < g->off[x + 1]; e++) {
            uint32_t v = g->nbr[e];
            if (in_omega[v]) {
                mask[moff[x] + 1 + (e - g->off[x])] = 0;
                touch(v);
                mask[moff[v] + 1 + g->rev[e]] = 0;
                resched.push_back(v);
                outer_dirty.push_back(v);
            } else {
                is_outer = true;
            }
        }
        outer[x] = is_outer ? 1 : 0;
        if (is_outer)
            new_outer.push_back(x);
    }

    /// Recomputes the outer flag of every old-`Ω` node next to something just added.
    void refresh_outer() {
        for (uint32_t v : outer_dirty) {
            bool is_outer = false;
            for (uint32_t e = g->off[v]; e < g->off[v + 1]; e++) {
                if (!in_omega[g->nbr[e]]) {
                    is_outer = true;
                    break;
                }
            }
            touch(v);
            outer[v] = is_outer ? 1 : 0;
        }
        outer_dirty.clear();
    }

    /// Reschedules the events of old-`Ω` nodes whose adjacency just grew. A node of a shrinking
    /// region is left alone: shrinking regions carry no node events, and the region reschedules all
    /// of its nodes when it stops shrinking.
    void apply_reschedules() {
        for (uint32_t v : resched) {
            pm::DetectorNode& nd = node(v);
            if (nd.region_that_arrived_top != nullptr && nd.region_that_arrived_top->radius.is_shrinking())
                continue;
            mwpm.flooder.reschedule_events_at_detector_node(nd);
        }
        resched.clear();
    }

    /// Puts every touched node back to the baseline. Regions must already have been extracted.
    void undo() {
        for (uint32_t u : touched) {
            std::copy(base_mask.begin() + moff[u], base_mask.begin() + moff[u + 1], mask.begin() + moff[u]);
            pm::DetectorNode& nd = node(u);
            nd.neighbor_weights[0] = base_slot_w[u];
            nd.neighbor_observables[0] = base_slot_obs[u];
            slot_kind[u] = base_slot_kind[u];
            in_omega[u] = base_in_omega[u];
            outer[u] = base_outer[u];
            touched_flag[u] = 0;
        }
        touched.clear();
        virtual_slots.clear();
        resched.clear();
        outer_dirty.clear();
        omega_size = base_omega_size;
        mwpm.flooder.outer_layer = nullptr;
        mwpm.flooder.outer_layer_hit = nullptr;
    }
};

/// A detection event injected at the current clock, with an even `y`-intercept (deviation 6).
void inject_uncapped(pm::Mwpm& mwpm, uint32_t node_index) {
    pm::GraphFlooder& flooder = mwpm.flooder;
    pm::GraphFillRegion* region = flooder.region_arena.alloc_default_constructed();
    pm::cumulative_time_int t = flooder.queue.cur_time;
    region->radius = (t & 1) ? pm::VaryingCT::growing_value_at_time(1, t)
                             : pm::VaryingCT::growing_varying_with_zero_distance_at_time(t);
    region->alt_tree_node = mwpm.node_arena.alloc_constructed(region);
    flooder.do_region_created_at_empty_detector_node(*region, flooder.graph.nodes[node_index]);
}

/// Runs a fallback instance until it settles or `(t1)` fires. Returns the claimed outer-layer node,
/// or `NONE` once the queue has drained.
uint32_t run_g_until_trigger(pm::Mwpm& mwpm, uint64_t& events) {
    pm::GraphFlooder& flooder = mwpm.flooder;
    while (true) {
        pm::MwpmEvent event = flooder.run_until_next_mwpm_notification();
        if (flooder.outer_layer_hit != nullptr) {
            uint32_t hit = (uint32_t)(flooder.outer_layer_hit - flooder.graph.nodes.data());
            flooder.outer_layer_hit = nullptr;
            return hit;
        }
        if (event.event_type == pm::NO_EVENT)
            return NONE;
        mwpm.process_event(event);
        events++;
    }
}

/// Is `region` (top level) matched to the boundary?
inline bool matched_to_boundary(const pm::GraphFillRegion* region) {
    return region != nullptr && region->alt_tree_node == nullptr && region->match.region == nullptr &&
           region->match.edge.loc_from != nullptr;
}

/// Tears down every region and alternating-tree node owned by `defects` (H indices of one
/// component, or of one subtree of it) in an instance whose edges to the rest are masked, then
/// drains the queue of their events. Every object is returned to its arena, so the instance's
/// liveness lists — which `fusion::run_until_settled` reads — stay exact.
struct Teardown {
    std::vector<pm::GraphFillRegion*> regions;
    std::vector<pm::AltTreeNode*> alts;
    uint64_t stray_events{0};

    void run(pm::Mwpm& mwpm, const uint32_t* defects, size_t count) {
        pm::MatchingGraph& graph = mwpm.flooder.graph;
        regions.clear();
        alts.clear();
        for (size_t i = 0; i < count; i++) {
            pm::GraphFillRegion* r = graph.nodes[defects[i]].region_that_arrived;
            while (r != nullptr) {
                regions.push_back(r);
                r = r->blossom_parent;
            }
        }
        std::sort(regions.begin(), regions.end());
        regions.erase(std::unique(regions.begin(), regions.end()), regions.end());
        for (pm::GraphFillRegion* r : regions)
            if (r->alt_tree_node != nullptr)
                alts.push_back(r->alt_tree_node);
        std::sort(alts.begin(), alts.end());
        alts.erase(std::unique(alts.begin(), alts.end()), alts.end());
        // Invalidate first, drain second, free last: a queued event names its object, and must find
        // a tracker that refuses it rather than a slot the arena has already handed out again.
        for (pm::GraphFillRegion* r : regions) {
            r->shrink_event_tracker.set_no_desired_event();
            for (pm::DetectorNode* nd : r->shell_area)
                nd->reset();
        }
        for (size_t i = 0; i < count; i++)
            graph.nodes[defects[i]].reset();
        while (true) {
            pm::FloodCheckEvent ev = mwpm.flooder.dequeue_valid();
            if (ev.tentative_event_type == pm::NO_FLOOD_CHECK_EVENT)
                break;
            stray_events++;
        }
        for (pm::AltTreeNode* a : alts)
            mwpm.node_arena.del(a);
        for (pm::GraphFillRegion* r : regions)
            mwpm.flooder.region_arena.del(r);
    }
};

// ---------------------------------------------------------------------------------------------
// B3. A cube decoded on its own, with its seam edges turned into boundary edges.
// ---------------------------------------------------------------------------------------------

struct CubePm {
    pm::Mwpm mwpm;
    uint32_t lo{0};
    uint32_t hi{0};
    std::vector<uint64_t> dets;
    std::vector<uint8_t> obs;

    void build(const StripGraph& g, uint32_t c) {
        lo = g.cube_off[c];
        hi = g.cube_off[c + 1];
        size_t n = hi - lo;
        pm::MatchingGraph mg(n, g.num_observables, g.normalising_constant);
        std::vector<pm::weight_int> bw(n, WEIGHT_INF);
        std::vector<pm::obs_int> bo(n, 0);
        auto obs_list = [&](pm::obs_int m) {
            std::vector<size_t> out;
            for (size_t i = 0; i < 64; i++)
                if (m >> i & 1)
                    out.push_back(i);
            return out;
        };
        for (uint32_t u = lo; u < hi; u++) {
            if (g.has_bnd[u]) {
                bw[u - lo] = g.bnd_w[u];
                bo[u - lo] = g.bnd_obs[u];
            }
            for (uint32_t e = g.off[u]; e < g.off[u + 1]; e++) {
                uint32_t v = g.nbr[e];
                if (v >= lo && v < hi) {
                    if (v > u)
                        mg.add_edge(u - lo, v - lo, (pm::signed_weight_int)g.w[e], obs_list(g.obs[e]));
                } else if (g.w[e] < bw[u - lo]) {
                    bw[u - lo] = g.w[e];
                    bo[u - lo] = g.obs[e];
                }
            }
        }
        for (size_t i = 0; i < n; i++)
            if (bw[i] != WEIGHT_INF)
                mg.add_boundary_edge(i, (pm::signed_weight_int)bw[i], obs_list(bo[i]));
        mwpm = pm::Mwpm(pm::GraphFlooder(std::move(mg)));
        obs.assign(std::max<size_t>(1, g.num_observables), 0);
    }
};

// ---------------------------------------------------------------------------------------------
// §8.3. Jobs and the simulated schedule.
// ---------------------------------------------------------------------------------------------

enum class Kind : uint8_t {
    ACCEL,
    MGR_UF,
    MGR_TREE,
    MGR_SCATTER,
    LEAF,
    FUSE,
    G_SEAM_CHECK,
    SEAM_FUSE,
    EXTRACT,
    FB_PHASE1,
    FB_OBSTACLES,
    P2_RELEASE,
    P2_CLOSURE,
    P2_RESUME,
    P2_EXTEND,
    G_COMBINE,
    DISPATCH,
};

const char* kind_name(Kind k) {
    switch (k) {
        case Kind::ACCEL:
            return "accel";
        case Kind::MGR_UF:
            return "mgr_uf";
        case Kind::MGR_TREE:
            return "mgr_tree";
        case Kind::MGR_SCATTER:
            return "mgr_scatter";
        case Kind::LEAF:
            return "leaf";
        case Kind::FUSE:
            return "fuse";
        case Kind::G_SEAM_CHECK:
            return "g_seam_check";
        case Kind::SEAM_FUSE:
            return "seam_fuse";
        case Kind::EXTRACT:
            return "extract";
        case Kind::FB_PHASE1:
            return "fb_phase1";
        case Kind::FB_OBSTACLES:
            return "fb_obstacles";
        case Kind::P2_RELEASE:
            return "p2_release";
        case Kind::P2_CLOSURE:
            return "p2_closure";
        case Kind::P2_RESUME:
            return "p2_resume";
        case Kind::P2_EXTEND:
            return "p2_extend";
        case Kind::G_COMBINE:
            return "g_combine";
        case Kind::DISPATCH:
            return "dispatch";
    }
    return "?";
}

enum class CoreClass : uint8_t { ACC, MGR, SOLVER, FB, GM };

enum JobStatus : uint8_t { ST_COMPLETE = 0, ST_TRUNCATED, ST_EVICTED, ST_T1, ST_T2 };

const char* status_name(uint8_t s) {
    switch (s) {
        case ST_COMPLETE:
            return "COMPLETE";
        case ST_TRUNCATED:
            return "TRUNCATED";
        case ST_EVICTED:
            return "EVICTED";
        case ST_T1:
            return "TRIGGER_T1";
        case ST_T2:
            return "TRIGGER_T2";
    }
    return "?";
}

enum Determined : uint8_t { BY_NOTHING = 0, BY_DEP, BY_CORE, BY_DISPATCH };

struct Job {
    Kind kind{Kind::LEAF};
    CoreClass cls{CoreClass::SOLVER};
    int32_t cube{-1};
    int32_t seam{-1};
    int32_t comp{-1};
    uint32_t dep_begin{0};
    uint32_t dep_end{0};
    /// For a fusion, a seam fusion or an extraction: the jobs whose cores it may run on (the one
    /// that finished last wins).
    uint32_t bind_begin{0};
    uint32_t bind_end{0};
    /// For a seam fusion: `(cube, defects)` of everything it fuses, so the link charge can be taken
    /// against whichever side does not host it.
    uint32_t moved_begin{0};
    uint32_t moved_end{0};
    uint32_t moved_total{0};
    /// For phase 2: items transferred between cores, charged `--link-ns` each.
    uint64_t link_items{0};
    uint64_t measured{0};
    uint64_t priority{0};
    uint32_t n_items{0};
    uint8_t status{ST_COMPLETE};
    /// Cube index of the dispatching manager, `num_cubes` for the global manager, `-1` for none.
    int32_t dispatcher{-1};
    /// Logged but not scheduled: `fb_obstacles` is a sub-interval of `fb_phase1`.
    bool informational{false};
    /// Filled by the scheduler.
    uint64_t ready{0};
    uint32_t ready_from{NONE};
    uint64_t start{0};
    uint64_t end{0};
    uint64_t delay{0};
    int32_t core{-1};
    uint32_t pred{NONE};
    uint8_t determined{BY_NOTHING};
    uint32_t pending{0};
    uint32_t moved_defects{0};
};

// ---------------------------------------------------------------------------------------------
// The per-shot record and the per-cell aggregate.
// ---------------------------------------------------------------------------------------------

struct ShotRow {
    uint64_t shot{0};
    uint32_t n_def{0};
    uint32_t n_def_cube_min{0};
    uint32_t n_def_cube_max{0};
    uint32_t n_edges{0};
    uint32_t n_cross_edges{0};
    uint32_t n_seam_components{0};
    uint32_t n_spanning{0};
    uint32_t max_span{0};
    uint32_t seam_depth{0};
    uint32_t n_cubes_truncated{0};
    uint32_t n_trunc_components{0};
    /// 0 none, 1 leaf, 2 fuse, 3 seam.
    uint8_t trunc_where{0};
    uint32_t n_groups{0};
    uint32_t omega_nodes{0};
    uint32_t omega_cubes{0};
    uint32_t n_t1{0};
    uint32_t n_t2{0};
    uint32_t n_closure_steps{0};
    uint32_t n_reset_components{0};
    uint32_t n_reset_defects{0};
    uint32_t n_obstacles{0};
    uint32_t n_promotions{0};
    uint32_t n_handoff{0};
    uint32_t handoff_nodes{0};
    uint32_t n_released_g{0};
    uint32_t n_released_h{0};
    uint32_t n_moved_defects{0};
    uint8_t p2_conflict{0};
    uint32_t n_pieces{0};
    uint32_t n_fusions{0};
    uint32_t n_seam_fusions{0};
    uint64_t t_spec_ns{0};
    uint64_t t_fb_max_ns{0};
    uint64_t t_p2_ns{0};
    uint64_t t_sys_ns{0};
    uint64_t t_strip_pm_ns{0};
    uint64_t t_global_escalate_ns{0};
    uint64_t t_cube_pm_max_ns{0};
    uint64_t cp_accel{0}, cp_mgr{0}, cp_leaf{0}, cp_fuse{0}, cp_seam{0}, cp_extract{0}, cp_gm{0}, cp_fb{0}, cp_p2{0},
        cp_link{0}, cp_idle{0};
    pm::total_weight_int weight{0};
    pm::total_weight_int ref_weight{0};
    pm::obs_int obs_xor{0};
    pm::obs_int ref_obs{0};
};

const char* trunc_where_name(uint8_t w) {
    switch (w) {
        case 1:
            return "leaf";
        case 2:
            return "fuse";
        case 3:
            return "seam";
    }
    return "";
}

struct PhaseTwoRow {
    std::string cubes;
    uint64_t start{0};
    uint64_t end{0};
    uint32_t omega_initial{0};
    uint32_t omega_final{0};
    uint32_t omega_cubes{0};
    uint32_t n_t1{0};
    uint32_t n_t2{0};
    uint32_t n_closure_steps{0};
    uint32_t n_reset_components{0};
    uint32_t n_reset_defects{0};
    uint32_t n_handoff{0};
    double max_y_over_t{0};
    uint64_t resume_events{0};
    std::vector<uint32_t> job_ids;
};

/// One `t_*` column of a cell: every shot's value, and whether that shot truncated.
struct Series {
    std::vector<long long> values;
    std::vector<uint8_t> truncated;
    void clear() {
        values.clear();
        truncated.clear();
    }
    void add(uint64_t v, bool t) {
        values.push_back((long long)v);
        truncated.push_back(t ? 1 : 0);
    }
};

struct CellAggregate {
    size_t N{0}, k{0}, d{0};
    double p{0}, T{0}, W{0};
    int orientation{1};
    uint64_t shots{0};
    uint64_t shots_truncated{0};
    uint64_t cubes_truncated{0};
    uint64_t trunc_components{0};
    uint64_t trunc_leaf{0}, trunc_fuse{0}, trunc_seam{0};
    Series t_spec, t_fb_max, t_p2, t_sys, t_strip_pm, t_global_escalate, t_cube_pm_max;
    __int128 cp_sum[11] = {0};
    std::vector<uint64_t> seam_depth_hist, omega_cubes_hist, triggers_hist, reset_defects_hist;
    uint64_t handoff_events{0};
    uint64_t shots_with_handoff{0};
    uint64_t omega_cubes_gt3{0};
    uint64_t p2_conflicts{0};
    uint64_t verify_shots{0};
    uint64_t weight_mismatches{0};
    uint64_t obs_mismatches{0};
    int64_t first_weight_mismatch_shot{-1};
    uint64_t determinism_mismatches{0};
    int determinism_check{-1};
    uint64_t n_spanning_shots{0};
    uint64_t seam_component_sum{0};
    uint64_t spanning_component_sum{0};
    uint64_t obstacles_sum{0};
    uint64_t promotions_sum{0};
    uint64_t moved_sum{0};
    uint64_t released_g_sum{0};
    uint64_t t1_sum{0}, t2_sum{0}, closure_sum{0};
    uint64_t stray_events{0};
    /// Debug-twin invariants (§9.1).
    uint64_t dbg_invariant_violations{0};
    uint64_t dbg_cap_violations{0};
    uint64_t dbg_injection_violations{0};
    uint64_t dbg_feasibility_violations{0};
    uint64_t dbg_g_unsettled{0};
    /// §9.2.
    std::vector<uint64_t> convexity_pairs, convexity_mismatches;
    std::string label;

    void reset(size_t hist_size) {
        seam_depth_hist.assign(hist_size, 0);
        omega_cubes_hist.assign(hist_size, 0);
        triggers_hist.assign(4096, 0);
        reset_defects_hist.assign(4096, 0);
    }
};

// ---------------------------------------------------------------------------------------------
// Everything built once per `(N, d, p, orientation)`.
// ---------------------------------------------------------------------------------------------

struct Corpus {
    LatticeSurgeryCircuit lc;
    stim::DetectorErrorModel dem;
    std::unique_ptr<BallDecoder> decoder;
    StripGraph g;
    double unit{0};
    double w_p{0};
    /// B1's instance.
    pm::Mwpm stock;
    std::vector<std::unique_ptr<CubePm>> cube_pm;
    /// The speculative path's shared `H` instance, sized for the whole strip.
    std::unique_ptr<BallMwpm> shared;
    fus::FusionInstance fusion;
    std::vector<std::unique_ptr<FallbackInstance>> fallback;
    std::vector<CubeRegion> regions;
};

// ---------------------------------------------------------------------------------------------
// One shot.
// ---------------------------------------------------------------------------------------------

struct QuotientAdj {
    uint32_t other;
    uint32_t count;
};

struct MergeCandidate {
    uint32_t count;
    uint32_t lo;
    uint32_t hi;
    bool operator<(const MergeCandidate& other) const {
        if (count != other.count)
            return count > other.count;
        if (lo != other.lo)
            return lo > other.lo;
        return hi > other.hi;
    }
};

struct Group {
    std::vector<uint32_t> cubes;
    uint32_t c0{0};
    std::vector<uint32_t> unfinished_comps;
};

struct ShotRunner {
    const Options* options{nullptr};
    Corpus* corpus{nullptr};
    size_t N{0};
    size_t k{0};
    horizon_int T_int{0};
    int64_t W_int{0};
    int64_t P_int{0};
    bool profile{false};
    bool sort_by_weight{false};
    double accel_ns{0};
    double link_ns{0};

    inline uint64_t tick() const {
        return profile ? hires_now_ns() : 0;
    }

    // ---- the shot's input
    const BallGraph* h{nullptr};
    uint32_t n{0};
    std::vector<uint32_t> hoff;  // N + 1
    std::vector<uint32_t> det_to_h;
    std::vector<uint32_t> cube_of_h;
    std::vector<std::vector<uint32_t>> cube_edges;
    std::vector<std::vector<uint32_t>> seam_edges;
    std::vector<uint32_t> bnd_lo, bnd_hi;
    std::vector<pm::obs_int> edge_obs;
    std::vector<uint32_t> identity;

    // ---- components of the strip's H (bookkeeping, untimed)
    std::vector<uint32_t> comp_parent;
    std::vector<uint32_t> comp_of;
    std::vector<uint32_t> comp_begin;
    std::vector<uint32_t> comp_members;
    std::vector<uint32_t> comp_cube_lo, comp_cube_hi;
    std::vector<int32_t> comp_claimed;
    std::vector<uint8_t> comp_unfinished;
    std::vector<uint8_t> comp_truncated;
    std::vector<uint8_t> comp_needs_h;
    std::vector<std::vector<uint32_t>> comp_jobs;
    uint32_t n_comps{0};

    // ---- §4.1 the per-cube cut (one StripCut, cubes using disjoint ranges)
    std::vector<uint32_t> parent, set_size;
    std::vector<pm::weight_int> w_boundary, bnd;
    std::vector<pm::obs_int> obs_boundary;
    std::vector<uint8_t> via_dummy;
    std::vector<uint32_t> piece_of;
    std::vector<uint32_t> piece_members;
    std::vector<uint32_t> piece_begin, piece_end, piece_cube;
    std::vector<uint32_t> piece_edge_begin, piece_edge_end;
    std::vector<std::vector<uint32_t>> cube_piece_edges;
    std::vector<std::vector<uint32_t>> cube_refused;
    std::vector<std::vector<uint32_t>> cube_cross;
    std::vector<uint32_t> cross_begin, cross_end;
    std::vector<uint32_t> node_degree;
    std::vector<uint32_t> mask_local;
    std::vector<uint64_t> root_stamp;
    std::vector<uint32_t> piece_index_of_root;
    uint64_t epoch{0};
    std::vector<uint32_t> pbase, pcount, tbase, tcount;
    uint32_t n_pieces_total{0};
    uint32_t n_tree_nodes{0};
    std::vector<uint32_t> fill_cursor;
    // quotient + greedy scratch
    std::vector<uint64_t> qkey;
    std::vector<uint32_t> qcount;
    std::vector<uint32_t> qslots;
    size_t qcap{0};
    std::vector<std::vector<QuotientAdj>> adj;
    std::vector<uint32_t> adj_touched;
    std::vector<uint64_t> mark_stamp;
    std::vector<uint32_t> mark_pos;
    uint64_t mark_epoch{0};
    std::vector<MergeCandidate> heap;
    // tree
    std::vector<uint32_t> tree_parent, tree_cost, tree_depth, tree_cube, tree_piece, tuf;
    std::vector<uint8_t> tree_alive, tree_dead;
    std::vector<uint32_t> child_begin, child_end, child_pool;
    std::vector<uint32_t> interior_begin, interior_end, interior_pool;
    std::vector<uint32_t> job_of_node;
    std::vector<uint64_t> defect_stamp;
    uint64_t defect_epoch{0};
    std::vector<uint32_t> cross_index_lo, cross_index_hi;
    std::vector<uint8_t> built;
    fus::PieceBuildScratch build_scratch;
    fus::PieceBuild build_input;
    uint32_t last_n{0};
    std::vector<uint8_t> cube_truncated;

    // ---- §4.3 the global seam check
    std::vector<uint32_t> pf_parent;
    struct SeamRec {
        uint32_t root;
        uint32_t seam;
        uint32_t edge;
    };
    std::vector<SeamRec> seam_recs;
    struct SeamJob {
        uint32_t comp;
        uint32_t seam;
        uint32_t begin;
        uint32_t end;
        bool first_in_chain;
    };
    std::vector<SeamJob> seam_jobs;
    std::vector<uint32_t> seam_job_edges;

    // ---- jobs
    std::vector<Job> jobs;
    std::vector<uint32_t> dep_pool;
    std::vector<uint32_t> bind_pool;
    std::vector<std::pair<uint32_t, uint32_t>> moved_pool;
    uint32_t job_g_seam_check{NONE};
    std::vector<uint32_t> job_mgr_scatter, job_mgr_uf, job_fb_phase1;
    std::vector<uint32_t> combine_deps;

    // ---- phase 2
    std::vector<uint8_t> is_g;
    std::vector<uint32_t> g_defects;
    std::vector<uint32_t> walk_stack;
    std::vector<uint32_t> scratch;
    std::vector<uint32_t> scratch2;
    std::vector<uint32_t> new_outer;
    std::vector<uint32_t> pending;
    std::vector<uint32_t> p1_defects;
    std::vector<uint32_t> promoted_outer;
    std::vector<uint8_t> cube_seen;
    Teardown teardown;
    Dijkstra dij, dij2;
    std::vector<Group> groups;
    std::vector<PhaseTwoRow> p2_rows;
    std::vector<uint8_t> p_of_group;
    struct ObstacleRec {
        uint32_t b;
        uint32_t inner_begin, inner_end;
        uint32_t shell_begin, shell_end;
    };
    std::vector<ObstacleRec> obstacle_recs;
    std::vector<uint32_t> obstacle_inner;
    std::vector<std::pair<uint32_t, int64_t>> obstacle_shell;
    std::vector<uint64_t> fb_phase1_ticks;

    // ---- outputs of the shot
    ShotRow row;
    pm::MatchingResult result;
    bool conflict{false};
    // debug counters
    uint64_t dbg_invariant_violations{0};
    uint64_t dbg_cap_violations{0};
    uint64_t dbg_injection_violations{0};
    uint64_t dbg_feasibility_violations{0};
    uint64_t dbg_g_unsettled{0};

    void configure(size_t n_max) {
        size_t tree_max = 3 * n_max + 2;
        det_to_h.assign(corpus->g.n, NONE);
        cube_of_h.assign(n_max, 0);
        comp_parent.assign(n_max, 0);
        comp_of.assign(n_max, 0);
        comp_begin.assign(n_max + 1, 0);
        comp_members.assign(n_max, 0);
        comp_cube_lo.assign(n_max, 0);
        comp_cube_hi.assign(n_max, 0);
        comp_claimed.assign(n_max, -1);
        comp_unfinished.assign(n_max, 0);
        comp_truncated.assign(n_max, 0);
        comp_needs_h.assign(n_max, 0);
        comp_jobs.assign(n_max, {});
        parent.assign(n_max, 0);
        set_size.assign(n_max, 0);
        w_boundary.assign(n_max, WEIGHT_INF);
        bnd.assign(n_max, WEIGHT_INF);
        obs_boundary.assign(n_max, 0);
        via_dummy.assign(n_max, 0);
        piece_of.assign(n_max, 0);
        piece_members.assign(n_max, 0);
        piece_begin.assign(n_max + 1, 0);
        piece_end.assign(n_max + 1, 0);
        piece_cube.assign(n_max + 1, 0);
        piece_edge_begin.assign(n_max + 1, 0);
        piece_edge_end.assign(n_max + 1, 0);
        cross_begin.assign(n_max, 0);
        cross_end.assign(n_max, 0);
        node_degree.assign(n_max, 0);
        mask_local.assign(n_max, 0);
        root_stamp.assign(n_max, 0);
        piece_index_of_root.assign(n_max, 0);
        fill_cursor.assign(tree_max + 1, 0);
        adj.assign(tree_max, {});
        mark_stamp.assign(tree_max, 0);
        mark_pos.assign(tree_max, 0);
        tree_parent.assign(tree_max, NONE);
        tree_cost.assign(tree_max, 0);
        tree_depth.assign(tree_max, 0);
        tree_cube.assign(tree_max, 0);
        tree_piece.assign(tree_max, NONE);
        tuf.assign(tree_max, 0);
        tree_alive.assign(tree_max, 0);
        tree_dead.assign(tree_max, 0);
        child_begin.assign(tree_max, 0);
        child_end.assign(tree_max, 0);
        interior_begin.assign(tree_max, 0);
        interior_end.assign(tree_max, 0);
        job_of_node.assign(tree_max, NONE);
        defect_stamp.assign(n_max, 0);
        built.assign(n_max, 0);
        build_scratch.reserve(n_max, 4 * INITIAL_EDGE_CAPACITY);
        pf_parent.assign(n_max + 1, 0);
        is_g.assign(n_max, 0);
        size_quotient_table(INITIAL_EDGE_CAPACITY);
        dij.configure(corpus->g.n);
        dij2.configure(corpus->g.n);
        cube_seen.assign(corpus->g.num_cubes, 0);
    }

    void size_quotient_table(size_t m) {
        size_t want = 4;
        while (want < 4 * (m + 1))
            want <<= 1;
        if (want > qcap) {
            qcap = want;
            qkey.assign(qcap, UINT64_MAX);
            qcount.assign(qcap, 0);
        }
    }

    inline uint32_t find(uint32_t v) {
        while (parent[v] != v) {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        return v;
    }
    inline uint32_t find_tree(uint32_t t) {
        while (tuf[t] != t) {
            tuf[t] = tuf[tuf[t]];
            t = tuf[t];
        }
        return t;
    }
    inline uint32_t find_pf(uint32_t p) {
        while (pf_parent[p] != p) {
            pf_parent[p] = pf_parent[pf_parent[p]];
            p = pf_parent[p];
        }
        return p;
    }
    inline uint32_t find_comp(uint32_t v) {
        while (comp_parent[v] != v) {
            comp_parent[v] = comp_parent[comp_parent[v]];
            v = comp_parent[v];
        }
        return v;
    }
    inline uint32_t tree_of_piece(uint32_t p) const {
        uint32_t c = piece_cube[p];
        return tbase[c] + (p - pbase[c]);
    }

    pm::Mwpm& hmwpm() {
        return corpus->shared->mwpm;
    }

    // -----------------------------------------------------------------------------------------
    // Jobs.
    // -----------------------------------------------------------------------------------------

    uint32_t add_job(Kind kind, CoreClass cls, int32_t cube, const std::vector<uint32_t>& deps) {
        Job job;
        job.kind = kind;
        job.cls = cls;
        job.cube = cube;
        job.dep_begin = (uint32_t)dep_pool.size();
        for (uint32_t d : deps)
            if (d != NONE)
                dep_pool.push_back(d);
        job.dep_end = (uint32_t)dep_pool.size();
        job.bind_begin = job.bind_end = (uint32_t)bind_pool.size();
        job.moved_begin = job.moved_end = (uint32_t)moved_pool.size();
        job.priority = UINT64_MAX;
        jobs.push_back(job);
        return (uint32_t)jobs.size() - 1;
    }
    void set_bind(uint32_t id, const std::vector<uint32_t>& cands) {
        jobs[id].bind_begin = (uint32_t)bind_pool.size();
        for (uint32_t c : cands)
            bind_pool.push_back(c);
        jobs[id].bind_end = (uint32_t)bind_pool.size();
    }

    // -----------------------------------------------------------------------------------------
    // The input, outside every timed region (§3: the accelerator stand-in).
    // -----------------------------------------------------------------------------------------

    pm::obs_int boundary_obs_mask(uint64_t det) const {
        const BallTables& t = corpus->decoder->tables;
        pm::obs_int mask = 0;
        for (uint64_t at = t.bcost_mask_offsets[det]; at < t.bcost_mask_offsets[det + 1]; at++)
            mask ^= (pm::obs_int)1 << t.bcost_mask_ids[at];
        return mask;
    }
    pm::obs_int edge_obs_mask(uint64_t entry) const {
        const BallTables& t = corpus->decoder->tables;
        pm::obs_int mask = 0;
        for (uint64_t at = t.ball_mask_offsets[entry]; at < t.ball_mask_offsets[entry + 1]; at++)
            mask ^= (pm::obs_int)1 << t.ball_mask_ids[at];
        return mask;
    }

    void prepare_input(const std::vector<uint64_t>& shot) {
        BallDecoder& decoder = *corpus->decoder;
        decoder.compute_seeded_detection_events(shot, decoder.seeded_scratch);
        build_ball_graph(decoder.tables, decoder.seeded_scratch, decoder.horizon, decoder.arena, decoder.config.mode,
                         nullptr, nullptr);
        h = &decoder.arena.graph;
        n = (uint32_t)h->num_nodes();
        const StripGraph& g = corpus->g;
        hoff.assign(N + 1, 0);
        for (uint32_t v = 0; v < n; v++) {
            uint32_t det = (uint32_t)h->h_to_det[v];
            det_to_h[det] = v;
            cube_of_h[v] = g.cube[det];
            hoff[cube_of_h[v] + 1]++;
        }
        for (size_t c = 0; c < N; c++)
            hoff[c + 1] += hoff[c];
        cube_edges.resize(N);
        seam_edges.resize(N);
        cube_refused.resize(N);
        cube_piece_edges.resize(N);
        cube_cross.resize(N);
        for (size_t c = 0; c < N; c++) {
            cube_edges[c].clear();
            seam_edges[c].clear();
        }
        for (uint32_t e = 0; e < (uint32_t)h->edges.size(); e++) {
            const BallGraphEdge& edge = h->edges[e];
            uint32_t ci = cube_of_h[edge.i];
            uint32_t cj = cube_of_h[edge.j];
            if (ci == cj) {
                cube_edges[ci].push_back(e);
            } else {
                if (ci + 1 != cj && cj + 1 != ci)
                    throw std::runtime_error("an H edge joins two non-adjacent cubes; 2T exceeds a cube's width");
                seam_edges[std::min(ci, cj)].push_back(e);
            }
        }
        if (sort_by_weight) {
            for (size_t c = 0; c < N; c++)
                std::sort(cube_edges[c].begin(), cube_edges[c].end(), [&](uint32_t a, uint32_t b) {
                    const BallGraphEdge& x = h->edges[a];
                    const BallGraphEdge& y = h->edges[b];
                    if (x.w_int != y.w_int)
                        return x.w_int < y.w_int;
                    if (x.i != y.i)
                        return x.i < y.i;
                    return x.j < y.j;
                });
        }
        bnd_lo.assign(N, 0);
        bnd_hi.assign(N, 0);
        {
            size_t at = 0;
            for (size_t c = 0; c < N; c++) {
                while (at < h->boundary_edges.size() && h->boundary_edges[at].i < hoff[c])
                    at++;
                bnd_lo[c] = (uint32_t)at;
                while (at < h->boundary_edges.size() && h->boundary_edges[at].i < hoff[c + 1])
                    at++;
                bnd_hi[c] = (uint32_t)at;
            }
        }
        for (const BallBoundaryEdge& edge : h->boundary_edges)
            obs_boundary[edge.i] = boundary_obs_mask(edge.det);
        edge_obs.resize(h->edges.size());
        for (size_t e = 0; e < h->edges.size(); e++)
            edge_obs[e] = edge_obs_mask(h->edges[e].entry);
        if (identity.size() < h->edges.size() + 1) {
            size_t old = identity.size();
            identity.resize(std::max<size_t>(h->edges.size() + 1, 2 * old));
            for (size_t i = old; i < identity.size(); i++)
                identity[i] = (uint32_t)i;
        }
        cross_index_lo.resize(std::max<size_t>(1, h->edges.size()));
        cross_index_hi.resize(std::max<size_t>(1, h->edges.size()));

        // The strip's components, for claims and labels. Bookkeeping only; nothing on a timed path
        // reads it except as a name.
        for (uint32_t v = 0; v < n; v++)
            comp_parent[v] = v;
        for (const BallGraphEdge& edge : h->edges) {
            uint32_t a = find_comp(edge.i);
            uint32_t b = find_comp(edge.j);
            if (a != b)
                comp_parent[std::max(a, b)] = std::min(a, b);
        }
        n_comps = 0;
        for (uint32_t v = 0; v < n; v++) {
            uint32_t r = find_comp(v);
            if (r == v)
                comp_of[v] = n_comps++;
            else
                comp_of[v] = comp_of[r];
        }
        for (uint32_t c = 0; c <= n_comps; c++)
            comp_begin[c] = 0;
        for (uint32_t v = 0; v < n; v++)
            comp_begin[comp_of[v] + 1]++;
        for (uint32_t c = 0; c < n_comps; c++)
            comp_begin[c + 1] += comp_begin[c];
        for (uint32_t c = 0; c < n_comps; c++) {
            fill_cursor[c] = comp_begin[c];
            comp_cube_lo[c] = UINT32_MAX;
            comp_cube_hi[c] = 0;
            comp_claimed[c] = -1;
            comp_unfinished[c] = 0;
            comp_truncated[c] = 0;
            comp_needs_h[c] = 0;
            comp_jobs[c].clear();
        }
        for (uint32_t v = 0; v < n; v++) {
            uint32_t c = comp_of[v];
            comp_members[fill_cursor[c]++] = v;
            comp_cube_lo[c] = std::min(comp_cube_lo[c], cube_of_h[v]);
            comp_cube_hi[c] = std::max(comp_cube_hi[c], cube_of_h[v]);
        }
    }

    void mark_all_unbuilt(uint32_t upto) {
        pm::MatchingGraph& graph = hmwpm().flooder.graph;
        for (uint32_t v = 0; v < upto; v++) {
            pm::DetectorNode& node = graph.nodes[v];
            node.neighbors.clear();
            node.neighbor_weights.clear();
            node.neighbor_observables.clear();
            node.neighbor_implied_weights.clear();
            built[v] = 0;
        }
    }

    // -----------------------------------------------------------------------------------------
    // §4.1. The per-cube manager: Steps 1-3 of §4.3 of the paper on `H_c`, with the cross-seam
    // edges treated as Step 1 crossing edges.
    // -----------------------------------------------------------------------------------------

    void mgr_uf(uint32_t c) {
        uint32_t lo = hoff[c];
        uint32_t hi = hoff[c + 1];
        uint32_t nc = hi - lo;
        uint32_t threshold = (uint32_t)((nc + k - 1) / k);
        if (threshold == 0)
            threshold = 1;
        for (uint32_t v = lo; v < hi; v++) {
            parent[v] = v;
            set_size[v] = 1;
            w_boundary[v] = WEIGHT_INF;
            bnd[v] = WEIGHT_INF;
            via_dummy[v] = 0;
        }
        for (uint32_t at = bnd_lo[c]; at < bnd_hi[c]; at++) {
            const BallBoundaryEdge& edge = h->boundary_edges[at];
            w_boundary[edge.i] = edge.w_int;
            bnd[edge.i] = edge.w_int;
        }
        std::vector<uint32_t>& refused = cube_refused[c];
        refused.clear();
        for (uint32_t e : cube_edges[c]) {
            const BallGraphEdge& edge = h->edges[e];
            uint32_t a = find(edge.i);
            uint32_t b = find(edge.j);
            if (a == b)
                continue;
            if ((uint64_t)set_size[a] + set_size[b] <= threshold) {
                if (set_size[a] < set_size[b] || (set_size[a] == set_size[b] && b < a))
                    std::swap(a, b);
                parent[b] = a;
                set_size[a] += set_size[b];
                continue;
            }
            refused.push_back(e);
            pm::weight_int half = fus::even_dummy_weight(edge.w_int);
            if (half < bnd[edge.i]) {
                bnd[edge.i] = half;
                via_dummy[edge.i] = 1;
            }
            if (half < bnd[edge.j]) {
                bnd[edge.j] = half;
                via_dummy[edge.j] = 1;
            }
        }
        // The cross-seam H-edges: a dummy on this cube's own endpoint, and only on it.
        for (int s : {(int)c - 1, (int)c}) {
            if (s < 0 || s + 1 >= (int)N)
                continue;
            for (uint32_t e : seam_edges[s]) {
                const BallGraphEdge& edge = h->edges[e];
                uint32_t u = cube_of_h[edge.i] == c ? edge.i : edge.j;
                pm::weight_int half = fus::even_dummy_weight(edge.w_int);
                if (half < bnd[u]) {
                    bnd[u] = half;
                    via_dummy[u] = 1;
                }
            }
        }
        // Pieces, numbered in ascending minimum member, after every earlier cube's.
        epoch++;
        uint32_t base = n_pieces_total;
        pbase[c] = base;
        uint32_t count = 0;
        for (uint32_t v = lo; v < hi; v++) {
            uint32_t root = find(v);
            if (root_stamp[root] != epoch) {
                root_stamp[root] = epoch;
                piece_index_of_root[root] = base + count;
                count++;
            }
            piece_of[v] = piece_index_of_root[root];
        }
        pcount[c] = count;
        for (uint32_t p = base; p <= base + count; p++)
            piece_begin[p] = 0;
        for (uint32_t v = lo; v < hi; v++)
            piece_begin[piece_of[v] + 1]++;
        piece_begin[base] = lo;
        for (uint32_t p = base; p < base + count; p++) {
            piece_begin[p + 1] += piece_begin[p];
            fill_cursor[p] = piece_begin[p];
            piece_cube[p] = c;
        }
        for (uint32_t p = base; p < base + count; p++)
            piece_end[p] = piece_begin[p + 1];
        for (uint32_t v = lo; v < hi; v++)
            piece_members[fill_cursor[piece_of[v]]++] = v;
        // The quotient graph on the final piece ids.
        size_quotient_table(refused.size());
        qslots.clear();
        for (uint32_t e : refused) {
            const BallGraphEdge& edge = h->edges[e];
            uint32_t pa = piece_of[edge.i];
            uint32_t pb = piece_of[edge.j];
            uint32_t a = std::min(pa, pb);
            uint32_t b = std::max(pa, pb);
            uint64_t key = ((uint64_t)a << 32) | b;
            size_t slot = (size_t)((key * 1099511628211ull) >> 24) & (qcap - 1);
            while (qkey[slot] != UINT64_MAX && qkey[slot] != key)
                slot = (slot + 1) & (qcap - 1);
            if (qkey[slot] == UINT64_MAX) {
                qkey[slot] = key;
                qcount[slot] = 0;
                qslots.push_back((uint32_t)slot);
            }
            qcount[slot]++;
        }
        n_pieces_total += count;
    }

    /// Step 2: cheapest-boundary-first fusion tree over the cube's pieces, then each refused edge
    /// assigned to the lowest common ancestor of its two leaves.
    void mgr_tree(uint32_t c) {
        uint32_t base = tbase[c] = n_tree_nodes;
        uint32_t P = pcount[c];
        uint32_t next = base + P;
        for (uint32_t i = 0; i < P; i++) {
            uint32_t t = base + i;
            tree_parent[t] = NONE;
            tree_cost[t] = 0;
            tree_alive[t] = 1;
            tree_dead[t] = 0;
            tree_cube[t] = c;
            tree_piece[t] = pbase[c] + i;
            tuf[t] = t;
            child_begin[t] = child_end[t] = (uint32_t)child_pool.size();
            interior_begin[t] = interior_end[t] = 0;
        }
        for (uint32_t t : adj_touched)
            adj[t].clear();
        adj_touched.clear();
        heap.clear();
        for (uint32_t slot : qslots) {
            uint64_t key = qkey[slot];
            uint32_t a = base + ((uint32_t)(key >> 32) - pbase[c]);
            uint32_t b = base + ((uint32_t)(key & 0xffffffffu) - pbase[c]);
            uint32_t count = qcount[slot];
            if (adj[a].empty())
                adj_touched.push_back(a);
            if (adj[b].empty())
                adj_touched.push_back(b);
            adj[a].push_back(QuotientAdj{b, count});
            adj[b].push_back(QuotientAdj{a, count});
            heap.push_back(MergeCandidate{count, a, b});
            qkey[slot] = UINT64_MAX;
        }
        qslots.clear();
        std::make_heap(heap.begin(), heap.end());
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end());
            MergeCandidate cand = heap.back();
            heap.pop_back();
            if (!tree_alive[cand.lo] || !tree_alive[cand.hi])
                continue;
            uint32_t a = cand.lo;
            uint32_t b = cand.hi;
            uint32_t p = next++;
            tree_parent[p] = NONE;
            tree_parent[a] = p;
            tree_parent[b] = p;
            tree_cost[p] = cand.count;
            tree_alive[a] = 0;
            tree_alive[b] = 0;
            tree_alive[p] = 1;
            tree_dead[p] = 0;
            tree_cube[p] = c;
            tree_piece[p] = NONE;
            tuf[p] = p;
            tuf[a] = p;
            tuf[b] = p;
            child_begin[p] = (uint32_t)child_pool.size();
            child_pool.push_back(a);
            child_pool.push_back(b);
            child_end[p] = (uint32_t)child_pool.size();
            std::vector<QuotientAdj>& big = adj[a].size() >= adj[b].size() ? adj[a] : adj[b];
            std::vector<QuotientAdj>& small = adj[a].size() >= adj[b].size() ? adj[b] : adj[a];
            if (adj[p].empty())
                adj_touched.push_back(p);
            adj[p].swap(big);
            mark_epoch++;
            size_t write = 0;
            for (size_t i = 0; i < adj[p].size(); i++) {
                uint32_t other = find_tree(adj[p][i].other);
                if (other == p)
                    continue;
                if (mark_stamp[other] == mark_epoch) {
                    adj[p][mark_pos[other]].count += adj[p][i].count;
                    continue;
                }
                mark_stamp[other] = mark_epoch;
                mark_pos[other] = (uint32_t)write;
                adj[p][write++] = QuotientAdj{other, adj[p][i].count};
            }
            adj[p].resize(write);
            for (const QuotientAdj& entry : small) {
                uint32_t other = find_tree(entry.other);
                if (other == p)
                    continue;
                if (mark_stamp[other] == mark_epoch) {
                    adj[p][mark_pos[other]].count += entry.count;
                    continue;
                }
                mark_stamp[other] = mark_epoch;
                mark_pos[other] = (uint32_t)adj[p].size();
                adj[p].push_back(QuotientAdj{other, entry.count});
            }
            big.clear();
            small.clear();
            for (const QuotientAdj& entry : adj[p]) {
                heap.push_back(MergeCandidate{entry.count, std::min(p, entry.other), std::max(p, entry.other)});
                std::push_heap(heap.begin(), heap.end());
            }
        }
        tcount[c] = next - base;
        n_tree_nodes = next;
        for (uint32_t t = base; t < next; t++)
            tree_depth[t] = 0;
        for (uint32_t t = next; t-- > base + P;) {
            uint32_t depth = tree_depth[t] + 1;
            for (uint32_t at = child_begin[t]; at < child_end[t]; at++)
                tree_depth[child_pool[at]] = depth;
        }
        auto lca = [&](uint32_t x, uint32_t y) -> uint32_t {
            while (x != NONE && y != NONE && tree_depth[x] > tree_depth[y])
                x = tree_parent[x];
            while (x != NONE && y != NONE && tree_depth[y] > tree_depth[x])
                y = tree_parent[y];
            while (x != NONE && y != NONE && x != y) {
                x = tree_parent[x];
                y = tree_parent[y];
            }
            if (x == NONE || y == NONE || x != y)
                throw std::runtime_error("fusion tree: a crossing edge's pieces have no common ancestor");
            return x;
        };
        const std::vector<uint32_t>& refused = cube_refused[c];
        uint32_t pool_base = (uint32_t)interior_pool.size();
        for (uint32_t t = base; t <= next; t++)
            fill_cursor[t - base] = 0;
        for (uint32_t e : refused) {
            uint32_t t = lca(tree_of_piece(piece_of[h->edges[e].i]), tree_of_piece(piece_of[h->edges[e].j]));
            fill_cursor[t - base + 1]++;
        }
        for (uint32_t t = base; t < next; t++)
            fill_cursor[t - base + 1] += fill_cursor[t - base];
        for (uint32_t t = base; t < next; t++) {
            interior_begin[t] = pool_base + fill_cursor[t - base];
            interior_end[t] = interior_begin[t];
        }
        interior_pool.resize(pool_base + refused.size());
        for (uint32_t e : refused) {
            uint32_t t = lca(tree_of_piece(piece_of[h->edges[e].i]), tree_of_piece(piece_of[h->edges[e].j]));
            interior_pool[interior_end[t]++] = e;
        }
    }

    /// Step 3: the cube's edge list grouped by piece, each defect's degree, crossing-edge lists and
    /// local mask offsets. Reads no ball table and touches no solver node.
    void mgr_scatter(uint32_t c) {
        uint32_t lo = hoff[c];
        uint32_t hi = hoff[c + 1];
        uint32_t pb = pbase[c];
        uint32_t P = pcount[c];
        for (uint32_t v = lo; v < hi; v++)
            node_degree[v] = bnd[v] == WEIGHT_INF ? 0 : 1;
        for (uint32_t p = pb; p <= pb + P; p++)
            fill_cursor[p - pb] = 0;
        for (uint32_t e : cube_edges[c]) {
            const BallGraphEdge& edge = h->edges[e];
            uint32_t a = piece_of[edge.i];
            uint32_t b = piece_of[edge.j];
            fill_cursor[a - pb + 1]++;
            if (a != b)
                fill_cursor[b - pb + 1]++;
            node_degree[edge.i]++;
            node_degree[edge.j]++;
        }
        for (int s : {(int)c - 1, (int)c}) {
            if (s < 0 || s + 1 >= (int)N)
                continue;
            for (uint32_t e : seam_edges[s]) {
                const BallGraphEdge& edge = h->edges[e];
                uint32_t u = cube_of_h[edge.i] == c ? edge.i : edge.j;
                fill_cursor[piece_of[u] - pb + 1]++;
                node_degree[u]++;
            }
        }
        for (uint32_t i = 0; i < P; i++)
            fill_cursor[i + 1] += fill_cursor[i];
        std::vector<uint32_t>& pe = cube_piece_edges[c];
        pe.resize(fill_cursor[P]);
        for (uint32_t i = 0; i < P; i++) {
            piece_edge_begin[pb + i] = fill_cursor[i];
            piece_edge_end[pb + i] = fill_cursor[i];
        }
        for (uint32_t e : cube_edges[c]) {
            const BallGraphEdge& edge = h->edges[e];
            uint32_t a = piece_of[edge.i];
            uint32_t b = piece_of[edge.j];
            pe[piece_edge_end[a]++] = e;
            if (a != b)
                pe[piece_edge_end[b]++] = e;
        }
        for (int s : {(int)c - 1, (int)c}) {
            if (s < 0 || s + 1 >= (int)N)
                continue;
            for (uint32_t e : seam_edges[s]) {
                const BallGraphEdge& edge = h->edges[e];
                uint32_t u = cube_of_h[edge.i] == c ? edge.i : edge.j;
                pe[piece_edge_end[piece_of[u]]++] = e;
            }
        }
        // The crossing edges at each defect: intra-cube refused edges and seam edges.
        for (uint32_t v = lo; v < hi; v++)
            cross_begin[v] = 0;
        std::vector<uint32_t>& cr = cube_cross[c];
        uint32_t total = 0;
        auto count_at = [&](uint32_t u) {
            cross_end[u]++;
        };
        for (uint32_t v = lo; v < hi; v++)
            cross_end[v] = 0;
        for (uint32_t e : cube_refused[c]) {
            count_at(h->edges[e].i);
            count_at(h->edges[e].j);
        }
        for (int s : {(int)c - 1, (int)c}) {
            if (s < 0 || s + 1 >= (int)N)
                continue;
            for (uint32_t e : seam_edges[s]) {
                const BallGraphEdge& edge = h->edges[e];
                count_at(cube_of_h[edge.i] == c ? edge.i : edge.j);
            }
        }
        for (uint32_t v = lo; v < hi; v++) {
            cross_begin[v] = total;
            total += cross_end[v];
            cross_end[v] = cross_begin[v];
        }
        cr.resize(total);
        for (uint32_t e : cube_refused[c]) {
            cr[cross_end[h->edges[e].i]++] = e;
            cr[cross_end[h->edges[e].j]++] = e;
        }
        for (int s : {(int)c - 1, (int)c}) {
            if (s < 0 || s + 1 >= (int)N)
                continue;
            for (uint32_t e : seam_edges[s]) {
                const BallGraphEdge& edge = h->edges[e];
                uint32_t u = cube_of_h[edge.i] == c ? edge.i : edge.j;
                cr[cross_end[u]++] = e;
            }
        }
        // Local mask offsets; the per-cube base is a fixed arena address in a real system.
        uint32_t running = 0;
        for (uint32_t v = lo; v < hi; v++) {
            mask_local[v] = running;
            running += node_degree[v];
        }
        for (uint32_t t = tbase[c]; t < tbase[c] + tcount[c]; t++)
            tuf[t] = t;
    }

    /// The mask arena and the build's view of the shot, rebased across cubes. Untimed: these are
    /// fixed arena addresses, not work.
    void bind_shot() {
        fus::FusionInstance& fusion = corpus->fusion;
        // Exactly `n + 1`: `EdgeMask::total_entries` reads `offsets.back()`, so a longer vector would
        // report a stale total and the build scratch would be sized from it. `attach` reserved
        // `n_max + 1`, so this never reallocates under the flooder's pointer.
        fusion.mask.offsets.resize((size_t)n + 1);
        uint32_t base = 0;
        for (size_t c = 0; c < N; c++) {
            uint32_t cube_total = 0;
            for (uint32_t v = hoff[c]; v < hoff[c + 1]; v++) {
                fusion.mask.offsets[v] = base + mask_local[v];
                cube_total = mask_local[v] + node_degree[v];
            }
            base += cube_total;
        }
        fusion.mask.offsets[n] = base;
        size_t bits_capacity = fusion.mask.bits.capacity();
        fusion.mask.bits.resize(std::max<size_t>(base, 1));
        if (fusion.mask.bits.capacity() != bits_capacity)
            fusion.mask.grow_events++;
        hmwpm().flooder.edge_mask = fusion.mask.bits.data();
        hmwpm().flooder.edge_mask_offsets = fusion.mask.offsets.data();
        build_scratch.size_for_shot(fusion.mask.total_entries(), n);
        build_input.edges = h->edges.data();
        build_input.edge_obs = edge_obs.data();
        build_input.piece_of = piece_of.data();
        build_input.boundary_weight = bnd.data();
        build_input.boundary_is_dummy = via_dummy.data();
        build_input.boundary_obs = obs_boundary.data();
        build_input.refused_slot = identity.data();
        build_input.cross_index_lo = cross_index_lo.data();
        build_input.cross_index_hi = cross_index_hi.data();
        build_input.built = built.data();
    }

    // -----------------------------------------------------------------------------------------
    // §4.3. The global manager's seam check: join the cubes' pieces across seams and find the
    // seam components and the spanning chains.
    // -----------------------------------------------------------------------------------------

    void g_seam_check() {
        for (uint32_t p = 0; p < n_pieces_total; p++)
            pf_parent[p] = p;
        auto unite = [&](uint32_t a, uint32_t b) {
            a = find_pf(a);
            b = find_pf(b);
            if (a != b)
                pf_parent[std::max(a, b)] = std::min(a, b);
        };
        for (size_t c = 0; c < N; c++)
            for (uint32_t e : cube_refused[c])
                unite(piece_of[h->edges[e].i], piece_of[h->edges[e].j]);
        seam_recs.clear();
        for (size_t s = 0; s + 1 < N; s++)
            for (uint32_t e : seam_edges[s])
                unite(piece_of[h->edges[e].i], piece_of[h->edges[e].j]);
        for (size_t s = 0; s + 1 < N; s++)
            for (uint32_t e : seam_edges[s])
                seam_recs.push_back(SeamRec{find_pf(piece_of[h->edges[e].i]), (uint32_t)s, e});
        std::sort(seam_recs.begin(), seam_recs.end(), [](const SeamRec& a, const SeamRec& b) {
            if (a.root != b.root)
                return a.root < b.root;
            if (a.seam != b.seam)
                return a.seam < b.seam;
            return a.edge < b.edge;
        });
        seam_jobs.clear();
        seam_job_edges.clear();
        row.n_seam_components = 0;
        row.n_spanning = 0;
        row.max_span = 0;
        size_t at = 0;
        while (at < seam_recs.size()) {
            size_t end = at;
            uint32_t span = 0;
            bool first = true;
            while (end < seam_recs.size() && seam_recs[end].root == seam_recs[at].root) {
                size_t seam_end = end;
                while (seam_end < seam_recs.size() && seam_recs[seam_end].root == seam_recs[at].root &&
                       seam_recs[seam_end].seam == seam_recs[end].seam)
                    seam_end++;
                SeamJob job;
                job.comp = comp_of[h->edges[seam_recs[end].edge].i];
                job.seam = seam_recs[end].seam;
                job.begin = (uint32_t)seam_job_edges.size();
                for (size_t i = end; i < seam_end; i++)
                    seam_job_edges.push_back(seam_recs[i].edge);
                job.end = (uint32_t)seam_job_edges.size();
                job.first_in_chain = first;
                first = false;
                seam_jobs.push_back(job);
                span++;
                end = seam_end;
            }
            row.n_seam_components++;
            if (span >= 2)
                row.n_spanning++;
            row.max_span = std::max(row.max_span, span);
            at = end;
        }
        row.seam_depth = row.max_span;
    }

    // -----------------------------------------------------------------------------------------
    // §4.4. A fusion, intra-cube or across a seam: one generic routine over a node's children.
    // -----------------------------------------------------------------------------------------

    void collect_defects(uint32_t node, std::vector<uint32_t>& out) {
        walk_stack.clear();
        walk_stack.push_back(node);
        while (!walk_stack.empty()) {
            uint32_t t = walk_stack.back();
            walk_stack.pop_back();
            if (tree_piece[t] != NONE) {
                uint32_t p = tree_piece[t];
                for (uint32_t i = piece_begin[p]; i < piece_end[p]; i++)
                    out.push_back(piece_members[i]);
                continue;
            }
            for (uint32_t at = child_end[t]; at-- > child_begin[t];)
                walk_stack.push_back(child_pool[at]);
        }
    }

    TimelineStatus fuse(uint32_t p) {
        pm::Mwpm& mwpm = hmwpm();
        fus::FusionInstance& fusion = corpus->fusion;
        fus::align_clock_for_new_growth(mwpm);
        for (uint32_t at = child_begin[p]; at < child_end[p]; at++)
            tuf[child_pool[at]] = p;
        for (uint32_t at = interior_begin[p]; at < interior_end[p]; at++) {
            uint32_t e = interior_pool[at];
            const BallGraphEdge& edge = h->edges[e];
            fus::unmask_crossing_edge(fusion, edge.i, cross_index_lo[e], edge.j, cross_index_hi[e]);
        }
        defect_epoch++;
        for (uint32_t at = interior_begin[p]; at < interior_end[p]; at++) {
            const BallGraphEdge& edge = h->edges[interior_pool[at]];
            for (uint32_t u : {edge.i, edge.j}) {
                if (defect_stamp[u] == defect_epoch || !via_dummy[u])
                    continue;
                defect_stamp[u] = defect_epoch;
                const std::vector<uint32_t>& cr = cube_cross[cube_of_h[u]];
                pm::weight_int wmin = WEIGHT_INF;
                for (uint32_t at2 = cross_begin[u]; at2 < cross_end[u]; at2++) {
                    const BallGraphEdge& other = h->edges[cr[at2]];
                    uint32_t far = other.i == u ? other.j : other.i;
                    if (find_tree(tree_of_piece(piece_of[far])) == p)
                        continue;
                    wmin = std::min(wmin, other.w_int);
                }
                pm::weight_int half = wmin == WEIGHT_INF ? WEIGHT_INF : fus::even_dummy_weight(wmin);
                pm::weight_int want = std::min(w_boundary[u], half);
                bool dummy = want < w_boundary[u];
                pm::weight_int previous = bnd[u];
                bool was_dummy = via_dummy[u] != 0;
                bnd[u] = want;
                via_dummy[u] = dummy ? 1 : 0;
                if (want == previous && dummy == was_dummy)
                    continue;
                pm::weight_int weight = want == WEIGHT_INF ? fusion.unreachable_boundary_weight() : want;
                pm::obs_int obs = (dummy || want == WEIGHT_INF) ? (pm::obs_int)0 : obs_boundary[u];
                fus::set_boundary_and_release(mwpm, u, weight, obs);
            }
        }
        return fus::run_until_settled(mwpm);
    }

    /// §4.5. A truncated job's regions leave the shared instance, so no later job processes them.
    void evict(uint32_t node) {
        scratch.clear();
        collect_defects(node, scratch);
        teardown.run(hmwpm(), scratch.data(), scratch.size());
        tree_dead[node] = 1;
        for (uint32_t v : scratch) {
            cube_truncated[cube_of_h[v]] = 1;
            comp_truncated[comp_of[v]] = 1;
            comp_unfinished[comp_of[v]] = 1;
        }
    }

    void check_feasibility_after(uint32_t node) {
#ifndef NDEBUG
        scratch2.clear();
        collect_defects(node, scratch2);
        if (!fus::feasibility_holds(corpus->fusion, scratch2.data(), scratch2.size(), nullptr))
            dbg_feasibility_violations++;
#else
        (void)node;
#endif
    }

    // -----------------------------------------------------------------------------------------
    // §5.1. Phase 1, on a cube's fallback core: obstacles, promotions, the stock solve.
    // -----------------------------------------------------------------------------------------

    struct PhaseOne {
        uint32_t job{NONE};
        uint64_t ticks{0};
        uint64_t obstacle_ticks{0};
        uint32_t obstacles{0};
        uint32_t promotions{0};
        uint32_t handoff_nodes{0};
    };

    /// Neighbour defects (det ids) of the scanned cubes that `pred` accepts.
    template <typename Pred>
    void scan_neighbour_defects(uint32_t c, const CubeRegion& r, std::vector<uint32_t>& out, const Pred& pred) {
        for (uint32_t cc = r.reach_lo; cc <= r.reach_hi; cc++) {
            if (cc == c)
                continue;
            for (uint32_t v = hoff[cc]; v < hoff[cc + 1]; v++) {
                uint32_t det = (uint32_t)h->h_to_det[v];
                if (pred(det))
                    out.push_back(det);
            }
        }
    }

    /// Extends `F`'s territory by `W` around `seeds` (plus every direct neighbour of each seed, so
    /// that an extension always makes progress), returning the nodes added in `added`.
    void extend_territory(FallbackInstance& F, const uint32_t* seeds, size_t count, std::vector<uint32_t>& added,
                          std::vector<uint32_t>& outer_out) {
        const StripGraph& g = corpus->g;
        added.clear();
        dij.run(g, seeds, count, W_int, [](uint32_t) { return true; }, [&](uint32_t v, int64_t) {
            if (!F.in_omega[v])
                added.push_back(v);
        });
        for (size_t i = 0; i < count; i++)
            for (uint32_t e = g.off[seeds[i]]; e < g.off[seeds[i] + 1]; e++)
                if (!F.in_omega[g.nbr[e]] && dij.dist[g.nbr[e]] == DIST_INF)
                    added.push_back(g.nbr[e]);
        std::sort(added.begin(), added.end());
        added.erase(std::unique(added.begin(), added.end()), added.end());
        for (uint32_t v : added)
            F.add_node(v, outer_out);
        F.refresh_outer();
    }

    PhaseOne run_phase1(uint32_t c) {
        PhaseOne out;
        FallbackInstance& F = *corpus->fallback[c];
        const CubeRegion& r = corpus->regions[c];
        const StripGraph& g = corpus->g;
        pm::Mwpm& mwpm = F.mwpm;
        mwpm.flooder.outer_layer = nullptr;
        mwpm.flooder.queue.cur_time = 0;
        p1_defects.clear();
        promoted_outer.clear();

        uint64_t t0 = tick();
        // ---- obstacles (§5.1, deviation 1)
        //
        // Every candidate — a defect outside the territory and strictly within `T` of it — gets
        // exactly one Dijkstra, out to `T` plus the longest edge. What that Dijkstra finds is a
        // property of `G` and the syndrome alone, so it is kept: the defects inside the candidate's
        // `T`-ball (which decide a promotion) and its shell nodes with their distances (which decide
        // the slots). A promotion changes only which nodes are in the territory, so the records are
        // re-read rather than re-run, and the slots are written once, from the final territory.
        for (uint32_t v = hoff[c]; v < hoff[c + 1]; v++)
            p1_defects.push_back((uint32_t)h->h_to_det[v]);
        scan_neighbour_defects(c, r, p1_defects, [&](uint32_t det) { return F.in_omega[det] != 0; });
        obstacle_recs.clear();
        obstacle_inner.clear();
        obstacle_shell.clear();
        scratch.clear();
        scan_neighbour_defects(c, r, scratch, [&](uint32_t det) { return r.in_band[det] && !F.in_omega[det]; });
        auto record = [&](uint32_t b) {
            ObstacleRec rec;
            rec.b = b;
            rec.inner_begin = (uint32_t)obstacle_inner.size();
            rec.shell_begin = (uint32_t)obstacle_shell.size();
            dij.run(g, &b, 1, T_int + (int64_t)g.w_max, [](uint32_t) { return true; }, [&](uint32_t v, int64_t d) {
                if (d < T_int && v != b && det_to_h[v] != NONE)
                    obstacle_inner.push_back(v);
            });
            for (uint32_t x : dij.touched) {
                int64_t dx = dij.dist[x];
                if (dx < T_int)
                    continue;
                for (uint32_t e = g.off[x]; e < g.off[x + 1]; e++) {
                    if (dij.dist[g.nbr[e]] < T_int) {
                        obstacle_shell.push_back({x, dx});
                        break;
                    }
                }
            }
            rec.inner_end = (uint32_t)obstacle_inner.size();
            rec.shell_end = (uint32_t)obstacle_shell.size();
            obstacle_recs.push_back(rec);
        };
        for (uint32_t b : scratch)
            record(b);
        while (true) {
            // A candidate is promoted when a defect of the territory lies inside its `T`-ball.
            std::vector<uint32_t>& promote = scratch2;
            promote.clear();
            for (const ObstacleRec& rec : obstacle_recs) {
                if (F.in_omega[rec.b])
                    continue;
                for (uint32_t at = rec.inner_begin; at < rec.inner_end; at++) {
                    if (F.in_omega[obstacle_inner[at]]) {
                        promote.push_back(rec.b);
                        break;
                    }
                }
            }
            if (promote.empty())
                break;
            out.promotions += (uint32_t)promote.size();
            std::vector<uint32_t> added;
            std::vector<uint32_t> outer_now;
            extend_territory(F, promote.data(), promote.size(), added, outer_now);
            for (uint32_t v : added) {
                if (det_to_h[v] != NONE)
                    p1_defects.push_back(v);
                if (!r.in_p[v])
                    out.handoff_nodes++;
            }
            // The outer layer moved: it carries the weight-0 slot, and a node that stopped being
            // outer gets its real slot back.
            for (uint32_t v : F.touched) {
                if (!F.in_omega[v])
                    continue;
                if (F.outer[v] && F.slot_kind[v] != SLOT_OUTER)
                    F.set_slot(v, 0, 0, SLOT_OUTER);
                else if (!F.outer[v] && F.slot_kind[v] == SLOT_OUTER)
                    F.restore_real_slot(v);
            }
            for (uint32_t v : outer_now)
                promoted_outer.push_back(v);
            F.resched.clear();
            // New candidates: defects within `T` of the part of the outer layer that just appeared.
            scratch.clear();
            dij2.run(g, outer_now.data(), outer_now.size(), T_int - 1, [&](uint32_t v) { return !F.in_omega[v]; },
                     [&](uint32_t v, int64_t) {
                         if (!F.in_omega[v] && det_to_h[v] != NONE)
                             scratch.push_back(v);
                     });
            for (uint32_t b : scratch) {
                bool known = false;
                for (const ObstacleRec& rec : obstacle_recs)
                    known = known || rec.b == b;
                if (!known)
                    record(b);
            }
        }
        out.obstacles = 0;
        for (const ObstacleRec& rec : obstacle_recs) {
            if (F.in_omega[rec.b])
                continue;
            out.obstacles++;
            for (uint32_t at = rec.shell_begin; at < rec.shell_end; at++) {
                auto [x, dx] = obstacle_shell[at];
                if (!F.in_omega[x])
                    continue;
                // `d(x, b) - T` rounded down to even; a smaller slot already there only helps.
                pm::weight_int wt = fus::even_dummy_weight((pm::weight_int)(2 * (dx - T_int)));
                if (F.mask[F.moff[x]] == 0 && F.node(x).neighbor_weights[0] <= wt)
                    continue;
                F.set_slot(x, wt, 0, SLOT_OBST);
                F.virtual_slots.push_back(x);
            }
        }
        uint64_t t1 = tick();
        // ---- the stock solve on G_c ∪ B_c (and whatever promotions added)
        for (uint32_t det : p1_defects)
            inject_uncapped(mwpm, det);
        uint64_t events = 0;
        uint32_t hit = run_g_until_trigger(mwpm, events);
        (void)hit;
        uint64_t t2 = tick();
        if (!mwpm.node_arena.live.empty())
            dbg_g_unsettled++;
        out.ticks = t2 - t0;
        out.obstacle_ticks = t1 - t0;
        out.job = add_job(Kind::FB_PHASE1, CoreClass::FB, (int32_t)c, {});
        jobs[out.job].measured = out.ticks;
        jobs[out.job].n_items = (uint32_t)p1_defects.size();
        uint32_t ob = add_job(Kind::FB_OBSTACLES, CoreClass::FB, (int32_t)c, {});
        jobs[ob].measured = out.obstacle_ticks;
        jobs[ob].n_items = out.obstacles;
        jobs[ob].informational = true;
        return out;
    }

    /// Hands a phase-1 instance back: extracts (untimed) and undoes every edit.
    void discard_fallback(FallbackInstance& F, const std::vector<uint32_t>& defects) {
        if (!F.mwpm.node_arena.live.empty()) {
            F.mwpm.reset();
        } else {
            fus::extract_root(F.mwpm, defects.data(), defects.size());
            reset_for_next_shot(F.mwpm);
        }
        F.undo();
    }

    // -----------------------------------------------------------------------------------------
    // §5.2. Phase 2.
    // -----------------------------------------------------------------------------------------

    int32_t current_group{-1};

    /// Does adding `x` to group `gi`'s territory meet another group?
    bool node_conflicts(uint32_t gi, uint32_t x) {
        for (uint32_t other = 0; other < groups.size(); other++) {
            if (other == gi)
                continue;
            if (corpus->fallback[groups[other].c0]->in_omega[x])
                return true;
        }
        return false;
    }

    /// Claims component `comp` for group `gi`: tears its `H` state down and records the jobs the
    /// claiming phase-2 job must wait for. Returns false on a conflict.
    bool claim(uint32_t gi, uint32_t comp, std::vector<uint32_t>& deps, uint64_t& moved, PhaseTwoRow& p2) {
        if (comp_claimed[comp] == (int32_t)gi)
            return true;
        if (comp_claimed[comp] >= 0)
            return false;
        comp_claimed[comp] = (int32_t)gi;
        p2.n_reset_components++;
        uint32_t size = comp_begin[comp + 1] - comp_begin[comp];
        p2.n_reset_defects += size;
        const Group& grp = groups[gi];
        for (uint32_t at = comp_begin[comp]; at < comp_begin[comp + 1]; at++) {
            uint32_t v = comp_members[at];
            if (std::find(grp.cubes.begin(), grp.cubes.end(), cube_of_h[v]) == grp.cubes.end())
                moved++;
        }
        teardown.run(hmwpm(), comp_members.data() + comp_begin[comp], size);
        for (uint32_t j : comp_jobs[comp]) {
            deps.push_back(j);
            if (jobs[j].status == ST_COMPLETE)
                jobs[j].status = ST_EVICTED;
        }
        comp_needs_h[comp] = 1;
        return true;
    }

    /// H-only members of `comp`: the ones outside `Ω`.
    void h_only_members(uint32_t comp, std::vector<uint32_t>& out) {
        out.clear();
        for (uint32_t at = comp_begin[comp]; at < comp_begin[comp + 1]; at++)
            if (!is_g[comp_members[at]])
                out.push_back(comp_members[at]);
    }

    /// Re-solves a claimed component's `H`-only members at `Y = 0` in the speculative instance,
    /// with every edge into `Ω` masked. Returns the members whose dual reached `T` on `(t2)`.
    bool solve_h_only(uint32_t comp, std::vector<uint32_t>& capped) {
        pm::Mwpm& mwpm = hmwpm();
        fus::FusionInstance& fusion = corpus->fusion;
        h_only_members(comp, scratch2);
        comp_needs_h[comp] = 0;
        if (scratch2.empty())
            return true;
        pm::MatchingGraph& graph = mwpm.flooder.graph;
        for (uint32_t u : scratch2) {
            pm::DetectorNode& node = graph.nodes[u];
            for (size_t kk = 0; kk < node.neighbors.size(); kk++) {
                if (node.neighbors[kk] == nullptr) {
                    bool real = w_boundary[u] != WEIGHT_INF;
                    node.neighbor_weights[0] = real ? w_boundary[u] : fusion.unreachable_boundary_weight();
                    node.neighbor_observables[0] = real ? obs_boundary[u] : 0;
                    fusion.mask.set(u, 0, 0);
                    continue;
                }
                uint32_t v = (uint32_t)(node.neighbors[kk] - graph.nodes.data());
                fusion.mask.set(u, (uint32_t)kk, is_g[v] ? 1 : 0);
            }
        }
        TimelineStatus status = fus::solve_piece(mwpm, T_int, scratch2.data(), scratch2.size());
        if (status == TimelineStatus::COMPLETE) {
#ifndef NDEBUG
            pm::cumulative_time_int now = mwpm.flooder.queue.cur_time;
            for (uint32_t u : scratch2)
                if (graph.nodes[u].local_radius().get_distance_at_time(now) > T_int)
                    dbg_cap_violations++;
#endif
            return true;
        }
        capped.clear();
        pm::cumulative_time_int now = mwpm.flooder.queue.cur_time;
        for (uint32_t u : scratch2)
            if (graph.nodes[u].local_radius().get_distance_at_time(now) >= T_int)
                capped.push_back((uint32_t)h->h_to_det[u]);
        if (capped.empty())
            for (uint32_t u : scratch2)
                capped.push_back((uint32_t)h->h_to_det[u]);
        teardown.run(mwpm, scratch2.data(), scratch2.size());
        comp_needs_h[comp] = 1;
        return false;
    }

    /// Debug twin: §5.2's invariant, every `H`-only defect at least `T` from `Ω`.
    void check_invariant(uint32_t gi, FallbackInstance& F) {
#ifndef NDEBUG
        const StripGraph& g = corpus->g;
        for (uint32_t comp = 0; comp < n_comps; comp++) {
            if (comp_claimed[comp] != (int32_t)gi)
                continue;
            for (uint32_t at = comp_begin[comp]; at < comp_begin[comp + 1]; at++) {
                uint32_t v = comp_members[at];
                if (is_g[v])
                    continue;
                uint32_t det = (uint32_t)h->h_to_det[v];
                bool bad = false;
                dij2.run(g, &det, 1, T_int - 1, [](uint32_t) { return true; }, [&](uint32_t x, int64_t) {
                    if (F.in_omega[x])
                        bad = true;
                });
                if (bad)
                    dbg_invariant_violations++;
            }
        }
#else
        (void)gi;
        (void)F;
#endif
    }

    /// One group's phase 2, from its lowest cube's phase-1 state. Returns false on a conflict.
    bool run_group(uint32_t gi, PhaseOne& p1, PhaseTwoRow& p2) {
        Group& grp = groups[gi];
        FallbackInstance& F = *corpus->fallback[grp.c0];
        const StripGraph& g = corpus->g;
        pm::Mwpm& mwpm = F.mwpm;
        current_group = (int32_t)gi;
        p2.cubes.clear();
        for (size_t i = 0; i < grp.cubes.size(); i++)
            p2.cubes += (i ? ";" : "") + std::to_string(grp.cubes[i]);
        g_defects.clear();
        for (uint32_t det : p1_defects) {
            g_defects.push_back(det);
            is_g[det_to_h[det]] = 1;
        }
        std::fill(cube_seen.begin(), cube_seen.end(), 0);
        for (size_t cc = 0; cc < N; cc++)
            if (corpus->regions[grp.c0].s_cubes[cc])
                cube_seen[cc] = 1;
        p2.omega_initial = (uint32_t)F.omega_size;
        uint32_t prev_job = p1.job;
        bool ok = true;

        auto group_p = [&](uint32_t x) {
            for (uint32_t c : grp.cubes)
                if (corpus->regions[c].in_p[x])
                    return true;
            return false;
        };
        // Adds nodes around `seeds` to Ω; returns false on a conflict.
        std::vector<uint32_t> added;
        auto extend = [&](const std::vector<uint32_t>& seeds, uint32_t& handoff_nodes) -> bool {
            extend_territory(F, seeds.data(), seeds.size(), added, new_outer);
            uint32_t outside = 0;
            for (uint32_t v : added) {
                if (node_conflicts(gi, v))
                    return false;
                cube_seen[g.cube[v]] = 1;
                if (!group_p(v))
                    outside++;
            }
            if (outside > 0) {
                p2.n_handoff++;
                handoff_nodes += outside;
            }
            return true;
        };

        // ---- p2_release
        std::vector<uint32_t> deps = {prev_job};
        uint64_t moved = 0;
        uint64_t r0 = tick();
        fus::align_clock_for_new_growth(mwpm);
        // Release (§5.2, G side): a region tight against an obstacle shell or the outer virtual
        // boundary becomes a root again with its radius kept, and every such slot goes back to what
        // the final problem has there — the real code boundary, or nothing.
        std::vector<pm::GraphFillRegion*> to_release;
        auto scan_virtual = [&](uint32_t x) {
            if (F.slot_kind[x] != SLOT_OBST && F.slot_kind[x] != SLOT_OUTER)
                return;
            pm::DetectorNode& nd = F.node(x);
            pm::GraphFillRegion* top = nd.region_that_arrived_top;
            if (matched_to_boundary(top) && nd.local_radius().get_distance_at_time(mwpm.flooder.queue.cur_time) ==
                                                (pm::cumulative_time_int)nd.neighbor_weights[0])
                to_release.push_back(top);
            F.restore_real_slot(x);
            // Only an occupied node's events depend on its boundary slot.
            if (top != nullptr && !top->radius.is_shrinking())
                mwpm.flooder.reschedule_events_at_detector_node(nd);
        };
        for (uint32_t x : F.virtual_slots)
            scan_virtual(x);
        for (uint32_t x : corpus->regions[grp.c0].outer_nodes)
            scan_virtual(x);
        for (uint32_t x : promoted_outer)
            scan_virtual(x);
        for (uint32_t x : F.touched)
            scan_virtual(x);
        // Deduplicated in discovery order, not by address: the order regions start growing in decides
        // the queue's tie-breaks, so an address order would make the trajectory depend on the heap.
        {
            size_t write = 0;
            for (size_t i = 0; i < to_release.size(); i++)
                if (std::find(to_release.begin(), to_release.begin() + write, to_release[i]) ==
                    to_release.begin() + write)
                    to_release[write++] = to_release[i];
            to_release.resize(write);
        }
        for (pm::GraphFillRegion* region : to_release) {
            region->match.clear();
            region->alt_tree_node = mwpm.node_arena.alloc_constructed(region);
            mwpm.flooder.set_region_growing(*region);
        }
        row.n_released_g += (uint32_t)to_release.size();
        // Claims: every component with a defect in Ω, and the group's unfinished components.
        for (uint32_t v = 0; v < n && ok; v++)
            if (F.in_omega[h->h_to_det[v]])
                ok = claim(gi, comp_of[v], deps, moved, p2);
        for (uint32_t comp : grp.unfinished_comps)
            if (ok)
                ok = claim(gi, comp, deps, moved, p2);
        // `(t1)` owed by phase 1: outer-layer nodes a region already owns.
        pending.clear();
        auto owed = [&](uint32_t x) {
            if (F.in_omega[x] && F.outer[x] && F.node(x).region_that_arrived != nullptr)
                pending.push_back(x);
        };
        for (uint32_t x : corpus->regions[grp.c0].outer_nodes)
            owed(x);
        for (uint32_t x : F.touched)
            owed(x);
        std::sort(pending.begin(), pending.end());
        pending.erase(std::unique(pending.begin(), pending.end()), pending.end());
        uint64_t r1 = tick();
        uint32_t jr = add_job(Kind::P2_RELEASE, CoreClass::FB, (int32_t)grp.c0, deps);
        jobs[jr].measured = r1 - r0;
        jobs[jr].n_items = (uint32_t)to_release.size();
        jobs[jr].link_items = moved;
        row.n_moved_defects += (uint32_t)moved;
        p2.job_ids.push_back(jr);
        prev_job = jr;
        if (!ok)
            return false;

        mwpm.flooder.outer_layer = F.outer.data();
        mwpm.flooder.outer_layer_hit = nullptr;
        bool first_closure = true;
        bool need_closure = true;
        std::vector<uint32_t> capped;
        uint32_t handoff_nodes = 0;

        while (true) {
            // ---- p2_extend
            if (!pending.empty()) {
                uint64_t e0 = tick();
                new_outer.clear();
                uint32_t before = (uint32_t)F.omega_size;
                uint32_t hn = 0;
                ok = extend(pending, hn);
                if (ok)
                    F.apply_reschedules();
                uint64_t e1 = tick();
                uint32_t je = add_job(Kind::P2_EXTEND, CoreClass::FB, (int32_t)grp.c0, {prev_job});
                jobs[je].measured = e1 - e0;
                jobs[je].n_items = (uint32_t)F.omega_size - before;
                jobs[je].link_items = hn;
                handoff_nodes += hn;
                p2.job_ids.push_back(je);
                prev_job = je;
                pending.clear();
                need_closure = true;
                if (!ok)
                    return false;
            }
            // ---- p2_closure
            if (need_closure) {
                deps.assign(1, prev_job);
                moved = 0;
                uint64_t c0t = tick();
                uint32_t steps = 0;
                // Sources: the static band on the first closure (plus anything promotions or an
                // extension added), and only the new outer layer after that.
                std::vector<uint32_t> found;
                if (first_closure) {
                    const CubeRegion& r = corpus->regions[grp.c0];
                    scan_neighbour_defects(grp.c0, r, found,
                                           [&](uint32_t det) { return r.in_band[det] && !F.in_omega[det]; });
                    for (uint32_t x : promoted_outer)
                        new_outer.push_back(x);
                }
                while (ok) {
                    if (!new_outer.empty()) {
                        std::vector<uint32_t> sources;
                        sources.swap(new_outer);
                        std::sort(sources.begin(), sources.end());
                        sources.erase(std::unique(sources.begin(), sources.end()), sources.end());
                        dij2.run(g, sources.data(), sources.size(), T_int - 1,
                                 [&](uint32_t v) { return !F.in_omega[v]; }, [&](uint32_t v, int64_t) {
                                     if (!F.in_omega[v] && det_to_h[v] != NONE)
                                         found.push_back(v);
                                 });
                    }
                    // Defects that landed inside Ω with the last extension.
                    for (uint32_t v : added)
                        if (det_to_h[v] != NONE && !is_g[det_to_h[v]])
                            ok = ok && claim(gi, comp_of[det_to_h[v]], deps, moved, p2);
                    added.clear();
                    std::sort(found.begin(), found.end());
                    found.erase(std::unique(found.begin(), found.end()), found.end());
                    std::vector<uint32_t> seeds;
                    for (uint32_t x : found) {
                        if (F.in_omega[x])
                            continue;
                        ok = ok && claim(gi, comp_of[det_to_h[x]], deps, moved, p2);
                        seeds.push_back(x);
                    }
                    found.clear();
                    if (!ok || seeds.empty())
                        break;
                    steps += (uint32_t)seeds.size();
                    uint32_t hn = 0;
                    ok = extend(seeds, hn);
                    handoff_nodes += hn;
                    moved += hn;
                }
                first_closure = false;
                need_closure = false;
                if (ok) {
                    F.apply_reschedules();
                    // Injection: every claimed member now inside Ω that is not yet a G defect.
                    for (uint32_t comp = 0; comp < n_comps; comp++) {
                        if (comp_claimed[comp] != (int32_t)gi)
                            continue;
                        bool changed = false;
                        for (uint32_t at = comp_begin[comp]; at < comp_begin[comp + 1]; at++) {
                            uint32_t v = comp_members[at];
                            uint32_t det = (uint32_t)h->h_to_det[v];
                            if (is_g[v] || !F.in_omega[det])
                                continue;
#ifndef NDEBUG
                            if (F.node(det).region_that_arrived != nullptr)
                                dbg_injection_violations++;
#endif
                            if (F.node(det).region_that_arrived != nullptr)
                                throw std::runtime_error("phase 2: a defect entered G territory on a covered node");
                            is_g[v] = 1;
                            g_defects.push_back(det);
                            inject_uncapped(mwpm, det);
                            if (F.outer[det])
                                pending.push_back(det);
                            changed = true;
                        }
                        if (changed && !comp_needs_h[comp]) {
                            // Its H-only remainder was solved without these members; redo it.
                            h_only_members(comp, scratch);
                            teardown.run(hmwpm(), scratch.data(), scratch.size());
                            comp_needs_h[comp] = 1;
                        }
                    }
                    mwpm.flooder.outer_layer_hit = nullptr;
                    check_invariant(gi, F);
                }
                uint64_t c1t = tick();
                uint32_t jc = add_job(Kind::P2_CLOSURE, CoreClass::FB, (int32_t)grp.c0, deps);
                jobs[jc].measured = c1t - c0t;
                jobs[jc].n_items = steps;
                jobs[jc].link_items = moved;
                row.n_moved_defects += (uint32_t)moved;
                p2.n_closure_steps += steps;
                p2.job_ids.push_back(jc);
                prev_job = jc;
                if (!ok)
                    return false;
                if (!pending.empty())
                    continue;
            }
            // ---- p2_resume
            uint64_t s0 = tick();
            uint8_t status = ST_COMPLETE;
            uint64_t events = 0;
            uint32_t hit = run_g_until_trigger(mwpm, events);
            p2.resume_events += events;
            if (hit != NONE) {
                status = ST_T1;
                pending.assign(1, hit);
                p2.n_t1++;
            } else {
                for (uint32_t comp = 0; comp < n_comps; comp++) {
                    if (comp_claimed[comp] != (int32_t)gi || !comp_needs_h[comp])
                        continue;
                    if (!solve_h_only(comp, capped)) {
                        status = ST_T2;
                        pending = capped;
                        p2.n_t2++;
                        break;
                    }
                }
            }
            uint64_t s1 = tick();
            uint32_t js = add_job(Kind::P2_RESUME, CoreClass::FB, (int32_t)grp.c0, {prev_job});
            jobs[js].measured = s1 - s0;
            jobs[js].status = status;
            jobs[js].n_items = (uint32_t)events;
            p2.job_ids.push_back(js);
            prev_job = js;
            if (status == ST_COMPLETE)
                break;
        }
        if (!mwpm.node_arena.live.empty())
            dbg_g_unsettled++;

        // ---- the group's extraction: G territory, then the H-only members of its components.
        p2.max_y_over_t = 0;
        {
            pm::cumulative_time_int now = hmwpm().flooder.queue.cur_time;
            for (uint32_t comp = 0; comp < n_comps; comp++) {
                if (comp_claimed[comp] != (int32_t)gi)
                    continue;
                for (uint32_t at = comp_begin[comp]; at < comp_begin[comp + 1]; at++) {
                    uint32_t v = comp_members[at];
                    if (is_g[v])
                        continue;
                    double y = (double)hmwpm().flooder.graph.nodes[v].local_radius().get_distance_at_time(now) /
                               (double)T_int;
                    p2.max_y_over_t = std::max(p2.max_y_over_t, y);
                }
            }
        }
        uint64_t x0 = tick();
        pm::MatchingResult part = fus::extract_root(mwpm, g_defects.data(), g_defects.size());
        for (uint32_t comp = 0; comp < n_comps; comp++) {
            if (comp_claimed[comp] != (int32_t)gi)
                continue;
            h_only_members(comp, scratch);
            part += fus::extract_root(hmwpm(), scratch.data(), scratch.size());
        }
        uint64_t x1 = tick();
        result += part;
        uint32_t jx = add_job(Kind::EXTRACT, CoreClass::FB, (int32_t)grp.c0, {prev_job});
        jobs[jx].measured = x1 - x0;
        jobs[jx].n_items = (uint32_t)g_defects.size();
        combine_deps.push_back(jx);
        p2.omega_final = (uint32_t)F.omega_size;
        uint32_t cubes = 0;
        for (size_t cc = 0; cc < N; cc++)
            cubes += cube_seen[cc];
        p2.omega_cubes = cubes;
        row.handoff_nodes += handoff_nodes;
        reset_for_next_shot(mwpm);
        for (uint32_t det : g_defects)
            is_g[det_to_h[det]] = 0;
        return true;
    }

    // -----------------------------------------------------------------------------------------
    // One shot, start to finish.
    // -----------------------------------------------------------------------------------------

    void run_shot(const std::vector<uint64_t>& shot) {
        Corpus& C = *corpus;
        row = ShotRow();
        result = pm::MatchingResult();
        conflict = false;
        jobs.clear();
        dep_pool.clear();
        bind_pool.clear();
        moved_pool.clear();
        child_pool.clear();
        interior_pool.clear();
        combine_deps.clear();
        groups.clear();
        p2_rows.clear();
        teardown.stray_events = 0;

        prepare_input(shot);
        mark_all_unbuilt(std::max(n, last_n));
        last_n = n;
        // `is_g` is per H index and a conflict abandons a group before its own clean-up runs, so it
        // is cleared here, for every shot, rather than trusted to have been cleared at the end of
        // the last one.
        std::fill(is_g.begin(), is_g.begin() + n, 0);

        row.n_def = n;
        row.n_edges = (uint32_t)h->edges.size();
        row.n_def_cube_min = UINT32_MAX;
        for (size_t c = 0; c < N; c++) {
            uint32_t nc = hoff[c + 1] - hoff[c];
            row.n_def_cube_min = std::min(row.n_def_cube_min, nc);
            row.n_def_cube_max = std::max(row.n_def_cube_max, nc);
        }
        for (size_t s = 0; s + 1 < N; s++)
            row.n_cross_edges += (uint32_t)seam_edges[s].size();
        cube_truncated.assign(N, 0);
        pbase.assign(N, 0);
        pcount.assign(N, 0);
        tbase.assign(N, 0);
        tcount.assign(N, 0);
        n_pieces_total = 0;
        n_tree_nodes = 0;
        job_mgr_scatter.assign(N, NONE);
        job_mgr_uf.assign(N, NONE);
        job_fb_phase1.assign(N, NONE);

        // ---- the per-cube managers
        for (uint32_t c = 0; c < N; c++) {
            uint32_t ja = add_job(Kind::ACCEL, CoreClass::ACC, (int32_t)c, {});
            jobs[ja].measured = (uint64_t)std::llround(accel_ns);
            uint64_t t0 = tick();
            mgr_uf(c);
            uint64_t t1 = tick();
            mgr_tree(c);
            uint64_t t2 = tick();
            mgr_scatter(c);
            uint64_t t3 = tick();
            uint32_t ju = add_job(Kind::MGR_UF, CoreClass::MGR, (int32_t)c, {ja});
            jobs[ju].measured = t1 - t0;
            jobs[ju].n_items = hoff[c + 1] - hoff[c];
            uint32_t jt = add_job(Kind::MGR_TREE, CoreClass::MGR, (int32_t)c, {ju});
            jobs[jt].measured = t2 - t1;
            jobs[jt].n_items = pcount[c];
            uint32_t js = add_job(Kind::MGR_SCATTER, CoreClass::MGR, (int32_t)c, {jt});
            jobs[js].measured = t3 - t2;
            jobs[js].n_items = (uint32_t)cube_edges[c].size();
            job_mgr_uf[c] = ju;
            job_mgr_scatter[c] = js;
            row.n_pieces += pcount[c];
            row.n_fusions += tcount[c] - pcount[c];
        }
        bind_shot();

        // ---- the global seam check
        {
            uint64_t t0 = tick();
            g_seam_check();
            uint64_t t1 = tick();
            job_g_seam_check = add_job(Kind::G_SEAM_CHECK, CoreClass::GM, -1, job_mgr_uf);
            jobs[job_g_seam_check].measured = t1 - t0;
            jobs[job_g_seam_check].n_items = (uint32_t)seam_recs.size();
        }

        // ---- leaves
        uint8_t first_trunc = 0;
        for (uint32_t c = 0; c < N; c++) {
            for (uint32_t i = 0; i < pcount[c]; i++) {
                uint32_t p = pbase[c] + i;
                uint32_t t = tbase[c] + i;
                uint32_t size = piece_end[p] - piece_begin[p];
                uint64_t t0 = tick();
                fus::build_piece(corpus->fusion, build_scratch, build_input, p, piece_members.data() + piece_begin[p],
                                 size, cube_piece_edges[c].data() + piece_edge_begin[p],
                                 piece_edge_end[p] - piece_edge_begin[p]);
                TimelineStatus status =
                    fus::solve_piece(hmwpm(), T_int, piece_members.data() + piece_begin[p], size);
                if (status == TimelineStatus::TRUNCATED)
                    evict(t);
                uint64_t t1 = tick();
                uint32_t j = add_job(Kind::LEAF, CoreClass::SOLVER, (int32_t)c, {job_mgr_scatter[c]});
                jobs[j].measured = t1 - t0;
                jobs[j].priority = size;
                jobs[j].n_items = size;
                jobs[j].dispatcher = (int32_t)c;
                jobs[j].comp = (int32_t)comp_of[piece_members[piece_begin[p]]];
                jobs[j].status = status == TimelineStatus::TRUNCATED ? ST_TRUNCATED : ST_COMPLETE;
                job_of_node[t] = j;
                if (status == TimelineStatus::TRUNCATED && first_trunc == 0)
                    first_trunc = 1;
                if (status == TimelineStatus::COMPLETE)
                    check_feasibility_after(t);
            }
        }
        // ---- intra-cube fusions, children before parents
        for (uint32_t c = 0; c < N; c++) {
            for (uint32_t t = tbase[c] + pcount[c]; t < tbase[c] + tcount[c]; t++) {
                bool dead_child = false;
                std::vector<uint32_t> deps;
                for (uint32_t at = child_begin[t]; at < child_end[t]; at++) {
                    uint32_t ch = child_pool[at];
                    if (tree_dead[ch] || job_of_node[ch] == NONE)
                        dead_child = true;
                    deps.push_back(job_of_node[ch]);
                }
                if (dead_child) {
                    tree_dead[t] = 1;
                    job_of_node[t] = NONE;
                    scratch.clear();
                    collect_defects(t, scratch);
                    comp_unfinished[comp_of[scratch[0]]] = 1;
                    continue;
                }
                uint64_t t0 = tick();
                TimelineStatus status = fuse(t);
                if (status == TimelineStatus::TRUNCATED)
                    evict(t);
                uint64_t t1 = tick();
                uint32_t j = add_job(Kind::FUSE, CoreClass::SOLVER, (int32_t)c, deps);
                set_bind(j, deps);
                jobs[j].measured = t1 - t0;
                jobs[j].priority = interior_end[t] - interior_begin[t];
                jobs[j].n_items = interior_end[t] - interior_begin[t];
                jobs[j].dispatcher = (int32_t)c;
                scratch.clear();
                collect_defects(t, scratch);
                jobs[j].comp = (int32_t)comp_of[scratch[0]];
                jobs[j].status = status == TimelineStatus::TRUNCATED ? ST_TRUNCATED : ST_COMPLETE;
                job_of_node[t] = j;
                if (status == TimelineStatus::TRUNCATED && first_trunc == 0)
                    first_trunc = 2;
                if (status == TimelineStatus::COMPLETE)
                    check_feasibility_after(t);
            }
        }
        // ---- seam fusions (§4.3): one layer, spanning chains in seam order
        {
            uint32_t prev_in_chain = NONE;
            bool chain_broken = false;
            for (const SeamJob& sj : seam_jobs) {
                if (sj.first_in_chain) {
                    prev_in_chain = NONE;
                    chain_broken = false;
                }
                uint32_t comp = sj.comp;
                bool claimed_away = false;
                for (uint32_t cc = comp_cube_lo[comp]; cc <= comp_cube_hi[comp]; cc++)
                    if (cube_truncated[cc])
                        claimed_away = true;
                // Children: the current tops of every endpoint of this seam's edges.
                scratch2.clear();
                for (uint32_t at = sj.begin; at < sj.end; at++) {
                    const BallGraphEdge& edge = h->edges[seam_job_edges[at]];
                    scratch2.push_back(find_tree(tree_of_piece(piece_of[edge.i])));
                    scratch2.push_back(find_tree(tree_of_piece(piece_of[edge.j])));
                }
                std::sort(scratch2.begin(), scratch2.end());
                scratch2.erase(std::unique(scratch2.begin(), scratch2.end()), scratch2.end());
                bool dead_child = false;
                for (uint32_t ch : scratch2)
                    if (tree_dead[ch] || job_of_node[ch] == NONE)
                        dead_child = true;
                if (claimed_away || chain_broken || dead_child) {
                    chain_broken = true;
                    comp_unfinished[comp] = 1;
                    continue;
                }
                uint32_t p = n_tree_nodes++;
                tree_parent[p] = NONE;
                tree_cost[p] = sj.end - sj.begin;
                tree_alive[p] = 1;
                tree_dead[p] = 0;
                tree_cube[p] = sj.seam;
                tree_piece[p] = NONE;
                tuf[p] = p;
                child_begin[p] = (uint32_t)child_pool.size();
                std::vector<uint32_t> deps = {job_g_seam_check};
                std::vector<uint32_t> bind;
                for (uint32_t ch : scratch2) {
                    child_pool.push_back(ch);
                    tree_parent[ch] = p;
                    tree_alive[ch] = 0;
                    deps.push_back(job_of_node[ch]);
                    bind.push_back(job_of_node[ch]);
                }
                child_end[p] = (uint32_t)child_pool.size();
                if (prev_in_chain != NONE)
                    deps.push_back(prev_in_chain);
                interior_begin[p] = (uint32_t)interior_pool.size();
                for (uint32_t at = sj.begin; at < sj.end; at++)
                    interior_pool.push_back(seam_job_edges[at]);
                interior_end[p] = (uint32_t)interior_pool.size();
                uint64_t t0 = tick();
                TimelineStatus status = fuse(p);
                if (status == TimelineStatus::TRUNCATED)
                    evict(p);
                uint64_t t1 = tick();
                uint32_t j = add_job(Kind::SEAM_FUSE, CoreClass::SOLVER, (int32_t)sj.seam, deps);
                set_bind(j, bind);
                jobs[j].seam = (int32_t)sj.seam;
                jobs[j].comp = (int32_t)comp;
                jobs[j].measured = t1 - t0;
                jobs[j].priority = sj.end - sj.begin;
                jobs[j].n_items = sj.end - sj.begin;
                jobs[j].dispatcher = (int32_t)N;
                jobs[j].status = status == TimelineStatus::TRUNCATED ? ST_TRUNCATED : ST_COMPLETE;
                // Per-cube defect counts, for the link charge against whichever side hosts it.
                scratch.clear();
                collect_defects(p, scratch);
                jobs[j].moved_begin = (uint32_t)moved_pool.size();
                for (uint32_t cc = comp_cube_lo[comp]; cc <= comp_cube_hi[comp]; cc++) {
                    uint32_t cnt = 0;
                    for (uint32_t v : scratch)
                        cnt += cube_of_h[v] == cc ? 1 : 0;
                    if (cnt)
                        moved_pool.push_back({cc, cnt});
                }
                jobs[j].moved_end = (uint32_t)moved_pool.size();
                jobs[j].moved_total = (uint32_t)scratch.size();
                job_of_node[p] = j;
                row.n_seam_fusions++;
                prev_in_chain = j;
                if (status == TimelineStatus::TRUNCATED) {
                    chain_broken = true;
                    if (first_trunc == 0)
                        first_trunc = 3;
                } else {
                    check_feasibility_after(p);
                }
            }
        }
        row.trunc_where = first_trunc;
        // Job lists per component, for the claims.
        for (uint32_t j = 0; j < jobs.size(); j++)
            if (jobs[j].comp >= 0)
                comp_jobs[jobs[j].comp].push_back(j);
        for (uint32_t c = 0; c < N; c++)
            row.n_cubes_truncated += cube_truncated[c];
        for (uint32_t comp = 0; comp < n_comps; comp++)
            row.n_trunc_components += comp_truncated[comp];

        // ---- groups (deviation 5)
        form_groups();
        row.n_groups = (uint32_t)groups.size();

        // ---- B1, B3
        uint64_t b1 = 0;
        pm::total_weight_int ref_weight = 0;
        pm::obs_int ref_obs = 0;
        {
            std::vector<uint8_t> obs(std::max<size_t>(1, C.g.num_observables), 0);
            uint64_t t0 = tick();
            pm::decode_detection_events(C.stock, shot, obs.data(), ref_weight, false);
            uint64_t t1 = tick();
            b1 = t1 - t0;
            for (size_t o = 0; o < C.g.num_observables; o++)
                ref_obs |= (pm::obs_int)(obs[o] & 1) << o;
        }
        uint64_t b3 = 0;
        for (uint32_t c = 0; c < N; c++) {
            CubePm& cp = *C.cube_pm[c];
            cp.dets.clear();
            for (uint32_t v = hoff[c]; v < hoff[c + 1]; v++)
                cp.dets.push_back(h->h_to_det[v] - cp.lo);
            pm::total_weight_int wgt = 0;
            uint64_t t0 = tick();
            pm::decode_detection_events(cp.mwpm, cp.dets, cp.obs.data(), wgt, false);
            uint64_t t1 = tick();
            b3 = std::max(b3, t1 - t0);
        }

        // ---- phase 1 on every cube that does not head a group (its state is not needed)
        std::vector<uint8_t> heads(N, 0);
        for (const Group& grp : groups)
            heads[grp.c0] = 1;
        fb_phase1_ticks.assign(N, 0);
        for (uint32_t c = 0; c < N; c++) {
            if (heads[c])
                continue;
            PhaseOne p1 = run_phase1(c);
            job_fb_phase1[c] = p1.job;
            fb_phase1_ticks[c] = p1.ticks;
            row.n_obstacles += p1.obstacles;
            row.n_promotions += p1.promotions;
            discard_fallback(*C.fallback[c], p1_defects);
        }
        // ---- phase 1 then phase 2, per group
        for (uint32_t gi = 0; gi < groups.size() && !conflict; gi++) {
            PhaseOne p1 = run_phase1(groups[gi].c0);
            job_fb_phase1[groups[gi].c0] = p1.job;
            fb_phase1_ticks[groups[gi].c0] = p1.ticks;
            row.n_obstacles += p1.obstacles;
            row.n_promotions += p1.promotions;
            PhaseTwoRow p2;
            bool ok = run_group(gi, p1, p2);
            p2_rows.push_back(p2);
            if (!ok)
                conflict = true;
        }
        if (conflict) {
            for (const Group& grp : groups) {
                corpus->fallback[grp.c0]->mwpm.reset();
            }
        }
        for (const Group& grp : groups)
            corpus->fallback[grp.c0]->undo();

        // ---- extraction of every unclaimed root, then the combine
        if (!conflict) {
            for (uint32_t comp = 0; comp < n_comps; comp++) {
                if (comp_claimed[comp] >= 0)
                    continue;
                uint32_t t = tree_of_piece(piece_of[comp_members[comp_begin[comp]]]);
                while (tree_parent[t] != NONE)
                    t = tree_parent[t];
                uint32_t root_job = job_of_node[t];
                uint64_t t0 = tick();
                scratch.clear();
                collect_defects(t, scratch);
                result += fus::extract_root(hmwpm(), scratch.data(), scratch.size());
                uint64_t t1 = tick();
                uint32_t j = add_job(Kind::EXTRACT, CoreClass::SOLVER, (int32_t)tree_cube[t], {root_job});
                set_bind(j, {root_job});
                jobs[j].measured = t1 - t0;
                jobs[j].comp = (int32_t)comp;
                jobs[j].n_items = (uint32_t)scratch.size();
                combine_deps.push_back(j);
            }
        }
        {
            uint64_t t0 = tick();
            pm::MatchingResult combined = result;
            uint64_t t1 = tick();
            uint32_t j = add_job(Kind::G_COMBINE, CoreClass::GM, -1, combine_deps);
            jobs[j].measured = t1 - t0;
            jobs[j].n_items = (uint32_t)combine_deps.size();
            result = combined;
        }

        // ---- teardown of the H instance, outside every region
        if (conflict || !hmwpm().node_arena.live.empty() || !hmwpm().flooder.region_arena.live.empty())
            abandon_shot(hmwpm());
        else
            reset_for_next_shot(hmwpm());
        for (uint32_t v = 0; v < n; v++)
            det_to_h[h->h_to_det[v]] = NONE;

        // ---- the schedule
        schedule();

        // ---- the row
        uint32_t combine_job = NONE;
        for (uint32_t j = 0; j < jobs.size(); j++)
            if (jobs[j].kind == Kind::G_COMBINE)
                combine_job = j;
        row.t_sys_ns = jobs[combine_job].end;
        walk_critical_path(combine_job);
        for (uint32_t c = 0; c < N; c++)
            row.t_fb_max_ns = std::max(row.t_fb_max_ns, fb_phase1_ticks[c]);
        uint64_t p2_start = UINT64_MAX;
        uint64_t p2_end = 0;
        for (const Job& job : jobs) {
            if (job.kind == Kind::P2_RELEASE || job.kind == Kind::P2_CLOSURE || job.kind == Kind::P2_RESUME ||
                job.kind == Kind::P2_EXTEND) {
                if (job.kind == Kind::P2_RELEASE)
                    p2_start = std::min(p2_start, job.start);
                p2_end = std::max(p2_end, job.end);
            }
            if (job.kind == Kind::SEAM_FUSE)
                row.n_moved_defects += job.moved_defects;
        }
        row.t_p2_ns = p2_start == UINT64_MAX ? 0 : p2_end - p2_start;
        bool truncated = row.n_cubes_truncated > 0 || !groups.empty();
        row.t_spec_ns = truncated ? 0 : row.t_sys_ns;
        row.t_strip_pm_ns = b1;
        row.t_global_escalate_ns = truncated ? b1 : row.t_spec_ns;
        row.t_cube_pm_max_ns = b3;
        row.ref_weight = ref_weight;
        row.ref_obs = ref_obs;
        if (conflict) {
            row.p2_conflict = 1;
            row.t_sys_ns = b1;
            row.weight = ref_weight;
            row.obs_xor = ref_obs;
        } else {
            row.weight = result.weight;
            row.obs_xor = result.obs_mask;
        }
        for (const PhaseTwoRow& p2 : p2_rows) {
            row.omega_nodes += p2.omega_final;
            row.omega_cubes = std::max(row.omega_cubes, p2.omega_cubes);
            row.n_t1 += p2.n_t1;
            row.n_t2 += p2.n_t2;
            row.n_closure_steps += p2.n_closure_steps;
            row.n_reset_components += p2.n_reset_components;
            row.n_reset_defects += p2.n_reset_defects;
            row.n_handoff += p2.n_handoff;
        }
    }

    /// Deviation 5: adjacent truncated cubes, and cubes sharing an unfinished component, form one
    /// group headed by its lowest cube.
    void form_groups() {
        std::vector<uint32_t> cp(N);
        for (uint32_t c = 0; c < N; c++)
            cp[c] = c;
        std::function<uint32_t(uint32_t)> f = [&](uint32_t x) {
            while (cp[x] != x)
                x = cp[x] = cp[cp[x]];
            return x;
        };
        auto unite = [&](uint32_t a, uint32_t b) {
            a = f(a);
            b = f(b);
            if (a != b)
                cp[std::max(a, b)] = std::min(a, b);
        };
        for (uint32_t c = 0; c + 1 < N; c++)
            if (cube_truncated[c] && cube_truncated[c + 1])
                unite(c, c + 1);
        for (uint32_t comp = 0; comp < n_comps; comp++) {
            if (!comp_unfinished[comp])
                continue;
            uint32_t first = NONE;
            for (uint32_t cc = comp_cube_lo[comp]; cc <= comp_cube_hi[comp]; cc++) {
                if (!cube_truncated[cc])
                    continue;
                if (first == NONE)
                    first = cc;
                else
                    unite(first, cc);
            }
            if (first == NONE) {
                std::ostringstream why;
                why << "an unfinished component touches no truncated cube: comp " << comp << " cubes ["
                    << comp_cube_lo[comp] << ", " << comp_cube_hi[comp] << "] size "
                    << comp_begin[comp + 1] - comp_begin[comp] << " truncated_flag " << (int)comp_truncated[comp]
                    << " cube_truncated=";
                for (uint32_t c = 0; c < N; c++)
                    why << (int)cube_truncated[c];
                throw std::logic_error(why.str());
            }
        }
        std::vector<int32_t> group_of(N, -1);
        for (uint32_t c = 0; c < N; c++) {
            if (!cube_truncated[c])
                continue;
            uint32_t root = f(c);
            if (group_of[root] < 0) {
                group_of[root] = (int32_t)groups.size();
                Group g;
                g.c0 = c;
                groups.push_back(g);
            }
            groups[group_of[root]].cubes.push_back(c);
            group_of[c] = group_of[root];
        }
        for (uint32_t comp = 0; comp < n_comps; comp++) {
            if (!comp_unfinished[comp])
                continue;
            for (uint32_t cc = comp_cube_lo[comp]; cc <= comp_cube_hi[comp]; cc++) {
                if (cube_truncated[cc]) {
                    groups[group_of[cc]].unfinished_comps.push_back(comp);
                    break;
                }
            }
        }
    }

    // -----------------------------------------------------------------------------------------
    // §8.2. The list scheduler over `N (k + 2) + 1` cores, and the critical path.
    // -----------------------------------------------------------------------------------------

    inline size_t core_base(uint32_t c) const {
        return (size_t)c * (k + 3);
    }
    inline size_t core_acc(uint32_t c) const {
        return core_base(c);
    }
    inline size_t core_mgr(uint32_t c) const {
        return core_base(c) + 1;
    }
    inline size_t core_solver(uint32_t c, size_t i) const {
        return core_base(c) + 2 + i;
    }
    inline size_t core_fb(uint32_t c) const {
        return core_base(c) + 2 + k;
    }
    inline size_t core_gm() const {
        return N * (k + 3);
    }
    std::string core_label(int32_t core) const {
        if (core < 0)
            return "";
        if ((size_t)core == core_gm())
            return "gm";
        size_t c = (size_t)core / (k + 3);
        size_t r = (size_t)core % (k + 3);
        std::string base = "c" + std::to_string(c) + ".";
        if (r == 0)
            return base + "acc";
        if (r == 1)
            return base + "mgr";
        if (r == k + 2)
            return base + "fb";
        return base + "s" + std::to_string(r - 2);
    }

    std::vector<uint64_t> core_free;
    std::vector<uint32_t> core_last;
    std::vector<uint32_t> dependents_off, dependents;
    std::vector<std::vector<uint32_t>> ready_heaps;

    bool heap_less(uint32_t a, uint32_t b) const {
        // `std::push_heap` keeps the max on top, so this orders *later* things lower.
        const Job& x = jobs[a];
        const Job& y = jobs[b];
        if (x.ready != y.ready)
            return x.ready > y.ready;
        if (x.priority != y.priority)
            return x.priority < y.priority;
        return a > b;
    }

    void schedule() {
        size_t num_cores = core_gm() + 1;
        core_free.assign(num_cores, 0);
        core_last.assign(num_cores, NONE);
        size_t real_jobs = jobs.size();
        dependents_off.assign(real_jobs + 1, 0);
        for (uint32_t j = 0; j < real_jobs; j++) {
            jobs[j].pending = jobs[j].dep_end - jobs[j].dep_begin;
            jobs[j].ready = 0;
            jobs[j].ready_from = NONE;
            for (uint32_t at = jobs[j].dep_begin; at < jobs[j].dep_end; at++)
                dependents_off[dep_pool[at] + 1]++;
        }
        for (size_t j = 0; j < real_jobs; j++)
            dependents_off[j + 1] += dependents_off[j];
        dependents.assign(dependents_off[real_jobs], 0);
        {
            std::vector<uint32_t> cur(dependents_off.begin(), dependents_off.end() - 1);
            for (uint32_t j = 0; j < real_jobs; j++)
                for (uint32_t at = jobs[j].dep_begin; at < jobs[j].dep_end; at++)
                    dependents[cur[dep_pool[at]]++] = j;
        }
        ready_heaps.assign(N + 2, {});
        auto heap_of = [&](const Job& job) -> size_t {
            return job.dispatcher < 0 ? N + 1 : (size_t)job.dispatcher;
        };
        auto cmp = [this](uint32_t a, uint32_t b) {
            return heap_less(a, b);
        };
        auto push = [&](uint32_t j) {
            auto& hp = ready_heaps[heap_of(jobs[j])];
            hp.push_back(j);
            std::push_heap(hp.begin(), hp.end(), cmp);
        };
        for (uint32_t j = 0; j < real_jobs; j++)
            if (jobs[j].pending == 0 && !jobs[j].informational)
                push(j);
        while (true) {
            size_t best = SIZE_MAX;
            for (size_t hp = 0; hp < ready_heaps.size(); hp++) {
                if (ready_heaps[hp].empty())
                    continue;
                if (best == SIZE_MAX || heap_less(ready_heaps[best].front(), ready_heaps[hp].front()))
                    best = hp;
            }
            if (best == SIZE_MAX)
                break;
            bool dispatched = best != N + 1;
            uint64_t d0 = dispatched ? tick() : 0;
            auto& hp = ready_heaps[best];
            std::pop_heap(hp.begin(), hp.end(), cmp);
            uint32_t id = hp.back();
            hp.pop_back();
            size_t core = pick_core(id);
            uint64_t d1 = dispatched ? tick() : 0;
            Job& job = jobs[id];
            uint64_t start = job.ready;
            job.pred = job.ready_from;
            job.determined = job.ready_from == NONE ? BY_NOTHING : BY_DEP;
            if (core_free[core] > start) {
                start = core_free[core];
                job.pred = core_last[core];
                job.determined = BY_CORE;
            }
            if (dispatched) {
                size_t dcore = best < N ? core_mgr((uint32_t)best) : core_gm();
                Job disp;
                disp.kind = Kind::DISPATCH;
                disp.cls = best < N ? CoreClass::MGR : CoreClass::GM;
                disp.cube = best < N ? (int32_t)best : -1;
                disp.measured = d1 - d0;
                disp.core = (int32_t)dcore;
                disp.start = job.ready;
                disp.pred = job.ready_from;
                disp.determined = job.ready_from == NONE ? BY_NOTHING : BY_DEP;
                if (core_free[dcore] > disp.start) {
                    disp.start = core_free[dcore];
                    disp.pred = core_last[dcore];
                    disp.determined = BY_CORE;
                }
                disp.end = disp.start + disp.measured;
                jobs.push_back(disp);
                uint32_t did = (uint32_t)jobs.size() - 1;
                core_free[dcore] = jobs[did].end;
                core_last[dcore] = did;
                Job& j2 = jobs[id];
                if (jobs[did].end > start) {
                    start = jobs[did].end;
                    j2.pred = did;
                    j2.determined = BY_DISPATCH;
                }
            }
            Job& j3 = jobs[id];
            uint64_t delay = 0;
            if (j3.kind == Kind::SEAM_FUSE) {
                uint32_t host_cube = (uint32_t)((size_t)core / (k + 3));
                uint32_t in_host = 0;
                for (uint32_t at = j3.moved_begin; at < j3.moved_end; at++)
                    if (moved_pool[at].first == host_cube)
                        in_host = moved_pool[at].second;
                j3.moved_defects = j3.moved_total - in_host;
                delay = (uint64_t)std::llround(link_ns * (double)j3.moved_defects);
            } else if (j3.link_items) {
                delay = (uint64_t)std::llround(link_ns * (double)j3.link_items);
            }
            j3.delay = delay;
            j3.core = (int32_t)core;
            // A seam fusion belongs to whichever side hosts it; the `seam` column keeps the seam.
            if (j3.kind == Kind::SEAM_FUSE)
                j3.cube = (int32_t)(core / (k + 3));
            j3.start = start + delay;
            j3.end = j3.start + j3.measured;
            core_free[core] = j3.end;
            core_last[core] = id;
            for (uint32_t at = dependents_off[id]; at < dependents_off[id + 1]; at++) {
                uint32_t d = dependents[at];
                Job& next = jobs[d];
                if (j3.end > next.ready || next.ready_from == NONE) {
                    if (j3.end >= next.ready) {
                        next.ready = j3.end;
                        next.ready_from = id;
                    }
                }
                if (--next.pending == 0)
                    push(d);
            }
        }
        // `fb_obstacles` sits at the front of its `fb_phase1`.
        for (uint32_t j = 0; j < real_jobs; j++) {
            if (!jobs[j].informational)
                continue;
            for (uint32_t i = 0; i < real_jobs; i++) {
                if (jobs[i].kind == Kind::FB_PHASE1 && jobs[i].cube == jobs[j].cube) {
                    jobs[j].core = jobs[i].core;
                    jobs[j].start = jobs[i].start;
                    jobs[j].end = jobs[i].start + jobs[j].measured;
                }
            }
        }
    }

    size_t pick_core(uint32_t id) {
        const Job& job = jobs[id];
        uint32_t c = job.cube < 0 ? 0 : (uint32_t)job.cube;
        switch (job.cls) {
            case CoreClass::ACC:
                return core_acc(c);
            case CoreClass::MGR:
                return core_mgr(c);
            case CoreClass::FB:
                return core_fb(c);
            case CoreClass::GM:
                return core_gm();
            case CoreClass::SOLVER:
                break;
        }
        if (job.bind_begin != job.bind_end) {
            // The core of the candidate that finished last; ties to the later job.
            uint32_t best = NONE;
            for (uint32_t at = job.bind_begin; at < job.bind_end; at++) {
                uint32_t cand = bind_pool[at];
                if (cand == NONE)
                    continue;
                if (best == NONE || jobs[cand].end > jobs[best].end ||
                    (jobs[cand].end == jobs[best].end && cand > best))
                    best = cand;
            }
            if (best != NONE && jobs[best].core >= 0)
                return (size_t)jobs[best].core;
        }
        size_t core = core_solver(c, 0);
        for (size_t i = 1; i < k; i++)
            if (core_free[core_solver(c, i)] < core_free[core])
                core = core_solver(c, i);
        return core;
    }

    void walk_critical_path(uint32_t from) {
        uint32_t at = from;
        size_t guard = 0;
        while (at != NONE && guard++ <= jobs.size() + 2) {
            const Job& job = jobs[at];
            uint64_t* bucket = &row.cp_idle;
            switch (job.kind) {
                case Kind::ACCEL:
                    bucket = &row.cp_accel;
                    break;
                case Kind::MGR_UF:
                case Kind::MGR_TREE:
                case Kind::MGR_SCATTER:
                    bucket = &row.cp_mgr;
                    break;
                case Kind::DISPATCH:
                    bucket = job.cls == CoreClass::GM ? &row.cp_gm : &row.cp_mgr;
                    break;
                case Kind::LEAF:
                    bucket = &row.cp_leaf;
                    break;
                case Kind::FUSE:
                    bucket = &row.cp_fuse;
                    break;
                case Kind::SEAM_FUSE:
                    bucket = &row.cp_seam;
                    break;
                case Kind::EXTRACT:
                    bucket = &row.cp_extract;
                    break;
                case Kind::G_SEAM_CHECK:
                case Kind::G_COMBINE:
                    bucket = &row.cp_gm;
                    break;
                case Kind::FB_PHASE1:
                case Kind::FB_OBSTACLES:
                    bucket = &row.cp_fb;
                    break;
                case Kind::P2_RELEASE:
                case Kind::P2_CLOSURE:
                case Kind::P2_RESUME:
                case Kind::P2_EXTEND:
                    bucket = &row.cp_p2;
                    break;
            }
            *bucket += job.measured;
            row.cp_link += job.delay;
            uint64_t base = job.start - job.delay;
            uint64_t prev_end = job.pred == NONE ? 0 : jobs[job.pred].end;
            if (base > prev_end)
                row.cp_idle += base - prev_end;
            at = job.pred;
        }
    }
};

// ---------------------------------------------------------------------------------------------
// Output.
// ---------------------------------------------------------------------------------------------

void write_shot_header(std::ofstream& out) {
    out << "shot,n_def,n_def_cube_min,n_def_cube_max,n_edges,n_cross_edges,n_seam_components,n_spanning,max_span,"
           "seam_depth,n_cubes_truncated,n_trunc_components,trunc_where,n_groups,omega_nodes,omega_cubes,n_t1,n_t2,"
           "n_closure_steps,n_reset_components,n_reset_defects,n_obstacles,n_promotions,n_handoff,handoff_nodes,"
           "n_released_g,n_released_h,n_moved_defects,p2_conflict,n_pieces,n_fusions,n_seam_fusions,t_spec_ns,"
           "t_fb_max_ns,t_p2_ns,t_sys_ns,t_strip_pm_ns,t_global_escalate_ns,t_cube_pm_max_ns,cp_accel,cp_mgr,"
           "cp_leaf,cp_fuse,cp_seam,cp_extract,cp_gm,cp_fb,cp_p2,cp_link,cp_idle,weight,ref_weight,obs_xor,"
           "ref_obs\n";
}

void write_shot_row(std::ofstream& out, const ShotRow& r) {
    out << r.shot << ',' << r.n_def << ',' << r.n_def_cube_min << ',' << r.n_def_cube_max << ',' << r.n_edges << ','
        << r.n_cross_edges << ',' << r.n_seam_components << ',' << r.n_spanning << ',' << r.max_span << ','
        << r.seam_depth << ',' << r.n_cubes_truncated << ',' << r.n_trunc_components << ','
        << trunc_where_name(r.trunc_where) << ',' << r.n_groups << ',' << r.omega_nodes << ',' << r.omega_cubes << ','
        << r.n_t1 << ',' << r.n_t2 << ',' << r.n_closure_steps << ',' << r.n_reset_components << ','
        << r.n_reset_defects << ',' << r.n_obstacles << ',' << r.n_promotions << ',' << r.n_handoff << ','
        << r.handoff_nodes << ',' << r.n_released_g << ',' << r.n_released_h << ',' << r.n_moved_defects << ','
        << (int)r.p2_conflict << ',' << r.n_pieces << ',' << r.n_fusions << ',' << r.n_seam_fusions << ','
        << r.t_spec_ns << ',' << r.t_fb_max_ns << ',' << r.t_p2_ns << ',' << r.t_sys_ns << ',' << r.t_strip_pm_ns
        << ',' << r.t_global_escalate_ns << ',' << r.t_cube_pm_max_ns << ',' << r.cp_accel << ',' << r.cp_mgr << ','
        << r.cp_leaf << ',' << r.cp_fuse << ',' << r.cp_seam << ',' << r.cp_extract << ',' << r.cp_gm << ','
        << r.cp_fb << ',' << r.cp_p2 << ',' << r.cp_link << ',' << r.cp_idle << ',' << r.weight << ','
        << r.ref_weight << ',' << r.obs_xor << ',' << r.ref_obs << '\n';
}

void write_jobs(std::ofstream& out, const ShotRunner& runner, uint64_t shot) {
    for (uint32_t j = 0; j < runner.jobs.size(); j++) {
        const Job& job = runner.jobs[j];
        if (job.kind == Kind::DISPATCH)
            continue;
        out << shot << ',' << j << ',' << kind_name(job.kind) << ',' << job.cube << ',' << job.seam << ',' << job.comp
            << ',' << runner.core_label(job.core) << ',';
        for (uint32_t at = job.dep_begin; at < job.dep_end; at++)
            out << (at > job.dep_begin ? ";" : "") << runner.dep_pool[at];
        out << ',' << job.start << ',' << job.end << ',' << job.measured << ',' << job.n_items << ','
            << status_name(job.status) << '\n';
    }
}

void write_phase2(std::ofstream& out, const ShotRunner& runner, uint64_t shot) {
    for (const PhaseTwoRow& p2 : runner.p2_rows) {
        uint64_t start = UINT64_MAX;
        uint64_t end = 0;
        for (uint32_t j : p2.job_ids) {
            start = std::min(start, runner.jobs[j].start);
            end = std::max(end, runner.jobs[j].end);
        }
        if (start == UINT64_MAX)
            start = 0;
        char y[32];
        std::snprintf(y, sizeof(y), "%.6f", p2.max_y_over_t);
        out << shot << ',' << p2.cubes << ',' << start << ',' << end << ',' << p2.omega_initial << ','
            << p2.omega_final << ',' << p2.omega_cubes << ',' << p2.n_t1 << ',' << p2.n_t2 << ','
            << p2.n_closure_steps << ',' << p2.n_reset_components << ',' << p2.n_reset_defects << ','
            << p2.n_handoff << ',' << y << ',' << p2.resume_events << '\n';
    }
}

struct Summary {
    double mean{0}, p50{0}, p99{0}, p99_9{0};
    long long max{0};
    uint64_t count{0};
};

Summary summarise(const Series& s, int which) {
    std::vector<long long> v;
    for (size_t i = 0; i < s.values.size(); i++) {
        if (which == 1 && s.truncated[i])
            continue;
        if (which == 2 && !s.truncated[i])
            continue;
        v.push_back(s.values[i]);
    }
    Summary out;
    out.count = v.size();
    if (v.empty())
        return out;
    long double sum = 0;
    for (long long x : v) {
        sum += x;
        out.max = std::max(out.max, x);
    }
    out.mean = (double)(sum / v.size());
    out.p50 = percentile(v, 0.5);
    out.p99 = percentile(v, 0.99);
    out.p99_9 = percentile(v, 0.999);
    return out;
}

void write_summary(std::ofstream& out, const char* name, const Series& s, bool last) {
    out << "      \"" << name << "\": {";
    const char* subsets[3] = {"all", "nontruncated", "truncated"};
    for (int w = 0; w < 3; w++) {
        Summary sm = summarise(s, w);
        out << "\"" << subsets[w] << "\": {\"count\": " << sm.count << ", \"mean\": " << fmt_g(sm.mean)
            << ", \"p50\": " << fmt_g(sm.p50) << ", \"p99\": " << fmt_g(sm.p99) << ", \"p99_9\": " << fmt_g(sm.p99_9)
            << ", \"max\": " << sm.max << "}" << (w < 2 ? ", " : "");
    }
    out << "}" << (last ? "" : ",") << "\n";
}

void write_hist(std::ofstream& out, const char* name, const std::vector<uint64_t>& hist) {
    out << "      \"" << name << "\": {";
    bool first = true;
    for (size_t i = 0; i < hist.size(); i++) {
        if (!hist[i])
            continue;
        out << (first ? "" : ", ") << "\"" << i << "\": " << hist[i];
        first = false;
    }
    out << "},\n";
}

void write_aggregate(std::ofstream& out, const Options& options, const std::vector<CellAggregate>& cells) {
    static const char* cp_names[11] = {"accel", "mgr",  "leaf", "fuse", "seam", "extract",
                                       "gm",    "fb",   "p2",   "link", "idle"};
    out << "{\n";
    out << "  \"binary\": \"lattice_surgery_profiler\",\n";
    out << "  \"git_hash\": \"" << SPECMATCHING_GIT_HASH << "\",\n";
    out << "  \"stim_version\": \"" << SPECMATCHING_STIM_VERSION << "\",\n";
    out << "  \"model\": \"lattice_surgery_k_core_model\",\n";
    out << "  \"tick_unit\": \"" << hires_timer_name() << ", uncorrected; see run.log timer_overhead_ns\",\n";
    out << "  \"seed\": " << options.seed << ",\n";
    out << "  \"seam_gap\": " << options.seam_gap << ",\n";
    out << "  \"merge_init\": " << (options.merge_init ? "true" : "false") << ",\n";
    out << "  \"accel_ns\": " << fmt_g(options.accel_ns) << ",\n";
    out << "  \"link_ns\": " << fmt_g(options.link_ns) << ",\n";
    out << "  \"sort_edges\": " << (options.sort_edges ? "true" : "false") << ",\n";
    out << "  \"verify\": " << (options.verify ? "true" : "false") << ",\n";
    out << "  \"t_spec_over\": \"non-truncated shots only (t_spec_ns is 0 on a truncated row)\",\n";
    out << "  \"cells\": [\n";
    for (size_t i = 0; i < cells.size(); i++) {
        const CellAggregate& c = cells[i];
        double shots = (double)std::max<uint64_t>(1, c.shots);
        out << "    {\n";
        out << "      \"cell\": \"" << c.label << "\",\n";
        out << "      \"N\": " << c.N << ", \"k\": " << c.k << ", \"d\": " << c.d << ", \"p\": " << fmt_g(c.p)
            << ", \"T\": " << fmt_g(c.T) << ", \"W\": " << fmt_g(c.W) << ", \"orientation\": \""
            << (c.orientation ? "long" : "short") << "\",\n";
        out << "      \"cores\": " << c.N * (c.k + 2) + 1 << ",\n";
        out << "      \"shots\": " << c.shots << ", \"shots_truncated\": " << c.shots_truncated
            << ", \"strip_trunc_rate\": " << fmt_g((double)c.shots_truncated / shots)
            << ", \"cube_trunc_rate\": " << fmt_g((double)c.cubes_truncated / (shots * (double)c.N))
            << ", \"trunc_components_mean\": " << fmt_g((double)c.trunc_components / shots)
            << ", \"first_trunc\": {\"leaf\": " << c.trunc_leaf << ", \"fuse\": " << c.trunc_fuse
            << ", \"seam\": " << c.trunc_seam << "},\n";
        out << "      \"q_reference\": \"compare cube_trunc_rate with an N=1 cell at the same (d, p, T, k): the same"
               " generator's d x d x d memory experiment\",\n";
        write_summary(out, "t_spec_ns", c.t_spec, false);
        write_summary(out, "t_fb_max_ns", c.t_fb_max, false);
        write_summary(out, "t_p2_ns", c.t_p2, false);
        write_summary(out, "t_sys_ns", c.t_sys, false);
        write_summary(out, "t_strip_pm_ns", c.t_strip_pm, false);
        write_summary(out, "t_global_escalate_ns", c.t_global_escalate, false);
        write_summary(out, "t_cube_pm_max_ns", c.t_cube_pm_max, false);
        double sys_mean = summarise(c.t_sys, 0).mean;
        out << "      \"speedup_vs_strip_pm\": "
            << fmt_g(sys_mean > 0 ? summarise(c.t_strip_pm, 0).mean / sys_mean : 0)
            << ", \"speedup_vs_global_escalate\": "
            << fmt_g(sys_mean > 0 ? summarise(c.t_global_escalate, 0).mean / sys_mean : 0) << ",\n";
        out << "      \"critical_path_means\": {";
        for (int b = 0; b < 11; b++)
            out << (b ? ", " : "") << "\"" << cp_names[b] << "\": " << fmt_g((double)c.cp_sum[b] / shots);
        out << "},\n";
        out << "      \"means\": {\"n_seam_components\": " << fmt_g((double)c.seam_component_sum / shots)
            << ", \"n_spanning\": " << fmt_g((double)c.spanning_component_sum / shots)
            << ", \"n_obstacles\": " << fmt_g((double)c.obstacles_sum / shots)
            << ", \"n_promotions\": " << fmt_g((double)c.promotions_sum / shots)
            << ", \"n_moved_defects\": " << fmt_g((double)c.moved_sum / shots)
            << ", \"n_released_g\": " << fmt_g((double)c.released_g_sum / shots)
            << ", \"n_t1\": " << fmt_g((double)c.t1_sum / shots) << ", \"n_t2\": " << fmt_g((double)c.t2_sum / shots)
            << ", \"n_closure_steps\": " << fmt_g((double)c.closure_sum / shots) << "},\n";
        out << "      \"shots_with_spanning_chain\": " << c.n_spanning_shots << ",\n";
        write_hist(out, "seam_depth_hist", c.seam_depth_hist);
        write_hist(out, "omega_cubes_hist", c.omega_cubes_hist);
        write_hist(out, "n_t1_plus_t2_hist", c.triggers_hist);
        write_hist(out, "n_reset_defects_hist", c.reset_defects_hist);
        out << "      \"handoff_events\": " << c.handoff_events << ", \"shots_with_handoff\": " << c.shots_with_handoff
            << ", \"omega_cubes_gt3\": " << c.omega_cubes_gt3 << ", \"p2_conflicts\": " << c.p2_conflicts << ",\n";
        out << "      \"verify_shots\": " << c.verify_shots << ", \"weight_mismatches\": " << c.weight_mismatches
            << ", \"obs_mismatches\": " << c.obs_mismatches
            << ", \"first_weight_mismatch_shot\": " << c.first_weight_mismatch_shot << ",\n";
        out << "      \"determinism_check\": " << c.determinism_check
            << ", \"determinism_mismatches\": " << c.determinism_mismatches << ",\n";
        out << "      \"stray_events_after_teardown\": " << c.stray_events << ",\n";
        out << "      \"debug\": {\"invariant_violations\": " << c.dbg_invariant_violations
            << ", \"cap_violations\": " << c.dbg_cap_violations
            << ", \"injection_violations\": " << c.dbg_injection_violations
            << ", \"feasibility_violations\": " << c.dbg_feasibility_violations
            << ", \"g_unsettled\": " << c.dbg_g_unsettled << "},\n";
        out << "      \"convexity\": [";
        for (size_t cc = 0; cc < c.convexity_pairs.size(); cc++)
            out << (cc ? ", " : "") << "{\"cube\": " << cc << ", \"pairs\": " << c.convexity_pairs[cc]
                << ", \"mismatches\": " << c.convexity_mismatches[cc] << "}";
        out << "]\n";
        out << "    }" << (i + 1 == cells.size() ? "" : ",") << "\n";
    }
    out << "  ]\n}\n";
}

uint64_t row_digest(const ShotRow& r) {
    // Everything a function of the shot alone: no tick, and nothing the schedule decides.
    Digest d;
    d.add(r.n_def);
    d.add(r.n_edges);
    d.add(r.n_cross_edges);
    d.add(r.n_seam_components);
    d.add(r.n_spanning);
    d.add(r.seam_depth);
    d.add(r.n_cubes_truncated);
    d.add(r.n_trunc_components);
    d.add(r.n_groups);
    d.add(r.omega_nodes);
    d.add(r.omega_cubes);
    d.add(r.n_t1);
    d.add(r.n_t2);
    d.add(r.n_closure_steps);
    d.add(r.n_reset_components);
    d.add(r.n_reset_defects);
    d.add(r.n_obstacles);
    d.add(r.n_promotions);
    d.add(r.n_handoff);
    d.add(r.n_released_g);
    d.add(r.n_pieces);
    d.add(r.n_fusions);
    d.add(r.n_seam_fusions);
    d.add(r.weight);
    d.add(r.obs_xor);
    return d.value;
}

// ---------------------------------------------------------------------------------------------
// The driver.
// ---------------------------------------------------------------------------------------------

void build_corpus(const Options& options, size_t N, size_t d, double p, int orientation, double T_max, Corpus& C,
                  std::ofstream& log) {
    LatticeSurgeryParams params;
    params.distance = d;
    params.patches = N;
    params.seam_gap = options.seam_gap;
    params.rounds = d;
    params.noise = p;
    params.merge_init = options.merge_init;
    params.long_orientation = orientation != 0;
    C.lc = generate_lattice_surgery_circuit(params);
    C.dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(C.lc.circuit, false, true, false, 0, false, false);

    pm::Mwpm probe = pm::detector_error_model_to_mwpm(C.dem, NUM_DISTINCT_WEIGHTS, false);
    C.unit = edge_weight_units(probe.flooder.graph);
    C.w_p = std::log((1 - p) / p);

    BallConfig config;
    config.T = T_max * C.unit;
    config.ball.T_max = T_max * C.unit;
    config.ball.R = 2.0 * config.ball.T_max;
    config.collect_component_stats = false;
    config.collect_harvest_diagnostics = false;
    config.collect_structural_counters = false;
    config.verify_component_decomposition = false;
    config.verify_against_g = false;
    config.ball.validate();
    C.decoder = std::make_unique<BallDecoder>(BallDecoder::from_detector_error_model(C.dem, config, NUM_DISTINCT_WEIGHTS));
    if (C.decoder->dem_has_negative_weights)
        throw std::runtime_error("the strip's DEM has a negative-weight edge; this binary assumes none");
    C.g.build(C.decoder->g_mwpm.flooder.graph, C.lc);
    C.stock = pm::detector_error_model_to_mwpm(C.dem, NUM_DISTINCT_WEIGHTS, false);
    C.cube_pm.clear();
    for (uint32_t c = 0; c < N; c++) {
        C.cube_pm.push_back(std::make_unique<CubePm>());
        C.cube_pm.back()->build(C.g, c);
    }
    size_t n_max = C.g.n;
    C.shared = std::make_unique<BallMwpm>();
    C.shared->configure(C.g.num_observables, C.g.normalising_constant);
    {
        BallGraph full;
        full.h_to_det.resize(n_max, 0);
        C.shared->rebuild(C.decoder->tables, full, C.decoder->arena, nullptr);
        C.shared->used_nodes = 0;
        C.shared->grow_events = 0;
    }
    C.fusion = fus::FusionInstance();
    C.fusion.mask.reserve(n_max, 4 * INITIAL_EDGE_CAPACITY);
    C.fallback.clear();
    for (uint32_t c = 0; c < N; c++) {
        C.fallback.push_back(std::make_unique<FallbackInstance>());
        C.fallback.back()->build(C.g, c);
    }
    C.regions.assign(N, CubeRegion());

    const BallStats& bs = C.decoder->tables.stats;
    log << "\ncorpus N=" << N << " d=" << d << " p=" << p << " orientation=" << (orientation ? "long" : "short")
        << " seam_gap=" << options.seam_gap << " merge_init=" << (options.merge_init ? "on" : "off") << "\n";
    log << "  generator_call=" << C.lc.describe() << "\n";
    log << "  strip width=" << C.lc.width << " height=" << C.lc.height << " detectors=" << C.g.n
        << " observables=" << C.g.num_observables << " seam_plaquettes=" << C.lc.seam_plaquettes
        << " removed_first_round_detectors=" << C.lc.removed_first_round_detectors << "\n";
    log << "  cube_rule=detector of plaquette (i, j) belongs to the cube of data column max(i, 0); a routing"
           " column belongs to the cube on its left; detector ids are cube-major\n";
    log << "  cube_detectors=";
    for (uint32_t c = 0; c < N; c++)
        log << (c ? "," : "") << C.g.cube_off[c + 1] - C.g.cube_off[c];
    log << "\n";
    log << "  w_p=ln((1-p)/p)=" << C.w_p << " edge_weight_unit(median lattice edge, DEM units)=" << C.unit
        << " normalising_constant=" << C.g.normalising_constant << "\n";
    log << "  w_min=" << C.g.w_min << " w_max=" << C.g.w_max << " (time units; in DEM units "
        << C.g.w_min / C.g.normalising_constant << ", " << C.g.w_max / C.g.normalising_constant
        << ") distinct_weights=" << C.g.n_distinct << "\n";
    log << "  ball_tables T_max=" << config.ball.T_max << " R=" << config.ball.R << " (DEM units) entries="
        << bs.total_entries << " mean_ball=" << bs.mean_ball_size << " bytes=" << bs.bytes_total
        << " compile_s=" << bs.compile_wall_seconds << " (compiled on the FULL strip G, §1.3)\n";
    log.flush();
}

void run_all(const Options& options, std::ofstream& log, bool interactive) {
    std::vector<CellAggregate> cells;
    size_t total_cells = options.patches.size() * options.distances.size() * options.error_rates.size() *
                         options.orientations.size() * options.horizons.size() * options.slab_widths.size() *
                         options.core_counts.size();
    size_t cell_index = 0;
    double T_max = *std::max_element(options.horizons.begin(), options.horizons.end());
    std::string error;

    auto write_agg = [&]() {
        std::ofstream agg(options.out_dir + "/agg.json");
        write_aggregate(agg, options, cells);
    };

    for (size_t N : options.patches) {
        for (size_t d : options.distances) {
            for (double p : options.error_rates) {
                for (int orientation : options.orientations) {
                    char setup_label[160];
                    std::snprintf(setup_label, sizeof(setup_label), "N=%zu d=%zu p=%g %s", N, d, p,
                                  orientation ? "long" : "short");
                    std::fprintf(stderr, interactive ? "\r%s  building circuit, DEM, ball tables...\x1b[K"
                                                     : "%s  building circuit, DEM, ball tables...\n",
                                 setup_label);
                    Corpus C;
                    build_corpus(options, N, d, p, orientation, T_max, C, log);
                    ShotSampler sampler;
                    sampler.circuit = C.lc.circuit;
                    sampler.num_detectors = C.lc.num_detectors;

                    for (double T : options.horizons) {
                        C.decoder->set_horizon(T * C.unit);
                        C.fusion.attach(C.shared->mwpm, C.g.n, C.decoder->horizon);
                        for (double W : options.slab_widths) {
                            double preload = options.preload_width > 0 ? options.preload_width : std::max(2 * W, W + 2 * T);
                            int64_t W_int = to_time_units(W * C.unit, C.g.normalising_constant);
                            int64_t P_int = to_time_units(preload * C.unit, C.g.normalising_constant);
                            int64_t T_int = C.decoder->horizon;
                            if (P_int < W_int + T_int)
                                throw std::invalid_argument("--preload-width must be at least W + T");
                            Dijkstra a, b;
                            a.configure(C.g.n);
                            b.configure(C.g.n);
                            std::vector<uint64_t> conv_pairs(N), conv_mism(N);
                            bool convex_ok = true;
                            for (uint32_t c = 0; c < N; c++) {
                                compute_cube_region(C.g, c, W_int, T_int, P_int, a, C.regions[c]);
                                check_convexity(C.g, 2 * T_int + W_int, options.convexity_sources, a, b, C.regions[c]);
                                C.fallback[c]->set_baseline(C.regions[c]);
                                conv_pairs[c] = C.regions[c].convexity_pairs;
                                conv_mism[c] = C.regions[c].convexity_mismatches;
                                if (C.regions[c].convexity_mismatches)
                                    convex_ok = false;
                                log << "  T=" << T << " W=" << W << " cube " << c << ": |S|=" << C.regions[c].s_nodes.size()
                                    << " |B|=" << C.regions[c].b_nodes << " |outer|=" << C.regions[c].outer_nodes.size()
                                    << " |P|=" << C.regions[c].p_nodes << " |band|=" << C.regions[c].band_nodes
                                    << " face_slab=" << C.regions[c].face_slab.size()
                                    << " convexity: sources=" << C.regions[c].convexity_sources_used
                                    << " pairs=" << C.regions[c].convexity_pairs
                                    << " mismatches=" << C.regions[c].convexity_mismatches << "\n";
                            }
                            log << "  T=" << T << " T_weight_units=" << T * C.unit << " T_int=" << T_int << " W=" << W
                                << " W_weight_units=" << W * C.unit << " W_int=" << W_int << " preload=" << preload
                                << " preload_weight_units=" << preload * C.unit << " preload_int=" << P_int
                                << " convexity=" << (convex_ok ? "pass" : "FAIL") << "\n";
                            log.flush();
                            if (!convex_ok && !options.allow_nonconvex) {
                                throw std::runtime_error(
                                    "geodesic convexity check failed (§9.2); see run.log, or pass --allow-nonconvex");
                            }

                            for (size_t k : options.core_counts) {
                                cell_index++;
                                std::ostringstream lab;
                                lab << "N=" << N << ",k=" << k << ",d=" << d << ",p=" << fmt_g(p) << ",T=" << fmt_g(T)
                                    << ",W=" << fmt_g(W) << "," << (orientation ? "long" : "short");
                                std::string label = lab.str();
                                std::ostringstream suffix;
                                suffix << "N" << N << "_k" << k << "_d" << d << "_p" << fmt_g(p) << "_T" << fmt_g(T)
                                       << "_W" << fmt_g(W) << "_" << (orientation ? "long" : "short") << ".csv";
                                std::ofstream shots_csv(options.out_dir + "/shots_" + suffix.str());
                                std::ofstream jobs_csv;
                                std::ofstream p2_csv(options.out_dir + "/phase2_" + suffix.str());
                                if (!shots_csv.is_open() || !p2_csv.is_open())
                                    throw std::runtime_error("could not open the cell's output files");
                                std::string header_comment = "# N=" + std::to_string(N) + ",k=" + std::to_string(k) +
                                                             ",d=" + std::to_string(d) + ",p=" + fmt_g(p) +
                                                             ",T=" + fmt_g(T) + ",W=" + fmt_g(W) +
                                                             ",orientation=" + (orientation ? "long" : "short") +
                                                             ",seam_gap=" + std::to_string(options.seam_gap) +
                                                             ",merge_init=" + (options.merge_init ? "on" : "off") +
                                                             ",seed=" + std::to_string(options.seed) +
                                                             ",tick_unit=" + hires_timer_name() + ",uncorrected=1\n";
                                shots_csv << header_comment;
                                write_shot_header(shots_csv);
                                p2_csv << header_comment;
                                p2_csv << "shot,cubes,start_ns,end_ns,omega_nodes_initial,omega_nodes_final,omega_cubes,"
                                          "n_t1,n_t2,n_closure_steps,n_reset_components,n_reset_defects,n_handoff,"
                                          "max_Y_over_T,resume_events\n";
                                if (options.jobs_log != JobsLog::NONE) {
                                    jobs_csv.open(options.out_dir + "/jobs_" + suffix.str());
                                    jobs_csv << header_comment;
                                    jobs_csv << "shot,job_id,kind,cube,seam,component,core,deps,start_ns,end_ns,"
                                                "measured_ns,n_items,status\n";
                                }

                                CellAggregate agg;
                                agg.N = N;
                                agg.k = k;
                                agg.d = d;
                                agg.p = p;
                                agg.T = T;
                                agg.W = W;
                                agg.orientation = orientation;
                                agg.label = label;
                                agg.reset(N + 2);
                                agg.convexity_pairs = conv_pairs;
                                agg.convexity_mismatches = conv_mism;

                                ShotRunner runner;
                                runner.options = &options;
                                runner.corpus = &C;
                                runner.N = N;
                                runner.k = k;
                                runner.T_int = T_int;
                                runner.W_int = W_int;
                                runner.P_int = P_int;
                                runner.sort_by_weight = options.sort_edges;
                                runner.accel_ns = options.accel_ns;
                                runner.link_ns = options.link_ns;
                                runner.configure(C.g.n);

                                log << "  cell " << label << " started=" << stamp() << "\n";
                                log.flush();
                                sampler.rng = std::mt19937_64(options.seed);
                                size_t done = 0;
                                size_t total = options.warmup + options.shots;
                                draw_progress(interactive, label, cell_index, total_cells, 0, total);
                                uint64_t determinism_mismatches = 0;
                                while (done < total) {
                                    size_t batch = std::min(SHOT_BATCH, total - done);
                                    sampler.sample(batch);
                                    for (size_t i = 0; i < batch; i++) {
                                        size_t index = done + i;
                                        bool measured = index >= options.warmup;
                                        runner.profile = measured;
                                        runner.run_shot(sampler.shots[i]);
                                        if (!measured)
                                            continue;
                                        uint64_t shot_index = index - options.warmup;
                                        ShotRow r = runner.row;
                                        r.shot = shot_index;
                                        if (options.check_determinism) {
                                            uint64_t first = row_digest(r);
                                            std::vector<Job> saved_jobs = runner.jobs;
                                            std::vector<uint32_t> saved_deps = runner.dep_pool;
                                            std::vector<PhaseTwoRow> saved_p2 = runner.p2_rows;
                                            runner.run_shot(sampler.shots[i]);
                                            if (row_digest(runner.row) != first)
                                                determinism_mismatches++;
                                            runner.jobs = saved_jobs;
                                            runner.dep_pool = saved_deps;
                                            runner.p2_rows = saved_p2;
                                        }
                                        write_shot_row(shots_csv, r);
                                        bool truncated = r.n_cubes_truncated > 0;
                                        if (options.jobs_log == JobsLog::ALL ||
                                            (options.jobs_log == JobsLog::SAMPLE &&
                                             (shot_index < JOBS_SAMPLE_FIRST || truncated)))
                                            write_jobs(jobs_csv, runner, shot_index);
                                        write_phase2(p2_csv, runner, shot_index);

                                        agg.shots++;
                                        agg.shots_truncated += truncated;
                                        agg.cubes_truncated += r.n_cubes_truncated;
                                        agg.trunc_components += r.n_trunc_components;
                                        agg.trunc_leaf += r.trunc_where == 1;
                                        agg.trunc_fuse += r.trunc_where == 2;
                                        agg.trunc_seam += r.trunc_where == 3;
                                        if (!truncated)
                                            agg.t_spec.add(r.t_spec_ns, false);
                                        agg.t_fb_max.add(r.t_fb_max_ns, truncated);
                                        if (truncated)
                                            agg.t_p2.add(r.t_p2_ns, true);
                                        agg.t_sys.add(r.t_sys_ns, truncated);
                                        agg.t_strip_pm.add(r.t_strip_pm_ns, truncated);
                                        agg.t_global_escalate.add(r.t_global_escalate_ns, truncated);
                                        agg.t_cube_pm_max.add(r.t_cube_pm_max_ns, truncated);
                                        uint64_t cps[11] = {r.cp_accel, r.cp_mgr, r.cp_leaf,    r.cp_fuse,
                                                            r.cp_seam,  r.cp_extract, r.cp_gm, r.cp_fb,
                                                            r.cp_p2,    r.cp_link,    r.cp_idle};
                                        for (int b = 0; b < 11; b++)
                                            agg.cp_sum[b] += cps[b];
                                        if (r.seam_depth < agg.seam_depth_hist.size())
                                            agg.seam_depth_hist[r.seam_depth]++;
                                        if (r.omega_cubes < agg.omega_cubes_hist.size())
                                            agg.omega_cubes_hist[r.omega_cubes]++;
                                        uint32_t trig = r.n_t1 + r.n_t2;
                                        agg.triggers_hist[std::min<size_t>(trig, agg.triggers_hist.size() - 1)]++;
                                        agg.reset_defects_hist[std::min<size_t>(r.n_reset_defects,
                                                                                agg.reset_defects_hist.size() - 1)]++;
                                        agg.handoff_events += r.n_handoff;
                                        agg.shots_with_handoff += r.n_handoff > 0;
                                        agg.omega_cubes_gt3 += r.omega_cubes > 3;
                                        agg.p2_conflicts += r.p2_conflict;
                                        agg.n_spanning_shots += r.n_spanning > 0;
                                        agg.seam_component_sum += r.n_seam_components;
                                        agg.spanning_component_sum += r.n_spanning;
                                        agg.obstacles_sum += r.n_obstacles;
                                        agg.promotions_sum += r.n_promotions;
                                        agg.moved_sum += r.n_moved_defects;
                                        agg.released_g_sum += r.n_released_g;
                                        agg.t1_sum += r.n_t1;
                                        agg.t2_sum += r.n_t2;
                                        agg.closure_sum += r.n_closure_steps;
                                        agg.stray_events += runner.teardown.stray_events;
                                        if (options.verify) {
                                            agg.verify_shots++;
                                            if (r.weight != r.ref_weight) {
                                                agg.weight_mismatches++;
                                                if (agg.first_weight_mismatch_shot < 0) {
                                                    agg.first_weight_mismatch_shot = (int64_t)shot_index;
                                                    log << "  VERIFY weight mismatch at shot " << shot_index
                                                        << ": weight=" << r.weight << " ref=" << r.ref_weight
                                                        << " truncated_cubes=" << r.n_cubes_truncated
                                                        << " groups=" << r.n_groups << "\n";
                                                }
                                            }
                                            if (r.obs_xor != r.ref_obs)
                                                agg.obs_mismatches++;
                                        }
                                    }
                                    done += batch;
                                    draw_progress(interactive, label, cell_index, total_cells, done, total);
                                }
                                if (interactive)
                                    std::fprintf(stderr, "\r\x1b[K");
                                agg.dbg_invariant_violations = runner.dbg_invariant_violations;
                                agg.dbg_cap_violations = runner.dbg_cap_violations;
                                agg.dbg_injection_violations = runner.dbg_injection_violations;
                                agg.dbg_feasibility_violations = runner.dbg_feasibility_violations;
                                agg.dbg_g_unsettled = runner.dbg_g_unsettled;
                                if (options.check_determinism) {
                                    agg.determinism_mismatches = determinism_mismatches;
                                    agg.determinism_check = determinism_mismatches == 0 ? 1 : 0;
                                }
                                log << "  cell " << label << " finished=" << stamp() << " shots=" << agg.shots
                                    << " shots_truncated=" << agg.shots_truncated
                                    << " cubes_truncated=" << agg.cubes_truncated << " p2_conflicts=" << agg.p2_conflicts
                                    << " handoff_events=" << agg.handoff_events;
                                if (options.verify)
                                    log << " verify_shots=" << agg.verify_shots
                                        << " weight_mismatches=" << agg.weight_mismatches
                                        << " obs_mismatches=" << agg.obs_mismatches;
                                if (options.check_determinism)
                                    log << " determinism_check="
                                        << (agg.determinism_check == 1 ? "pass" : "FAIL");
                                log << " debug(invariant=" << agg.dbg_invariant_violations
                                    << ",cap=" << agg.dbg_cap_violations
                                    << ",injection=" << agg.dbg_injection_violations
                                    << ",feasibility=" << agg.dbg_feasibility_violations
                                    << ",g_unsettled=" << agg.dbg_g_unsettled << ")"
                                    << " stray_events=" << agg.stray_events << "\n";
                                log.flush();
                                cells.push_back(std::move(agg));
                                write_agg();
                                const CellAggregate& last = cells.back();
                                if (last.dbg_invariant_violations || last.dbg_cap_violations ||
                                    last.dbg_injection_violations || last.dbg_feasibility_violations ||
                                    last.dbg_g_unsettled)
                                    throw std::runtime_error("a §9.1 debug invariant failed at " + label);
                                if (options.check_determinism && last.determinism_check != 1)
                                    throw std::runtime_error("--check-determinism failed at " + label);
                            }
                        }
                    }
                }
            }
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        options = parse_options(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }
    std::error_code dir_error;
    std::filesystem::create_directories(options.out_dir, dir_error);
    if (dir_error) {
        std::cerr << "error: could not create " << options.out_dir << ": " << dir_error.message() << "\n";
        return 1;
    }
    std::ofstream log(options.out_dir + "/run.log");
    if (!log.is_open()) {
        std::cerr << "error: could not open " << options.out_dir << "/run.log\n";
        return 1;
    }
    double overhead_ns = hires_timer_overhead_ns();
    bool thread_scoped = hires_timer_is_thread_scoped();

    log << "lattice_surgery_profiler - speculative decoding of a 1xN lattice-surgery strip with depth-1 seam fusion"
           " and local fallback resolution (docs/lattice_surgery_profiler.md)\n";
    log << "git_hash=" << SPECMATCHING_GIT_HASH << "\n";
    log << "build_flags=" << SPECMATCHING_BUILD_FLAGS << "\n";
    log << "stim_version=" << SPECMATCHING_STIM_VERSION << "\n";
    log << "command_line=" << options.command_line << "\n";
#if defined(__APPLE__)
    log << "platform=macos\n";
#else
    log << "platform=linux\n";
#endif
    log << "started=" << stamp() << "\n";
    log << "hires_timer_name=" << hires_timer_name() << "\n";
    log << "timer_is_thread_scoped=" << (thread_scoped ? 1 : 0) << "\n";
    log << "timer_overhead_ns=" << overhead_ns << " (10^6 back-to-back hires_now_ns() calls)\n";
    log << "rows_are_uncorrected=1\n";
    log << "threads=1 (no thread is ever spawned)\n";
    log << "model=lattice_surgery_k_core_model: every unit of work is executed on the profiling thread in a fixed"
           " dependency-respecting order and timed; a list scheduler then places the measurements on N(k+2)+1"
           " simulated cores (per cube: manager, k solvers, fallback; plus the global manager). t_sys is the makespan"
           " to g_combine. The scheduler's own dispatches (ready-heap pop + core pick) are timed as they run and"
           " charged to the dispatching manager.\n";
    log << "modelled_cores=N*(k+2)+1 (accelerators are a constant --accel-ns, on a per-cube pseudo-core)\n";
    log << "accel_ns=" << options.accel_ns << " link_ns=" << options.link_ns
        << " (link charged per transferred defect on seam fusions and phase-2 claims, per node on hand-offs)\n";
    log << "transfer_cost_modelled=" << (options.link_ns > 0 ? 1 : 0)
        << " warm_cache_caveat=1 evicting_truncated_components=1\n";
    log << "timed_region_input=compute_seeded_detection_events + build_ball_graph on the FULL strip's ball tables,"
           " the per-cube/per-seam edge split, the edge-list sort (--sort), the per-edge observable masks, the"
           " components (bookkeeping), the instance resets and output writing are all OUTSIDE every timed region\n";
    log << "edge_order=" << (options.sort_edges ? "sorted (--sort)" : "input_order") << "\n";
    log << "fallback=every cube's fallback core runs fb_phase1 from t=0 on every shot; its time is logged even when"
           " unused. Phase-1 territory Omega_0 = G_c u B_c (B_c: nodes within W of G_c), outer layer a weight-0"
           " virtual boundary; obstacles are defects OUTSIDE Omega_0 whose T-ball reaches it, each setting"
           " d(x,b)-T rounded down to even on its shell nodes inside Omega_0 (min with what is there).\n";
    log << "spare_pool=not used (deviation 3): H-only defects are re-solved in the speculative H instance\n";
    log << "baselines=B1 strip_pm: stock pymatching on the full strip G on one core (t_strip_pm_ns, ref_weight);"
           " B2 global_escalate: any cube truncated ? B1 : t_spec; B3 cube_pm_max: max over cubes of stock pymatching"
           " on G_c with seam edges turned into boundary edges (not exact; latency reference only)\n";
    log << "t_spec_ns=t_sys on a shot with no truncation, 0 otherwise; agg.json takes it over non-truncated shots\n";
    log << "t_p2_ns=start of the first p2_release to the end of the last p2_* job\n";
    log << "verify=" << (options.verify ? 1 : 0)
        << " (weight == B1's weight on every shot is the exactness test; mismatches are counted, not fatal; the"
           " observable mask is compared and counted separately, since degenerate matchings share a weight)\n";
    log << "check_determinism=" << (options.check_determinism ? 1 : 0)
        << " (every shot decoded twice; structural columns, weight and observables compared)\n";
    log << "jobs_log=" << (options.jobs_log == JobsLog::ALL ? "all" : options.jobs_log == JobsLog::SAMPLE ? "sample" : "none")
        << " (sample: the first " << JOBS_SAMPLE_FIRST << " shots and every shot that truncates)\n";
    log << "convexity_check=per (T, W) and cube, face-slab pairs within 2T+W, distance over P_c vs over G; sources"
           " capped at "
        << options.convexity_sources << " per cube (0 = all); allow_nonconvex=" << (options.allow_nonconvex ? 1 : 0)
        << ". Exactness here does NOT depend on it: (t1) fires before any region leaves Omega, so a pair whose"
           " geodesic leaves Omega is feasible by the outer-layer argument.\n";
    log << "merge_init_note=with --merge-init on, the decoded type's seam plaquettes have no first-round detector, so"
           " the window's first round is a time boundary along the seam with no pre-merge history behind it; the"
           " logical error rates of such a window are not those of a surgery inside a longer computation, and this"
           " binary does not report them\n";
    log << "vendored_edits=graph_flooder.{h,cc}: a masked boundary slot is invisible to the scan; the outer_layer flag"
           " and (t1) hook (outer_layer / outer_layer_hit / note_node_claimed). See docs/upstream_edits.md.\n";
    log << "deviation_1=slab defects are phase-1 defects; obstacles are the defects outside Omega_0 whose T-ball"
           " reaches it; an obstacle whose ball contains a phase-1 defect is promoted (Omega_0 grows by W around it)\n";
    log << "deviation_2=the outer layer of Omega_0 is a weight-0 virtual boundary, which gives Y(a) <= d(a, x) for every"
           " outer node x - the same inequality (t1) maintains in phase 2\n";
    log << "deviation_3=H-only defects are re-solved in the speculative H instance; the cross H-edges into Omega are"
           " masked; section 10 items 3 and 5 are not needed\n";
    log << "deviation_4=phase 2 claims every component with a defect in or within T of Omega AND every unfinished"
           " component; a seam component with a truncated cube's defects is not handed to a seam fusion\n";
    log << "deviation_5=adjacent truncated cubes, and cubes sharing an unfinished component, form one group run in"
           " the lowest cube's instance; the others' unfinished components enter through claims and (t2); two groups"
           " whose territories or claims meet are a counted conflict, resolved with B1's answer and latency\n";
    log << "deviation_6=a defect injected into G territory at an odd clock starts at radius 1 (even y-intercept)\n";
    log << "deviation_7=T, W and the preload width are in lattice edge weights (median discretised edge of G), as"
           " --T is in every other binary in this directory; w_p is logged per corpus\n";
    log << "deviation_8=time-ordered list schedule: per-dispatcher ready heaps ordered by (ready, larger cost, id)\n";
    log << "deviation_9=p2_resume runs G until settled or (t1), then each H solve until complete or (t2); first"
           " trigger ends the segment; resume_events counts the G side's matcher events\n";
    log << "deviation_10=no H region is released in place (n_released_h = 0): every component with a dummy toward a"
           " phase-2 territory is claimed and reset\n";
    log << "q_reference=an N=1 cell at the same (d, p, T, k) is this generator's d x d x d memory experiment; its"
           " cube_trunc_rate is the q to compare against\n";
    log << "seed=" << options.seed << " (shot i of a corpus is the same across the T, W and k sweeps)\n";
    log << "shots_per_cell=" << options.shots << " warmup_shots_per_cell=" << options.warmup << "\n";
    if (!thread_scoped)
        log << "WARNING: no thread-scoped clock; every tick is an upper bound only.\n";
    log.flush();

    bool interactive = stderr_is_terminal();
    bool ok = true;
    try {
        run_all(options, log, interactive);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "\n");
        std::cerr << "error: " << error.what() << "\n";
        log << "ERROR=" << error.what() << "\n";
        ok = false;
    }
    log << "\nfinished=" << stamp() << "\n";
    log << (ok ? "status=ok\n" : "status=incomplete (see ERROR above)\n");
    std::printf("wrote %s/{shots_*.csv,jobs_*.csv,phase2_*.csv,agg.json,run.log}\n", options.out_dir.c_str());
    return ok ? 0 : 1;
}
