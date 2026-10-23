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

/// `load_balancer_profiler` — the k-core critical-path latency of the decomposed `H` solve.
///
/// The component profiler priced one component of `H`; this binary asks what a machine with `k`
/// solver cores would pay for a whole shot once those components are spread over them. Per shot it
/// unions the shot's `H` into components, assigns the components to `k` cores by LPT under a
/// size-only cost model, builds and solves one sub-`H` per core, and logs
///
/// ```
///   sparse_k = uf + balance + max_c (build_c + solve_c + extract_c) + combine
///   system   = escalated ? fallback : min(fallback, sparse_k)
/// ```
///
/// `system` is defined by that and by nothing else (hard constraint 6), and every table derived
/// from `shots.csv` computes it from those two lines.
///
/// Beside the model, and without changing it, the CPU cost of producing the model's input is timed
/// too, so the whole software pipeline of a shot has a distribution:
///
/// ```
///   pipeline        = seed + hbuild + sparse_k
///   system_pipeline = escalated ? fallback : min(fallback, pipeline)
/// ```
///
/// `sparse_k` and `system` keep assuming the edge list arrives for free (from hardware);
/// `pipeline` and `system_pipeline` are what the same shot costs when `H` is built on this CPU.
///
/// **Exactness is untouched.** The one edit under `src/specmatching/` is the `ResetPolicy` parameter
/// on `Harvester::extract_only_to_*` described below; it defaults to the previous behaviour, changes
/// no production call site, and moves no work — it only lets a caller choose where the teardown
/// runs. The solve, the extraction and the escalating set are exactly as they were.
///
/// What changes between `k = 1` and `k = 8` is *which instance solves which
/// set of components*, and the escalation predicate is invariant under that: a core's timeline is
/// the disjoint union of its components' timelines, so the core truncates at `T` iff at least one
/// of its components does, and `any core TRUNCATED` is the same predicate as `any component
/// TRUNCATED`. `--verify` checks that against the monolithic per-component path shot by shot.
///
/// ## Execution model — `k_core_critical_path_model`
///
/// **Single thread; no thread is ever spawned** (hard constraint 2). The `k` loads are solved one
/// after another on the profiling thread, each on its own pre-allocated `BallMwpm` instance, each
/// bracketed by the thread-CPU timer of `hires_timer.h`; the per-core figure reported is the **max**
/// over loads, never the sum. Two assumptions come with every number:
///
///  1. Cores share nothing and start together, and no cost is modelled for moving a core's edge
///     slice to it. The `balance` region does the scatter on the profiling thread, so the scatter is
///     paid once serially rather than `k` times in parallel; that is conservative.
///  2. The fallback core and the `k` solver cores start at the same instant, when the edge list is
///     ready.
///
/// And the warm-cache caveat travels with the data: a load solved immediately after the previous
/// load sees a warmer cache than a real core would.
///
/// ## What is timed
///
/// Per shot, in this order, each region bracketed by one timer pair. Nothing else in the shot is
/// timed. The edge list — `compute_seeded_detection_events` plus `build_ball_graph` — is **the
/// input** of the `k`-core model: its two regions are timed on the CPU and reported in their own
/// columns, but they are in no term of `sparse_k` or `system`.
///
/// ```
///   [ seed ]      compute_seeded_detection_events: the negative-weight preamble on G, sorted dets
///   [ hbuild ]    build_ball_graph: the ball intersection that writes H, then its canonical sort.
///                 Split inside the call by BallGraphTiming into hbuild_intersect and
///                 hbuild_finalize, on the same thread-scoped clock; the split does not sum to
///                 hbuild (the call's setup and the inner timer's own reads are the gap)
///   [ uf ]        union-find over H's defect-defect edges; comp id per defect; sizes
///   [ balance ]   cost lookup, LPT assignment, scatter of defects and edges into per-core slices
///   for c in [0, k):                                     (a core with no components is skipped)
///     [ build_c ]    sub-H_c from core c's slice, then BallMwpm::rebuild, then the dets fill
///     [ solve_c ]    process_timeline_until_horizon(instance_c, dets_c, T_int)
///     [ extract_c ]  extract_only_to_obs(instance_c, dets_c)      -- COMPLETE cores only
///     abandon_shot(instance_c)                          -- TRUNCATED cores only     NOT timed
///   [ combine ]   XOR of the k partial obs masks, sum of the k partial weights, OR of the flags
///   [ fallback ]  pm::decode_detection_events on G from the full syndrome
///     the monolithic cross-check under --verify                                     NOT timed
/// ```
///
/// Three region definitions are stated here rather than left to be inferred, because a table read
/// under the wrong one is wrong in a way that looks plausible:
///
///  - **`build_c` contains `BallMwpm::rebuild`.** The design's `build_c` is "edges and boundary
///    edges copied with `H`'s weights, masks and cutoffs"; the *masks* are written by `rebuild`
///    (`mask_of` over the ball tables), not by the sub-graph copy, which carries only the canonical
///    ball `entry`. And a core that did not rebuild has no instance to solve on, so leaving the
///    rebuild in no region at all would make the critical path omit work the core provably does.
///    The dets fill — `dets_c = 0..n_c-1` — is in `build_c` for the same reason.
///  - **`extract_c` is the extraction alone — `reset_for_next_shot` is outside it.** Extraction
///    produces the observable mask and the weight, which are the decoder's answer; the reset
///    prepares the instance for the *next* shot, and nothing about this shot's answer depends on it,
///    so a pipelined decoder that double-buffered its instances would not pay it here. The two used
///    to be inseparable — `extract_only_impl` ended in the reset unconditionally — so this is a
///    **decoder change**: `Harvester::extract_only_to_*` now takes a `ResetPolicy`, defaulting to the
///    old behaviour so every production call site is untouched, and this profiler passes
///    `CALLER_RESETS` and runs the reset itself on the far side of the second clock read. The reset
///    still happens, in the same place in the same order; it is charged to no column.
///  - **`abandon_shot` is outside every region**, on the `TRUNCATED` branch, as the design asks. It
///    is `Mwpm::reset`, which sweeps the instance's whole node vector — and these instances are
///    sized for the full `H` — so a region containing it would be priced by `n_max` rather than by
///    the load.
///
/// A core with **no components assigned** does no work and opens no timer: its three regions are
/// `0`, and it contributes nothing to the makespan, which is what an idle core contributes in the
/// system being modelled. `run.log` records this as `empty_core_policy`.
///
/// Timer reads per **measured** shot: `3 (seed, hbuild) + 3 inside build_ball_graph + 2 (uf) +
/// 2 (balance) + 4 per non-empty core + 2 (combine) + 2 (fallback)`. Four per core rather than six because consecutive regions share their boundary
/// read; a `TRUNCATED` core reads three, since it has no extraction.
///
/// **A warm-up shot is not profiled.** It runs the identical workload — same builds, same solves,
/// same extractions, same fallback, same order, same buffers — and reads no clock of its own:
/// `run_shot` is templated on the switch, so the reads compile out of the warm-up instantiation
/// rather than being taken and discarded. The one exception is the three reads `build_ball_graph`
/// takes inside `hbuild`, which a warm-up shot also takes, so that the measured path is warm. Rows are stored **uncorrected** —
/// `timer_overhead_ns` goes to `run.log` and the analysis subtracts it, because a correction baked
/// into a row cannot be undone by a reader who disagrees with it.
///
/// ## `crit_core` is the measured argmax
///
/// §3.3 of the design defines `c* = argmax_c (build_c + solve_c + extract_c)`, ties to the lowest
/// index, and §7.1 defines `ticks_core_max` as that core's three regions summed. That is the `max_c`
/// of the latency model, so it is a function of the measured ticks and **not** deterministic across
/// runs — which contradicts §8's list, where `crit_core` sits in the deterministic set. §3.3 wins,
/// because §0's model is hard constraint 6 and `ticks_core_max` has to be the max for `sparse_k` to
/// mean what §0 says it means.
///
/// So `--check-determinism` diffs the columns that really are deterministic — `n_def`, `n_comps`,
/// `largest`, `escalated`, `pred_max`, `pred_sum`, `imbalance` — and, instead of the three tick-
/// derived `crit_*` columns, the **whole LPT assignment**: `pred_cost[c]`, `sizes[c]` and
/// `n_comps[c]` for every core. That is the deterministic content §8 is about, and checking all `k`
/// cores is a stronger statement than checking the one that happened to be slowest.
///
/// ## Units
///
/// `--T` is in **multiples of one lattice edge weight**, as in `sparse_graph_stats` and
/// `component_profiler` — the normalisation under which the M1 structural results are
/// `d`-independent, and the one every horizon in `benchmarks/spec_matching/` is quoted in. The
/// design's §4 quotes `--T` in raw DEM float weight units; that reading is a deviation recorded
/// rather than silently taken, and it is the same one the two companions record, for the same
/// reason: 2 raw DEM float units is a fraction of an edge at `p = 1e-3` and would truncate every
/// component. All three readings — `T` (edge weights), `T_weight_units` (DEM float) and `T_int`
/// (the flooder's integer time units, read off the decoder rather than reconverted) — go into
/// `run.log`.
///
/// The file also lives in `benchmarks/spec_matching/` rather than the design's
/// `benchmarks/two_phase/`, which is this tree's name for that corpus; the generator call is the
/// corpus's own (see `generator_call` for its two inherited deviations).
///
/// ## Usage
///
///   load_balancer_profiler [--d 5,7,9] [--p 1e-3,5e-4] [--T 1.5,2] [--k 1,2,4,8] [--alpha 1]
///                          [--shots N] [--warmup 1000] [--seed S] [--verify]
///                          [--check-determinism] [--out DIR]
///
/// Writes `shots.csv`, `agg.json` and `run.log` to `--out`. The agent runs no campaigns; smoke runs
/// at `--shots 100 --d 5` exist only to check that the output parses, and their numbers appear in
/// no document.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "benchmarks/spec_matching/hires_timer.h"
#include "benchmarks/spec_matching/profiler_util.h"
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

namespace {

/// Shots sampled per batch, as in the two companions. The stream is a function of the seed and the
/// batch sizes alone, and both are identical across the `T` and `k` sweeps, which is what makes shot
/// `i` of `(d, p)` the same shot in every cell.
constexpr size_t SHOT_BATCH = 4096;

/// The fixed row arena. One row per shot, reserved once at startup and never grown, drained to
/// `shots.csv` between batches and at the end of every cell — from between-shot positions only,
/// never from inside a timed region. Overflow is impossible in practice (a batch is 4096 rows
/// against a flush threshold of 65536) but is handled without allocating anyway: the arena stops
/// taking rows and `agg.json` records `raw_truncated_at_shot`.
constexpr size_t ROW_ARENA_CAPACITY = 1u << 17;
constexpr size_t ROW_ARENA_FLUSH_AT = ROW_ARENA_CAPACITY / 2;

/// Initial per-core edge-slice and sub-`H` edge capacity. `|E(H)|` has no static bound — it is
/// `O(n^2)` in the worst case — so this is a starting point that the warm-up grows to the corpus's
/// high-water mark, exactly as `BallGraphArena` and `BallMwpm` do. `run.log` reports
/// `no_growth_after_warmup` per cell, which is how invariant 11 is stated in this tree: a cell whose
/// warm-up was too short for its corpus says so in its own artifact rather than being asserted away.
constexpr size_t INITIAL_EDGE_SLICE_CAPACITY = 4096;

struct Options {
    std::vector<size_t> distances = {5};
    std::vector<double> error_rates = {1e-3};
    std::vector<double> horizons = {1.5, 2.0};
    std::vector<size_t> core_counts = {1, 2, 4, 8};
    /// The balancer's cost exponent, `cost(s) = s^alpha`. The design's §3.2 specifies a default of
    /// 1.1, fitted to per-component solve time; this tree defaults to **1.0** — a linear cost model —
    /// on instruction. `run.log` records whichever value a run actually used.
    double alpha = 1.0;
    size_t shots = 10000;
    size_t warmup = 1000;
    uint64_t seed = 20260906;
    bool verify = false;
    bool check_determinism = false;
    std::string out_dir = "benchmarks/spec_matching/results/load_balancer_profiler";
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; i++) {
        std::string flag = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument("missing value for " + flag);
            return argv[++i];
        };
        if (flag == "--d" || flag == "--distances") {
            options.distances = parse_list<size_t>(next());
        } else if (flag == "--p" || flag == "--error-rates") {
            options.error_rates = parse_list<double>(next());
        } else if (flag == "--T" || flag == "--horizons") {
            options.horizons = parse_list<double>(next());
        } else if (flag == "--k" || flag == "--cores") {
            options.core_counts = parse_list<size_t>(next());
        } else if (flag == "--alpha") {
            options.alpha = std::stod(next());
        } else if (flag == "--shots") {
            options.shots = std::stoul(next());
        } else if (flag == "--warmup") {
            options.warmup = std::stoul(next());
        } else if (flag == "--seed") {
            options.seed = std::stoull(next());
        } else if (flag == "--verify") {
            options.verify = true;
        } else if (flag == "--check-determinism") {
            options.check_determinism = true;
        } else if (flag == "--out" || flag == "--out-dir") {
            options.out_dir = next();
        } else {
            throw std::invalid_argument("unrecognised flag " + flag);
        }
    }
    if (options.distances.empty() || options.error_rates.empty() || options.horizons.empty() ||
        options.core_counts.empty())
        throw std::invalid_argument("--d, --p, --T and --k each need at least one value");
    if (options.shots == 0)
        throw std::invalid_argument("--shots must be positive");
    for (size_t distance : options.distances) {
        if (distance < 3 || distance % 2 == 0)
            throw std::invalid_argument("--d values must be odd and at least 3");
    }
    for (double noise : options.error_rates) {
        if (!(noise > 0) || !(noise < 1))
            throw std::invalid_argument("--p values must lie in (0, 1)");
    }
    for (double horizon : options.horizons) {
        if (!(horizon > 0))
            throw std::invalid_argument("--T values must be positive");
    }
    for (size_t k : options.core_counts) {
        if (k == 0)
            throw std::invalid_argument("--k values must be positive");
    }
    if (!(options.alpha > 0))
        throw std::invalid_argument("--alpha must be positive");
    // Descending `k` would make the "same escalating set across k" cross-check read oddly (it
    // compares against the first `k` of the cell), and a duplicate `k` would write two identical
    // cells. Neither is an error worth failing on, so both are normalised here and recorded.
    std::sort(options.core_counts.begin(), options.core_counts.end());
    options.core_counts.erase(
        std::unique(options.core_counts.begin(), options.core_counts.end()), options.core_counts.end());
    return options;
}

/// `%.10g` — enough to round-trip the balancer's doubles for a reader, deterministic for a fixed
/// binary, and short enough that a `pred_sum` column does not dominate the file.
std::string fmt_g(double value) {
    char buffer[40];
    std::snprintf(buffer, sizeof(buffer), "%.10g", value);
    return std::string(buffer);
}

/// `__int128` has no stream inserter. The 128-bit accumulators of the aggregate are printed as JSON
/// *strings*, because a sum of squares of nanosecond ticks over 10^6 shots leaves the range a JSON
/// number is guaranteed to survive.
std::string fmt_i128(__int128 value) {
    if (value == 0)
        return "0";
    bool negative = value < 0;
    unsigned __int128 magnitude = negative ? (unsigned __int128)(-value) : (unsigned __int128)value;
    char buffer[48];
    size_t at = sizeof(buffer);
    while (magnitude != 0) {
        buffer[--at] = (char)('0' + (int)(magnitude % 10));
        magnitude /= 10;
    }
    std::string out = negative ? "-" : "";
    out.append(buffer + at, sizeof(buffer) - at);
    return out;
}

/// FNV-1a over the deterministic columns of a cell. Two runs of the same `(seed, d, p, T, k, alpha)`
/// must agree on it; see the header for what is in the set and why `crit_*` is not.
struct Digest {
    uint64_t value{1469598103934665603ull};

    void bytes(const void* pointer, size_t count) {
        const uint8_t* at = (const uint8_t*)pointer;
        for (size_t i = 0; i < count; i++) {
            value ^= at[i];
            value *= 1099511628211ull;
        }
    }
    template <typename T>
    void add(const T& item) {
        bytes(&item, sizeof(T));
    }
};

/// One shot's row of `shots.csv`.
///
/// `d`, `p`, `T`, `k` and `alpha` are cell constants and are stamped on at flush time rather than
/// stored per row.
struct Row {
    uint64_t shot;
    uint32_t n_def;
    uint32_t n_comps;
    uint32_t largest;
    uint8_t escalated;
    double pred_max;
    double pred_sum;
    double imbalance;
    uint32_t crit_core;
    uint32_t crit_size;
    uint32_t crit_ncomps;
    uint64_t uf;
    uint64_t balance;
    uint64_t build_crit;
    uint64_t solve_crit;
    uint64_t extract_crit;
    uint64_t core_max;
    uint64_t combine;
    uint64_t sparse_k;
    uint64_t fallback;
    uint64_t system;
    uint64_t seed;
    uint64_t hbuild;
    uint64_t hbuild_intersect;
    uint64_t hbuild_finalize;
    uint64_t pipeline;
    uint64_t system_pipeline;
};

struct CellLabels {
    size_t distance{0};
    double noise{0};
    double horizon{0};
    size_t cores{0};
    double alpha{0};
};

struct RowArena {
    std::vector<Row> rows;
    CellLabels labels;
    std::ofstream* out{nullptr};
    uint64_t rows_this_cell{0};
    uint64_t rows_total{0};
    /// `-1` unless the arena refused a row this cell; then the 0-based shot index of the first one
    /// it refused. Reported in `agg.json`, never worked around by allocating.
    int64_t truncated_at_shot{-1};

    void configure(std::ofstream& stream) {
        out = &stream;
        rows.reserve(ROW_ARENA_CAPACITY);
    }

    inline void append(const Row& row) {
        if (rows.size() >= ROW_ARENA_CAPACITY) {
            if (truncated_at_shot < 0)
                truncated_at_shot = (int64_t)row.shot;
            return;
        }
        rows.push_back(row);
    }

    /// Called from between-shot positions only.
    void flush() {
        for (const Row& row : rows) {
            *out << labels.distance << "," << fmt_g(labels.noise) << "," << fmt_g(labels.horizon) << ","
                 << labels.cores << "," << fmt_g(labels.alpha) << "," << row.shot << "," << row.n_def << ","
                 << row.n_comps << "," << row.largest << "," << (int)row.escalated << "," << fmt_g(row.pred_max) << ","
                 << fmt_g(row.pred_sum) << "," << fmt_g(row.imbalance) << "," << row.crit_core << "," << row.crit_size
                 << "," << row.crit_ncomps << "," << row.uf << "," << row.balance << "," << row.build_crit << ","
                 << row.solve_crit << "," << row.extract_crit << "," << row.core_max << "," << row.combine << ","
                 << row.sparse_k << "," << row.fallback << "," << row.system << "," << row.seed << ","
                 << row.hbuild << "," << row.hbuild_intersect << "," << row.hbuild_finalize << "," << row.pipeline
                 << "," << row.system_pipeline << "\n";
        }
        rows_this_cell += rows.size();
        rows_total += rows.size();
        rows.clear();
        out->flush();
    }
};

/// The three latency series of §0, accumulated in 128-bit integers with the raw tick samples kept
/// for the percentiles. No floating point touches the accumulation path; the means and the speedup
/// are computed once, at print time.
struct TickSeries {
    __int128 sum{0};
    __int128 sum_sq{0};
    uint64_t max{0};
    std::vector<long long> samples;

    void reserve(size_t count) {
        samples.reserve(count);
    }
    void clear() {
        sum = 0;
        sum_sq = 0;
        max = 0;
        samples.clear();
    }
    void add(uint64_t ticks) {
        sum += (__int128)ticks;
        sum_sq += (__int128)ticks * (__int128)ticks;
        max = std::max(max, ticks);
        samples.push_back((long long)ticks);
    }
    double mean() const {
        return samples.empty() ? 0.0 : (double)sum / (double)samples.size();
    }
};

/// One `(d, p, T, k)` cell of `agg.json`.
struct CellAggregate {
    size_t distance{0};
    double noise{0};
    double horizon{0};
    size_t cores{0};
    double alpha{0};
    uint64_t shots{0};
    uint64_t escalations{0};

    TickSeries fallback;
    TickSeries sparse_k;
    TickSeries system;
    /// The CPU input pipeline: `H`'s construction alone, and the two whole-pipeline series.
    TickSeries hbuild;
    TickSeries pipeline;
    TickSeries system_pipeline;

    __int128 uf_sum{0};
    __int128 balance_sum{0};
    __int128 build_crit_sum{0};
    __int128 solve_crit_sum{0};
    __int128 extract_crit_sum{0};
    __int128 combine_sum{0};
    __int128 seed_sum{0};
    __int128 hbuild_intersect_sum{0};
    __int128 hbuild_finalize_sum{0};

    double imbalance_sum{0};

    /// Indexed by value; only the non-zero entries are written out.
    std::vector<uint64_t> crit_ncomps_hist;
    std::vector<uint64_t> largest_hist;

    int64_t raw_truncated_at_shot{-1};
    /// `1` / `0` / `-1` for pass / fail / not run.
    int determinism_check{-1};
    /// Same, for "the escalating set is the one the first `k` of this `(d, p, T)` produced".
    int escalation_matches_first_k{-1};

    void reset(size_t histogram_size) {
        shots = 0;
        escalations = 0;
        fallback.clear();
        sparse_k.clear();
        system.clear();
        hbuild.clear();
        pipeline.clear();
        system_pipeline.clear();
        seed_sum = 0;
        hbuild_intersect_sum = 0;
        hbuild_finalize_sum = 0;
        uf_sum = 0;
        balance_sum = 0;
        build_crit_sum = 0;
        solve_crit_sum = 0;
        extract_crit_sum = 0;
        combine_sum = 0;
        imbalance_sum = 0;
        crit_ncomps_hist.assign(histogram_size, 0);
        largest_hist.assign(histogram_size, 0);
        raw_truncated_at_shot = -1;
        determinism_check = -1;
        escalation_matches_first_k = -1;
    }
};

/// One of the `k` emulated solver cores: its instance, its slice of the shot, and the three ticks it
/// contributed. Constructed once per `(d, p)` and reused for every shot of every `(T, k)` cell.
struct SolverCore {
    BallMwpm instance;
    Harvester harvester;
    /// The core's sub-`H`, rewritten in `build_c`. Its own buffer rather than a shared one, so a
    /// core's build writes into memory no other core touches — which is what the model assumes.
    BallGraph sub;
    /// `0..n_c-1`: the sub-`H`'s nodes are exactly its detection events, by construction.
    std::vector<uint64_t> dets;
    HarvestResult part;

    /// The `balance` scatter's output: indices into `H`'s node, edge and boundary-edge lists,
    /// ascending in all three.
    std::vector<uint32_t> slice_nodes;
    std::vector<uint32_t> slice_edges;
    std::vector<uint32_t> slice_boundary;

    uint32_t ncomps{0};
    double pred_cost{0};

    uint64_t build_ticks{0};
    uint64_t solve_ticks{0};
    uint64_t extract_ticks{0};
    bool truncated{false};

    void reserve(size_t n_max) {
        sub.h_to_det.reserve(n_max);
        sub.edges.reserve(INITIAL_EDGE_SLICE_CAPACITY);
        sub.boundary_edges.reserve(n_max);
        dets.reserve(n_max);
        slice_nodes.reserve(n_max);
        slice_edges.reserve(INITIAL_EDGE_SLICE_CAPACITY);
        slice_boundary.reserve(n_max);
    }

    /// Total capacity across the core's own buffers, so that "did a measured shot allocate?" can be
    /// answered by comparing one number before and after.
    uint64_t capacity_signature() const {
        return (uint64_t)sub.h_to_det.capacity() + sub.edges.capacity() + sub.boundary_edges.capacity() +
               dets.capacity() + slice_nodes.capacity() + slice_edges.capacity() + slice_boundary.capacity() +
               instance.grow_events;
    }
};

/// The serial front half of a shot: the union-find of §3.1 and the LPT balancer of §3.2, over fixed
/// arenas sized to `n_max` at cell setup.
struct Balancer {
    /// Union-find over `H`'s defect-defect edges. Weighted union (by set size) with path halving.
    std::vector<uint32_t> parent;
    std::vector<uint32_t> set_size;

    /// Component index of each `H` node, and the per-component tables. Components are numbered in
    /// **ascending minimum `H`-node index**, which is ascending minimum detector id because
    /// `h_to_det` ascends — the tie-break §3.2's LPT order needs.
    std::vector<uint32_t> comp_id;
    std::vector<uint32_t> comp_size;
    std::vector<uint32_t> comp_min;
    std::vector<uint32_t> core_of_comp;

    /// `comp_index_of_root[r]` is valid iff `root_stamp[r] == epoch`, so nothing has to be cleared
    /// between shots. `epoch` is bumped once per shot, inside the `uf` region.
    std::vector<uint32_t> comp_index_of_root;
    std::vector<uint64_t> root_stamp;
    uint64_t epoch{0};

    /// Counting sort by component size: `order` is the components in descending size, ties ascending
    /// minimum defect id. `counter` is zero on entry to every shot and left zero on exit — only the
    /// sizes actually present are touched, so no `O(n_max)` clear is paid per shot.
    std::vector<uint32_t> order;
    std::vector<uint32_t> counter;

    /// `cost_tab[s] = s^alpha`, filled once per `(d, p)` with `pow`. No `pow` call per shot
    /// (hard constraint 4).
    std::vector<double> cost_tab;

    /// Position of each `H` node within its core's sub-`H`. Written in `build_c`; a node belongs to
    /// exactly one core, so the cores never collide in it.
    std::vector<uint32_t> local_index;

    uint32_t n_comps{0};
    uint32_t largest{0};

    void configure(size_t n_max, double alpha) {
        parent.assign(n_max, 0);
        set_size.assign(n_max, 0);
        comp_id.assign(n_max, 0);
        comp_size.assign(n_max, 0);
        comp_min.assign(n_max, 0);
        core_of_comp.assign(n_max, 0);
        comp_index_of_root.assign(n_max, 0);
        root_stamp.assign(n_max, 0);
        order.assign(n_max, 0);
        counter.assign(n_max + 2, 0);
        local_index.assign(n_max, 0);
        epoch = 0;
        cost_tab.assign(n_max + 1, 0.0);
        for (size_t s = 1; s <= n_max; s++)
            cost_tab[s] = std::pow((double)s, alpha);
    }

    inline uint32_t find(uint32_t node) {
        while (parent[node] != node) {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        return node;
    }
};

/// Everything a shot measures, kept in one place so that the row writer, the aggregate and the
/// determinism digest all read the same values rather than three recomputations of them.
struct ShotMeasurement {
    uint32_t n_def{0};
    uint32_t n_comps{0};
    uint32_t largest{0};
    uint8_t escalated{0};
    double pred_max{0};
    double pred_sum{0};
    double imbalance{0};
    uint32_t crit_core{0};
    uint32_t crit_size{0};
    uint32_t crit_ncomps{0};
    uint64_t uf{0};
    uint64_t balance{0};
    uint64_t build_crit{0};
    uint64_t solve_crit{0};
    uint64_t extract_crit{0};
    uint64_t core_max{0};
    uint64_t combine{0};
    uint64_t sparse_k{0};
    uint64_t fallback{0};
    uint64_t system{0};
    uint64_t seed{0};
    uint64_t hbuild{0};
    uint64_t hbuild_intersect{0};
    uint64_t hbuild_finalize{0};
    uint64_t pipeline{0};
    uint64_t system_pipeline{0};
    /// The combined Phase-1 answer, for `--verify`. Not a column of `shots.csv`.
    pm::obs_int obs_mask{0};
    pm::total_weight_int weight{0};
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

constexpr size_t PROGRESS_BAR_WIDTH = 24;

/// The campaign's progress: a bar over the current cell's shots, and the cell's place in the grid.
/// No ETA and no shots-per-second — every number on the line is a ratio of counts the loop already
/// keeps, so drawing it reads no clock. Drawn between batches, to stderr, from outside every timed
/// region.
void draw_progress(
    bool interactive, const char* label, size_t cell_index, size_t total_cells, size_t done, size_t cell_shots) {
    double fraction = cell_shots ? (double)done / (double)cell_shots : 1.0;
    if (!interactive) {
        std::fprintf(
            stderr, "[%zu/%zu] %s  %5.1f%%  %zu/%zu shots\n", cell_index, total_cells, label, 100.0 * fraction, done,
            cell_shots);
        std::fflush(stderr);
        return;
    }
    char bar[PROGRESS_BAR_WIDTH + 1];
    size_t filled = (size_t)(fraction * (double)PROGRESS_BAR_WIDTH);
    for (size_t i = 0; i < PROGRESS_BAR_WIDTH; i++)
        bar[i] = i < filled ? '#' : '.';
    bar[PROGRESS_BAR_WIDTH] = '\0';
    std::fprintf(
        stderr, "\r[%zu/%zu] %s  [%s] %5.1f%%  %zu/%zu shots\x1b[K", cell_index, total_cells, label, bar,
        100.0 * fraction, done, cell_shots);
    std::fflush(stderr);
}

void draw_setup(bool interactive, const char* label, size_t cell_index, size_t total_cells) {
    std::fprintf(
        stderr,
        interactive ? "\r[%zu/%zu] %s  building DEM, ball tables and solver instances...\x1b[K"
                    : "[%zu/%zu] %s  building DEM, ball tables and solver instances...\n",
        cell_index, total_cells, label);
    std::fflush(stderr);
}

/// The generator call, as one string, so `run.log` records what `(d, p)` actually meant.
///
/// It is `profiler_util.h`'s `ShotSampler::make` — the same call the rest of the
/// `benchmarks/spec_matching/` corpus is generated by, which is the requirement: `(d, p)` has to
/// mean here what it means in every table beside it.
///
/// Two deviations from the design's prose, recorded rather than silently taken, and both inherited
/// from that shared call: the basis is `rotated_memory_x`, not `_z`; and three noise parameters are
/// set, not four (`before_round_data_depolarization` is left at 0).
std::string generator_call(size_t distance, size_t rounds, double noise) {
    std::ostringstream out;
    out << "stim::CircuitGenParameters(rounds=" << rounds << ", distance=" << distance
        << ", task=\"rotated_memory_x\"); after_clifford_depolarization=" << noise
        << "; after_reset_flip_probability=" << noise << "; before_measure_flip_probability=" << noise
        << "; before_round_data_depolarization=0 (unset)"
        << "; stim::generate_surface_code_circuit(...).circuit"
        << "; stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false)";
    return out.str();
}

/// Bytes an instance's matching graph is holding: the node pool plus the four parallel adjacency
/// arrays of every node, taken from the containers' own capacities rather than estimated. Read once
/// per cell, outside every timed region, so that `run.log` can state what `k` instances cost.
uint64_t instance_bytes(const BallMwpm& core) {
    const pm::MatchingGraph& graph = core.mwpm.flooder.graph;
    // Element sizes are taken from the containers themselves rather than named, so that widening
    // `pm::weight_int` or reshaping an implied-weight rule cannot silently make this a guess.
    uint64_t total = (uint64_t)graph.nodes.capacity() * sizeof(decltype(graph.nodes)::value_type);
    for (const pm::DetectorNode& node : graph.nodes) {
        total += (uint64_t)node.neighbors.capacity() * sizeof(decltype(node.neighbors)::value_type);
        total += (uint64_t)node.neighbor_weights.capacity() * sizeof(decltype(node.neighbor_weights)::value_type);
        total +=
            (uint64_t)node.neighbor_observables.capacity() * sizeof(decltype(node.neighbor_observables)::value_type);
        total += (uint64_t)node.neighbor_implied_weights.capacity() *
                 sizeof(decltype(node.neighbor_implied_weights)::value_type);
        for (const auto& implied : node.neighbor_implied_weights)
            total += (uint64_t)implied.capacity() * sizeof(typename std::decay_t<decltype(implied)>::value_type);
    }
    return total;
}

std::string stamp() {
    std::time_t now = std::time(nullptr);
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&now));
    return std::string(buffer);
}

/// The whole of one `(d, p, T, k)` cell's per-shot work, and the only place a clock is read.
struct CellRunner {
    BallDecoder* decoder{nullptr};
    pm::Mwpm* fallback_instance{nullptr};
    Balancer* balancer{nullptr};
    std::vector<SolverCore>* cores{nullptr};
    size_t k{0};

    /// The fallback's output buffer. Zeroed before the region, never inside it.
    std::vector<uint8_t> obs_bytes;
    /// LPT accumulators, `k_max` long.
    std::vector<double> core_load;

    void configure(size_t num_observables, size_t k_max) {
        obs_bytes.assign(std::max<size_t>(1, num_observables), 0);
        core_load.assign(k_max, 0.0);
    }

    /// One shot, start to finish. `out` is filled with everything the row, the aggregate and the
    /// digest need; nothing else is kept.
    ///
    /// `Profile` is the warm-up switch, and it is a template parameter rather than a flag so that a
    /// measured shot carries **no** branch the profiled path did not already have. A warm-up shot
    /// runs the identical workload — the same builds, the same solves, the same extractions, the
    /// same fallback, in the same order, touching the same buffers — and reads no clock at all: the
    /// `tick()` calls below compile away entirely, so warm-up costs no timer reads and produces no
    /// measurement. Only measured shots are profiled.
    template <bool Profile>
    void run_shot(const std::vector<uint64_t>& shot, ShotMeasurement& out) {
        Balancer& bal = *balancer;
        std::vector<SolverCore>& core = *cores;
        auto tick = []() -> uint64_t {
            if constexpr (Profile)
                return hires_now_ns();
            else
                return 0;
        };

        // ---- The input. `compute_seeded_detection_events` and `build_ball_graph` are the ball-table
        // generator's edge list, and hard constraint 3 keeps producing it out of `sparse_k`. It is
        // still timed here, on the CPU, into its own columns: `seed`, `hbuild`, and `hbuild`'s
        // two-stage split, which `build_ball_graph` measures itself when handed a `BallGraphTiming`.
        // A warm-up shot hands it one too and discards it: the timing branch inside the call is then
        // warm when the first measured shot takes it, where with `nullptr` that shot paid a cold
        // first pass (~40 us at d = 5). These three reads are the only ones a warm-up shot takes.
        BallGraphTiming hbuild_timing;
        uint64_t seed_started = tick();
        decoder->compute_seeded_detection_events(shot, decoder->seeded_scratch);
        uint64_t seed_ended = tick();
        build_ball_graph(
            decoder->tables,
            decoder->seeded_scratch,
            decoder->horizon,
            decoder->arena,
            decoder->config.mode,
            &hbuild_timing,
            nullptr);
        uint64_t hbuild_ended = tick();
        const BallGraph& h = decoder->arena.graph;
        uint32_t n = (uint32_t)h.num_nodes();

        // ---- §3.1. Union-find over `H`'s defect-defect edges, and only those: the boundary is not a
        // node of `H`, so two defects that both reach it are not thereby connected, and a defect with
        // only boundary edges is a size-1 component.
        uint64_t uf_started = tick();
        for (uint32_t i = 0; i < n; i++) {
            bal.parent[i] = i;
            bal.set_size[i] = 1;
        }
        for (const BallGraphEdge& edge : h.edges) {
            uint32_t a = bal.find(edge.i);
            uint32_t b = bal.find(edge.j);
            if (a == b)
                continue;
            // Weighted union: the smaller set is linked under the larger, ties to the smaller index
            // so that the link is a function of `H` alone.
            if (bal.set_size[a] < bal.set_size[b] || (bal.set_size[a] == bal.set_size[b] && b < a))
                std::swap(a, b);
            bal.parent[b] = a;
            bal.set_size[a] += bal.set_size[b];
        }
        // Components numbered in ascending minimum member. Iterating `v` upwards, the first `v` of a
        // set *is* its minimum, so the numbering falls out of the pass with no sort; that is §3.2's
        // "ties by ascending minimum defect id", since `h_to_det` ascends with the node index.
        bal.epoch++;
        bal.n_comps = 0;
        bal.largest = 0;
        for (uint32_t v = 0; v < n; v++) {
            uint32_t root = bal.find(v);
            if (bal.root_stamp[root] != bal.epoch) {
                bal.root_stamp[root] = bal.epoch;
                bal.comp_index_of_root[root] = bal.n_comps;
                bal.comp_min[bal.n_comps] = v;
                bal.comp_size[bal.n_comps] = 0;
                bal.n_comps++;
            }
            uint32_t c = bal.comp_index_of_root[root];
            bal.comp_id[v] = c;
            bal.comp_size[c]++;
            bal.largest = std::max(bal.largest, bal.comp_size[c]);
        }
        uint64_t uf_ended = tick();

        // ---- §3.2. Cost lookup, LPT assignment, scatter. The balancer's *statistics* are read after
        // the region closes (hard constraint 3); only the accumulators are written inside it.
        uint64_t balance_started = tick();
        for (size_t c = 0; c < k; c++) {
            core_load[c] = 0.0;
            core[c].ncomps = 0;
            core[c].slice_nodes.clear();
            core[c].slice_edges.clear();
            core[c].slice_boundary.clear();
        }
        if (bal.n_comps > 0) {
            // Counting sort on size into the fixed `counter` arena: descending size, ties ascending
            // minimum defect id (the components are already numbered in that order).
            for (uint32_t c = 0; c < bal.n_comps; c++)
                bal.counter[bal.comp_size[c]]++;
            uint32_t running = 0;
            for (uint32_t s = bal.largest; s >= 1; s--) {
                uint32_t count = bal.counter[s];
                bal.counter[s] = running;
                running += count;
            }
            for (uint32_t c = 0; c < bal.n_comps; c++)
                bal.order[bal.counter[bal.comp_size[c]]++] = c;
            // Left zero for the next shot. The clear has to cover **every** `s` in `[1, largest]`,
            // not just the sizes this shot actually had: the prefix scan above writes a cursor into
            // `counter[s]` for every `s` it walks, absent sizes included, and an absent size's cursor
            // is nonzero whenever a larger component exists. Clearing only the present sizes leaves
            // that residue behind, and the next shot's `counter[...]++` then starts from it and
            // writes past the end of `order`. This is `O(largest)`, the same order as the prefix
            // scan it undoes, so it costs the sort nothing asymptotically.
            for (uint32_t s = 1; s <= bal.largest; s++)
                bal.counter[s] = 0;

            // LPT: longest first, each to the least-loaded core, ties to the lowest core index.
            for (uint32_t position = 0; position < bal.n_comps; position++) {
                uint32_t c = bal.order[position];
                size_t best = 0;
                for (size_t candidate = 1; candidate < k; candidate++) {
                    if (core_load[candidate] < core_load[best])
                        best = candidate;
                }
                bal.core_of_comp[c] = (uint32_t)best;
                core_load[best] += bal.cost_tab[bal.comp_size[c]];
                core[best].ncomps++;
            }

            // Scatter: one pass over the defects, one over the edges, one over the boundary edges.
            // All three of `H`'s lists are ascending and both endpoints of an `H` edge are in one
            // component, so every slice comes out ascending too — which is what lets `build_c`
            // relabel by position and inherit `H`'s `(i, j)` sort for free.
            for (uint32_t v = 0; v < n; v++)
                core[bal.core_of_comp[bal.comp_id[v]]].slice_nodes.push_back(v);
            for (uint32_t e = 0; e < (uint32_t)h.edges.size(); e++)
                core[bal.core_of_comp[bal.comp_id[h.edges[e].i]]].slice_edges.push_back(e);
            for (uint32_t b = 0; b < (uint32_t)h.boundary_edges.size(); b++)
                core[bal.core_of_comp[bal.comp_id[h.boundary_edges[b].i]]].slice_boundary.push_back(b);
        }
        uint64_t balance_ended = tick();

        // ---- §3.2's recorded-but-untimed statistics, read off the accumulators now that the region
        // has closed.
        out.pred_max = 0;
        out.pred_sum = 0;
        for (size_t c = 0; c < k; c++) {
            core[c].pred_cost = core_load[c];
            out.pred_max = std::max(out.pred_max, core_load[c]);
            out.pred_sum += core_load[c];
        }
        out.imbalance = out.pred_sum > 0 ? (double)k * out.pred_max / out.pred_sum : 1.0;

        // ---- §3.3. The `k` loads, solved one after another, each on its own instance. The figure
        // that goes into the model is the **max** over them, never the sum.
        out.crit_core = 0;
        out.core_max = 0;
        for (size_t c = 0; c < k; c++) {
            SolverCore& load = core[c];
            load.build_ticks = 0;
            load.solve_ticks = 0;
            load.extract_ticks = 0;
            load.truncated = false;
            // An idle core does nothing and is charged nothing, not even a timer read.
            if (load.slice_nodes.empty()) {
                load.part.committed = pm::MatchingResult();
                continue;
            }

            // ---- `build_c`. Sub-`H_c` from the slice, then the instance on it, then the dets.
            uint64_t build_started = tick();
            BallGraph& sub = load.sub;
            sub.clear();
            uint32_t count = (uint32_t)load.slice_nodes.size();
            for (uint32_t j = 0; j < count; j++) {
                uint32_t v = load.slice_nodes[j];
                bal.local_index[v] = j;
                sub.h_to_det.push_back(h.h_to_det[v]);
            }
            for (uint32_t index : load.slice_edges) {
                const BallGraphEdge& edge = h.edges[index];
                sub.edges.push_back(
                    BallGraphEdge{bal.local_index[edge.i], bal.local_index[edge.j], edge.w_int, edge.entry});
            }
            for (uint32_t index : load.slice_boundary) {
                const BallBoundaryEdge& edge = h.boundary_edges[index];
                sub.boundary_edges.push_back(BallBoundaryEdge{bal.local_index[edge.i], edge.w_int, edge.det});
            }
            load.instance.rebuild(decoder->tables, sub, decoder->arena, nullptr);
            // `dets` always holds `0, 1, ..., size - 1`, so shrinking is a truncation and growing is
            // an append: no `O(n_c)` refill, and nothing allocates once the high-water mark is set.
            if (load.dets.size() > count) {
                load.dets.resize(count);
            } else {
                while (load.dets.size() < count)
                    load.dets.push_back((uint64_t)load.dets.size());
            }
            uint64_t build_ended = tick();

            // ---- `solve_c`. Truncated sparse blossom on the core's whole load, which is the
            // disjoint union of its components — the same event sequence, component by component,
            // that each would have produced alone (§1's independence property).
            TimelineStatus status = process_timeline_until_horizon(load.instance.mwpm, load.dets, decoder->horizon);
            uint64_t solve_ended = tick();

            // ---- `extract_c`, on a completed core only, and **the extraction alone**: the
            // observable mask and the weight, which are the decoder's answer for this core.
            //
            // `ResetPolicy::CALLER_RESETS` keeps `reset_for_next_shot` out of the region. That reset
            // prepares the instance for the *next* shot and nothing about this shot's answer depends
            // on it, so a pipelined decoder that double-buffered its instances would not pay it on
            // this critical path. It still runs, in the same place in the same order, immediately
            // below — it is simply on the far side of the second clock read.
            uint64_t extract_ended = solve_ended;
            if (status == TimelineStatus::COMPLETE) {
                load.part = load.harvester.extract_only_to_obs(
                    load.instance.mwpm, load.dets, ResetPolicy::CALLER_RESETS);
                extract_ended = tick();
                // The debt `CALLER_RESETS` incurs, paid before anything else touches the instance.
                reset_for_next_shot(load.instance.mwpm);
            } else {
                // Outside every region, as the design asks: `abandon_shot` is `Mwpm::reset`, which
                // sweeps the instance's entire node vector, and these instances are sized for the
                // full `H`.
                abandon_shot(load.instance.mwpm);
                load.part.committed = pm::MatchingResult();
                load.truncated = true;
            }

            load.build_ticks = build_ended - build_started;
            load.solve_ticks = solve_ended - build_ended;
            load.extract_ticks = extract_ended - solve_ended;
            uint64_t total = load.build_ticks + load.solve_ticks + load.extract_ticks;
            if (total > out.core_max) {
                out.core_max = total;
                out.crit_core = (uint32_t)c;
            }
        }

        // ---- `combine`. A `k`-fold reduction, unconditionally over all `k`: a truncated or idle
        // core's partial was zeroed above, so it contributes the identity to both halves.
        uint64_t combine_started = tick();
        pm::obs_int obs_mask = 0;
        pm::total_weight_int weight = 0;
        bool escalated = false;
        for (size_t c = 0; c < k; c++) {
            obs_mask ^= core[c].part.committed.obs_mask;
            weight += core[c].part.committed.weight;
            escalated |= core[c].truncated;
        }
        uint64_t combine_ended = tick();

        // ---- `fallback`. The stock decode on `G` from the **full** syndrome, on its own instance,
        // exactly the region the previous latency profiler timed. Run on every shot, not only on
        // escalating ones, so that `min(fallback, sparse_k)` is available on every shot.
        std::fill(obs_bytes.begin(), obs_bytes.end(), (uint8_t)0);
        pm::total_weight_int fallback_weight = 0;
        uint64_t fallback_started = tick();
        pm::decode_detection_events(*fallback_instance, shot, obs_bytes.data(), fallback_weight, false);
        uint64_t fallback_ended = tick();

        out.n_def = n;
        out.n_comps = bal.n_comps;
        out.largest = bal.largest;
        out.escalated = escalated ? 1 : 0;
        out.crit_size = (uint32_t)core[out.crit_core].slice_nodes.size();
        out.crit_ncomps = core[out.crit_core].ncomps;
        out.uf = uf_ended - uf_started;
        out.balance = balance_ended - balance_started;
        out.build_crit = core[out.crit_core].build_ticks;
        out.solve_crit = core[out.crit_core].solve_ticks;
        out.extract_crit = core[out.crit_core].extract_ticks;
        out.combine = combine_ended - combine_started;
        out.sparse_k = out.uf + out.balance + out.core_max + out.combine;
        out.fallback = fallback_ended - fallback_started;
        out.system = out.escalated ? out.fallback : std::min(out.fallback, out.sparse_k);
        out.seed = seed_ended - seed_started;
        out.hbuild = hbuild_ended - seed_ended;
        out.hbuild_intersect = (uint64_t)hbuild_timing.intersect_ns;
        out.hbuild_finalize = (uint64_t)hbuild_timing.finalize_ns;
        out.pipeline = out.seed + out.hbuild + out.sparse_k;
        out.system_pipeline = out.escalated ? out.fallback : std::min(out.fallback, out.pipeline);
        out.obs_mask = obs_mask;
        out.weight = weight;
    }
};

/// §6. The monolithic per-component path on the same shot, untimed, against the `k`-core one.
///
/// Deviation from the design's §6, recorded rather than silently taken: the combined `obs` and
/// `weight` are compared **on non-escalating shots only**. On an escalating shot the two paths
/// legitimately differ — a core whose timeline truncates abandons the components it *did* complete,
/// while the monolithic path extracts each completed component individually — and both discard
/// Phase 1 in full, so the difference is not observable in the decoder's output and is not a bug.
/// The escalation flag, which is what the exactness claim rests on, is compared on every shot.
void verify_shot(BallDecoder& decoder, const std::vector<uint64_t>& shot, uint64_t shot_index, const ShotMeasurement& m) {
    Phase1Outcome reference = decoder.decode_phase1_production(shot);
    bool reference_escalated = reference.status == TimelineStatus::TRUNCATED;
    if (reference_escalated != (m.escalated != 0)) {
        throw std::runtime_error(
            "verify: escalation disagrees with the monolithic per-component path at shot " +
            std::to_string(shot_index));
    }
    if (m.escalated)
        return;
    if (reference.harvest.committed.obs_mask != m.obs_mask) {
        throw std::runtime_error(
            "verify: combined observable mask disagrees with the monolithic per-component path at shot " +
            std::to_string(shot_index));
    }
    if (reference.harvest.committed.weight != m.weight) {
        throw std::runtime_error(
            "verify: combined weight disagrees with the monolithic per-component path at shot " +
            std::to_string(shot_index));
    }
}

/// Total capacity across every buffer a measured shot could make grow: the cores' own arenas and
/// instances, and the `H` builder's arena. Invariant 11 as this tree states it — if this number
/// moves after the warm-up, a measured shot allocated.
uint64_t capacity_signature(const CellRunner& runner) {
    uint64_t total = runner.decoder->arena.grow_events;
    for (const SolverCore& core : *runner.cores)
        total += core.capacity_signature();
    return total;
}

/// One `(d, p, T, k)` cell: warm-up, then `--shots` measured shots.
///
/// `record` says whether the pass fills `arena` and `aggregate`; the `--check-determinism` second
/// pass runs the identical path with it off, so the two digests compare like with like.
Digest run_cell(
    const Options& options,
    CellRunner& runner,
    ShotSampler& sampler,
    RowArena& arena,
    CellAggregate& aggregate,
    Digest* escalation_digest,
    uint64_t* capacity_after_warmup,
    bool record,
    bool interactive,
    const char* label,
    size_t cell_index,
    size_t total_cells) {
    Digest digest;
    Digest escalation;

    sampler.rng = std::mt19937_64(options.seed);
    size_t done = 0;
    size_t cell_total = options.warmup + options.shots;
    ShotMeasurement measurement;

    // Read before the first shot so that `--warmup 0` compares against the as-configured state
    // rather than against zero.
    if (capacity_after_warmup != nullptr)
        *capacity_after_warmup = capacity_signature(runner);

    draw_progress(interactive, label, cell_index, total_cells, done, cell_total);
    while (done < cell_total) {
        size_t batch = std::min(SHOT_BATCH, cell_total - done);
        sampler.sample(batch);
        for (size_t i = 0; i < batch; i++) {
            size_t index = done + i;
            bool measured = index >= options.warmup;
            // Warm-up shots run the same workload and are not profiled.
            if (measured)
                runner.run_shot<true>(sampler.shots[i], measurement);
            else
                runner.run_shot<false>(sampler.shots[i], measurement);
            if (!measured) {
                if (capacity_after_warmup != nullptr && index + 1 == options.warmup)
                    *capacity_after_warmup = capacity_signature(runner);
                continue;
            }
            uint64_t shot_index = (uint64_t)(index - options.warmup);
            if (options.verify)
                verify_shot(*runner.decoder, sampler.shots[i], shot_index, measurement);

            // The deterministic set: everything that is a function of `(seed, d, p, T, k, alpha)`
            // and of nothing else. The `crit_*` columns are excluded and the whole LPT assignment is
            // included in their place; see the header.
            digest.add(shot_index);
            digest.add(measurement.n_def);
            digest.add(measurement.n_comps);
            digest.add(measurement.largest);
            digest.add(measurement.escalated);
            digest.add(measurement.pred_max);
            digest.add(measurement.pred_sum);
            digest.add(measurement.imbalance);
            for (size_t c = 0; c < runner.k; c++) {
                const SolverCore& load = (*runner.cores)[c];
                uint32_t size = (uint32_t)load.slice_nodes.size();
                digest.add(size);
                digest.add(load.ncomps);
                digest.add(load.pred_cost);
            }
            escalation.add(shot_index);
            escalation.add(measurement.escalated);

            if (!record)
                continue;

            Row row;
            row.shot = shot_index;
            row.n_def = measurement.n_def;
            row.n_comps = measurement.n_comps;
            row.largest = measurement.largest;
            row.escalated = measurement.escalated;
            row.pred_max = measurement.pred_max;
            row.pred_sum = measurement.pred_sum;
            row.imbalance = measurement.imbalance;
            row.crit_core = measurement.crit_core;
            row.crit_size = measurement.crit_size;
            row.crit_ncomps = measurement.crit_ncomps;
            row.uf = measurement.uf;
            row.balance = measurement.balance;
            row.build_crit = measurement.build_crit;
            row.solve_crit = measurement.solve_crit;
            row.extract_crit = measurement.extract_crit;
            row.core_max = measurement.core_max;
            row.combine = measurement.combine;
            row.sparse_k = measurement.sparse_k;
            row.fallback = measurement.fallback;
            row.system = measurement.system;
            row.seed = measurement.seed;
            row.hbuild = measurement.hbuild;
            row.hbuild_intersect = measurement.hbuild_intersect;
            row.hbuild_finalize = measurement.hbuild_finalize;
            row.pipeline = measurement.pipeline;
            row.system_pipeline = measurement.system_pipeline;
            arena.append(row);

            aggregate.shots++;
            aggregate.escalations += measurement.escalated;
            aggregate.fallback.add(measurement.fallback);
            aggregate.sparse_k.add(measurement.sparse_k);
            aggregate.system.add(measurement.system);
            aggregate.hbuild.add(measurement.hbuild);
            aggregate.pipeline.add(measurement.pipeline);
            aggregate.system_pipeline.add(measurement.system_pipeline);
            aggregate.seed_sum += (__int128)measurement.seed;
            aggregate.hbuild_intersect_sum += (__int128)measurement.hbuild_intersect;
            aggregate.hbuild_finalize_sum += (__int128)measurement.hbuild_finalize;
            aggregate.uf_sum += (__int128)measurement.uf;
            aggregate.balance_sum += (__int128)measurement.balance;
            aggregate.build_crit_sum += (__int128)measurement.build_crit;
            aggregate.solve_crit_sum += (__int128)measurement.solve_crit;
            aggregate.extract_crit_sum += (__int128)measurement.extract_crit;
            aggregate.combine_sum += (__int128)measurement.combine;
            aggregate.imbalance_sum += measurement.imbalance;
            if (measurement.crit_ncomps < aggregate.crit_ncomps_hist.size())
                aggregate.crit_ncomps_hist[measurement.crit_ncomps]++;
            if (measurement.largest < aggregate.largest_hist.size())
                aggregate.largest_hist[measurement.largest]++;
        }
        done += batch;
        // Between batches, from a between-shot position: no timed region is open, so neither the
        // flush nor the progress line lands in any row.
        if (record && arena.rows.size() >= ROW_ARENA_FLUSH_AT)
            arena.flush();
        draw_progress(interactive, label, cell_index, total_cells, done, cell_total);
    }
    if (record)
        arena.flush();
    if (interactive)
        std::fprintf(stderr, "\r\x1b[K");
    if (escalation_digest != nullptr)
        *escalation_digest = escalation;
    return digest;
}

void write_series(std::ofstream& out, const char* name, const TickSeries& series, const char* trailing) {
    std::vector<long long> samples = series.samples;
    out << "      \"" << name << "\": {\"mean\": " << fmt_g(series.mean()) << ", \"max\": " << series.max
        << ", \"p50\": " << fmt_g(percentile(samples, 0.50)) << ", \"p99\": " << fmt_g(percentile(samples, 0.99))
        << ", \"p99_9\": " << fmt_g(percentile(samples, 0.999)) << ", \"sum\": \"" << fmt_i128(series.sum)
        << "\", \"sum_sq\": \"" << fmt_i128(series.sum_sq) << "\"}" << trailing << "\n";
}

void write_histogram(std::ofstream& out, const char* name, const std::vector<uint64_t>& histogram, const char* trailing) {
    out << "      \"" << name << "\": {";
    bool first = true;
    for (size_t value = 0; value < histogram.size(); value++) {
        if (histogram[value] == 0)
            continue;
        if (!first)
            out << ", ";
        out << "\"" << value << "\": " << histogram[value];
        first = false;
    }
    out << "}" << trailing << "\n";
}

double mean_of(__int128 sum, uint64_t count) {
    return count == 0 ? 0.0 : (double)sum / (double)count;
}

void write_aggregate(std::ofstream& out, const Options& options, const std::vector<CellAggregate>& cells) {
    out << "{\n";
    out << "  \"binary\": \"load_balancer_profiler\",\n";
    out << "  \"git_hash\": \"" << SPECMATCHING_GIT_HASH << "\",\n";
    out << "  \"stim_version\": \"" << SPECMATCHING_STIM_VERSION << "\",\n";
    out << "  \"model\": \"k_core_critical_path_model\",\n";
    out << "  \"tick_unit\": \"" << hires_timer_name() << ", uncorrected; see run.log timer_overhead_ns\",\n";
    out << "  \"seed\": " << options.seed << ",\n";
    out << "  \"alpha\": " << fmt_g(options.alpha) << ",\n";
    out << "  \"verify\": " << (options.verify ? "true" : "false") << ",\n";
    out << "  \"cells\": [\n";
    for (size_t i = 0; i < cells.size(); i++) {
        const CellAggregate& cell = cells[i];
        std::ostringstream key;
        key << "d=" << cell.distance << ",p=" << fmt_g(cell.noise) << ",T=" << fmt_g(cell.horizon)
            << ",k=" << cell.cores;
        double fallback_mean = cell.fallback.mean();
        double system_mean = cell.system.mean();
        out << "    {\n";
        out << "      \"cell\": \"" << key.str() << "\",\n";
        out << "      \"d\": " << cell.distance << ", \"p\": " << fmt_g(cell.noise)
            << ", \"T\": " << fmt_g(cell.horizon) << ", \"k\": " << cell.cores << ", \"alpha\": "
            << fmt_g(cell.alpha) << ",\n";
        out << "      \"shots\": " << cell.shots << ", \"escalations\": " << cell.escalations << ",\n";
        write_series(out, "fallback", cell.fallback, ",");
        write_series(out, "sparse_k", cell.sparse_k, ",");
        write_series(out, "system", cell.system, ",");
        write_series(out, "hbuild", cell.hbuild, ",");
        write_series(out, "pipeline", cell.pipeline, ",");
        write_series(out, "system_pipeline", cell.system_pipeline, ",");
        // The input's regions, kept out of `module_means` so that those six still sum to `sparse_k`.
        out << "      \"input_means\": {\"seed\": " << fmt_g(mean_of(cell.seed_sum, cell.shots))
            << ", \"hbuild\": " << fmt_g(cell.hbuild.mean())
            << ", \"hbuild_intersect\": " << fmt_g(mean_of(cell.hbuild_intersect_sum, cell.shots))
            << ", \"hbuild_finalize\": " << fmt_g(mean_of(cell.hbuild_finalize_sum, cell.shots)) << "},\n";
        out << "      \"module_means\": {\"uf\": " << fmt_g(mean_of(cell.uf_sum, cell.shots))
            << ", \"balance\": " << fmt_g(mean_of(cell.balance_sum, cell.shots))
            << ", \"build_crit\": " << fmt_g(mean_of(cell.build_crit_sum, cell.shots))
            << ", \"solve_crit\": " << fmt_g(mean_of(cell.solve_crit_sum, cell.shots))
            << ", \"extract_crit\": " << fmt_g(mean_of(cell.extract_crit_sum, cell.shots))
            << ", \"combine\": " << fmt_g(mean_of(cell.combine_sum, cell.shots)) << "},\n";
        double imbalance_mean = cell.shots == 0 ? 0.0 : cell.imbalance_sum / (double)cell.shots;
        out << "      \"imbalance_mean\": " << fmt_g(imbalance_mean) << ",\n";
        // The design lists "mean imbalance" and "mean pred_max / pred_sum * k" separately; by §3.2's
        // definition of `imbalance` they are the same quantity, so the second key is emitted as an
        // alias rather than as a second, subtly different number.
        out << "      \"pred_max_over_pred_sum_times_k_mean\": " << fmt_g(imbalance_mean)
            << ",\n";
        out << "      \"speedup_mean\": " << fmt_g(system_mean > 0 ? fallback_mean / system_mean : 0.0) << ",\n";
        double system_pipeline_mean = cell.system_pipeline.mean();
        out << "      \"speedup_pipeline_mean\": "
            << fmt_g(system_pipeline_mean > 0 ? fallback_mean / system_pipeline_mean : 0.0) << ",\n";
        write_histogram(out, "crit_ncomps_hist", cell.crit_ncomps_hist, ",");
        write_histogram(out, "largest_hist", cell.largest_hist, ",");
        out << "      \"raw_truncated_at_shot\": " << cell.raw_truncated_at_shot << ",\n";
        out << "      \"determinism_check\": " << cell.determinism_check << ",\n";
        out << "      \"escalation_matches_first_k\": " << cell.escalation_matches_first_k << "\n";
        out << "    }" << (i + 1 == cells.size() ? "" : ",") << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

/// One `(d, p)` of the grid: circuit, DEM, ball tables, `k_max` solver instances and the fallback
/// instance, then the `T` x `k` sweep over the same shots.
void run_dp_cell(
    const Options& options,
    size_t distance,
    double noise,
    RowArena& arena,
    std::vector<CellAggregate>& aggregates,
    std::ofstream& log,
    bool interactive,
    size_t& cell_index,
    size_t total_cells) {
    char setup_label[160];
    std::snprintf(setup_label, sizeof(setup_label), "d=%zu p=%g", distance, noise);
    draw_setup(interactive, setup_label, cell_index + 1, total_cells);

    // ---- §5 step 1. Circuit -> DEM, through the corpus's own generator call.
    ShotSampler sampler = ShotSampler::make(distance, distance, noise, options.seed);
    log << "\ncorpus d=" << distance << " rounds=" << distance << " p=" << noise << "\n";
    log << "  generator_call=" << generator_call(distance, distance, noise) << "\n";

    // ---- §5 step 2. Ball tables once, at `T_max = max(T list)` and `R = 2 * T_max`.
    pm::Mwpm probe = pm::detector_error_model_to_mwpm(sampler.dem, NUM_DISTINCT_WEIGHTS, false);
    double unit = edge_weight_units(probe.flooder.graph);
    double max_horizon = *std::max_element(options.horizons.begin(), options.horizons.end());

    BallConfig ball_config;
    ball_config.T = max_horizon * unit;
    ball_config.ball.T_max = max_horizon * unit;
    ball_config.ball.R = 2.0 * ball_config.ball.T_max;
    ball_config.collect_component_stats = false;
    ball_config.collect_harvest_diagnostics = false;
    ball_config.collect_structural_counters = false;
    ball_config.verify_component_decomposition = false;
    ball_config.verify_against_g = false;
    ball_config.ball.validate();
    if (!(ball_config.ball.R >= 2.0 * ball_config.ball.T_max))
        throw std::runtime_error("the compiled ball radius is below 2 * T_max");

    BallDecoder decoder = BallDecoder::from_detector_error_model(sampler.dem, ball_config, NUM_DISTINCT_WEIGHTS);
    decoder.harvester.collect_diagnostics = false;
    decoder.harvester.use_legacy_enumeration = false;

    // The fallback core's own instance, outside the `k`. Separate from `decoder.g_mwpm`, which the
    // §M2.1 preamble and the `--verify` path use, so that neither is warmed by the other.
    pm::Mwpm fallback_instance = pm::detector_error_model_to_mwpm(sampler.dem, NUM_DISTINCT_WEIGHTS, false);

    const pm::MatchingGraph& g_graph = decoder.g_mwpm.flooder.graph;
    size_t n_max = g_graph.nodes.size();
    size_t num_observables = g_graph.num_observables;
    size_t k_max = options.core_counts.back();

    log << "  dem_has_negative_weights=" << (decoder.dem_has_negative_weights ? 1 : 0)
        << " edge_weight_unit=" << unit << " T_max=" << ball_config.ball.T_max << " R=" << ball_config.ball.R
        << " n_max=" << n_max << " num_observables=" << num_observables
        << " ball_graph_build_mode=" << (decoder.config.mode == BallGraphBuildMode::BITSET ? "BITSET" : "SCAN") << "\n";

    // ---- §5 step 4. Everything allocated here, before the first measured shot.
    Balancer balancer;
    balancer.configure(n_max, options.alpha);

    std::vector<SolverCore> cores(k_max);
    uint64_t instance_memory = 0;
    {
        // §3.3: each of the `k` instances is sized for the **full** `H` — the worst case is one
        // giant component landing on one core — and is built once, here. `rebuild` sizes the node
        // pool from the graph it is given, so a node-only graph of `n_max` nodes is what sizes it.
        BallGraph full;
        full.h_to_det.resize(n_max, 0);
        for (SolverCore& core : cores) {
            core.instance.configure(num_observables, g_graph.normalising_constant);
            core.harvester.collect_diagnostics = false;
            core.harvester.use_legacy_enumeration = false;
            core.reserve(n_max);
            core.instance.rebuild(decoder.tables, full, decoder.arena, nullptr);
            // The sizing rebuild leaves `used_nodes = n_max`, which would make every later rebuild
            // clear the whole pool instead of the nodes the previous load actually used — an
            // `O(n_max)` term inside `build_c` that no real core pays. The pool keeps its size;
            // only the "how much is stale" cursor is reset. Nodes above a load's own `n_c` are
            // never reachable from it, so the adjacency left in them is inert.
            core.instance.used_nodes = 0;
            core.instance.grow_events = 0;
            instance_memory += instance_bytes(core.instance);
        }
    }
    log << "  solver_instances=" << k_max << " node_capacity_each=" << n_max
        << " instance_memory_total_bytes=" << instance_memory
        << " (node pools and adjacency capacities, summed over the k instances, right after sizing)\n";
    log.flush();

    CellRunner runner;
    runner.decoder = &decoder;
    runner.fallback_instance = &fallback_instance;
    runner.balancer = &balancer;
    runner.cores = &cores;
    runner.configure(num_observables, k_max);

    for (double horizon : options.horizons) {
        decoder.set_horizon(horizon * unit);
        // "Identical escalating set across `k`" (§2, §8) is checked against the first `k` of this
        // `(d, p, T)`, under `--check-determinism`.
        Digest first_k_escalation;
        bool have_first_k = false;

        for (size_t k : options.core_counts) {
            cell_index++;
            char label[160];
            std::snprintf(label, sizeof(label), "d=%zu p=%g T=%g k=%zu", distance, noise, horizon, k);

            runner.k = k;
            arena.labels = CellLabels{distance, noise, horizon, k, options.alpha};
            arena.rows_this_cell = 0;
            arena.truncated_at_shot = -1;

            aggregates.emplace_back();
            CellAggregate& aggregate = aggregates.back();
            aggregate.distance = distance;
            aggregate.noise = noise;
            aggregate.horizon = horizon;
            aggregate.cores = k;
            aggregate.alpha = options.alpha;
            aggregate.reset(n_max + 1);
            aggregate.fallback.reserve(options.shots);
            aggregate.sparse_k.reserve(options.shots);
            aggregate.system.reserve(options.shots);
            aggregate.hbuild.reserve(options.shots);
            aggregate.pipeline.reserve(options.shots);
            aggregate.system_pipeline.reserve(options.shots);

            log << "  cell T=" << horizon << " k=" << k << " T_weight_units=" << horizon * unit
                << " T_int=" << (int64_t)decoder.horizon << " started=" << stamp() << "\n";
            log.flush();

            Digest escalation;
            uint64_t capacity_at_warmup_end = 0;
            Digest digest = run_cell(
                options, runner, sampler, arena, aggregate, &escalation, &capacity_at_warmup_end, true, interactive,
                label, cell_index, total_cells);
            uint64_t capacity_at_end = capacity_signature(runner);

            aggregate.raw_truncated_at_shot = arena.truncated_at_shot;

            if (options.check_determinism) {
                Digest replay = run_cell(
                    options, runner, sampler, arena, aggregate, nullptr, nullptr, false, interactive, label, cell_index,
                    total_cells);
                aggregate.determinism_check = replay.value == digest.value ? 1 : 0;
                if (have_first_k)
                    aggregate.escalation_matches_first_k = escalation.value == first_k_escalation.value ? 1 : 0;
                else
                    aggregate.escalation_matches_first_k = 1;
            }
            if (!have_first_k) {
                first_k_escalation = escalation;
                have_first_k = true;
            }

            log << "  cell T=" << horizon << " k=" << k << " finished=" << stamp() << " rows=" << arena.rows_this_cell
                << " escalations=" << aggregate.escalations
                // Invariant 11: after the warm-up no buffer grew, so no measured shot allocated.
                // Reported rather than asserted, because a cell whose warm-up was too short for its
                // corpus should say so in its own artifact.
                << " no_growth_after_warmup=" << (capacity_at_end == capacity_at_warmup_end ? 1 : 0);
            if (options.check_determinism) {
                log << " determinism_check=" << (aggregate.determinism_check == 1 ? "pass" : "FAIL")
                    << " escalation_matches_first_k="
                    << (aggregate.escalation_matches_first_k == 1 ? "pass" : "FAIL");
            }
            if (arena.truncated_at_shot >= 0)
                log << " raw_truncated_at_shot=" << arena.truncated_at_shot;
            log << "\n";
            log.flush();

            if (options.check_determinism && aggregate.determinism_check != 1)
                throw std::runtime_error("--check-determinism: the non-tick columns differ between two runs of " +
                                         std::string(label));
            if (options.check_determinism && aggregate.escalation_matches_first_k != 1)
                throw std::runtime_error("--check-determinism: the escalating set changed with k at " +
                                         std::string(label));
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

    std::ofstream shots_csv(options.out_dir + "/shots.csv");
    std::ofstream log(options.out_dir + "/run.log");
    if (!shots_csv.is_open() || !log.is_open()) {
        std::cerr << "error: could not open the output files under " << options.out_dir << "\n";
        return 1;
    }
    shots_csv << "d,p,T,k,alpha,shot,n_def,n_comps,largest,escalated,"
                 "pred_max,pred_sum,imbalance,crit_core,crit_size,crit_ncomps,"
                 "ticks_uf,ticks_balance,ticks_build_crit,ticks_solve_crit,ticks_extract_crit,"
                 "ticks_core_max,ticks_combine,ticks_sparse_k,ticks_fallback,ticks_system,"
                 "ticks_seed,ticks_hbuild,ticks_hbuild_intersect,ticks_hbuild_finalize,"
                 "ticks_pipeline,ticks_system_pipeline\n";

    RowArena arena;
    arena.configure(shots_csv);

    // ---- The startup calibration, before any decoding, so it measures the clock rather than the
    // clock plus whatever the decoder left in the caches.
    double overhead_ns = hires_timer_overhead_ns();
    bool thread_scoped = hires_timer_is_thread_scoped();

    log << "load_balancer_profiler - k-core critical-path latency of the decomposed H solve\n";
    log << "git_hash=" << SPECMATCHING_GIT_HASH << "\n";
    log << "build_flags=" << SPECMATCHING_BUILD_FLAGS << "\n";
    log << "stim_version=" << SPECMATCHING_STIM_VERSION << "\n";
#if defined(__APPLE__)
    log << "platform=macos\n";
#else
    log << "platform=linux\n";
#endif
    log << "hires_timer_name=" << hires_timer_name() << "\n";
    log << "timer_is_thread_scoped=" << (thread_scoped ? 1 : 0) << "\n";
    log << "timer_overhead_ns=" << overhead_ns << " (10^6 back-to-back hires_now_ns() calls)\n";
    log << "rows_are_uncorrected=1 (every ticks_* column is a raw delta of this clock; the analysis"
           " subtracts timer_overhead_ns, and one tick is one nanosecond as this backend reports it)\n";
    log << "threads=1 (no thread is ever spawned)\n";
    log << "k_core_critical_path_model=1: the k loads are solved one after another on the profiling"
           " thread, each on its own pre-allocated Mwpm instance, and the per-core figure is the MAX"
           " over loads, never the sum. Assumption 1: cores share nothing and start together, and no"
           " cost is modelled for moving a core's edge slice to it - the balance region does the"
           " scatter once, serially, on the profiling thread, which is conservative. Assumption 2:"
           " the fallback core and the k solver cores start at the same instant, when the edge list"
           " is ready.\n";
    log << "warm_cache_caveat=1: a load solved immediately after the previous load sees a warmer"
           " cache than a real core would; this travels with every number in shots.csv and"
           " agg.json.\n";
    log << "timed_region_input=compute_seeded_detection_events + build_ball_graph are the INPUT of the"
           " k-core model and are in no term of sparse_k or system; they are timed on the CPU into their"
           " own columns (see timed_region_seed, timed_region_hbuild)\n";
    log << "timed_region_seed=compute_seeded_detection_events: the negative-weight preamble on G and"
           " the sorted detection-event list\n";
    log << "timed_region_hbuild=build_ball_graph end to end, on the CPU. hbuild_intersect and"
           " hbuild_finalize are the split build_ball_graph measures itself (BallGraphTiming, same"
           " thread-scoped clock): the ball intersection, then the canonical edge sort and arena"
           " cleanup. They do not sum to hbuild; the gap is the call's setup and the inner timer's"
           " three reads, which hbuild includes\n";
    log << "pipeline=seed + hbuild + sparse_k; system_pipeline=escalated ? fallback :"
           " min(fallback, pipeline). sparse_k and system are unchanged and still exclude the input\n";
    log << "timed_region_uf=union-find over H's defect-defect edges, component id per defect,"
           " component sizes\n";
    log << "timed_region_balance=cost lookup, LPT assignment, and the scatter of defects, edges and"
           " boundary edges into the per-core slices; the balancer statistics (pred_cost, pred_max,"
           " pred_sum, imbalance, sizes, n_comps) are read AFTER the region closes\n";
    log << "timed_region_build_c=sub-H_c construction from the core's slice (nodes relabelled in"
           " ascending full-H order, edges and boundary edges copied with H's weights and cutoffs),"
           " then BallMwpm::rebuild (which is where H's observable masks are written), then the"
           " detection-event fill dets_c = 0..n_c-1\n";
    log << "timed_region_solve_c=process_timeline_until_horizon(instance_c, dets_c, T_int) ONLY\n";
    log << "timed_region_extract_c=extract_only_to_obs(instance_c, dets_c, CALLER_RESETS) on a"
           " COMPLETE core only - the EXTRACTION ALONE (the observable mask and the weight)."
           " reset_for_next_shot is OUTSIDE this region: it prepares the instance for the next shot"
           " and this shot's answer does not depend on it, so a pipelined decoder that"
           " double-buffered its instances would not pay it on this critical path. It still runs, in"
           " the same place in the same order, immediately after the region closes.\n";
    log << "timed_region_combine=XOR of the k partial observable masks, sum of the k partial"
           " weights, OR of the k truncation flags\n";
    log << "timed_region_fallback=pm::decode_detection_events on G from the FULL syndrome, on its"
           " own instance, run on every shot (not only escalating ones)\n";
    log << "untimed=reset_for_next_shot after extraction on a COMPLETE core; abandon_shot"
           " (Mwpm::reset) on a TRUNCATED core; the fallback's obs buffer fill and weight reset; the"
           " --verify cross-check\n";
    log << "decoder_change=Harvester::extract_only_to_obs / extract_only_to_match_edges take a"
           " ResetPolicy, defaulting to RESET_BEFORE_RETURN (the previous behaviour, so every"
           " production call site is unchanged). This binary passes CALLER_RESETS and runs"
           " reset_for_next_shot itself, outside the timed region. Nothing else under"
           " src/specmatching/ is edited; the solve, the extraction and the escalating set are"
           " unchanged.\n";
    log << "empty_core_policy=a core with no components assigned does no work and opens no timer;"
           " its build/solve/extract ticks are 0 and it contributes nothing to the makespan\n";
    log << "timer_reads_per_measured_shot=2 (uf) + 2 (balance) + 4 per non-empty core + 2 (combine)"
           " + 2 (fallback); four per core because consecutive regions share their boundary read,"
           " and three on a TRUNCATED core, which has no extraction\n";
    log << "warmup_is_not_profiled=1: a warm-up shot runs the identical workload - same builds,"
           " solves, extractions and fallback, in the same order, on the same buffers - and reads no"
           " clock at all. run_shot is templated on the switch, so the reads compile out of the"
           " warm-up instantiation rather than being taken and thrown away.\n";
    log << "sparse_k=uf + balance + max_c(build_c + solve_c + extract_c) + combine\n";
    log << "system=escalated ? fallback : min(fallback, sparse_k); nothing else in this artifact is"
           " called system\n";
    log << "escalation_predicate=any core's timeline is TRUNCATED, which equals any component"
           " TRUNCATED because a core's timeline is the disjoint union of its components'\n";
    log << "cost_model=cost(s) = s^alpha, a function of component SIZE ONLY - the status of a"
           " component is not known before it is solved. cost_tab[s] is filled once per (d, p) with"
           " pow; there is no pow call per shot.\n";
    log << "alpha=" << options.alpha << "\n";
    log << "balancer=LPT (longest-processing-time-first): components in descending size, ties by"
           " ascending minimum defect id, each to the core with the smallest accumulated cost, ties"
           " to the lowest core index. Counting sort on size into a fixed arena; least-loaded core"
           " by linear scan over k.\n";
    log << "crit_core=argmax_c (build_c + solve_c + extract_c), ties to the lowest index - a"
           " function of the MEASURED ticks, so crit_core, crit_size and crit_ncomps are NOT in the"
           " deterministic set. --check-determinism diffs n_def, n_comps, largest, escalated,"
           " pred_max, pred_sum, imbalance and the whole LPT assignment (sizes[c], n_comps[c],"
           " pred_cost[c] for every core) in their place.\n";
    log << "verify=" << (options.verify ? 1 : 0)
        << " (per shot, the monolithic per-component path untimed; the escalation flag is compared on"
           " every shot, the combined obs mask and weight on non-escalating shots only - on an"
           " escalating shot a core abandons the components it did complete while the monolithic path"
           " extracts each of them, and both paths discard Phase 1 in full)\n";
    log << "check_determinism=" << (options.check_determinism ? 1 : 0) << "\n";
    log << "seed=" << options.seed << " (shot i of (d, p) is identical across the T and k sweeps)\n";
    log << "shots_per_cell=" << options.shots << " (after warm-up)\n";
    log << "warmup_shots_per_cell=" << options.warmup << " (run and discarded)\n";
    log << "T_unit=one lattice edge weight (the median discretised edge of G, in DEM float units);"
           " each cell below records T, T_weight_units and T_int\n";
    log << "decoder_otherwise_unchanged=1 (apart from the ResetPolicy parameter recorded above,"
           " nothing under src/specmatching/ is edited; the vendored solver, escalate.{h,cc} and the"
           " sub-H construction are used as they are)\n";

    if (!thread_scoped) {
        log << "WARNING: no thread-scoped clock was available, so every tick on this run includes"
               " time the thread spent off the CPU and is an upper bound only.\n";
        std::fprintf(
            stderr,
            "\n  WARNING: latencies are being measured with a WALL CLOCK, backend '%s'.\n"
            "  Every tick includes time the thread spent off the CPU and is an upper bound only.\n\n",
            hires_timer_name());
        std::fflush(stderr);
    }
    log.flush();

    bool interactive = stderr_is_terminal();
    size_t total_cells = options.distances.size() * options.error_rates.size() * options.horizons.size() *
                         options.core_counts.size();
    size_t cell_index = 0;
    bool ok = true;
    std::vector<CellAggregate> aggregates;
    aggregates.reserve(total_cells);

    for (size_t distance : options.distances) {
        if (!ok)
            break;
        for (double noise : options.error_rates) {
            try {
                run_dp_cell(options, distance, noise, arena, aggregates, log, interactive, cell_index, total_cells);
            } catch (const std::exception& error) {
                std::fprintf(stderr, "\n");
                std::cerr << "error at d=" << distance << " p=" << noise << ": " << error.what() << "\n";
                log << "  ERROR=" << error.what() << "\n";
                ok = false;
                break;
            }
        }
    }

    {
        std::ofstream agg(options.out_dir + "/agg.json");
        if (agg.is_open()) {
            write_aggregate(agg, options, aggregates);
        } else {
            std::cerr << "warning: could not open " << options.out_dir << "/agg.json for writing\n";
            log << "WARNING: could not open agg.json for writing\n";
        }
    }

    log << "\nfinished=" << stamp() << "\n";
    log << "rows_total=" << arena.rows_total << "\n";
    log << (ok ? "status=ok\n" : "status=incomplete (see ERROR above)\n");

    std::printf("wrote %s/{shots.csv,agg.json,run.log}\n", options.out_dir.c_str());
    return ok ? 0 : 1;
}
