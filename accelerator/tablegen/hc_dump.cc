// Table compiler and vector generator for the h_construct accelerator.
//
//     hc_dump <out-dir> [--d 5] [--rounds 5] [--p 0.005] [--shots 20] [--keep 0] [--seed 7] [--lanes 64]
//                       [--time-coord -1] [--no-prune] [--no-dedup] [--stats-only]
//
// Compiles the software ball tables (`compile_ball_tables`, T = T_max = 2 lattice edge weights,
// R = 2T) for a generated rotated-surface-code memory experiment and re-encodes them into the two
// tables the accelerator reads:
//
//   * ball table (`ball.hex`): one fixed-width word per *unique* ball row. An entry is
//     `{w, delta}` with `delta = target - source`, so a row no longer names its source and two
//     detectors whose balls are time shifts of each other hold the same row. Rows are shared only
//     between detectors with the same space coordinate (every detector coordinate except the time
//     one): two detectors at different space positions never share a row, even when their
//     `{w, delta}` lists happen to coincide, so the hardware stays reprogrammable for DEMs whose
//     space-like weights vary within a layer. The lists are compared exactly, so nothing about the
//     DEM is assumed: a DEM that is not time invariant just yields more rows.
//   * boundary table (`bnd.hex`): one word per detector, `{row, hb, bw, size}` -- the row the
//     detector's ball lives in, the boundary edge (`hb` set when `bcost <= T`), and the row length.
//
// With pruning on (the default), entries with `target <= source` are dropped offline: H keeps each
// pair once, from its lower endpoint, and that test depends only on the row. The lanes still test
// the sign of `delta`, so an unpruned table (`--no-prune`) is run by the same RTL, and padding
// entries are written with `delta = -1` so they never pass.
//
// `--shots` is the size of the sampled batch and `--keep N` (N > 0) keeps only its first N shots,
// which is how the checked-in vectors were drawn (the first 20 of a 60-shot batch).
//
// Outputs, all in <out-dir>: ball.hex, bnd.hex, params.vh, syn_<k>.txt (`M` then the `M` ascending
// detector ids), ref_<k>.txt (the software H from `build_ball_graph(..., SCAN)` on the uncompressed
// tables, as `E u v w` / `B u w`), and stats.txt (one `key=value` line, also printed).
//
// `--stats-only` skips the vector files and only writes stats.txt. Its `model_*` fields are the
// per-shot cycle model of the RTL evaluated on the sampled shots for every lane count in
// {16, 32, 64, 128}, both for the compressed tables and for the legacy layout (one full ball per
// detector, absolute targets), so a single run reports the change in table size and in latency.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "specmatching/spec_matching/manifold/ball_graph.h"
#include "specmatching/spec_matching/manifold/ball_tables.h"
#include "tests/spec_matching/spec_matching_test_util.h"

namespace {

using pm::spec_matching::BallGraphArena;
using pm::spec_matching::BallGraphBuildMode;
using pm::spec_matching::BallParams;
using pm::spec_matching::BallTables;

/// Integer weight width. Every stored weight is <= 2T (entries) or <= T (boundary edges), so it is
/// sized from 2T in run(), with the 16 bits of the original design as the floor.
size_t WEIGHT_W = 16;
const size_t MODEL_LANES[] = {16, 32, 64, 128};

struct Options {
    std::string out_dir;
    size_t d = 5;
    size_t rounds = 5;
    double p = 0.005;
    size_t shots = 20;
    size_t keep = 0;  // 0: keep every sampled shot
    uint64_t seed = 7;
    size_t lanes = 64;
    long time_coord = -1;  // index of the time coordinate; negative counts from the end
    bool prune = true;
    bool dedup = true;
    bool stats_only = false;
};

size_t clog2(uint64_t x) {
    size_t b = 0;
    while (((uint64_t)1 << b) < x)
        b++;
    return b;
}

size_t ceil_div(size_t a, size_t b) {
    return (a + b - 1) / b;
}

/// A fixed-width little-endian bit vector, written out as one `$readmemh` line.
struct Word {
    std::vector<uint64_t> limbs;
    size_t width;

    explicit Word(size_t width) : limbs(ceil_div(width, 64) + 1, 0), width(width) {
    }

    void set(size_t pos, size_t len, uint64_t value) {
        for (size_t b = 0; b < len; b++) {
            if ((value >> b) & 1)
                limbs[(pos + b) / 64] |= (uint64_t)1 << ((pos + b) % 64);
        }
    }

    std::string hex() const {
        static const char* digits = "0123456789abcdef";
        size_t n = ceil_div(width, 4);
        std::string s(n, '0');
        for (size_t k = 0; k < n; k++) {
            size_t pos = 4 * k;
            unsigned nib = (unsigned)((limbs[pos / 64] >> (pos % 64)) & 0xf);
            s[n - 1 - k] = digits[nib];
        }
        return s;
    }
};

using Entry = std::pair<int64_t, uint32_t>;  // {delta, w}

/// The accelerator's view of one ball table compile.
struct Compressed {
    std::vector<std::vector<Entry>> rows;  // unique rows
    std::vector<uint32_t> det_row;         // row of each detector (0 when its row is empty)
    std::vector<uint32_t> det_size;        // entries in that row
    std::vector<uint8_t> det_hb;
    std::vector<uint32_t> det_bw;
    size_t max_row{0};
    size_t entries_stored{0};
    size_t entries_logical{0};  // sum of det_size, i.e. before dedup
    int64_t min_delta{-1};      // -1: the padding entry has to be representable too
    int64_t max_delta{0};
    size_t num_space_positions{0};
};

FILE* open_or_die(const std::string& path, const char* mode) {
    FILE* f = std::fopen(path.c_str(), mode);
    if (f == nullptr)
        throw std::runtime_error("cannot open " + path);
    return f;
}

Compressed compress(
    const BallTables& tables,
    const stim::DetectorErrorModel& dem,
    pm::cumulative_time_int two_t,
    pm::cumulative_time_int horizon,
    const Options& opt) {
    size_t n = tables.num_nodes;
    Compressed c;
    c.det_row.assign(n, 0);
    c.det_size.assign(n, 0);
    c.det_hb.assign(n, 0);
    c.det_bw.assign(n, 0);

    std::set<uint64_t> all_dets;
    for (uint64_t k = 0; k < n; k++)
        all_dets.insert(k);
    std::map<uint64_t, std::vector<double>> coords = dem.get_detector_coordinates(all_dets);

    // A detector without coordinates gets its own space key, so it never shares a row.
    using SpaceKey = std::pair<int64_t, std::vector<double>>;
    std::map<std::pair<SpaceKey, std::vector<Entry>>, uint32_t> index;
    std::set<SpaceKey> spaces;

    for (uint64_t det = 0; det < n; det++) {
        std::vector<Entry> row;
        for (uint64_t e = tables.ball_begin(det); e < tables.ball_end(det); e++) {
            pm::weight_int w = tables.ball_w_int[e];
            if ((pm::cumulative_time_int)w > two_t)
                break;  // entries are sorted by (w, target)
            uint32_t target = tables.ball_target[e];
            if (opt.prune && target <= det)
                continue;
            if ((uint64_t)w >= ((uint64_t)1 << WEIGHT_W))
                throw std::runtime_error("ball weight does not fit WEIGHT_W");
            row.emplace_back((int64_t)target - (int64_t)det, (uint32_t)w);
        }

        if (tables.has_bcost[det] && (pm::cumulative_time_int)tables.bcost_w_int[det] <= horizon) {
            if ((uint64_t)tables.bcost_w_int[det] >= ((uint64_t)1 << WEIGHT_W))
                throw std::runtime_error("boundary weight does not fit WEIGHT_W");
            c.det_hb[det] = 1;
            c.det_bw[det] = (uint32_t)tables.bcost_w_int[det];
        }

        SpaceKey space;
        const std::vector<double>& xyz = coords[det];
        if (xyz.empty()) {
            space = {(int64_t)det, {}};
        } else {
            long t = opt.time_coord < 0 ? (long)xyz.size() + opt.time_coord : opt.time_coord;
            if (t < 0 || t >= (long)xyz.size())
                throw std::runtime_error("--time-coord is outside detector " + std::to_string(det) + "'s coordinates");
            space.first = -1;
            for (long k = 0; k < (long)xyz.size(); k++) {
                if (k != t)
                    space.second.push_back(xyz[k]);
            }
        }
        spaces.insert(space);

        c.det_size[det] = (uint32_t)row.size();
        c.entries_logical += row.size();
        c.max_row = std::max(c.max_row, row.size());
        for (const Entry& en : row) {
            c.min_delta = std::min(c.min_delta, en.first);
            c.max_delta = std::max(c.max_delta, en.first);
        }
        if (row.empty())
            continue;  // no ball read is issued for an empty row, so it needs no word

        auto key = std::make_pair(space, row);
        auto it = opt.dedup ? index.find(key) : index.end();
        if (it != index.end()) {
            c.det_row[det] = it->second;
        } else {
            uint32_t id = (uint32_t)c.rows.size();
            c.rows.push_back(row);
            c.entries_stored += row.size();
            c.det_row[det] = id;
            if (opt.dedup)
                index.emplace(std::move(key), id);
        }
    }
    c.num_space_positions = spaces.size();
    return c;
}

/// Smallest two's complement width holding every value in [lo, hi].
size_t signed_width(int64_t lo, int64_t hi) {
    size_t w = 1;
    while (lo < -((int64_t)1 << (w - 1)) || hi > ((int64_t)1 << (w - 1)) - 1)
        w++;
    return w;
}

/// The RTL's per-shot cycle count, from the iStart edge to the edge oDone is seen on, inclusive --
/// what sim_main.cc reports.
///
/// A register-level replay of the control path, not a closed form: the address generator's
/// slots (one per unstalled S_SCAN cycle, `max(1, chunks)` per defect, a chunk issued in each slot
/// of a non-empty ball), the one-cycle SRAM/sideband stage, the batch FIFO and its `cnt >= 3`
/// stall, the lane pipeline of depth `1 + clog2(cap + 1)` frozen while the collector is not ready,
/// the collector draining one hit per cycle, and the done FSM's quiet test. Only the hit count of
/// each batch is carried, since that is all the timing depends on. It is the RTL's sum
/// `7 + clog2(cap + 1) + fetch + drain` except where scan slots hide under a drain (e.g. a trailing
/// empty ball, which pruning makes common), and run_crosscheck.sh checks it cycle-exact.
///
/// `cap` is the syndrome capacity NUM_SYND of the instance; `cap = 0` uses the per-shot `M`, which
/// is how the latency sweep reports it (a fixed-capacity instance pays the difference). `fetch`
/// (scan slots) and `drain` (hits beyond the first per batch) are reported alongside.
struct ShotModel {
    uint64_t cycles{0};
    uint64_t fetch{0};
    uint64_t drain{0};
};

template <typename RowOf>
ShotModel model_shot(
    const std::vector<uint64_t>& shot,
    const std::vector<uint32_t>& is_defect,
    const std::vector<uint8_t>& hb,
    size_t lanes,
    size_t cap,
    RowOf row_of) {
    ShotModel m;
    size_t num = shot.size();
    if (num == 0)
        return m;

    // Per defect, the hit count of each of its chunks.
    std::vector<std::vector<uint32_t>> hits(num);
    for (size_t i = 0; i < num; i++) {
        uint64_t det = shot[i];
        const std::vector<uint32_t>& targets = row_of(det);  // absolute targets, in row order
        size_t chunks = ceil_div(targets.size(), lanes);
        m.fetch += std::max<size_t>(1, chunks);
        for (size_t ch = 0; ch < chunks; ch++) {
            uint32_t h = 0;
            for (size_t k = ch * lanes; k < std::min(targets.size(), (ch + 1) * lanes); k++) {
                uint32_t t = targets[k];
                if (t > det && is_defect[t])
                    h++;
            }
            hits[i].push_back(h);
            if (h > 1)
                m.drain += h - 1;
        }
    }

    enum { S_IDLE, S_NODE, S_HDR, S_SCAN };
    const size_t lat = 1 + clog2((cap > 0 ? cap : num) + 1);
    const int fifo_high = 3;  // hc_fifo: FIFO_DEPTH - 1

    int ag = S_IDLE;  // hc_addr_gen
    size_t i = 0, ni = 0, ccnt = 0;
    bool rd_v = false;  // hc_fifo's SRAM-aligned sideband stage
    uint32_t rd_hits = 0;
    std::vector<uint32_t> fifo;  // hit count of each buffered batch
    size_t fifo_head = 0;
    std::vector<int64_t> pipe(lat, -1);  // lane pipeline, -1 = bubble
    uint32_t pend = 0;                   // hc_collect: hits still to emit
    bool bnd_v = false;                  // hc_bnd_check's oBndValid
    bool running = false;                // hc_fifo S_RUN

    for (uint64_t t = 0;; t++) {
        size_t cnt = fifo.size() - fifo_head;
        bool start = (t == 0);
        bool go = !running && start;
        bool stall = (int)cnt >= fifo_high;
        bool empty_c = (ccnt == 0);
        bool last_c = empty_c || ccnt == 1;
        bool step = ag == S_SCAN && !stall;
        bool load = ag == S_HDR || (step && last_c && i != num - 1);
        bool bnd_rd = ag == S_NODE || (load && ni != num);
        bool ch_valid = step && !empty_c;
        uint32_t ch_hits = ch_valid ? hits[i][hits[i].size() - ccnt] : 0;
        bool ready = pend <= 1;
        bool take = ready && pipe[lat - 1] >= 0;
        bool bnd_busy = load || bnd_v;
        bool quiet = ag == S_IDLE && !rd_v && cnt == 0 && pend == 0 && !bnd_busy &&
                     std::all_of(pipe.begin(), pipe.end(), [](int64_t v) { return v < 0; });

        // ---- registers
        size_t load_det = (ag == S_HDR) ? 0 : i + 1;
        switch (ag) {
            case S_IDLE:
                if (go) {
                    i = 0;
                    ni = 0;
                    ag = S_NODE;
                }
                break;
            case S_NODE:
                ag = S_HDR;
                break;
            case S_HDR:
                ag = S_SCAN;
                break;
            default:
                if (!stall) {
                    ccnt--;
                    if (last_c) {
                        if (i == num - 1)
                            ag = S_IDLE;
                        else
                            i++;
                    }
                }
        }
        if (bnd_rd)
            ni++;
        if (load)
            ccnt = hits[load_det].size();
        bnd_v = load && hb[shot[load_det]];

        int64_t head = cnt > 0 ? (int64_t)fifo[fifo_head] : -1;
        if (ready && cnt > 0)
            fifo_head++;
        if (rd_v)
            fifo.push_back(rd_hits);
        rd_v = ch_valid;
        rd_hits = ch_hits;

        uint32_t out_hits = pipe[lat - 1] >= 0 ? (uint32_t)pipe[lat - 1] : 0;
        if (ready) {
            for (size_t k = lat - 1; k > 0; k--)
                pipe[k] = pipe[k - 1];
            pipe[0] = head;
        }
        pend = take ? out_hits : (pend > 0 ? pend - 1 : 0);

        bool done_now = running && quiet;
        if (go)
            running = true;
        if (done_now) {
            m.cycles = t + 1;
            return m;
        }
    }
}

void parse_args(int argc, char** argv, Options& opt) {
    if (argc < 2)
        throw std::invalid_argument("usage: hc_dump <out-dir> [options]");
    opt.out_dir = argv[1];
    for (int k = 2; k < argc; k++) {
        std::string a = argv[k];
        auto next = [&]() -> std::string {
            if (k + 1 >= argc)
                throw std::invalid_argument("missing value for " + a);
            return argv[++k];
        };
        if (a == "--d")
            opt.d = std::stoul(next());
        else if (a == "--rounds")
            opt.rounds = std::stoul(next());
        else if (a == "--p")
            opt.p = std::stod(next());
        else if (a == "--shots")
            opt.shots = std::stoul(next());
        else if (a == "--keep")
            opt.keep = std::stoul(next());
        else if (a == "--seed")
            opt.seed = std::stoull(next());
        else if (a == "--lanes")
            opt.lanes = std::stoul(next());
        else if (a == "--time-coord")
            opt.time_coord = std::stol(next());
        else if (a == "--no-prune")
            opt.prune = false;
        else if (a == "--no-dedup")
            opt.dedup = false;
        else if (a == "--stats-only")
            opt.stats_only = true;
        else
            throw std::invalid_argument("unknown option " + a);
    }
    if (opt.lanes < 2 || (opt.lanes & (opt.lanes - 1)) != 0)
        throw std::invalid_argument("--lanes must be a power of two >= 2");
}

int run(int argc, char** argv) {
    Options opt;
    parse_args(argc, argv, opt);

    pm::spec_matching::test::Corpus corpus =
        pm::spec_matching::test::generate_surface_code_corpus(opt.d, opt.rounds, opt.p, opt.shots, opt.seed);
    if (opt.keep > 0 && opt.keep < corpus.shots.size())
        corpus.shots.resize(opt.keep);
    pm::Mwpm mwpm = corpus.to_mwpm();
    const pm::MatchingGraph& graph = mwpm.flooder.graph;
    if (!mwpm.flooder.negative_weight_detection_events.empty() || mwpm.flooder.negative_weight_sum != 0)
        throw std::runtime_error("negative weights: the shots would need the preamble applied");

    double unit = pm::spec_matching::test::edge_weight_units(graph);
    BallParams params;
    params.T_max = 2.0 * unit;
    params.R = 4.0 * unit;
    BallTables tables = pm::spec_matching::compile_ball_tables(mwpm, params);
    pm::cumulative_time_int horizon = tables.t_max_int;
    pm::cumulative_time_int two_t = 2 * horizon;
    size_t n = tables.num_nodes;
    WEIGHT_W = std::max<size_t>(16, clog2((uint64_t)two_t + 1));

    Compressed c = compress(tables, corpus.dem, two_t, horizon, opt);

    // ---- geometry
    size_t log_num_det = std::max<size_t>(1, clog2(n));
    size_t delta_w = signed_width(c.min_delta, c.max_delta);
    if (delta_w > log_num_det)
        log_num_det = delta_w;  // the lane sign-extends delta to the detector id width
    size_t num_rows = std::max<size_t>(1, c.rows.size());
    size_t log_num_rows = std::max<size_t>(1, clog2(num_rows));
    size_t log_max_ball = std::max<size_t>(1, clog2(c.max_row + 1));
    size_t lane_w = WEIGHT_W + delta_w;
    size_t chunks_per_ball = std::max<size_t>(1, ceil_div(c.max_row, opt.lanes));
    size_t chunk_w = opt.lanes * lane_w;
    size_t ball_width = chunks_per_ball * chunk_w;
    size_t bnd_width = log_num_rows + 1 + WEIGHT_W + log_max_ball;

    // ---- legacy layout, for the comparison: one full ball per detector, absolute targets
    size_t legacy_max = 0;
    for (uint64_t det = 0; det < n; det++) {
        size_t k = 0;
        for (uint64_t e = tables.ball_begin(det); e < tables.ball_end(det); e++)
            if ((pm::cumulative_time_int)tables.ball_w_int[e] <= two_t)
                k++;
        legacy_max = std::max(legacy_max, k);
    }
    size_t legacy_log_num_det = std::max<size_t>(1, clog2(n));
    size_t legacy_lane_w = WEIGHT_W + legacy_log_num_det;
    size_t legacy_bnd_width = 1 + WEIGHT_W + std::max<size_t>(1, clog2(legacy_max + 1));

    // ---- shots, and the cycle model over them
    std::vector<std::vector<uint32_t>> new_targets(n), legacy_targets(n);
    for (uint64_t det = 0; det < n; det++) {
        if (c.det_size[det] > 0)
            for (const Entry& en : c.rows[c.det_row[det]])
                new_targets[det].push_back((uint32_t)((int64_t)det + en.first));
        for (uint64_t e = tables.ball_begin(det); e < tables.ball_end(det); e++)
            if ((pm::cumulative_time_int)tables.ball_w_int[e] <= two_t)
                legacy_targets[det].push_back(tables.ball_target[e]);
    }
    std::vector<uint32_t> is_defect(n, 0);
    size_t max_m = 0;
    double sum_m = 0;
    for (const auto& shot : corpus.shots)
        max_m = std::max(max_m, shot.size());
    size_t num_synd = std::max<size_t>(1, max_m);
    std::vector<double> sum_new(4, 0), sum_legacy(4, 0), sum_new_fetch(4, 0), sum_legacy_fetch(4, 0),
        sum_new_drain(4, 0), sum_legacy_drain(4, 0);
    std::vector<uint64_t> max_new(4, 0), max_legacy(4, 0);
    std::vector<uint64_t> shot_cycles;  // at --lanes, compressed tables, capacity = num_synd
    for (const auto& shot : corpus.shots) {
        sum_m += (double)shot.size();
        for (uint64_t det : shot)
            is_defect[det] = 1;
        for (size_t li = 0; li < 4; li++) {
            size_t L = MODEL_LANES[li];
            ShotModel a = model_shot(shot, is_defect, c.det_hb, L, 0, [&](uint64_t det) -> const std::vector<uint32_t>& {
                return new_targets[det];
            });
            ShotModel b = model_shot(shot, is_defect, c.det_hb, L, 0, [&](uint64_t det) -> const std::vector<uint32_t>& {
                return legacy_targets[det];
            });
            sum_new[li] += (double)a.cycles;
            sum_legacy[li] += (double)b.cycles;
            sum_new_fetch[li] += (double)a.fetch;
            sum_legacy_fetch[li] += (double)b.fetch;
            sum_new_drain[li] += (double)a.drain;
            sum_legacy_drain[li] += (double)b.drain;
            max_new[li] = std::max(max_new[li], a.cycles);
            max_legacy[li] = std::max(max_legacy[li], b.cycles);
        }
        shot_cycles.push_back(
            model_shot(shot, is_defect, c.det_hb, opt.lanes, num_synd, [&](uint64_t det) -> const std::vector<uint32_t>& {
                return new_targets[det];
            }).cycles);
        for (uint64_t det : shot)
            is_defect[det] = 0;
    }

    // ---- stats
    std::ostringstream st;
    double shots = (double)std::max<size_t>(1, corpus.shots.size());
    st << "d=" << opt.d << " rounds=" << opt.rounds << " p=" << opt.p << " shots=" << corpus.shots.size()
       << " seed=" << opt.seed << " prune=" << opt.prune << " dedup=" << opt.dedup << " lanes=" << opt.lanes << " weight_w=" << WEIGHT_W
       << " num_dets=" << n << " space_positions=" << c.num_space_positions << " rows=" << c.rows.size()
       << " entries_logical=" << c.entries_logical << " entries_stored=" << c.entries_stored
       << " max_row=" << c.max_row << " delta_min=" << c.min_delta << " delta_max=" << c.max_delta
       << " delta_w=" << delta_w << " lane_w=" << lane_w << " chunks_per_ball=" << chunks_per_ball
       << " ball_width=" << ball_width << " bnd_width=" << bnd_width
       << " ball_bits=" << (uint64_t)c.rows.size() * ball_width << " bnd_bits=" << (uint64_t)n * bnd_width
       << " legacy_max_ball=" << legacy_max << " legacy_lane_w=" << legacy_lane_w
       << " legacy_chunks_per_ball=" << std::max<size_t>(1, ceil_div(legacy_max, opt.lanes))
       << " legacy_ball_bits="
       << (uint64_t)n * std::max<size_t>(1, ceil_div(legacy_max, opt.lanes)) * opt.lanes * legacy_lane_w
       << " legacy_bnd_bits=" << (uint64_t)n * legacy_bnd_width << " mean_m=" << sum_m / shots
       << " max_m=" << max_m;
    for (size_t li = 0; li < 4; li++) {
        size_t L = MODEL_LANES[li];
        st << " model_L" << L << "_mean=" << sum_new[li] / shots << " model_L" << L << "_fetch=" << sum_new_fetch[li] / shots
           << " model_L" << L << "_drain=" << sum_new_drain[li] / shots << " model_L" << L << "_max=" << max_new[li]
           << " legacy_L" << L << "_mean=" << sum_legacy[li] / shots << " legacy_L" << L
           << "_fetch=" << sum_legacy_fetch[li] / shots << " legacy_L" << L << "_drain=" << sum_legacy_drain[li] / shots
           << " legacy_L" << L << "_max=" << max_legacy[li];
    }
    std::string stats = st.str();
    std::printf("%s\n", stats.c_str());
    FILE* sf = open_or_die(opt.out_dir + "/stats.txt", "w");
    std::fprintf(sf, "%s\n", stats.c_str());
    std::fclose(sf);
    if (opt.stats_only)
        return 0;

    // ---- ball.hex: unique rows, padded to the SRAM depth with all-padding words
    uint64_t pad_delta = ((uint64_t)1 << delta_w) - 1;  // -1
    size_t ball_depth = (size_t)1 << log_num_rows;
    FILE* bf = open_or_die(opt.out_dir + "/ball.hex", "w");
    for (size_t r = 0; r < ball_depth; r++) {
        Word w(ball_width);
        for (size_t k = 0; k < chunks_per_ball * opt.lanes; k++) {
            size_t pos = k * lane_w;
            if (r < c.rows.size() && k < c.rows[r].size()) {
                const Entry& en = c.rows[r][k];
                w.set(pos, delta_w, (uint64_t)en.first & pad_delta);
                w.set(pos + delta_w, WEIGHT_W, en.second);
            } else {
                w.set(pos, delta_w, pad_delta);
            }
        }
        std::fprintf(bf, "%s\n", w.hex().c_str());
    }
    std::fclose(bf);

    // ---- bnd.hex: {row, hb, bw, size} per detector, padded to the SRAM depth with zeros
    size_t node_depth = (size_t)1 << log_num_det;
    FILE* nf = open_or_die(opt.out_dir + "/bnd.hex", "w");
    for (size_t det = 0; det < node_depth; det++) {
        Word w(bnd_width);
        if (det < n) {
            w.set(0, log_max_ball, c.det_size[det]);
            w.set(log_max_ball, WEIGHT_W, c.det_bw[det]);
            w.set(log_max_ball + WEIGHT_W, 1, c.det_hb[det]);
            w.set(log_max_ball + WEIGHT_W + 1, log_num_rows, c.det_row[det]);
        }
        std::fprintf(nf, "%s\n", w.hex().c_str());
    }
    std::fclose(nf);

    // ---- params.vh
    FILE* pf = open_or_die(opt.out_dir + "/params.vh", "w");
    std::fprintf(pf, "// Generated by accelerator/tablegen/hc_dump.cc -- do not edit.\n");
    std::fprintf(
        pf,
        "// Table: generate_surface_code_corpus(%zu, %zu, %g, %zu, %llu), T = T_max = 2 units, R = 4 units, "
        "prune=%d, dedup=%d.\n",
        opt.d, opt.rounds, opt.p, opt.shots, (unsigned long long)opt.seed, (int)opt.prune, (int)opt.dedup);
    std::fprintf(pf, "`ifndef HC_PARAMS_VH\n`define HC_PARAMS_VH\n");
    auto def = [&](const char* k, uint64_t v) { std::fprintf(pf, "`define %s %llu\n", k, (unsigned long long)v); };
    def("HC_NUM_SYND", num_synd);
    def("HC_LOG_NUM_DET", log_num_det);
    def("HC_LOG_NUM_ROWS", log_num_rows);
    def("HC_LOG_MAX_BALL", log_max_ball);
    def("HC_WEIGHT_W", WEIGHT_W);
    def("HC_DELTA_W", delta_w);
    def("HC_NUM_LANES", opt.lanes);
    def("HC_LANE_W", lane_w);
    def("HC_CHUNK_W", chunk_w);
    def("HC_CHUNKS_PER_BALL", chunks_per_ball);
    def("HC_NUM_DETS", n);
    def("HC_NUM_ROWS", c.rows.size());
    def("HC_BND_DEPTH", node_depth);
    def("HC_BALL_DEPTH", ball_depth);
    def("HC_BND_WIDTH", bnd_width);
    def("HC_BALL_WIDTH", ball_width);
    def("HC_NUM_SHOTS", corpus.shots.size());
    def("HC_MAX_BALL", c.max_row);
    std::fprintf(pf, "`define HC_BALL_HEX \"ball.hex\"\n`define HC_BND_HEX \"bnd.hex\"\n`endif\n");
    std::fclose(pf);

    // ---- per shot: syndrome, reference H, model cycles
    BallGraphArena arena;
    FILE* mf = open_or_die(opt.out_dir + "/model_cycles.txt", "w");
    for (size_t k = 0; k < corpus.shots.size(); k++) {
        const std::vector<uint64_t>& shot = corpus.shots[k];
        FILE* f = open_or_die(opt.out_dir + "/syn_" + std::to_string(k) + ".txt", "w");
        std::fprintf(f, "%zu\n", shot.size());
        for (uint64_t det : shot)
            std::fprintf(f, "%llu\n", (unsigned long long)det);
        std::fclose(f);

        pm::spec_matching::build_ball_graph(tables, shot, (pm::spec_matching::horizon_int)horizon, arena, BallGraphBuildMode::SCAN);
        f = open_or_die(opt.out_dir + "/ref_" + std::to_string(k) + ".txt", "w");
        for (const auto& b : arena.graph.boundary_edges)
            std::fprintf(f, "B %u %llu\n", b.i, (unsigned long long)b.w_int);
        for (const auto& e : arena.graph.edges)
            std::fprintf(f, "E %u %u %llu\n", e.i, e.j, (unsigned long long)e.w_int);
        std::fclose(f);
        std::fprintf(mf, "%zu %llu\n", k, (unsigned long long)shot_cycles[k]);
    }
    std::fclose(mf);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "hc_dump: %s\n", ex.what());
        return 2;
    }
}
