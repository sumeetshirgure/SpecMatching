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

/// §2.1 of `docs/lattice_surgery_profiler.md` — the `1×N` merged-strip circuit.
///
/// `stim::generate_surface_code_circuit` makes square patches only, so this builds the rotated
/// surface code on a `d × L` rectangle itself, with the noise channels of `profiler_util.h`'s
/// `generate` (after-Clifford depolarisation, after-reset flip and before-measure flip, all at `p`).
///
/// ## Geometry
///
/// Data qubit `(i, j)`, column `i ∈ [0, L)` along the strip, row `j ∈ [0, d)`. Plaquette `(i, j)`
/// for `i ∈ [-1, L-1]`, `j ∈ [-1, d-1]` covers the in-range data of `(i, j), (i+1, j), (i, j+1),
/// (i+1, j+1)`; its type is the checkerboard parity of `i + j`. Weight-2 plaquettes on the top and
/// bottom edges are kept only for one type (`TB`) and on the left and right ends only for the other
/// (`LR`); corners are dropped. A data qubit on the top row is then in two `TB` plaquettes and one
/// `LR` plaquette, so **the `LR` type's matching graph has its boundary along the long sides**.
///
///  - `--orientation long`: the decoded (X) type is `LR`. Its boundaries run along the strip, so
///    every detector is within about `d/2` of one.
///  - `--orientation short`: the decoded type is `TB`. Its boundaries are the two short ends, and a
///    detector in the middle of the strip can be `~Nd/2` away from both.
///
/// The memory basis is X throughout (data reset with `RX`, measured with `MX`), as in the corpus's
/// `rotated_memory_x`. Observable 0 is the logical X string that commutes with every Z stabilizer:
/// row `j = 0` for `long`, column `i = 0` for `short`.
///
/// ## What differs from the memory corpus, deliberately
///
///  - **Only the decoded type's detectors are declared.** The doc's orientation argument is about
///    "the boundaries the decoded detector type sees", and declaring the other type too would put
///    a second, oppositely oriented matching graph in every cell. A consequence: a data `Y` error is
///    graphlike here (it flips X detectors only), where in the corpus it flips four detectors of two
///    types and `iter_detector_error_model_edges` drops it. The DEM call is the corpus's own
///    (`decompose_errors = false`).
///  - **Detectors are declared cube-major** (cube, then round, then plaquette), after the last
///    measurement, so that the detector ids of one cube are contiguous. The node order only breaks
///    ties; it is what lets every per-cube structure in the profiler be a contiguous range.
///
/// ## Merge initialisation (`--merge-init on`)
///
/// The **seam** plaquettes of the decoded type — those whose support spans two patches, or touches
/// a routing column — get no detector in their first round, because in a merge their first outcome
/// is random. The product of those first-round outcomes is observable 1, the joint logical
/// measurement. With `--merge-init off` every plaquette gets its first-round detector and there is
/// no observable 1: a plain rectangular memory experiment, kept as a control.
///
/// ## CNOT schedule
///
/// Four layers. Both types start at corner `(i, j)` and end at `(i+1, j+1)`; they differ in the
/// order of the middle two, which is what makes the measurements commute. The Z ancilla's last two
/// targets are its hook pair; they are placed **parallel** to the decoded type's logical error
/// string (horizontal for `long`, vertical for `short`), so a hook does not shorten it.

#ifndef SPECMATCHING_BENCHMARKS_SPEC_MATCHING_LATTICE_SURGERY_CIRCUIT_H
#define SPECMATCHING_BENCHMARKS_SPEC_MATCHING_LATTICE_SURGERY_CIRCUIT_H

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "stim.h"

namespace pm {
namespace spec_matching {
namespace profiler {

struct LatticeSurgeryParams {
    size_t distance{5};
    size_t patches{2};
    /// 0: patches abut, `L = N d`. 1: one routing column between neighbours, `L = N d + N - 1`.
    size_t seam_gap{0};
    size_t rounds{5};
    double noise{1e-3};
    bool merge_init{true};
    /// `true` for `--orientation long`.
    bool long_orientation{true};
};

struct LatticeSurgeryCircuit {
    LatticeSurgeryParams params;
    stim::Circuit circuit;
    size_t width{0};
    size_t height{0};
    size_t num_detectors{0};
    size_t num_observables{0};
    /// Cube of every detector; ids of one cube are contiguous, `[cube_offsets[c], cube_offsets[c+1])`.
    std::vector<uint32_t> detector_cube;
    std::vector<uint32_t> cube_offsets;
    /// Seam plaquettes of the decoded type, and how many first-round detectors merge-init removed.
    size_t seam_plaquettes{0};
    size_t removed_first_round_detectors{0};

    /// The cube a data column belongs to. A routing column belongs to the cube on its left.
    static size_t cube_of_column(const LatticeSurgeryParams& p, long column) {
        if (column < 0)
            column = 0;
        size_t period = p.distance + p.seam_gap;
        size_t cube = (size_t)column / period;
        return std::min(cube, p.patches - 1);
    }
    static bool is_routing_column(const LatticeSurgeryParams& p, long column) {
        if (p.seam_gap == 0 || column < 0)
            return false;
        size_t period = p.distance + p.seam_gap;
        return (size_t)column % period >= p.distance;
    }

    /// The generator call as one line, for `run.log`.
    std::string describe() const {
        std::ostringstream out;
        out << "lattice_surgery_circuit(distance=" << params.distance << ", patches=" << params.patches
            << ", seam_gap=" << params.seam_gap << ", rounds=" << params.rounds << ", width=" << width
            << ", height=" << height << ", orientation=" << (params.long_orientation ? "long" : "short")
            << ", merge_init=" << (params.merge_init ? "on" : "off") << ", basis=X, detectors=decoded type only"
            << ", cnot_order=hook-parallel-to-logical-error"
            << "); after_clifford_depolarization=" << params.noise
            << "; after_reset_flip_probability=" << params.noise
            << "; before_measure_flip_probability=" << params.noise
            << "; before_round_data_depolarization=0 (unset)"
            << "; stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false)";
        return out.str();
    }
};

inline LatticeSurgeryCircuit generate_lattice_surgery_circuit(const LatticeSurgeryParams& p) {
    if (p.distance < 3 || p.distance % 2 == 0)
        throw std::invalid_argument("lattice surgery circuit: distance must be odd and at least 3");
    if (p.patches == 0)
        throw std::invalid_argument("lattice surgery circuit: need at least one patch");
    if (p.seam_gap > 1)
        throw std::invalid_argument("lattice surgery circuit: seam gap must be 0 or 1");
    if (p.rounds == 0)
        throw std::invalid_argument("lattice surgery circuit: need at least one round");

    LatticeSurgeryCircuit out;
    out.params = p;
    const long W = (long)(p.patches * p.distance + (p.patches - 1) * p.seam_gap);
    const long H = (long)p.distance;
    out.width = (size_t)W;
    out.height = (size_t)H;

    // Types: TB is the type kept on the top and bottom edges, LR the one kept on the ends.
    const int TB = 0;
    const int LR = 1;
    const int x_type = p.long_orientation ? LR : TB;

    auto data_index = [&](long i, long j) -> uint32_t {
        return (uint32_t)(j * W + i);
    };
    auto parity = [](long v) -> int {
        return (int)(((v % 2) + 2) % 2);
    };

    struct Plaquette {
        long i;
        long j;
        int type;
        bool is_x;
        bool seam;
        uint32_t cube;
        uint32_t ancilla;
        /// The four corners in CNOT-layer order; `UINT32_MAX` where the corner is off the patch.
        uint32_t order[4];
        std::vector<uint32_t> support;
    };
    std::vector<Plaquette> plaquettes;
    uint32_t next_qubit = (uint32_t)(W * H);

    for (long j = -1; j <= H - 1; j++) {
        for (long i = -1; i <= W - 1; i++) {
            bool top_bottom = (j == -1 || j == H - 1);
            bool left_right = (i == -1 || i == W - 1);
            if (top_bottom && left_right)
                continue;
            int type = parity(i + j);
            if (top_bottom && type != TB)
                continue;
            if (left_right && type != LR)
                continue;
            Plaquette q;
            q.i = i;
            q.j = j;
            q.type = type;
            q.is_x = type == x_type;
            // Corners c00 = (i, j), c10 = (i+1, j), c01 = (i, j+1), c11 = (i+1, j+1).
            uint32_t corner[4];
            long ci[4] = {i, i + 1, i, i + 1};
            long cj[4] = {j, j, j + 1, j + 1};
            for (int k = 0; k < 4; k++) {
                bool inside = ci[k] >= 0 && ci[k] < W && cj[k] >= 0 && cj[k] < H;
                corner[k] = inside ? data_index(ci[k], cj[k]) : UINT32_MAX;
                if (inside)
                    q.support.push_back(corner[k]);
            }
            // Z hook (its last two targets) parallel to the decoded type's logical error string:
            // horizontal (c01, c11) for `long`, vertical (c10, c11) for `short`. X uses the other
            // middle order.
            bool z_horizontal_hook = p.long_orientation;
            bool use_horizontal_last_pair = q.is_x ? !z_horizontal_hook : z_horizontal_hook;
            if (use_horizontal_last_pair) {
                q.order[0] = corner[0];
                q.order[1] = corner[1];
                q.order[2] = corner[2];
                q.order[3] = corner[3];
            } else {
                q.order[0] = corner[0];
                q.order[1] = corner[2];
                q.order[2] = corner[1];
                q.order[3] = corner[3];
            }
            // Columns spanned: i and i+1, clipped.
            long lo = std::max<long>(i, 0);
            long hi = std::min<long>(i + 1, W - 1);
            bool spans_two = LatticeSurgeryCircuit::cube_of_column(p, lo) != LatticeSurgeryCircuit::cube_of_column(p, hi) ||
                             (p.seam_gap == 1 && (LatticeSurgeryCircuit::is_routing_column(p, lo) ||
                                                  LatticeSurgeryCircuit::is_routing_column(p, hi)));
            q.seam = p.patches > 1 && spans_two;
            q.cube = (uint32_t)LatticeSurgeryCircuit::cube_of_column(p, i);
            q.ancilla = next_qubit++;
            plaquettes.push_back(std::move(q));
        }
    }
    if ((long)plaquettes.size() != W * H - 1)
        throw std::logic_error(
            "lattice surgery circuit: expected W*H-1 = " + std::to_string(W * H - 1) + " stabilizers, built " +
            std::to_string(plaquettes.size()));

    std::vector<uint32_t> data;
    for (long j = 0; j < H; j++)
        for (long i = 0; i < W; i++)
            data.push_back(data_index(i, j));
    std::vector<uint32_t> ancillas, x_ancillas;
    for (const Plaquette& q : plaquettes) {
        ancillas.push_back(q.ancilla);
        if (q.is_x)
            x_ancillas.push_back(q.ancilla);
    }

    stim::Circuit& c = out.circuit;
    auto op = [&](const char* name, const std::vector<uint32_t>& targets, double arg = -1) {
        if (targets.empty())
            return;
        std::ostringstream text;
        text << name;
        if (arg >= 0)
            text << "(" << arg << ")";
        for (uint32_t t : targets)
            text << " " << t;
        text << "\n";
        c.append_from_text(text.str());
    };
    for (long j = 0; j < H; j++)
        for (long i = 0; i < W; i++) {
            std::ostringstream text;
            text << "QUBIT_COORDS(" << 2 * i + 1 << ", " << 2 * j + 1 << ") " << data_index(i, j) << "\n";
            c.append_from_text(text.str());
        }
    for (const Plaquette& q : plaquettes) {
        std::ostringstream text;
        text << "QUBIT_COORDS(" << 2 * q.i + 2 << ", " << 2 * q.j + 2 << ") " << q.ancilla << "\n";
        c.append_from_text(text.str());
    }

    // Measurement record bookkeeping: absolute index of every measurement.
    size_t measurement_count = 0;
    // per round, per plaquette index -> absolute measurement index
    std::vector<std::vector<size_t>> stab_meas(p.rounds, std::vector<size_t>(plaquettes.size(), 0));
    std::vector<size_t> data_meas(data.size(), 0);

    const double noise = p.noise;
    op("RX", data);
    op("Z_ERROR", data, noise);
    for (size_t r = 0; r < p.rounds; r++) {
        op("R", ancillas);
        op("X_ERROR", ancillas, noise);
        op("H", x_ancillas);
        op("DEPOLARIZE1", x_ancillas, noise);
        for (int layer = 0; layer < 4; layer++) {
            std::vector<uint32_t> pairs;
            for (const Plaquette& q : plaquettes) {
                uint32_t target = q.order[layer];
                if (target == UINT32_MAX)
                    continue;
                if (q.is_x) {
                    pairs.push_back(q.ancilla);
                    pairs.push_back(target);
                } else {
                    pairs.push_back(target);
                    pairs.push_back(q.ancilla);
                }
            }
            op("CX", pairs);
            op("DEPOLARIZE2", pairs, noise);
        }
        op("H", x_ancillas);
        op("DEPOLARIZE1", x_ancillas, noise);
        op("X_ERROR", ancillas, noise);
        op("M", ancillas);
        for (size_t k = 0; k < plaquettes.size(); k++)
            stab_meas[r][k] = measurement_count + k;
        measurement_count += plaquettes.size();
        c.append_from_text("TICK\n");
    }
    op("Z_ERROR", data, noise);
    op("MX", data);
    for (size_t k = 0; k < data.size(); k++)
        data_meas[k] = measurement_count + k;
    measurement_count += data.size();

    auto rec = [&](size_t absolute) -> std::string {
        return "rec[-" + std::to_string(measurement_count - absolute) + "]";
    };

    // Detectors, cube-major.
    out.cube_offsets.assign(p.patches + 1, 0);
    for (uint32_t cube = 0; cube < (uint32_t)p.patches; cube++) {
        out.cube_offsets[cube] = (uint32_t)out.detector_cube.size();
        for (size_t r = 0; r <= p.rounds; r++) {
            for (size_t k = 0; k < plaquettes.size(); k++) {
                const Plaquette& q = plaquettes[k];
                if (!q.is_x || q.cube != cube)
                    continue;
                std::ostringstream text;
                text << "DETECTOR(" << 2 * q.i + 2 << ", " << 2 * q.j + 2 << ", " << r << ")";
                if (r == 0) {
                    if (p.merge_init && q.seam) {
                        out.removed_first_round_detectors++;
                        continue;
                    }
                    text << " " << rec(stab_meas[0][k]);
                } else if (r < p.rounds) {
                    text << " " << rec(stab_meas[r][k]) << " " << rec(stab_meas[r - 1][k]);
                } else {
                    // Final round: the product of the data X measurements on the support, against
                    // the last stabilizer measurement.
                    text << " " << rec(stab_meas[p.rounds - 1][k]);
                    for (uint32_t d : q.support)
                        text << " " << rec(data_meas[d]);
                }
                text << "\n";
                c.append_from_text(text.str());
                out.detector_cube.push_back(cube);
            }
        }
    }
    out.cube_offsets[p.patches] = (uint32_t)out.detector_cube.size();

    // Observable 0: the logical X string.
    {
        std::ostringstream text;
        text << "OBSERVABLE_INCLUDE(0)";
        if (p.long_orientation) {
            for (long i = 0; i < W; i++)
                text << " " << rec(data_meas[data_index(i, 0)]);
        } else {
            for (long j = 0; j < H; j++)
                text << " " << rec(data_meas[data_index(0, j)]);
        }
        text << "\n";
        c.append_from_text(text.str());
    }
    out.num_observables = 1;
    for (const Plaquette& q : plaquettes)
        out.seam_plaquettes += (q.is_x && q.seam) ? 1 : 0;
    if (p.merge_init && out.seam_plaquettes > 0) {
        std::ostringstream text;
        text << "OBSERVABLE_INCLUDE(1)";
        for (size_t k = 0; k < plaquettes.size(); k++) {
            if (plaquettes[k].is_x && plaquettes[k].seam)
                text << " " << rec(stab_meas[0][k]);
        }
        text << "\n";
        c.append_from_text(text.str());
        out.num_observables = 2;
    }
    out.num_detectors = out.detector_cube.size();
    if (out.circuit.count_detectors() != out.num_detectors)
        throw std::logic_error("lattice surgery circuit: detector count mismatch");
    return out;
}

}  // namespace profiler
}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_BENCHMARKS_SPEC_MATCHING_LATTICE_SURGERY_CIRCUIT_H
