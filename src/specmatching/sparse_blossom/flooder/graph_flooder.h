// Copyright 2022 PyMatching Contributors
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

#ifndef SPECMATCHING2_GRAPH_FLOODER_H
#define SPECMATCHING2_GRAPH_FLOODER_H

#include <queue>

#include "specmatching/sparse_blossom/arena.h"
#include "specmatching/sparse_blossom/flooder/graph.h"
#include "specmatching/sparse_blossom/flooder/graph_fill_region.h"
#include "specmatching/sparse_blossom/flooder_matcher_interop/mwpm_event.h"
#include "specmatching/sparse_blossom/flooder_matcher_interop/region_edge.h"
#include "specmatching/sparse_blossom/tracker/flood_check_event.h"
#include "specmatching/sparse_blossom/tracker/radix_heap_queue.h"

namespace pm {

/// Counts blossom formations, for §M2.9.6 measurement 3.
///
/// §M2.9.4's eager-cached-base option is worth a field and a mutation site only if formations are
/// much rarer than base readouts, and that ratio is not observable after the fact: a blossom that
/// forms and then shatters again leaves nothing behind to count. Hence a counter at the formation
/// site. Unlike `pm::horizon_gate_stats` this one is live in release builds too, because the
/// measurement it feeds is a release measurement; the cost is one increment per blossom, against a
/// formation that already sweeps the blossom's whole area to reschedule events. Not thread safe.
struct BlossomFormationStats {
    uint64_t formations{0};

    void clear() {
        formations = 0;
    }
};

/// Process-wide counter for blossom formations. Not thread safe.
inline BlossomFormationStats blossom_formation_stats{};

struct GraphFlooder {
    /// The graph of detector nodes that is being flooded.
    MatchingGraph graph;
    /// Tracks the next thing that will occur as flooding proceeds.
    /// The events are "tentative" because processing an event may remove another,
    /// for example if a region stops growing due to colliding with another region
    /// then the growing region will no longer reach other nodes even if those
    /// events were scheduled in this queue before the region collision was processed.
    ///
    /// Events are ordered by time; by when they will occur in a timeline.
    radix_heap_queue<false> queue;

    Arena<GraphFillRegion> region_arena;

    std::vector<CompressedEdge> match_edges;

    /// These are the detection events that would occur if an error occurred on every edge that has a negative weight.
    /// Stored as a sorted vector of indices of detection events.
    std::vector<uint64_t> negative_weight_detection_events;
    /// These are the observables that would be flipped if an error occurred on every edge that has a negative weight.
    /// Stored as a sorted vector of indices of observables.
    std::vector<size_t> negative_weight_observables;
    /// Observable mask corresponding to the observables that would be flipped if an error occurred on every edge that
    /// has a negative weight. Only used for fewer than 64 (=sizeof(pm::obs_int)*8) observables.
    pm::obs_int negative_weight_obs_mask;
    /// The sum of the edge weights of all edges with negative edge weights.
    pm::total_weight_int negative_weight_sum;

    /// Truncation horizon `T`, in cumulative time units (M1).
    ///
    /// Node/grow events that would occur strictly after this time are never inserted into the
    /// queue, and the timeline loop stops rather than advancing past it. The sentinel
    /// `pm::NO_HORIZON` disables truncation entirely and is the only value the stock decode path
    /// ever sees; it is set (and restored) by `process_timeline_until_horizon`.
    cumulative_time_int horizon{NO_HORIZON};

    /// §4.1 of the fusion design — the edge mask of the divide-and-conquer solve.
    ///
    /// One byte per *directed* adjacency entry, laid out CSR-style: node `i`'s entries start at
    /// `edge_mask_offsets[i]` and run parallel to `nodes[i].neighbors`. A non-zero byte hides that
    /// half-edge from the neighbour scan, so no collision event is ever generated across it and no
    /// region can grow through it. Both directions of an undirected edge must carry the same bit;
    /// nothing here enforces that, because the only writer is the fusion driver, which owns both.
    ///
    /// `nullptr` — the default, and the only value the stock decode path ever sees — makes the scan
    /// the code it was before this field existed: the two scan bodies are templated on whether a
    /// mask is present and the choice is made **once per node event**, not once per neighbour, so an
    /// unmasked instance pays one perfectly-predicted branch per event and no extra loads at all.
    /// `DetectorNode` is not widened by a byte, which is why the mask lives here and not there.
    const uint8_t* edge_mask{nullptr};
    const uint32_t* edge_mask_offsets{nullptr};

    /// §4.2 — the per-region dual cap `T`, in cumulative time units, or `NO_HORIZON` for "no cap".
    ///
    /// The global `horizon` above bounds the *clock*; on a monolithic run started at time 0 that
    /// bounds every dual too, because every dual is then at most the elapsed time. The fusion driver
    /// solves many pieces one after another on one instance, so its clock has long since left any
    /// single defect's dual behind and the bound has to be checked where it lives: on each growing
    /// region. A region is capped at `radius == dual_cap - inner_max`, where `inner_max` is the
    /// largest frozen (blossom-nested) part among its member defects; reaching that means some
    /// member's dual is exactly `dual_cap`, which is legal, and the event is scheduled for the first
    /// instant *past* it, which is not.
    ///
    /// The cap event rides the region's `shrink_event_tracker`, which is idle for the whole time a
    /// region is growing — `set_region_growing` used to call `set_no_desired_event()` on it and now
    /// schedules the cap there instead. That is what keeps `GraphFillRegion` at 128 bytes.
    cumulative_time_int dual_cap{NO_HORIZON};
    /// Set when a growing region's cap event fired, i.e. some defect's dual would have exceeded
    /// `dual_cap`. `run_until_next_mwpm_notification` stops and reports `NO_EVENT` as soon as it is
    /// set; the driver reads it, reports `TRUNCATED` and clears it. Never set while
    /// `dual_cap == NO_HORIZON`.
    bool dual_cap_hit{false};

    /// `lattice_surgery_profiler` §7.5 — the `(t1)` trigger of the phase-2 hybrid resume.
    ///
    /// One byte per node, indexed like `graph.nodes`: non-zero marks a node on the outer layer of the
    /// full-resolution territory `Ω`, i.e. a node of `Ω` with a `G`-edge leaving it. The moment a
    /// region claims such a node — by arriving at it, or by being created on it — `outer_layer_hit`
    /// is set to that node and `run_until_next_mwpm_notification` stops exactly as it does for a
    /// fired dual cap, so the driver can grow `Ω` before any part of an edge leaving it is covered.
    ///
    /// `nullptr` — the default, and the only value the stock decode path ever sees — costs one
    /// predicted branch per claimed node and nothing else.
    const uint8_t* outer_layer{nullptr};
    /// The outer-layer node whose claim stopped the timeline, or `nullptr`. Never set while
    /// `outer_layer == nullptr`; the driver reads it and clears it.
    DetectorNode* outer_layer_hit{nullptr};

    GraphFlooder();
    explicit GraphFlooder(MatchingGraph graph);
    GraphFlooder(GraphFlooder&&) noexcept;
    MwpmEvent run_until_next_mwpm_notification();
    void set_region_growing(pm::GraphFillRegion& region);
    void set_region_frozen(pm::GraphFillRegion& region);
    void set_region_shrinking(pm::GraphFillRegion& region);
    GraphFillRegion* create_blossom(std::vector<RegionEdge>& contained_regions);
    void schedule_tentative_shrink_event(GraphFillRegion& region);
    /// §4.2. Puts `region`'s dual-cap event into the queue, given the largest frozen part among its
    /// member defects. Only legal while `region.radius` is growing and `dual_cap != NO_HORIZON`.
    /// Exposed because the fusion driver creates detection events at a non-zero clock itself and has
    /// to cap them exactly as `set_region_growing` would.
    void schedule_dual_cap_event(GraphFillRegion& region, cumulative_time_int inner_max);
    /// §4.2. Reports a growing region that is about to stop growing while its own dual-cap event is
    /// already due, which means its radius is one tick past the cap.
    ///
    /// `schedule_dual_cap_event` puts the cap at the first instant the region would grow *past* it,
    /// so that a collision at a dual of exactly `dual_cap` — which is legal, and is how a `2T` edge
    /// goes tight — still runs. The queue gives no order among events that share a timestamp,
    /// though, so an event at the cap's own time can be dequeued first; if it matches or blossoms
    /// the region, `set_region_frozen` clears the tracker slot the cap was riding and the cap is
    /// lost. The timeline then finishes with that region's defects one unit over `dual_cap`. This is
    /// the tie-break: a cap that was due when growth ended counts as fired.
    void note_growth_ending(GraphFillRegion& region);
    /// `lattice_surgery_profiler` §7.5: records a claim of `node` if it is on the outer layer.
    inline void note_node_claimed(DetectorNode& node) {
        if (outer_layer != nullptr && outer_layer[(size_t)(&node - graph.nodes.data())])
            outer_layer_hit = &node;
    }
    void reschedule_events_at_detector_node(DetectorNode& detector_node);
    void do_region_created_at_empty_detector_node(GraphFillRegion& region, DetectorNode& detector_node);
    void do_region_arriving_at_empty_detector_node(
        GraphFillRegion& region, DetectorNode& empty_node, const DetectorNode& from_node, size_t from_to_empty_index);
    MwpmEvent do_region_shrinking(GraphFillRegion& shrinking_region);
    pm::MwpmEvent do_neighbor_interaction(DetectorNode& src, size_t src_to_dst_index, DetectorNode& dst);
    pm::MwpmEvent do_region_hit_boundary_interaction(DetectorNode& node);
    static MwpmEvent do_degenerate_implosion(const GraphFillRegion& region);
    static MwpmEvent do_blossom_shattering(GraphFillRegion& region);
    bool dequeue_decision(pm::FloodCheckEvent ev);
    std::pair<size_t, pm::cumulative_time_int> find_next_event_at_node_returning_neighbor_index_and_time(
        const DetectorNode& detector_node) const;
    pm::MwpmEvent do_look_at_node_event(DetectorNode& node);

    pm::FloodCheckEvent dequeue_valid();
    pm::MwpmEvent process_tentative_event_returning_mwpm_event(FloodCheckEvent tentative_event);

    void sync_negative_weight_observables_and_detection_events();
};

}  // namespace pm

#endif  // SPECMATCHING2_GRAPH_FLOODER_H
