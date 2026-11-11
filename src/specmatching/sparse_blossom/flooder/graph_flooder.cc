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

#include "specmatching/sparse_blossom/flooder/graph_flooder.h"

#include <algorithm>

#include "specmatching/sparse_blossom/flooder/graph.h"
#include "specmatching/sparse_blossom/flooder/graph_fill_region.h"
#include "specmatching/sparse_blossom/flooder_matcher_interop/varying.h"
#include "specmatching/sparse_blossom/matcher/alternating_tree.h"

using namespace pm;

GraphFlooder::GraphFlooder(MatchingGraph graph)
    : graph(std::move(graph)), negative_weight_obs_mask(0), negative_weight_sum(0) {
}

GraphFlooder::GraphFlooder(GraphFlooder &&flooder) noexcept
    : graph(std::move(flooder.graph)),
      queue(std::move(flooder.queue)),
      region_arena(std::move(flooder.region_arena)),
      match_edges(std::move(flooder.match_edges)),
      negative_weight_detection_events(std::move(flooder.negative_weight_detection_events)),
      negative_weight_observables(std::move(flooder.negative_weight_observables)),
      negative_weight_obs_mask(flooder.negative_weight_obs_mask),
      negative_weight_sum(flooder.negative_weight_sum),
      horizon(flooder.horizon),
      edge_mask(flooder.edge_mask),
      edge_mask_offsets(flooder.edge_mask_offsets),
      dual_cap(flooder.dual_cap),
      dual_cap_hit(flooder.dual_cap_hit),
      outer_layer(flooder.outer_layer),
      outer_layer_hit(flooder.outer_layer_hit) {
}

void GraphFlooder::do_region_created_at_empty_detector_node(GraphFillRegion &region, DetectorNode &detector_node) {
    detector_node.reached_from_source = &detector_node;
    detector_node.radius_of_arrival = 0;
    detector_node.region_that_arrived = &region;
    detector_node.region_that_arrived_top = &region;
    detector_node.wrapped_radius_cached = 0;
    region.shell_area.push_back(&detector_node);
    reschedule_events_at_detector_node(detector_node);
    note_node_claimed(detector_node);
}

/// `Masked` selects between the stock scan and the §4.1 masked one. With `Masked == false` the
/// compiler sees the original body verbatim — `mask` is never named — so an unmasked instance pays
/// nothing per neighbour.
template <bool Masked>
std::pair<size_t, cumulative_time_int> find_next_event_at_node_not_occupied_by_growing_top_region(
    const DetectorNode &detector_node, VaryingCT rad1, const uint8_t *mask) {
    cumulative_time_int best_time = std::numeric_limits<cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;

    size_t start = 0;
    if (!detector_node.neighbors.empty() && detector_node.neighbors[0] == nullptr)
        start++;

    // Handle non-boundary neighbors.
    for (size_t i = start; i < detector_node.neighbors.size(); i++) {
        if constexpr (Masked) {
            // §4.1: a masked half-edge is invisible. No collision across it, and no growth through
            // it, which is what makes a cut of `H` a virtual boundary rather than a real one.
            if (mask[i])
                continue;
        }
        auto weight = detector_node.neighbor_weights[i];

        auto neighbor = detector_node.neighbors[i];

        auto rad2 = neighbor->local_radius();

        if (rad2.is_growing()) {
            auto collision_time = weight - rad1.y_intercept() - rad2.y_intercept();
            if (collision_time < best_time) {
                best_time = collision_time;
                best_neighbor = i;
            }
        }
    }
    return {best_neighbor, best_time};
}

template <bool Masked>
std::pair<size_t, cumulative_time_int> find_next_event_at_node_occupied_by_growing_top_region(
    const DetectorNode &detector_node, const VaryingCT &rad1, const uint8_t *mask) {
    cumulative_time_int best_time = std::numeric_limits<cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;
    size_t start = 0;
    if (!detector_node.neighbors.empty() && detector_node.neighbors[0] == nullptr) {
        // Growing towards boundary. The fusion driver never masks a boundary half-edge — a piece's
        // dummy boundary *is* a boundary edge, and the cut it stands for is carried by the masked
        // defect-defect edges instead — but `lattice_surgery_profiler`'s fallback instance gives
        // every node a boundary slot and masks the ones that hold no boundary this shot, so a masked
        // slot is skipped here like any other masked half-edge.
        bool visible = true;
        if constexpr (Masked)
            visible = mask[0] == 0;
        if (visible) {
            auto weight = detector_node.neighbor_weights[0];
            auto collision_time = weight - rad1.y_intercept();
            if (collision_time < best_time) {
                best_time = collision_time;
                best_neighbor = 0;
            }
        }
        start++;
    }

    // Handle non-boundary neighbors.
    for (size_t i = start; i < detector_node.neighbors.size(); i++) {
        if constexpr (Masked) {
            if (mask[i])
                continue;
        }
        auto weight = detector_node.neighbor_weights[i];

        auto neighbor = detector_node.neighbors[i];
        if (detector_node.has_same_owner_as(*neighbor)) {
            continue;
        }
        auto rad2 = neighbor->local_radius();
        if (rad2.is_shrinking()) {
            continue;
        }

        auto collision_time = weight - rad1.y_intercept() - rad2.y_intercept();
        if (rad2.is_growing()) {
            collision_time >>= 1;
        }
        if (collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = i;
        }
    }
    return {best_neighbor, best_time};
}

std::pair<size_t, cumulative_time_int> GraphFlooder::find_next_event_at_node_returning_neighbor_index_and_time(
    const DetectorNode &detector_node) const {
    auto rad1 = detector_node.local_radius();

    // One branch per node event, taken on a field that is `nullptr` for the whole lifetime of every
    // instance but the fusion driver's. See `GraphFlooder::edge_mask`.
    if (edge_mask == nullptr) {
        if (rad1.is_growing()) {
            return find_next_event_at_node_occupied_by_growing_top_region<false>(detector_node, rad1, nullptr);
        } else {
            return find_next_event_at_node_not_occupied_by_growing_top_region<false>(detector_node, rad1, nullptr);
        }
    }
    const uint8_t *mask = edge_mask + edge_mask_offsets[(size_t)(&detector_node - graph.nodes.data())];
    if (rad1.is_growing()) {
        return find_next_event_at_node_occupied_by_growing_top_region<true>(detector_node, rad1, mask);
    } else {
        return find_next_event_at_node_not_occupied_by_growing_top_region<true>(detector_node, rad1, mask);
    }
}

void GraphFlooder::reschedule_events_at_detector_node(DetectorNode &detector_node) {
    auto x = find_next_event_at_node_returning_neighbor_index_and_time(detector_node);
    if (x.first == SIZE_MAX) {
        detector_node.node_event_tracker.set_no_desired_event();
    } else {
        // GATED: this is a growth/collision event. If it lies past the horizon it will never be
        // reached within the truncated timeline, so it must not enter the queue at all.
        detector_node.node_event_tracker.set_desired_event(
            {
                &detector_node,
                cyclic_time_int{x.second},
            },
            queue,
            x.second,
            horizon);
    }
}

void GraphFlooder::schedule_tentative_shrink_event(GraphFillRegion &region) {
    cumulative_time_int t;
    if (region.shell_area.empty()) {
        t = region.radius.time_of_x_intercept_for_shrinking();
    } else {
        t = region.shell_area.back()->local_radius().time_of_x_intercept_for_shrinking();
    }
    // EXEMPT: shrink events are never gated. A shrinking region's shell-area bookkeeping has to be
    // able to complete, or the region state is inconsistent at harvest time.
    region.shrink_event_tracker.set_desired_event(
        {
            &region,
            cyclic_time_int{t},
        },
        queue);
}

void GraphFlooder::schedule_dual_cap_event(GraphFillRegion &region, cumulative_time_int inner_max) {
    // While growing, `radius(t) = t + y_intercept`, so the radius reaches its cap `dual_cap -
    // inner_max` at `t_cap`. Reaching it is legal — a dual of exactly `T` is what makes a `2T` edge
    // tight, and `H` has those — so the event goes one tick later, at the first instant the region
    // would grow *past* the cap. Everything scheduled at `t_cap` therefore still runs.
    cumulative_time_int t_cap = (dual_cap - inner_max) - region.radius.y_intercept() + 1;
    // A region released at a radius already at its cap is capped immediately rather than in the
    // past; the queue is monotone and will not take an event cycle-before the present.
    if (t_cap < queue.cur_time)
        t_cap = queue.cur_time;
    // EXEMPT from the horizon gate, like the shrink event whose tracker slot this borrows: the cap
    // *is* this instance's truncation test, and gating it against a horizon it replaces would be
    // circular.
    region.shrink_event_tracker.set_desired_event(
        {
            &region,
            cyclic_time_int{t_cap},
        },
        queue);
}

void GraphFlooder::do_region_arriving_at_empty_detector_node(
    GraphFillRegion &region, DetectorNode &empty_node, const DetectorNode &from_node, size_t from_to_empty_index) {
    empty_node.observables_crossed_from_source =
        (from_node.observables_crossed_from_source ^ from_node.neighbor_observables[from_to_empty_index]);
    empty_node.reached_from_source = from_node.reached_from_source;
    empty_node.radius_of_arrival = region.radius.get_distance_at_time(queue.cur_time);
    empty_node.region_that_arrived = &region;
    empty_node.region_that_arrived_top = region.blossom_parent_top;
    empty_node.wrapped_radius_cached = empty_node.compute_wrapped_radius();
    region.shell_area.push_back(&empty_node);
    reschedule_events_at_detector_node(empty_node);
    note_node_claimed(empty_node);
}

MwpmEvent GraphFlooder::do_region_shrinking(GraphFillRegion &region) {
    // §4.2. A region that is not shrinking has no business having a shrink event: while it grows,
    // its tracker slot holds the dual-cap event instead. So this *is* the cap firing — one of its
    // member defects would now have a dual past `dual_cap`. Stop, and let the driver escalate.
    if (dual_cap != NO_HORIZON && !region.radius.is_shrinking()) {
        dual_cap_hit = true;
        return MwpmEvent::no_event();
    }
    if (region.shell_area.empty()) {
        return do_blossom_shattering(region);
    } else if (region.shell_area.size() == 1 && region.blossom_children.empty()) {
        return do_degenerate_implosion(region);
    } else {
        auto leaving_node = region.shell_area.back();
        region.shell_area.pop_back();
        leaving_node->region_that_arrived = nullptr;
        leaving_node->region_that_arrived_top = nullptr;
        leaving_node->wrapped_radius_cached = 0;
        leaving_node->reached_from_source = nullptr;
        leaving_node->radius_of_arrival = 0;
        leaving_node->observables_crossed_from_source = 0;
        reschedule_events_at_detector_node(*leaving_node);
        schedule_tentative_shrink_event(region);
        return MwpmEvent::no_event();
    }
}

MwpmEvent GraphFlooder::do_neighbor_interaction(DetectorNode &src, size_t src_to_dst_index, DetectorNode &dst) {
    // First check if one region is moving into an empty location
    if (src.region_that_arrived && !dst.region_that_arrived) {
        do_region_arriving_at_empty_detector_node(*src.region_that_arrived_top, dst, src, src_to_dst_index);
        return MwpmEvent::no_event();
    } else if (dst.region_that_arrived && !src.region_that_arrived) {
        do_region_arriving_at_empty_detector_node(*dst.region_that_arrived_top, src, dst, dst.index_of_neighbor(&src));
        return MwpmEvent::no_event();
    } else {
        // Two regions colliding
        return RegionHitRegionEventData{
            src.region_that_arrived_top,
            dst.region_that_arrived_top,
            CompressedEdge{
                src.reached_from_source,
                dst.reached_from_source,
                src.observables_crossed_from_source ^ dst.observables_crossed_from_source ^
                    src.neighbor_observables[src_to_dst_index]},
        };
    }
}

MwpmEvent GraphFlooder::do_region_hit_boundary_interaction(DetectorNode &node) {
    return RegionHitBoundaryEventData{
        node.region_that_arrived_top,
        CompressedEdge{
            node.reached_from_source, nullptr, node.observables_crossed_from_source ^ node.neighbor_observables[0]}};
}

MwpmEvent GraphFlooder::do_degenerate_implosion(const GraphFillRegion &region) {
    return RegionHitRegionEventData{
        region.alt_tree_node->parent.alt_tree_node->outer_region,
        region.alt_tree_node->outer_region,
        CompressedEdge{
            region.alt_tree_node->parent.edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.obs_mask ^ region.alt_tree_node->parent.edge.obs_mask}};
}

MwpmEvent GraphFlooder::do_blossom_shattering(GraphFillRegion &region) {
    return BlossomShatterEventData{
        &region,
        region.alt_tree_node->parent.edge.loc_from->heir_region_on_shatter(),
        region.alt_tree_node->inner_to_outer_edge.loc_from->heir_region_on_shatter()};
}

GraphFillRegion *GraphFlooder::create_blossom(std::vector<RegionEdge> &contained_regions) {
    // §M2.9.6 measurement 3. The one place a blossom comes into existence.
    blossom_formation_stats.formations++;
    auto blossom_region = region_arena.alloc_default_constructed();
    blossom_region->radius = VaryingCT::growing_varying_with_zero_distance_at_time(queue.cur_time);
    blossom_region->blossom_children = std::move(contained_regions);
    for (auto &region_edge : blossom_region->blossom_children) {
        // Freezing a child here bypasses `set_region_frozen`, so the cap tie-break is taken by hand.
        note_growth_ending(*region_edge.region);
        region_edge.region->radius = region_edge.region->radius.then_frozen_at_time(queue.cur_time);
        region_edge.region->wrap_into_blossom(blossom_region);
        region_edge.region->shrink_event_tracker.set_no_desired_event();
    }

    // §4.2: the blossom's cap needs the largest frozen part among its member defects, and that is
    // exactly `max wrapped_radius_cached` over the defects of its area — the children's radii were
    // just frozen into the wrapped radii above. The sweep the formation already does is where it is
    // taken from, so the cap costs a comparison per node and no second traversal.
    cumulative_time_int inner_max = 0;
    bool capped = dual_cap != NO_HORIZON;
    blossom_region->do_op_for_each_node_in_total_area([this, capped, &inner_max](DetectorNode *n) {
        if (capped && n->reached_from_source == n)
            inner_max = std::max(inner_max, (cumulative_time_int)n->wrapped_radius_cached);
        reschedule_events_at_detector_node(*n);
    });
    if (capped)
        schedule_dual_cap_event(*blossom_region, inner_max);

    return blossom_region;
}

bool GraphFlooder::dequeue_decision(FloodCheckEvent ev) {
    switch (ev.tentative_event_type) {
        case FloodCheckEventType::LOOK_AT_NODE: {
            auto &node = *ev.data_look_at_node;
            return node.node_event_tracker.dequeue_decision(ev, queue);
        }
        case FloodCheckEventType::LOOK_AT_SHRINKING_REGION: {
            auto &region = *ev.data_look_at_shrinking_region;
            return region.shrink_event_tracker.dequeue_decision(ev, queue);
        }
        case FloodCheckEventType::NO_FLOOD_CHECK_EVENT:
            return true;
        default:
            throw std::invalid_argument("Unrecognized event type.");
    }
}

FloodCheckEvent GraphFlooder::dequeue_valid() {
    while (true) {
        FloodCheckEvent ev = queue.dequeue();
        if (dequeue_decision(ev)) {
            return ev;
        }
    }
}

void GraphFlooder::set_region_growing(GraphFillRegion &region) {
    region.radius = region.radius.then_growing_at_time(queue.cur_time);

    // No shrinking event can occur while growing.
    region.shrink_event_tracker.set_no_desired_event();

    // Node events can occur while growing, and events in the queue may occur sooner than
    // previously scheduled. Therefore, we must reschedule all the nodes.
    //
    // §4.2 rides along: a growing region's `shrink_event_tracker` slot, cleared just above, now
    // carries its dual cap instead, and the cap's `inner_max` is read off the same sweep.
    cumulative_time_int inner_max = 0;
    bool capped = dual_cap != NO_HORIZON;
    region.do_op_for_each_node_in_total_area([this, capped, &inner_max](DetectorNode *n) {
        if (capped && n->reached_from_source == n)
            inner_max = std::max(inner_max, (cumulative_time_int)n->wrapped_radius_cached);
        reschedule_events_at_detector_node(*n);
    });
    if (capped)
        schedule_dual_cap_event(region, inner_max);
}

void GraphFlooder::note_growth_ending(GraphFillRegion &region) {
    if (dual_cap == NO_HORIZON || !region.radius.is_growing())
        return;
    if (!region.shrink_event_tracker.has_desired_time)
        return;
    // Monotone queue: `cur_time` never passes an event that is still in it, so "due" is equality.
    // The comparison is cyclic because that is the form the tracker stores.
    if (!(region.shrink_event_tracker.desired_time == cyclic_time_int{queue.cur_time}))
        return;
    dual_cap_hit = true;
}

void GraphFlooder::set_region_frozen(GraphFillRegion &region) {
    note_growth_ending(region);
    bool was_shrinking = region.radius.is_shrinking();
    region.radius = region.radius.then_frozen_at_time(queue.cur_time);

    // No shrinking event can occur while frozen.
    region.shrink_event_tracker.set_no_desired_event();

    // Node events can occur while frozen, from other regions growing into this one.
    // However, those events can only be sooner than the currently scheduled events
    // if the region was previously shrinking (as opposed to growing).
    if (was_shrinking) {
        region.do_op_for_each_node_in_total_area([&](DetectorNode *n) {
            reschedule_events_at_detector_node(*n);
        });
    }
}

void GraphFlooder::set_region_shrinking(GraphFillRegion &region) {
    note_growth_ending(region);
    region.radius = region.radius.then_shrinking_at_time(queue.cur_time);

    // Shrinking events can now occur.
    schedule_tentative_shrink_event(region);

    // No node events can occur while shrinking.
    region.do_op_for_each_node_in_total_area([&](DetectorNode *n) {
        n->node_event_tracker.set_no_desired_event();
    });
}

MwpmEvent GraphFlooder::do_look_at_node_event(DetectorNode &node) {
    auto next = find_next_event_at_node_returning_neighbor_index_and_time(node);
    if (next.second == queue.cur_time) {
        // Need to revisit this node immediately after the mwpm event is handled. There may be an event to handle
        // along another edge, or even along the same edge at a later time (e.g. if this event isn't the *first* event
        // for the neighbor along the current edge when the neighbor is being rescheduled).
        // GATED: a node event, like every other. The gate never actually fires here, because the
        // event is scheduled for the present and the present is never past the horizon; the
        // argument is passed anyway so that no node-event scheduling site is un-gated.
        node.node_event_tracker.set_desired_event(
            {
                &node,
                cyclic_time_int{queue.cur_time},
            },
            queue,
            queue.cur_time,
            horizon);

        if (node.neighbors[next.first] == nullptr) {
            return do_region_hit_boundary_interaction(node);
        }
        auto &neighbor = *node.neighbors[next.first];
        return do_neighbor_interaction(node, next.first, neighbor);
    } else if (next.first != SIZE_MAX) {
        // Need to revisit this node at a later time.
        // GATED: a growth/collision event past the horizon is unreachable in the truncated
        // timeline, so it must not enter the queue at all.
        node.node_event_tracker.set_desired_event(
            {
                &node,
                cyclic_time_int{next.second},
            },
            queue,
            next.second,
            horizon);
    }

    return MwpmEvent::no_event();
}

MwpmEvent GraphFlooder::process_tentative_event_returning_mwpm_event(FloodCheckEvent tentative_event) {
    switch (tentative_event.tentative_event_type) {
        case LOOK_AT_NODE: {
            return do_look_at_node_event(*tentative_event.data_look_at_node);
        }
        case LOOK_AT_SHRINKING_REGION: {
            return do_region_shrinking(*tentative_event.data_look_at_shrinking_region);
        }
        default:
            throw std::invalid_argument("Unknown tentative event type.");
    }
}

MwpmEvent GraphFlooder::run_until_next_mwpm_notification() {
    while (true) {
        // A cap that `note_growth_ending` recorded during the caller's `process_event` is read here,
        // before another event is dequeued: the timeline is already past `dual_cap` and everything
        // after this point would be solving a different problem. An outer-layer claim stops it the
        // same way; `outer_layer_hit` is `nullptr` on every instance but the lattice-surgery
        // fallback's.
        if (dual_cap_hit || outer_layer_hit != nullptr) {
            return MwpmEvent::no_event();
        }
        FloodCheckEvent tentative_event = dequeue_valid();
        if (tentative_event.tentative_event_type == NO_FLOOD_CHECK_EVENT) {
            return MwpmEvent::no_event();
        }
        // Truncation stop (M1). Dequeueing advances `queue.cur_time` to the event's time. Node
        // events past the horizon were never enqueued, but shrink events are exempt from the gate,
        // so a surviving shrink event is the one thing that can carry the timeline past `T`.
        // Report "nothing further happens" instead of processing it. With the sentinel horizon
        // this comparison is never true and the stock path is unchanged.
        if (queue.cur_time > horizon) {
            return MwpmEvent::no_event();
        }
        MwpmEvent notification = process_tentative_event_returning_mwpm_event(tentative_event);
        // §4.2's truncation stop. A fired dual cap means some defect's dual would pass `dual_cap`,
        // so the timeline must not go on: report "nothing further happens" and let the driver
        // read the flag. With `dual_cap == NO_HORIZON` this is never set and the stock path is
        // unchanged. The `(t1)` trigger likewise: the claim has been made, nothing past it has.
        if (dual_cap_hit || outer_layer_hit != nullptr) {
            return MwpmEvent::no_event();
        }
        if (notification.event_type != NO_EVENT) {
            return notification;
        }
    }
}

void GraphFlooder::sync_negative_weight_observables_and_detection_events() {
    /// Move set of negative weight detection events into a sorted vector, for faster processing during decoding
    negative_weight_detection_events.clear();
    negative_weight_detection_events.reserve(graph.negative_weight_detection_events_set.size());
    for (auto &det : graph.negative_weight_detection_events_set)
        negative_weight_detection_events.push_back(det);
    /// Move set of negative weight observable indices into a sorted vector, for faster processing during decoding
    negative_weight_observables.clear();
    negative_weight_observables.reserve(graph.negative_weight_observables_set.size());
    for (auto &obs : graph.negative_weight_observables_set)
        negative_weight_observables.push_back(obs);
    negative_weight_obs_mask = 0;
    if (graph.num_observables <= sizeof(pm::obs_int) * 8) {
        // Compute the observable mask for 64 or fewer observables
        for (auto &i : negative_weight_observables)
            negative_weight_obs_mask ^= (pm::obs_int)1 << i;
    }
    negative_weight_sum = graph.negative_weight_sum;
}

GraphFlooder::GraphFlooder() : negative_weight_obs_mask(0), negative_weight_sum(0) {
}
