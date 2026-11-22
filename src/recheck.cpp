#include "recheck.h"
#include "engine_state.h"
#include "message_updates.h"
#include "visits.h"
namespace pap_internal::search {
int run_recheck_wave(RecheckContext& context, int round) {
  auto& complete_recheck_enabled = context.complete_recheck_enabled;
  auto& pending_recheck_count = context.pending_recheck_count;
  auto& graph_revision = context.graph_revision;
  auto& parents = context.parents;
  auto& expected_year = context.expected_year;
  auto& id_order_rank = context.id_order_rank;
  auto& inputs = context.inputs;
  auto& final_focal_order = context.final_focal_order;
  auto& final_active_input = context.final_active_input;
  auto& dam_candidates = context.dam_candidates;
  auto& sire_candidates = context.sire_candidates;
  auto& message_revision = context.message_revision;
  auto& compiled_sexp = context.compiled_sexp;
  auto& record_typed_ages = context.record_typed_ages;
  auto& record_complete_evaluations = context.record_complete_evaluations;
  auto& marking_round = context.marking_round;
  auto& mark_complete_rechecks_stale = context.mark_complete_rechecks_stale;
  auto& defer_intermediate_publication = context.defer_intermediate_publication;
  auto& deferred_publication_phases = context.deferred_publication_phases;
  auto& publish = context.publish;
  auto& publication = context.publication;
  auto& recheck_pending = context.recheck_pending;
  auto& last_result = context.last_result;
  auto& recheck_move_blocks = context.recheck_move_blocks;
  auto& recheck_focal_blocks = context.recheck_focal_blocks;
  auto& recheck_configuration_diagnostic_blocks = context.recheck_configuration_diagnostic_blocks;
  auto& recheck_graph_state_diagnostic_blocks = context.recheck_graph_state_diagnostic_blocks;
  auto& recheck_relationship_provenance_blocks = context.recheck_relationship_provenance_blocks;
  auto& recheck_phase_round = context.recheck_phase_round;
  auto& recheck_phase_targets = context.recheck_phase_targets;
  auto& recheck_phase_changes = context.recheck_phase_changes;
  auto& recheck_phase_pending_after = context.recheck_phase_pending_after;
  auto& recheck_phase_seconds = context.recheck_phase_seconds;
  auto& recheck_phase_setup_seconds = context.recheck_phase_setup_seconds;
  auto& recheck_phase_pass_seconds = context.recheck_phase_pass_seconds;
  auto& recheck_phase_publication_seconds = context.recheck_phase_publication_seconds;
    if (!complete_recheck_enabled || pending_recheck_count() == 0) return 0;
    const auto started = std::chrono::steady_clock::now();
    const int graph_before_recheck = graph_revision;
    IntegerVector order = focal_order(parents, expected_year, id_order_rank);
    std::vector<int> targets;
    for (int at = 0; at < order.size(); ++at)
      if (recheck_pending[static_cast<std::size_t>(order[at])] &&
          parents(order[at] - 1, 0) > 0 && parents(order[at] - 1, 1) > 0)
        targets.push_back(order[at]);
    if (targets.empty()) return 0;
    inputs["parents"] = parents;
    inputs["focal_order"] = order;
    final_focal_order = order;
    final_active_input = refresh_active_candidates(
      inputs, parents, dam_candidates, sire_candidates);
    inputs["start_graph_revision"] = graph_revision;
    inputs["message_revision"] = message_revision;
    inputs["recheck_only"] = true;
    inputs["progress_parentage_round"] = round;
    inputs["recheck_focals"] = wrap(targets);
    const auto pass_started = std::chrono::steady_clock::now();
    const List recheck_input = inputs;
    List recheck_result = run_visits(
      compiled_sexp, recheck_input);
    record_typed_ages(recheck_result);
    const double pass_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - pass_started).count();
    const DataFrame recheck_moves = as<DataFrame>(recheck_result["accepted_moves"]);
    const DataFrame recheck_summary = as<DataFrame>(recheck_result["focal_summary"]);
    parents = clone(as<IntegerMatrix>(recheck_result["final_parents"]));
    graph_revision = as<int>(recheck_result["final_graph_revision"]);
    record_complete_evaluations(recheck_summary);
    marking_round = round;
    mark_complete_rechecks_stale(recheck_moves, graph_revision);

    std::vector<int> dirty_nodes;
    if (as<bool>(inputs["adaptive_selective_publisher"])) {
      if (!recheck_result.containsElementNamed("dirty_union_node_ids"))
        engine_stop("adaptive selective publisher did not receive recheck dirty state");
      IntegerVector dirty = recheck_result["dirty_union_node_ids"];
      dirty_nodes.assign(dirty.begin(), dirty.end());
    }
    double publication_seconds = 0.0;
    if (recheck_moves.nrows() > 0) {
      DirtyRegion dirty_region;
      const IntegerVector recheck_moved_children = recheck_moves["child"];
      dirty_region.changed_nodes.assign(
        recheck_moved_children.begin(), recheck_moved_children.end());
      dirty_region.publication_nodes = dirty_nodes;
      dirty_region.graph_revision_before = graph_before_recheck;
      dirty_region.graph_revision_after = graph_revision;
      dirty_region.normalize();
      if (defer_intermediate_publication) {
        ++deferred_publication_phases;
      } else {
        publish("post-complete-recheck", round, true, true, dirty_region);
        publication_seconds = publication.seconds + publication.dispatch_seconds;
      }
    }
    recheck_move_blocks.push_back(recheck_moves);
    recheck_focal_blocks.push_back(recheck_summary);
    if (as<bool>(inputs["configuration_diagnostics"]))
      recheck_configuration_diagnostic_blocks.push_back(
        recheck_result["configuration_diagnostics"]);
    if (as<bool>(inputs["configuration_diagnostics"]))
      recheck_graph_state_diagnostic_blocks.push_back(
        recheck_result["diagnostic_graph_states"]);
    if (as<bool>(inputs["configuration_diagnostics"]))
      recheck_relationship_provenance_blocks.push_back(
        recheck_result["relationship_provenance"]);
    recheck_phase_round.push_back(round);
    recheck_phase_targets.push_back(static_cast<int>(targets.size()));
    recheck_phase_changes.push_back(recheck_moves.nrows());
    recheck_phase_pending_after.push_back(pending_recheck_count());
    const double phase_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
    recheck_phase_seconds.push_back(phase_seconds);
    recheck_phase_setup_seconds.push_back(std::chrono::duration<double>(
      pass_started - started).count());
    recheck_phase_pass_seconds.push_back(pass_seconds);
    recheck_phase_publication_seconds.push_back(publication_seconds);
    inputs["recheck_only"] = false;
    inputs["progress_parentage_round"] = round;
    inputs["recheck_focals"] = IntegerVector();
    last_result = recheck_result;
    return static_cast<int>(recheck_moves.nrows());
}
}
