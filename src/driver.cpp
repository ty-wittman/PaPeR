#include "pap_removed_path.h"
#include "input_checks.h"
#include "typed_age.h"
#include <atomic>
// Typed-parent search driver and R entry points.

#include <Rcpp.h>
#include <array>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <deque>
#include <functional>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <limits>
#include <list>
#include <mutex>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <stdexcept>
#include <condition_variable>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#include "latent_mendelian_support.h"
#include "pedigree_graph_snapshot.h"
#include "probability_log.h"
#include "partition_totals.h"
#include "score_types.h"
#include "predictive_models.h"
#include "panel_rejection.h"
#include "panel_staging.h"
#include "support_rejection.h"
#include "joint_support.h"
#include "initial_pair_check.h"
#include "predictive_error_bounds.h"
#include "result_storage.h"
#include "scorer.h"

using namespace Rcpp;
using pap_internal::DirtyRegion;
using pap_internal::PedigreeGraphSnapshot;
using pap_internal::kNoEvidenceRow;
using pap_internal::union_nodes;

#include "scoring.h"
#include "decisions.h"
#include "message_updates.h"
#include "visits.h"
#include "recheck.h"
#include "graph_traversal.h"
using namespace pap_internal::search;

// Visits share the compiled locus tables and current pedigree inputs.

List visit_candidates_cpp(SEXP compiled_sexp, List inputs) {
  return run_visits(compiled_sexp, inputs);
}

// [[Rcpp::export]]
List publish_messages_cpp(SEXP compiled_sexp,
                                      IntegerMatrix parents,
                                      int graph_revision,
                                      int message_revision,
                                      int locus_parallel_workers = 1) {
  XPtr<EngineState> pointer(compiled_sexp);
  if (pointer.get() == nullptr)
    engine_stop("native Stage-1 publisher compiled snapshot pointer is null");
  EngineState& compiled = *pointer;
  if (parents.nrow() != compiled.n || parents.ncol() != 2 ||
      graph_revision < 0 || message_revision < 0 ||
      locus_parallel_workers != 1)
    engine_stop("native Stage-1 publisher input is invalid");
  compiled.message_revision = message_revision;

  const MessagePublication publication = publish_messages_full(
    compiled, parents, graph_revision);
  List messages(compiled.loci), transmitted(compiled.loci);
  for (int locus_index = 0; locus_index < compiled.loci; ++locus_index) {
    const LocusTables& locus = compiled.locus[
      static_cast<std::size_t>(locus_index)];
    NumericMatrix message(locus.G, compiled.n);
    NumericMatrix transmission(locus.k, compiled.n);
    std::copy(locus.message_ln.begin(), locus.message_ln.end(), message.begin());
    std::copy(locus.transmitted_message_ln.begin(),
              locus.transmitted_message_ln.end(), transmission.begin());
    messages[locus_index] = message;
    transmitted[locus_index] = transmission;
  }
  return List::create(
    _["schema_version"] = "simp-native-stage1-publication-v1",
    _["messages_ln"] = messages,
    _["transmitted_messages_ln"] = transmitted,
    _["node_likelihood_ln"] = publication.node_likelihood_ln,
    _["locus_likelihood_ln"] = locus_likelihood_matrix(
      publication, compiled.loci, compiled.n),
    _["total_likelihood_ln"] = publication.total_likelihood_ln,
    _["message_evidence_signature"] = publication.evidence_signature,
    _["publication_work"] = publication_work_list(publication),
    _["publication_seconds"] = publication.seconds,
    _["graph_revision"] = graph_revision,
    _["message_revision"] = message_revision);
}

// Initial CheckPedigree is a support check on supplied complete pairs, not a
// second chooser. Its four counterfactuals own the same frozen bounded D/C;
// only the focal birth changes, and every missing mate is a fresh founder.
// Reuse the production Q constructor and its ordered marker/status reduction.
static std::vector<ConfigurationScore> check_pair_raw_common(
    SearchState& state,
    std::vector<std::unique_ptr<LocusWorkspace>>& workspaces,
    int child, int dam, int sire, std::string& unavailability_reason) { pap_internal::removed_path("check_pair_raw_common"); }

// With the public canonical transmission/emission primitives, a common
// homozygote witnesses positive support on every legal typed DAG.  The
// publisher also constructs IBD bases in probability space, so finite log p
// alone is insufficient: certify its actual rounded multiplication order.
// An inconclusive certificate retains publications; it is never rejection.
static bool publication_positive_support_certified(
    const EngineState& compiled) { pap_internal::removed_path("publication_positive_support_certified"); }

// One external invocation owns every fixed-50 search round.  The accepted
// pass scorer is reused internally; the static candidate bundle is updated
// only for native parent/revision/focal-order fields between passes.
List search_rounds_cpp(SEXP compiled_sexp, List inputs) {
  XPtr<EngineState> pointer(compiled_sexp);
  if (pointer.get() == nullptr)
    engine_stop("native Stage-1 compiled snapshot pointer is null");
  EngineState& compiled = *pointer;
  if (!inputs.containsElementNamed("parents") ||
      !inputs.containsElementNamed("expected_birth_year") ||
      !inputs.containsElementNamed("id_labels"))
    engine_stop("native Stage-1 input lacks parents, expected years, or IDs");
  const std::string engine = as<std::string>(inputs["search_loop_engine"]);
  if (engine != "bounded_messages")
    engine_stop("native Stage-1 engine accepts bounded_messages only");
  IntegerMatrix parents = clone(as<IntegerMatrix>(inputs["parents"]));
  const NumericVector expected_year = inputs["expected_birth_year"];
  const IntegerVector id_order_rank = inputs["id_order_rank"];
  const CharacterVector ids = inputs["id_labels"];
  const int max_rounds = as<int>(inputs["max_rounds"]);
  const double decrease_tolerance = as<double>(
    inputs["likelihood_decrease_tolerance"]);
  const bool run_check = as<bool>(inputs["run_check_pedigree"]);
  const std::string complete_recheck_mode =
    inputs.containsElementNamed("complete_recheck") ?
      as<std::string>(inputs["complete_recheck"]) : "none";
  const std::string complete_recheck_scope =
    inputs.containsElementNamed("complete_recheck_scope") ?
      as<std::string>(inputs["complete_recheck_scope"]) : "all_complete";
  if (complete_recheck_mode != "none" && complete_recheck_mode != "end_round" &&
      complete_recheck_mode != "convergence" &&
      complete_recheck_mode != "first_dirty_wave_convergence")
    engine_stop("complete_recheck must be none, end_round, convergence, or first_dirty_wave_convergence");
  if (complete_recheck_scope != "all_complete" && complete_recheck_scope != "local")
    engine_stop("complete_recheck_scope must be all_complete or local");
  // Local scope: an accepted move reopens a complete pair only when the focal
  // or its current parents lie within complete_recheck_radius parent-child
  // links of the moved child or its old or new parents, or one of its
  // candidates lies within complete_recheck_candidate_radius (0: the
  // candidate itself was moved or gained or lost a child).
  const int complete_recheck_radius =
    inputs.containsElementNamed("complete_recheck_radius") ?
      as<int>(inputs["complete_recheck_radius"]) : 2;
  const int complete_recheck_candidate_radius =
    inputs.containsElementNamed("complete_recheck_candidate_radius") ?
      as<int>(inputs["complete_recheck_candidate_radius"]) : 0;
  if (complete_recheck_radius < 0 || complete_recheck_candidate_radius < 0)
    engine_stop("complete_recheck radii must be nonnegative");
  // With after_evaluation, a move reopens a complete pair only when it was
  // accepted after that pair's own last evaluation (each move's revision,
  // not the pass's final one); moves the evaluation already saw do not.
  const bool complete_recheck_after_evaluation =
    inputs.containsElementNamed("complete_recheck_after_evaluation") &&
    as<bool>(inputs["complete_recheck_after_evaluation"]);
  if (complete_recheck_after_evaluation && complete_recheck_scope != "local")
    engine_stop("complete_recheck_after_evaluation requires the local recheck scope");
  const bool complete_recheck_enabled = complete_recheck_mode != "none";
  const auto publication_flag = [&](const char* name, bool fallback = false) {
    return inputs.containsElementNamed(name) ? as<bool>(inputs[name]) : fallback;
  };
  const bool defer_publication_requested = publication_flag(
    "defer_intermediate_publication");
  bool defer_intermediate_publication = false;
  std::string publication_deferral_guard = "not_requested";
  if (defer_publication_requested) {
    // The current production registry path constructs a non-null owned common
    // boundary for EVERY required side. A failed boundary stays unavailable;
    // it cannot fall back to frozen message marginals. The empty-registry seam
    // and arbitrary observers preserve their original publication contract.
    publication_deferral_guard = "ownPhi_common_boundary";
    const int mode = inputs.containsElementNamed("relationship_mode") ?
      as<int>(inputs["relationship_mode"]) : 0;
    if ((mode != 0 && mode != 1) ||
        !inputs.containsElementNamed("full_relationship_registry") ||
        as<List>(inputs["full_relationship_registry"]).size() == 0)
      publication_deferral_guard = "legacy_or_missing_relationship_registry";
    for (const char* name : {"progress_callback"})
      if (inputs.containsElementNamed(name) &&
          static_cast<SEXP>(inputs[name]) != R_NilValue)
        publication_deferral_guard = "observer_requires_current_publication";
    for (const char* name : {"configuration_diagnostics",
        "configuration_locus_diagnostics",
        "benchmark_baseline_mode",
        "private_test_relationship_comparison_exception",
        "private_test_raw_comparison_exception",
        "private_test_configuration_post_score_exception",
        "raw_common_background_legacy_missing_mates"})
      if (publication_flag(name))
        publication_deferral_guard = "reference_or_observer_route";
    if (!publication_flag("relationship_common_pair_fast", true) ||
        !publication_flag("relationship_common_scaled_fast", true))
      publication_deferral_guard = "reference_or_observer_route";
    if (inputs.containsElementNamed("raw_common_observation_offsets") &&
        Rf_length(inputs["raw_common_observation_offsets"]) != 0)
      publication_deferral_guard = "reference_or_observer_route";
    if (inputs.containsElementNamed("replay_moves") &&
        as<IntegerMatrix>(inputs["replay_moves"]).nrow() != 0)
      publication_deferral_guard = "reference_or_observer_route";
    if (publication_deferral_guard == "ownPhi_common_boundary" &&
        !publication_positive_support_certified(compiled))
      publication_deferral_guard = "publication_positive_support_not_certified";
    if (publication_flag("message_scoring") ||
        publication_flag("stored_message_scoring"))
      publication_deferral_guard = "message_scoring_requires_publication";
    defer_intermediate_publication =
      publication_deferral_guard == "ownPhi_common_boundary";
  }
  int deferred_publication_phases = 0;
  double deferred_final_publication_seconds = 0.0;
  double deferred_final_publication_wall_seconds = 0.0;
  int graph_revision = as<int>(inputs["start_graph_revision"]);
  int message_revision = as<int>(inputs["message_revision"]);
  if (expected_year.size() != compiled.n || ids.size() != compiled.n || max_rounds < 1)
    engine_stop("native Stage-1 static input dimensions are invalid");
  const CandidateSide dam_candidates(inputs, "dam", compiled.n);
  const CandidateSide sire_candidates(inputs, "sire", compiled.n);
  bool supplied_links_present = false;
  for (int child = 0; child < compiled.n; ++child)
    supplied_links_present = supplied_links_present ||
      parents(child, 0) > 0 || parents(child, 1) > 0;

  std::vector<std::string> publication_reason;
  std::vector<int> publication_round, publication_graph, publication_message;
  std::vector<double> publication_total, publication_seconds,
    publication_dispatch_seconds, publication_dirty_nodes,
    publication_total_nodes, publication_nodes_visited,
    publication_nodes_recomputed, publication_message_vectors,
    publication_genotype_elements, publication_transmitted_elements,
    publication_likelihood_cells, publication_aggregate_updates,
    publication_signature_global_scans, publication_copies,
    publication_allocations, publication_finalization_seconds,
    publication_worker_seconds_max, publication_worker_seconds_mean,
    publication_worker_imbalance;
  std::vector<int> publication_requested_workers,
    publication_effective_workers, publication_worker_thread_count,
    publication_job_count, publication_worker_r_api_calls;
  std::vector<std::string> publication_parent, publication_evidence;
  std::vector<std::string> publication_mode;
  MessagePublication publication;
  int last_message_graph_revision = graph_revision;
  auto publish = [&](const char* reason, int round, bool advance,
                     bool selective_allowed,
                     const DirtyRegion& dirty_region) {
    if (advance) ++message_revision;
    // Original publication owns the lease. Conservatively invalidate retained
    // directed cuts at republication; graph moves use their broader family
    // dependency envelope rather than the forward-only descendant dirty set.
    compiled.message_revision = message_revision;
    const auto dispatch_started = std::chrono::steady_clock::now();
    const std::vector<int>& dirty_nodes = dirty_region.publication_nodes;
    const double dirty_fraction = compiled.n > 0 ?
      static_cast<double>(dirty_nodes.size()) / compiled.n : 1.0;
    const bool use_selective = selective_allowed &&
      as<bool>(inputs["adaptive_selective_publisher"]) &&
      dirty_nodes.size() < static_cast<std::size_t>(compiled.n) &&
      dirty_fraction < as<double>(inputs["adaptive_selective_threshold"]);
    if (use_selective) {
      publish_messages_selective(
        compiled, parents, graph_revision, publication, dirty_region);
      publication.requested_workers = 1;
      publication.effective_publication_workers = 1;
      publication.worker_thread_count = 0;
      publication.publication_job_count = 0;
      publication.worker_r_api_calls = 0;
      publication.worker_seconds.clear();
      publication.locus_owner.clear();
      publication.finalization_seconds = 0.0;
    } else {
      publication = publish_messages_full(
        compiled, parents, graph_revision);
    }
    publication.dispatch_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - dispatch_started).count() -
      publication.seconds;
    if (publication.dispatch_seconds < 0.0) publication.dispatch_seconds = 0.0;
    publication_reason.push_back(reason);
    publication_round.push_back(round);
    publication_graph.push_back(graph_revision);
    publication_message.push_back(message_revision);
    publication_total.push_back(publication.total_likelihood_ln);
    publication_seconds.push_back(publication.seconds);
    publication_dispatch_seconds.push_back(publication.dispatch_seconds);
    const bool dirty_state_available = selective_allowed &&
      as<bool>(inputs["adaptive_selective_publisher"]);
    publication_dirty_nodes.push_back(static_cast<double>(
      dirty_state_available ? dirty_nodes.size() : compiled.n));
    publication_total_nodes.push_back(static_cast<double>(compiled.n));
    publication_nodes_visited.push_back(static_cast<double>(
      publication.nodes_visited));
    publication_nodes_recomputed.push_back(static_cast<double>(
      publication.nodes_recomputed));
    publication_message_vectors.push_back(static_cast<double>(
      publication.message_vectors_rebuilt));
    publication_genotype_elements.push_back(static_cast<double>(
      publication.genotype_state_elements));
    publication_transmitted_elements.push_back(static_cast<double>(
      publication.transmitted_state_elements));
    publication_likelihood_cells.push_back(static_cast<double>(
      publication.likelihood_cells));
    publication_aggregate_updates.push_back(static_cast<double>(
      publication.aggregate_updates));
    publication_signature_global_scans.push_back(static_cast<double>(
      publication.signature_global_scans));
    publication_copies.push_back(static_cast<double>(publication.copies));
    publication_allocations.push_back(static_cast<double>(
      publication.allocations));
    publication_requested_workers.push_back(publication.requested_workers);
    publication_effective_workers.push_back(
      publication.effective_publication_workers);
    publication_worker_thread_count.push_back(
      publication.worker_thread_count);
    publication_job_count.push_back(publication.publication_job_count);
    publication_worker_r_api_calls.push_back(publication.worker_r_api_calls);
    publication_finalization_seconds.push_back(
      publication.finalization_seconds);
    const double worker_total = std::accumulate(
      publication.worker_seconds.begin(), publication.worker_seconds.end(), 0.0);
    const double worker_mean = publication.worker_seconds.empty() ? 0.0 :
      worker_total / static_cast<double>(publication.worker_seconds.size());
    const double worker_max = publication.worker_seconds.empty() ? 0.0 :
      *std::max_element(publication.worker_seconds.begin(),
                        publication.worker_seconds.end());
    publication_worker_seconds_max.push_back(worker_max);
    publication_worker_seconds_mean.push_back(worker_mean);
    publication_worker_imbalance.push_back(
      worker_mean > 0.0 ? worker_max / worker_mean : 0.0);
    publication_mode.push_back(use_selective ? "selective" : "full");
    publication_parent.push_back(parent_signature(ids, parents));
    publication_evidence.push_back(publication.evidence_signature);
    last_message_graph_revision = graph_revision;
  };

  DirtyRegion initial_region;
  initial_region.graph_revision_before = graph_revision;
  initial_region.graph_revision_after = graph_revision;
  publish("initial", 0, false, false, initial_region);
  std::vector<int> check_child, check_old_dam, check_old_sire, check_pass;
  std::vector<int> check_available, check_proven_rejection, check_hard_rejection;
  std::vector<int> check_graph_revision;
  std::vector<double> check_support, check_raw_null, check_raw_dam,
    check_raw_sire, check_raw_pair, check_rejection_upper;
  std::vector<std::string> check_action, check_unavailability_reason;
  const bool joint_initial_check=pap_internal::native_demand::live&&
    inputs.containsElementNamed("pedigree_likelihood")&&as<bool>(inputs["pedigree_likelihood"])&&
    (!inputs.containsElementNamed("likelihood_check_pedigree_joint")||as<bool>(inputs["likelihood_check_pedigree_joint"]));
  std::vector<double> check_identity_null,check_identity_dam,check_identity_sire,check_identity_pair;
  if (run_check && supplied_links_present) {
    // One owner for the entire check wave. Graph/role metadata is updated by
    // the same authoritative mutation used by search, not copied per pair.
    SearchInputs check_input(inputs, compiled.n);
    check_input.start_graph_revision = graph_revision;
    check_input.message_revision = message_revision;
    check_input.lazy_independent_gates = false;
    check_input.raw_common_background_legacy_missing_mates = false;
    // Integrate this scientific correction with retained arithmetic first.
    // Optional accelerated Q requires its existing complete-panel refinement
    // contract before it may make CheckPedigree's >= threshold decision.
    check_input.raw_scaled_fast = false;
    check_input.validate();
    SearchState check_state(compiled, check_input, parents);
    apply_raw_provider_settings(check_input, check_state.pass);
    std::vector<std::unique_ptr<LocusWorkspace>> check_workspaces;
    check_workspaces.reserve(static_cast<std::size_t>(compiled.loci));
    for (int locus = 0; locus < compiled.loci; ++locus)
      check_workspaces.emplace_back(new LocusWorkspace(check_state.pass, locus));
    for (int child = 1; child <= compiled.n; ++child) {
      const int dam = parents(child - 1, 0);
      const int sire = parents(child - 1, 1);
      if (dam == 0 || sire == 0) continue;
      // Preserve the original check's explicit structural rejection scope.
      // New role, age, OH or ME policy is outside this scientific work unit.
      const bool hard_rejected = dam < 0 || dam > compiled.n ||
        sire < 0 || sire > compiled.n || dam == sire ||
        dam == child || sire == child;
      std::string unavailability_reason;
      namespace check = pap_internal::initial_pair_check;
      std::array<check::Panel,4> Q;
      check::Result checked;
      std::array<double,4> identity{{NA_REAL,NA_REAL,NA_REAL,NA_REAL}};
      if(joint_initial_check){
        try{
          check_state.prepare_visit_ages(child);
          pap_internal::typed_age::Scope age_scope(&check_state.visit_age_bounds);
          auto* live=pap_internal::native_demand::live;
          // The initial check precedes ordinary visits. Establish exactly the
          // same live owner and focal cavity before asking any numerical question.
          live->begin_visit(child,check_state.pass.dam,check_state.pass.sire,
            check_state.graph_revision,false,inputs);
          const auto gate=pap_internal::joint_support_gate::evaluate(
            dam,sire,false,false,check_input.Tassign_ln,[&](int d,int s){
              const int at=(d?1:0)+(s?2:0);
              try{identity[at]=live->joint_support(d,s);}
              catch(const std::exception& error){if(unavailability_reason.empty())unavailability_reason=error.what();}
              return identity[at];
            },check_input.support_available_baselines);
          checked.support_ln=gate.support[2];checked.support_available=gate.available;
          checked.support_pass=gate.pass&&!hard_rejected;
          checked.proven_rejection=hard_rejected||gate.rejected;
          for(double value:gate.support)if(!std::isnan(value)&&
            (std::isnan(checked.rejection_upper_ln)||value<checked.rejection_upper_ln))checked.rejection_upper_ln=value;
        }catch(const std::exception& error){
          unavailability_reason=error.what();checked.proven_rejection=hard_rejected;
        }
      }else{
        const auto scores=check_pair_raw_common(
          check_state,check_workspaces,child,dam,sire,unavailability_reason);
        for(std::size_t i=0;i<scores.size();++i){
          Q[i].ln=scores[i].objective_ln;
          Q[i].available=scores[i].status==0&&scores[i].raw_comparison_available;
          Q[i].complete=scores[i].raw_loci_accumulated==static_cast<std::size_t>(compiled.loci);
        }
        checked=check::evaluate(Q,check_input.Tassign_ln,hard_rejected);
      }
      check_unavailability_reason.push_back(std::move(unavailability_reason));
      const double support = checked.support_ln;
      const bool support_available = checked.support_available;
      const bool support_pass = checked.support_pass;
      const bool proven_rejection = checked.proven_rejection;
      check_child.push_back(child);
      check_old_dam.push_back(dam);
      check_old_sire.push_back(sire);
      check_graph_revision.push_back(check_state.graph_revision);
      check_support.push_back(support);
      check_rejection_upper.push_back(checked.rejection_upper_ln);
      check_raw_null.push_back(joint_initial_check?NA_REAL:Q[0].ln);
      check_raw_dam.push_back(joint_initial_check?NA_REAL:Q[1].ln);
      check_raw_sire.push_back(joint_initial_check?NA_REAL:Q[2].ln);
      check_raw_pair.push_back(joint_initial_check?NA_REAL:Q[3].ln);
      check_identity_null.push_back(identity[0]);check_identity_dam.push_back(identity[1]);
      check_identity_sire.push_back(identity[2]);check_identity_pair.push_back(identity[3]);
      check_available.push_back(support_available ? 1 : 0);
      check_proven_rejection.push_back(proven_rejection ? 1 : 0);
      check_hard_rejection.push_back(hard_rejected ? 1 : 0);
      check_pass.push_back(proven_rejection ? 0 :
        (support_available ? (support_pass ? 1 : 0) : NA_INTEGER));
      check_action.push_back(proven_rejection ? "drop_complete_pair" :
        (support_available ? "retained" : "retained_unavailable"));
      if (proven_rejection) {
        check_state.apply_move(child, 0, 0);
        parents(child - 1, 0) = 0;
        parents(child - 1, 1) = 0;
        graph_revision = check_state.graph_revision;
      }
    }
    DirtyRegion check_region;
    check_region.graph_revision_before = initial_region.graph_revision_after;
    check_region.graph_revision_after = graph_revision;
    publish("post-CheckPedigree", 0, true, false, check_region);
  }

  const double initial_likelihood = publication.total_likelihood_ln;
  std::vector<double> likelihood_history(1, initial_likelihood);
  std::vector<int> likelihood_history_current(1, 1);
  std::vector<List> move_blocks, focal_blocks, configuration_diagnostic_blocks,
    configuration_locus_diagnostic_blocks, graph_state_diagnostic_blocks,
    relationship_provenance_blocks;
  std::vector<int> round_number, round_changes, round_recheck_changes,
    round_pending, round_search_residual, round_links,
    round_graph_start, round_graph_end, round_message_start,
    round_message_end, round_publication, round_source_converged,
    round_converged, round_decrease, round_likelihood_current,
    round_likelihood_graph, round_publication_deferred;
  std::vector<double> round_total, round_delta, round_pass_seconds,
    round_publication_seconds, round_wall_seconds, round_cpu_seconds;
  std::vector<std::string> round_publication_reason;
  std::string status = "production_round_limit_reached";
  bool converged = false;
  List last_result;
  List final_active_input;
  IntegerVector final_focal_order;
  // Complete-pair rechecks are native passes over the same materialized
  // fixed-50 arrays. A focal is pending only when an accepted move has made
  // its previously evaluated complete state stale.
  std::vector<int> last_complete_evaluation_revision(
    static_cast<std::size_t>(compiled.n + 1), NA_INTEGER);
  std::vector<char> recheck_pending(static_cast<std::size_t>(compiled.n + 1), 0);
  std::vector<List> recheck_move_blocks, recheck_focal_blocks,
    recheck_configuration_diagnostic_blocks, recheck_graph_state_diagnostic_blocks,
    recheck_relationship_provenance_blocks;
  std::vector<int> recheck_phase_round, recheck_phase_targets,
    recheck_phase_changes, recheck_phase_pending_after;
  std::vector<double> recheck_phase_seconds, recheck_phase_setup_seconds,
    recheck_phase_pass_seconds, recheck_phase_publication_seconds;
  bool first_dirty_wave_done = false;
  if (complete_recheck_enabled)
    for (int focal = 1; focal <= compiled.n; ++focal)
      if (parents(focal - 1, 0) > 0 && parents(focal - 1, 1) > 0)
        last_complete_evaluation_revision[static_cast<std::size_t>(focal)] =
          graph_revision;
  // Diagnostic (opt-in): complete pairs listed in
  // complete_recheck_initial_focals start pending, so the first recheck wave
  // reopens them from the supplied starting state.
  if (complete_recheck_enabled &&
      inputs.containsElementNamed("complete_recheck_initial_focals")) {
    const IntegerVector initial = inputs["complete_recheck_initial_focals"];
    for (int focal : initial) {
      if (focal < 1 || focal > compiled.n)
        engine_stop("complete_recheck_initial_focals is outside the native ID range");
      recheck_pending[static_cast<std::size_t>(focal)] = 1;
    }
  }
  // Opt-in: a complete pair whose last evaluation had every margin (identity,
  // pair, dam and sire support) above complete_recheck_confident_ln is not
  // reopened by later moves.
  const double complete_recheck_confident_ln =
    inputs.containsElementNamed("complete_recheck_confident_ln") ?
      as<double>(inputs["complete_recheck_confident_ln"]) : R_PosInf;
  // complete_recheck_confident_relationship (default TRUE) also requires the
  // full relationship contrast of the pair above the threshold.
  const bool complete_recheck_confident_relationship =
    !inputs.containsElementNamed("complete_recheck_confident_relationship") ||
    as<bool>(inputs["complete_recheck_confident_relationship"]);
  // The relationship contrast (PO against its best alternative) is on a
  // smaller scale than the supports; it may take its own threshold.
  const double complete_recheck_confident_relationship_ln =
    inputs.containsElementNamed("complete_recheck_confident_relationship_ln") ?
      as<double>(inputs["complete_recheck_confident_relationship_ln"]) : complete_recheck_confident_ln;
  std::vector<double> last_complete_margin(static_cast<std::size_t>(compiled.n + 1), R_NegInf);
  std::vector<std::array<double, 5>> last_complete_margins(static_cast<std::size_t>(compiled.n + 1),
    std::array<double, 5>{{NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL}});
  // Every reopening withheld by the confidence rule, for audit.
  std::vector<int> skip_round, skip_focal, skip_evaluated, skip_move;
  std::vector<double> skip_identity, skip_pair, skip_dam, skip_sire, skip_relationship;
  int marking_round = 0;
  // complete_recheck_confident_no_child_change (opt-in): the exemption also
  // requires that no accepted move since the pair's evaluation had the focal
  // as an old or new parent (it gained or lost a child, or a child's co-parent
  // changed).
  const bool complete_recheck_confident_no_child_change =
    inputs.containsElementNamed("complete_recheck_confident_no_child_change") &&
    as<bool>(inputs["complete_recheck_confident_no_child_change"]);
  const bool ancestral_rechecks=inputs.containsElementNamed("likelihood_evidence")?
    as<std::string>(inputs["likelihood_evidence"])=="ancestral":
    (inputs.containsElementNamed("pedigree_likelihood")&&as<bool>(inputs["pedigree_likelihood"]));
  const bool track_child_change=complete_recheck_confident_no_child_change&&!ancestral_rechecks;
  std::vector<int> last_child_change_revision(track_child_change?static_cast<std::size_t>(compiled.n + 1):0, -1);
  int child_change_vetoes = 0;
  auto record_complete_evaluations = [&](const DataFrame& summary) {
    if (!complete_recheck_enabled) return;
    const IntegerVector evaluated = summary["focal"];
    const IntegerVector evaluated_graph = summary["graph_revision_after"];
    const IntegerVector skipped = summary["complete_pair_skipped"];
    const NumericVector identity = summary["identity_margin_ln"], pair = summary["pair_support_ln"],
      dam_support = summary["dam_support_ln"], sire_support = summary["sire_support_ln"],
      relationship = summary["full_relationship_ln"];
    for (int index = 0; index < evaluated.size(); ++index) {
      const int focal = evaluated[index];
      // Skipping an incumbent pair is not a scientific evaluation. A newly
      // installed pair was evaluated at its own accepted revision, but later
      // mutations in this pass can make that question stale again.
      if (skipped[index]) continue;
      recheck_pending[static_cast<std::size_t>(focal)] = 0;
      if (parents(focal - 1, 0) > 0 && parents(focal - 1, 1) > 0) {
        last_complete_evaluation_revision[static_cast<std::size_t>(focal)] =
          evaluated_graph[index];
        double margin = R_PosInf;
        for (double x : {identity[index], pair[index], dam_support[index], sire_support[index]})
          margin = ISNAN(x) ? R_NegInf : std::min(margin, x);
        last_complete_margin[static_cast<std::size_t>(focal)] = margin;
        last_complete_margins[static_cast<std::size_t>(focal)] = {{identity[index], pair[index],
          dam_support[index], sire_support[index], relationship[index]}};
      }
    }
  };
  auto confident_complete = [&](int focal, int move_revision) {
    if (!(last_complete_margin[static_cast<std::size_t>(focal)] > complete_recheck_confident_ln)) return false;
    const auto& m = last_complete_margins[static_cast<std::size_t>(focal)];
    if (complete_recheck_confident_relationship && !(m[4] > complete_recheck_confident_relationship_ln)) return false;
    if (track_child_change &&
        last_child_change_revision[static_cast<std::size_t>(focal)] >
          last_complete_evaluation_revision[static_cast<std::size_t>(focal)]) { ++child_change_vetoes; return false; }
    skip_round.push_back(marking_round); skip_focal.push_back(focal);
    skip_evaluated.push_back(last_complete_evaluation_revision[static_cast<std::size_t>(focal)]);
    skip_move.push_back(move_revision);
    skip_identity.push_back(m[0]); skip_pair.push_back(m[1]); skip_dam.push_back(m[2]);
    skip_sire.push_back(m[3]); skip_relationship.push_back(m[4]);
    return true;
  };
  List typed_age_blocks;
  std::vector<int> latest_typed_age_revision(compiled.n+1,-1);
  auto record_typed_ages=[&](const List& result){
    if(result.containsElementNamed("typed_age_propagation"))typed_age_blocks.push_back(result["typed_age_propagation"]);
    if(result.containsElementNamed("typed_age_changed_revisions")){
      const IntegerVector changed=result["typed_age_changed_revisions"];
      for(int v=1;v<=compiled.n;++v)latest_typed_age_revision[v]=std::max(latest_typed_age_revision[v],int(changed[v]));}
  };
  RecheckGraph recheck_graph;
  auto mark_complete_rechecks_stale = [&](const DataFrame& moves,
                                           int move_graph_revision) {
    if (!complete_recheck_enabled || moves.nrows() == 0) return;
    recheck_graph.synchronize(parents, move_graph_revision);
    // Date feasibility is separate from genotype evidence. Descendant moves
    // can change an ancestor's typed constraint component, hence applicability of a
    // question whose genotype boundary itself did not change. Exact dates
    // never enter this additional reopening path.
    if(std::any_of(latest_typed_age_revision.begin(),latest_typed_age_revision.end(),[](int x){return x>=0;})){
      auto latest=latest_typed_age_revision;std::deque<int> todo;
      for(int v=1;v<=compiled.n;++v)if(latest[v]>=0)todo.push_back(v);
      recheck_graph.spread_descendant_revisions(latest, std::move(todo));
      for(int focal=1;focal<=compiled.n;++focal)if(parents(focal-1,0)&&parents(focal-1,1)){
        int newest=latest[focal];
        for(const CandidateSide* side:{&dam_candidates,&sire_candidates})
          for(int at=side->offsets[focal-1];at<side->offsets[focal];++at)
            if(side->permanent[at]||side->ids[at]==parents(focal-1,side==&dam_candidates?0:1))newest=std::max(newest,latest[side->ids[at]]);
        if(newest>last_complete_evaluation_revision[focal]&&!confident_complete(focal,newest))recheck_pending[focal]=1;
      }
    }
    if(ancestral_rechecks) {
      // A moved birth changes that child's ancestry and its descendants only.
      // Reopen a complete question only if an endpoint/candidate ancestry is newer.
      std::vector<int> latest(compiled.n+1,-1);std::deque<int> queue;
      const IntegerVector child=moves["child"],revision=moves["graph_revision_after"];
      for(int j=0;j<child.size();++j){latest[child[j]]=std::max(latest[child[j]],int(revision[j]));queue.push_back(child[j]);}
      recheck_graph.spread_descendant_revisions(latest, std::move(queue));
      for(int focal=1;focal<=compiled.n;++focal){
        if(!parents(focal-1,0)||!parents(focal-1,1))continue;
        int newest=std::max({latest[focal],latest[parents(focal-1,0)],latest[parents(focal-1,1)]});
        for(const CandidateSide* side:{&dam_candidates,&sire_candidates})
          for(int at=side->offsets[focal-1];at<side->offsets[focal];++at)
            if(side->permanent[at]||side->ids[at]==parents(focal-1,side==&dam_candidates?0:1))
              newest=std::max(newest,latest[side->ids[at]]);
        if(newest>last_complete_evaluation_revision[focal]&&!confident_complete(focal,newest))recheck_pending[focal]=1;
      }
      return;
    }
    if(track_child_change){
      const IntegerVector revision = moves["graph_revision_after"];
      for (const char* column : {"old_dam", "old_sire", "new_dam", "new_sire"}) {
        const IntegerVector ids = moves[column];
        for (int i = 0; i < ids.size(); ++i)
          if (ids[i] > 0 && ids[i] <= compiled.n)
            last_child_change_revision[ids[i]] = std::max(last_child_change_revision[ids[i]], static_cast<int>(revision[i]));
      }
    }
    std::vector<int> near;   // distance to a moved node, or -1
    if (complete_recheck_scope == "local") {
      const int n = compiled.n;
      std::vector<int> distance(static_cast<std::size_t>(n + 1), -1);
      std::deque<int> queue;
      for (const char* column : {"child", "old_dam", "old_sire", "new_dam", "new_sire"}) {
        const IntegerVector ids = moves[column];
        for (int v : ids)
          if (v > 0 && v <= n && distance[v] < 0) { distance[v] = 0; queue.push_back(v); }
      }
      const int reach = std::max(complete_recheck_radius, complete_recheck_candidate_radius);
      recheck_graph.spread_distances(distance, std::move(queue), reach);
      near = std::move(distance);
    }
    if (complete_recheck_after_evaluation) {
      // Latest move revision within each radius of every individual.
      const int n = compiled.n;
      std::vector<int> latest(static_cast<std::size_t>(n + 1), -1);
      const IntegerVector revision = moves["graph_revision_after"];
      for (const char* column : {"child", "old_dam", "old_sire", "new_dam", "new_sire"}) {
        const IntegerVector ids = moves[column];
        for (int i = 0; i < ids.size(); ++i)
          if (ids[i] > 0 && ids[i] <= n) latest[ids[i]] = std::max(latest[ids[i]], static_cast<int>(revision[i]));
      }
      const std::vector<int> own = recheck_graph.spread_neighborhood_max(latest, complete_recheck_radius);
      const std::vector<int> candidate = recheck_graph.spread_neighborhood_max(latest, complete_recheck_candidate_radius);
      for (int focal = 1; focal <= n; ++focal) {
        if (parents(focal - 1, 0) == 0 || parents(focal - 1, 1) == 0) continue;
        int newest = std::max({own[focal], own[parents(focal - 1, 0)], own[parents(focal - 1, 1)]});
        for (const CandidateSide* side : {&dam_candidates, &sire_candidates})
          for (int at = side->offsets[focal - 1]; at < side->offsets[focal]; ++at) {
            const int id = side->ids[static_cast<std::size_t>(at)];
            if (id > 0 && id <= n) newest = std::max(newest, candidate[id]);
          }
        if (newest > last_complete_evaluation_revision[static_cast<std::size_t>(focal)] && !confident_complete(focal, newest))
          recheck_pending[static_cast<std::size_t>(focal)] = 1;
      }
      return;
    }
    auto within = [&](int v, int radius) { return v > 0 && near[v] >= 0 && near[v] <= radius; };
    auto touched = [&](int focal) {
      if (near.empty()) return true;
      if (within(focal, complete_recheck_radius) || within(parents(focal - 1, 0), complete_recheck_radius) ||
          within(parents(focal - 1, 1), complete_recheck_radius)) return true;
      for (const CandidateSide* side : {&dam_candidates, &sire_candidates})
        for (int at = side->offsets[focal - 1]; at < side->offsets[focal]; ++at) {
          const int candidate = side->ids[static_cast<std::size_t>(at)];
          if (candidate <= compiled.n && within(candidate, complete_recheck_candidate_radius)) return true;
        }
      return false;
    };
    for (int focal = 1; focal <= compiled.n; ++focal) {
      if (parents(focal - 1, 0) == 0 || parents(focal - 1, 1) == 0) continue;
      if (!touched(focal)) continue;
      // The contracted all-complete envelope includes moved focals whose
      // actual evaluation preceded a later accepted family/graph change.
      if (last_complete_evaluation_revision[static_cast<std::size_t>(focal)] <
            move_graph_revision && !confident_complete(focal, move_graph_revision))
        recheck_pending[static_cast<std::size_t>(focal)] = 1;
    }
  };
  auto pending_recheck_count = [&]() {
    int count = 0;
    for (int focal = 1; focal <= compiled.n; ++focal)
      if (recheck_pending[static_cast<std::size_t>(focal)] &&
          parents(focal - 1, 0) > 0 && parents(focal - 1, 1) > 0)
        ++count;
    return count;
  };
  // Search convergence is a logical fixed point. Publisher likelihood,
  // revisions, and timing are diagnostic telemetry only. Pending complete
  // obligations distinguish two otherwise identical pedigrees.
  auto logical_state_signature = [&]() {
    std::ostringstream out;
    out << parent_signature(ids, parents) << "|pending=";
    for (int focal = 1; focal <= compiled.n; ++focal)
      out << (recheck_pending[static_cast<std::size_t>(focal)] ? '1' : '0');
    return out.str();
  };
  std::set<std::string> seen_logical_states;
  seen_logical_states.insert(logical_state_signature());
  int logical_cycle_events = 0;
  RecheckContext recheck_context{
    complete_recheck_enabled,
    pending_recheck_count,
    graph_revision,
    parents,
    expected_year,
    id_order_rank,
    inputs,
    final_focal_order,
    final_active_input,
    dam_candidates,
    sire_candidates,
    message_revision,
    compiled_sexp,
    record_typed_ages,
    record_complete_evaluations,
    marking_round,
    mark_complete_rechecks_stale,
    defer_intermediate_publication,
    deferred_publication_phases,
    publish,
    publication,
    recheck_pending,
    last_result,
    recheck_move_blocks,
    recheck_focal_blocks,
    recheck_configuration_diagnostic_blocks,
    recheck_graph_state_diagnostic_blocks,
    recheck_relationship_provenance_blocks,
    recheck_phase_round,
    recheck_phase_targets,
    recheck_phase_changes,
    recheck_phase_pending_after,
    recheck_phase_seconds,
    recheck_phase_setup_seconds,
    recheck_phase_pass_seconds,
    recheck_phase_publication_seconds
  };

  for (int round = 1; round <= max_rounds; ++round) {
    const auto round_started = std::chrono::steady_clock::now();
    const auto round_cpu_started = std::clock();
    const int before_graph = graph_revision;
    const int before_message = message_revision;
    const bool before_likelihood_current = graph_revision == last_message_graph_revision;
    const double before_likelihood = publication.total_likelihood_ln;
    inputs["parents"] = parents;
    inputs["recheck_only"] = false;
    inputs["progress_parentage_round"] = round;
    inputs["recheck_focals"] = IntegerVector();
    final_focal_order = focal_order(
      parents, expected_year, id_order_rank);
    inputs["focal_order"] = final_focal_order;
    final_active_input = refresh_active_candidates(
      inputs, parents, dam_candidates, sire_candidates);
    inputs["start_graph_revision"] = graph_revision;
    inputs["message_revision"] = message_revision;
    const auto pass_started = std::chrono::steady_clock::now();
    const List parentage_input = inputs;
    last_result = run_visits(
      compiled_sexp, parentage_input);
    record_typed_ages(last_result);
    const double pass_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - pass_started).count();
    parents = clone(as<IntegerMatrix>(last_result["final_parents"]));
    graph_revision = as<int>(last_result["final_graph_revision"]);
    DataFrame moves = as<DataFrame>(last_result["accepted_moves"]);
    int changes = moves.nrows();
    record_complete_evaluations(as<DataFrame>(last_result["focal_summary"]));
    marking_round = round;
    mark_complete_rechecks_stale(moves, graph_revision);

    std::vector<int> dirty_nodes;
    if (as<bool>(inputs["adaptive_selective_publisher"])) {
      if (!last_result.containsElementNamed("dirty_union_node_ids"))
        engine_stop("adaptive selective publisher did not receive pass dirty state");
      IntegerVector dirty = last_result["dirty_union_node_ids"];
      dirty_nodes.assign(dirty.begin(), dirty.end());
    }
    move_blocks.push_back(moves);
    focal_blocks.push_back(last_result["focal_summary"]);
    if (as<bool>(inputs["configuration_diagnostics"]))
      configuration_diagnostic_blocks.push_back(
        last_result["configuration_diagnostics"]);
    if (as<bool>(inputs["configuration_diagnostics"]))
      graph_state_diagnostic_blocks.push_back(
        last_result["diagnostic_graph_states"]);
    if (as<bool>(inputs["configuration_diagnostics"]))
      relationship_provenance_blocks.push_back(
        last_result["relationship_provenance"]);
    if (as<bool>(inputs["configuration_locus_diagnostics"]))
      configuration_locus_diagnostic_blocks.push_back(
        last_result["configuration_locus_diagnostics"]);
    double publish_seconds = 0.0;
    if (changes == 0) {
      // Reuse is recorded in the round ledger, not as a publication event.
      // This matches the accepted R lifecycle's message-revision history.
    } else {
      DirtyRegion dirty_region;
      const IntegerVector moved_children = moves["child"];
      dirty_region.changed_nodes.assign(
        moved_children.begin(), moved_children.end());
      dirty_region.publication_nodes = dirty_nodes;
      dirty_region.graph_revision_before = before_graph;
      dirty_region.graph_revision_after = graph_revision;
      dirty_region.normalize();
      if (defer_intermediate_publication) {
        ++deferred_publication_phases;
      } else {
        publish("post-parentage-pass", round, true, true, dirty_region);
        publish_seconds = publication.seconds + publication.dispatch_seconds;
      }
    }
    const bool likelihood_current = graph_revision == last_message_graph_revision;
    const double current_likelihood = likelihood_current ?
      publication.total_likelihood_ln : NA_REAL;
    const double delta = likelihood_current && before_likelihood_current ?
      publication.total_likelihood_ln - before_likelihood : NA_REAL;
    const bool decrease = delta < -decrease_tolerance;
    // A bounded-cavity publication can move in either direction. Accepted
    // parentage moves remain valid; only the assignment/search state controls
    // convergence.
    const bool source_converged = changes == 0;
    int links = 0;
    for (int i = 0; i < parents.nrow(); ++i)
      links += (parents(i, 0) > 0 ? 1 : 0) + (parents(i, 1) > 0 ? 1 : 0);
    round_number.push_back(round);
    round_changes.push_back(changes);
    round_recheck_changes.push_back(0);
    round_pending.push_back(0);
    round_search_residual.push_back(changes);
    round_links.push_back(links);
    round_graph_start.push_back(before_graph);
    round_graph_end.push_back(graph_revision);
    round_message_start.push_back(before_message);
    round_message_end.push_back(message_revision);
    round_publication.push_back(changes > 0 && !defer_intermediate_publication ? 1 : 0);
    round_publication_deferred.push_back(changes > 0 && defer_intermediate_publication ? 1 : 0);
    round_source_converged.push_back(source_converged ? 1 : 0);
    bool round_is_converged = source_converged;
    round_decrease.push_back(likelihood_current && before_likelihood_current ?
      (decrease ? 1 : 0) : NA_INTEGER);
    round_likelihood_current.push_back(likelihood_current ? 1 : 0);
    round_likelihood_graph.push_back(last_message_graph_revision);
    round_total.push_back(current_likelihood);
    round_delta.push_back(delta);
    round_pass_seconds.push_back(pass_seconds);
    round_publication_seconds.push_back(publish_seconds);
    round_publication_reason.push_back(changes > 0 ?
      (defer_intermediate_publication ? "ownPhi-publication-deferred" : "post-parentage-pass") :
      (likelihood_current ? "unchanged-state-reuse" : "unpublished-state-reuse"));
    round_wall_seconds.push_back(std::chrono::duration<double>(
      std::chrono::steady_clock::now() - round_started).count());
    const auto round_cpu_finished = std::clock();
    round_cpu_seconds.push_back(round_cpu_started == std::clock_t(-1) ||
      round_cpu_finished == std::clock_t(-1) ? NA_REAL :
      double(round_cpu_finished - round_cpu_started) / CLOCKS_PER_SEC);
    likelihood_history.push_back(current_likelihood);
    likelihood_history_current.push_back(likelihood_current ? 1 : 0);
    if (!complete_recheck_enabled) {
      const int residual = changes;
      round_search_residual.back() = residual;
      if (residual == 0) {
        status = "converged";
        converged = true;
        round_is_converged = true;
        round_converged.push_back(1);
        break;
      }
      if (!seen_logical_states.insert(logical_state_signature()).second) {
        ++logical_cycle_events;
        status = "cycle_detected";
        round_converged.push_back(0);
        break;
      }
      round_converged.push_back(0);
      continue;
    }
    int accepted_rechecks = 0;
    const bool h1_intermediate_wave =
      complete_recheck_mode == "first_dirty_wave_convergence" &&
      !first_dirty_wave_done && !source_converged &&
      pending_recheck_count() > 0;
    const bool run_convergence_flush = source_converged;
    if (complete_recheck_mode == "end_round" || h1_intermediate_wave ||
        run_convergence_flush) {
      accepted_rechecks = run_recheck_wave(recheck_context, round);
      if (h1_intermediate_wave) first_dirty_wave_done = true;
    }
    const int pending_after = pending_recheck_count();
    const int residual = changes + accepted_rechecks + pending_after;
    round_recheck_changes.back() = accepted_rechecks;
    round_pending.back() = pending_after;
    round_search_residual.back() = residual;
    round_is_converged = residual == 0;
    round_converged.push_back(round_is_converged ? 1 : 0);
    if (round_is_converged) {
      status = "converged";
      converged = true;
      break;
    }
    if (!seen_logical_states.insert(logical_state_signature()).second) {
      ++logical_cycle_events;
      status = "cycle_detected";
      break;
    }
  }

  // Final output/hydration/FinalLLR always receive a physical publication for
  // the final typed graph. Revisions advance only here, never at deferred waves.
  if (defer_intermediate_publication && graph_revision != last_message_graph_revision) {
    const auto started = std::chrono::steady_clock::now();
    DirtyRegion final_region;
    final_region.graph_revision_before = last_message_graph_revision;
    final_region.graph_revision_after = graph_revision;
    publish("final-ownPhi-publication", round_number.back(), true, false, final_region);
    deferred_final_publication_seconds = publication.seconds + publication.dispatch_seconds;
    deferred_final_publication_wall_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  }

  IntegerVector final_role = last_result["final_inferred_role"];
  List architecture = List::create(
    _["native_engine_calls"] = 1L,
    _["r_round_driver_calls_in_production"] = 0L,
    _["r_move_replay_calls_in_production"] = 0L,
    _["r_search_publications_after_engine_launch"] = 0L,
    _["r_repeated_pass_input_materializations"] = 0L,
    _["repeated_pass_input_materializations"] = 0L,
    _["residual_entry_attempts"] = 0L,
    _["recursive_entry_attempts"] = 0L,
    _["automatic_fallbacks"] = 0L,
    _["defer_intermediate_publication_requested"] = defer_publication_requested,
    _["defer_intermediate_publication_enabled"] = defer_intermediate_publication,
    _["publication_deferral_guard"] = publication_deferral_guard,
    _["deferred_publication_phases"] = deferred_publication_phases,
    _["deferred_final_publication_seconds"] = deferred_final_publication_seconds,
    _["deferred_final_publication_wall_seconds"] = deferred_final_publication_wall_seconds,
    _["final_message_graph_revision"] = last_message_graph_revision,
    _["publication_revision_contract"] = "advance_only_on_actual_publication",
    _["complete_recheck"] = complete_recheck_mode,
    _["complete_recheck_scope"] = complete_recheck_scope,
    _["complete_recheck_obligation_contract"] =
      "actual non-skipped evaluation revision; dependency marks drive pending complete obligations",
    _["first_dirty_wave_done"] = first_dirty_wave_done,
    _["logical_cycle_events"] = logical_cycle_events,
    _["logical_cycle_signature"] =
      "parent matrix plus pending complete obligations");
  DataFrame publications = DataFrame::create(
    _["publication_reason"] = publication_reason,
    _["graph_revision"] = publication_graph,
    _["message_revision"] = publication_message,
    _["total_likelihood_ln"] = publication_total,
    _["parent_signature"] = publication_parent,
    _["message_evidence_signature"] = publication_evidence,
    _["total_likelihood_seconds"] = NumericVector(publication_seconds.size(), 0.0),
    _["publication_core_seconds"] = publication_seconds,
    _["publication_seconds"] = publication_seconds,
    _["publication_mode"] = publication_mode,
    _["dirty_nodes"] = publication_dirty_nodes,
    _["total_nodes"] = publication_total_nodes,
    _["dirty_fraction"] = [&]() {
      NumericVector out(publication_dirty_nodes.size());
      for (R_xlen_t i = 0; i < out.size(); ++i)
        out[i] = publication_total_nodes[i] > 0.0 ?
          publication_dirty_nodes[i] / publication_total_nodes[i] : 0.0;
      return out;
    }(),
    _["dispatch_seconds"] = publication_dispatch_seconds,
    _["nodes_visited"] = publication_nodes_visited,
    _["nodes_recomputed"] = publication_nodes_recomputed,
    _["message_vectors_rebuilt"] = publication_message_vectors,
    _["genotype_state_elements"] = publication_genotype_elements,
    _["transmitted_state_elements"] = publication_transmitted_elements,
    _["likelihood_cells"] = publication_likelihood_cells,
    _["aggregate_updates"] = publication_aggregate_updates,
    _["signature_global_scans"] = publication_signature_global_scans,
    _["copies"] = publication_copies,
    _["allocations"] = publication_allocations,
    _["requested_workers"] = publication_requested_workers,
    _["effective_publication_workers"] = publication_effective_workers,
    _["worker_thread_count"] = publication_worker_thread_count,
    _["publication_job_count"] = publication_job_count,
    _["worker_r_api_calls"] = publication_worker_r_api_calls,
    _["worker_rcpp_api_calls"] = IntegerVector(
      publication_worker_r_api_calls.size(), 0),
    _["finalization_seconds"] = publication_finalization_seconds,
    _["worker_seconds_max"] = publication_worker_seconds_max,
    _["worker_seconds_mean"] = publication_worker_seconds_mean,
    _["worker_imbalance_max_over_mean"] = publication_worker_imbalance,
    _["assignment_policy"] = CharacterVector(
      publication_mode.size(),
      "deterministic_lpt_G2_descending_original_locus_index_lightest_worker"),
    _["reduction_policy"] = CharacterVector(
      publication_mode.size(), "serial_original_locus_order"),
    _["round"] = publication_round);
  DataFrame rounds = DataFrame::create(
    _["round"] = round_number,
    _["accepted_changes"] = round_changes,
    _["accepted_recheck_changes"] = round_recheck_changes,
    _["pending_obligations"] = round_pending,
    _["search_residual_count"] = round_search_residual,
    _["total_assigned_links"] = round_links,
    _["graph_revision_start"] = round_graph_start,
    _["graph_revision_end"] = round_graph_end,
    _["message_revision_start"] = round_message_start,
    _["message_revision_end"] = round_message_end,
    _["publication_performed"] = round_publication,
    _["publication_reason"] = round_publication_reason,
    _["total_likelihood_ln"] = round_total,
    _["total_likelihood_delta_ln"] = round_delta,
    _["total_likelihood_current"] = round_likelihood_current,
    _["likelihood_graph_revision"] = round_likelihood_graph,
    _["publication_deferred"] = round_publication_deferred,
    _["unexplained_likelihood_decrease"] = round_decrease,
    _["source_converged"] = round_source_converged,
    _["recheck_converged"] = round_converged,
    _["parentage_pass_seconds"] = round_pass_seconds,
    _["message_publication_seconds"] = round_publication_seconds,
    _["round_wall_seconds"] = round_wall_seconds,
    _["round_cpu_seconds"] = round_cpu_seconds,
    _["focal_memo_requests"] = NumericVector(round_number.size()),
    _["focal_memo_hits"] = NumericVector(round_number.size()),
    _["focal_memo_admissions"] = NumericVector(round_number.size()));
  DataFrame recheck_phases = DataFrame::create(
    _["round"] = recheck_phase_round,
    _["target_count"] = recheck_phase_targets,
    _["accepted_moves"] = recheck_phase_changes,
    _["pending_after"] = recheck_phase_pending_after,
    _["seconds"] = recheck_phase_seconds,
    _["setup_seconds"] = recheck_phase_setup_seconds,
    _["pass_seconds"] = recheck_phase_pass_seconds,
    _["publication_seconds"] = recheck_phase_publication_seconds);
  DataFrame checks = DataFrame::create(
    _["child"] = check_child, _["old_dam"] = check_old_dam,
    _["old_sire"] = check_old_sire,
    _["check_unavailability_reason"] = check_unavailability_reason,
    _["pair_support_ln"] = check_support,
    _["pair_support_available"] = LogicalVector(check_available.begin(), check_available.end()),
    _["pair_support_proven_rejection"] = LogicalVector(check_proven_rejection.begin(), check_proven_rejection.end()),
    _["pair_support_rejection_upper_ln"] = check_rejection_upper,
    _["pair_hard_rejection"] = LogicalVector(check_hard_rejection.begin(), check_hard_rejection.end()),
    _["raw_null_ln"] = check_raw_null,
    _["raw_dam_only_ln"] = check_raw_dam,
    _["raw_sire_only_ln"] = check_raw_sire,
    _["raw_pair_ln"] = check_raw_pair,
    _["identity_null_ln"] = check_identity_null,
    _["identity_dam_only_ln"] = check_identity_dam,
    _["identity_sire_only_ln"] = check_identity_sire,
    _["identity_pair_ln"] = check_identity_pair,
    _["graph_revision"] = check_graph_revision,
    _["support_boundary"] = CharacterVector(check_child.size(),
      joint_initial_check?"configuration_local_joint_identity_strict_support":
        "frozen_bounded_ownPhi_predictive_fresh_missing_mates"),
    _["active_pair_check_pass"] = check_pass,
    _["individual_side_check_active"] = LogicalVector(check_child.size(), false),
    _["action"] = check_action);

  return List::create(
    _["pair_rescue"] = pair_rescue_output(compiled),
    _["raw_joint_census"] = raw_joint_output(compiled),
    _["schema_version"] = "simp-native-stage1-engine-v1",
    _["status"] = status, _["converged"] = converged,
    _["rounds"] = static_cast<int>(round_number.size()),
    _["final_parents"] = parents,
    _["complete_trio_record_state"] = List::create(
      _["schema_version"] = "complete-trio-dependency-v1",
      _["enabled"] = complete_recheck_enabled && ancestral_rechecks &&
        complete_recheck_confident_ln == R_PosInf,
      _["converged"] = converged,
      _["final_graph_revision"] = graph_revision,
      _["final_parents"] = parents,
      _["last_evaluated_revision"] = wrap(last_complete_evaluation_revision),
      _["recheck_pending"] = IntegerVector(recheck_pending.begin(), recheck_pending.end())),
    _["final_inferred_role"] = final_role,
    _["final_graph_revision"] = graph_revision,
    _["final_message_revision"] = message_revision,
    _["final_focal_order"] = final_focal_order,
    _["final_messages_ln"] = message_matrices(compiled),
    _["final_locus_likelihood_ln"] = locus_likelihood_matrix(
      publication, compiled.loci, compiled.n),
    _["final_message_evidence_signature"] = publication.evidence_signature,
    _["node_likelihood_ln"] = publication.node_likelihood_ln,
    _["total_likelihood_ln"] = publication.total_likelihood_ln,
    _["initial_total_likelihood_ln"] = initial_likelihood,
    _["likelihood_history"] = likelihood_history,
    _["likelihood_history_current"] = likelihood_history_current,
    _["deferred_final_publication_seconds"] = deferred_final_publication_seconds,
    _["deferred_final_publication_wall_seconds"] = deferred_final_publication_wall_seconds,
    _["typed_age_propagation_blocks"] = typed_age_blocks,
    _["accepted_move_blocks"] = move_blocks,
    _["accepted_recheck_move_blocks"] = recheck_move_blocks,
    _["focal_summary_blocks"] = focal_blocks,
    _["recheck_focal_summary_blocks"] = recheck_focal_blocks,
    _["recheck_configuration_diagnostics_blocks"] =
      recheck_configuration_diagnostic_blocks,
    _["recheck_graph_state_diagnostic_blocks"] =
      recheck_graph_state_diagnostic_blocks,
    _["recheck_relationship_provenance_blocks"] =
      recheck_relationship_provenance_blocks,
    _["configuration_diagnostics_blocks"] = configuration_diagnostic_blocks,
    _["graph_state_diagnostic_blocks"] = graph_state_diagnostic_blocks,
    _["relationship_provenance_blocks"] = relationship_provenance_blocks,
    _["configuration_locus_diagnostics_blocks"] =
      configuration_locus_diagnostic_blocks,
    _["round_summary"] = rounds,
    _["recheck_phases"] = recheck_phases,
    _["complete_recheck_child_change_vetoes"] = child_change_vetoes,
    _["complete_recheck_confident_skips"] = DataFrame::create(
      _["round"] = wrap(skip_round), _["focal"] = wrap(skip_focal),
      _["evaluated_revision"] = wrap(skip_evaluated), _["move_revision"] = wrap(skip_move),
      _["identity_margin_ln"] = wrap(skip_identity), _["pair_support_ln"] = wrap(skip_pair),
      _["dam_support_ln"] = wrap(skip_dam), _["sire_support_ln"] = wrap(skip_sire),
      _["full_relationship_ln"] = wrap(skip_relationship)),
    _["message_revision_history"] = publications,
    _["check_pedigree"] = checks,
    _["last_pass"] = last_result,

    _["architecture_telemetry"] = architecture,
    _["final_active_input"] = final_active_input,
    _["compiled_snapshot_persistent_bytes"] = compiled.bytes,
    _["persistent_bytes"] = compiled.bytes);
}
