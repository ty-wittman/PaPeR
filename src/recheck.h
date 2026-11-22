#pragma once
#include <Rcpp.h>
#include <functional>
#include <vector>
namespace pap_internal { struct DirtyRegion; }
namespace pap_internal::search {
struct CandidateSide;
struct MessagePublication;
struct RecheckContext {
  bool complete_recheck_enabled;
  std::function<int()> pending_recheck_count;
  int& graph_revision;
  Rcpp::IntegerMatrix& parents;
  const Rcpp::NumericVector& expected_year;
  const Rcpp::IntegerVector& id_order_rank;
  Rcpp::List& inputs;
  Rcpp::IntegerVector& final_focal_order;
  Rcpp::List& final_active_input;
  const CandidateSide& dam_candidates;
  const CandidateSide& sire_candidates;
  int& message_revision;
  SEXP compiled_sexp;
  std::function<void(const Rcpp::List&)> record_typed_ages;
  std::function<void(const Rcpp::DataFrame&)> record_complete_evaluations;
  int& marking_round;
  std::function<void(const Rcpp::DataFrame&, int)> mark_complete_rechecks_stale;
  bool defer_intermediate_publication;
  int& deferred_publication_phases;
  std::function<void(const char*, int, bool, bool, const DirtyRegion&)> publish;
  MessagePublication& publication;
  std::vector<char>& recheck_pending;
  Rcpp::List& last_result;
  std::vector<Rcpp::List>& recheck_move_blocks;
  std::vector<Rcpp::List>& recheck_focal_blocks;
  std::vector<Rcpp::List>& recheck_configuration_diagnostic_blocks;
  std::vector<Rcpp::List>& recheck_graph_state_diagnostic_blocks;
  std::vector<Rcpp::List>& recheck_relationship_provenance_blocks;
  std::vector<int>& recheck_phase_round;
  std::vector<int>& recheck_phase_targets;
  std::vector<int>& recheck_phase_changes;
  std::vector<int>& recheck_phase_pending_after;
  std::vector<double>& recheck_phase_seconds;
  std::vector<double>& recheck_phase_setup_seconds;
  std::vector<double>& recheck_phase_pass_seconds;
  std::vector<double>& recheck_phase_publication_seconds;
};
int run_recheck_wave(RecheckContext& context, int round);
}
