#pragma once
#include "engine_state.h"
#include "execution_config.h"
namespace pap_internal::search {

struct SearchInputs {
 std::shared_ptr<const pap_internal::ExecutionConfig> config;
  int n = 0;
  std::vector<int> focal_order;
  // Private frozen-state probes: explicit visits evolve ONE live graph. The
  // public focal_order remains a permutation; no state is republished between.
  std::vector<int> diagnostic_visit_order, diagnostic_visit_recheck_order;
  std::vector<int> dam_offsets;
  std::vector<int> sire_offsets;
  std::vector<int> dam_ids;
  std::vector<int> sire_ids;
  std::vector<double> dam_quick;
  std::vector<double> sire_quick;
  std::vector<int> id_order_rank;
  std::vector<int> sex;
  std::vector<int> birth_lo;
  std::vector<int> birth_hi;
  std::vector<int> year_last;
  std::vector<double> min_age;
  std::vector<double> max_age;
  double Tassign_ln = 0.0;
  double identity_margin_ln = 0.0;
  double tol = 0.0;
  double tie_tolerance = 1e-6;
  int max_mendelian_errors = 0;
  // Opt-in unconditional trio screen (sequoia CalcP2 analogue): a parent pair
  // whose trio likelihood ratio falls below this natural-log threshold is
  // screened out like a Mendelian-error failure. -Inf disables it.
  double trio_lr_filter_ln = R_NegInf;
  double trio_lr_margin_ln = 0.0;
  bool trio_lr_margin_recheck_only = false;
  double duo_lr_margin_ln = 0.0;
  bool pair_rescue_enabled=false,pair_rescue_shadow=false,pair_rescue_diagnostics=false,pair_rescue_support_only=false,pair_rescue_calibration=false;
  bool complete_identity_competitors=false;
  bool completion_rank_only=true;
  double pair_rescue_margin_ln=0,pair_rescue_trio_ln=0;
  int pair_rescue_top_k=50;
  double gate_window_ln = 0.0;
  double gate_window_base = 0.0;
  bool gate_window_full = false;
  bool duo_lr_margin_recheck_only = false;
  int start_graph_revision = 0;
  int message_revision = 0;
  bool fixed50 = true;

  bool configuration_diagnostics = false;
  bool configuration_locus_diagnostics = false;
  // Private bounded trajectory capture; these controls never select science.
  int diagnostic_start_position = 1;
  int diagnostic_end_position = 0;
  std::vector<int> initial_dirty_union_node_ids;
  bool recheck_only = false;
  bool private_frozen_decision_audit = false;
  bool private_decision_only_diagnostics = false;
  std::vector<int> recheck_focals;
  int configuration_locus_focal_position = 0;
  // Private one-focal fault-injection seam; public inference never sets it.
  bool private_test_relationship_comparison_exception = false;

  // Saved-state reference observer retains the discarded trial for diagnosis.
  // Public production never supplies these private benchmark controls.
  bool private_complete_message_trial = false;

  bool private_test_raw_comparison_exception = false;
  bool private_test_configuration_post_score_exception = false;
  bool lazy_independent_gates = true;
  // Optional work partition only; known-negative semantics are always shared
  // with the exhaustive corrected reducer. Exhaustive observers stay OFF.

  bool incumbent_full_scoring = false;
  bool joint_identity_objective = false;
  bool joint_support_gates = true;
  bool identity_prescreen = false;
  bool skip_unused_raw = false;
  bool identity_union = false;
  bool raw_support_observer = false;
  bool support_available_baselines = true;
  bool null_joint_identity = true;

  bool raw_joint_census = false;

  bool raw_common_background_legacy_missing_mates = false;

  std::vector<double> raw_common_observation_offsets;

  bool raw_scaled_fast = false;
  bool message_scoring = false;
  bool stored_message_scoring = false;

  bool adaptive_selective_publisher = false;
  double adaptive_selective_threshold = 0.0;
  int progress_parentage_round = 0;
  int progress_locus_interval = 5;
  double progress_heartbeat_seconds = 15.0;
  bool progress_locus_heartbeats = true;
  SEXP progress_callback = R_NilValue;
  SEXP visit_progress_callback = R_NilValue;
  bool benchmark_baseline_mode = false;

  bool benchmark_opt5_diag_mode = false;

  bool phase_timing = false;
  std::string search_loop_engine = "visited_cavity";
  // 0 = the standard relationship reduction; 1 = that same reduction
  // plus reviewed full-only compound relationship members.
  int relationship_mode = 0;
  std::vector<std::array<double, 9>> full_relationship_delta;
  std::vector<unsigned char> full_relationship_orientation_mask;
  std::vector<int> full_relationship_parent_numerator;
  std::vector<int> full_relationship_denominator_alternative;
  std::map<std::string, int> full_relationship_index;
  std::array<int, 9> full_contextual_relationship_index;
  bool bounded_kernel_profile = false;
  int focal_start_position = 1;
  int focal_end_position = 0;
  IntegerMatrix replay_moves;

  bool elide_raw() const;

  static IntegerVector int_field(const List& inputs, const char* name);

  static NumericVector numeric_field(const List& inputs, const char* name);

  SearchInputs(const List& inputs, int n_);

  bool is_recheck_focal(int focal) const;

  void validate() const;
};

}
