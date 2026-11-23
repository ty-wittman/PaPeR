#pragma once
#include "search_inputs.h"
#include "locus_workspace.h"
#include "genotype_checks.h"
namespace pap_internal::search {

static constexpr int kRelationshipCount = 6;
static constexpr int kRelationshipMemberCount =
  2 * kRelationshipCount;
static constexpr int kFullParentNumeratorCount = 3;
static constexpr int kFullContextualRelationshipCount = 9;
static constexpr int kFullContextualOrientationCount = 2;
static constexpr int kFullContextualMemberCount =
  2 * kFullContextualRelationshipCount * kFullContextualOrientationCount;
static constexpr unsigned char FULL_ORIENTATION_FORWARD = 1;
static constexpr unsigned char FULL_ORIENTATION_REVERSE = 2;

static constexpr unsigned char QUESTION_PO_RAW_CONDITIONAL = 0;
static constexpr unsigned char QUESTION_PO_CONTEXT_U_RELATIVE = 1;
static constexpr unsigned char QUESTION_PO_CONTEXT_DISPLAY_REFERENCE = 2;

static constexpr unsigned char QUESTION_PO_CONTEXT_ZERO_U_REFUSAL = 3;

struct SideTotals {
  std::array<double, kRelationshipMemberCount> value;
  std::array<unsigned char, kRelationshipMemberCount> seen;
  // Contextual named motifs retain their raw partition functions and every
  // biological orientation through the complete locus sum. Orientation
  // selection is performed once at the panel boundary, never per locus.
  std::array<double, kRelationshipMemberCount> relationship_partition_value;
  std::array<unsigned char, kRelationshipMemberCount>
    relationship_partition_seen;
  std::array<double, 2 * 22> relationship_orientation_partition_value;
  std::array<unsigned char, 2 * 22> relationship_orientation_partition_seen;
  std::array<double, 2> relationship_po_partition_value;
  std::array<unsigned char, 2> relationship_po_partition_seen;
  // Direct PO is a member of a tested-side request, not a configuration-wide
  // scalar.  For a two-parent proposal, the two requests retain different
  // external evidence boundaries and must remain available through reduction.
  std::array<double, 2> relationship_po_member_value;
  std::array<unsigned char, 2> relationship_po_member_seen;
  std::array<unsigned char, 2> relationship_po_member_scale;
  // PO is the existing bounded configuration score.  HS_PO and PO_HA are
  // aliases of that same score when the proposed-parent motif is present;
  // they never invoke a second parentage likelihood implementation.
  std::array<double, kFullParentNumeratorCount> parent_value;
  std::array<unsigned char, kFullParentNumeratorCount> parent_seen;
  // Contextual denominator members retain candidate side and certified
  // orientation through the complete locus sum.  No maximum or gate is taken
  // at this checkpoint.
  std::array<double, kFullContextualMemberCount> contextual_value;
  std::array<unsigned char, kFullContextualMemberCount> contextual_seen;
  // One structural descriptor per candidate side, shared by every locus in
  // the configuration. Evidence vectors remain locus-local workspace values.
  std::array<RelationshipContext, 2> relationship_context;
  std::array<unsigned char, 2> relationship_context_seen;
  std::array<unsigned char, 2> relationship_common_supported;
  // Conservative panel uncertainty relative to the retained log kernel.
  // This selects numerical recomputation only; it is not a gate tolerance.
  std::array<double, 2> relationship_common_error_bound;
  // A relationship member is accumulated in the serial locus reducer, once
  // for every locus of this exact configuration.  Retain the count in
  // diagnostic builds so provenance captures can distinguish that path from
  // a future cached or pre-configuration member.
  std::array<int, 2> relationship_loci_accumulated;
  // An at-least-once presence flag is not a complete fixed-panel score.
  std::size_t relationship_expected_loci = 0U;
  bool relationship_panel_length_valid = true;
  std::array<int, 2> relationship_po_loci_accumulated{{0, 0}};
  std::array<int, kRelationshipMemberCount> relationship_partition_loci_accumulated{};
  std::array<int, 2 * 22> relationship_orientation_loci_accumulated{};
  // Full likelihoods retain literal hypotheses through the whole-panel sum.
  // Storage is allocated only for active full-mode sides, with one bounded
  // compiled request; it never holds a candidate-times-locus factor graph.
  std::array<unsigned char, 2> full_common_requested;
  std::array<std::shared_ptr<pap_internal::side_likelihood::FullPanelTotals>, 2>
    full_common_panel;
  std::array<std::shared_ptr<std::array<pap_internal::side_likelihood::PanelCell, 25>>, 2>
    full_common_relationship_panel;
  // Optional exact numerical prefix storage, admitted once per focal under
  // a hard aggregate byte budget. No graph or marker factor is duplicated.
  std::array<std::vector<pap_internal::side_staging::Prefix>, 2>
    relationship_stage_prefix;
  SideTotals();
};

std::size_t contextual_member_offset(
    int side, int relationship, int orientation);

struct SideReduction {
  // Unavailable comparisons are not negative relationship evidence. Carry
  // that distinction through the caller rather than encoding it as a gate.
  bool comparison_available = true;
  bool operational_complete = false;
  bool alternative_complete = false;
  double operational_ln = NA_REAL;
  double alternative_ln = NA_REAL;
  double relationship_ln = NA_REAL;
  std::vector<int> operational_ties;
  std::vector<int> alternative_ties;
  // Diagnostics/test-only visibility into the side-local comparisons.  The
  // existing configuration gate continues to consume one scalar reduction;
  // its numerator is selected from the same side as that scalar alternative.
  std::array<double, 2> side_operational_ln{{NA_REAL, NA_REAL}};
  std::array<double, 2> side_alternative_ln{{NA_REAL, NA_REAL}};
  std::array<double, 2> side_relationship_ln{{NA_REAL, NA_REAL}};
  // This may be a certified bound from a complete parent panel plus known
  // nonparent panels even when other alternatives remain unavailable. It
  // never claims that side_relationship_ln is known in that case.
  std::array<double, 2> side_relationship_upper_ln{{NA_REAL, NA_REAL}};
  std::array<unsigned char, 2> side_complete{{0, 0}};
  std::array<double, 2> side_reference_uncertainty{{0.0, 0.0}};
  int selected_side = -1;
};

bool full_relationship_known_rejection(
    const SideReduction& reduction, double Tassign_ln);

SideReduction full_reduce_members(
    const SideTotals* totals, int nparents,
    double null_objective_ln);

struct ConfigurationScore {
  int dam = 0;
  int sire = 0;
  int nparents = 0;
  int status = 0;
  int reason = 0;
  int route = ROUTE_NOT_APPLICABLE;
  int trio_me = -1;
  double objective_ln = R_NegInf;
  double raw_reference_error = 0.0;
  bool raw_message_approximate = false;
  // Ordered finite-marker accumulation can overflow without any model zero
  // or missing marker. A later exact zero may still settle that panel.
  bool raw_accumulation_overflow = false;
  bool raw_comparison_available = true;
  // Raw staging changes only evaluation order; these values remain the
  // original complete ordered Q panel on the frozen candidate background.
  bool raw_precomputed = false;
  bool raw_numerical_cached = false;
  bool raw_support_rejected = false;
  bool incumbent_dam = false, incumbent_sire = false;
  bool full_identity = false;
  bool joint_identity = false;
  double joint_identity_ln = NA_REAL;
  std::size_t raw_loci_accumulated = 0U;
  ::relationship_topology::configuration::PredictiveStatus raw_status =
    ::relationship_topology::configuration::PredictiveStatus::DefinedFinite;
  double dam_support_ln = R_NegInf;
  double sire_support_ln = R_NegInf;
  double pair_support_ln = R_NegInf;
  double identity_margin_ln = R_NegInf;
  bool hard_feasible = false;
  bool hard_valid = false;
  bool trio_me_pass = true;
  bool relationship_pass = false;
  bool relationship_comparison_available = true;
  bool support_comparison_available = true;
  bool independently_rejected = false;
  // Full-mode relationship evidence is an additional gate.  Simp rows retain
  // the default true value and therefore execute the byte-for-byte chooser
  // path that preceded this field.
  bool full_relationship_gate_pass = true;
  bool full_proven_rejected = false;
  // Execution-state flags are not likelihoods and never enter numerical reuse.
  bool full_parent_staged = false;
  bool full_stage_proven_rejected = false;
  bool relationship_prefix_staged = false;
  bool relationship_stage_proven_rejected = false;
  bool gate_deferred = false;
  bool staged_parent_precomputed=false;
  double full_relationship_ln = NA_REAL;
  double full_operational_ln = NA_REAL;
  double full_alternative_ln = NA_REAL;
  bool assignment_confident = false;
  bool eligible_for_retention = false;
  bool eligible_for_addition = false;
  bool eligible_for_replacement = false;
  bool eligible_for_removal = false;
  bool eligible = false;
  bool near_tie = false;
  double cheap_score_ln = R_NegInf;
  // Null in simp mode.  Full mode points to one external accumulator owned by
  // the focal pass, so simp neither allocates nor evaluates full members.
  SideTotals* full_relationship_members = nullptr;
  std::shared_ptr<const SideReduction> native_full_reduction;
};

double decision_objective(const ConfigurationScore& row);

void set_incumbent_policy(const SearchInputs& input,
    const ConfigurationScore& current, ConfigurationScore& row);

struct SideContribution {
  std::array<double, kRelationshipMemberCount> relationship_value;
  std::array<unsigned char, kRelationshipMemberCount> relationship_seen;
  std::array<double, kRelationshipMemberCount> relationship_partition_value;
  std::array<unsigned char, kRelationshipMemberCount>
    relationship_partition_seen;
  std::array<double, 2 * 22> relationship_orientation_partition_value;
  std::array<unsigned char, 2 * 22>
    relationship_orientation_partition_seen;
  std::array<double, 2> relationship_po_partition_value;
  std::array<unsigned char, 2> relationship_po_partition_seen;
  std::array<double, kFullParentNumeratorCount> parent_value;
  std::array<unsigned char, kFullParentNumeratorCount> parent_seen;
  std::array<double, kFullContextualMemberCount> contextual_value;
  std::array<unsigned char, kFullContextualMemberCount> contextual_seen;
  std::array<RelationshipContext, 2> relationship_context;
  std::array<unsigned char, 2> relationship_context_seen;
  std::array<unsigned char, 2> relationship_common_supported;
  std::array<double, 2> relationship_common_error_bound;
  std::array<int, 2> relationship_loci_accumulated;
  std::array<unsigned char, 2> full_common_requested;
  std::array<std::vector<pap_internal::side_likelihood::PartitionValue>, 2>
    full_common_relationship_values;
  std::array<std::vector<pap_internal::side_likelihood::PartitionValue>, 2>
    full_common_values;
  std::size_t full_common_expected_loci = 0U;
  bool full_parent_stage_only = false;
  bool full_alternative_stage_only = false;
  pap_internal::side_staging::Phase relationship_stage_phase =
    pap_internal::side_staging::Phase::complete;
  int relationship_stage_marker = -1;

  SideContribution();
};

struct ConfigurationResult {
  bool scored = false;
  int status = 0;
  int reason = 0;
  int route = ROUTE_NOT_APPLICABLE;
  double locus_objective_ln = R_NegInf;
  double raw_reference_error = 0.0;
  bool raw_message_approximate = false;
  bool raw_comparison_available = true;
  bool raw_only = false;

  // Deferred is distinct from a scored marker or an unavailable calculation.
  bool raw_marker_materialized = true;
  ::relationship_topology::configuration::PredictiveStatus raw_status =
    ::relationship_topology::configuration::PredictiveStatus::DefinedFinite;
  SideContribution full;
};

struct LocusResult {
  int locus_index = 0;
  std::vector<ConfigurationResult> configurations;

  explicit LocusResult(std::size_t configuration_count = 0);

  void reset(int locus, std::size_t configuration_count);
};

struct SingleConfigurationIndices {
  std::vector<int> dam_only_by_parent;
  std::vector<int> sire_only_by_parent;

  int dam_only(int parent) const;

  int sire_only(int parent) const;
};

SingleConfigurationIndices single_configuration_indices(
    const std::vector<ConfigurationScore>& configurations);

bool full_relationship_gate_for_configuration(
    const std::vector<ConfigurationScore>& configurations,
    const std::vector<SideReduction>& reduction,
    int configuration_index, double Tassign_ln);

struct AcceptedMove {
  int child = 0;
  int old_dam = 0;
  int old_sire = 0;
  int new_dam = 0;
  int new_sire = 0;
  double current_objective_ln = R_NegInf;
  double selected_objective_ln = R_NegInf;
  double objective_delta_ln = R_NegInf;
  double dam_support_ln = R_NegInf;
  double sire_support_ln = R_NegInf;
  double pair_support_ln = R_NegInf;
  double identity_margin_ln = R_NegInf;
  int trio_me = -1;
  int graph_revision_before = 0;
  int graph_revision_after = 0;
};

struct VisitSummary {
  int focal = 0;
  int starting_dam = 0;
  int starting_sire = 0;
  int selected_dam = 0;
  int selected_sire = 0;
  double current_objective_ln = R_NegInf;
  double selected_objective_ln = R_NegInf;
  double dam_support_ln = R_NegInf;
  double sire_support_ln = R_NegInf;
  double pair_support_ln = R_NegInf;
  double identity_margin_ln = R_NegInf;
  double full_relationship_ln = NA_REAL;
  int trio_me = -1;
  int trio_me_pass = 1;
  int near_tie = 0;
  int accepted = 0;
  int rejection_reason = 1;
  int graph_revision_before = 0;
  int graph_revision_after = 0;
  int dam_candidate_count = 0;
  int sire_candidate_count = 0;
  int configurations_total = 0;
  int single_configurations_total = 0;
  int pair_configurations_total = 0;
  int single_hard_feasible = 0;
  int pair_hard_feasible = 0;
  int single_trio_me_pass = 0;
  int pair_trio_me_pass = 0;
  int single_relationship_pass = 0;
  int pair_relationship_pass = 0;
  int single_eligible = 0;
  int pair_eligible = 0;
  int configurations_eligible = 0;
  int current_hard_valid = 0;
  int current_assignment_confident = 0;
  int current_retention_eligible = 0;
  int recheck_evaluation = 0;
  int complete_pair_skipped = 0;
  int message_scoring_model = 0;
  int raw_message_approximate = 0;
  double focal_wall_seconds = 0.0;
  double chooser_seconds = 0.0;
  double bounded_score_calls = 0.0, bounded_cavity_requests = 0.0;
  double bounded_cavity_computations = 0.0, bounded_cavity_cache_hits = 0.0;
  double bounded_family_requests = 0.0, bounded_family_computations = 0.0;
  double bounded_family_cache_hits = 0.0, bounded_child_contractions = 0.0;
  double bounded_message_reads = 0.0, bounded_operation_count = 0.0;
  double recursive_visited_cavity_calls = 0.0, recursive_family_to_cavity_calls = 0.0;
};

void merge_visit_counters(
    EngineCounters& target, const EngineCounters& source);

struct RelationshipProvenanceDiagnostic {
  int focal = 0;
  int candidate = 0;
  int tested_parent_side = 0;
  int proposed_dam = 0;
  int proposed_sire = 0;
  int focal_dam = 0;
  int focal_sire = 0;
  int candidate_dam = 0;
  int candidate_sire = 0;
  int fixed_other_parent = 0;
  int has_relevant_context = 0;
  int focal_endpoint_excluded_child = 0;
  int candidate_endpoint_excluded_child = 0;
  int fs_orientation_mask = 0;
  int hs_orientation_mask = 0;
  int gp_orientation_mask = 0;
  int fa_orientation_mask = 0;
  int ha_orientation_mask = 0;
  int loci_accumulated = 0;
  std::uint64_t context_signature = 0;
  std::array<double, kRelationshipCount> partition_ln;
  double proposed_po_partition_ln = R_NegInf;
  std::array<double, 22> orientation_partition_ln;

  RelationshipProvenanceDiagnostic();
};

template<class Cells>
static List audit_operand_panel(
    const pap_internal::side_likelihood::Boundary& boundary,
    const Cells& cells, std::size_t expected, std::size_t loci) { pap_internal::removed_path("whole_audit_operand_panel"); }

struct ConfigurationDiagnostic {
  double identity_objective_ln=NA_REAL;
  int focal_position = 0;
  int focal = 0;
  int graph_revision = 0;
  int configuration_position = 0;
  int dam = 0;
  int sire = 0;
  int nparents = 0;
  int status = 0;
  int reason = 0;
  int route = ROUTE_NOT_APPLICABLE;
  int trio_me = -1;
  double objective_ln = R_NegInf;
  double raw_reference_error = 0.0;
  std::array<double, 2> side_reference_uncertainty{{0.0, 0.0}};
  int raw_comparison_available = 1;
  ::relationship_topology::configuration::PredictiveStatus raw_status =
    ::relationship_topology::configuration::PredictiveStatus::DefinedFinite;
  double dam_support_ln = R_NegInf;
  double sire_support_ln = R_NegInf;
  double pair_support_ln = R_NegInf;
  double identity_margin_ln = R_NegInf;
  int hard_feasible = 0;
  int hard_valid = 0;
  int trio_me_pass = 1;
  int relationship_pass = 0;
  int relationship_comparison_available = 1;
  int full_relationship_gate_pass = 1;
  double full_relationship_ln = NA_REAL;
  double full_operational_ln = NA_REAL;
  double full_alternative_ln = NA_REAL;
  int full_common_boundary = 0;
  int full_limiting_side = -1;
  std::array<double, 2> full_side_relationship_ln{{NA_REAL, NA_REAL}};
  std::string full_parent_model;
  std::string full_alternative_model;
  int assignment_confident = 0;
  int eligible_for_retention = 0;
  int eligible_for_addition = 0;
  int eligible_for_replacement = 0;
  int eligible_for_removal = 0;
  int eligible = 0;
  int near_tie = 0;
  int is_current = 0;
  int is_selected = 0;
  int focal_rejection_reason = 0;
  std::array<double, kRelationshipMemberCount> full_relationship_member_ln;
  std::array<double, kFullParentNumeratorCount> full_parent_numerator_ln;
  std::array<double, kFullContextualMemberCount> full_contextual_member_ln;
  int current_dam_before_test = 0;
  int current_sire_before_test = 0;
  std::array<RelationshipProvenanceDiagnostic, 2> relationship_provenance;
  std::array<unsigned char, 2> relationship_provenance_seen;

  ConfigurationDiagnostic();
};

struct GraphStateDiagnostic {
  int graph_revision = 0;
  std::vector<int> dam;
  std::vector<int> sire;
};

struct ConfigurationLocusDiagnostic {
  int focal_position = 0;
  int focal = 0;
  int graph_revision = 0;
  int configuration_position = 0;
  int locus = 0;
  double objective_ln = R_NegInf;
};

void report_visits(const SearchInputs& input, std::size_t done,
                                std::size_t total, std::size_t moves);

bool is_na(int value);

int assigned_links(const VisitGraph& pass);

struct SearchState {
  const SearchInputs& input;
  std::unique_ptr<pap_internal::typed_age::State> typed_ages;
  pap_internal::typed_age::Bounds visit_age_bounds;
  mutable double age_proposal_seconds=0, age_proposals=0, age_rejections=0;
  double age_move_seconds=0, age_moves=0, age_move_nodes=0, age_move_relaxations=0;
  double age_visit_seconds=0;
  std::vector<int> age_changed_revision;

  void prepare_visit_ages(int child);
  VisitGraph pass;
  std::vector<std::array<int, 2>> role_counts;
  std::vector<char> publication_dirty_union;
  int graph_revision = 0;

  SearchState(EngineState& compiled, const SearchInputs& input_,
                 const IntegerMatrix& parents);

  PedigreeGraphSnapshot graph_snapshot() const;

  bool initial_graph_valid() const;

  bool reaches_focal(int child, int proposed_dam, int proposed_sire,
                     int start) const;

  int proposal_reason(int child, int proposed_dam, int proposed_sire) const;

  DirtyRegion apply_move(int child, int new_dam, int new_sire);
};

struct DeterministicLocusAssignments {
  std::vector<std::vector<int>> by_worker;
  std::vector<int> owner_by_locus;
  std::vector<int> local_position_by_locus;
  std::vector<double> weight_by_locus;
};

void apply_raw_provider_settings(
    const SearchInputs& input, VisitGraph& pass);

struct LocusContext {
  EvaluationState evaluation_state;
  VisitGraph pass;
  // Declared before workspaces so every borrowing workspace dies first.
  std::vector<int> owned_loci;
  std::vector<std::unique_ptr<LocusWorkspace>> workspaces;
  std::vector<LocusResult> locus_results;
  std::size_t result_local_first = 0U, result_local_last = 0U;

  LocusContext(EngineState& compiled,
                          const VisitGraph& authoritative_pass,
                          const PedigreeGraphSnapshot& graph,
                          const SearchInputs& input,
                          const std::vector<int>& owned_loci_);

  void assert_graph_agrees(const PedigreeGraphSnapshot& graph,
                           const char* boundary) const;

  void synchronize_focal(const PedigreeGraphSnapshot& graph, int focal);

  void synchronize_graph(const PedigreeGraphSnapshot& graph,
                         const DirtyRegion& dirty);



};

double log_contrast(double focal, double alternative);

double configuration_contrast(const ConfigurationScore& a,const ConfigurationScore& b,bool decision=false);

bool objective_tie(double first, double second, double tolerance);

bool identity_competitor(const ConfigurationScore& configuration);

std::vector<double> identity_margin_values(
    const std::vector<ConfigurationScore>& configurations,
    bool independently_supported_only = false);

int id_rank(const SearchInputs& input, int id);

std::shared_ptr<const pap_internal::predictive_likelihood::RawBoundary>
build_raw_common_boundary(EngineState&, const SearchState& state,
    int child, const std::vector<ConfigurationScore>& configurations);

struct RelationshipContexts {
  std::array<RelationshipContext, 2> value;
  std::array<unsigned char, 2> seen;
  std::shared_ptr<const pap_internal::predictive_likelihood::RawBoundary> raw_boundary;
  // Transient immutable execution mask; never retained by the numerical cache.
  std::size_t raw_request_index = 0U;

  RelationshipContexts();
};

RelationshipContexts build_side_contexts(
    const SearchState& state, int child,
    const ConfigurationScore& configuration,
    const std::vector<ConfigurationScore>&);

void finalize_contextual_side_members(
    SideTotals* totals);

void align_contextual_operational_display_reference(
    SideTotals* totals, double display_reference_ln);

std::size_t contextual_member_offset(
    int side, int relationship, int orientation);

bool raw_stage_needed(const ConfigurationScore& configuration);

bool relationship_stage_needed(const ConfigurationScore& configuration);

void mark_raw_support_rejections(const SearchInputs& input,
    std::size_t expected_loci, std::vector<ConfigurationScore>& configurations);

bool window_requires_unknown_identity(const SearchInputs& input,
    const std::vector<ConfigurationScore>& configurations);

void restore_cached_raw_parent_alias(const SearchInputs& input,
    std::size_t loci, const ConfigurationScore& configuration,
    SideTotals* totals);

bool evaluate_configuration_locus(
    const SearchInputs& input, VisitGraph& pass,
    LocusWorkspace& workspace, std::size_t locus_index,
    int child, const ConfigurationScore& configuration,
    const RelationshipContexts& contexts,
    ConfigurationResult& result,
    bool raw_only = false);

void reduce_configuration_locus(
    const ConfigurationResult& result,
    ConfigurationScore& configuration);

}
