#include "scoring.h"
namespace pap_internal::search {
static double relationship_partition_contrast(double parent, double alternative);
static double relationship_upper_bound(double parent, double alternative,
                                              double uncertainty);
static pap_internal::panel_rejection::Certificate
common_full_rejection_certificate(
    const RelationshipContext& context,
    const std::array<pap_internal::side_likelihood::PanelCell, 25>& ordinary,
    const pap_internal::side_likelihood::FullPanelTotals& full,
    double Tassign_ln);
static SideReduction reduce_common_full_members(
    const SideTotals& totals, int nparents,
    double display_reference_ln);
static bool full_relationship_passes_threshold(
    const SideReduction& reduction, double Tassign_ln);

static bool edge_feasible(const SearchInputs& input, int child,
                                int parent, int side);
static PedigreeGraphSnapshot make_graph_snapshot(const VisitGraph& pass,
                                                  int graph_revision);

static std::vector<int> native_parent_columns(
    const PedigreeGraphSnapshot& graph);
static void copy_evaluation_settings(
    const VisitGraph& source, VisitGraph& target);
static DeterministicLocusAssignments deterministic_lpt_locus_assignments(
    const EngineState& compiled, int effective_workers);
static SideParentSlots side_parent_slots(
    const VisitGraph& pass, int individual);
static int side_context_parent(
    const VisitGraph& pass, int individual, int side,
    int focal, int candidate);
static SideParentSlots side_context_parent_slots(
    const VisitGraph& pass, int individual, int focal, int candidate);
static bool side_edge_allowed(
    const VisitGraph& pass, const LifeHistoryView* input,
    int child, int proposed_parent, int side,
    int focal, int candidate);
static RelationshipContext
build_side_relationship_context(
    const VisitGraph& pass, int focal, int candidate,
    const LifeHistoryView* = nullptr,
    int proposed_parent_side = 0, int proposed_dam = 0,
    int proposed_sire = 0,
    const std::vector<double>* = nullptr,
    std::uint64_t = 2000000U,
    bool include_full = false,
    const std::shared_ptr<const pap_internal::side_likelihood::Boundary>&
      prepared_boundary = nullptr);
static void evaluate_shared_side_members(
    const SearchInputs& input, LocusWorkspace& workspace, int child,
    const ConfigurationScore& configuration,
    const RelationshipContexts& contexts,
    SideContribution& contribution, bool alternatives_only = false, bool relationship_only = false);
static void present_common_full_members(
    SideTotals& totals, double display_reference_ln);

static std::array<double, 9> reverse_delta(
    const std::array<double, 9>& forward);
static bool can_defer_cached_raw_marker(const SearchInputs& input,
    const VisitGraph& pass, const ConfigurationScore& configuration,
    const RelationshipContexts& contexts, bool raw_only);
static void reduce_shared_side_members(
    const SideContribution& contribution,
    SideTotals& totals);

static double relationship_partition_contrast(double parent, double alternative) {
  if (!finite_or_neg_inf(parent) || !finite_or_neg_inf(alternative) ||
      (is_negative_infinity_ln(parent) && is_negative_infinity_ln(alternative)))
    return NA_REAL;
  return parent - alternative;
}

static double relationship_upper_bound(double parent, double alternative,
                                              double uncertainty) { pap_internal::removed_path("whole_relationship_upper_bound"); }

bool full_relationship_known_rejection(
    const SideReduction& reduction, double Tassign_ln) {
  for (double upper : reduction.side_relationship_upper_ln)
    if (!std::isnan(upper) && upper <= Tassign_ln) return true;
  // Descriptor-only static reduction supplies its exact scalar explicitly.
  // Production side reductions always carry their own error-bound arrays.
  if (reduction.comparison_available && reduction.operational_complete &&
      reduction.alternative_complete && reduction.selected_side < 0 &&
      !std::isnan(reduction.relationship_ln))
    return reduction.relationship_ln <= Tassign_ln;
  return false;
}

static pap_internal::panel_rejection::Certificate
common_full_rejection_certificate(
    const RelationshipContext& context,
    const std::array<pap_internal::side_likelihood::PanelCell, 25>& ordinary,
    const pap_internal::side_likelihood::FullPanelTotals& full,
    double Tassign_ln) {
  namespace proof = pap_internal::panel_rejection;
  proof::Request request;
  request.authority = {context.common_boundary.get(), context.common_panel_owner,
    context.common_graph_revision, context.common_message_revision,
    context.proposed_parent_side};
  request.expected_loci = full.expected_loci;
  request.Tassign_ln = Tassign_ln;
  // PanelCell already covers per-model contractions and original-order
  // summation. A zero envelope means literal retained-reference arithmetic.
  request.retained_reference_arithmetic = true;
  for (const auto& cell : ordinary)
    if (cell.reference_error != 0.0) request.retained_reference_arithmetic = false;
  for (const auto& cell : full.cells)
    if (cell.reference_error != 0.0) request.retained_reference_arithmetic = false;
  proof::OwnedPanels panels{context.common_boundary.get(), &ordinary, &full,
    request.authority};
  return proof::certify(request, panels);
}

static SideReduction reduce_common_full_members(
    const SideTotals& totals, int nparents,
    double display_reference_ln) {
  using namespace pap_internal::side_likelihood;
  SideReduction out;
  out.comparison_available = false;
  const bool sire_single = nparents == 1 &&
    totals.full_common_requested[1] && !totals.full_common_requested[0];
  const int side_start = sire_single ? 1 : 0;
  const int side_end = nparents == 2 ? 2 : side_start + 1;
  bool all_sides_defined = true;
  std::array<AbsoluteSideReduction, 2> reduced;
  for (int side = side_start; side < side_end; ++side) {
    const std::size_t at = static_cast<std::size_t>(side);
    const auto& context = totals.relationship_context[at];
    // Preserve a valid whole-panel upper bound even when unrelated
    // alternatives prevent an exact complete side contrast. Threshold zero
    // selects no policy: the certificate's bound is compared to the current
    // actual Tassign only by the common known-rejection caller below.
    if (totals.full_common_requested[at] && context.common_boundary &&
        context.common_boundary->full_boundary && totals.full_common_panel[at] &&
        totals.full_common_relationship_panel[at]) {
      const auto certificate = common_full_rejection_certificate(context,
        *totals.full_common_relationship_panel[at], *totals.full_common_panel[at], 0.0);
      if (certificate.parent_maximum_complete && !std::isnan(certificate.upper_ln))
        out.side_relationship_upper_ln[at] = certificate.upper_ln;
    }
    if (!totals.full_common_requested[at] || !totals.relationship_common_supported[at] ||
        !context.common_boundary || !context.common_boundary->full_boundary ||
        !totals.full_common_panel[at] || !totals.full_common_relationship_panel[at]) {
      all_sides_defined = false;
      continue;
    }
    reduced[at] = reduce_full_side(*context.common_boundary,
      *totals.full_common_relationship_panel[at], *totals.full_common_panel[at]);
    const auto& value = reduced[at];
    if (value.state != ContrastState::defined) {
      all_sides_defined = false;
      continue;
    }
    // Side-owned absolute partitions determine the contrast. A finite null
    // reference is a display shift only and need not itself be available.
    const auto display = [&](double partition) {
      return R_finite(display_reference_ln) && value.u_relative_available ?
        (is_negative_infinity_ln(partition) ? R_NegInf :
         display_reference_ln + (partition - value.own_u_ln)) : partition;
    };
    out.side_complete[at] = 1;
    out.side_operational_ln[at] = display(value.parent_ln);
    out.side_alternative_ln[at] = display(value.alternative_ln);
    out.side_relationship_ln[at] = value.relationship_ln;
    out.side_relationship_upper_ln[at] = value.relationship_upper_ln;
    if (R_finite(value.relationship_ln))
      out.side_reference_uncertainty[at] = std::max(
        std::fabs(value.relationship_ln - value.relationship_lower_ln),
        std::fabs(value.relationship_upper_ln - value.relationship_ln));
  }
  // Keep all genuine side results above even when the whole AND is unknown.
  // A separate rejection helper consumes only their certified upper bounds.
  if (!all_sides_defined) return out;
  int limiting_side = side_start;
  for (int side = side_start + 1; side < side_end; ++side)
    if (out.side_relationship_ln[static_cast<std::size_t>(side)] <
        out.side_relationship_ln[static_cast<std::size_t>(limiting_side)])
      limiting_side = side;
  const std::size_t at = static_cast<std::size_t>(limiting_side);
  out.comparison_available = true;
  out.operational_complete = out.alternative_complete = true;
  out.selected_side = limiting_side;
  out.operational_ln = out.side_operational_ln[at];
  out.alternative_ln = out.side_alternative_ln[at];
  out.relationship_ln = out.side_relationship_ln[at];
  for (std::size_t id : reduced[at].parent_ties)
    out.operational_ties.push_back(static_cast<int>(id));
  for (std::size_t id : reduced[at].alternative_ties)
    out.alternative_ties.push_back(static_cast<int>(id));
  return out;
}

SideReduction full_reduce_members(
    const SideTotals* totals, int nparents,
    double null_objective_ln) {
  SideReduction out;
  if (nparents == 0) {
    out.operational_complete = out.alternative_complete = true;
    out.operational_ln = out.alternative_ln = null_objective_ln;
    out.relationship_ln = relationship_partition_contrast(
      null_objective_ln, null_objective_ln);
    out.comparison_available = !std::isnan(out.relationship_ln);
    out.operational_ties.push_back(0);
    out.alternative_ties.push_back(0);
    return out;
  }
  if (totals == nullptr) {
    out.comparison_available = false;
    return out;
  }
  if (totals->full_common_requested[0] || totals->full_common_requested[1])
    return reduce_common_full_members(*totals, nparents, null_objective_ln);
  auto reduce = [](const std::vector<std::pair<int, double>>& values,
                   bool& complete, std::vector<int>& ties) {
    if (values.empty()) return NA_REAL;
    complete = true;
    double maximum = R_NegInf;
    for (const auto& value : values) {
      if (!finite_or_neg_inf(value.second)) {
        complete = false;
        return NA_REAL;
      }
      maximum = std::max(maximum, value.second);
    }
    for (const auto& value : values)
      if (value.second == maximum) ties.push_back(value.first);
    std::sort(ties.begin(), ties.end());
    return maximum;
  };
  std::vector<std::pair<int, double>> operational, alternative;
  std::array<std::vector<std::pair<int, double>>, 2> side_alternative;
  std::array<std::vector<int>, 2> side_alternative_ties;
  for (int relationship = 0; relationship < kFullParentNumeratorCount;
       ++relationship)
    if (totals->parent_seen[static_cast<std::size_t>(relationship)])
      operational.emplace_back(100 + relationship,
        totals->parent_value[static_cast<std::size_t>(relationship)]);
  // A single-parent configuration may occupy either typed side.  H_D0 uses
  // the dam cells and H_0S uses the sire cells; H_PP retains both.  The
  // accumulated U cell is always present for an activated candidate, so it is
  // the compact side discriminator without changing the member storage.
  const bool sire_single = nparents == 1 &&
    totals->seen[static_cast<std::size_t>(kRelationshipCount)];
  const int side_start = sire_single ? 1 : 0;
  const int side_end = nparents == 2 ? 2 : side_start + 1;
  for (int side = side_start; side < side_end; ++side)
    if (!totals->relationship_common_supported[static_cast<std::size_t>(side)]) {
      out.comparison_available = false;
      return out;  // Applicable but unscored: never use the raw objective.
    }
  for (int side = side_start; side < side_end; ++side)
    if (totals->relationship_po_member_scale[static_cast<std::size_t>(side)] ==
        QUESTION_PO_CONTEXT_ZERO_U_REFUSAL)
      out.comparison_available = false;
  for (int side = side_start; side < side_end; ++side) {
    for (int relationship = 0; relationship < kRelationshipCount;
         ++relationship) {
      const std::size_t index = static_cast<std::size_t>(
        side * kRelationshipCount + relationship);
      if (totals->seen[index]) {
        alternative.emplace_back(
          side * 1000 + relationship, totals->value[index]);
        side_alternative[static_cast<std::size_t>(side)].emplace_back(
          side * 1000 + relationship, totals->value[index]);
      }
    }
    for (int relationship = 0;
         relationship < kFullContextualRelationshipCount; ++relationship)
      for (int orientation = 0; orientation < kFullContextualOrientationCount;
           ++orientation) {
        const std::size_t index = contextual_member_offset(
          side, relationship, orientation);
        if (totals->contextual_seen[index]) {
          alternative.emplace_back(
            2000 + side * 100 + relationship * 2 + orientation,
            totals->contextual_value[index]);
          side_alternative[static_cast<std::size_t>(side)].emplace_back(
            2000 + side * 100 + relationship * 2 + orientation,
            totals->contextual_value[index]);
        }
      }
  }
  out.operational_ln = reduce(operational, out.operational_complete,
                               out.operational_ties);
  // Ordinary contextual alternatives use each side's own U partition; static
  // ordinary and additional compound alternatives use the frozen pair D9
  // basis. They share a numerical relative-member gauge, which alone does not
  // certify that the compound and contextual ordinary evidence boundaries
  // match. The static path supplies a raw conditional-child PO; production
  // contextual PO is converted to this display reference before reduction.
  // This null constant is a presentation reference for contextual members,
  // never a substitute for either side's own contextual U partition.
  const double alternative_relative_ln = reduce(
    alternative, out.alternative_complete, out.alternative_ties);
  if (out.alternative_complete)
    out.alternative_ln = is_negative_infinity_ln(null_objective_ln) ||
      is_negative_infinity_ln(alternative_relative_ln) ? R_NegInf :
      null_objective_ln + alternative_relative_ln;
  for (int side = side_start; side < side_end; ++side) {
    const std::size_t side_index = static_cast<std::size_t>(side);
    bool side_alternative_complete = false;
    const double side_relative_ln = reduce(
      side_alternative[side_index], side_alternative_complete,
      side_alternative_ties[side_index]);
    if (!side_alternative_complete ||
        !totals->relationship_po_member_seen[side_index])
      continue;
    out.side_complete[side_index] = 1;
    out.side_operational_ln[side_index] =
      totals->relationship_po_member_value[side_index];
    out.side_alternative_ln[side_index] =
      is_negative_infinity_ln(null_objective_ln) ||
      is_negative_infinity_ln(side_relative_ln) ? R_NegInf :
      null_objective_ln + side_relative_ln;
    out.side_relationship_ln[side_index] = relationship_partition_contrast(
      out.side_operational_ln[side_index], out.side_alternative_ln[side_index]);
    out.side_relationship_upper_ln[side_index] = out.side_relationship_ln[side_index];
    if (std::isnan(out.side_relationship_ln[side_index])) {
      out.side_complete[side_index] = 0;
      out.comparison_available = false;
    }
  }
  int limiting_side = -1;
  for (int side = side_start; side < side_end; ++side) {
    const std::size_t side_index = static_cast<std::size_t>(side);
    if (!out.side_complete[side_index]) continue;
    if (limiting_side < 0 || out.side_relationship_ln[side_index] <
        out.side_relationship_ln[static_cast<std::size_t>(limiting_side)])
      limiting_side = side;
  }
  if (limiting_side >= 0) {
    const std::size_t side_index = static_cast<std::size_t>(limiting_side);
    out.selected_side = limiting_side;
    out.operational_complete = true;
    out.alternative_complete = true;
    out.operational_ln = out.side_operational_ln[side_index];
    out.alternative_ln = out.side_alternative_ln[side_index];
    out.relationship_ln = out.side_relationship_ln[side_index];
    out.operational_ties.clear();
    out.operational_ties.push_back(100 + limiting_side);
    out.alternative_ties = side_alternative_ties[side_index];
    return out;
  }
  if (out.operational_complete && out.alternative_complete) {
    out.relationship_ln = relationship_partition_contrast(
      out.operational_ln, out.alternative_ln);
    if (std::isnan(out.relationship_ln)) out.comparison_available = false;
  }
  return out;
}

static bool full_relationship_passes_threshold(
    const SideReduction& reduction, double Tassign_ln) {
  return reduction.comparison_available && reduction.operational_complete &&
    reduction.alternative_complete &&
    reduction.relationship_ln > Tassign_ln;
}

double decision_objective(const ConfigurationScore& row) {
  if(row.joint_identity)return row.joint_identity_ln;
  return row.full_identity ? row.full_operational_ln : row.objective_ln;
}

void set_incumbent_policy(const SearchInputs& input,
    const ConfigurationScore& current, ConfigurationScore& row) {
  row.full_identity = input.incumbent_full_scoring &&
    input.relationship_mode == 1 && current.nparents > 0;
  row.joint_identity=input.joint_identity_objective;
  row.incumbent_dam = row.full_identity && row.dam > 0 && row.dam == current.dam;
  row.incumbent_sire = row.full_identity && row.sire > 0 && row.sire == current.sire;
}

SingleConfigurationIndices single_configuration_indices(
    const std::vector<ConfigurationScore>& configurations) {
  int maximum_parent = 0;
  for (const ConfigurationScore& configuration : configurations)
    maximum_parent = std::max(maximum_parent,
      std::max(configuration.dam, configuration.sire));
  SingleConfigurationIndices out;
  out.dam_only_by_parent.assign(static_cast<std::size_t>(maximum_parent + 1), -1);
  out.sire_only_by_parent.assign(static_cast<std::size_t>(maximum_parent + 1), -1);
  for (std::size_t i = 0; i < configurations.size(); ++i) {
    const ConfigurationScore& configuration = configurations[i];
    if (configuration.nparents != 1) continue;
    const int parent = configuration.dam > 0 ? configuration.dam :
      configuration.sire;
    std::vector<int>& index = configuration.dam > 0 ?
      out.dam_only_by_parent : out.sire_only_by_parent;
    if (index[static_cast<std::size_t>(parent)] >= 0)
      engine_stop("whole-pass configuration set duplicated a single-parent row");
    index[static_cast<std::size_t>(parent)] = static_cast<int>(i);
  }
  return out;
}

bool full_relationship_gate_for_configuration(
    const std::vector<ConfigurationScore>& configurations,
    const std::vector<SideReduction>& reduction,
    int configuration_index, double Tassign_ln) {
  if (configuration_index < 0 ||
      configuration_index >= static_cast<int>(configurations.size()) ||
      reduction.size() != configurations.size())
    engine_stop("native full relationship gate received invalid configuration rows");
  const ConfigurationScore& configuration = configurations[
    static_cast<std::size_t>(configuration_index)];
  if (configuration.nparents == 0) return true;
  const auto passes = [&](int index) {
    return index >= 0 && full_relationship_passes_threshold(
      reduction[static_cast<std::size_t>(index)], Tassign_ln);
  };
  return passes(configuration_index);
}

void merge_visit_counters(EngineCounters& target, const EngineCounters& source) {
  target.message_partition_unavailable += source.message_partition_unavailable;
}

void report_visits(const SearchInputs& input, std::size_t done,
                                std::size_t total, std::size_t moves) {
  if (input.visit_progress_callback == R_NilValue) return;
  Function callback(input.visit_progress_callback);
  callback(IntegerVector::create(input.progress_parentage_round,
    input.recheck_only ? 1 : 0, static_cast<int>(done), static_cast<int>(total),
    static_cast<int>(moves)));
}

bool is_na(int value) {
  return value == NA_INTEGER;
}

static bool edge_feasible(const SearchInputs& input, int child,
                                int parent, int side) {
  if (child < 1 || child > input.n || parent < 1 || parent > input.n ||
      child == parent) return false;
  const int parent_sex = input.sex[parent - 1];
  if (side == 1 && !is_na(parent_sex) && parent_sex == 2) return false;
  if (side == 2 && !is_na(parent_sex) && parent_sex == 1) return false;
  const int child_lo = input.birth_lo[child - 1];
  const int child_hi = input.birth_hi[child - 1];
  const int parent_lo = input.birth_lo[parent - 1];
  const int parent_hi = input.birth_hi[parent - 1];
  const int last = input.year_last[parent - 1];
  if (is_na(child_lo) || is_na(child_hi)) return true;
  if (is_na(parent_lo) || is_na(parent_hi))
    return is_na(last) || child_lo <= last;
  long long effective_child_hi = child_hi;
  if (!is_na(last)) effective_child_hi = std::min<long long>(
    effective_child_hi, last);
  if (effective_child_hi < child_lo) return false;
  const double min_age = input.min_age[side - 1];
  const double max_age = input.max_age[side - 1];
  const long long minimum_difference =
    static_cast<long long>(child_lo) - static_cast<long long>(parent_hi);
  const long long maximum_difference =
    effective_child_hi - static_cast<long long>(parent_lo);
  return maximum_difference >= static_cast<long long>(std::ceil(min_age)) &&
    minimum_difference <= static_cast<long long>(std::floor(max_age));
}

static PedigreeGraphSnapshot make_graph_snapshot(const VisitGraph& pass,
                                                  int graph_revision) {
  return PedigreeGraphSnapshot(
    pass.compiled.n, pass.dam, pass.sire, pass.reverse,
    pass.node_to_evidence_row, graph_revision);
}

int assigned_links(const VisitGraph& pass) {
  int out = 0;
  for (int child = 1; child <= pass.compiled.n; ++child) {
    if (pass.dam[child] > 0) ++out;
    if (pass.sire[child] > 0) ++out;
  }
  return out;
}

static std::vector<int> native_parent_columns(
    const PedigreeGraphSnapshot& graph) {
  const std::size_t n = static_cast<std::size_t>(graph.node_count());
  std::vector<int> out(n * 2U, 0);
  for (int child = 1; child <= graph.node_count(); ++child) {
    out[static_cast<std::size_t>(child - 1)] = graph.parent(child, 1);
    out[n + static_cast<std::size_t>(child - 1)] = graph.parent(child, 2);
  }
  return out;
}

static void copy_evaluation_settings(
    const VisitGraph& source, VisitGraph& target) {
  target.options = source.options;

  target.bounded_messages = source.bounded_messages;

  target.raw_scaled_fast = source.raw_scaled_fast;
  target.message_scoring = source.message_scoring;

  target.evidence_exclusions = source.evidence_exclusions;
}

static DeterministicLocusAssignments deterministic_lpt_locus_assignments(
    const EngineState& compiled, int effective_workers) {
  DeterministicLocusAssignments out;
  out.by_worker.resize(static_cast<std::size_t>(effective_workers));
  out.owner_by_locus.assign(static_cast<std::size_t>(compiled.loci), -1);
  out.local_position_by_locus.assign(
    static_cast<std::size_t>(compiled.loci), -1);
  out.weight_by_locus.assign(static_cast<std::size_t>(compiled.loci), 0.0);
  std::vector<int> locus_order(static_cast<std::size_t>(compiled.loci));
  for (int locus = 0; locus < compiled.loci; ++locus) {
    const double G = static_cast<double>(
      compiled.locus[static_cast<std::size_t>(locus)].G);
    out.weight_by_locus[static_cast<std::size_t>(locus)] = G * G;
    locus_order[static_cast<std::size_t>(locus)] = locus;
  }
  std::sort(locus_order.begin(), locus_order.end(), [&](int left, int right) {
    const double left_weight = out.weight_by_locus[
      static_cast<std::size_t>(left)];
    const double right_weight = out.weight_by_locus[
      static_cast<std::size_t>(right)];
    return left_weight == right_weight ? left < right : left_weight > right_weight;
  });
  std::vector<double> assigned_weight(
    static_cast<std::size_t>(effective_workers), 0.0);
  for (int locus : locus_order) {
    int worker_id = 0;
    for (int candidate = 1; candidate < effective_workers; ++candidate)
      if (assigned_weight[static_cast<std::size_t>(candidate)] <
          assigned_weight[static_cast<std::size_t>(worker_id)])
        worker_id = candidate;
    const std::size_t owner = static_cast<std::size_t>(worker_id);
    out.by_worker[owner].push_back(locus);
    out.owner_by_locus[static_cast<std::size_t>(locus)] = worker_id;
    assigned_weight[owner] += out.weight_by_locus[
      static_cast<std::size_t>(locus)];
  }
  for (std::vector<int>& assignment : out.by_worker)
    std::sort(assignment.begin(), assignment.end());
  for (int worker_id = 0; worker_id < effective_workers; ++worker_id)
    for (std::size_t local = 0;
         local < out.by_worker[static_cast<std::size_t>(worker_id)].size();
         ++local)
      out.local_position_by_locus[static_cast<std::size_t>(out.by_worker[
        static_cast<std::size_t>(worker_id)][local])] =
        static_cast<int>(local);
  return out;
}

void apply_raw_provider_settings(
    const SearchInputs& input, VisitGraph& pass) {
  pass.options = input.config->evaluation;

  pass.raw_scaled_fast = input.raw_scaled_fast;
  pass.message_scoring = input.message_scoring;

}

double log_contrast(double focal, double alternative) {
  if ((!R_finite(focal) && !is_negative_infinity_ln(focal)) ||
      (!R_finite(alternative) && !is_negative_infinity_ln(alternative)) ||
      (is_negative_infinity_ln(focal) && is_negative_infinity_ln(alternative)))
    return NA_REAL;
  if (is_negative_infinity_ln(focal)) return R_NegInf;
  return focal - alternative;
}

double configuration_contrast(const ConfigurationScore& a,const ConfigurationScore& b,bool decision) {
  if(decision&&(a.joint_identity||b.joint_identity))
    return log_contrast(decision_objective(a),decision_objective(b));
  return log_contrast(decision?decision_objective(a):a.objective_ln,decision?decision_objective(b):b.objective_ln);
}

bool objective_tie(double first, double second, double tolerance) {
  if (first == second) return true;
  return R_finite(first) && R_finite(second) &&
    std::fabs(first - second) <= tolerance;
}

bool identity_competitor(const ConfigurationScore& configuration) {
  if (configuration.gate_deferred && configuration.full_identity&&!configuration.joint_identity) return false;
  if (!configuration.hard_feasible || configuration.status != 0 ||
      !configuration.trio_me_pass ||
      configuration.independently_rejected ||
      (configuration.nparents > 0 &&
       (configuration.joint_identity ?
         is_negative_infinity_ln(configuration.joint_identity_ln) : configuration.full_identity ?
         is_negative_infinity_ln(configuration.full_operational_ln) :
         (configuration.raw_comparison_available && !R_finite(configuration.objective_ln)))))
    return false;
  // An unscored relationship comparison is uncertainty, not evidence that
  // this identity has been rejected. Keep it in competition. In particular,
  // never derive this mask from the identity margin that it will define.
  return configuration.nparents == 0 ||
    !configuration.relationship_comparison_available ||
    !configuration.support_comparison_available || configuration.relationship_pass;
}

std::vector<double> identity_margin_values(
    const std::vector<ConfigurationScore>& configurations,
    bool independently_supported_only ) {
  int best_index = -1;
  int second_index = -1;
  for (std::size_t i = 0; i < configurations.size(); ++i) {
    if (independently_supported_only &&
        !identity_competitor(configurations[i])) continue;
    if (!((configurations[i].joint_identity||configurations[i].full_identity) ? !std::isnan(decision_objective(configurations[i])) : configurations[i].raw_comparison_available))
      return std::vector<double>(configurations.size(), NA_REAL);
    const double objective = decision_objective(configurations[i]);
    if (best_index < 0 || objective > decision_objective(configurations[best_index])) {
      second_index = best_index;
      best_index = static_cast<int>(i);
    } else if (second_index < 0 ||
               objective > decision_objective(configurations[second_index])) {
      second_index = static_cast<int>(i);
    }
  }
  std::vector<double> out(configurations.size(), R_PosInf);
  for (std::size_t i = 0; i < configurations.size(); ++i) {
    const int alternative_index =
      static_cast<int>(i) == best_index ? second_index : best_index;
    if (alternative_index >= 0)
      out[i] = log_contrast(decision_objective(configurations[i]),
                               decision_objective(configurations[alternative_index]));
  }
  return out;
}

int id_rank(const SearchInputs& input, int id) {
  return id == 0 ? 0 : input.id_order_rank[id - 1];
}

static SideParentSlots side_parent_slots(
    const VisitGraph& pass, int individual) {
  SideParentSlots out;
  if (individual < 1 || individual > pass.compiled.n) return out;
  out.dam = pass.dam[static_cast<std::size_t>(individual)];
  out.sire = pass.sire[static_cast<std::size_t>(individual)];
  return out;
}

static int side_context_parent(
    const VisitGraph& pass, int individual, int side,
    int focal, int candidate) {
  if (individual < 1 || individual > pass.compiled.n) return 0;
  const int parent = pass.parent(individual, side);
  if ((individual == focal && parent == candidate) ||
      (individual == candidate && parent == focal))
    return 0;
  return parent;
}

static SideParentSlots side_context_parent_slots(
    const VisitGraph& pass, int individual, int focal, int candidate) {
  SideParentSlots out;
  out.dam = side_context_parent(
    pass, individual, 1, focal, candidate);
  out.sire = side_context_parent(
    pass, individual, 2, focal, candidate);
  return out;
}

static bool side_edge_allowed(
    const VisitGraph& pass, const LifeHistoryView* input,
    int child, int proposed_parent, int side,
    int focal, int candidate) { pap_internal::removed_path("whole_relationship_edge_allowed"); }

static RelationshipContext
build_side_relationship_context(
    const VisitGraph& pass, int focal, int candidate,
    const LifeHistoryView* ,
    int proposed_parent_side , int proposed_dam ,
    int proposed_sire ,
    const std::vector<double>* ,
    std::uint64_t ,
    bool include_full ,
    const std::shared_ptr<const pap_internal::side_likelihood::Boundary>&
      prepared_boundary ) {

  RelationshipContext out;
  out.common_panel_owner = &pass.compiled;
  out.common_graph_revision = static_cast<std::uint64_t>(pass.graph_revision);
  out.common_message_revision = static_cast<std::uint64_t>(pass.compiled.message_revision);
  out.focal = focal; out.candidate = candidate;
  out.proposed_parent_side = proposed_parent_side;
  out.proposed_dam = proposed_dam; out.proposed_sire = proposed_sire;
  if (proposed_parent_side < 0 || proposed_parent_side > 2 ||
      proposed_dam < 0 || proposed_sire < 0 ||
      (proposed_parent_side == 1 && proposed_dam != candidate) ||
      (proposed_parent_side == 2 && proposed_sire != candidate))
    engine_stop("ordinary relationship proposal descriptor is malformed");
  for (int side = 1; side <= 2; ++side) {
    if (pass.parent(focal, side) == candidate)
      out.candidate_endpoint_excluded_child = focal;
    if (pass.parent(candidate, side) == focal)
      out.focal_endpoint_excluded_child = candidate;
  }
  out.focal_parents = side_context_parent_slots(pass, focal, focal, candidate);
  if (proposed_parent_side != 0) {
    out.focal_parents = SideParentSlots();
    if (proposed_parent_side == 1) out.focal_parents.sire = proposed_sire;
    else out.focal_parents.dam = proposed_dam;
  }
  out.candidate_parents = side_context_parent_slots(pass, candidate, focal, candidate);
  out.focal_dam_parents = side_parent_slots(pass, out.focal_parents.dam);
  out.focal_sire_parents = side_parent_slots(pass, out.focal_parents.sire);
  out.candidate_dam_parents = side_parent_slots(pass, out.candidate_parents.dam);
  out.candidate_sire_parents = side_parent_slots(pass, out.candidate_parents.sire);
  out.has_relevant_context = out.focal_parents.has_any() || out.candidate_parents.has_any() ||
    out.focal_dam_parents.has_any() || out.focal_sire_parents.has_any() ||
    out.candidate_dam_parents.has_any() || out.candidate_sire_parents.has_any();
  if (proposed_parent_side != 0) {
    auto boundary = prepared_boundary;
    auto* live = pap_internal::native_demand::live;
    const pap_internal::native_demand::Tags tags{
      1, 0, focal, candidate, proposed_parent_side, proposed_dam, proposed_sire};
    if (!boundary && live && live->native_catalog)
      boundary = live->catalog(tags, include_full);
    if (!boundary)
      throw std::logic_error("Stage-1 relationship topology requires its live M_ped2 catalog");
    out.common_boundary = std::move(boundary);
    std::uint64_t signature = UINT64_C(1469598103934665603);
    const auto mix = [&](int value) {
      signature ^= static_cast<std::uint32_t>(value);
      signature *= UINT64_C(1099511628211);
    };
    mix(focal); mix(candidate); mix(proposed_parent_side); mix(proposed_dam); mix(proposed_sire);
    mix(pass.graph_revision);
    for (int id : out.common_boundary->evidence_ids) mix(id);
    for (const auto& model : out.common_boundary->models) {
      mix(static_cast<int>(model.topology.status));
      for (const auto& node : model.topology.nodes) {
        mix(node.typed_id); mix(node.has_birth); mix(node.parents[0]); mix(node.parents[1]);
      }
      if (!model.topology.ready()) continue;
      const int orientation = model.topology.orientation_index;
      if (model.topology.relationship == "FS") out.fs_orientation_mask = 1;
      else if (orientation < 0) continue;
      else if (orientation < 2) out.hs_orientation_mask |= 1U << orientation;
      else if (orientation < 10) out.gp_orientation_mask |= 1U << (orientation - 2);
      else if (orientation < 14) out.fa_orientation_mask |= 1U << (orientation - 10);
      else out.ha_orientation_mask |= 1U << (orientation - 14);
    }
    out.signature = signature;
  }
  return out;
}

std::shared_ptr<const pap_internal::predictive_likelihood::RawBoundary>
build_raw_common_boundary(EngineState&, const SearchState& state,
    int child, const std::vector<ConfigurationScore>& configurations) {
  auto* live = pap_internal::native_demand::live;
  if (!live || !live->native_catalog)
    throw std::logic_error("Stage-1 raw topology requires its live M_ped2 catalog");
  std::vector<std::array<int, 2>> proposals;
  std::vector<unsigned char> needed;
  proposals.reserve(configurations.size()); needed.reserve(configurations.size());
  for (const auto& c : configurations) {
    proposals.push_back({{c.dam, c.sire}});
    needed.push_back(!state.input.lazy_independent_gates ||
      state.input.configuration_diagnostics || state.input.configuration_locus_diagnostics ||
      c.nparents != 2 || c.trio_me_pass);
  }
  return live->raw_catalog(child, proposals, needed);
}

RelationshipContexts build_side_contexts(
    const SearchState& state, int child,
    const ConfigurationScore& configuration,
    const std::vector<ConfigurationScore>&) {

  RelationshipContexts out;
  if (state.input.full_relationship_delta.empty() ||
      configuration.full_relationship_members == nullptr) return out;
  auto* live = pap_internal::native_demand::live;
  if (!live || !live->native_catalog)
    throw std::logic_error("Stage-1 relationship topology requires its live M_ped2 catalog");
  const int candidate_by_side[] = {configuration.dam, configuration.sire};
  for (int side = 0; side < 2; ++side) {
    const int candidate = candidate_by_side[side];
    if (candidate == 0) continue;

    const pap_internal::native_demand::Tags tags{
      1, 0, child, candidate, side + 1, configuration.dam, configuration.sire};
    auto boundary = live->catalog(tags, state.input.relationship_mode == 1);
    if (!boundary) throw std::logic_error("M_ped2 returned no ordinary topology catalog");
    live->observe_catalog(tags, *boundary);

    out.value[side] = build_side_relationship_context(
      state.pass, child, candidate, nullptr, side + 1,
      configuration.dam, configuration.sire, nullptr, 2000000U,
      state.input.relationship_mode == 1, boundary);
    out.seen[static_cast<std::size_t>(side)] = 1;
  }
  return out;
}

static void evaluate_shared_side_members(
    const SearchInputs& input, LocusWorkspace& workspace, int child,
    const ConfigurationScore& configuration,
    const RelationshipContexts& contexts,
    SideContribution& contribution, bool alternatives_only , bool relationship_only ) {
  namespace stage = pap_internal::side_staging;
  const stage::Phase relationship_phase = configuration.relationship_prefix_staged ?
    (alternatives_only ? stage::Phase::remainder : stage::Phase::prefix) : stage::Phase::complete;
  contribution.relationship_stage_phase = relationship_phase;
  contribution.relationship_stage_marker = workspace.locus_index;

  if (input.full_relationship_delta.empty() ||
      configuration.full_relationship_members == nullptr) return;
  static constexpr std::array<const char*, kRelationshipCount>
    relationship_id{{"U", "FS", "HS", "GP", "FA", "HA"}};
  const int candidate_by_side[] = {configuration.dam, configuration.sire};
  for (int side = 0; side < 2; ++side) {
    int candidate = 0;
    const RelationshipContext* context_pointer = nullptr;
    {

      candidate = candidate_by_side[side];
      if (candidate == 0) continue;
      if (!contexts.seen[static_cast<std::size_t>(side)])
        engine_stop("whole locus result is missing ordinary context");
      context_pointer = &contexts.value[static_cast<std::size_t>(side)];
    }

    const RelationshipContext& context = *context_pointer;
    const bool prefix_phase = relationship_phase == stage::Phase::prefix;
    const bool remainder_phase = relationship_phase == stage::Phase::remainder;
    const stage::Prefix* prefix = nullptr;
    if (remainder_phase) {
      const auto& saved = configuration.full_relationship_members->relationship_stage_prefix[
        static_cast<std::size_t>(side)];
      if (workspace.locus_index >= 0 && static_cast<std::size_t>(workspace.locus_index) < saved.size())
        prefix = &saved[static_cast<std::size_t>(workspace.locus_index)];
    }
    const auto relationship_present = [&](int relationship) {
      if ((prefix_phase && relationship >= 2) || (remainder_phase && relationship < 2)) return false;
      if (context.common_boundary == nullptr || relationship == 0) return true;
      const auto& models = context.common_boundary->models;
      if (models.size() != 25U) return false;
      const int begin[] = {1, 2, 3, 5, 13, 17};
      const int end[] = {2, 3, 5, 13, 17, 25};
      for (int m = begin[relationship]; m < end[relationship]; ++m)
        if (context.common_boundary->native_support ? workspace.common_last_status[m] != 0 : models[m].topology.status != ::relationship_topology::TopologyStatus::Inapplicable)
          return true;
      return false;
    };
    {

      contribution.relationship_context[static_cast<std::size_t>(side)] = context;
      contribution.relationship_context_seen[static_cast<std::size_t>(side)] = 1;
      contribution.relationship_loci_accumulated[static_cast<std::size_t>(side)] =
        remainder_phase ? 0 : 1;
    }
    std::array<double, kRelationshipCount> value;
    std::array<double, kRelationshipCount> locus_partition;
    std::array<double, 22> locus_orientation_partition;
    {

      value.fill(R_NegInf);
      locus_partition.fill(R_NegInf);
      locus_orientation_partition.fill(R_NegInf);
    }
    {

        value = workspace.bounded_contextual_side_members(
          context, &locus_partition,
          &locus_orientation_partition, relationship_phase, prefix);
        if (context.common_boundary != nullptr) {
          contribution.relationship_common_supported[static_cast<std::size_t>(side)] =
            workspace.common_last_supported ? 1 : 0;
          double alternative_error = 0.0, maximum_score = 0.0;
          for (std::size_t model = 0; model < 25U; ++model) {
            if (model > 1U)
              alternative_error = std::max(alternative_error,
                workspace.common_last_error_bound[model]);
            if (R_finite(workspace.common_last_partition[model]))
              maximum_score = std::max(maximum_score,
                std::fabs(workspace.common_last_partition[model]));
          }
          const double kernel_error = workspace.common_last_error_bound[0] +
            2.0 * workspace.common_last_error_bound[1] + alternative_error;
          // Cover the original-order marker sum and side-local U display
          // subtractions as well as the individual model contractions.
          contribution.relationship_common_error_bound[static_cast<std::size_t>(side)] =
            kernel_error == 0.0 ? 0.0 : kernel_error +
              32.0 * std::numeric_limits<double>::epsilon() *
                (1.0 + maximum_score) * workspace.pass.compiled.loci;
        }

    }
    {

      for (int relationship = 0;
           relationship < kRelationshipCount; ++relationship) {
        const std::size_t index = static_cast<std::size_t>(
          side * kRelationshipCount + relationship);
        contribution.relationship_value[index] =
          value[static_cast<std::size_t>(relationship)];
        contribution.relationship_seen[index] = relationship_present(relationship) ? 1 : 0;
      }
    }
    if (context.has_relevant_context) {
      const double po_partition = [&]() {

        return workspace.common_last_partition[0];
      }();
      {

        contribution.relationship_po_partition_value[
          static_cast<std::size_t>(side)] = po_partition;
        contribution.relationship_po_partition_seen[
          static_cast<std::size_t>(side)] = remainder_phase ? 0 : 1;
      }
      {

        for (int relationship = 0;
             relationship < kRelationshipCount; ++relationship) {
          const std::size_t index = static_cast<std::size_t>(
            side * kRelationshipCount + relationship);
          const double locus_value = locus_partition[
            static_cast<std::size_t>(relationship)];
          contribution.relationship_partition_value[index] = locus_value;
          contribution.relationship_partition_seen[index] = relationship_present(relationship) ? 1 : 0;
        }
      }
      {

        for (int orientation = 0; orientation < 22; ++orientation) {
          const std::size_t index = static_cast<std::size_t>(
            side * 22 + orientation);
          const double locus_value = locus_orientation_partition[
            static_cast<std::size_t>(orientation)];
          contribution.relationship_orientation_partition_value[index] =
            locus_value;
          contribution.relationship_orientation_partition_seen[index] =
            prefix_phase ? 0 : context.common_boundary == nullptr ? 1 :
            context.common_boundary->models.size() == 25U &&
            (context.common_boundary->native_support ? workspace.common_last_status[orientation+3] != 0 : context.common_boundary->models[static_cast<std::size_t>(orientation + 3)]
              .topology.status != ::relationship_topology::TopologyStatus::Inapplicable) ? 1 : 0;
        }
      }
    }
    if ((configuration.relationship_prefix_staged ||
         contribution.full_common_requested[static_cast<std::size_t>(side)]) &&
        context.common_boundary) {
      using namespace pap_internal::side_likelihood;
      workspace.common_side_partitions(context, relationship_phase, prefix);
      auto& values = contribution.full_common_relationship_values[static_cast<std::size_t>(side)];
      values.resize(25U);
      for (std::size_t model = 0U; model < values.size(); ++model) {
        values[model].state = static_cast<PartitionState>(workspace.common_last_status[model]);
        values[model].log_partition = workspace.common_last_partition[model];
        values[model].reference_error = workspace.common_last_error_bound[model];
        if (values[model].state == PartitionState::unavailable)
          values[model].reason = workspace.common_last_reason[model];
      }
      if (!relationship_only && contribution.full_common_requested[static_cast<std::size_t>(side)] &&
          context.common_boundary->full_boundary)
      contribution.full_common_values[static_cast<std::size_t>(side)] =
        workspace.common_full_partitions(context, alternatives_only ?
          pap_internal::panel_staging::Phase::alternatives : configuration.full_parent_staged ?
          pap_internal::panel_staging::Phase::parents :
          pap_internal::panel_staging::Phase::complete);
    }
  }
}

void finalize_contextual_side_members(
    SideTotals* totals) {
  if (totals == nullptr) return;
  static constexpr std::array<int, 6> orientation_begin{{0, 0, 0, 2, 10, 14}};
  static constexpr std::array<int, 6> orientation_end{{0, 0, 2, 10, 14, 22}};
  for (int side = 0; side < 2; ++side) {
    if (!totals->relationship_context_seen[static_cast<std::size_t>(side)] ||
        !totals->relationship_context[static_cast<std::size_t>(side)]
           .has_relevant_context)
      continue;
    if (!totals->relationship_common_supported[static_cast<std::size_t>(side)])
      continue;
    const std::size_t base = static_cast<std::size_t>(
      side * kRelationshipCount);
    const double partition_u = totals->relationship_partition_value[base];
    const double legacy_u = totals->value[base];
    // The proposal-aware PO member is the operational anchor whenever this
    // side was built from a configuration request.  It has the same U
    // partition as the proposal-aware ordinary alternatives. This expression
    // produces a context-U-relative member, not a raw conditional-child PO.
    // Keep that provenance explicit until the reducer-boundary adapter below.
    const bool zero_u_structural_po_refusal =
      totals->relationship_po_partition_seen[static_cast<std::size_t>(side)] &&
      totals->relationship_partition_seen[base] &&
      is_negative_infinity_ln(partition_u) &&
      is_negative_infinity_ln(totals->relationship_po_partition_value[
        static_cast<std::size_t>(side)]);
    if (zero_u_structural_po_refusal ||
        (totals->relationship_po_partition_seen[static_cast<std::size_t>(side)] &&
         !is_negative_infinity_ln(partition_u) &&
         !is_negative_infinity_ln(legacy_u) &&
         (!is_negative_infinity_ln(totals->relationship_po_partition_value[
           static_cast<std::size_t>(side)]) ||
          (R_finite(partition_u) && R_finite(legacy_u))))) {
      // A scored PO structural zero on a finite U boundary is present,
      // not unavailable. Keep it -Inf so a raw objective or the other
      // parental side cannot replace this impossible contextual hypothesis.
      const double po_partition = totals->relationship_po_partition_value[
        static_cast<std::size_t>(side)];
      const double po_member = is_negative_infinity_ln(po_partition) ?
        R_NegInf : legacy_u + po_partition - partition_u;
      totals->relationship_po_member_value[static_cast<std::size_t>(side)] =
        po_member;
      totals->relationship_po_member_seen[static_cast<std::size_t>(side)] = 1;
      // U=0 prevents a PO/U ratio. This distinct refusal tag bypasses the
      // display adapter; -Inf is only the unchanged gate's rejection sentinel,
      // never a defined contextual PO/U log likelihood ratio for 0/0.
      totals->relationship_po_member_scale[static_cast<std::size_t>(side)] =
        zero_u_structural_po_refusal ? QUESTION_PO_CONTEXT_ZERO_U_REFUSAL :
          QUESTION_PO_CONTEXT_U_RELATIVE;
      // Keep the existing scalar as compatibility telemetry for consumers
      // outside reduction.  full_reduce_members() ignores it whenever
      // the selected alternative has a side-matched contextual PO member.
      if (side == 0 || !totals->relationship_context_seen[0]) {
        totals->parent_value[0] = po_member;
        totals->parent_seen[0] = 1;
        // These compound operational rows are certified aliases of the
        // direct parent topology.  Keep their mechanically shared member on
        // the same configuration-complete boundary as PO.
        for (int member = 1; member < kFullParentNumeratorCount; ++member)
          if (totals->parent_seen[static_cast<std::size_t>(member)])
            totals->parent_value[static_cast<std::size_t>(member)] = po_member;
      }
    }
    for (int relationship = 1;
         relationship < kRelationshipCount; ++relationship) {
      double relation_partition = R_NegInf;
      if (relationship == 1) {
        relation_partition = totals->relationship_partition_value[
          base + static_cast<std::size_t>(relationship)];
      } else {
        for (int orientation =
               orientation_begin[static_cast<std::size_t>(relationship)];
             orientation <
               orientation_end[static_cast<std::size_t>(relationship)];
             ++orientation) {
          const std::size_t index = static_cast<std::size_t>(
            side * 22 + orientation);
          if (totals->relationship_orientation_partition_seen[index])
            relation_partition = std::max(
              relation_partition,
              totals->relationship_orientation_partition_value[index]);
        }
      }
      totals->value[base + static_cast<std::size_t>(relationship)] =
        is_negative_infinity_ln(partition_u) ||
        is_negative_infinity_ln(relation_partition) ||
        is_negative_infinity_ln(legacy_u) ? R_NegInf :
        legacy_u + relation_partition - partition_u;
    }
  }
}

static void present_common_full_members(
    SideTotals& totals, double display_reference_ln) {
  using namespace pap_internal::side_likelihood;
  totals.parent_value.fill(NA_REAL);
  totals.parent_seen.fill(0);
  totals.contextual_value.fill(NA_REAL);
  totals.contextual_seen.fill(0);
  const int nparents = static_cast<int>(totals.full_common_requested[0] != 0) +
    static_cast<int>(totals.full_common_requested[1] != 0);
  const SideReduction reduction = reduce_common_full_members(
    totals, nparents, display_reference_ln);
  if (!reduction.comparison_available) return;
  static constexpr std::array<const char*, kFullContextualRelationshipCount>
    contextual_id{{"HS_HA", "HS_F1C", "HS_FA", "HS_GP", "HS_GGG", "FA_HA",
                   "FS_HA", "FS_H1C", "FS_PARENTS_FS"}};
  totals.parent_value.fill(R_NegInf);
  totals.contextual_value.fill(R_NegInf);
  for (int side = 0; side < 2; ++side) {
    const std::size_t at = static_cast<std::size_t>(side);
    if (!totals.full_common_requested[at]) continue;
    const auto& ordinary = *totals.relationship_context[at].common_boundary;
    const auto& full = *ordinary.full_boundary;
    const auto& panel = *totals.full_common_panel[at];
    const auto& relationship_panel = *totals.full_common_relationship_panel[at];
    const bool relative_defined = relationship_panel[1].state == PartitionState::positive;
    const double own_u = relationship_panel[1].log_partition;
    const auto display = [&](double partition) {
      return relative_defined && R_finite(display_reference_ln) ?
        (is_negative_infinity_ln(partition) ? R_NegInf :
         display_reference_ln + (partition - own_u)) : partition;
    };
    if (side == reduction.selected_side) {
      totals.parent_value[0] = display(relationship_panel[0].log_partition);
      totals.parent_seen[0] = 1;
    }
    for (const auto& descriptor : full.members) {
      if (descriptor.compiled_index >= panel.cells.size()) continue;
      const double partition = panel.cells[descriptor.compiled_index].log_partition;
      if (descriptor.parent_numerator) {
        if (side != reduction.selected_side) continue;
        const int member = descriptor.relationship == "HS_PO" ? 1 :
          descriptor.relationship == "PO_HA" ? 2 : -1;
        if (member < 0) continue;
        const std::size_t index = static_cast<std::size_t>(member);
        totals.parent_value[index] = std::max(totals.parent_value[index], display(partition));
        totals.parent_seen[index] = 1;
      } else {
        for (int member = 0; member < kFullContextualRelationshipCount; ++member) {
          if (descriptor.relationship != contextual_id[static_cast<std::size_t>(member)])
            continue;
          const std::size_t index = contextual_member_offset(
            side, member, descriptor.reverse ? 1 : 0);
          totals.contextual_value[index] = relative_defined ?
            std::max(totals.contextual_value[index],
              is_negative_infinity_ln(partition) ? R_NegInf : partition - own_u) : NA_REAL;
          totals.contextual_seen[index] = 1;
        }
      }
    }
  }
}

void align_contextual_operational_display_reference(
    SideTotals* totals, double display_reference_ln) {
  if (totals == nullptr) return;
  if (totals->full_common_requested[0] || totals->full_common_requested[1]) {
    present_common_full_members(*totals, display_reference_ln);
    return;
  }
  if (!R_finite(display_reference_ln)) return;
  for (int side = 0; side < 2; ++side) {
    const std::size_t side_index = static_cast<std::size_t>(side);
    if (!totals->relationship_po_member_seen[side_index] ||
        totals->relationship_po_member_scale[side_index] !=
          QUESTION_PO_CONTEXT_U_RELATIVE)
      continue;
    const std::size_t base = static_cast<std::size_t>(
      side * kRelationshipCount);
    const double legacy_u = totals->value[base];
    const double partition_u = totals->relationship_partition_value[base];
    const double partition_po = totals->relationship_po_partition_value[side_index];
    const bool structural_zero_po =
      is_negative_infinity_ln(partition_po);
    if (!R_finite(legacy_u) || !R_finite(partition_u) ||
        (!R_finite(partition_po) && !structural_zero_po))
      continue;
    const double relative_po = structural_zero_po ? R_NegInf :
      legacy_u + partition_po - partition_u;
    const double display_po = structural_zero_po ? R_NegInf :
      display_reference_ln + relative_po;
    if (!R_finite(display_po) && !structural_zero_po) continue;
    totals->relationship_po_member_value[side_index] = display_po;
    totals->relationship_po_member_scale[side_index] =
      QUESTION_PO_CONTEXT_DISPLAY_REFERENCE;
    // Preserve the exact side ownership and certified alias rules of the
    // finalizer. The underlying PO/HS_PO/PO_HA likelihoods are unchanged.
    if (side == 0 || !totals->relationship_context_seen[0]) {
      totals->parent_value[0] = display_po;
      for (int member = 1; member < kFullParentNumeratorCount; ++member)
        if (totals->parent_seen[static_cast<std::size_t>(member)])
          totals->parent_value[static_cast<std::size_t>(member)] = display_po;
    }
  }
}

std::size_t contextual_member_offset(
    int side, int relationship, int orientation) {
  return static_cast<std::size_t>(
    (side * kFullContextualRelationshipCount + relationship) *
      kFullContextualOrientationCount + orientation);
}

static std::array<double, 9> reverse_delta(
    const std::array<double, 9>& forward) { pap_internal::removed_path("whole_reverse_delta"); }

bool raw_stage_needed(const ConfigurationScore& configuration) {
  return configuration.status != 2 && (configuration.nparents < 2 ||
    (configuration.hard_feasible && configuration.status == 0 &&
     configuration.trio_me_pass));
}

bool relationship_stage_needed(const ConfigurationScore& configuration) {
  if (configuration.full_stage_proven_rejected || configuration.relationship_stage_proven_rejected ||
      configuration.gate_deferred) return false;
  return !configuration.raw_precomputed ||
    (configuration.nparents > 0 && !configuration.raw_support_rejected &&
     configuration.trio_me_pass);
}

void mark_raw_support_rejections(const SearchInputs& input,
    std::size_t expected_loci, std::vector<ConfigurationScore>& configurations) {
  namespace proof = pap_internal::support_rejection;
  const SingleConfigurationIndices singles =
    single_configuration_indices(configurations);
  int null_index = -1;
  for (std::size_t i = 0; i < configurations.size(); ++i)
    if (configurations[i].nparents == 0) { null_index = static_cast<int>(i); break; }
  if (null_index < 0) engine_stop("raw support staging lost its null row");
  const auto panel = [&](int index) {
    proof::Panel value;
    if (index < 0) return value;
    const ConfigurationScore& configuration = configurations[static_cast<std::size_t>(index)];
    value.ln = configuration.objective_ln;
    value.available = configuration.raw_comparison_available;
    value.complete = expected_loci > 0U && configuration.raw_precomputed &&
      configuration.raw_loci_accumulated == expected_loci && configuration.status != 2;
    return value;
  };
  const proof::Panel null = panel(null_index);
  // Gate window (opt-in, visits ranked by raw objective only): a row whose
  // raw objective trails the best by more than the window keeps no
  // relationship answer and cannot be assigned. Under the joint chooser it
  // triggers no identity work until the window is widened.
  // With gate_window_full the window also applies when an incumbent makes
  // identity use full objectives: the current row is always answered, and a
  // deferred row, having no full objective, leaves the identity competition.
  for (auto& configuration : configurations) configuration.gate_deferred = false;
  if (input.gate_window_ln > 0.0) {
    double best = R_NegInf;
    bool full = false, incumbent_dam = false, incumbent_sire = false;
    for (const auto& configuration : configurations) {
      full = full || configuration.full_identity;
      incumbent_dam = incumbent_dam || configuration.incumbent_dam;
      incumbent_sire = incumbent_sire || configuration.incumbent_sire;
      if (configuration.raw_comparison_available && R_finite(configuration.objective_ln))
        best = std::max(best, configuration.objective_ln);
    }
    if ((!full || input.gate_window_full) && R_finite(best))
    for (auto& configuration : configurations) {
      const bool current = (incumbent_dam ? configuration.incumbent_dam : configuration.dam == 0) &&
        (incumbent_sire ? configuration.incumbent_sire : configuration.sire == 0);
      if (configuration.nparents > 0 && configuration.raw_comparison_available && !(full && current) &&
          !(configuration.objective_ln >= best - input.gate_window_ln))
        configuration.gate_deferred = true;
    }
  }

  for (std::size_t i = 0; i < configurations.size(); ++i) {
    ConfigurationScore& configuration = configurations[i];
    configuration.raw_support_rejected = false;
    if (!configuration.hard_feasible || configuration.status != 0 ||
        configuration.nparents == 0) continue;
    // Gate order is numerical cost dispatch, with unchanged gate definitions.
    // Measured S7 SNP and multiallelic controls favor joint support before
    // their relationship suites. Other suites retain relationship-first work.
    // A window-deferred row cannot trigger either joint route here.
    if(auto* live=pap_internal::native_demand::live;live&&input.joint_support_gates){
      if(!configuration.gate_deferred&&configuration.trio_me_pass&&live->prefer_joint_support_first())
        configuration.raw_support_rejected=pap_internal::joint_support_gate::evaluate(
          configuration.dam,configuration.sire,configuration.incumbent_dam,
          configuration.incumbent_sire,input.Tassign_ln,
          [&](int d,int s){return live->joint_support(d,s);},input.support_available_baselines).rejected;
      continue;
    }
    const proof::Panel focal = panel(static_cast<int>(i));
    if (configuration.incumbent_dam || configuration.incumbent_sire) {
      const auto rejected = [&](const proof::Panel& denominator) {
        return focal.complete && focal.available && denominator.complete && denominator.available &&
          log_contrast(focal.ln, denominator.ln) <= input.Tassign_ln;
      };
      if (configuration.nparents == 2) {
        if (!configuration.incumbent_dam)
          configuration.raw_support_rejected = rejected(panel(singles.sire_only(configuration.sire)));
        if (!configuration.incumbent_sire)
          configuration.raw_support_rejected = configuration.raw_support_rejected ||
            rejected(panel(singles.dam_only(configuration.dam)));
      }
      continue;
    }
    configuration.raw_support_rejected = configuration.nparents == 1 ?
      proof::single(focal, null, configuration.dam > 0 ? 1 : 2, input.Tassign_ln).rejected :
      proof::pair(focal, panel(singles.dam_only(configuration.dam)),
        panel(singles.sire_only(configuration.sire)), null, input.Tassign_ln).rejected;
  }
}

bool window_requires_unknown_identity(const SearchInputs& input,
    const std::vector<ConfigurationScore>& configurations) {
  auto* live=pap_internal::native_demand::live;
  if(!live||!input.joint_identity_objective||input.gate_window_ln<=0)return false;
  const ConfigurationScore* best=nullptr;
  for(const auto& row:configurations)
    if(row.raw_comparison_available&&R_finite(row.objective_ln)&&
       (!best||row.objective_ln>best->objective_ln))best=&row;
  if(!best||best->nparents==0)return false;
  for(const auto& row:configurations){
    if(!row.gate_deferred||!row.hard_feasible||row.status!=0||!row.trio_me_pass||
       row.raw_support_rejected||row.full_stage_proven_rejected||
       row.relationship_stage_proven_rejected||row.full_proven_rejected)continue;
    if(row.native_full_reduction&&full_relationship_known_rejection(
         *row.native_full_reduction,input.Tassign_ln))continue;
    double exact=NA_REAL;
    if(!live->cached_joint_identity(row.dam,row.sire,exact)||std::isnan(exact))return true;
  }
  return false;
}

static bool can_defer_cached_raw_marker(const SearchInputs& input,
    const VisitGraph& pass, const ConfigurationScore& configuration,
    const RelationshipContexts& contexts, bool raw_only) {
  if (!pass.bounded_messages || raw_only || !configuration.raw_numerical_cached ||
      !configuration.raw_precomputed || configuration.raw_loci_accumulated !=
        static_cast<std::size_t>(pass.compiled.loci) || !configuration.raw_comparison_available ||
      input.full_relationship_delta.empty() || configuration.full_relationship_members == nullptr ||
      input.configuration_diagnostics || input.configuration_locus_diagnostics ||
      input.benchmark_baseline_mode ||
      input.private_test_relationship_comparison_exception || input.private_test_raw_comparison_exception ||
      input.private_test_configuration_post_score_exception) return false;
  for (std::size_t side = 0U; side != 2U; ++side)
    if ((side == 0U ? configuration.dam : configuration.sire) > 0 &&
        (!contexts.seen[side] || !contexts.value[side].common_boundary)) return false;
  return true;
}

void restore_cached_raw_parent_alias(const SearchInputs& input,
    std::size_t loci, const ConfigurationScore& configuration,
    SideTotals* totals) {
  namespace cfg = ::relationship_topology::configuration;
  if (!totals || input.relationship_mode == 1 || configuration.nparents == 0 ||
      !configuration.raw_precomputed ||
      configuration.raw_loci_accumulated != loci || !configuration.raw_comparison_available ||
      (configuration.raw_status != cfg::PredictiveStatus::DefinedFinite &&
       configuration.raw_status != cfg::PredictiveStatus::DefinedGeneticZero) ||
      !finite_or_neg_inf(configuration.objective_ln)) return;
  totals->parent_value[0] = configuration.objective_ln;
  totals->parent_seen[0] = 1U;
}

bool evaluate_configuration_locus(
    const SearchInputs& input, VisitGraph& pass,
    LocusWorkspace& workspace, std::size_t locus_index,
    int child, const ConfigurationScore& configuration,
    const RelationshipContexts& contexts,
    ConfigurationResult& result,
    bool raw_only ) {

  if(!raw_only && configuration.staged_parent_precomputed && pass.options.relationship_common_scaled_fast) {
    result=ConfigurationResult();return true;
  }
  NativeErrorScope native_errors;

  {

    result = ConfigurationResult();
    result.scored = true;
    result.raw_only = raw_only;
    result.full.full_parent_stage_only = !raw_only && configuration.full_parent_staged;
    if (pass.bounded_messages) {
      result.raw_comparison_available = false;
      result.raw_status = ::relationship_topology::configuration::PredictiveStatus::Unavailable;
      result.locus_objective_ln = NA_REAL;
    }
    // Raw-only work requests no relationship members. Both normal modes
    // still require the complete fixed marker panel for a relationship proof.
    if (!raw_only) result.full.full_common_expected_loci =
      static_cast<std::size_t>(pass.compiled.loci);
    if (!raw_only && input.relationship_mode == 1) {
      for (int side = 0; side < 2; ++side)
        result.full.full_common_requested[static_cast<std::size_t>(side)] =
          (side == 0 ? configuration.dam : configuration.sire) > 0 ? 1 : 0;
    }
  }

  try {
    {

    }
    int route = ROUTE_NOT_APPLICABLE;

    {

    }

    const double value = [&]() {

      namespace cfg = ::relationship_topology::configuration;
      if (input.elide_raw() || can_defer_cached_raw_marker(input, pass, configuration, contexts, raw_only)) {
        result.raw_marker_materialized = false;
        result.raw_comparison_available = false; // this marker was not scored
        result.raw_status = cfg::PredictiveStatus::Unavailable;
        route = configuration.route;
        return NA_REAL;
      }
      cfg::PredictiveResult raw;
      try {
        if (input.private_test_raw_comparison_exception)
          throw std::runtime_error("Private raw comparison fault");
        // The raw-only prepass owns the exact marker values. Workspace reset
        // deliberately retains this boundary-owned panel; a pointer/index
        // check prevents reuse outside the immutable focal request.
        if (configuration.raw_precomputed && !raw_only &&
            workspace.raw_last_boundary == contexts.raw_boundary &&
            contexts.raw_request_index < workspace.raw_last_values.size())
          raw = workspace.raw_last_values[contexts.raw_request_index];
        else raw = workspace.raw_common_predictive(
          contexts.raw_boundary, contexts.raw_request_index);
      } catch (const std::exception& error) { raw.reason = error.what(); }
      result.raw_status = raw.status;
      result.raw_reference_error = raw.reference_error;
      result.raw_message_approximate = raw.message_approximate;
      result.raw_comparison_available = raw.status == cfg::PredictiveStatus::DefinedFinite ||
        raw.status == cfg::PredictiveStatus::DefinedGeneticZero;
      route = configuration.nparents == 0 ? ROUTE_NOT_APPLICABLE : ROUTE_UNCONNECTED;
      return result.raw_comparison_available ? raw.ln : NA_REAL;
    }();

    result.locus_objective_ln = value;
    if (input.private_test_configuration_post_score_exception)
      throw std::runtime_error("Private post-score construction fault");
    // No other proposal consumes a pair's raw objective. Its independently
    // established ME gate settles eligibility before expensive
    // relationship construction. Null/single baselines remain available for
    // the support of every surviving pair. Exhaustive diagnostics score all.
    if (!raw_only && input.lazy_independent_gates && !input.configuration_diagnostics &&
        !input.configuration_locus_diagnostics &&
        (!configuration.trio_me_pass)) {
      result.locus_objective_ln = value;
      result.route = route;
      return true;
    }
    {

    }

    if (!raw_only) {

      try {
        if (input.private_test_relationship_comparison_exception &&
            (contexts.value[0].common_boundary != nullptr ||
             contexts.value[1].common_boundary != nullptr))
          throw std::runtime_error("private ordinary comparison fault");
        evaluate_shared_side_members(
          input, workspace, child, configuration, contexts, result.full);
      } catch (const std::exception&) {
        const bool common_request =
          contexts.value[0].common_boundary != nullptr ||
          contexts.value[1].common_boundary != nullptr;
        if (!common_request) throw;
        // Failure of the new comparison cannot invalidate an already-scored
        // raw objective or supply genetic evidence for incumbent removal.
        workspace.common_last_supported = false;
        for (int side = 0; side < 2; ++side)
          if ((side == 0 ? configuration.dam : configuration.sire) > 0) {
            result.full.relationship_common_supported[static_cast<std::size_t>(side)] = 0;
            // Preserve the requested side discriminator even if the exception
            // preceded its first member. No value is consumed while unsupported.
            result.full.relationship_seen[static_cast<std::size_t>(
              side * kRelationshipCount)] = 1;
          }

      }

      {

      }

      // Native full proposals use literal models on the same side boundary.
      // Frozen-pair marginals are not substitutes for contextual likelihoods.

      {

      }
    }
    {

      result.locus_objective_ln = value;
      result.route = route;
    }

    {

    }

    {

    }
    {

      if (result.raw_marker_materialized && result.raw_comparison_available && !finite_or_neg_inf(value)) {
        result.status = 2;
        result.reason = 4;
        result.locus_objective_ln = R_NegInf;
        return false;
      }
    }

    return true;
  } catch (const std::exception&) {
    if (pass.bounded_messages) {
      // Numerical/construction failure is an unavailable evaluation, not
      // biological invalidity. Preserve any independently completed raw
      // objective; incomplete relationship questions cannot justify removal.
      result.status = 0;
      result.reason = 4;
      result.full.relationship_common_supported.fill(0);

      return false;
    }
    result.status = 2;
    result.reason = 4;
    result.locus_objective_ln = R_NegInf;
    return false;
  }
}

static void reduce_shared_side_members(
    const SideContribution& contribution,
    SideTotals& totals) {
  if (totals.relationship_expected_loci == 0U)
    totals.relationship_expected_loci = contribution.full_common_expected_loci;
  else if (totals.relationship_expected_loci != contribution.full_common_expected_loci)
    totals.relationship_panel_length_valid = false;
  for (int side = 0; side < 2; ++side) {
    const std::size_t side_index = static_cast<std::size_t>(side);
    totals.full_common_requested[side_index] =
      totals.full_common_requested[side_index] ||
      contribution.full_common_requested[side_index];
    // A failed comparison may precede the first context contribution. Its
    // availability still belongs to the request and must survive reduction.
    // Inactive sides default to supported, so they remain neutral here.
    totals.relationship_common_supported[side_index] =
      totals.relationship_common_supported[side_index] &&
      contribution.relationship_common_supported[side_index];
    totals.relationship_common_error_bound[side_index] +=
      contribution.relationship_common_error_bound[side_index];
    if (contribution.relationship_context_seen[side_index]) {
      if (!totals.relationship_context_seen[side_index]) {
        totals.relationship_context[side_index] =
          contribution.relationship_context[side_index];
        totals.relationship_context_seen[side_index] = 1;
      }
      totals.relationship_loci_accumulated[side_index] +=
        contribution.relationship_loci_accumulated[side_index];
    }
    if (contribution.relationship_stage_phase !=
        pap_internal::side_staging::Phase::complete) {
      namespace stage = pap_internal::side_staging;
      auto& ordinary = totals.full_common_relationship_panel[side_index];
      if (contribution.relationship_context_seen[side_index] && ordinary) {
        if (!stage::append_locus(*ordinary, contribution.full_common_relationship_values[side_index],
              contribution.relationship_stage_phase, totals.relationship_expected_loci))
          totals.relationship_common_supported[side_index] = 0U;
        if (contribution.relationship_stage_phase == stage::Phase::prefix) {
          auto& prefix = totals.relationship_stage_prefix[side_index];
          const int marker = contribution.relationship_stage_marker;
          if (marker >= 0 && static_cast<std::size_t>(marker) < prefix.size())
            prefix[static_cast<std::size_t>(marker)] = stage::capture(
              totals.relationship_context[side_index].common_boundary.get(), marker,
              contribution.full_common_relationship_values[side_index],
              contribution.relationship_common_supported[side_index]);
          else totals.relationship_common_supported[side_index] = 0U;
        }
      }
    }
    if (contribution.full_common_requested[side_index]) {
      using namespace pap_internal::side_likelihood;
      const auto& context = totals.relationship_context[side_index];
      if (context.common_boundary && context.common_boundary->full_boundary) {
        const auto& boundary = *context.common_boundary->full_boundary;
        auto& panel = totals.full_common_panel[side_index];
        auto& relationship_panel = totals.full_common_relationship_panel[side_index];
        if (!panel) {
          panel = std::make_shared<FullPanelTotals>(boundary,
            contribution.full_common_expected_loci);
          if (!relationship_panel) relationship_panel = std::make_shared<std::array<PanelCell, 25>>();
        }
        const auto& contributed_context = contribution.relationship_context[side_index];
        if (!contribution.relationship_context_seen[side_index] ||
            !contributed_context.common_boundary ||
            contributed_context.common_boundary->full_boundary.get() != panel->boundary ||
            panel->expected_loci != contribution.full_common_expected_loci ||
            contribution.full_common_relationship_values[side_index].size() != 25U) {
          panel->failure = "Full locus is missing or belongs to a different side request";
        } else {
          pap_internal::panel_staging::append_locus(*panel,
            contribution.full_common_values[side_index], contribution.full_alternative_stage_only ?
              pap_internal::panel_staging::Phase::alternatives : contribution.full_parent_stage_only ?
              pap_internal::panel_staging::Phase::parents :
              pap_internal::panel_staging::Phase::complete);
          if (contribution.relationship_stage_phase ==
              pap_internal::side_staging::Phase::complete)
            for (std::size_t model = 0U; model < 25U; ++model)
              append_panel_cell((*relationship_panel)[model],
                contribution.full_common_relationship_values[side_index][model]);
        }
      }
    }
    for (int relationship = 0;
         relationship < kRelationshipCount; ++relationship) {
      const std::size_t index = static_cast<std::size_t>(
        side * kRelationshipCount + relationship);
      if (!contribution.relationship_seen[index]) continue;
      double& total = totals.value[index];
      if (!totals.seen[index]) {
        total = contribution.relationship_value[index];
        totals.seen[index] = 1;
      } else if (is_negative_infinity_ln(contribution.relationship_value[index])) {
        total = R_NegInf;
      } else {
        total += contribution.relationship_value[index];
      }
    }
    for (int relationship = 0;
         relationship < kRelationshipCount; ++relationship) {
      const std::size_t index = static_cast<std::size_t>(
        side * kRelationshipCount + relationship);
      if (!contribution.relationship_partition_seen[index]) continue;
      ++totals.relationship_partition_loci_accumulated[index];
      double& total = totals.relationship_partition_value[index];
      if (!totals.relationship_partition_seen[index]) {
        total = contribution.relationship_partition_value[index];
        totals.relationship_partition_seen[index] = 1;
      } else if (!finite_or_neg_inf(total) ||
                 !finite_or_neg_inf(contribution.relationship_partition_value[index])) {
        total = NA_REAL; // Unavailability is absorbing, never a genetic zero.
      } else if (is_negative_infinity_ln(
                   contribution.relationship_partition_value[index])) {
        total = R_NegInf;
      } else if (!is_negative_infinity_ln(total)) {
        total += contribution.relationship_partition_value[index];
      }
    }
    if (contribution.relationship_po_partition_seen[side_index]) {
      ++totals.relationship_po_loci_accumulated[side_index];
      double& total = totals.relationship_po_partition_value[side_index];
      if (!totals.relationship_po_partition_seen[side_index]) {
        total = contribution.relationship_po_partition_value[side_index];
        totals.relationship_po_partition_seen[side_index] = 1;
      } else if (!finite_or_neg_inf(total) ||
                 !finite_or_neg_inf(contribution.relationship_po_partition_value[side_index])) {
        total = NA_REAL;
      } else if (is_negative_infinity_ln(total) ||
                 is_negative_infinity_ln(
                   contribution.relationship_po_partition_value[side_index])) {
        total = R_NegInf;
      } else {
        total += contribution.relationship_po_partition_value[side_index];
      }
    }
    for (int orientation = 0; orientation < 22; ++orientation) {
      const std::size_t index = static_cast<std::size_t>(side * 22 + orientation);
      if (!contribution.relationship_orientation_partition_seen[index]) continue;
      ++totals.relationship_orientation_loci_accumulated[index];
      double& total = totals.relationship_orientation_partition_value[index];
      if (!totals.relationship_orientation_partition_seen[index]) {
        total = contribution.relationship_orientation_partition_value[index];
        totals.relationship_orientation_partition_seen[index] = 1;
      } else if (!finite_or_neg_inf(total) ||
                 !finite_or_neg_inf(contribution.relationship_orientation_partition_value[index])) {
        total = NA_REAL;
      } else if (is_negative_infinity_ln(
                   contribution.relationship_orientation_partition_value[index])) {
        total = R_NegInf;
      } else if (!is_negative_infinity_ln(total)) {
        total += contribution.relationship_orientation_partition_value[index];
      }
    }
  }
}

void reduce_configuration_locus(
    const ConfigurationResult& result,
    ConfigurationScore& configuration) {

  if (!result.scored) return;

  if (result.raw_only) ++configuration.raw_loci_accumulated;
  if (!result.raw_only && configuration.full_relationship_members != nullptr) {
    reduce_shared_side_members(
      result.full, *configuration.full_relationship_members);

  }
  if (result.status != 0) {
    configuration.status = result.status;
    configuration.reason = result.reason;
    configuration.objective_ln = R_NegInf;
    return;
  }
  // Phase two still reduces all actual relationship members above, while
  // the original raw marker values must not be accumulated a second time.
  configuration.raw_message_approximate = configuration.raw_message_approximate ||
    result.raw_message_approximate;
  if (configuration.raw_precomputed && !result.raw_only) return;
  if (!result.raw_comparison_available) {
    configuration.raw_comparison_available = false;
    configuration.raw_accumulation_overflow = false;
    // Undefined conditioning has precedence over unavailability for diagnosis;
    // both dominate a genetic zero at another marker in the panel.
    if (result.raw_status == ::relationship_topology::configuration::PredictiveStatus::ConditioningZero ||
        configuration.raw_status != ::relationship_topology::configuration::PredictiveStatus::ConditioningZero)
      configuration.raw_status = result.raw_status;
    configuration.objective_ln = NA_REAL;
    return;
  }
  using RawStatus = ::relationship_topology::configuration::PredictiveStatus;
  if (!configuration.raw_comparison_available) {
    if (!configuration.raw_accumulation_overflow ||
        result.raw_status != RawStatus::DefinedGeneticZero) return;
    // Every marker so far was scored with positive conditioning support.
    // Unlike an unavailable marker, finite-sum overflow cannot defeat this
    // exact zero certificate. This also makes zero status marker-order invariant.
    configuration.raw_comparison_available = true;
  }
  if (configuration.raw_status == RawStatus::DefinedGeneticZero ||
      result.raw_status == RawStatus::DefinedGeneticZero) {
    configuration.raw_status = RawStatus::DefinedGeneticZero;
    configuration.raw_accumulation_overflow = false;
    configuration.objective_ln = R_NegInf;
    configuration.raw_reference_error = 0.0; // exact structural panel zero
  } else {
    const bool approximate = configuration.raw_reference_error != 0.0 ||
      result.raw_reference_error != 0.0;
    configuration.objective_ln += result.locus_objective_ln;
    if (!R_finite(configuration.objective_ln)) {
      if (approximate) {
        // A nonfinite accelerated sum remains pending until the complete
        // retained panel is retried. It is not a structural genetic zero.
        configuration.raw_reference_error = R_PosInf;
      } else {
        // Finite reference marker scores can overflow their ordered sum.
        // Preserve numerical unavailability rather than inventing a zero.
        configuration.raw_comparison_available = false;
        configuration.raw_accumulation_overflow = true;
        configuration.raw_status = RawStatus::Unavailable;
        configuration.objective_ln = NA_REAL;
        configuration.raw_reference_error = 0.0;
      }
    } else {
      configuration.raw_reference_error = pap_internal::predictive_error_bounds::accumulated_error(
        configuration.objective_ln, configuration.raw_reference_error, result.raw_reference_error);
    }
  }
  configuration.route = result.route;
}

DataFrame raw_joint_output(const EngineState& compiled) {
  const char* names[]={"round","focal","graph_revision","recheck","rows","pairs",
    "raw_unavailable","joint_unavailable","max_abs_ln","contrast_range_ln",
    "worst_dam","worst_sire","raw_null_ln","seconds"};
  List columns;
  for(int k=0;k<14;++k){NumericVector x(compiled.raw_joint_checks.size());
    for(size_t j=0;j<compiled.raw_joint_checks.size();++j)x[j]=compiled.raw_joint_checks[j][k];
    columns[names[k]]=x;
  }
  columns.attr("class")="data.frame";
  columns.attr("row.names")=IntegerVector::create(NA_INTEGER,-int(compiled.raw_joint_checks.size()));
  return as<DataFrame>(columns);
}

List pair_rescue_output(const EngineState& compiled) {
  const auto& visits=compiled.pair_rescue_visits;List v;
  auto count=[&](const char* name,auto field){IntegerVector x(visits.size());
    for(std::size_t i=0;i<visits.size();++i)x[i]=field(visits[i]);v[name]=x;};
#define PAIR_COUNT(name) count(#name,[](const auto& r){return r.name;});
  PAIR_COUNT(round) PAIR_COUNT(focal) PAIR_COUNT(revision) PAIR_COUNT(recheck)
  PAIR_COUNT(singles) PAIR_COUNT(cache_hits) PAIR_COUNT(unknown) PAIR_COUNT(passing)
  PAIR_COUNT(support_shortcuts)
  PAIR_COUNT(retained) PAIR_COUNT(pairs) PAIR_COUNT(kept) PAIR_COUNT(failed_both)
  PAIR_COUNT(rank_pruned) PAIR_COUNT(trio_pruned)
  PAIR_COUNT(configurations_before) PAIR_COUNT(configurations_after) PAIR_COUNT(roots_before) PAIR_COUNT(roots_after)
  PAIR_COUNT(completion_round) PAIR_COUNT(completion_rows)
  PAIR_COUNT(rank_capped_dams) PAIR_COUNT(rank_capped_sires)
  PAIR_COUNT(completion_identity_requests) PAIR_COUNT(completion_identity_hits) PAIR_COUNT(completion_identity_computations)
  PAIR_COUNT(dam_floor_pass) PAIR_COUNT(sire_floor_pass)
#undef PAIR_COUNT
  NumericVector seconds(visits.size());for(std::size_t i=0;i<visits.size();++i)seconds[i]=visits[i].seconds;
  v["single_scoring_seconds"]=seconds;
  v.attr("class")="data.frame";v.attr("row.names")=IntegerVector::create(NA_INTEGER,-int(visits.size()));
  const auto& scores=compiled.pair_rescue_scores;List d;
  auto column=[&](const char* name,auto field){NumericVector x(scores.size());
    for(std::size_t i=0;i<scores.size();++i)x[i]=field(scores[i]);d[name]=x;};
#define PAIR_FIELD(name) column(#name,[](const auto& r){return r.name;});
  PAIR_FIELD(round) PAIR_FIELD(focal) PAIR_FIELD(revision) PAIR_FIELD(recheck) PAIR_FIELD(side)
#undef PAIR_FIELD
#define CANDIDATE_FIELD(name) column(#name,[](const auto& r){return r.candidate.name;});
  CANDIDATE_FIELD(id) CANDIDATE_FIELD(feasible) CANDIDATE_FIELD(incumbent)
  CANDIDATE_FIELD(passing) CANDIDATE_FIELD(retained)
#undef CANDIDATE_FIELD
#define SCORE_FIELD(name) column(#name,[](const auto& r){return r.candidate.score.name;});
  SCORE_FIELD(support) SCORE_FIELD(relationship) SCORE_FIELD(error) SCORE_FIELD(complete) SCORE_FIELD(cached)
  SCORE_FIELD(support_rejected)
  SCORE_FIELD(relationship_skipped) SCORE_FIELD(support_only)
#undef SCORE_FIELD
  d.attr("class")="data.frame";d.attr("row.names")=IntegerVector::create(NA_INTEGER,-int(scores.size()));
  const auto& pairs=compiled.pair_rescue_pairs;List p;
  auto pair_column=[&](const char* name,auto field){NumericVector x(pairs.size());
    for(std::size_t i=0;i<pairs.size();++i)x[i]=field(pairs[i]);p[name]=x;};
#define PAIR_DIAG(name) pair_column(#name,[](const auto& r){return r.name;});
  PAIR_DIAG(round) PAIR_DIAG(focal) PAIR_DIAG(revision) PAIR_DIAG(recheck)
  PAIR_DIAG(dam) PAIR_DIAG(sire) PAIR_DIAG(incumbent) PAIR_DIAG(hard_feasible)
  PAIR_DIAG(trio_me_pass) PAIR_DIAG(candidate_keep) PAIR_DIAG(trio_ln)
#undef PAIR_DIAG
  p.attr("class")="data.frame";p.attr("row.names")=IntegerVector::create(NA_INTEGER,-int(pairs.size()));
  const auto& required=compiled.pair_rescue_requirements;List q;
  auto requirement_column=[&](const char* name,auto field){NumericVector x(required.size());
    for(std::size_t i=0;i<required.size();++i)x[i]=field(required[i]);q[name]=x;};
#define REQUIREMENT_FIELD(name) requirement_column(#name,[](const auto& r){return r.name;});
  REQUIREMENT_FIELD(round) REQUIREMENT_FIELD(focal) REQUIREMENT_FIELD(revision)
  REQUIREMENT_FIELD(recheck) REQUIREMENT_FIELD(dam) REQUIREMENT_FIELD(sire)
  REQUIREMENT_FIELD(dam_rank) REQUIREMENT_FIELD(sire_rank) REQUIREMENT_FIELD(margin)
  REQUIREMENT_FIELD(current) REQUIREMENT_FIELD(selected) REQUIREMENT_FIELD(trio_pass)
#undef REQUIREMENT_FIELD
  q.attr("class")="data.frame";q.attr("row.names")=IntegerVector::create(NA_INTEGER,-int(required.size()));
  return List::create(_["visits"]=v,_["singles"]=d,_["pairs"]=p,_["requirements"]=q);
}

SideTotals::SideTotals() {
    value.fill(R_NegInf);
    seen.fill(0);
    relationship_partition_value.fill(R_NegInf);
    relationship_partition_seen.fill(0);
    relationship_orientation_partition_value.fill(R_NegInf);
    relationship_orientation_partition_seen.fill(0);
    relationship_po_partition_value.fill(R_NegInf);
    relationship_po_partition_seen.fill(0);
    relationship_po_member_value.fill(R_NegInf);
    relationship_po_member_seen.fill(0);
    relationship_po_member_scale.fill(QUESTION_PO_RAW_CONDITIONAL);
    parent_value.fill(R_NegInf);
    parent_seen.fill(0);
    contextual_value.fill(R_NegInf);
    contextual_seen.fill(0);
    relationship_context_seen.fill(0);
    relationship_common_supported.fill(1);
    relationship_common_error_bound.fill(0.0);
    relationship_loci_accumulated.fill(0);
    full_common_requested.fill(0);
  }

SideContribution::SideContribution() {
    relationship_value.fill(R_NegInf);
    relationship_seen.fill(0);
    relationship_partition_value.fill(R_NegInf);
    relationship_partition_seen.fill(0);
    relationship_orientation_partition_value.fill(R_NegInf);
    relationship_orientation_partition_seen.fill(0);
    relationship_po_partition_value.fill(R_NegInf);
    relationship_po_partition_seen.fill(0);
    parent_value.fill(R_NegInf);
    parent_seen.fill(0);
    contextual_value.fill(R_NegInf);
    contextual_seen.fill(0);
    relationship_context_seen.fill(0);
    relationship_common_supported.fill(1);
    relationship_common_error_bound.fill(0.0);
    relationship_loci_accumulated.fill(0);
    full_common_requested.fill(0);
  }

LocusResult::LocusResult(std::size_t configuration_count )
      : configurations(configuration_count) {}

auto LocusResult::reset(int locus, std::size_t configuration_count) -> void  {

    locus_index = locus;
    if (configurations.size() != configuration_count)
      configurations.resize(configuration_count);
    auto* live=pap_internal::native_demand::live;
    if(live && live->engine_result_reset) {
      // Unscored rows are never reduced. Both evaluators fully initialize a
      // selected row before marking it scored. Preserve the original slot
      // count, packet representation and evaluation/reduction order.
      for(auto& result:configurations) {
        result.scored=false;
        // Match the reference's release of owned payloads and boundaries;
        // only the fixed scalar arrays remain stale behind scored=false.
        for(auto& context:result.full.relationship_context)context.common_boundary.reset();
        for(auto& panel:result.full.full_common_relationship_values)
          std::vector<pap_internal::side_likelihood::PartitionValue>().swap(panel);
        for(auto& panel:result.full.full_common_values)
          std::vector<pap_internal::side_likelihood::PartitionValue>().swap(panel);
      }
      live->invalidated_result_rows+=configurations.size();
    } else for (ConfigurationResult& result : configurations)
      result = ConfigurationResult();
  }

auto SingleConfigurationIndices::dam_only(int parent) const -> int  {
    return parent > 0 && parent < static_cast<int>(dam_only_by_parent.size()) ?
      dam_only_by_parent[static_cast<std::size_t>(parent)] : -1;
  }

auto SingleConfigurationIndices::sire_only(int parent) const -> int  {
    return parent > 0 && parent < static_cast<int>(sire_only_by_parent.size()) ?
      sire_only_by_parent[static_cast<std::size_t>(parent)] : -1;
  }

RelationshipProvenanceDiagnostic::RelationshipProvenanceDiagnostic() {
    partition_ln.fill(R_NegInf);
    orientation_partition_ln.fill(R_NegInf);
  }

ConfigurationDiagnostic::ConfigurationDiagnostic() {
    full_relationship_member_ln.fill(R_NegInf);
    full_parent_numerator_ln.fill(R_NegInf);
    full_contextual_member_ln.fill(R_NegInf);
    relationship_provenance_seen.fill(0);
  }

auto SearchState::prepare_visit_ages(int child) -> void  {
    const auto start=std::chrono::steady_clock::now();
    auto cut=typed_ages->propose(child-1,{{-1,-1}});
    if(!cut.feasible)engine_stop("typed age source is infeasible");
    visit_age_bounds=typed_ages->bounds(&cut);
    age_visit_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  }

SearchState::SearchState(EngineState& compiled, const SearchInputs& input_,
                 const IntegerMatrix& parents)
      : input(input_),
        pass(compiled, input_.focal_order.front(), input_.start_graph_revision,
             parents),
        role_counts(static_cast<std::size_t>(input_.n + 1),
                    std::array<int, 2>{{0, 0}}),
        publication_dirty_union(
          input_.adaptive_selective_publisher ?
            static_cast<std::size_t>(input_.n + 1) : 0,
          0),
        graph_revision(input_.start_graph_revision) {

    pass.bounded_messages = input.search_loop_engine == "bounded_messages";

    for (int node : input_.initial_dirty_union_node_ids)
      publication_dirty_union[static_cast<std::size_t>(node)] = 1;

    for (int child = 1; child <= input.n; ++child) {
      if (pass.dam[child] > 0) ++role_counts[pass.dam[child]][0];
      if (pass.sire[child] > 0) ++role_counts[pass.sire[child]][1];
    }
    std::vector<std::array<int,2>> links(input.n);
    for(int v=0;v<input.n;++v)links[v]={{pass.dam[v+1]-1,pass.sire[v+1]-1}};
    typed_ages=std::make_unique<pap_internal::typed_age::State>(
      input.birth_lo,input.birth_hi,input.year_last,
      std::array<double,2>{{input.min_age[0],input.min_age[1]}},
      std::array<double,2>{{input.max_age[0],input.max_age[1]}},links);
    visit_age_bounds=typed_ages->bounds();
    age_changed_revision.assign(input.n+1,-1);
  }

auto SearchState::graph_snapshot() const -> PedigreeGraphSnapshot  {
    return make_graph_snapshot(pass, graph_revision);
  }

auto SearchState::initial_graph_valid() const -> bool  {
    if(!typed_ages->feasible)return false;
    std::vector<int> colour(static_cast<std::size_t>(input.n + 1), 0);
    std::function<bool(int)> visit = [&](int node) {
      if (colour[node] == 1) return true;
      if (colour[node] == 2) return false;
      colour[node] = 1;
      for (int parent : {pass.dam[node], pass.sire[node]})
        if (parent > 0 && visit(parent)) return true;
      colour[node] = 2;
      return false;
    };
    for (int node = 1; node <= input.n; ++node)
      if (visit(node)) return false;
    for (int child = 1; child <= input.n; ++child) {
      if (pass.dam[child] > 0 &&
          !edge_feasible(input, child, pass.dam[child], 1)) return false;
      if (pass.sire[child] > 0 &&
          !edge_feasible(input, child, pass.sire[child], 2)) return false;
    }
    for (int node = 1; node <= input.n; ++node) {
      if (role_counts[node][0] > 0 && role_counts[node][1] > 0) return false;
      const int sx = input.sex[node - 1];
      if (!is_na(sx) && ((sx == 2 && role_counts[node][0] > 0) ||
                               (sx == 1 && role_counts[node][1] > 0)))
        return false;
    }
    return true;
  }

auto SearchState::reaches_focal(int child, int proposed_dam, int proposed_sire,
                     int start) const -> bool  {
    std::vector<char> seen(static_cast<std::size_t>(input.n + 1), 0);
    std::queue<int> queue;
    queue.push(start);
    while (!queue.empty()) {
      const int node = queue.front();
      queue.pop();
      if (node == child) return true;
      if (node < 1 || node > input.n || seen[node]) continue;
      seen[node] = 1;
      const int dam = node == child ? proposed_dam : pass.dam[node];
      const int sire = node == child ? proposed_sire : pass.sire[node];
      if (dam > 0) queue.push(dam);
      if (sire > 0) queue.push(sire);
    }
    return false;
  }

auto SearchState::proposal_reason(int child, int proposed_dam, int proposed_sire) const -> int  {
    if (proposed_dam < 0 || proposed_dam > input.n ||
        proposed_sire < 0 || proposed_sire > input.n) return 1;
    if (proposed_dam > 0 && proposed_dam == proposed_sire) return 2;
    if (proposed_dam == child || proposed_sire == child) return 3;
    const int old_dam = pass.dam[child];
    const int old_sire = pass.sire[child];
    std::vector<int> affected;
    for (int parent : {old_dam, old_sire, proposed_dam, proposed_sire})
      if (parent > 0 && std::find(affected.begin(), affected.end(), parent) == affected.end())
        affected.push_back(parent);
    for (int node : affected) {
      int dam_count = role_counts[node][0];
      int sire_count = role_counts[node][1];
      if (old_dam == node) --dam_count;
      if (old_sire == node) --sire_count;
      if (proposed_dam == node) ++dam_count;
      if (proposed_sire == node) ++sire_count;
      if (dam_count > 0 && sire_count > 0) return 4;
      const int sx = input.sex[node - 1];
      if (!is_na(sx) && ((sx == 2 && dam_count > 0) ||
                               (sx == 1 && sire_count > 0))) return 4;
    }
    if (proposed_dam > 0 && proposed_dam != old_dam &&
        !edge_feasible(input, child, proposed_dam, 1)) return 5;
    if (proposed_sire > 0 && proposed_sire != old_sire &&
        !edge_feasible(input, child, proposed_sire, 2)) return 5;
    if (proposed_dam > 0 && reaches_focal(child, proposed_dam,
                                          proposed_sire, proposed_dam)) return 6;
    if (proposed_sire > 0 && reaches_focal(child, proposed_dam,
                                           proposed_sire, proposed_sire)) return 6;
    const auto started=std::chrono::steady_clock::now();++age_proposals;
    auto ages=typed_ages->propose(child-1,{{proposed_dam-1,proposed_sire-1}});
    age_proposal_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    if(!ages.feasible){++age_rejections;return 5;}
    return 0;
  }

auto SearchState::apply_move(int child, int new_dam, int new_sire) -> DirtyRegion  {
    const auto started = std::chrono::steady_clock::now();
    const int old_dam = pass.dam[child];
    const int old_sire = pass.sire[child];
    const PedigreeGraphSnapshot before = graph_snapshot();
    DirtyRegion dirty;
    dirty.changed_nodes = {child};
    dirty.graph_revision_before = graph_revision;
    if (!publication_dirty_union.empty())
      dirty.publication_nodes = before.descendants(std::vector<int>{child});
    const std::vector<int> old_component = before.affected_component(
      std::vector<int>{child, old_dam, old_sire});
    if (!publication_dirty_union.empty())
      for (int node : dirty.publication_nodes)
        publication_dirty_union[static_cast<std::size_t>(node)] = 1;
    for (int parent : {old_dam, old_sire}) if (parent > 0) {
      if (parent == old_dam) --role_counts[parent][0];
      if (parent == old_sire) --role_counts[parent][1];
    }
    if (new_dam > 0) ++role_counts[new_dam][0];
    if (new_sire > 0) ++role_counts[new_sire][1];
    const auto age_started=std::chrono::steady_clock::now();
    auto age_patch=typed_ages->propose(child-1,{{new_dam-1,new_sire-1}});
    if(!age_patch.feasible)engine_stop("accepted move has infeasible typed ages");
    // A changed difference constraint can alter joint pair feasibility even
    // when every marginal interval is unchanged. Version the non-exact
    // component, not only the endpoints whose numeric bounds tightened.
    for(int v:age_patch.nodes)if(!typed_ages->exact(v))age_changed_revision[v+1]=graph_revision+1;
    typed_ages->commit(child-1,{{new_dam-1,new_sire-1}},age_patch);
    ++age_moves;age_move_nodes+=age_patch.nodes.size();age_move_relaxations+=age_patch.relaxations;
    age_move_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-age_started).count();
    pass.apply_parent_move(child, new_dam, new_sire);
    ++graph_revision;
    pass.graph_revision = graph_revision;
    dirty.graph_revision_after = graph_revision;
    const PedigreeGraphSnapshot after = graph_snapshot();
    dirty.affected_nodes = union_nodes(
      old_component, after.affected_component(
        std::vector<int>{child, old_dam, old_sire, new_dam, new_sire}));
    dirty.normalize();
    if(pap_internal::input_checks::enabled()){if (!dirty.revision_advanced_once())engine_stop("typed parent move produced an invalid dirty-region revision");}

    return dirty;
  }

LocusContext::LocusContext(EngineState& compiled,
                          const VisitGraph& authoritative_pass,
                          const PedigreeGraphSnapshot& graph,
                          const SearchInputs& input,
                          const std::vector<int>& owned_loci_)
      : evaluation_state(),
        pass(compiled, authoritative_pass.focal,
             graph.graph_revision(),
             native_parent_columns(graph),
             evaluation_state.telemetry),
        owned_loci(owned_loci_),
        workspaces(),
        locus_results() {
    copy_evaluation_settings(authoritative_pass, pass);
    apply_raw_provider_settings(input, pass);

    workspaces.reserve(owned_loci.size());
    locus_results.reserve(owned_loci.size());
    for (int locus : owned_loci) {
      workspaces.emplace_back(new LocusWorkspace(
        pass, locus));
      locus_results.emplace_back();

    }
  }

auto LocusContext::assert_graph_agrees(const PedigreeGraphSnapshot& graph,
                           const char* boundary) const -> void  {
    if(!pap_internal::input_checks::enabled())return;

    if (pass.graph_revision != graph.graph_revision())
      engine_stop("whole locus worker graph revision mismatch at %s", boundary);
    for (int node = 1; node <= graph.node_count(); ++node)
      if (pass.dam[static_cast<std::size_t>(node)] != graph.parent(node, 1) ||
          pass.sire[static_cast<std::size_t>(node)] != graph.parent(node, 2))
        engine_stop("whole locus worker parent state mismatch at %s", boundary);
  }

auto LocusContext::synchronize_focal(const PedigreeGraphSnapshot& graph, int focal) -> void  {
    assert_graph_agrees(graph, "focal boundary");
    pass.set_focal(focal, graph.graph_revision());
    for(auto& workspace:workspaces) {
      // Reuse within one graph revision only; acquire still checks immutable
      // locus owner, exact p/T, neutral mask and the original bounded budget.
    }
  }

auto LocusContext::synchronize_graph(const PedigreeGraphSnapshot& graph,
                         const DirtyRegion& dirty) -> void  {
    if(pap_internal::input_checks::enabled()){
      if (pass.graph_revision != dirty.graph_revision_before || graph.graph_revision() != dirty.graph_revision_after || !dirty.revision_advanced_once())engine_stop("whole locus worker dirty-region revision mismatch");}
    for (int node : dirty.changed_nodes)
      pass.apply_parent_move(
        node, graph.parent(node, 1), graph.parent(node, 2));
    pass.graph_revision = graph.graph_revision();
    assert_graph_agrees(graph, "post-move boundary");
  }

RelationshipContexts::RelationshipContexts() { seen.fill(0); }

VisitScorer::VisitScorer(EngineState& compiled,
                       const SearchState& authoritative,
                       const SearchInputs& input) {
    DeterministicLocusAssignments schedule =
      deterministic_lpt_locus_assignments(compiled, effective_workers);
    assignments = std::move(schedule.by_worker);
    owner_by_locus = std::move(schedule.owner_by_locus);
    local_position_by_locus = std::move(schedule.local_position_by_locus);
    weight_by_locus = std::move(schedule.weight_by_locus);

    contexts.reserve(static_cast<std::size_t>(effective_workers));
    for (int worker_id = 0; worker_id < effective_workers; ++worker_id) {
      contexts.emplace_back(new LocusContext(
        compiled, authoritative.pass, authoritative.graph_snapshot(), input,
        assignments[static_cast<std::size_t>(worker_id)]));

    }
  }

auto VisitScorer::result_metadata_bytes() const -> std::size_t  {
    namespace wave = pap_internal::result_storage;
    std::size_t bytes = sizeof(*this);
    bytes = wave::add(bytes, wave::multiply(contexts.capacity(), sizeof(contexts[0])));
    bytes = wave::add(bytes, wave::multiply(assignments.capacity(), sizeof(assignments[0])));
    for (const auto& owned : assignments)
      bytes = wave::add(bytes, wave::multiply(owned.capacity(), sizeof(int)));
    bytes = wave::add(bytes, wave::multiply(owner_by_locus.capacity(), sizeof(int)));
    bytes = wave::add(bytes, wave::multiply(local_position_by_locus.capacity(), sizeof(int)));
    bytes = wave::add(bytes, wave::multiply(weight_by_locus.capacity(), sizeof(double)));
    for (const auto& context : contexts) {
      bytes = wave::add(bytes, sizeof(*context));
      bytes = wave::add(bytes, wave::multiply(context->owned_loci.capacity(), sizeof(int)));
      bytes = wave::add(bytes, wave::multiply(context->workspaces.capacity(), sizeof(context->workspaces[0])));
      bytes = wave::add(bytes, wave::multiply(context->locus_results.capacity(), sizeof(LocusResult)));
    }
    return bytes;
  }

auto VisitScorer::capture_result_storage(const SearchInputs& input, bool completed_wave ) -> void  {
    namespace wave = pap_internal::result_storage;
    using pap_internal::side_likelihood::PartitionValue;
    std::size_t metadata = result_metadata_bytes(), payload = 0U;
    std::size_t raw_metadata = 0U, raw_payload = 0U;
    for (const auto& context : contexts) {
      for (const auto& result : context->locus_results) {
        if (result.configurations.capacity() > 0U) metadata = wave::add(metadata, 2U * sizeof(void*));
        payload = wave::add(payload, wave::multiply(result.configurations.capacity(), sizeof(ConfigurationResult)));
        for (const auto& row : result.configurations)
          for (std::size_t side = 0U; side != 2U; ++side) {
            const std::vector<PartitionValue>* panels[] = {
              &row.full.full_common_relationship_values[side], &row.full.full_common_values[side]};
            for (const auto* panel : panels) {
              if (panel->capacity() > 0U) metadata = wave::add(metadata, 2U * sizeof(void*));
              payload = wave::add(payload, wave::multiply(panel->capacity(), sizeof(PartitionValue)));
              for (const auto& value : *panel) {
                const auto bytes = wave::heap_string_bytes(value.reason);
                payload = wave::add(payload, bytes);
                if (bytes > 0U) metadata = wave::add(metadata, 2U * sizeof(void*));
              }
            }
          }
      }
      if (!completed_wave) for (const auto& workspace : context->workspaces) {
        raw_metadata = wave::add(raw_metadata, sizeof(workspace->raw_last_values) + sizeof(workspace->raw_last_boundary));
        if (workspace->raw_last_values.capacity() > 0U) raw_metadata = wave::add(raw_metadata, 2U * sizeof(void*));
        raw_payload = wave::add(raw_payload, wave::multiply(workspace->raw_last_values.capacity(),
          sizeof(::relationship_topology::configuration::PredictiveResult)));
        for (const auto& value : workspace->raw_last_values) {
          const auto bytes = wave::heap_string_bytes(value.reason);
          raw_payload = wave::add(raw_payload, bytes);
          if (bytes > 0U) raw_metadata = wave::add(raw_metadata, 2U * sizeof(void*));
        }
      }
    }
    auto& telemetry = primary().pass.compiled.telemetry;

    const std::size_t total = wave::add(metadata, payload);

  }

auto VisitScorer::score_raw_panel(const SearchInputs& input, int child,
      const std::vector<RelationshipContexts>& side_contexts,
      std::vector<ConfigurationScore>& configurations,
      const std::function<void(int)>& progress ) -> void  {
    namespace cfg = ::relationship_topology::configuration;
    if (input.elide_raw()) {
      for (auto& configuration : configurations) {
        // This is a completed dispatch phase, not a computed Q panel.
        configuration.raw_precomputed = true;
        configuration.raw_numerical_cached = false;
        configuration.raw_loci_accumulated = 0U;
        configuration.raw_reference_error = 0.0;
        configuration.raw_message_approximate = false;
        configuration.raw_accumulation_overflow = false;
        configuration.raw_comparison_available = false;
        configuration.raw_status = cfg::PredictiveStatus::Unavailable;
        configuration.objective_ln = NA_REAL;
        configuration.raw_support_rejected = false;
      }
      mark_raw_support_rejections(input, owner_by_locus.size(), configurations);
      return;
    }
    bool pending_raw = false;
    for (ConfigurationScore& configuration : configurations) {
      configuration.raw_precomputed = true;
      configuration.raw_support_rejected = false;
      if (configuration.raw_numerical_cached) continue;
      pending_raw = pending_raw || raw_stage_needed(configuration);
      configuration.raw_loci_accumulated = 0U;
      configuration.raw_reference_error = 0.0;
      configuration.raw_message_approximate = false;
      configuration.raw_accumulation_overflow = false;
      configuration.raw_comparison_available = raw_stage_needed(configuration);
      configuration.raw_status = configuration.raw_comparison_available ?
        cfg::PredictiveStatus::DefinedFinite : cfg::PredictiveStatus::Unavailable;
      configuration.objective_ln = configuration.raw_comparison_available ? 0.0 : NA_REAL;
    }
    if (!pending_raw) {
      // All required numerical panels are complete. Still redo support/margins
      // on the fresh current configuration universe before any authorization.
      refine_raw_panels(input, side_contexts, configurations);
      mark_raw_support_rejections(input, owner_by_locus.size(), configurations);
      return;
    }
    const double declines_before=primary().pass.telemetry.message_partition_unavailable;
    const auto failed=[&]() {return input.message_scoring && !input.private_complete_message_trial &&
      primary().pass.telemetry.message_partition_unavailable > declines_before;};
    {
      LocusContext& context = primary();
      LocusResult& result = context.locus_results.front();
      const bool compact=pap_internal::native_demand::live&&pap_internal::native_demand::live->engine_compact_results&&!input.configuration_diagnostics&&!input.configuration_locus_diagnostics&&!input.message_scoring;
      for (int locus = 0; locus < static_cast<int>(owner_by_locus.size()); ++locus) {
        const std::size_t local = static_cast<std::size_t>(
          local_position_by_locus[static_cast<std::size_t>(locus)]);
        LocusWorkspace& workspace = *context.workspaces[local];
        result.reset(locus, compact?1:configurations.size());
        if(compact)pap_internal::native_demand::live->compact_packets_avoided+=configurations.size()-1;
        for (std::size_t i = 0; i < configurations.size(); ++i) {
          ConfigurationScore& configuration = configurations[i];
          if (configuration.raw_numerical_cached || !raw_stage_needed(configuration)) continue;
          evaluate_configuration_locus(input, context.pass, workspace,
            static_cast<std::size_t>(locus), child, configuration, side_contexts[i],
            result.configurations[compact?0:i], true);
          if (failed()) return;
          reduce_configuration_locus(result.configurations[compact?0:i], configuration);
        }

        if (progress) progress(locus + 1);
      }
    }
    // No running-marker LLR is consulted. The proof consumes complete Q
    // panels, including all null/single counterfactuals and status precedence.
    capture_result_storage(input);
    refine_raw_panels(input, side_contexts, configurations);
    mark_raw_support_rejections(input, owner_by_locus.size(), configurations);
  }

auto VisitScorer::refine_raw_panels(const SearchInputs& input,
      const std::vector<RelationshipContexts>& side_contexts,
      std::vector<ConfigurationScore>& configurations) -> void  {
    namespace guard = pap_internal::predictive_error_bounds;
    namespace cfg = ::relationship_topology::configuration;
    auto& compiled=primary().pass.compiled;
    bool any = false;
    for (const auto& value : configurations)
      any = any || (value.raw_reference_error != 0.0 && value.raw_comparison_available);
    if (!any) return;
    const auto started = input.bounded_kernel_profile ? std::chrono::steady_clock::now() :
      std::chrono::steady_clock::time_point();
    std::vector<unsigned char> mask;
    try {
      std::vector<guard::Estimate> estimates;
      estimates.reserve(configurations.size());
      for (const auto& value : configurations) estimates.push_back({value.objective_ln,
        value.raw_reference_error, value.raw_comparison_available &&
          value.raw_loci_accumulated == owner_by_locus.size() && value.status == 0});

      mask = guard::refinement_mask(estimates, {input.Tassign_ln, input.identity_margin_ln,
        input.tol, input.tie_tolerance, -input.tie_tolerance});
    } catch (const std::bad_alloc&) {
      // Declining an optional guard optimization does not decline scientific
      // scoring. Recompute every uncertain complete Q with retained arithmetic.
    } catch (const std::length_error&) {}
    bool refined = false;
    for (std::size_t i = 0U; i < configurations.size(); ++i) {
      auto& value = configurations[i];
      const bool pending = value.raw_reference_error != 0.0 && value.raw_comparison_available &&
        value.raw_loci_accumulated == owner_by_locus.size() && value.status == 0;
      if (!pending || (!mask.empty() && !mask[i] && R_finite(value.objective_ln))) continue;
      ConfigurationScore replacement = value;
      replacement.objective_ln = 0.0;
      replacement.raw_reference_error = 0.0;
      replacement.raw_message_approximate = false;
      replacement.raw_accumulation_overflow = false;
      replacement.raw_loci_accumulated = 0U;
      replacement.raw_comparison_available = true;
      replacement.raw_status = cfg::PredictiveStatus::DefinedFinite;
      replacement.raw_support_rejected = false;
      for (int locus = 0; locus < static_cast<int>(owner_by_locus.size()); ++locus) {
        const int owner = owner_by_locus[static_cast<std::size_t>(locus)];
        const int local = local_position_by_locus[static_cast<std::size_t>(locus)];
        auto& context = *contexts[static_cast<std::size_t>(owner)];
        auto& workspace = *context.workspaces[static_cast<std::size_t>(local)];

        cfg::PredictiveResult scored;
        try {
          scored = workspace.raw_common_predictive(side_contexts[i].raw_boundary,
            side_contexts[i].raw_request_index, true);
        } catch (const std::exception& error) {
          // An unsuccessful optional reference refinement is unavailable
          // evidence, never a genetic zero or a hard-invalid incumbent.
          scored.reason = error.what();
        }
        ConfigurationResult contribution;
        contribution.scored = contribution.raw_only = true;
        contribution.raw_status = scored.status;
        contribution.raw_reference_error = scored.reference_error;
        contribution.raw_message_approximate = scored.message_approximate;
        contribution.raw_comparison_available = scored.status == cfg::PredictiveStatus::DefinedFinite ||
          scored.status == cfg::PredictiveStatus::DefinedGeneticZero;
        contribution.locus_objective_ln = contribution.raw_comparison_available ? scored.ln : NA_REAL;
        contribution.route = value.route;
        reduce_configuration_locus(contribution, replacement);

      }
      value.objective_ln = replacement.objective_ln;
      value.raw_reference_error = replacement.raw_reference_error;
      value.raw_message_approximate = replacement.raw_message_approximate;
      value.raw_accumulation_overflow = replacement.raw_accumulation_overflow;
      value.raw_status = replacement.raw_status;
      value.raw_comparison_available = replacement.raw_comparison_available;
      value.raw_loci_accumulated = replacement.raw_loci_accumulated;
      value.raw_support_rejected = false;

      refined = true;
    }

    // The finalized rows below are the only rows that may enter the memo.

  }

auto VisitScorer::configure_raw_arithmetic(bool fast) -> void  {
    for (auto& context : contexts) context->pass.raw_scaled_fast = fast;
  }

auto VisitScorer::configure_message_scoring(bool selected, bool published) -> void  {
    unavailable_message_trial.store(false,std::memory_order_relaxed);
    for (auto& context : contexts) {
      context->pass.message_scoring = selected;

    }
  }

auto VisitScorer::probe_required_raw_backend(const SearchInputs& input,
      const std::vector<RelationshipContexts>& side_contexts,
      const std::vector<ConfigurationScore>& configurations) -> void  {
    if (!input.message_scoring || !input.stored_message_scoring ||
        input.private_complete_message_trial || owner_by_locus.empty()) return;
    std::size_t row=0;
    while(row<configurations.size() && (configurations[row].raw_numerical_cached ||
      !raw_stage_needed(configurations[row]) || !side_contexts[row].raw_boundary ||
      !side_contexts[row].raw_boundary->ready())) ++row;
    if(row==configurations.size()) return;
    int locus=0;
    for(int i=1;i<static_cast<int>(owner_by_locus.size());++i)
      if(primary().pass.compiled.locus[static_cast<std::size_t>(i)].k >
         primary().pass.compiled.locus[static_cast<std::size_t>(locus)].k) locus=i;
    auto& context=*contexts[static_cast<std::size_t>(owner_by_locus[static_cast<std::size_t>(locus)])];
    auto& workspace=*context.workspaces[static_cast<std::size_t>(local_position_by_locus[static_cast<std::size_t>(locus)])];

    const auto started=std::chrono::steady_clock::now();
    try {
      (void)workspace.raw_common_predictive(side_contexts[row].raw_boundary,
        side_contexts[row].raw_request_index);
    } catch(const std::exception&) {
      // An optional probe failure cannot manufacture a scientific refusal.
      // Retain any genuine recorded partition decline; otherwise the original
      // raw caller rebinds and applies its established exception semantics.

    }

  }

auto VisitScorer::message_trial_aborted() const -> bool  {
    return unavailable_message_trial.load(std::memory_order_relaxed);
  }

auto VisitScorer::refine_relationship_panels(
      const SearchInputs& input, int child,
      const std::vector<RelationshipContexts>& side_contexts,
      std::vector<ConfigurationScore>& configurations) -> void  {
    if (!input.config->evaluation.relationship_common_scaled_fast && !input.config->evaluation.relationship_common_pair_fast) return;
    int null_index = -1;
    for (std::size_t i = 0; i < configurations.size(); ++i)
      if (configurations[i].nparents == 0) {
        null_index = static_cast<int>(i);
        break;
      }
    if (null_index < 0) return;
    const double null_objective = input.elide_raw() ? 0.0 : configurations[null_index].objective_ln;
    bool refined = false;
    for (std::size_t i = 0; i < configurations.size(); ++i) {
      ConfigurationScore& configuration = configurations[i];
      if (!configuration.hard_feasible || configuration.status != 0 ||
          configuration.nparents == 0 || configuration.raw_support_rejected ||
          !relationship_stage_needed(configuration) ||
          configuration.full_relationship_members == nullptr) continue;
      SideTotals trial = *configuration.full_relationship_members;
      double uncertainty = std::max(
        trial.relationship_common_error_bound[0], trial.relationship_common_error_bound[1]);
      finalize_contextual_side_members(&trial);
      align_contextual_operational_display_reference(&trial, null_objective);
      const SideReduction reduction = full_reduce_members(
        &trial, configuration.nparents, null_objective);
      uncertainty = std::max(uncertainty, std::max(
        reduction.side_reference_uncertainty[0], reduction.side_reference_uncertainty[1]));

      const bool no_refinement=!(uncertainty>0.0)||!reduction.comparison_available||!R_finite(reduction.relationship_ln)||std::fabs(reduction.relationship_ln-input.Tassign_ln)>uncertainty;
      if(auto* live=pap_internal::native_demand::live;live&&live->engine_reduction_reuse&&no_refinement){
        *configuration.full_relationship_members=std::move(trial);
        configuration.native_full_reduction=std::make_shared<const SideReduction>(reduction);
      }
      if (!(uncertainty > 0.0)) continue;
      if (!reduction.comparison_available || !R_finite(reduction.relationship_ln) ||
          std::fabs(reduction.relationship_ln - input.Tassign_ln) > uncertainty)
        continue;
      SideTotals replacement;
      ConfigurationScore reference = configuration;
      reference.full_relationship_members = &replacement;
      restore_cached_raw_parent_alias(input, owner_by_locus.size(), configuration, &replacement);
      reference.objective_ln = 0.0;
      for (int locus = 0; locus < static_cast<int>(owner_by_locus.size()); ++locus) {
        const int owner = owner_by_locus[static_cast<std::size_t>(locus)];
        const int local = local_position_by_locus[static_cast<std::size_t>(locus)];
        LocusContext& context = *contexts[static_cast<std::size_t>(owner)];
        const bool saved_scaled = context.pass.options.relationship_common_scaled_fast;
        const bool saved_pair = context.pass.options.relationship_common_pair_fast;
        context.pass.options.relationship_common_scaled_fast = false;
        context.pass.options.relationship_common_pair_fast = false;
        ConfigurationResult result;
        evaluate_configuration_locus(
          input, context.pass, *context.workspaces[static_cast<std::size_t>(local)],
          static_cast<std::size_t>(locus), child, reference, side_contexts[i],
          result);
        context.pass.options.relationship_common_scaled_fast = saved_scaled;
        context.pass.options.relationship_common_pair_fast = saved_pair;
        reduce_configuration_locus(result, reference);

      }
      // Raw objectives and their identity/conditional competitors retain
      // their original ordered values; only the relationship payload changes.
      *configuration.full_relationship_members = std::move(replacement);
      refined = true;
    }

  }

auto VisitScorer::provenance() const -> List  {
    IntegerVector owners(owner_by_locus.size());
    NumericVector weights(weight_by_locus.size());
    for (std::size_t locus = 0; locus < owner_by_locus.size(); ++locus) {
      owners[locus] = owner_by_locus[locus] + 1;
      weights[locus] = weight_by_locus[locus];
    }
    return List::create(
      _["requested_workers"] = requested_workers,
      _["effective_workers"] = effective_workers,
      _["worker_thread_count"] = worker_thread_count,
      _["locus_owner"] = owners,
      _["locus_weight"] = weights,
      _["assignment_policy"] =
        "deterministic_lpt_G2_descending_original_locus_index_lightest_worker",
      _["reduction_policy"] = "serial_original_locus_order",
      _["worker_r_api_calls"] = 0);
  }

}
