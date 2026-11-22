#ifndef PAP_PANEL_REJECTION_H
#define PAP_PANEL_REJECTION_H

// Whole-panel negative proof for one already declared full-mode side. This
// helper authorizes no addition, replacement, retention, or removal. It proves
// only that the current strict relationship gate must reject, even if some
// NONPARENT likelihoods have not been evaluated. Every applicable direct or
// compound parent must be covered on its own same-D/C complete marker panel.
// It never treats an unknown parent as ordinary PO, takes per-marker maxima,
// ignores a parental orientation, or calls an unavailable rival rejected.
#include "score_rejection.h"
#include "partition_totals.h"

namespace pap_internal {
namespace panel_rejection {

using Outcome = score_rejection::Outcome;
enum class Reason : unsigned char {
  none, different_owner, invalid_threshold, incomplete_parent_maximum,
  parent_edge_ownership, invalid_numerical_bound, unavailable_u_integrity,
  zero_u_integrity, all_zero, rejection_witness, threshold_overlap
};

struct OwnedPanels {
  const side_likelihood::Boundary* ordinary = nullptr;
  const std::array<side_likelihood::PanelCell, 25>* relationship_cells = nullptr;
  const side_likelihood::FullPanelTotals* full = nullptr;
  // Exact common-boundary/panel/revision ownership supplied at scoring time.
  // This tag is not a score cache key; the caller owns immutable inputs.
  score_rejection::Boundary authority;
};

struct Request {
  score_rejection::Boundary authority;
  std::size_t expected_loci = 0U;
  double Tassign_ln = std::numeric_limits<double>::quiet_NaN();
  // True only when every parent and the selected rejection witness are
  // their retained ordered outputs. Other alternatives are unnecessary to
  // the proof, irrespective of whether they have been scored.
  bool retained_reference_arithmetic = false;
  double additional_error_bound_ln = 0.0;
};

struct Certificate {
  Outcome outcome = Outcome::continue_complete;
  Reason reason = Reason::none;
  bool parent_maximum_complete = false;
  bool alternatives_complete = false;
  std::size_t applicable_parents = 0U, available_parents = 0U;
  std::size_t applicable_alternatives = 0U, available_alternatives = 0U;
  std::size_t unscored_alternatives = 0U, unavailable_alternatives = 0U;
  std::size_t witness_id = std::numeric_limits<std::size_t>::max();
  double parent_max_ln = std::numeric_limits<double>::quiet_NaN();
  double parent_upper_ln = std::numeric_limits<double>::quiet_NaN();
  double witness_ln = std::numeric_limits<double>::quiet_NaN();
  double witness_lower_ln = std::numeric_limits<double>::quiet_NaN();
  // Bound on the COMPLETE parent-minus-nonparent maximum contrast. This is
  // not an exact LLR when alternatives_complete is false. Never display it
  // as a likelihood evaluation of skipped hypotheses.
  double upper_ln = std::numeric_limits<double>::quiet_NaN();
};

namespace detail {
inline bool parent_edge(const ::relationship_topology::TopologyModel& model, int side) {
  if (!model.ready() || (side != 1 && side != 2) || model.focal < 0 || model.candidate < 0 ||
      static_cast<std::size_t>(model.focal) >= model.nodes.size() ||
      static_cast<std::size_t>(model.candidate) >= model.nodes.size()) return false;
  const auto& child = model.nodes[static_cast<std::size_t>(model.focal)];
  return child.has_birth && child.parents[static_cast<std::size_t>(side - 1)] == model.candidate;
}
inline bool ready_cell(const side_likelihood::PanelCell& cell,
                       std::size_t loci) {
  return side_likelihood::full_scoring_detail::complete_cell(
    cell, loci, ::relationship_topology::TopologyStatus::Ready);
}
inline double upper(const side_likelihood::PanelCell& cell) {
  if (cell.state == side_likelihood::PartitionState::zero) return -std::numeric_limits<double>::infinity();
  const double value = cell.log_partition + cell.reference_error;
  return cell.reference_error == 0.0 ? value :
    std::nextafter(value, std::numeric_limits<double>::infinity());
}
inline double lower(const side_likelihood::PanelCell& cell) {
  if (cell.state == side_likelihood::PartitionState::zero) return -std::numeric_limits<double>::infinity();
  const double value = cell.log_partition - cell.reference_error;
  return cell.reference_error == 0.0 ? value :
    std::nextafter(value, -std::numeric_limits<double>::infinity());
}
}  // namespace detail

inline Certificate certify(const Request& request, const OwnedPanels& panels) {
  using namespace side_likelihood;
  Certificate out;
  const auto decline = [&](Reason reason) { out.reason = reason; return out; };
  if (!panels.ordinary || !panels.relationship_cells || !panels.full ||
      !request.authority.same(panels.authority) ||
      request.authority.request != panels.ordinary || request.expected_loci == 0U)
    return decline(Reason::different_owner);
  const auto& ordinary = *panels.ordinary;
  const auto& relationship_cells = *panels.relationship_cells;
  const auto& full = *panels.full;
  if (!full.boundary || ordinary.full_boundary.get() != full.boundary ||
      !ordinary.failure.empty() || !full.boundary->failure.empty() ||
      !full.failure.empty() || ordinary.models.size() != 25U ||
      ordinary.models[0].topology.relationship != "PO" || ordinary.models[1].topology.relationship != "U" ||
      ordinary.evidence_ids != full.boundary->compiled.evidence_ids ||
      full.expected_loci != request.expected_loci ||
      full.cells.size() != full.boundary->compiled.models.size() ||
      full.boundary->tested_side != request.authority.side)
    return decline(Reason::different_owner);
  if (!std::isfinite(request.Tassign_ln)) return decline(Reason::invalid_threshold);
  if (!std::isfinite(request.additional_error_bound_ln) || request.additional_error_bound_ln < 0.0)
    return decline(Reason::invalid_numerical_bound);

  double parent_value = -std::numeric_limits<double>::infinity();
  double parent_upper = parent_value;
  double best_alternative_lower = parent_value;
  double best_alternative_value = parent_value;
  bool any_error = request.additional_error_bound_ln != 0.0;
  bool parental_failure = false;
  Reason parental_reason = Reason::incomplete_parent_maximum;
  const auto consume = [&](const Model* model, ::relationship_topology::TopologyStatus applicability,
      const PanelCell* cell, bool is_parent, std::size_t id) {
    if (cell) applicability = full_scoring_detail::applicability(*cell,applicability,ordinary.native_support);
    if (applicability == ::relationship_topology::TopologyStatus::Inapplicable) return;
    if (is_parent) {
      ++out.applicable_parents;
      if (applicability != ::relationship_topology::TopologyStatus::Ready || !model || (!ordinary.native_support && !model->supported) ||
          !cell || !detail::ready_cell(*cell, request.expected_loci)) {
        parental_failure = true;
        return;
      }
      if (!detail::parent_edge(model->topology, request.authority.side)) {
        parental_failure = true;
        parental_reason = Reason::parent_edge_ownership;
        return;
      }
      ++out.available_parents;
      any_error = any_error || cell->reference_error != 0.0;
      parent_value = std::max(parent_value, cell->log_partition);
      parent_upper = std::max(parent_upper, detail::upper(*cell));
    } else {
      ++out.applicable_alternatives;
      if (applicability != ::relationship_topology::TopologyStatus::Ready || !model || (!ordinary.native_support && !model->supported) ||
          !cell || !detail::ready_cell(*cell, request.expected_loci)) {
        if (cell && cell->loci == 0U) ++out.unscored_alternatives;
        else ++out.unavailable_alternatives;
        return;
      }
      ++out.available_alternatives;
      const double lower = detail::lower(*cell);
      if (lower > best_alternative_lower) {
        best_alternative_lower = lower;
        best_alternative_value = cell->log_partition;
        out.witness_id = id;
      }
    }
  };
  for (std::size_t i = 0U; i < ordinary.models.size(); ++i)
    consume(&ordinary.models[i], ordinary.models[i].topology.status, &relationship_cells[i], i == 0U, i);
  for (std::size_t i = 0U; i < full.boundary->members.size(); ++i) {
    const auto& descriptor = full.boundary->members[i];
    const auto index = descriptor.compiled_index;
    const Model* model = index < full.boundary->compiled.models.size() ?
      &full.boundary->compiled.models[index] : nullptr;
    const PanelCell* cell = index < full.cells.size() ? &full.cells[index] : nullptr;
    consume(model, descriptor.topology_status, cell, descriptor.parent_numerator, 25U + i);
  }
  out.parent_maximum_complete = !parental_failure && out.applicable_parents != 0U &&
    out.available_parents == out.applicable_parents;
  out.alternatives_complete = out.available_alternatives == out.applicable_alternatives;
  if (!out.parent_maximum_complete) return decline(parental_reason);
  out.parent_max_ln = parent_value;
  out.parent_upper_ln = parent_upper;
  out.witness_ln = best_alternative_value;
  out.witness_lower_ln = best_alternative_lower;
  const double zero = -std::numeric_limits<double>::infinity();

  if (parent_value != zero) {
    const auto& u = relationship_cells[1];
    // Retain the current compound support-integrity guard. A missing U does
    // not establish that positive parent support is compatible with its
    // declared common-graph boundary. Zero parents need no U-relative ratio.
    if (!detail::ready_cell(u, request.expected_loci)) return decline(Reason::unavailable_u_integrity);
    if (u.state == PartitionState::zero) return decline(Reason::zero_u_integrity);
  }
  if (parent_value == zero && best_alternative_lower == zero) {
    if (out.alternatives_complete) {
      out.outcome = Outcome::undefined;
      out.reason = Reason::all_zero;
    }
    return out;
  }
  if (best_alternative_lower == zero) return out;
  if (parent_value == zero) {
    out.outcome = Outcome::rejected;
    out.reason = Reason::rejection_witness;
    out.upper_ln = zero;
    return out;
  }
  const auto& witness = out.witness_id < 25U ? relationship_cells[out.witness_id] :
    full.cells[full.boundary->members[out.witness_id - 25U].compiled_index];
  any_error = any_error || witness.reference_error != 0.0;
  if (request.retained_reference_arithmetic && any_error)
    return decline(Reason::invalid_numerical_bound);
  const double upper = parent_upper - best_alternative_lower;
  if (!std::isfinite(upper)) return decline(Reason::invalid_numerical_bound);
  out.upper_ln = request.retained_reference_arithmetic ? upper :
    std::nextafter(upper + request.additional_error_bound_ln,
                   std::numeric_limits<double>::infinity());
  if (out.upper_ln <= request.Tassign_ln) {
    out.outcome = Outcome::rejected;
    out.reason = Reason::rejection_witness;
  } else if (parent_value - best_alternative_value <= request.Tassign_ln) {
    out.outcome = Outcome::needs_reference;
    out.reason = Reason::threshold_overlap;
  }
  return out;
}

}  // namespace panel_rejection
}  // namespace pap_internal
#endif
