#ifndef PAP_PARTITION_TOTALS_H
#define PAP_PARTITION_TOTALS_H
#include "input_checks.h"
#include "score_types.h"

// Full-mode result aggregation for the live pedigree likelihood scorer. Topology catalog
// metadata and partition result values live in topology_types.h.

namespace pap_internal {
namespace side_likelihood {

inline bool valid_partition_value(const PartitionValue& value) {
  if (value.state == PartitionState::positive)
    return std::isfinite(value.log_partition) &&
      std::isfinite(value.reference_error) && value.reference_error >= 0.0;
  if (value.state == PartitionState::zero)
    return value.log_partition == -std::numeric_limits<double>::infinity() &&
      value.reference_error == 0.0;
  return value.state == PartitionState::inapplicable ||
         value.state == PartitionState::unavailable;
}

// Serial locus-order accumulation. A scored zero is absorbing. Unavailable
// is absorbing too, but is never manufactured into zero or inapplicability.
inline bool append_panel_cell(PanelCell& cell, const PartitionValue& value) {
  const bool first = cell.loci == 0U;
  ++cell.loci;
  if (!valid_partition_value(value)) {
    cell.state = PartitionState::unavailable;
    cell.log_partition = std::numeric_limits<double>::quiet_NaN();
    cell.reason = "Invalid raw partition or reference error";
    return false;
  }
  if (first) {
    cell.state = value.state;
    cell.log_partition = value.state == PartitionState::inapplicable ?
      -std::numeric_limits<double>::infinity() : value.log_partition;
    cell.reference_error = value.reference_error;
    cell.reason = value.reason;
    return cell.state != PartitionState::unavailable;
  }
  if (cell.state == PartitionState::unavailable) return false;
  if (value.state == PartitionState::unavailable ||
      ((cell.state == PartitionState::inapplicable) !=
       (value.state == PartitionState::inapplicable))) {
    cell.state = PartitionState::unavailable;
    cell.log_partition = std::numeric_limits<double>::quiet_NaN();
    cell.reason = value.state == PartitionState::unavailable ? value.reason :
      "A fixed topology changed applicability across loci";
    return false;
  }
  if (cell.state == PartitionState::inapplicable) return true;
  if (cell.state == PartitionState::zero || value.state == PartitionState::zero) {
    cell.state = PartitionState::zero;
    cell.log_partition = -std::numeric_limits<double>::infinity();
    cell.reference_error = 0.0;
    return true;
  }
  const double previous = cell.log_partition;
  cell.log_partition += value.log_partition;
  // The panel reference is this exact same serial order applied to the own
  // retained per-locus outputs. If every input equals that reference, the
  // accumulated output does too; do not create a spurious refinement loop.
  if (cell.reference_error != 0.0 || value.reference_error != 0.0) {
    const double rounding = 4.0 * std::numeric_limits<double>::epsilon() *
      (std::fabs(previous) + std::fabs(value.log_partition) + std::fabs(cell.log_partition));
    cell.reference_error = std::nextafter(
      cell.reference_error + value.reference_error + rounding,
      std::numeric_limits<double>::infinity());
  }
  if (!std::isfinite(cell.log_partition) || !std::isfinite(cell.reference_error)) {
    cell.state = PartitionState::unavailable;
    cell.log_partition = std::numeric_limits<double>::quiet_NaN();
    cell.reason = "Panel partition or reference-error accumulation overflow";
    return false;
  }
  return true;
}

inline bool append_full_locus(FullPanelTotals& panel,
                             const std::vector<PartitionValue>& values) {
  if(input_checks::enabled()){

  if (!panel.boundary || values.size() != panel.cells.size() ||
      panel.loci >= panel.expected_loci) {

    panel.failure = "Full locus contribution does not match its fixed request or panel length";
    return false;
  }
  }
  ++panel.loci;
  bool available = true;
  for (std::size_t i = 0U; i < values.size(); ++i) {
    const auto& model = panel.boundary->compiled.models[i];
    PartitionValue value = values[i];
    if (!panel.boundary->compiled.native_support && (!model.supported || !model.topology.ready() ||
        value.state == PartitionState::inapplicable)) {
      value.state = PartitionState::unavailable;
      value.reason = !model.supported ? model.reason :
        "A compiled biologically Ready full model was not scored";
    }
    available = append_panel_cell(panel.cells[i], value) && available;
  }
  return available;
}

enum class ContrastState { defined, undefined, unavailable };
struct AbsoluteSideReduction {
  ContrastState state = ContrastState::unavailable;
  bool parent_complete = false;
  bool alternative_complete = false;
  double parent_ln = std::numeric_limits<double>::quiet_NaN();
  double alternative_ln = std::numeric_limits<double>::quiet_NaN();
  double own_u_ln = std::numeric_limits<double>::quiet_NaN();
  bool u_relative_available = false;
  double relationship_ln = std::numeric_limits<double>::quiet_NaN();
  double relationship_lower_ln = std::numeric_limits<double>::quiet_NaN();
  double relationship_upper_ln = std::numeric_limits<double>::quiet_NaN();
  // Ordinary ties use their existing 0..24 model IDs. Full ties use
  // 25 + descriptor index, retaining every literal direction/slot variant.
  std::vector<std::size_t> parent_ties, alternative_ties;
  std::string reason;
};

namespace full_scoring_detail {
struct Maximum {
  double value = -std::numeric_limits<double>::infinity();
  double lower = -std::numeric_limits<double>::infinity();
  double upper = -std::numeric_limits<double>::infinity();
  bool seen = false;
  std::vector<std::size_t> ties;
  void add(const PanelCell& cell, std::size_t id) {
    if (cell.state == PartitionState::inapplicable) return;
    seen = true;
    if (cell.log_partition > value) {
      value = cell.log_partition;
      ties.clear();
    }
    if (cell.log_partition == value) ties.push_back(id);
    if (cell.state == PartitionState::zero) return;
    const double raw_lower = cell.log_partition - cell.reference_error;
    const double raw_upper = cell.log_partition + cell.reference_error;
    lower = std::max(lower, cell.reference_error == 0.0 ? raw_lower :
      std::nextafter(raw_lower, -std::numeric_limits<double>::infinity()));
    upper = std::max(upper, cell.reference_error == 0.0 ? raw_upper :
      std::nextafter(raw_upper, std::numeric_limits<double>::infinity()));
  }
};
inline ::relationship_topology::TopologyStatus applicability(const PanelCell& cell,
    ::relationship_topology::TopologyStatus original, bool native) {
  return native ? (cell.state == PartitionState::inapplicable ?
    ::relationship_topology::TopologyStatus::Inapplicable : ::relationship_topology::TopologyStatus::Ready) : original;
}
inline bool complete_cell(const PanelCell& cell, std::size_t loci,
                          ::relationship_topology::TopologyStatus applicability) {
  if (cell.loci != loci || loci == 0U) return false;
  if (applicability == ::relationship_topology::TopologyStatus::Inapplicable)
    return cell.state == PartitionState::inapplicable;
  if (applicability != ::relationship_topology::TopologyStatus::Ready) return false;
  return (cell.state == PartitionState::positive || cell.state == PartitionState::zero) &&
    valid_partition_value({cell.state, cell.log_partition, cell.reference_error,
                           false, cell.reason});
}
}  // namespace full_scoring_detail

// Maxima are taken only after every fixed literal topology has accumulated
// the whole panel. A side's raw partitions cannot be borrowed by another
// side, and no marginal alias/frozen-context member is admitted here.
inline AbsoluteSideReduction reduce_full_side(const Boundary& ordinary,
    const std::array<PanelCell, 25>& relationship_panel, const FullPanelTotals& full) {
  AbsoluteSideReduction out;
  if (!full.boundary || ordinary.full_boundary.get() != full.boundary ||
      !ordinary.failure.empty() || !full.boundary->failure.empty() ||
      !full.failure.empty() || ordinary.models.size() != 25U ||
      ordinary.evidence_ids != full.boundary->compiled.evidence_ids ||
      full.loci != full.expected_loci || full.expected_loci == 0U ||
      full.cells.size() != full.boundary->compiled.models.size()) {
    out.reason = "Full/ordinary panel is incomplete or not owned by one common side boundary";
    return out;
  }
  if (ordinary.models[0].topology.relationship != "PO" ||
      ordinary.models[1].topology.relationship != "U") {
    out.reason = "Full reduction lost the authoritative ordinary PO/U model order";
    return out;
  }
  full_scoring_detail::Maximum parent, alternative;
  for (std::size_t i = 0U; i < relationship_panel.size(); ++i) {
    const auto& model = ordinary.models[i];
    if (!full_scoring_detail::complete_cell(relationship_panel[i], full.expected_loci,
                                          full_scoring_detail::applicability(relationship_panel[i],model.topology.status,ordinary.native_support)) ||
        (!ordinary.native_support && model.topology.ready() && !model.supported)) {
      out.reason = "An applicable ordinary model is unscored on this full panel";
      return out;
    }
    (i == 0U ? parent : alternative).add(relationship_panel[i], i);
  }
  for (std::size_t i = 0U; i < full.boundary->members.size(); ++i) {
    const auto& descriptor = full.boundary->members[i];
    if (descriptor.topology_status == ::relationship_topology::TopologyStatus::Inapplicable) continue;
    if (descriptor.topology_status != ::relationship_topology::TopologyStatus::Ready ||
        descriptor.compiled_index >= full.cells.size()) {
      out.reason = "An applicable full hypothesis lacks a supported compiled model";
      return out;
    }
    const auto& model = full.boundary->compiled.models[descriptor.compiled_index];
    const auto& cell = full.cells[descriptor.compiled_index];
    if ((!ordinary.native_support && !model.supported) || !full_scoring_detail::complete_cell(
        cell, full.expected_loci, full_scoring_detail::applicability(cell,descriptor.topology_status,ordinary.native_support))) {
      out.reason = "An applicable full model is unscored on this panel";
      return out;
    }
    (descriptor.parent_numerator ? parent : alternative).add(cell, 25U + i);
  }
  out.parent_complete = true;
  out.alternative_complete = true;
  out.parent_ln = parent.value;
  out.alternative_ln = alternative.value;
  out.own_u_ln = relationship_panel[1].log_partition;
  out.u_relative_available = relationship_panel[1].state == PartitionState::positive;
  out.parent_ties = std::move(parent.ties);
  out.alternative_ties = std::move(alternative.ties);
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  // Zero U cannot define a U-relative display/reference. Positive parent
  // support alongside zero U violates the common-graph support contract;
  // do not convert that construction failure into infinite parent support.
  if (out.parent_ln != negative_infinity && out.own_u_ln == negative_infinity) {
    out.reason = "Positive full parent hypothesis with zero common U support";
    return out;
  }
  if (out.parent_ln == negative_infinity && out.alternative_ln == negative_infinity) {
    out.state = ContrastState::undefined;
    out.reason = "All competing parent and nonparent hypotheses have scored zero likelihood";
    return out;
  }
  out.state = ContrastState::defined;
  out.relationship_ln = out.parent_ln - out.alternative_ln;
  out.relationship_lower_ln = parent.lower - alternative.upper;
  out.relationship_upper_ln = parent.upper - alternative.lower;
  return out;
}


}  // namespace side_likelihood
}  // namespace pap_internal

#endif
