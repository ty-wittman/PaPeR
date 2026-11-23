#ifndef PAP_PANEL_STAGING_H
#define PAP_PANEL_STAGING_H

// Execution partitioning of the 25-member relationship catalog. Every
// model keeps its own fixed marker sum; no hypothesis, likelihood or gate
// is defined here. Compact prefix values are exact saved scorer outputs.
#include "partition_totals.h"

namespace pap_internal {
namespace side_staging {
namespace b = side_likelihood;

enum class Phase : unsigned char { complete, prefix, remainder };
inline bool selected(std::size_t model, Phase phase) {
  return model < 25U && (phase == Phase::complete ||
    (phase == Phase::prefix ? model < 3U : model >= 3U));
}
inline bool valid_catalog(const b::Boundary& boundary) {
  return boundary.failure.empty() && boundary.models.size() == 25U &&
    boundary.models[0].topology.relationship == "PO" &&
    boundary.models[1].topology.relationship == "U" &&
    boundary.models[2].topology.relationship == "FS";
}

// No factors, messages, strings or per-model graphs are retained. Unsupported
// prefix operands use the original full evaluation on completion instead of
// caching their refusal provenance. This cannot make them a scored zero.
struct Prefix {
  const b::Boundary* owner = nullptr;
  int marker = -1;
  std::array<double, 3> partition{};
  std::array<double, 3> error{};
  std::array<unsigned char, 3> status{{3U, 3U, 3U}};
  bool supported = false;
  bool present = false;

  bool reusable(const b::Boundary* requested_owner, int requested_marker) const {
    if (!present || owner != requested_owner || marker != requested_marker ||
        !requested_owner || !valid_catalog(*requested_owner)) return false;
    for (std::size_t model = 0U; model != 3U; ++model) {
      if (status[model] > 2U || !std::isfinite(error[model]) || error[model] < 0.0)
        return false;
      if (status[model] == 1U && !std::isfinite(partition[model])) return false;
      if (status[model] != 1U &&
          !(std::isinf(partition[model]) && partition[model] < 0.0)) return false;
      if (status[model] == 2U && error[model] != 0.0) return false;
    }
    return true;
  }
};

inline Prefix capture(const b::Boundary* owner, int marker,
                      const std::vector<b::PartitionValue>& values,
                      bool supported) {
  Prefix out;
  if (!owner || values.size() != 25U || !valid_catalog(*owner)) return out;
  out.owner = owner; out.marker = marker; out.supported = supported;
  out.present = true;
  for (std::size_t model = 0U; model != 3U; ++model) {
    out.partition[model] = values[model].log_partition;
    out.error[model] = values[model].reference_error;
    out.status[model] = static_cast<unsigned char>(values[model].state);
  }
  return out;
}

// Validate the stage/count protocol before touching any scalar. Complete
// and remainder preserve the original per-model append_panel_cell arithmetic.
// A legitimate unavailable value is appended with its original provenance;
// the return value describes count/owner layout, not likelihood availability.
inline bool append_locus(std::array<b::PanelCell, 25>& panel,
                         const std::vector<b::PartitionValue>& values,
                         Phase phase, std::size_t expected_loci) {
  if (values.size() != panel.size() || expected_loci == 0U) return false;
  std::size_t previous = 0U;
  bool seen = false;
  for (std::size_t model = 0U; model != panel.size(); ++model) {
    if (phase == Phase::remainder && model < 3U &&
        panel[model].loci != expected_loci) return false;
    if (!selected(model, phase)) continue;
    if (panel[model].loci >= expected_loci ||
        (seen && panel[model].loci != previous)) return false;
    previous = panel[model].loci; seen = true;
  }
  for (std::size_t model = 0U; model != panel.size(); ++model)
    if (selected(model, phase)) b::append_panel_cell(panel[model], values[model]);
  return true;
}

}  // namespace side_staging
}  // namespace pap_internal


// Work partitioning only. The caller scores the SAME immutable literal models
// and keeps marker order. This helper calculates no likelihood or threshold.
#include "partition_totals.h"
#include "input_checks.h"

namespace pap_internal {
namespace panel_staging {
namespace b = side_likelihood;

enum class Phase : unsigned char { complete, parents, alternatives };

inline bool valid_map(const b::FullBoundary& boundary) {
  if(!input_checks::enabled())return true;

  if (boundary.compiled_to_member.size() != boundary.compiled.models.size()){return false;}
  for (std::size_t model = 0U; model < boundary.compiled.models.size(); ++model) {
    const auto member = boundary.compiled_to_member[model];
    if (member >= boundary.members.size() || boundary.members[member].compiled_index != model)
      {return false;}
  }
  return true;
}

inline bool selected_unchecked(const b::FullBoundary& boundary, std::size_t model, Phase phase) {
  return phase==Phase::complete||boundary.members[boundary.compiled_to_member[model]].parent_numerator==(phase==Phase::parents);
}
inline bool selected(const b::FullBoundary& boundary, std::size_t model, Phase phase) {
  if(!input_checks::enabled())return selected_unchecked(boundary,model,phase);

  if (phase == Phase::complete) return model < boundary.compiled.models.size();
  if (model >= boundary.compiled_to_member.size()) return false;
  const auto member = boundary.compiled_to_member[model];
  if (member >= boundary.members.size() || boundary.members[member].compiled_index != model) return false;
  return boundary.members[member].parent_numerator == (phase == Phase::parents);
}

// The original append contract covers per-model unavailable/zero arithmetic.
// Parents do not increment the complete-catalog marker count. Alternatives
// finish one complete marker only after every parent already owns ALL markers.
inline bool append_locus(b::FullPanelTotals& panel,
                         const std::vector<b::PartitionValue>& values,
                         Phase phase) {
  if (phase == Phase::complete) return b::append_full_locus(panel, values);
  const auto refuse = [&](const char* reason) { panel.failure = reason; return false; };
  if(input_checks::enabled()){
  const bool map_ok=panel.boundary&&valid_map(*panel.boundary);

  if (!panel.boundary || values.size() != panel.cells.size() ||
      panel.cells.size() != panel.boundary->compiled.models.size() ||
      !map_ok || panel.expected_loci == 0U)
    return refuse("Staged full contribution is not owned by its immutable catalog");
  if (phase == Phase::parents) {
    if (panel.loci != 0U)
      return refuse("Staged parent panel follows an already completed catalog");
    std::size_t first_count = panel.expected_loci;
    bool seen = false;
    for (std::size_t model = 0U; model < values.size(); ++model) {
      if (!selected_unchecked(*panel.boundary, model, phase)) continue;
      if (panel.cells[model].loci >= panel.expected_loci)
        return refuse("Staged parent model received a duplicate marker");
      if (seen && panel.cells[model].loci != first_count)
        return refuse("Staged parent models have inconsistent ordered panels");
      first_count = panel.cells[model].loci;
      seen = true;
    }
  } else {
    if (panel.loci >= panel.expected_loci)
      return refuse("Staged full alternatives received a duplicate marker");
    for (std::size_t model = 0U; model < values.size(); ++model) {
      const bool parent = selected_unchecked(*panel.boundary, model, Phase::parents);
      if (parent && panel.cells[model].loci != panel.expected_loci)
        return refuse("Staged full alternatives preceded a complete parent panel");
      if (!parent && panel.cells[model].loci != panel.loci)
        return refuse("Staged alternatives have inconsistent ordered panels");
    }
  }
  }
  if(phase==Phase::alternatives)++panel.loci;
  bool available = true;
  for (std::size_t model_index = 0U; model_index < values.size(); ++model_index) {
    if (!selected(*panel.boundary, model_index, phase)) continue;
    const auto& model = panel.boundary->compiled.models[model_index];
    auto value = values[model_index];
    // Identical scientific availability normalization to append_full_locus.
    if (!panel.boundary->compiled.native_support && (!model.supported || !model.topology.ready() || value.state == b::PartitionState::inapplicable)) {
      value.state = b::PartitionState::unavailable;
      value.reason = !model.supported ? model.reason :
        "A compiled biologically Ready full model was not scored";
    }
    available = b::append_panel_cell(panel.cells[model_index], value) && available;
  }
  return available;
}

}  // namespace panel_staging
}  // namespace pap_internal
#endif
