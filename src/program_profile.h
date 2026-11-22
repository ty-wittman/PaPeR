#pragma once
#include <cstdint>

namespace pap_internal::pedigree_likelihood {
// Stable identifier used by the optional exact-resource execution receipt.
inline std::uint64_t program_profile_identity(const Program& p) {
  std::uint64_t h = 1469598103934665603ULL;
  auto add = [&](std::uint64_t x) { h = (h ^ x) * 1099511628211ULL; };
  add(p.G); add(p.k); add(p.model_variables); add(p.model_births);
  add(p.model_pieces); add(p.retained_axes); add(p.arena); add(p.result);
  for (const auto& slot : p.slots) { add(slot.kind); add(slot.index); add(slot.offset); }
  for (const auto& step : p.steps) {
    add(step.variable); add(step.cells); add(step.dim); add(step.select); add(step.select_role);
    for (int x : step.inputs) add(x);
    for (int x : step.dims) add(x);
    for (auto x : step.axis0) add(x);
    for (auto x : step.strides) add(x);
  }
  return h;
}
}
