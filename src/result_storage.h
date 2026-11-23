#ifndef PAP_RESULT_STORAGE_H
#define PAP_RESULT_STORAGE_H

// Result ownership/scheduling only. It never changes numerical arithmetic,
// candidate availability, per-locus ownership, or scientific reductions.
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace pap_internal {
namespace result_storage {
inline std::size_t add(std::size_t a, std::size_t b) {
  return b > std::numeric_limits<std::size_t>::max() - a ?
    std::numeric_limits<std::size_t>::max() : a + b;
}
inline std::size_t multiply(std::size_t a, std::size_t b) {
  return a != 0U && b > std::numeric_limits<std::size_t>::max() / a ?
    std::numeric_limits<std::size_t>::max() : a * b;
}
inline std::size_t heap_string_bytes(const std::string& value) {
  // A small-string buffer lies inside its owner. Its bytes are already in
  // sizeof(PartitionValue); count only an external capacity allocation.
  const auto data = reinterpret_cast<std::uintptr_t>(value.data());
  const auto begin = reinterpret_cast<std::uintptr_t>(&value);
  const auto end = begin + sizeof(value);
  return data >= begin && data < end ? 0U : add(value.capacity(), 1U);
}
struct Plan {
  std::size_t width = 1U, metadata = 0U, per_locus = 0U;
  std::size_t estimated_wave = 0U, estimated_all_loci = 0U;
  bool irreducible_budget_exceeded = false;
};
inline Plan plan(std::size_t loci, std::size_t per_locus,
                 std::size_t metadata, std::size_t budget) {
  Plan out;
  out.metadata = metadata; out.per_locus = per_locus;
  out.estimated_all_loci = add(metadata, multiply(loci, per_locus));
  if (loci == 0U) { out.width = 0U; out.estimated_wave = metadata; return out; }
  const std::size_t remaining = budget > metadata ? budget - metadata : 0U;
  out.width = budget == 0U || per_locus == 0U ? loci :
    std::min(loci, std::max(std::size_t(1U), remaining / per_locus));
  out.estimated_wave = add(metadata, multiply(out.width, per_locus));
  out.irreducible_budget_exceeded = budget != 0U &&
    add(metadata, per_locus) > budget;
  return out;
}
} // namespace result_storage
} // namespace pap_internal
#endif
