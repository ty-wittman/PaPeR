#pragma once
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace pap_internal::pedigree_likelihood::live_arena {
struct Stats {
  std::size_t packed_peak_cells = 0;
  std::size_t concurrent_live_peak_cells = 0;
  std::size_t originally_requested_cells = 0;
  std::size_t free_fragments_at_end = 0;
  bool within_limit = false;
  int first_exceeding_step = -1;
};

// All region sizes and offsets are expressed in double elements, not bytes.
// Critical invariant: the output of a step cannot alias ANY input of that
// same step. Input regions are released only after the step has completed.
// This is a static offset compiler; it never allocates tables at runtime.
template <class Program>
Stats pack_program(Program& p, std::size_t limit_cells) {
  const auto npos = std::numeric_limits<std::size_t>::max();
  const int steps = static_cast<int>(p.steps.size());
  const int slots = static_cast<int>(p.slots.size());
  std::vector<int> last_use(slots, -1);
  std::vector<std::size_t> length(slots, 0), offset(slots, npos);
  std::vector<unsigned char> active(slots, 0);
  for (int i = 0; i < steps; ++i) {
    for (int s : p.steps[i].inputs) {
      if (s < 0 || s >= slots) throw std::logic_error("invalid input slot");
      last_use[s] = i;
    }
  }
  for (int s : p.scalars) {
    if (s < 0 || s >= slots) throw std::logic_error("invalid scalar slot");
    last_use[s] = steps;
  }
  if (p.result >= 0) {
    if (p.result >= slots) throw std::logic_error("invalid result slot");
    last_use[p.result] = steps;
  }

  std::map<std::size_t, std::size_t> free_by_offset;
  std::size_t high_water = 0, live = 0, peak_live = 0;
  const auto add_checked = [](std::size_t a, std::size_t b) -> std::size_t {
    if (b > std::numeric_limits<std::size_t>::max() - a)
      throw std::overflow_error("arena offsets overflow");
    return a + b;
  };
  const auto allocate = [&](std::size_t count) -> std::size_t {
    if (!count) throw std::logic_error("zero length region");
    auto best = free_by_offset.end();
    for (auto it = free_by_offset.begin(); it != free_by_offset.end(); ++it)
      if (it->second >= count && (best == free_by_offset.end() ||
          it->second < best->second ||
          (it->second == best->second && it->first < best->first))) best = it;
    std::size_t at;
    if (best != free_by_offset.end()) {
      at = best->first;
      const auto rem = best->second - count;
      free_by_offset.erase(best);
      if (rem) free_by_offset.emplace(add_checked(at, count), rem);
    } else {
      // An unused interval touching the high-water mark may be extended
      // rather than abandoning that tail and appending the entire request.
      auto tail = free_by_offset.empty() ? free_by_offset.end() :
                  std::prev(free_by_offset.end());
      if (tail != free_by_offset.end() &&
          add_checked(tail->first, tail->second) == high_water) {
        at = tail->first;
        const auto existing = tail->second;
        free_by_offset.erase(tail);
        high_water = add_checked(high_water, count - existing);
      } else {
        at = high_water;
        high_water = add_checked(high_water, count);
      }
    }
    live = add_checked(live, count);
    peak_live = std::max(peak_live, live);
    return at;
  };
  const auto release = [&](std::size_t at, std::size_t count) {
    if (at == npos || !count || live < count) throw std::logic_error("invalid release");
    live -= count;
    auto next = free_by_offset.lower_bound(at);
    if (next != free_by_offset.end() && add_checked(at, count) > next->first)
      throw std::logic_error("overlapping free region");
    if (next != free_by_offset.begin()) {
      auto prev = std::prev(next);
      if (add_checked(prev->first, prev->second) > at)
        throw std::logic_error("overlapping free region");
      if (add_checked(prev->first, prev->second) == at) {
        at = prev->first;
        count = add_checked(count, prev->second);
        free_by_offset.erase(prev);
      }
    }
    if (next != free_by_offset.end() && add_checked(at, count) == next->first) {
      count = add_checked(count, next->second);
      free_by_offset.erase(next);
    }
    free_by_offset.emplace(at, count);
  };

  for (int s = 0; s < slots; ++s) {
    if (p.slots[s].kind != Program::UnaryTable) continue;
    const int v = p.slots[s].index;
    if (v < 0 || static_cast<std::size_t>(v) >= p.unary_offset.size())
      throw std::logic_error("invalid unary variable");
    length[s] = static_cast<std::size_t>(p.G);
    offset[s] = allocate(length[s]);
    active[s] = 1;
    p.slots[s].offset = offset[s];
    p.unary_offset[static_cast<std::size_t>(v)] = offset[s];
  }
  int first_exceeding_step = -1;
  for (int i = 0; i < steps; ++i) {
    auto& step = p.steps[i];
    const int out = step.out;
    if (out < 0 || out >= slots || p.slots[out].kind != Program::Produced || active[out])
      throw std::logic_error("invalid produced slot");
    for (int s : step.inputs) {
      if ((p.slots[s].kind == Program::UnaryTable ||
           p.slots[s].kind == Program::Produced) && !active[s])
        throw std::logic_error("consumed internal factor not live");
    }
    length[out] = step.cells;
    offset[out] = allocate(step.cells);  // allocate BEFORE any input is freed
    active[out] = 1;
    if (high_water > limit_cells && first_exceeding_step < 0) first_exceeding_step = i;
    step.offset = offset[out];
    p.slots[out].offset = offset[out];
    for (int s : step.inputs) {
      if (last_use[s] != i || !active[s]) continue;
      release(offset[s], length[s]);
      active[s] = 0;
    }
    if (last_use[out] < 0) { // a dead, unreferenced produced factor
      release(offset[out], length[out]);
      active[out] = 0;
    }
  }
  Stats result;
  result.originally_requested_cells = p.arena;
  result.packed_peak_cells = high_water;
  result.concurrent_live_peak_cells = peak_live;
  result.free_fragments_at_end = free_by_offset.size();
  result.within_limit = high_water <= limit_cells;
  result.first_exceeding_step = first_exceeding_step;
  p.arena = high_water;
  return result;
}
} // namespace pap_internal::pedigree_likelihood::live_arena
