// Native functions that no supported input reaches have their bodies replaced
// by removed_path(). Reaching one records its name, which the R side reports
// as an internal error after the native call, even when a caller catches the
// exception.
#pragma once
#include <atomic>
#include <stdexcept>
#include <string>

namespace pap_internal {
inline std::atomic<const char*>& removed_path_hit() {
  static std::atomic<const char*> hit{nullptr};
  return hit;
}
[[noreturn]] inline void removed_path(const char* name) {
  const char* none = nullptr;
  removed_path_hit().compare_exchange_strong(none, name);
  throw std::logic_error(std::string("PaPeR: removed code path: ") + name);
}
}  // namespace pap_internal
