#ifndef PAP_INITIAL_PAIR_CHECK_H
#define PAP_INITIAL_PAIR_CHECK_H

// Inclusive support check for one supplied pair. These are four COMPLETE
// marker panels on one frozen bounded predictive background: null,D,S,pair.
// This helper neither scores hypotheses nor proposes assignments or moves.
#include "support_rejection.h"
#include <array>

namespace pap_internal {
namespace initial_pair_check {

using Panel = support_rejection::Panel;

struct Result {
  double support_ln = std::numeric_limits<double>::quiet_NaN();
  double rejection_upper_ln = std::numeric_limits<double>::quiet_NaN();
  bool support_available = false;
  bool support_pass = false;
  bool proven_rejection = false;
};

inline Result evaluate(const std::array<Panel, 4>& Q, double threshold,
                       bool hard_rejected = false) {
  namespace raw = support_rejection;
  Result result;
  result.proven_rejection = hard_rejected;
  if (!std::isfinite(threshold) &&
      threshold != -std::numeric_limits<double>::infinity()) return result;
  bool all_defined = true;
  for (const Panel& panel : Q) all_defined = all_defined && raw::defined(panel);
  if (all_defined) {
    const Panel baseline{std::max(Q[0].ln, std::max(Q[1].ln, Q[2].ln)), true, true};
    result.support_ln = raw::contrast(Q[3], baseline);
    result.support_available = !std::isnan(result.support_ln);
    result.support_pass = result.support_available && result.support_ln >= threshold;
  }
  // Any complete baseline is a lower bound on max(null,D,S), hence a
  // corresponding upper bound on pair support. Missing rivals can only make
  // support smaller; they cannot invalidate an already proven rejection.
  // Equality is NOT a rejection because CheckPedigree preserves >=.
  for (std::size_t i = 0U; i != 3U; ++i) {
    const double upper = raw::contrast(Q[3], Q[i]);
    if (!std::isnan(upper) && (std::isnan(result.rejection_upper_ln) ||
        upper < result.rejection_upper_ln)) result.rejection_upper_ln = upper;
  }
  result.proven_rejection = result.proven_rejection ||
    (!std::isnan(result.rejection_upper_ln) && result.rejection_upper_ln < threshold);
  return result;
}

} // namespace initial_pair_check
} // namespace pap_internal
#endif
