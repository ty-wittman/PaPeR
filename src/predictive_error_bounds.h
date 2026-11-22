#ifndef PAP_PREDICTIVE_ERROR_BOUNDS_H
#define PAP_PREDICTIVE_ERROR_BOUNDS_H

// Numerical refinement only. Every bound is relative to the retained ordered
// scorer on the SAME declared model, panel and primitives. This helper changes
// neither an evidence boundary nor a production confidence threshold. It can
// request reference Q values; it never authorizes or rejects a parent.
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace pap_internal {
namespace predictive_error_bounds {

namespace detail {
inline double upward(long double value) {
  if (!(value >= 0.0L) || !std::isfinite(value) ||
      value > std::numeric_limits<double>::max())
    return std::numeric_limits<double>::infinity();
  if (value == 0.0L) return 0.0;
  return std::nextafter(static_cast<double>(value),
                        std::numeric_limits<double>::infinity());
}
inline bool valid_error(double value) {
  return std::isfinite(value) && value >= 0.0;
}
inline double operation_error(double value, double first_error, double second_error) {
  if (!valid_error(first_error) || !valid_error(second_error))
    return std::numeric_limits<double>::infinity();
  if (first_error == 0.0 && second_error == 0.0) return 0.0;
  if (!std::numeric_limits<double>::is_iec559 || !std::isfinite(value))
    return std::numeric_limits<double>::infinity();
  const long double error = static_cast<long double>(first_error) + second_error;
  // Bound both binary64 roundings (fast and reference), including gradual
  // underflow. Four epsilons conservatively cover the 1/(1-eps) envelope.
  const long double epsilon = std::numeric_limits<double>::epsilon();
  return upward(error + 4.0L * epsilon *
    (std::abs(static_cast<long double>(value)) + error) +
    4.0L * std::numeric_limits<double>::denorm_min());
}
} // namespace detail

inline double predictive_error(double q, double numerator_error,
                               double denominator_error) {
  return detail::operation_error(q, numerator_error, denominator_error);
}

// Call with the newly computed sum and the two operand error bounds. Marker
// accumulation order is unchanged. An exact reference operand pair yields
// zero error, so reference controls remain exact rather than accumulating an
// unnecessary mathematical-roundoff interval around their own arithmetic.
inline double accumulated_error(double new_sum, double previous_error,
                                double marker_error) {
  return detail::operation_error(new_sum, previous_error, marker_error);
}

struct Estimate {
  double value = std::numeric_limits<double>::quiet_NaN();
  double reference_error = 0.0;
  bool complete_available = false;
};

inline bool potentially_approximate(const Estimate& value) {
  return value.complete_available && std::isfinite(value.value) &&
    value.reference_error != 0.0;
}

// All current raw consumers use differences, or a difference from a maximum.
// Preserving their pairwise threshold comparisons preserves these consumers
// for EVERY subset of rows surviving independent relationship gates, including
// incumbent retention/replacement. Include Tassign, identity margin, tol,
// and +/-tie_tolerance as supplied by the current chooser. Include zero when
// a consumer actually uses a strict-zero decision. A max support > Tassign
// is equivalent to all corresponding pairwise contrasts > Tassign; exact
// ordering inside a chooser's existing tie region need not be refined.
// Both orientations
// of a contrast are tested; thresholds are never enlarged by an error bound.
// A marked approximate row must be recomputed over its COMPLETE fixed marker
// panel before support rejection, chooser eligibility or mutation. Do not
// replace its score by an interval endpoint or a favorable bound.
inline std::vector<unsigned char> refinement_mask(
    const std::vector<Estimate>& values, const std::vector<double>& boundaries) {
  std::vector<unsigned char> mask(values.size(), 0U);
  bool any = false, valid_boundaries = true, valid_errors = true;
  double largest_error = 0.0;
  long double largest_magnitude = 0.0L;
  for (double threshold : boundaries) valid_boundaries &= std::isfinite(threshold);
  for (std::size_t at = 0U; at < values.size(); ++at) {
    any = any || potentially_approximate(values[at]);
    if (values[at].complete_available && std::isfinite(values[at].value)) {
      valid_errors &= detail::valid_error(values[at].reference_error);
      largest_error = std::max(largest_error, values[at].reference_error);
      largest_magnitude = std::max(largest_magnitude,
        std::abs(static_cast<long double>(values[at].value)));
    }
  }
  if (!any) return mask;
  const long double error = static_cast<long double>(largest_error) * 2.0L;
  const double contrast_error = detail::upward(error +
    4.0L * std::numeric_limits<double>::epsilon() *
      (2.0L * largest_magnitude + error) +
    4.0L * std::numeric_limits<double>::denorm_min());
  if (!valid_boundaries || !valid_errors || !std::isfinite(contrast_error)) {
    for (std::size_t at = 0U; at < values.size(); ++at)
      mask[at] = potentially_approximate(values[at]) ? 1U : 0U;
    return mask;
  }
  std::vector<std::pair<double, std::size_t>> sorted;
  sorted.reserve(values.size());
  for (std::size_t at = 0U; at < values.size(); ++at)
    if (values[at].complete_available && std::isfinite(values[at].value))
      sorted.emplace_back(values[at].value, at);
  std::sort(sorted.begin(), sorted.end());
  // Range marks avoid enumerating every pair of equal/near-equal rows.
  // Storage is O(C), time O(K*C*log C), including a 2,500-row tie cluster.
  std::vector<long long> ranges(sorted.size() + 1U, 0);
  for (std::size_t first = 0U; first < sorted.size(); ++first) {
    const std::size_t owner = sorted[first].second;
    if (!potentially_approximate(values[owner])) continue;
    for (double threshold : boundaries) for (int orientation : {-1, 1}) {
      const long double target = static_cast<long double>(sorted[first].first) -
        static_cast<long double>(orientation) * threshold;
      const auto left = std::lower_bound(sorted.begin(), sorted.end(), target - contrast_error,
        [](const std::pair<double, std::size_t>& entry, long double bound) {
          return static_cast<long double>(entry.first) < bound;
        });
      const auto right = std::upper_bound(left, sorted.end(), target + contrast_error,
        [](long double bound, const std::pair<double, std::size_t>& entry) {
          return bound < static_cast<long double>(entry.first);
        });
      if (left == right || (right - left == 1 && left->second == owner)) continue;
      mask[owner] = 1U;
      ++ranges[static_cast<std::size_t>(left - sorted.begin())];
      --ranges[static_cast<std::size_t>(right - sorted.begin())];
    }
  }
  long long covered = 0;
  for (std::size_t at = 0U; at < sorted.size(); ++at) {
    covered += ranges[at];
    if (covered > 0 && potentially_approximate(values[sorted[at].second]))
      mask[sorted[at].second] = 1U;
  }
  return mask;
}

} // namespace predictive_error_bounds
} // namespace pap_internal
#endif
