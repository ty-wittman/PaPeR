#include "pap_removed_path.h"
#ifndef PAP_SCORE_REJECTION_H
#define PAP_SCORE_REJECTION_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

// Optional execution-order proof only. This header calculates no likelihood,
// constructs no relationship graph, and defines no assignment threshold. The
// caller owns each request and panel for the lifetime of these descriptors.
namespace pap_internal {
namespace score_rejection {

enum class State : unsigned char {
  unscored, scored_positive, scored_zero, inapplicable, unavailable
};

struct Boundary {
  const void* request = nullptr;  // Exact owned common-D/C request, not its hash.
  const void* panel = nullptr;    // Exact compiled observation/marker panel.
  std::uint64_t graph_revision = 0;
  std::uint64_t message_revision = 0;
  int side = 0;                  // 1=dam request, 2=sire request.

  bool valid() const {
    return request != nullptr && panel != nullptr && (side == 1 || side == 2);
  }
  bool same(const Boundary& other) const {
    return valid() && other.valid() && request == other.request &&
      panel == other.panel && graph_revision == other.graph_revision &&
      message_revision == other.message_revision && side == other.side;
  }
};

struct Operand {
  Boundary boundary;
  State state = State::unscored;
  double partition_ln = std::numeric_limits<double>::quiet_NaN();
  double absolute_error_bound_ln = 0.0;
  std::size_t loci_accumulated = 0;
};

enum class Outcome : unsigned char {
  continue_complete, rejected, needs_reference, undefined
};

enum class Reason : unsigned char {
  none, incumbent, incomplete_operational_maximum, incomplete_panel,
  different_boundary, invalid_operand, undefined_u, undefined_contrast,
  invalid_numerical_bound, rejection_witness, threshold_overlap
};

struct Request {
  bool incumbent = false;
  // True only after the complete operational maximum has been established.
  // Direct PO is sufficient for the ordinary common-D/C simp side. An
  // unscored compound parental numerator must not be presumed <= direct PO.
  bool operational_maximum_certified = false;
  bool retained_reference_arithmetic = false;
  std::size_t expected_loci = 0;
  double Tassign_ln = std::numeric_limits<double>::quiet_NaN();
  // The caller uses the exact existing side reducer's expression/order,
  // including its side-local U and finite display reference. Do not replace
  // this value by a cross-side score, a per-locus maximum, or a raw objective.
  double contrast_ln = std::numeric_limits<double>::quiet_NaN();
  // Covers panel summation and display/reducer arithmetic not already in the
  // operands' error envelopes, relative to the retained reference decision.
  double additional_error_bound_ln = 0.0;
};

struct Certificate {
  Outcome outcome = Outcome::continue_complete;
  Reason reason = Reason::none;
  double witness_contrast_ln = std::numeric_limits<double>::quiet_NaN();
  double lower_ln = std::numeric_limits<double>::quiet_NaN();
  double upper_ln = std::numeric_limits<double>::quiet_NaN();
};

inline bool negative_infinity(double x) {
  return std::isinf(x) && x < 0.0;
}

inline bool valid_scored_operand(const Operand& operand) { pap_internal::removed_path("valid_scored_operand"); }

inline Certificate certify(const Request& request, const Operand& po,
                           const Operand& u, const Operand& alternative) {
  Certificate out;
  const auto decline = [&](Reason reason) {
    out.reason = reason;
    return out;
  };
  if (request.incumbent) return decline(Reason::incumbent);
  if (!request.operational_maximum_certified)
    return decline(Reason::incomplete_operational_maximum);
  if (request.expected_loci == 0 ||
      po.loci_accumulated != request.expected_loci ||
      u.loci_accumulated != request.expected_loci ||
      alternative.loci_accumulated != request.expected_loci)
    return decline(Reason::incomplete_panel);
  if (!po.boundary.same(u.boundary) ||
      !po.boundary.same(alternative.boundary))
    return decline(Reason::different_boundary);
  if (!valid_scored_operand(po) || !valid_scored_operand(u) ||
      !valid_scored_operand(alternative))
    return decline(Reason::invalid_operand);
  // A partial request cannot manufacture a defined PO/U contrast from U=0.
  if (u.state != State::scored_positive) {
    out.outcome = Outcome::undefined;
    out.reason = Reason::undefined_u;
    return out;
  }
  if (po.state == State::scored_zero &&
      alternative.state == State::scored_zero) {
    out.outcome = Outcome::undefined;
    out.reason = Reason::undefined_contrast;
    return out;
  }
  if (!std::isfinite(request.Tassign_ln))
    return decline(Reason::invalid_operand);
  if (po.state == State::scored_zero) {
    if (!negative_infinity(request.contrast_ln))
      return decline(Reason::invalid_operand);
    out.outcome = Outcome::rejected;
    out.reason = Reason::rejection_witness;
    out.witness_contrast_ln = out.lower_ln = out.upper_ln =
      -std::numeric_limits<double>::infinity();
    return out;
  }
  if (alternative.state == State::scored_zero)
    return decline(Reason::none);  // This alternative cannot reject finite PO.
  if (!std::isfinite(request.contrast_ln))
    return decline(Reason::invalid_operand);
  if (!std::isfinite(request.additional_error_bound_ln) ||
      request.additional_error_bound_ln < 0.0)
    return decline(Reason::invalid_numerical_bound);
  double error = 0.0;
  const double error_terms[] = {
    po.absolute_error_bound_ln, u.absolute_error_bound_ln,
    u.absolute_error_bound_ln, alternative.absolute_error_bound_ln,
    request.additional_error_bound_ln};
  for (double term : error_terms) {
    if (term == 0.0) continue;
    error = std::nextafter(error + term,
      std::numeric_limits<double>::infinity());
    if (!std::isfinite(error)) return decline(Reason::invalid_numerical_bound);
  }
  if (request.retained_reference_arithmetic && error != 0.0)
    return decline(Reason::invalid_numerical_bound);
  out.witness_contrast_ln = request.contrast_ln;
  if (request.retained_reference_arithmetic) {
    out.lower_ln = out.upper_ln = request.contrast_ln;
  } else {
    // An outward step also encloses rounding of the envelope and endpoints.
    error = std::nextafter(error, std::numeric_limits<double>::infinity());
    out.lower_ln = std::nextafter(request.contrast_ln - error,
      -std::numeric_limits<double>::infinity());
    out.upper_ln = std::nextafter(request.contrast_ln + error,
      std::numeric_limits<double>::infinity());
  }
  if (out.upper_ln <= request.Tassign_ln) {
    out.outcome = Outcome::rejected;
    out.reason = Reason::rejection_witness;
  } else if (out.lower_ln <= request.Tassign_ln) {
    out.outcome = Outcome::needs_reference;
    out.reason = Reason::threshold_overlap;
  }
  return out;
}

}  // namespace score_rejection
}  // namespace pap_internal

#endif
