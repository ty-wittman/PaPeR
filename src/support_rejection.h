#ifndef PAP_SUPPORT_REJECTION_H
#define PAP_SUPPORT_REJECTION_H

// Negative proof for a strict raw conditional-support AND gate. Each input is
// the complete ordered marker panel on the same frozen candidate background.
// An unavailable panel or zero/zero contrast is unknown, never a rejection.
// This helper neither changes thresholds nor certifies a relationship class.
#include <algorithm>
#include <cmath>
#include <limits>

namespace pap_internal {
namespace support_rejection {

struct Panel {
  double ln = std::numeric_limits<double>::quiet_NaN();
  bool available = false;
  bool complete = false;
};

struct Proof {
  bool rejected = false;
  bool all_defined = false;
  double dam_ln = std::numeric_limits<double>::quiet_NaN();
  double sire_ln = std::numeric_limits<double>::quiet_NaN();
  double pair_ln = std::numeric_limits<double>::quiet_NaN();
};

inline bool defined(const Panel& panel) {
  return panel.available && panel.complete &&
    (std::isfinite(panel.ln) || panel.ln == -std::numeric_limits<double>::infinity());
}

inline double contrast(const Panel& numerator, const Panel& denominator) {
  const double zero = -std::numeric_limits<double>::infinity();
  if (!defined(numerator) || !defined(denominator) ||
      (numerator.ln == zero && denominator.ln == zero))
    return std::numeric_limits<double>::quiet_NaN();
  return numerator.ln == zero ? zero : numerator.ln - denominator.ln;
}

inline Proof single(const Panel& parent, const Panel& null, int side,
                    double threshold) {
  Proof proof;
  if ((side != 1 && side != 2) || !std::isfinite(threshold)) return proof;
  const double value = contrast(parent, null);
  if (side == 1) proof.dam_ln = value; else proof.sire_ln = value;
  proof.all_defined = !std::isnan(value);
  proof.rejected = proof.all_defined && value <= threshold;
  return proof;
}

inline Proof pair(const Panel& parents, const Panel& dam, const Panel& sire,
                  const Panel& null, double threshold) {
  Proof proof;
  if (!std::isfinite(threshold)) return proof;
  proof.dam_ln = contrast(parents, sire);
  proof.sire_ln = contrast(parents, dam);
  if (defined(null) && defined(dam) && defined(sire)) {
    Panel baseline;
    baseline.ln = std::max(null.ln, std::max(dam.ln, sire.ln));
    baseline.available = baseline.complete = true;
    proof.pair_ln = contrast(parents, baseline);
  }
  const double support[] = {proof.dam_ln, proof.sire_ln, proof.pair_ln};
  proof.all_defined = true;
  for (double value : support) {
    proof.all_defined = proof.all_defined && !std::isnan(value);
    proof.rejected = proof.rejected || (!std::isnan(value) && value <= threshold);
  }
  return proof;
}

} // namespace support_rejection
} // namespace pap_internal
#endif
