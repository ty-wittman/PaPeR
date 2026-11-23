#ifndef PAP_PROBABILITY_LOG_H
#define PAP_PROBABILITY_LOG_H

#include <Rcpp.h>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>

// This header is also used by the scheduled native-worker evaluator.  Its
// numerical-contract failures therefore stay native exceptions; the Rcpp
// entry boundary translates them on the calling thread after every worker has
// joined.
[[noreturn]] static inline void probability_log_contract_error(
    const char* context, const char* requirement) {
  throw std::runtime_error(std::string(context) + " " + requirement);
}

// Stage-1 probability/log boundary.  Missing observations are resolved before
// this helper is reached; zero here is therefore a genuine structural zero.
static inline double log_probability_ln(double probability, const char* context) {
  if (!std::isfinite(probability) || probability < 0.0)
    probability_log_contract_error(context, "must be nonnegative finite");
  return probability == 0.0 ? R_NegInf : std::log(probability);
}

static inline bool is_negative_infinity_ln(double value) {
  return std::isinf(value) && value < 0.0;
}

static inline double log_likelihood_contrast_ln(double focal_ln,
                                                double alternative_ln,
                                                const char* context) {
  const bool focal_ok = std::isfinite(focal_ln) || is_negative_infinity_ln(focal_ln);
  const bool alternative_ok = std::isfinite(alternative_ln) ||
    is_negative_infinity_ln(alternative_ln);
  if (!focal_ok || !alternative_ok)
    probability_log_contract_error(context, "must be finite or -Inf");
  // Both hypotheses assign zero probability: the contrast is undefined.
  if (is_negative_infinity_ln(focal_ln) && is_negative_infinity_ln(alternative_ln))
    return NA_REAL;
  if (is_negative_infinity_ln(focal_ln)) return R_NegInf;
  if (is_negative_infinity_ln(alternative_ln)) return R_PosInf;
  return focal_ln - alternative_ln;
}


class LogSumExpAccumulator {
 public:
  LogSumExpAccumulator() : max_ln_(R_NegInf), sum_scaled_(0.0) {}

  void add_log(double value_ln, const char* context) {
    if (!(std::isfinite(value_ln) || is_negative_infinity_ln(value_ln)))
      probability_log_contract_error(context, "must be finite or -Inf");
    if (is_negative_infinity_ln(value_ln)) return;
    if (is_negative_infinity_ln(max_ln_)) {
      max_ln_ = value_ln;
      sum_scaled_ = 1.0;
    } else if (value_ln > max_ln_) {
      sum_scaled_ = sum_scaled_ * std::exp(max_ln_ - value_ln) + 1.0;
      max_ln_ = value_ln;
    } else {
      sum_scaled_ += std::exp(value_ln - max_ln_);
    }
  }

  // Exact hot-kernel variant for values whose finite-or--Inf contract was
  // validated at the owning snapshot/buffer boundary. The accumulator update
  // order and floating-point operations are identical to add_log().
  void add_log_prevalidated(double value_ln) {
    if (value_ln == R_NegInf) return;
    if (max_ln_ == R_NegInf) {
      max_ln_ = value_ln;
      sum_scaled_ = 1.0;
    } else if (value_ln > max_ln_) {
      sum_scaled_ = sum_scaled_ * std::exp(max_ln_ - value_ln) + 1.0;
      max_ln_ = value_ln;
    } else {
      sum_scaled_ += std::exp(value_ln - max_ln_);
    }
  }

  double value() const {
    return is_negative_infinity_ln(max_ln_) ? R_NegInf :
      max_ln_ + std::log(sum_scaled_);
  }

 private:
  double max_ln_;
  double sum_scaled_;
};

#endif
