#include "evidence_views.h"
namespace pap_internal::search {

thread_local bool native_error_scope = false;

bool finite_or_neg_inf(double value) {
  return std::isfinite(value) || is_negative_infinity_ln(value);
}

double log_add(double first, double second) {
  if (is_negative_infinity_ln(first)) return second;
  if (is_negative_infinity_ln(second)) return first;
  const double high = std::max(first, second);
  const double low = std::min(first, second);
  return high + std::log1p(std::exp(low - high));
}

double log_sum(const std::vector<double>& values, bool sorted,
                      const char* context) {
  if (values.empty()) return R_NegInf;
  std::vector<double> ordered;
  const std::vector<double>* source = &values;
  if (sorted) {
    ordered = values;
    std::sort(ordered.begin(), ordered.end());
    source = &ordered;
  }
  double maximum = R_NegInf;
  for (std::size_t index = 0; index < source->size(); ++index) {
    const double value = (*source)[index];
    if (!finite_or_neg_inf(value))
      engine_stop("%s must be finite or -Inf (index=%d, value=%g)", context,
           static_cast<int>(index), value);
    if (value > maximum) maximum = value;
  }
  if (is_negative_infinity_ln(maximum)) return R_NegInf;
  double total = 0.0;
  for (double value : *source)
    if (!is_negative_infinity_ln(value)) total += std::exp(value - maximum);
  return maximum + std::log(total);
}

int genotype_index_zero(int first, int second, int k) {
  if (first > second) std::swap(first, second);
  return first * k - first * (first - 1) / 2 + (second - first);
}

void normalize_in_place(std::vector<double>& value) {
  const double z = log_sum(value, false, "simp normalized message");
  if (is_negative_infinity_ln(z)) return;
  for (double& x : value) x -= z;
}

NativeErrorScope::NativeErrorScope() : previous(native_error_scope) { native_error_scope = true; }

NativeErrorScope::~NativeErrorScope() { native_error_scope = previous; }





















auto SideParentSlots::has_dam() const -> bool  { return dam > 0; }

auto SideParentSlots::has_sire() const -> bool  { return sire > 0; }

auto SideParentSlots::has_any() const -> bool  { return has_dam() || has_sire(); }

















}
