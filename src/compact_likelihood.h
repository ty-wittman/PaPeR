#ifndef PAP_COMPACT_LIKELIHOOD_H
#define PAP_COMPACT_LIKELIHOOD_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "probability_log.h"

// Unordered genotype indexing: (1,1), (1,2), ..., (1,k), (2,2), ... .
static inline int compact_geno_index(int a, int b, int k) {
  if (a > b) { int swap = a; a = b; b = swap; }
  const int offset = (a - 1) * k - (a - 1) * (a - 2) / 2;
  return offset + (b - a);
}

static inline int compact_geno_count(int k) {
  return k * (k + 1) / 2;
}

struct CompactEmissionKernel {
  int k;
  double Err;
  double w_Het;
  double w_Hom;
  double w_Far;
  double Z_hom;
  double Z_het;
  int model;
  double Dropout;
  double Foreign;
  Rcpp::NumericVector q;
};

static const int COMPACT_MODEL_LEGACY = 0;
static const int COMPACT_MODEL_DROPOUT_FOREIGN = 1;

static inline CompactEmissionKernel compact_make_emission_kernel(
    int k, double Err, const Rcpp::NumericVector& profile) {
  if (k < 1) Rcpp::stop("allele count must be positive");
  if (!R_finite(Err) || Err <= 0.0 || Err > 1.0)
    Rcpp::stop("Err must be strictly positive and finite, at most 1; supply a small positive value");
  if (profile.size() != 3) Rcpp::stop("ErrProfile must contain exactly three weights");
  for (int q = 0; q < 3; ++q)
    if (!R_finite(profile[q]) || profile[q] < 0.0)
      Rcpp::stop("ErrProfile weights must be finite and nonnegative");

  CompactEmissionKernel out = {k, Err, profile[0], profile[1], profile[2], 0.0, 0.0};
  const double km1 = static_cast<double>(k - 1);
  out.Z_hom = km1 * out.w_Het + km1 * static_cast<double>(k) / 2.0 * out.w_Far;
  const double km2 = static_cast<double>(k - 2);
  out.Z_het = 2.0 * out.w_Hom + 2.0 * km2 * out.w_Het +
    km2 * km1 / 2.0 * out.w_Far;
  if (Err > 0.0) {
    if (k == 1) Rcpp::stop("a monomorphic locus cannot represent positive Err");
    if (!(out.Z_hom > 0.0) || !(out.Z_het > 0.0) ||
        !R_finite(out.Z_hom) || !R_finite(out.Z_het))
      Rcpp::stop("ErrProfile is non-normalizable for at least one true genotype row");
  }
  return out;
}

static inline CompactEmissionKernel compact_make_dropout_foreign_kernel(
    int k, double Dropout, double Foreign, const Rcpp::NumericVector& q) {
  if (k < 1) Rcpp::stop("allele count must be positive");
  if (!R_finite(Dropout) || Dropout <= 0.0 || Dropout > 1.0)
    Rcpp::stop("Dropout must be strictly positive and finite, at most 1; supply a small positive value");
  if (!R_finite(Foreign) || Foreign <= 0.0 || Foreign > 1.0)
    Rcpp::stop("Foreign must be strictly positive and finite, at most 1; supply a small positive value");
  if (q.size() != k)
    Rcpp::stop("ForeignQ must match the allele count");
  double total = 0.0;
  for (int a = 0; a < k; ++a) {
    if (!R_finite(q[a]) || q[a] < 0.0)
      Rcpp::stop("ForeignQ must be finite and nonnegative");
    total += q[a];
  }
  if (!R_finite(total) || total <= 0.0 || std::fabs(total - 1.0) > 1e-12)
    Rcpp::stop("ForeignQ must be normalized to sum to one");
  CompactEmissionKernel out;
  out.k = k;
  out.Err = 0.0;
  out.w_Het = 0.0;
  out.w_Hom = 0.0;
  out.w_Far = 0.0;
  out.Z_hom = 0.0;
  out.Z_het = 0.0;
  out.model = COMPACT_MODEL_DROPOUT_FOREIGN;
  out.Dropout = Dropout;
  out.Foreign = Foreign;
  out.q = q;
  return out;
}

static inline int compact_shared_allele_copies(int true_a, int true_b,
                                                int observed_a, int observed_b) {
  if (true_a == true_b)
    return (observed_a == true_a) + (observed_b == true_a);
  return (observed_a == true_a || observed_b == true_a) +
    (observed_a == true_b || observed_b == true_b);
}

static inline double compact_foreign_q_probability(
    const CompactEmissionKernel& kernel, int truth, int observed) {
  return (truth == observed ? 1.0 - kernel.Foreign : 0.0) +
    kernel.Foreign * kernel.q[observed - 1];
}

static inline double compact_foreign_residual_probability(
    const CompactEmissionKernel& kernel, int true_a, int true_b,
    int observed_a, int observed_b) {
  const double q_au = compact_foreign_q_probability(kernel, true_a, observed_a);
  const double q_av = compact_foreign_q_probability(kernel, true_a, observed_b);
  const double q_bu = compact_foreign_q_probability(kernel, true_b, observed_a);
  const double q_bv = compact_foreign_q_probability(kernel, true_b, observed_b);
  if (observed_a == observed_b) return q_au * q_bu;
  return q_au * q_bv + q_av * q_bu;
}

static inline double compact_emission_probability(
    const CompactEmissionKernel& kernel, int true_a, int true_b,
    int observed_a, int observed_b) {
  if (true_a > true_b) { int swap = true_a; true_a = true_b; true_b = swap; }
  if (observed_a > observed_b) {
    int swap = observed_a; observed_a = observed_b; observed_b = swap;
  }
  if (kernel.model == COMPACT_MODEL_DROPOUT_FOREIGN) {
    const double residual = compact_foreign_residual_probability(
      kernel, true_a, true_b, observed_a, observed_b);
    if (true_a == true_b) return residual;
    return (1.0 - kernel.Dropout) * residual +
      0.5 * kernel.Dropout * compact_foreign_residual_probability(
        kernel, true_a, true_a, observed_a, observed_b) +
      0.5 * kernel.Dropout * compact_foreign_residual_probability(
        kernel, true_b, true_b, observed_a, observed_b);
  }
  if (true_a == observed_a && true_b == observed_b) return 1.0 - kernel.Err;
  if (kernel.Err == 0.0) return 0.0;
  const int overlap = compact_shared_allele_copies(
    true_a, true_b, observed_a, observed_b);
  const double weight = overlap == 0 ? kernel.w_Far :
    (observed_a == observed_b ? kernel.w_Hom : kernel.w_Het);
  const double normalizer = true_a == true_b ? kernel.Z_hom : kernel.Z_het;
  return kernel.Err * weight / normalizer;
}

static inline void compact_validate_locus_table(const Rcpp::List& tab) {
  if (!tab.containsElementNamed("schema_version"))
    Rcpp::stop("locus table is missing schema_version; rebuild with build_tables()");
  const std::string schema = Rcpp::as<std::string>(tab["schema_version"]);
  const bool v1 = schema == "compact-locus-v1";
  const bool v2 = schema == "compact-locus-v2";
  if (!v1 && !v2)
    Rcpp::stop("locus table must use compact-locus-v1 or compact-locus-v2; rebuild with build_tables()");
  for (const char* name : {"S", "YS", "X", "Delta", "relnames"})
    if (tab.containsElementNamed(name))
      Rcpp::stop("locus table retains obsolete dense field '%s'; rebuild with build_tables()", name);
  for (const char* name : {"schema_id", "k", "G", "geno", "p", "pi", "mObs",
                           "p_oh_true_po", "p_oh_unrelated", "p_me_true_trio"})
    if (!tab.containsElementNamed(name))
      Rcpp::stop("compact locus table is missing required field '%s'", name);
  const int schema_id = Rcpp::as<int>(tab["schema_id"]);
  if ((v1 && schema_id != 1) || (v2 && schema_id != 2))
    Rcpp::stop("locus table schema_id does not match schema_version");
  const auto positive_rate = [&](const char* name) {
    SEXP value = tab[name];
    if ((TYPEOF(value) != REALSXP && TYPEOF(value) != INTSXP) || Rf_xlength(value) != 1)
      Rcpp::stop("%s must be one positive panel-wide scalar; supply a small positive value", name);
    const double rate = Rcpp::as<double>(value);
    if (!R_finite(rate) || rate <= 0.0 || rate > 1.0)
      Rcpp::stop("%s must be strictly positive and finite, at most 1; supply a small positive value", name);
  };
  if (v1) {
    for (const char* name : {"Err", "ErrProfile", "Z_hom", "Z_het"})
      if (!tab.containsElementNamed(name))
        Rcpp::stop("compact locus table is missing required field '%s'", name);
    for (const char* name : {"ErrorModel", "Dropout", "Foreign", "ForeignQ",
                             "ForeignQMode"})
      if (tab.containsElementNamed(name))
        Rcpp::stop("legacy compact locus table contains new-model field '%s'", name);
    positive_rate("Err");
  } else {
    for (const char* name : {"ErrorModel", "Dropout", "Foreign", "ForeignQ",
                             "ForeignQMode"})
      if (!tab.containsElementNamed(name))
        Rcpp::stop("dropout_foreign compact locus table is missing required field '%s'", name);
    if (Rcpp::as<std::string>(tab["ErrorModel"]) != "dropout_foreign")
      Rcpp::stop("compact-locus-v2 requires ErrorModel='dropout_foreign'");
    for (const char* name : {"Err", "ErrProfile", "Z_hom", "Z_het"})
      if (tab.containsElementNamed(name))
        Rcpp::stop("dropout_foreign compact locus table contains legacy field '%s'", name);
    const std::string q_mode = Rcpp::as<std::string>(tab["ForeignQMode"]);
    if (q_mode != "p" && q_mode != "supplied")
      Rcpp::stop("ForeignQMode must be 'p' or 'supplied'");
    positive_rate("Dropout");
    positive_rate("Foreign");
  }
}

static inline CompactEmissionKernel compact_kernel_from_table(const Rcpp::List& tab) {
  compact_validate_locus_table(tab);
  const std::string schema = Rcpp::as<std::string>(tab["schema_version"]);
  if (schema == "compact-locus-v2") {
    return compact_make_dropout_foreign_kernel(
      Rcpp::as<int>(tab["k"]), Rcpp::as<double>(tab["Dropout"]),
      Rcpp::as<double>(tab["Foreign"]),
      Rcpp::as<Rcpp::NumericVector>(tab["ForeignQ"]));
  }
  const CompactEmissionKernel kernel = compact_make_emission_kernel(
    Rcpp::as<int>(tab["k"]), Rcpp::as<double>(tab["Err"]),
    Rcpp::as<Rcpp::NumericVector>(tab["ErrProfile"]));
  const double stored_hom = Rcpp::as<double>(tab["Z_hom"]);
  const double stored_het = Rcpp::as<double>(tab["Z_het"]);
  if (std::fabs(stored_hom - kernel.Z_hom) > 1e-12 ||
      std::fabs(stored_het - kernel.Z_het) > 1e-12)
    Rcpp::stop("locus table has inconsistent error-kernel normalizers");
  return kernel;
}

static inline double compact_founder_probability(const Rcpp::NumericVector& p,
                                                 int a, int b) {
  return a == b ? p[a - 1] * p[a - 1] : 2.0 * p[a - 1] * p[b - 1];
}

struct CompactObservedDescriptor {
  std::vector<double> d;
  std::vector<double> r;
  std::vector<double> log_d;
  std::vector<double> log_r;
  double h;
  double m;
  double log_h;
  double log_m;
};

// Construct the observed-state descriptor in O(k).  The three-weight error
// kernel means that, for a fixed transmitted allele, all non-exact heterozygous
// truths fall into one shared-allele or no-shared-allele class.  Summing their
// allele-frequency masses directly avoids an O(k^2) truth-by-observation loop.
static inline CompactObservedDescriptor compact_observed_descriptor_from_parts(
    const int k, const Rcpp::NumericVector& p,
    const Rcpp::IntegerMatrix& geno, const CompactEmissionKernel& kernel,
    int observed) {
  CompactObservedDescriptor out;
  out.d.assign(k + 1, 1.0);
  out.r.assign(k + 1, 1.0);
  out.log_d.assign(k + 1, 0.0);
  out.log_r.assign(k + 1, 0.0);
  out.h = 1.0;
  out.m = 1.0;
  out.log_h = 0.0;
  out.log_m = 0.0;
  if (observed < 0) return out;
  const int G = compact_geno_count(k);
  if (observed >= G) Rcpp::stop("observed genotype index is out of range");
  const int observed_a = geno(observed, 0), observed_b = geno(observed, 1);
  const bool observed_hom = observed_a == observed_b;
  if (kernel.model == COMPACT_MODEL_DROPOUT_FOREIGN) {
    // For a fixed transmitted allele a, sum the residual observation kernel
    // over the other founder allele in O(k).  These two single-allele masses
    // are sufficient for S_ab[a]; S_bb is shared by every a.
    double s_u = 0.0, s_v = 0.0, s_bb = 0.0;
    for (int b = 1; b <= k; ++b) {
      const double p_b = p[b - 1];
      const double q_bu = compact_foreign_q_probability(kernel, b, observed_a);
      const double q_bv = compact_foreign_q_probability(kernel, b, observed_b);
      s_u += p_b * q_bu;
      s_v += p_b * q_bv;
      s_bb += p_b * (observed_hom ? q_bu * q_bu : 2.0 * q_bu * q_bv);
    }
    out.h = 0.0;
    out.m = 0.0;
    for (int a = 1; a <= k; ++a) {
      const double q_au = compact_foreign_q_probability(
        kernel, a, observed_a);
      const double q_av = compact_foreign_q_probability(
        kernel, a, observed_b);
      const double s_ab = observed_hom ? q_au * s_u :
        q_au * s_v + q_av * s_u;
      out.d[a] = compact_foreign_residual_probability(
        kernel, a, a, observed_a, observed_b);
      out.r[a] = (1.0 - kernel.Dropout) * s_ab +
        0.5 * kernel.Dropout * out.d[a] +
        0.5 * kernel.Dropout * s_bb;
      out.log_d[a] = out.d[a] == 0.0 ? R_NegInf : std::log(out.d[a]);
      out.log_r[a] = out.r[a] == 0.0 ? R_NegInf : std::log(out.r[a]);
      out.h += p[a - 1] * out.d[a];
      out.m += p[a - 1] * out.r[a];
    }
    out.log_h = out.h == 0.0 ? R_NegInf : std::log(out.h);
    out.log_m = out.m == 0.0 ? R_NegInf : std::log(out.m);
    return out;
  }
  const double shared_het_error = kernel.Err == 0.0 ? 0.0 :
    kernel.Err * (observed_hom ? kernel.w_Hom : kernel.w_Het) / kernel.Z_het;
  const double far_het_error = kernel.Err == 0.0 ? 0.0 :
    kernel.Err * kernel.w_Far / kernel.Z_het;
  out.h = 0.0;
  out.m = 0.0;
  for (int a = 1; a <= k; ++a) {
    out.d[a] = compact_emission_probability(kernel, a, a, observed_a, observed_b);
    out.log_d[a] = out.d[a] == 0.0 ? R_NegInf : std::log(out.d[a]);
    const double p_a = p[a - 1];
    const bool observed_contains_a = a == observed_a || a == observed_b;
    double hetero_mass = 0.0;
    if (observed_contains_a) {
      const double exact_mass = !observed_hom ?
        p[(a == observed_a ? observed_b : observed_a) - 1] : 0.0;
      hetero_mass = exact_mass * (1.0 - kernel.Err) +
        (1.0 - p_a - exact_mass) * shared_het_error;
    } else {
      const double shared_mass = p[observed_a - 1] +
        (observed_hom ? 0.0 : p[observed_b - 1]);
      hetero_mass = shared_mass * shared_het_error +
        (1.0 - p_a - shared_mass) * far_het_error;
    }
    out.r[a] = p_a * out.d[a] + hetero_mass;
    if (out.r[a] > 0.0) {
      out.log_r[a] = std::log(out.r[a]);
    } else {
      LogSumExpAccumulator r_log_sum;
      for (int b = 1; b <= k; ++b) {
        const double emission = compact_emission_probability(
          kernel, a, b, observed_a, observed_b);
        if (emission > 0.0)
          r_log_sum.add_log(std::log(p[b - 1]) + std::log(emission),
                            "observed descriptor probability");
      }
      out.log_r[a] = r_log_sum.value();
    }
    out.h += p[a - 1] * out.d[a];
    out.m += p[a - 1] * out.r[a];
  }
  if (out.h > 0.0) {
    out.log_h = std::log(out.h);
  } else {
    LogSumExpAccumulator h_log_sum;
    for (int a = 1; a <= k; ++a)
      if (!is_negative_infinity_ln(out.log_d[a]))
        h_log_sum.add_log(std::log(p[a - 1]) + out.log_d[a],
                          "homozygote descriptor probability");
    out.log_h = h_log_sum.value();
  }
  if (out.m > 0.0) {
    out.log_m = std::log(out.m);
  } else {
    LogSumExpAccumulator m_log_sum;
    for (int a = 1; a <= k; ++a)
      if (!is_negative_infinity_ln(out.log_r[a]))
        m_log_sum.add_log(std::log(p[a - 1]) + out.log_r[a],
                           "observed descriptor probability");
    out.log_m = m_log_sum.value();
  }
  return out;
}

static inline CompactObservedDescriptor compact_observed_descriptor(
    const Rcpp::List& tab, int observed) {
  const int k = Rcpp::as<int>(tab["k"]);
  const Rcpp::NumericVector p = tab["p"];
  const Rcpp::IntegerMatrix geno = tab["geno"];
  const CompactEmissionKernel kernel = compact_kernel_from_table(tab);
  return compact_observed_descriptor_from_parts(k, p, geno, kernel, observed);
}

static inline double compact_d7_probability(const Rcpp::List& tab, int left,
                                            int right) {
  if (left < 0 && right < 0) return 1.0;
  const int k = Rcpp::as<int>(tab["k"]);
  const Rcpp::NumericVector p = tab["p"];
  const Rcpp::IntegerMatrix geno = tab["geno"];
  const CompactEmissionKernel kernel = compact_kernel_from_table(tab);
  const int G = Rcpp::as<int>(tab["G"]);
  if (left >= G || right >= G) Rcpp::stop("observed genotype index is out of range");
  double out = 0.0;
  for (int a = 1; a <= k; ++a) for (int b = 1; b <= k; ++b) {
    const double left_emission = left < 0 ? 1.0 : compact_emission_probability(
      kernel, a, b, geno(left, 0), geno(left, 1));
    const double right_emission = right < 0 ? 1.0 : compact_emission_probability(
      kernel, a, b, geno(right, 0), geno(right, 1));
    out += p[a - 1] * p[b - 1] * left_emission * right_emission;
  }
  return out;
}

static inline std::vector<double> compact_pair_state_probabilities(
    const Rcpp::List& tab, int left, int right) {
  const int k = Rcpp::as<int>(tab["k"]);
  const Rcpp::NumericVector p = tab["p"];
  const CompactObservedDescriptor i = compact_observed_descriptor(tab, left);
  const CompactObservedDescriptor j = compact_observed_descriptor(tab, right);
  std::vector<double> out(9, 0.0);
  for (int a = 1; a <= k; ++a) {
    const double pa = p[a - 1];
    out[0] += pa * i.d[a] * j.d[a];
    out[2] += pa * i.d[a] * j.r[a];
    out[4] += pa * i.r[a] * j.d[a];
    out[7] += pa * i.r[a] * j.r[a];
  }
  out[1] = i.h * j.h;
  out[3] = i.h * j.m;
  out[5] = i.m * j.h;
  out[6] = compact_d7_probability(tab, left, right);
  out[8] = i.m * j.m;
  return out;
}

#endif
