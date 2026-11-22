// ============================================================================
// Compact Stage-1 likelihood kernel.  Static locus tables are O(G), never
// O(G^2): relationship coefficients live once per run and observed pair cells
// are evaluated directly from the nine condensed identity-state formulas.
// ============================================================================

#include <Rcpp.h>
#include <algorithm>
#include <cfloat>
#include <map>
#include <vector>
#include "latent_mendelian_support.h"
#include "probability_log.h"
#include "compact_likelihood.h"
using namespace Rcpp;

static inline int observed_genotype_index(int a, int b, int k) {
  return (a <= 0 || b <= 0) ? -1 : compact_geno_index(a, b, k);
}


static double search_logsumexp_values(std::vector<double>& values,
                                    bool sorted, const char* context) {
  if (values.empty()) return R_NegInf;
  if (sorted) std::sort(values.begin(), values.end());
  double maximum = R_NegInf;
  for (const double value : values) {
    if (!(R_finite(value) || is_negative_infinity_ln(value)))
      stop("%s must be finite or -Inf", context);
    if (value > maximum) maximum = value;
  }
  if (is_negative_infinity_ln(maximum)) return R_NegInf;
  double sum = 0.0;
  for (const double value : values)
    if (!is_negative_infinity_ln(value)) sum += std::exp(value - maximum);
  return maximum + std::log(sum);
}

// Parent-joint construction used by the finalized simp boundary recursion.
// Keeping this bounded G x G operation beside the existing child kernels
// removes repeated R outer/matrix/log-normalization allocations while
// preserving the same natural-log equation and column-major traversal.
// [[Rcpp::export]]
NumericMatrix parent_joint_log_cpp(
    NumericVector prior_ln, NumericVector log_p, IntegerMatrix child_index,
    NumericVector first_ln, NumericVector second_ln,
    double k0, double k1, double k2) {
  const int G = prior_ln.size();
  const int k = log_p.size();
  if (k < 1 || G != k * (k + 1) / 2)
    stop("simp parent joint prior dimensions are invalid");
  if (child_index.nrow() != k || child_index.ncol() != k)
    stop("simp parent joint child index dimensions are invalid");
  if (first_ln.size() != G || second_ln.size() != G)
    stop("simp parent joint messages must have length G");
  if (!R_finite(k0) || !R_finite(k1) || !R_finite(k2) ||
      k0 < 0.0 || k1 < 0.0 || k2 < 0.0)
    stop("simp parent joint IBD weights must be nonnegative finite");
  for (int g = 0; g < G; ++g) {
    if (!(R_finite(prior_ln[g]) || is_negative_infinity_ln(prior_ln[g])) ||
        !(R_finite(first_ln[g]) || is_negative_infinity_ln(first_ln[g])) ||
        !(R_finite(second_ln[g]) || is_negative_infinity_ln(second_ln[g])))
      stop("simp parent joint inputs must be finite or -Inf");
  }
  if (k1 == 0.0 && k2 == 0.0) {
    NumericMatrix out(G, G);
    for (int right = 0; right < G; ++right)
      for (int left = 0; left < G; ++left)
        out(left, right) = first_ln[left] + second_ln[right];
    return out;
  }
  std::vector<double> prior(static_cast<std::size_t>(G));
  for (int g = 0; g < G; ++g) prior[g] = std::exp(prior_ln[g]);
  std::vector<double> probability(static_cast<std::size_t>(G) * G, 0.0);
  for (int right = 0; right < G; ++right)
    for (int left = 0; left < G; ++left)
      probability[left + right * G] = k0 * prior[left] * prior[right];
  if (k1 > 0.0) {
    std::vector<double> allele(static_cast<std::size_t>(k));
    for (int a = 0; a < k; ++a) allele[a] = std::exp(log_p[a]);
    for (int shared = 0; shared < k; ++shared)
      for (int left = 0; left < k; ++left)
        for (int right = 0; right < k; ++right) {
          const int g1 = child_index(std::min(shared, left),
                                      std::max(shared, left)) - 1;
          const int g2 = child_index(std::min(shared, right),
                                      std::max(shared, right)) - 1;
          probability[g1 + g2 * G] +=
            k1 * allele[shared] * allele[left] * allele[right];
        }
  }
  if (k2 > 0.0)
    for (int g = 0; g < G; ++g)
      probability[g + g * G] += k2 * prior[g];
  std::vector<double> values(static_cast<std::size_t>(G) * G, R_NegInf);
  for (int right = 0; right < G; ++right)
    for (int left = 0; left < G; ++left) {
      const double mass = probability[left + right * G];
      if (mass > 0.0)
        values[left + right * G] = std::log(mass) +
          first_ln[left] - prior_ln[left] +
          second_ln[right] - prior_ln[right];
    }
  const double normalizer = search_logsumexp_values(
    values, false, "simp parent joint normalization");
  if (!is_negative_infinity_ln(normalizer))
    for (double& value : values) value -= normalizer;
  NumericMatrix out(G, G);
  std::copy(values.begin(), values.end(), out.begin());
  return out;
}


// Exact OH masses for ordered observed pairs.  For each first observed state o,
// D[a, o] is the conditional mass of a second observed genotype that is
// allele-disjoint from o given transmitted allele a.  Summing P(o, o2) over
// ordered pairs is correct because a PO pair has labelled parent/offspring
// observations; the unrelated calculation uses the same ordered convention.
// No G x G observed-pair matrix is constructed.
static std::vector<double> compact_oh_probabilities(
    const int k, const NumericVector& p, const IntegerMatrix& geno,
    const CompactEmissionKernel& kernel) {
  const int G = compact_geno_count(k);
  std::vector<CompactObservedDescriptor> desc(G);
  std::vector<double> containment((k + 1) * (k + 1), 0.0);
  std::vector<double> marginal_containment(k + 1, 0.0);
  const auto c_at = [&containment, k](int a, int u) -> double& {
    return containment[a * (k + 1) + u];
  };
  for (int observed = 0; observed < G; ++observed) {
    desc[observed] = compact_observed_descriptor_from_parts(
      k, p, geno, kernel, observed);
    const int u = geno(observed, 0), v = geno(observed, 1);
    for (int a = 1; a <= k; ++a) {
      c_at(a, u) += desc[observed].r[a];
      if (v != u) c_at(a, v) += desc[observed].r[a];
    }
    marginal_containment[u] += desc[observed].m;
    if (v != u) marginal_containment[v] += desc[observed].m;
  }
  double p_true_po = 0.0, p_unrelated = 0.0;
  for (int observed = 0; observed < G; ++observed) {
    const int u = geno(observed, 0), v = geno(observed, 1);
    const double disjoint_unrelated = v == u ?
      1.0 - marginal_containment[u] :
      1.0 - marginal_containment[u] - marginal_containment[v] + desc[observed].m;
    p_unrelated += desc[observed].m * disjoint_unrelated;
    for (int a = 1; a <= k; ++a) {
      const double disjoint = v == u ? 1.0 - c_at(a, u) :
        1.0 - c_at(a, u) - c_at(a, v) + desc[observed].r[a];
      p_true_po += p[a - 1] * desc[observed].r[a] * disjoint;
    }
  }
  const auto clamp_probability = [](double x, const char* name) {
    const double tolerance = 128.0 * DBL_EPSILON;
    if (x < -tolerance || x > 1.0 + tolerance || !R_finite(x))
      stop("compact %s OH probability is invalid", name);
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
  };
  return std::vector<double>{clamp_probability(p_true_po, "true-PO"),
                             clamp_probability(p_unrelated, "unrelated")};
}

// Exact probability of an observed Mendelian-incompatible trio conditional on
// all three observations being nonmissing.  The two transmitted alleles are
// summed directly, while each nontransmitted parental homologue is integrated
// through the compact r_o[a] descriptors.  This retains no observed-trio
// table, G x G matrix, or G^3 object: the temporary work is O(G k) for the
// descriptors/containment masses and O(G k^2) for the final contraction.
static double compact_true_trio_me_probability(
    const int k, const NumericVector& p, const IntegerMatrix& geno,
    const CompactEmissionKernel& kernel) {
  const int G = compact_geno_count(k);
  std::vector<CompactObservedDescriptor> desc(G);
  std::vector<double> containment((k + 1) * (k + 1), 0.0);
  const auto c_at = [&containment, k](int a, int u) -> double& {
    return containment[a * (k + 1) + u];
  };
  for (int observed = 0; observed < G; ++observed) {
    desc[observed] = compact_observed_descriptor_from_parts(
      k, p, geno, kernel, observed);
    const int u = geno(observed, 0), v = geno(observed, 1);
    for (int a = 1; a <= k; ++a) {
      c_at(a, u) += desc[observed].r[a];
      if (v != u) c_at(a, v) += desc[observed].r[a];
    }
  }

  long double compatible = 0.0L;
  for (int a = 1; a <= k; ++a) for (int b = 1; b <= k; ++b) {
    long double observed_compatible = 0.0L;
    for (int observed = 0; observed < G; ++observed) {
      const int x = geno(observed, 0), y = geno(observed, 1);
      const double emission = compact_emission_probability(
        kernel, a, b, x, y);
      if (emission == 0.0) continue;
      double K = 0.0;
      if (x == y) {
        K = c_at(a, x) * c_at(b, x);
      } else {
        K = c_at(a, x) * c_at(b, y) +
          c_at(a, y) * c_at(b, x) -
          desc[observed].r[a] * desc[observed].r[b];
      }
      const double tolerance = 1024.0 * DBL_EPSILON;
      if (!R_finite(K) || K < -tolerance || K > 1.0 + tolerance)
        stop("compact true-trio compatibility probability is invalid");
      if (K < 0.0) K = 0.0;
      if (K > 1.0) K = 1.0;
      observed_compatible += static_cast<long double>(emission) *
        static_cast<long double>(K);
    }
    compatible += static_cast<long double>(p[a - 1]) *
      static_cast<long double>(p[b - 1]) * observed_compatible;
  }

  const long double tolerance = 1024.0L *
    static_cast<long double>(DBL_EPSILON);
  if (!std::isfinite(compatible) || compatible < -tolerance ||
      compatible > 1.0L + tolerance)
    stop("compact true-trio compatible probability is invalid");
  if (compatible < 0.0L) compatible = 0.0L;
  if (compatible > 1.0L) compatible = 1.0L;
  long double result = 1.0L - compatible;
  if (result < 0.0L && result > -tolerance) result = 0.0L;
  if (result > 1.0L && result < 1.0L + tolerance) result = 1.0L;
  if (result < 0.0L || result > 1.0L || !std::isfinite(result))
    stop("compact true-trio Mendelian-error probability is invalid");
  return static_cast<double>(result);
}

// [[Rcpp::export]]
List build_locus_table_cpp(NumericVector p, double Err, NumericVector profile) {
  const int k = p.size();
  const int G = compact_geno_count(k);
  if (k < 1 || p.size() < 1) stop("allele frequency vector must be nonempty");
  double p_total = 0.0;
  for (int a = 0; a < k; ++a) {
    if (!R_finite(p[a]) || p[a] < 0.0) stop("allele frequencies must be nonnegative finite");
    p_total += p[a];
  }
  if (!(p_total > 0.0) || !R_finite(p_total)) stop("allele frequencies must have positive finite total");
  for (int a = 0; a < k; ++a) p[a] /= p_total;

  const CompactEmissionKernel kernel = compact_make_emission_kernel(k, Err, profile);
  IntegerMatrix geno(G, 2);
  NumericVector pi(G), mObs(G);
  for (int a = 1; a <= k; ++a) for (int b = a; b <= k; ++b) {
    const int g = compact_geno_index(a, b, k);
    geno(g, 0) = a;
    geno(g, 1) = b;
    pi[g] = compact_founder_probability(p, a, b);
  }
  for (int observed = 0; observed < G; ++observed)
    mObs[observed] = compact_observed_descriptor_from_parts(
      k, p, geno, kernel, observed).m;
  const std::vector<double> oh = compact_oh_probabilities(k, p, geno, kernel);
  const double p_me_true_trio = compact_true_trio_me_probability(
    k, p, geno, kernel);

  return List::create(
    _["schema_version"] = "compact-locus-v1",
    _["schema_id"] = 1L,
    _["k"] = k,
    _["G"] = G,
    _["geno"] = geno,
    _["p"] = p,
    _["pi"] = pi,
    _["mObs"] = mObs,
    _["Err"] = Err,
    _["ErrProfile"] = profile,
    _["Z_hom"] = kernel.Z_hom,
    _["Z_het"] = kernel.Z_het,
    _["p_oh_true_po"] = oh[0],
    _["p_oh_unrelated"] = oh[1],
    _["p_me_true_trio"] = p_me_true_trio);
}

// [[Rcpp::export]]
List build_locus_table_dropout_foreign_cpp(NumericVector p, double Dropout,
                                           double Foreign, NumericVector q) {
  const int k = p.size();
  const int G = compact_geno_count(k);
  if (k < 1 || p.size() < 1) stop("allele frequency vector must be nonempty");
  double p_total = 0.0;
  for (int a = 0; a < k; ++a) {
    if (!R_finite(p[a]) || p[a] < 0.0)
      stop("allele frequencies must be nonnegative finite");
    p_total += p[a];
  }
  if (!(p_total > 0.0) || !R_finite(p_total))
    stop("allele frequencies must have positive finite total");
  for (int a = 0; a < k; ++a) p[a] /= p_total;

  const CompactEmissionKernel kernel = compact_make_dropout_foreign_kernel(
    k, Dropout, Foreign, q);
  IntegerMatrix geno(G, 2);
  NumericVector pi(G), mObs(G);
  for (int a = 1; a <= k; ++a) for (int b = a; b <= k; ++b) {
    const int g = compact_geno_index(a, b, k);
    geno(g, 0) = a;
    geno(g, 1) = b;
    pi[g] = compact_founder_probability(p, a, b);
  }
  for (int observed = 0; observed < G; ++observed)
    mObs[observed] = compact_observed_descriptor_from_parts(
      k, p, geno, kernel, observed).m;
  const std::vector<double> oh = compact_oh_probabilities(k, p, geno, kernel);
  const double p_me_true_trio = compact_true_trio_me_probability(
    k, p, geno, kernel);

  return List::create(
    _["schema_version"] = "compact-locus-v2",
    _["schema_id"] = 2L,
    _["k"] = k,
    _["G"] = G,
    _["geno"] = geno,
    _["p"] = p,
    _["pi"] = pi,
    _["mObs"] = mObs,
    _["ErrorModel"] = "dropout_foreign",
    _["Dropout"] = Dropout,
    _["Foreign"] = Foreign,
    _["ForeignQ"] = q,
    _["ForeignQMode"] = "supplied",
    _["p_oh_true_po"] = oh[0],
    _["p_oh_unrelated"] = oh[1],
    _["p_me_true_trio"] = p_me_true_trio);
}

// [[Rcpp::export]]
NumericVector emission_probabilities_cpp(List tab, IntegerVector true_index,
                                         IntegerVector observed_index) {
  if (true_index.size() != observed_index.size())
    stop("true_index and observed_index must have equal length");
  compact_validate_locus_table(tab);
  const int G = tab["G"];
  const IntegerMatrix geno = tab["geno"];
  const CompactEmissionKernel kernel = compact_kernel_from_table(tab);
  NumericVector out(true_index.size());
  for (int q = 0; q < true_index.size(); ++q) {
    if (true_index[q] < 0 || true_index[q] >= G ||
        observed_index[q] < 0 || observed_index[q] >= G)
      stop("emission genotype index is out of range");
    out[q] = compact_emission_probability(
      kernel, geno(true_index[q], 0), geno(true_index[q], 1),
      geno(observed_index[q], 0), geno(observed_index[q], 1));
  }
  return out;
}


static NumericMatrix legacy_delta() {
  NumericMatrix out(6, 9);
  const double values[6][9] = {
    {0, 0, 0, 0, 0, 0, 0.0, 0.0, 1.0},
    {0, 0, 0, 0, 0, 0, 0.0, 1.0, 0.0},
    {0, 0, 0, 0, 0, 0, 0.25, 0.5, 0.25},
    {0, 0, 0, 0, 0, 0, 0.0, 0.5, 0.5},
    {0, 0, 0, 0, 0, 0, 0.0, 0.25, 0.75},
    {0, 0, 0, 0, 0, 0, 1.0, 0.0, 0.0}
  };
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 9; ++j) out(i, j) = values[i][j];
  out.attr("dimnames") = List::create(
    CharacterVector::create("U", "PO", "FS", "DEG2", "HA", "MZ"), R_NilValue);
  return out;
}

static inline double mixed_pair_probability(const std::vector<double>& states,
                                            const NumericMatrix& delta, int row) {
  double out = 0.0;
  for (int state = 0; state < 9; ++state) out += delta(row, state) * states[state];
  return out;
}

// [[Rcpp::export]]
NumericVector calc_pair_LL_cpp(List tables, IntegerMatrix G, int i, int j) {
  if (tables.size() < 1) stop("at least one compact locus table is required");
  if (i < 0 || j < 0 || i >= G.nrow() || j >= G.nrow()) stop("pair row is out of range");
  const NumericMatrix delta = legacy_delta();
  NumericVector out(delta.nrow());
  out.attr("names") = rownames(delta);
  for (int l = 0; l < tables.size(); ++l) {
    const List tab = tables[l];
    compact_validate_locus_table(tab);
    const int k = tab["k"];
    const int left = observed_genotype_index(G(i, 2 * l), G(i, 2 * l + 1), k);
    const int right = observed_genotype_index(G(j, 2 * l), G(j, 2 * l + 1), k);
    const std::vector<double> states = compact_pair_state_probabilities(tab, left, right);
    for (int row = 0; row < delta.nrow(); ++row)
      out[row] += log_probability_ln(mixed_pair_probability(states, delta, row),
                                     "pair likelihood probability");
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix score_basis_pairs_cpp(List tables, IntegerMatrix obsidx,
                                    IntegerVector i, IntegerVector j,
                                    NumericMatrix Delta) {
  if (i.size() != j.size()) stop("i and j must have equal length");
  if (Delta.ncol() != 9) stop("Delta must have exactly nine columns");
  if (obsidx.ncol() != tables.size()) stop("obsidx must have one column per locus");
  // Deduplicate labelled individual pairs first, then deduplicate the observed
  // descriptors and ordered observed-state pairs within each locus.  Ordered
  // keys are necessary because D3/D5 and D4/D6 exchange on reversal.
  std::map<std::pair<int, int>, int> individual_pair_index;
  std::vector<int> requested_to_unique(i.size());
  std::vector<std::pair<int, int> > unique_individual_pairs;
  for (int q = 0; q < i.size(); ++q)
    if (i[q] < 0 || j[q] < 0 || i[q] >= obsidx.nrow() || j[q] >= obsidx.nrow())
      stop("pair row is out of range");
  for (int q = 0; q < i.size(); ++q) {
    const std::pair<int, int> key(i[q], j[q]);
    std::map<std::pair<int, int>, int>::iterator found =
      individual_pair_index.find(key);
    if (found == individual_pair_index.end()) {
      const int index = unique_individual_pairs.size();
      individual_pair_index[key] = index;
      unique_individual_pairs.push_back(key);
      requested_to_unique[q] = index;
    } else {
      requested_to_unique[q] = found->second;
    }
  }
  NumericMatrix unique_out(unique_individual_pairs.size(), Delta.nrow());
  double unique_descriptors = 0.0, unique_state_pairs = 0.0;
  for (int l = 0; l < tables.size(); ++l) {
    const List tab = tables[l];
    compact_validate_locus_table(tab);
    const int k = tab["k"];
    const NumericVector p = tab["p"];
    const IntegerMatrix geno = tab["geno"];
    const CompactEmissionKernel kernel = compact_kernel_from_table(tab);
    std::map<int, CompactObservedDescriptor> descriptors;
    std::map<std::pair<int, int>, std::vector<double> > state_pairs;
    for (int q = 0; q < static_cast<int>(unique_individual_pairs.size()); ++q) {
      const int left = obsidx(unique_individual_pairs[q].first, l);
      const int right = obsidx(unique_individual_pairs[q].second, l);
      if (descriptors.find(left) == descriptors.end())
        descriptors[left] = compact_observed_descriptor_from_parts(
          k, p, geno, kernel, left);
      if (descriptors.find(right) == descriptors.end())
        descriptors[right] = compact_observed_descriptor_from_parts(
          k, p, geno, kernel, right);
      const std::pair<int, int> pair_key(left, right);
      std::map<std::pair<int, int>, std::vector<double> >::iterator pair_found =
        state_pairs.find(pair_key);
      if (pair_found == state_pairs.end()) {
        const CompactObservedDescriptor& x = descriptors.find(left)->second;
        const CompactObservedDescriptor& y = descriptors.find(right)->second;
        std::vector<double> states(9, 0.0);
        for (int a = 1; a <= k; ++a) {
          const double pa = p[a - 1];
          states[0] += pa * x.d[a] * y.d[a];
          states[2] += pa * x.d[a] * y.r[a];
          states[4] += pa * x.r[a] * y.d[a];
          states[7] += pa * x.r[a] * y.r[a];
        }
        states[1] = x.h * y.h;
        states[3] = x.h * y.m;
        states[5] = x.m * y.h;
        states[6] = compact_d7_probability(tab, left, right);
        states[8] = x.m * y.m;
        pair_found = state_pairs.insert(std::make_pair(pair_key, states)).first;
      }
      const std::vector<double>& states = pair_found->second;
      for (int row = 0; row < Delta.nrow(); ++row)
        unique_out(q, row) += log_probability_ln(
          mixed_pair_probability(states, Delta, row),
          "pair-basis likelihood probability");
    }
    unique_descriptors += descriptors.size();
    unique_state_pairs += state_pairs.size();
  }
  NumericMatrix out(i.size(), Delta.nrow());
  for (int q = 0; q < i.size(); ++q)
    for (int row = 0; row < Delta.nrow(); ++row)
      out(q, row) = unique_out(requested_to_unique[q], row);
  out.attr("compact_batch_counters") = List::create(
    _["requested_pairs"] = static_cast<double>(i.size()),
    _["unique_individual_pairs"] = static_cast<double>(unique_individual_pairs.size()),
    _["requested_pair_locus_cells"] =
      static_cast<double>(i.size()) * static_cast<double>(tables.size()),
    _["unique_observed_descriptors"] = unique_descriptors,
    _["unique_ordered_observed_state_pairs"] = unique_state_pairs,
    _["reused_pair_locus_cells"] =
      static_cast<double>(i.size()) * static_cast<double>(tables.size()) - unique_state_pairs,
    _["route"] = "compiled_batch");
  return out;
}

struct CompactTransmission {
  std::vector<double> u;
  std::vector<double> log_u;
  double m;
  double log_m;
};

static inline double compact_log_founder_probability(
    const NumericVector& p, int a, int b) {
  if (a == b) return 2.0 * std::log(p[a - 1]);
  return std::log(2.0) + std::log(p[a - 1]) + std::log(p[b - 1]);
}

static CompactTransmission compact_transmitted_mass(const List& tab, int observed) {
  const int k = tab["k"];
  const NumericVector p = tab["p"];
  CompactTransmission out;
  out.u.assign(k + 1, 0.0);
  out.log_u.assign(k + 1, R_NegInf);
  if (observed < 0) {
    for (int a = 1; a <= k; ++a) {
      out.u[a] = p[a - 1];
      out.log_u[a] = std::log(p[a - 1]);
    }
    out.m = 1.0;
    out.log_m = 0.0;
    return out;
  }
  const int G = tab["G"];
  if (observed >= G) stop("observed genotype index is out of range");
  const IntegerMatrix geno = tab["geno"];
  const CompactEmissionKernel kernel = compact_kernel_from_table(tab);
  std::vector<LogSumExpAccumulator> allele_log_mass(k + 1);
  for (int g = 0; g < G; ++g) {
    const double emission = compact_emission_probability(
      kernel, geno(g, 0), geno(g, 1), geno(observed, 0), geno(observed, 1));
    if (emission == 0.0) continue;
    const double log_mass = compact_log_founder_probability(
      p, geno(g, 0), geno(g, 1)) + std::log(emission);
    allele_log_mass[geno(g, 0)].add_log(log_mass - std::log(2.0),
                                        "transmitted allele probability");
    allele_log_mass[geno(g, 1)].add_log(log_mass - std::log(2.0),
                                        "transmitted allele probability");
  }
  LogSumExpAccumulator total_log_mass;
  for (int a = 1; a <= k; ++a) {
    out.log_u[a] = allele_log_mass[a].value();
    out.u[a] = is_negative_infinity_ln(out.log_u[a]) ? 0.0 :
      std::exp(out.log_u[a]);
    total_log_mass.add_log(out.log_u[a], "transmitted mass");
  }
  out.log_m = total_log_mass.value();
  out.m = is_negative_infinity_ln(out.log_m) ? 0.0 : std::exp(out.log_m);
  return out;
}

static inline double compact_offspring_emission(const List& tab,
                                                 const CompactEmissionKernel& kernel,
                                                 const IntegerMatrix& geno,
                                                 int observed, int a, int b) {
  return observed < 0 ? 1.0 : compact_emission_probability(
    kernel, a, b, geno(observed, 0), geno(observed, 1));
}

static inline double compact_log_dot_transmission(
    const CompactTransmission& x, const CompactObservedDescriptor& o, int k) {
  LogSumExpAccumulator out;
  for (int a = 1; a <= k; ++a) {
    if (is_negative_infinity_ln(x.log_u[a]) ||
        is_negative_infinity_ln(o.log_r[a])) continue;
    out.add_log(x.log_u[a] + o.log_r[a],
                "one-parent trio likelihood probability");
  }
  return out.value();
}

static inline double compact_log_two_parent_probability(
    const CompactTransmission& dam, const CompactTransmission& sire,
    const List& tab, const CompactEmissionKernel& kernel,
    const IntegerMatrix& geno, int observed, int k) {
  LogSumExpAccumulator out;
  for (int a = 1; a <= k; ++a) for (int b = 1; b <= k; ++b) {
    if (is_negative_infinity_ln(dam.log_u[a]) ||
        is_negative_infinity_ln(sire.log_u[b])) continue;
    const double emission = compact_offspring_emission(
      tab, kernel, geno, observed, a, b);
    if (emission == 0.0) continue;
    out.add_log(dam.log_u[a] + sire.log_u[b] + std::log(emission),
                "two-parent trio likelihood probability");
  }
  return out.value();
}

static NumericVector compact_trio_one(List tables, IntegerMatrix G, int o, int d, int s) {
  if (o < 0 || d < 0 || s < 0 || o >= G.nrow() || d >= G.nrow() || s >= G.nrow())
    stop("trio row is out of range");
  double ll_pp = 0.0, ll_d0 = 0.0, ll_0s = 0.0, ll_u = 0.0;
  for (int l = 0; l < tables.size(); ++l) {
    const List tab = tables[l];
    compact_validate_locus_table(tab);
    const int k = tab["k"];
    const int oo = observed_genotype_index(G(o, 2 * l), G(o, 2 * l + 1), k);
    const int od = observed_genotype_index(G(d, 2 * l), G(d, 2 * l + 1), k);
    const int os = observed_genotype_index(G(s, 2 * l), G(s, 2 * l + 1), k);
    const CompactObservedDescriptor child = compact_observed_descriptor(tab, oo);
    const CompactTransmission dam = compact_transmitted_mass(tab, od);
    const CompactTransmission sire = compact_transmitted_mass(tab, os);
    const IntegerMatrix geno = tab["geno"];
    const CompactEmissionKernel kernel = compact_kernel_from_table(tab);
    const double h_pp_ln = compact_log_two_parent_probability(
      dam, sire, tab, kernel, geno, oo, k);
    const double h_d0_ln = compact_log_dot_transmission(dam, child, k) +
      sire.log_m;
    const double h_0s_ln = compact_log_dot_transmission(sire, child, k) +
      dam.log_m;
    const double h_u_ln = child.log_m + dam.log_m + sire.log_m;
    ll_pp += h_pp_ln;
    ll_d0 += h_d0_ln;
    ll_0s += h_0s_ln;
    ll_u += h_u_ln;
  }
  return NumericVector::create(_["H_PP"] = ll_pp, _["H_D0"] = ll_d0,
                               _["H_0S"] = ll_0s, _["H_U"] = ll_u);
}

// [[Rcpp::export]]
NumericVector calc_trio_LL_cpp(List tables, IntegerMatrix G, int o, int d, int s) {
  return compact_trio_one(tables, G, o, d, s);
}


// Truncated Poisson-binomial count distribution per row, for the error-rate
// fit. Same recurrence and operation order as the former R loop; an uncalled
// locus with a finite probability leaves the row unchanged, so it is skipped.
// [[Rcpp::export]]
NumericMatrix error_count_pmf_cpp(NumericVector probabilities, LogicalMatrix called,
                                  int ceiling) {
  const int n = called.nrow(), L = called.ncol(), K = ceiling + 1;
  if (probabilities.size() != L || ceiling < 0) stop("error_count_pmf_cpp: bad dimensions");
  NumericMatrix out(n, K);
  std::vector<double> dp(K);
  for (int i = 0; i < n; ++i) {
    std::fill(dp.begin(), dp.end(), 0.0); dp[0] = 1.0;
    for (int l = 0; l < L; ++l) {
      const int c = called(i, l);
      if (c == NA_LOGICAL) stop("error_count_pmf_cpp: missing called flag");
      const double q = probabilities[l];
      if (!c && std::isfinite(q)) continue;
      const double p = q * (c ? 1.0 : 0.0), keep = 1 - p;
      for (int k = K - 1; k >= 1; --k) dp[k] = dp[k] * keep + dp[k - 1] * p;
      dp[0] = dp[0] * keep;
    }
    for (int k = 0; k < K; ++k) out(i, k) = dp[k];
  }
  return out;
}
