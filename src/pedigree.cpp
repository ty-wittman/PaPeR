// Typed-parent relationship screening and diagnostics.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include "genotype_compat.h"
#include "probability_log.h"
#include "compact_likelihood.h"
using namespace Rcpp;

static inline double compact_mixed_pair_probability(
    const std::vector<double>& states, const NumericMatrix& delta, int row) {
  double out = 0.0;
  for (int state = 0; state < 9; ++state) out += delta(row, state) * states[state];
  return out;
}

// Score a bounded all-pairs typed relationship screen. Delta contains the
// canonical nine-state Jacquard rows selected by R; no relationship name enters
// compiled code. Missing observations use exact all-ones-emission marginals.
// [[Rcpp::export]]
DataFrame screen_rel_pairs_cpp(List tables, IntegerMatrix G,
                               NumericMatrix Delta, IntegerVector focal,
                               IntegerVector alts, double thresh,
                               IntegerVector byear, int max_agediff,
                               int min_agediff, int oh_max = -1) {
  int n = G.nrow(), loci = tables.size();
  int nf = focal.size(), na = alts.size();
  if (Delta.ncol() != 9 || nf == 0 || na == 0)
    stop("validated Delta, focal, and alternative indices are required");
  for (int i = 0; i < nf; ++i)
    if (focal[i] < 0 || focal[i] >= Delta.nrow()) stop("focal index is out of range");
  for (int i = 0; i < na; ++i)
    if (alts[i] < 0 || alts[i] >= Delta.nrow()) stop("alternative index is out of range");

  for (int l = 0; l < loci; ++l) compact_validate_locus_table(tables[l]);

  bool use_age = max_agediff >= 0, use_oh = oh_max >= 0;
  std::vector<int> I, J, fw, aw; std::vector<double> llr_out;
  for (int i = 0; i < n - 1; ++i) for (int j = i + 1; j < n; ++j) {
    if (use_age && byear.size() == n && byear[i] != NA_INTEGER &&
        byear[j] != NA_INTEGER) {
      int gap = std::abs(byear[i] - byear[j]);
      if (gap > max_agediff || gap < min_agediff) continue;
    }
    if (use_oh) {
      int oh = 0; bool reject = false;
      for (int l = 0; l < loci; ++l) if (diploid_disjoint(
          G(i, 2 * l), G(i, 2 * l + 1), G(j, 2 * l), G(j, 2 * l + 1)))
        if (++oh > oh_max) { reject = true; break; }
      if (reject) continue;
    }
    std::vector<double> f(nf, 0.0), a(na, 0.0);
    for (int l = 0; l < loci; ++l) {
      List tab = tables[l];
      int k = tab["k"];
      int oi = (G(i, 2 * l) <= 0 || G(i, 2 * l + 1) <= 0) ? -1 :
        compact_geno_index(G(i, 2 * l), G(i, 2 * l + 1), k);
      int oj = (G(j, 2 * l) <= 0 || G(j, 2 * l + 1) <= 0) ? -1 :
        compact_geno_index(G(j, 2 * l), G(j, 2 * l + 1), k);
      const std::vector<double> states = compact_pair_state_probabilities(tab, oi, oj);
      for (int q = 0; q < nf; ++q) f[q] += log_probability_ln(
        compact_mixed_pair_probability(states, Delta, focal[q]),
        "relationship-screen likelihood probability");
      for (int q = 0; q < na; ++q) a[q] += log_probability_ln(
        compact_mixed_pair_probability(states, Delta, alts[q]),
        "relationship-screen likelihood probability");
    }
    int best_f = 0, best_a = 0;
    for (int q = 1; q < nf; ++q) if (f[q] > f[best_f]) best_f = q;
    for (int q = 1; q < na; ++q) if (a[q] > a[best_a]) best_a = q;
    double llr = log_likelihood_contrast_ln(
      f[best_f], *std::max_element(a.begin(), a.end()),
      "relationship-screen likelihood contrast");
    if (llr > thresh) {
      I.push_back(i + 1); J.push_back(j + 1); llr_out.push_back(llr);
      fw.push_back(focal[best_f] + 1); aw.push_back(alts[best_a] + 1);
    }
  }
  return DataFrame::create(_["i"] = I, _["j"] = J, _["llr"] = llr_out,
                           _["focal_win"] = fw, _["alternative_win"] = aw);
}

static int count_mendelian_errors_one(
    const IntegerMatrix& G, int o, int d, int s) {
  const int loci = G.ncol() / 2;
  int me = 0;
  for (int l = 0; l < loci; ++l) {
    int ox = G(o, 2 * l), oy = G(o, 2 * l + 1);
    if (ox <= 0 || oy <= 0) continue;
    int dx = d >= 0 ? G(d, 2 * l) : 0, dy = d >= 0 ? G(d, 2 * l + 1) : 0;
    int sx = s >= 0 ? G(s, 2 * l) : 0, sy = s >= 0 ? G(s, 2 * l + 1) : 0;
    bool dok = d >= 0 && typed_diploid(dx, dy);
    bool sok = s >= 0 && typed_diploid(sx, sy);
    if (dok && sok) {
      bool dam_x = ox == dx || ox == dy, dam_y = oy == dx || oy == dy;
      bool sire_x = ox == sx || ox == sy, sire_y = oy == sx || oy == sy;
      if (!((dam_x && sire_y) || (dam_y && sire_x))) ++me;
    }
  }
  return me;
}

// Count typed OH and Mendelian incompatibilities for a final typed parent pair.
// Rows are zero-based; -1 denotes an unassigned parent.
// [[Rcpp::export]]
IntegerVector count_oh_me_trio_cpp(IntegerMatrix G, int o, int d, int s) {
  int loci = G.ncol() / 2, ohd = 0, ohs = 0;
  for (int l = 0; l < loci; ++l) {
    int ox = G(o, 2 * l), oy = G(o, 2 * l + 1);
    if (ox <= 0 || oy <= 0) continue;
    int dx = d >= 0 ? G(d, 2 * l) : 0, dy = d >= 0 ? G(d, 2 * l + 1) : 0;
    int sx = s >= 0 ? G(s, 2 * l) : 0, sy = s >= 0 ? G(s, 2 * l + 1) : 0;
    bool dok = d >= 0 && typed_diploid(dx, dy);
    bool sok = s >= 0 && typed_diploid(sx, sy);
    if (dok && diploid_disjoint(ox, oy, dx, dy)) ++ohd;
    if (sok && diploid_disjoint(ox, oy, sx, sy)) ++ohs;
  }
  return IntegerVector::create(_["OHdam"] = ohd, _["OHsire"] = ohs,
                               _["MEpair"] = count_mendelian_errors_one(G, o, d, s));
}

// Exact bounded batch for offspring x candidate x one known co-parent. Missing
// calls reduce the number of counted incompatibilities because a trio error is
// evaluated only when all three calls are typed.
// [[Rcpp::export]]
IntegerVector count_mendelian_errors_batch_cpp(IntegerMatrix G, int o,
                                               IntegerVector candidates,
                                               int co_parent) {
  if (o < 0 || o >= G.nrow()) stop("offspring index is out of range");
  if (co_parent < 0 || co_parent >= G.nrow())
    stop("co-parent index is out of range");
  IntegerVector out(candidates.size());
  for (int i = 0; i < candidates.size(); ++i) {
    const int candidate = candidates[i];
    if (candidate < 0 || candidate >= G.nrow())
      stop("candidate index is out of range");
    out[i] = count_mendelian_errors_one(G, o, candidate, co_parent);
  }
  return out;
}
