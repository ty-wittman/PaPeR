// ============================================================================
// exclusion.cpp  --  Multiallelic error-tolerant OH exclusion filter
//
// Generalized opposing homozygote (OH) for parent-offspring: a candidate parent
// must transmit one allele to the offspring, so their genotypes must SHARE at
// least one allele. If the observed genotypes are DISJOINT (share no allele),
// that is a Mendelian mismatch. Direct comparison needs at most four equality
// tests and is independent of the numerical size of an allele identity code.
//
// Dropout/miscall make a TRUE PO pair carry a few false mismatches, so the
// threshold T is calibrated (see po_mismatch_prob_cpp + R calibrate_maxmismatch)
// to retain ~all true PO while excluding unrelated pairs.
// ============================================================================

#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include <vector>
#include "genotype_compat.h"
#include "compact_likelihood.h"
using namespace Rcpp;


// Fuse exact contiguous birth-year/Year.last feasibility with the symmetric
// unordered-pair OH screen.  Sex remains deliberately outside this frontend:
// the returned direction mask is biological only and is filtered after scoring.
// [[Rcpp::export]]
DataFrame screen_candidate_pairs_sparse_cpp(
    IntegerMatrix G, IntegerVector birth_lo, IntegerVector birth_hi,
    IntegerVector year_last, IntegerVector min_age, IntegerVector max_age,
    int T) {
  const int n = G.nrow();
  if (birth_lo.size() != n || birth_hi.size() != n || year_last.size() != n)
    stop("birth-year and Year.last vectors must match G rows");
  if (min_age.size() != 2 || max_age.size() != 2)
    stop("parent-age bounds must contain maternal and paternal values");
  for (int role = 0; role < 2; ++role) {
    if (IntegerVector::is_na(min_age[role]) ||
        IntegerVector::is_na(max_age[role]) || max_age[role] < min_age[role])
      stop("parent-age bounds must be finite ordered whole-number values");
  }

  const int nLoci = G.ncol() / 2;
  const std::size_t reserve = std::min<std::size_t>(
      static_cast<std::size_t>(n), static_cast<std::size_t>(4096));
  std::vector<int> keep_left, keep_right, keep_direction, keep_oh, keep_loci;
  keep_left.reserve(reserve);
  keep_right.reserve(reserve);
  keep_direction.reserve(reserve);
  keep_oh.reserve(reserve);
  keep_loci.reserve(reserve);

  auto known_birth = [&](int i) {
    return !IntegerVector::is_na(birth_lo[i]) &&
      !IntegerVector::is_na(birth_hi[i]);
  };
  auto direction_ok = [&](int parent, int child, int role) {
    // An unknown child support imposes no age or Year.last restriction, which
    // is the same NULL-support branch used by the old R direction_ok().
    if (!known_birth(child)) return true;

    const int child_lo = birth_lo[child];
    const int child_hi = birth_hi[child];
    const int last = year_last[parent];
    if (!known_birth(parent)) {
      // Age is unrestricted; a finite Year.last only needs one supported
      // child year no later than that cutoff.
      return IntegerVector::is_na(last) || child_lo <= last;
    }

    int effective_child_hi = child_hi;
    if (!IntegerVector::is_na(last))
      effective_child_hi = std::min(effective_child_hi, last);
    if (effective_child_hi < child_lo) return false;

    const long long minimum_difference =
      static_cast<long long>(child_lo) - birth_hi[parent];
    const long long maximum_difference =
      static_cast<long long>(effective_child_hi) - birth_lo[parent];
    return maximum_difference >= min_age[role] &&
      minimum_difference <= max_age[role];
  };

  const double theoretical_pairs =
    static_cast<double>(n) * static_cast<double>(n - 1) / 2.0;
  double hard_feasible_pairs = 0.0;
  double oh_total_loci_examined = 0.0;
  double oh_rejected_loci_examined_total = 0.0;
  double hard_feasibility_seconds = 0.0;
  double oh_seconds = 0.0;
  double oh_rejected_pairs = 0.0;

  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const std::chrono::steady_clock::time_point hard_started =
        std::chrono::steady_clock::now();
      int direction = 0;
      // bits: i -> j maternal/paternal 1/2; j -> i maternal/paternal 4/8.
      if (direction_ok(i, j, 0)) direction |= 1;
      if (direction_ok(i, j, 1)) direction |= 2;
      if (direction_ok(j, i, 0)) direction |= 4;
      if (direction_ok(j, i, 1)) direction |= 8;
      hard_feasibility_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - hard_started).count();
      if (direction == 0) continue;
      ++hard_feasible_pairs;

      const std::chrono::steady_clock::time_point oh_started =
        std::chrono::steady_clock::now();
      int oh = 0;
      int loci = 0;
      for (int locus = 0; locus < nLoci; ++locus) {
        ++loci;
        if (diploid_disjoint(G(i, 2 * locus), G(i, 2 * locus + 1),
                             G(j, 2 * locus), G(j, 2 * locus + 1)) &&
            ++oh > T)
          break;
      }
      oh_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - oh_started).count();
      oh_total_loci_examined += loci;
      if (oh <= T) {
        keep_left.push_back(i);
        keep_right.push_back(j);
        keep_direction.push_back(direction);
        keep_oh.push_back(oh);
        keep_loci.push_back(loci);
      } else {
        ++oh_rejected_pairs;
        oh_rejected_loci_examined_total += loci;
      }
    }
  }

  DataFrame out = DataFrame::create(
    _["left"] = wrap(keep_left), _["right"] = wrap(keep_right),
    _["direction_mask"] = wrap(keep_direction), _["oh_count"] = wrap(keep_oh),
    _["loci_examined"] = wrap(keep_loci));
  out.attr("candidate_unordered_pairs_total") = theoretical_pairs;
  out.attr("candidate_hard_feasible_pairs") = hard_feasible_pairs;
  out.attr("candidate_oh_survivors") = static_cast<double>(keep_left.size());
  out.attr("oh_pairs_examined") = hard_feasible_pairs;
  out.attr("oh_pairs_early_rejected") = oh_rejected_pairs;
  out.attr("oh_total_loci_examined") = oh_total_loci_examined;
  out.attr("oh_rejected_pairs") = oh_rejected_pairs;
  out.attr("oh_rejected_loci_examined_total") =
    oh_rejected_loci_examined_total;
  out.attr("candidate_hard_feasibility_seconds") = hard_feasibility_seconds;
  out.attr("candidate_oh_seconds") = oh_seconds;
  return out;
}


// ============================================================================
// po_mismatch_prob_cpp: per-locus P(OH | true PO pair, both typed).
// Returns scalars constructed with the compact locus table.  No observed-pair
// enumeration occurs at calibration time.
// ============================================================================
// [[Rcpp::export]]
DataFrame po_mismatch_prob_cpp(List tables) {
  int nLoci = tables.size();
  NumericVector pfm(nLoci), pum(nLoci), pme(nLoci);
  for (int l = 0; l < nLoci; ++l) {
    List tab = tables[l];
    compact_validate_locus_table(tab);
    pfm[l] = as<double>(tab["p_oh_true_po"]);
    pum[l] = as<double>(tab["p_oh_unrelated"]);
    pme[l] = as<double>(tab["p_me_true_trio"]);
    if (!R_finite(pfm[l]) || !R_finite(pum[l]) || !R_finite(pme[l]) ||
        pfm[l] < 0.0 || pfm[l] > 1.0 || pum[l] < 0.0 || pum[l] > 1.0 ||
        pme[l] < 0.0 || pme[l] > 1.0)
      stop("compact locus table contains invalid screening probabilities");
  }
  return DataFrame::create(_["locus"] = seq_len(nLoci),
                           _["p_fm"] = pfm, _["p_um"] = pum,
                           _["p_me_true_trio"] = pme);
}
