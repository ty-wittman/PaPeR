#pragma once
#include "evidence_views.h"
namespace pap_internal::search {
struct LocusTables {
  int k = 0;
  int G = 0;
  std::vector<double> prior_ln;
  std::vector<double> log_p;
  // Rcpp matrices are column-major: genotype + allele * G.
  std::vector<double> transmission_ln;
  std::vector<int> child_index;
  // Genotype-state allele pairs in canonical one-based R order, stored as
  // zero-based C++ indices.  The compiled snapshot does not retain the
  // redundant G x n observation matrix after constructing child kernels.
  std::vector<int> genotype_first;
  std::vector<int> genotype_second;
  std::vector<double> message_ln;
  // k x k observation contractions for every child, stored child-major.
  std::vector<double> child_kernel_ln;
  // k-length transmitted masses for every published individual message.
  std::vector<double> transmitted_message_ln;

  double child_kernel(int child, int first, int second, int n) const;
  double child_emission(int child, int genotype) const;
};
}
