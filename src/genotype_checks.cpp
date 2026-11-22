#include "genotype_checks.h"
#include "visit_graph.h"
namespace pap_internal::search {

void validate_locus_matrix(const NumericMatrix& matrix, int rows,
                                  int columns, const char* name) {
  if (matrix.nrow() != rows || matrix.ncol() != columns)
    engine_stop("simp snapshot %s has incompatible dimensions", name);
  for (int j = 0; j < columns; ++j) for (int i = 0; i < rows; ++i)
    if (!finite_or_neg_inf(matrix(i, j)))
      engine_stop("simp snapshot %s must be finite or -Inf", name);
}

int count_trio_me(const EngineState& compiled,
                         int child, int dam, int sire) {
  if (compiled.genotype_alleles.empty()) return -1;
  const int loci = static_cast<int>(compiled.genotype_alleles.size() /
                                    (static_cast<std::size_t>(compiled.n) * 2));
  int errors = 0;
  for (int locus = 0; locus < loci; ++locus) {
    const int ox = compiled.genotype_alleles[
      child - 1 + (2 * locus) * compiled.n];
    const int oy = compiled.genotype_alleles[
      child - 1 + (2 * locus + 1) * compiled.n];
    if (ox <= 0 || oy <= 0) continue;
    const int dx = compiled.genotype_alleles[dam - 1 + (2 * locus) * compiled.n];
    const int dy = compiled.genotype_alleles[dam - 1 + (2 * locus + 1) * compiled.n];
    const int sx = compiled.genotype_alleles[sire - 1 + (2 * locus) * compiled.n];
    const int sy = compiled.genotype_alleles[sire - 1 + (2 * locus + 1) * compiled.n];
    if (dx <= 0 || dy <= 0 || sx <= 0 || sy <= 0) continue;
    const bool dam_x = ox == dx || ox == dy;
    const bool dam_y = oy == dx || oy == dy;
    const bool sire_x = ox == sx || ox == sy;
    const bool sire_y = oy == sx || oy == sy;
    if (!((dam_x && sire_y) || (dam_y && sire_x))) ++errors;
  }
  return errors;
}

double trio_llr_ln(const EngineState& compiled, int child, int dam, int sire) {
  double total = 0.0;
  std::vector<double> from_dam, from_sire;
  for (const auto& L : compiled.locus) {
    const int k = L.k, G = L.G;
    // An absent parent (id 0) transmits population gametes.
    auto scale = [&](int id) {
      if (id <= 0) return 0.0;
      double top = R_NegInf;
      for (int g = 0; g < G; ++g) top = std::max(top, L.child_emission(id - 1, g));
      return top;
    };
    double child_top = R_NegInf;
    for (int a = 0; a < k; ++a) for (int b = 0; b < k; ++b) child_top = std::max(child_top, L.child_kernel(child - 1, a, b, 0));
    const double dam_top = scale(dam), sire_top = scale(sire);
    if (!R_FINITE(child_top) || !R_FINITE(dam_top) || !R_FINITE(sire_top)) continue;
    // Gamete distribution of each parent from its own observation and prior.
    auto gametes = [&](int id, double top, std::vector<double>& out) {
      out.assign(k, 0.0);
      double z = 0.0;
      for (int g = 0; g < G; ++g) {
        const double w = std::exp(L.prior_ln[g] + (id > 0 ? L.child_emission(id - 1, g) : 0.0) - top);
        z += w;
        for (int a = 0; a < k; ++a) out[a] += w * std::exp(L.transmission_ln[g + a * G]);
      }
      return z;
    };
    const double zd = gametes(dam, dam_top, from_dam), zs = gametes(sire, sire_top, from_sire);
    double trio = 0.0;
    for (int a = 0; a < k; ++a) for (int b = 0; b < k; ++b)
      trio += from_dam[a] * from_sire[b] * std::exp(L.child_kernel(child - 1, a, b, 0) - child_top);
    double alone = 0.0;
    for (int g = 0; g < G; ++g) alone += std::exp(L.prior_ln[g] + L.child_emission(child - 1, g) - child_top);
    if (!(trio > 0.0)) return R_NegInf;
    total += std::log(trio) - std::log(zd) - std::log(zs) - std::log(alone);
  }
  return total;
}

}
