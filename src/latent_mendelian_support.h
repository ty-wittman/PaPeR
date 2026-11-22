#ifndef PAP_LATENT_MENDELIAN_SUPPORT_H
#define PAP_LATENT_MENDELIAN_SUPPORT_H

#include <algorithm>
#include <array>

struct LatentMendelianChildSupport {
  std::array<int, 4> genotype{{-1, -1, -1, -1}};
  int count = 0;
};

// Exact latent-child support for one fixed pair of latent true parent
// genotypes.  The allele arrays describe canonical genotype states; neither
// observed calls nor emission probabilities are inputs to this helper.
template <typename GenotypeIndex>
static inline LatentMendelianChildSupport latent_mendelian_child_support(
    int dam_genotype, int sire_genotype,
    const int* genotype_first, const int* genotype_second,
    GenotypeIndex genotype_index) {
  const int dam_first = genotype_first[dam_genotype];
  const int dam_second = genotype_second[dam_genotype];
  const int sire_first = genotype_first[sire_genotype];
  const int sire_second = genotype_second[sire_genotype];

  std::array<int, 4> gamete_paths{{
    genotype_index(dam_first, sire_first),
    genotype_index(dam_first, sire_second),
    genotype_index(dam_second, sire_first),
    genotype_index(dam_second, sire_second)
  }};
  std::sort(gamete_paths.begin(), gamete_paths.end());

  LatentMendelianChildSupport out;
  for (int genotype : gamete_paths) {
    if (out.count == 0 || out.genotype[out.count - 1] != genotype)
      out.genotype[out.count++] = genotype;
  }
  return out;
}

namespace pap_internal::search {
struct LocusTables;
LatentMendelianChildSupport locus_mendelian_child_support(
    const LocusTables& locus, int dam_genotype, int sire_genotype);
}

#endif
