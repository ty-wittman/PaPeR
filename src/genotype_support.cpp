#include "locus_tables.h"
#include "latent_mendelian_support.h"
namespace pap_internal::search {


LatentMendelianChildSupport locus_mendelian_child_support(
    const LocusTables& locus, int dam_genotype, int sire_genotype) {
  return latent_mendelian_child_support(
    dam_genotype, sire_genotype,
    locus.genotype_first.data(), locus.genotype_second.data(),
    [&locus](int first, int second) {
      return genotype_index_zero(first, second, locus.k);
    });
}
}
