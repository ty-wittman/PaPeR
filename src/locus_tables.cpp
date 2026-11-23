#include "locus_tables.h"
namespace pap_internal::search {


auto LocusTables::child_kernel(int child, int first, int second, int n) const -> double  {
    return child_kernel_ln[(child * k + first) * k + second];
  }

auto LocusTables::child_emission(int child, int genotype) const -> double  {
    return child_kernel(child, genotype_first[genotype],
                        genotype_second[genotype], 0);
  }
}
