#pragma once
#include "locus_tables.h"
namespace pap_internal::search {
// Unavailable partitions affect visit fallback and public availability.
struct EngineCounters {
  double message_partition_unavailable = 0.0;
};

static constexpr int relationship_family_factored_max_k = 9;
static constexpr int relationship_family_factored_max_k2 =
  relationship_family_factored_max_k * relationship_family_factored_max_k;

enum class SideFamilyPartitionRoute : unsigned char {
  dense = 0,
  factored_child1 = 1,
  factored_child2 = 2
};

}
