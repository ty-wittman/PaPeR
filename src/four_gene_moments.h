#ifndef PAP_FOUR_GENE_MOMENTS_H
#define PAP_FOUR_GENE_MOMENTS_H
#include "four_gene_identity.h"
#include "pair_moments.h"
namespace pap_internal::four_gene {
// Stored recipe type retained by the scorer's cache layout.
struct Recipe {Identity identity;std::array<double,9> weights{};};
// Distinguish declared structural zeros from exponentiation underflow before
// consulting any linear moment cache. Subnormal positives use the log route.
inline bool representable(const std::vector<double>& logs,const std::vector<double>& values){
 if(logs.size()!=values.size())return false;
 for(size_t i=0;i<logs.size();++i){if(logs[i]==-std::numeric_limits<double>::infinity()){if(values[i]!=0)return false;}
  else if(!std::isfinite(logs[i])||!std::isfinite(values[i])||values[i]<std::numeric_limits<double>::min())return false;}
 return true;
}

}
#endif
