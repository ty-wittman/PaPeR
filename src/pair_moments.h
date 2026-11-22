#ifndef PAP_PAIR_MOMENTS_H
#define PAP_PAIR_MOMENTS_H
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>
namespace pap_internal::pair_moments {
struct Unsafe:std::runtime_error {Unsafe():std::runtime_error("factorized pair linear range unavailable") {}};
inline double product(double a,double b) {
 double v=a*b;
 if(!std::isfinite(v)||(a>0&&b>0&&(v==0||v<std::numeric_limits<double>::min())))throw Unsafe();
 return v;
}
struct Moments {double mean=0,within=0;std::vector<double> allele,diagonal,matrix;};
// The caller owns this cache in one immutable, locus-specific EffectivePhi.
inline const Moments& get(std::optional<Moments>& cache,const std::vector<double>& phi,
 const std::vector<double>& p,const std::vector<double>& prior,
 const std::vector<std::array<int,2>>& alleles) {
 if(cache)return *cache;
 Moments out;out.allele.resize(p.size());
 for(size_t g=0;g<phi.size();++g){
  out.mean+=product(prior[g],phi[g]);auto [a,b]=alleles[g];
  out.allele[a]+=product(p[b],phi[g]);
  if(a!=b)out.allele[b]+=product(p[a],phi[g]);
 }
 if(!std::isfinite(out.mean))throw Unsafe();
 for(double v:out.allele)if(!std::isfinite(v))throw Unsafe();
 cache=std::move(out);return *cache;
}
// Lazy extension for within-individual IBD; ordinary pair work is unchanged.
inline const Moments& get_inbred(std::optional<Moments>& cache,const std::vector<double>& phi,
 const std::vector<double>& p,const std::vector<double>& prior,
 const std::vector<std::array<int,2>>& alleles) {
 get(cache,phi,p,prior,alleles);auto& out=*cache;if(!out.diagonal.empty())return out;
 std::vector<double> diagonal(p.size());double within=0;
 for(size_t g=0;g<phi.size();++g){auto [a,b]=alleles[g];if(a==b){diagonal[a]=phi[g];within+=product(p[a],phi[g]);}}
 if(!std::isfinite(within))throw Unsafe();out.within=within;out.diagonal=std::move(diagonal);return out;
}
inline double contract(const std::array<double,3>& ibd,const Moments* a,const Moments* b,
 const std::vector<double>& phi_a,const std::vector<double>& phi_b,
 const std::vector<double>& p,const std::vector<double>& prior) {
 double sum=0;
 if(ibd[0]>0)sum+=product(ibd[0],product(a->mean,b->mean));
 if(ibd[1]>0){double one=0;for(size_t x=0;x<p.size();++x)one+=product(p[x],product(a->allele[x],b->allele[x]));sum+=product(ibd[1],one);}
 if(ibd[2]>0){double two=0;for(size_t g=0;g<prior.size();++g)two+=product(prior[g],product(phi_a[g],phi_b[g]));sum+=product(ibd[2],two);}
 if(!std::isfinite(sum))throw Unsafe();return sum;
}
}
#endif
