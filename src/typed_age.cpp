#include <Rcpp.h>
#include "typed_age.h"
using namespace Rcpp;
namespace age=pap_internal::typed_age;

// [[Rcpp::export]]
SEXP typed_age_state_cpp(IntegerMatrix parents,IntegerVector lo,IntegerVector hi,
                         IntegerVector last,NumericVector minimum,NumericVector maximum) {
  std::vector<std::array<int,2>> links(parents.nrow());
  for(int v=0;v<parents.nrow();++v)for(int s=0;s<2;++s)links[v][s]=parents(v,s)-1;
  if(minimum.size()!=2||maximum.size()!=2)stop("typed age limits need two sides");
  XPtr<age::State> ptr(new age::State(as<std::vector<int>>(lo),as<std::vector<int>>(hi),
    as<std::vector<int>>(last),{{minimum[0],minimum[1]}},{{maximum[0],maximum[1]}},links),true);
  return ptr;
}

// [[Rcpp::export]]
List typed_age_query_cpp(SEXP pointer,int child,IntegerVector parents) {
  XPtr<age::State> ptr(pointer);
  age::Patch patch;
  if(child>0){if(parents.size()!=2)stop("typed age proposal needs two parent slots");patch=ptr->propose(child-1,{{parents[0]-1,parents[1]-1}});}
  auto b=ptr->bounds(child>0&&patch.feasible?&patch:nullptr);
  auto nodes=patch.nodes;for(int& v:nodes)++v;
  return List::create(_["feasible"]=ptr->feasible&&patch.feasible,
    _["birth_lo"]=b.lo,_["birth_hi"]=b.hi,_["component"]=nodes,
    _["relaxations"]=double(patch.relaxations));
}

#include "pap_removed_path.h"

// Name of the first removed code path reached since the last call, or NA;
// reading it resets the record.
// [[Rcpp::export]]
CharacterVector pap_removed_path_hit_cpp() {
  const char* hit = pap_internal::removed_path_hit().exchange(nullptr);
  return hit == nullptr ? CharacterVector::create(NA_STRING) : CharacterVector::create(hit);
}
