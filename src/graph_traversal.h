#pragma once
#include <Rcpp.h>
#include <deque>
#include <vector>
namespace pap_internal::search {
// Cache graph structure only; evidence, ages and reopening policies remain callers' decisions.
class RecheckGraph {
 int revision=-1;
 std::vector<std::vector<int>> children, adjacent;
 public:
 void synchronize(const Rcpp::IntegerMatrix& parents,int graph_revision);
 void spread_descendant_revisions(std::vector<int>& latest,std::deque<int> queue)const;
 void spread_distances(std::vector<int>& distance,std::deque<int> queue,int reach)const;
 std::vector<int> spread_neighborhood_max(std::vector<int> values,int radius)const;
};
}
