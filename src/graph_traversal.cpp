#include "graph_traversal.h"
#include <algorithm>
namespace pap_internal::search {
void RecheckGraph::synchronize(const Rcpp::IntegerMatrix& parents,int graph_revision) {
 if(revision==graph_revision)return;
 children.assign(parents.nrow()+1,{});adjacent.assign(parents.nrow()+1,{});
 for(int child=1;child<=parents.nrow();++child)
  for(int side=0;side<2;++side)if(int parent=parents(child-1,side)) {
   children[parent].push_back(child);
   adjacent[child].push_back(parent);adjacent[parent].push_back(child);
  }
 revision=graph_revision;
}
void RecheckGraph::spread_descendant_revisions(std::vector<int>& latest,std::deque<int> queue)const {
 while(!queue.empty()) {
  int node=queue.front();queue.pop_front();
  for(int child:children[node])if(latest[child]<latest[node]) {
   latest[child]=latest[node];queue.push_back(child);
  }
 }
}
void RecheckGraph::spread_distances(std::vector<int>& distance,std::deque<int> queue,int reach)const {
 while(!queue.empty()) {
  int node=queue.front();queue.pop_front();
  if(distance[node]>=reach)continue;
  for(int next:adjacent[node])if(distance[next]<0) {
   distance[next]=distance[node]+1;queue.push_back(next);
  }
 }
}
std::vector<int> RecheckGraph::spread_neighborhood_max(std::vector<int> values,int radius)const {
 for(int step=0;step<radius;++step) {
  auto next=values;
  for(std::size_t node=1;node<adjacent.size();++node)
   for(int neighbor:adjacent[node])next[node]=std::max(next[node],values[neighbor]);
  values.swap(next);
 }
 return values;
}
}
