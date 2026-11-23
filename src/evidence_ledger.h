#ifndef PAP_EVIDENCE_LEDGER_H
#define PAP_EVIDENCE_LEDGER_H
#include "pedigree_graph_snapshot.h"
#include <deque>
#include <set>
#include <vector>
#include <map>
#include <numeric>
#include <array>

namespace pap_internal::pedigree_message {
struct EvidenceLedger {
 int focal=0,revision=0;
 std::set<int> observed,structural,cut_births,core;
 std::set<std::pair<int,int>> split_parents;
 int variant=0,components=0,max_candidates_per_component=0;
 std::vector<int> component_candidate_counts;
 static std::set<int> family_core(const PedigreeGraphSnapshot& graph,int focal,const std::vector<int>& candidates){
  std::set<int> roots(candidates.begin(),candidates.end()),out;roots.insert(focal);roots.erase(0);
  for(int root:roots){out.insert(root);for(int side:{1,2})if(int p=graph.parent(root,side))out.insert(p);
   for(int child:graph.children(root)){out.insert(child);for(int side:{1,2})if(int p=graph.parent(child,side))out.insert(p);}
   int dam=graph.parent(root,1),sire=graph.parent(root,2);
   if(dam&&sire)for(int sib:graph.children(dam))if(graph.parent(sib,2)==sire)out.insert(sib);
  }return out;
 }
 void finish(const PedigreeGraphSnapshot& graph){
  structural=observed;cut_births.clear();
  for(int node:observed)if(node!=focal)for(int side:{1,2})if(int p=graph.parent(node,side))structural.insert(p);
  for(int node:structural)if(!observed.count(node)&&(graph.parent(node,1)||graph.parent(node,2)))cut_births.insert(node);
 }
 void far_tree(const PedigreeGraphSnapshot& graph,int distance,bool protect_core_families=true){
  auto protected_nodes=core,layer=core;
  for(int d=0;d<distance;++d){std::set<int> next;for(int v:layer){for(int side:{1,2})if(int p=graph.parent(v,side))if(protected_nodes.insert(p).second)next.insert(p);
   for(int child:graph.children(v))if(child!=focal&&protected_nodes.insert(child).second)next.insert(child);}layer=std::move(next);}
  std::map<int,int> index;for(int v:structural)index[v]=index.size();int n=index.size();
  std::vector<int> uf(2*n);std::iota(uf.begin(),uf.end(),0);
  auto find=[&](int v){while(uf[v]!=v){uf[v]=uf[uf[v]];v=uf[v];}return v;};
  auto join=[&](int a,int b){a=find(a);b=find(b);if(a==b)return false;uf[a]=b;return true;};
  auto protected_birth=[&](int child){
   // One generation from a registered role includes the whole incident birth,
   // including its co-parent. Never sever a core parent from its own child
   // merely because that co-parent lies outside the individual-node radius.
   bool touches_core=core.count(child)!=0;for(int side:{1,2})touches_core|=core.count(graph.parent(child,side))!=0;
   if(protect_core_families&&touches_core)return true;
   if(!protected_nodes.count(child))return false;for(int side:{1,2})if(int parent=graph.parent(child,side))if(!protected_nodes.count(parent))return false;return true;};
  for(int pass=0;pass<2;++pass)for(int child:observed)if(child!=focal&&protected_birth(child)==(pass==0)){
   int factor=n+index.at(child);join(factor,index.at(child));
   for(int side:{1,2})if(int p=graph.parent(child,side))if(!join(factor,index.at(p))&&pass==1)split_parents.emplace(child,side);
  }
 }
 void candidate_components(const PedigreeGraphSnapshot& graph,const std::vector<int>& candidates){
  std::map<int,int> labels;int label=0;for(int root:structural)if(!labels.count(root)){std::deque<int> q{root};labels[root]=label;
   while(!q.empty()){int v=q.front();q.pop_front();auto add=[&](int w){if(w&&structural.count(w)&&!labels.count(w)){labels[w]=label;q.push_back(w);}};
    if(v!=focal&&observed.count(v))for(int side:{1,2})if(!split_parents.count({v,side}))add(graph.parent(v,side));
    for(int child:graph.children(v))if(child!=focal&&observed.count(child)){bool linked=false;for(int side:{1,2})linked|=graph.parent(child,side)==v&&!split_parents.count({child,side});if(linked)add(child);}
   }++label;}
  components=label;component_candidate_counts.assign(label,0);std::set<int> unique(candidates.begin(),candidates.end());
  for(int v:unique)if(v&&labels.count(v))++component_candidate_counts[labels.at(v)];
  for(int count:component_candidate_counts)max_candidates_per_component=std::max(max_candidates_per_component,count);
 }
 // Recursive own-family evidence. The same closure recipe applies to every
 // individual; the focal's incoming birth is excluded from every cavity.
 static EvidenceLedger build(const PedigreeGraphSnapshot& graph,int focal,
                             const std::vector<int>& candidates) {
  EvidenceLedger out;out.focal=focal;out.revision=graph.graph_revision();
  out.observed.insert(focal);
  for(int side:{1,2})if(graph.parent(focal,side))out.observed.insert(graph.parent(focal,side));
  for(int child:graph.children(focal)){out.observed.insert(child);for(int side:{1,2})if(graph.parent(child,side))out.observed.insert(graph.parent(child,side));}
  auto layer=out.observed;
  for(int depth=0;depth<2;++depth){std::set<int>next;for(int node:layer)for(int side:{1,2}){int p=graph.parent(node,side);if(p&&out.observed.insert(p).second)next.insert(p);}layer=std::move(next);}
  for(int root:candidates)if(root){std::deque<int>todo{root};std::set<int>seen{root};
   while(!todo.empty()){int node=todo.front();todo.pop_front();out.observed.insert(node);
    if(node!=focal)for(int side:{1,2}){int p=graph.parent(node,side);if(p&&seen.insert(p).second)todo.push_back(p);}
    for(int child:graph.children(node))if(child!=focal&&seen.insert(child).second)todo.push_back(child);
   }
  }
  out.structural=out.observed;
  for(int node:out.observed)if(node!=focal)for(int side:{1,2})if(graph.parent(node,side))out.structural.insert(graph.parent(node,side));
  for(int node:out.structural)if(!out.observed.count(node)&&
    (graph.parent(node,1)||graph.parent(node,2)))out.cut_births.insert(node);
  return out;
 }
};
}
#endif
