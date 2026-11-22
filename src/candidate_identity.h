#ifndef PAP_CANDIDATE_IDENTITY_H
#define PAP_CANDIDATE_IDENTITY_H
#include "message_graph.h"
#include <unordered_map>
#include <unordered_set>
namespace pap_internal::candidate_identity {
struct VectorHash {
 size_t operator()(const std::vector<int>& x)const noexcept {
  size_t h=x.size();for(int v:x)h^=std::hash<int>{}(v)+0x9e3779b9U+(h<<6)+(h>>2);return h;
 }
};
using AtomKey=std::array<int,4>;
struct AtomHash {size_t operator()(const AtomKey& a)const noexcept {
 size_t h=0;for(int v:a)h^=std::hash<int>{}(v)+0x9e3779b9U+(h<<6)+(h>>2);return h;
}};
inline AtomKey atom(const pedigree_message::Atom& a,const std::vector<int>& ids){
 return {{int(a.kind),ids.at(a.child),a.kind==pedigree_message::Kind::birth?ids.at(a.dam):-1,
  a.kind==pedigree_message::Kind::birth?ids.at(a.sire):-1}};
}
struct CatalogIdentities {
 int next;std::unordered_map<std::string,int> latent;
 explicit CatalogIdentities(int typed):next(typed){}
 std::vector<int> nodes(const ::relationship_topology::TopologyModel& m){
  std::vector<int> ids;for(const auto& n:m.nodes){if(n.typed_id>0)ids.push_back(n.typed_id-1);
   else {auto [at,fresh]=latent.emplace(n.label,next);if(fresh)++next;ids.push_back(at->second);}}
  return ids;
 }
};
// A source-missing identity is certified by incidence, not a latent label.
// It must be an unobserved founder used by exactly one actual source birth.
}
#endif
