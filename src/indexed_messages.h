#ifndef PAP_INDEXED_MESSAGES_H
#define PAP_INDEXED_MESSAGES_H
// Included inside pedigree_message, after TreeStructure. Symbolic ownership is
// compiled once per graph revision and shared across every locus. Numerical
// messages are vectors indexed by individual/family incidence, not string keys.
struct IndexedForestPlan {
 std::shared_ptr<const TreeStructure> tree;
 std::vector<std::pair<int,int>> incidence; // family, individual
 std::map<std::pair<int,int>,size_t> slots;
 std::vector<std::string> identities;
 std::vector<std::vector<uint32_t>> to_dependencies,from_dependencies,local;
 std::vector<std::vector<size_t>> incident;
 std::vector<uint32_t> axes;
 IndexedForestPlan(const relationship_topology::TopologyModel& m,
                   std::shared_ptr<const TreeStructure> source,Store& symbols)
 :tree(std::move(source)),local(m.nodes.size()),incident(m.nodes.size()),axes(m.nodes.size()){
  for(int v:tree->variables)axes[v]=symbols.symbol("Axis:"+identity(m.nodes[v]));
  std::vector<std::vector<uint32_t>> births(tree->families.size());
  for(const auto& a:tree->atoms)if(a.kind!=Kind::birth)local[a.child].push_back(symbols.symbol(a.key));
  for(size_t f=0;f<tree->families.size();++f){const auto& family=tree->families[f];
   for(size_t at:family.births)births[f].push_back(symbols.symbol(tree->atoms[at].key));
   for(int v:family.members){size_t slot=incidence.size();incidence.emplace_back(f,v);slots.emplace(std::make_pair(int(f),v),slot);incident[v].push_back(slot);
    identities.push_back(encoded(identity(m.nodes[v]))+encoded(identity(m.nodes[family.dam]))+encoded(identity(m.nodes[family.sire])));}
  }
  to_dependencies.resize(incidence.size());from_dependencies.resize(incidence.size());
  std::vector<unsigned char> ready_to(incidence.size()),ready_from(incidence.size());
  auto canonical=[](std::vector<uint32_t>& x){std::sort(x.begin(),x.end());x.erase(std::unique(x.begin(),x.end()),x.end());};
  auto build=[&](auto&& self,size_t slot,bool incoming)->const std::vector<uint32_t>&{
   auto& ready=incoming?ready_from:ready_to;auto& result=incoming?from_dependencies[slot]:to_dependencies[slot];
   if(ready[slot])return result;ready[slot]=1;auto [f,v]=incidence[slot];
   if(incoming){result=births[f];for(int other:tree->families[f].members)if(other!=v){const auto& d=self(self,slots.at({f,other}),false);result.insert(result.end(),d.begin(),d.end());}}
   else{result=local[v];for(size_t other:incident[v])if(other!=slot){const auto& d=self(self,other,true);result.insert(result.end(),d.begin(),d.end());}}
   canonical(result);return result;
  };
  for(size_t slot=0;slot<incidence.size();++slot){build(build,slot,false);build(build,slot,true);}
 }
 size_t bytes()const{size_t n=sizeof(*this)+incidence.size()*sizeof(incidence[0])+slots.size()*64;
  for(const auto& id:identities)n+=sizeof(id)+id.capacity();
  for(const auto* vectors:{&to_dependencies,&from_dependencies,&local}){n+=vectors->size()*sizeof(std::vector<uint32_t>);for(const auto& x:*vectors)n+=x.capacity()*sizeof(uint32_t);}
  for(const auto& x:incident)n+=sizeof(x)+x.capacity()*sizeof(size_t);return n+axes.capacity()*sizeof(uint32_t);}
};
struct IndexedForestMessages {
 std::shared_ptr<const IndexedForestPlan> plan;
 std::vector<Cached> to,from;
 std::vector<unsigned char> to_ready,from_ready;
 size_t computed=0,reused=0,restored=0;
 explicit IndexedForestMessages(std::shared_ptr<const IndexedForestPlan> p):plan(std::move(p)),
  to(plan->incidence.size()),from(to.size()),to_ready(to.size()),from_ready(to.size()){}
 // Restore only unchanged directed dependency sets. A changed birth invalidates
 // exactly the directions whose branches include it, including across merges.
 size_t slot(int family,int individual)const{return plan->slots.at({family,individual});}
 size_t cells()const{size_t n=0;for(size_t j=0;j<to.size();++j){if(to_ready[j])n+=to[j].value.size();if(from_ready[j])n+=from[j].value.size();}return n;}
 size_t bytes()const{size_t n=sizeof(*this)+(to.size()+from.size())*sizeof(Cached)+to_ready.size()+from_ready.size();for(const auto* values:{&to,&from})for(const auto& x:*values)n+=x.value.capacity()*sizeof(double);return n;}
};
struct IndexedForestProjection {
 size_t epoch=0,component=SIZE_MAX;int individual=-1;
 bool valid=false,phi=false,prior=false;
 std::vector<size_t> incoming;
};
#endif
