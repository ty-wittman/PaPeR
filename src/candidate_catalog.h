#include "pap_removed_path.h"
#ifndef PAP_CANDIDATE_CATALOG_H
#define PAP_CANDIDATE_CATALOG_H
#include <Rcpp.h>
#include "input_checks.h"
#include "typed_age.h"
#include "relationship_topology.h"
#include "relationship_catalog.h"
#include "score_types.h"
#include "predictive_models.h"
#include "score_types.h"
#include "evidence_ledger.h"
#include "candidate_identity.h"
#include <map>
#include <memory>
namespace pap_internal::candidate_catalog {
namespace t=::relationship_topology;
namespace pm=pedigree_message;
namespace id=candidate_identity;
struct Model {t::TopologyModel topology;const t::TopologyModel* finished=nullptr;std::string name_override;const t::TopologyModel& source()const{return finished?*finished:topology;}std::set<int> core;size_t invariant_factors=0;bool shared_context=false;};
struct Request {int model,locus,kind,focal,candidate,side,dam,sire,term,hypothesis;};
inline std::vector<Request> requests(Rcpp::IntegerMatrix events){
 if(events.ncol()!=10)throw std::invalid_argument("M_ped events require ten integer columns");
 std::vector<Request> out;out.reserve(events.nrow());
 for(int i=0;i<events.nrow();++i)out.push_back({events(i,0)-1,events(i,1)-1,events(i,2),events(i,3),events(i,4),events(i,5),events(i,6),events(i,7),events(i,8),events(i,9)});
 return out;
}
struct Variant{int definition;bool reverse;unsigned bits;int order;};
// A reverse PO is an alternative to the tested candidate-parent edge.
// Keep its codes outside the fixed ordinary25 and registered compound catalog.
inline const t::full_catalog::Definition& reverse_po_definition(){
 static const t::full_catalog::Definition definition{
  "PO",false,true,true,{{"A",{{"B","Q"}}}}};
 return definition;
}
// Bind only ancestor roles visited by registered motifs. Parent-slot variants
// require both source sides; unbound roles remain latent in the motif builder.
template<class Parent> inline std::set<int> registered_ancestors_by(Parent source_parent,const std::set<int>& endpoints){ pap_internal::removed_path("registered_ancestors_by"); }
inline std::set<int> registered_ancestors(const PedigreeGraphSnapshot& graph,const std::set<int>& endpoints){return registered_ancestors_by([&](int child,int side){return graph.parent(child,side);},endpoints);}
inline std::set<int> registered_ancestors(const std::vector<std::array<int,2>>& source,const std::set<int>& endpoints){return registered_ancestors_by([&](int child,int side){return source.at(child-1)[side-1]+1;},endpoints);}
inline pm::EvidenceLedger make_ledger(const PedigreeGraphSnapshot& graph,int focal,const std::vector<int>& candidates,Rcpp::List input){
 auto flag=[&](const char* name){return input.containsElementNamed(name)&&Rcpp::as<bool>(input[name]);};
 int bounded=flag("pedigree_ledger_bounded"),d1=flag("pedigree_ledger_far_tree_d1"),d2=flag("pedigree_ledger_far_tree_d2");
 if(bounded+d1+d2>1)throw std::invalid_argument("ledger experiment variants are mutually exclusive");
 const int variant=bounded?1:d1?2:d2?3:0;
 const bool legacy_node_radius=flag("pedigree_ledger_far_tree_legacy_node_radius");
 std::set<std::pair<int,int>> retain_links;
 if(input.containsElementNamed("pedigree_ledger_far_tree_keep")){
  Rcpp::IntegerMatrix links=input["pedigree_ledger_far_tree_keep"];
  if(links.ncol()!=2||(!d1&&!d2))throw std::invalid_argument("selective far-tree links require C1/C2 and child/side columns");
  for(int row=0;row<links.nrow();++row){int child=links(row,0),side=links(row,1);
   if(child<1||child>graph.node_count()||side<1||side>2||!graph.parent(child,side))throw std::invalid_argument("invalid retained far-tree parent link");
   retain_links.emplace(child,side);
  }
 }
 // The memo and scorer request the same immutable focal context. Construct it
 // once, then copy its small structural census; no per-hypothesis closure.
 struct Last {int focal=-1,variant=-1;std::vector<int> candidates;std::vector<std::array<int,2>> parents;pm::EvidenceLedger ledger;std::set<std::pair<int,int>> retain_links;bool legacy_node_radius=false;};
 static thread_local Last last;
 std::vector<std::array<int,2>> parents;for(int v=1;v<=graph.node_count();++v)parents.push_back({graph.parent(v,1),graph.parent(v,2)});
 auto ordered=candidates;std::sort(ordered.begin(),ordered.end());ordered.erase(std::unique(ordered.begin(),ordered.end()),ordered.end());
 if(last.focal==focal&&last.variant==variant&&last.candidates==ordered&&last.parents==parents&&last.retain_links==retain_links&&last.legacy_node_radius==legacy_node_radius)return last.ledger;
 pm::EvidenceLedger out;if(!bounded)out=pm::EvidenceLedger::build(graph,focal,candidates);else{out.focal=focal;out.revision=graph.graph_revision();}out.variant=variant;
 if(out.variant){out.core=pm::EvidenceLedger::family_core(graph,focal,candidates);std::set<int> endpoints(candidates.begin(),candidates.end());endpoints.insert(focal);endpoints.erase(0);
  auto ancestors=registered_ancestors(graph,endpoints);if(!bounded)out.core=endpoints;out.core.insert(ancestors.begin(),ancestors.end());
  if(bounded){out.observed=out.core;out.finish(graph);}else {out.far_tree(graph,d1?1:2,!legacy_node_radius);for(const auto& link:retain_links)out.split_parents.erase(link);}
 }
 out.candidate_components(graph,candidates);last={focal,variant,std::move(ordered),std::move(parents),out,retain_links,legacy_node_radius};return out;
}

struct ImmutableInput {
 t::Life life;std::vector<int> doff,soff,di,si;
 std::vector<t::TopologyModel> simple=t::topology_detail::specifications();
 std::vector<t::full_catalog::Definition> full=t::full_catalog::definitions(false,true);
 std::map<int,Variant> variants;Rcpp::List primitives;
 static constexpr const char* names[]={"sex","birth_lo","birth_hi","year_last","min_age","max_age","dam_offsets","sire_offsets","dam_ids","sire_ids"};
 explicit ImmutableInput(Rcpp::List input):primitives(10){
  for(int i=0;i<10;++i)primitives[i]=input[names[i]];
  life.sex=Rcpp::as<std::vector<int>>(input["sex"]);
  life.birth_lo=Rcpp::as<std::vector<int>>(input["birth_lo"]);life.birth_hi=Rcpp::as<std::vector<int>>(input["birth_hi"]);life.year_last=Rcpp::as<std::vector<int>>(input["year_last"]);
  auto lo=Rcpp::as<std::vector<double>>(input["min_age"]),hi=Rcpp::as<std::vector<double>>(input["max_age"]);for(int s=0;s<2;++s){life.min_age[s]=lo.at(s);life.max_age[s]=hi.at(s);}
  doff=Rcpp::as<std::vector<int>>(input["dam_offsets"]);soff=Rcpp::as<std::vector<int>>(input["sire_offsets"]);di=Rcpp::as<std::vector<int>>(input["dam_ids"]);si=Rcpp::as<std::vector<int>>(input["sire_ids"]);
  int order=28;for(size_t v=0;v<full.size();++v)for(int reverse=0;reverse<=(full[v].reverse_allowed?1:0);++reverse)for(unsigned bits=0;bits<(1U<<full[v].births.size());++bits)variants.emplace(order,Variant{int(v),bool(reverse),bits,order}),++order;
 }
 bool matches(Rcpp::List input)const{for(int i=0;i<10;++i){SEXP now=input[names[i]],old=primitives[i];if(now!=old&&!R_compute_identical(now,old,0))return false;}return true;}
 static std::shared_ptr<const ImmutableInput> parse(Rcpp::List input){return std::make_shared<const ImmutableInput>(input);}
};
struct Visit {
 bool parentless_hspo=false;
 bool reuse_source=false;size_t source_builds=0,source_copies=0,raw_null_reuses=0;
 int focal;Rcpp::IntegerMatrix parents;
 std::shared_ptr<const ImmutableInput> immutable;
 const t::Life& life;typed_age::Bounds typed_bounds;std::vector<int> candidates,dam_candidates,sire_candidates;
 std::unique_ptr<PedigreeGraphSnapshot> snapshot;pm::EvidenceLedger ledger;
 std::vector<std::array<int,2>> source;
 const std::vector<t::TopologyModel>& simple;const std::vector<t::full_catalog::Definition>& full;const std::map<int,Variant>& variants;
 size_t snapshots=0,lifetimes=0,ledgers=0,contexts=0,models_built=0;
 // pedigree likelihood: questions use only local individuals; the rest enters as stored messages.
 bool local_inputs=false;
 bool local_visit_facts=false;
 std::vector<std::array<int,2>> local_role_counts;

 explicit Visit(Rcpp::List input,int f,std::shared_ptr<const ImmutableInput> shared={}):focal(f),parents(input["parents"]),immutable(shared?std::move(shared):ImmutableInput::parse(input)),life(immutable->life),simple(immutable->simple),full(immutable->full),variants(immutable->variants){
  parentless_hspo=input.containsElementNamed("pedigree_hypothesis_suite")&&Rcpp::as<std::string>(input["pedigree_hypothesis_suite"])=="S37";
  local_visit_facts=input.containsElementNamed("likelihood_identity_visit_facts")&&Rcpp::as<bool>(input["likelihood_identity_visit_facts"]);
  if(local_visit_facts)local_role_counts.resize(parents.nrow()+1);


  ++lifetimes;const auto& doff=immutable->doff;const auto& soff=immutable->soff;const auto& di=immutable->di;const auto& si=immutable->si;
  for(int i=doff[focal-1];i<doff[focal];++i)dam_candidates.push_back(di[i]);
  for(int i=soff[focal-1];i<soff[focal];++i)sire_candidates.push_back(si[i]);
  for(int s=0;s<2;++s)if(parents(focal-1,s))(s?sire_candidates:dam_candidates).push_back(parents(focal-1,s));
  candidates=dam_candidates;candidates.insert(candidates.end(),sire_candidates.begin(),sire_candidates.end());
  int n=parents.nrow();std::vector<int>d(n+1),s(n+1),e(n+1,-1);std::vector<std::vector<int>>children(n+1);source.resize(n);
  for(int v=1;v<=n;++v){d[v]=parents(v-1,0);s[v]=parents(v-1,1);e[v]=v-1;source[v-1]={{d[v]-1,s[v]-1}};
   for(int p:{d[v],s[v]})if(p)children[p].push_back(v);
   if(local_visit_facts){if(d[v])++local_role_counts[d[v]][0];if(s[v])++local_role_counts[s[v]][1];}}
  snapshot=std::make_unique<PedigreeGraphSnapshot>(n,d,s,children,e,0);++snapshots;
  // pedigree likelihood questions read local inputs and stored forest messages, never this
  // ledger; its whole-component closure costs O(N) per candidate per visit.
  if(!(input.containsElementNamed("pedigree_likelihood")&&Rcpp::as<bool>(input["pedigree_likelihood"]))){ledger=make_ledger(*snapshot,focal,candidates,input);++ledgers;}
  if(typed_age::active)typed_bounds=*typed_age::active;
  else {
   auto links=source;links[focal-1]={{-1,-1}};
   typed_age::State ages(life.birth_lo,life.birth_hi,life.year_last,life.min_age,life.max_age,links);
   if(!ages.feasible)throw std::runtime_error("typed age source is infeasible");
   typed_bounds=ages.bounds();
  }

 }
 t::Life lifetime_view()const {
  t::Life out;out.borrowed_sex=&life.sex;out.borrowed_birth_lo=&life.birth_lo;
  out.borrowed_birth_lo=&typed_bounds.lo;out.borrowed_birth_hi=&typed_bounds.hi;out.borrowed_year_last=&life.year_last;
  out.min_age=life.min_age;out.max_age=life.max_age;return out;
 }
 // pedigree likelihood input: the focal, the tested candidate and opposite parent with
 // their parents and grandparents (raw: the focal and every candidate). A
 // birth is kept only when all its known parents are present, so every model
 // birth is either a real birth among these individuals or the hypothesis.
 template<class Parent> t::Input local_input(int candidate,int side,int opposite,t::Life& context_life,Parent&& parent){
  std::set<int> ids{focal};
  if(candidate||opposite){
   std::set<int> layer{focal};if(candidate)layer.insert(candidate);if(opposite)layer.insert(opposite);
   ids.insert(layer.begin(),layer.end());
   for(int depth=0;depth<2;++depth){std::set<int> next;for(int v:layer)for(int s=0;s<2;++s)if(int p=parent(v,s)){if(ids.insert(p).second)next.insert(p);}layer=std::move(next);}
  }else ids.insert(candidates.begin(),candidates.end());
  ids.erase(0);
  t::Input out;out.max_typed_nodes=4096;out.max_nodes=8192;out.optional_source_birth_budget=0;out.life=&context_life;
  // The local map has at most the protected endpoint ancestry, independent
  // of pedigree size. The reference path retains its old N-sized allocation.
  std::vector<int> wide_index;if(!local_visit_facts)wide_index.assign(parents.nrow()+1,-1);
  std::map<int,int> small_index;
  auto index=[&](int v)->int&{return local_visit_facts?
    small_index.try_emplace(v,-1).first->second:wide_index[v];};
  for(int v:ids){index(v)=out.nodes.size();t::TopologyNode node;node.typed_id=v;node.observed=true;node.label="T"+std::to_string(v);out.nodes.push_back(node);out.retained_source_birth_ids.push_back(v);}
  for(int v:ids){auto& node=out.nodes[index(v)];bool complete=true,any=false;
   for(int s=0;s<2;++s)if(int p=parent(v,s)){any=true;if(index(p)<0)complete=false;}
   if(any&&complete){for(int s=0;s<2;++s)if(int p=parent(v,s))node.parents[s]=index(p);node.has_birth=true;}}
  std::vector<int> roles;
  if(!local_visit_facts){roles.resize(parents.nrow()+1);
    for(int v=1;v<=parents.nrow();++v)for(int s=0;s<2;++s)if(int p=parent(v,s))roles[p]|=1<<s;
  }
  auto role=[&](int id){
    if(!local_visit_facts)return roles[id];
    auto count=local_role_counts[id];
    // Only the focal and candidate births differ in input()'s cavity. Repair
    // role counts by those incidences, retaining roles from other offspring.
    for(int which=0;which<2;++which){int v=which?candidate:focal;
      if(v<=0||(which&&v==focal))continue;
      for(int s=0;s<2;++s){if(parents(v-1,s)==id)--count[s];if(parent(v,s)==id)++count[s];}
    }
    return int(count[0]>0)+2*int(count[1]>0);
  };
  for(int v:ids)if(t::topology_detail::life_sex_value(context_life,v)!=1&&t::topology_detail::life_sex_value(context_life,v)!=2){
    int r=role(v);if(r==1||r==2)context_life.sex_overrides.emplace_back(v,r);
  }
  out.focal=index(focal);out.candidate=candidate?index(candidate):-1;out.proposed_parent_side=side;return out;
 }
 t::Input input(int candidate,int side,int opposite,t::Life& context_life){
  auto parent=[&](int v,int s){if(v==focal)return side==0?0:s==side-1?0:opposite;
   int p=parents(v-1,s);return v==candidate&&p==focal?0:p;};
  if(local_inputs)return local_input(candidate,side,opposite,context_life,parent);
  std::set<int> ids=ledger.observed;ids.insert(focal);if(candidate)ids.insert(candidate);if(opposite)ids.insert(opposite);
  for(int v:ledger.observed)for(int s=0;s<2;++s)if(parent(v,s))ids.insert(parent(v,s));
  t::Input out;out.max_typed_nodes=4096;out.max_nodes=8192;out.optional_source_birth_budget=0;out.life=&context_life;
  std::vector<int> index(parents.nrow()+1,-1);for(int v:ids){index[v]=out.nodes.size();t::TopologyNode node;
   node.typed_id=v;node.observed=ledger.observed.count(v)||v==focal||v==candidate||v==opposite;node.label="T"+std::to_string(v);out.nodes.push_back(node);out.retained_source_birth_ids.push_back(v);}
  for(int v:ids){auto& node=out.nodes[index[v]];if(ledger.observed.count(v))for(int s=0;s<2;++s){int p=parent(v,s);if(p&&index[p]>=0)node.parents[s]=index[p];}
   node.has_birth=node.parents[0]>=0||node.parents[1]>=0;}
  for(const auto& edge:ledger.split_parents){int child=index.at(edge.first),side=edge.second-1;
   // Copies are fresh per source incidence and shared by every hypothesis.
   t::TopologyNode copy;copy.label="LedgerCopy_"+std::to_string(edge.first)+"_"+std::to_string(edge.second);
   out.nodes[child].parents[side]=out.nodes.size();out.nodes.push_back(std::move(copy));
  }
  std::vector<int> roles(parents.nrow()+1);for(int v=1;v<=parents.nrow();++v)for(int s=0;s<2;++s){int p=parent(v,s);if(p)roles[p]|=1<<s;}
  for(int v:ids)if(t::topology_detail::life_sex_value(context_life,v)!=1&&t::topology_detail::life_sex_value(context_life,v)!=2&&(roles[v]==1||roles[v]==2))context_life.sex_overrides.emplace_back(v,roles[v]);
  out.focal=index[focal];out.candidate=candidate?index[candidate]:-1;out.proposed_parent_side=side;return out;
 }
 t::TopologyModel build(const t::Input& input,int code,const t::topology_detail::Builder* source=nullptr){
  ++models_built;
  if(code==25||code==26)return t::full_catalog::detail::build(
   input,reverse_po_definition(),true,unsigned(code-25),code,source).topology;
  if(code>=28){auto it=variants.find(code);if(it==variants.end())throw std::invalid_argument("unknown integer full hypothesis");
   const auto& v=it->second;return t::full_catalog::detail::build(input,full[v.definition],v.reverse,v.bits,v.order,source,parentless_hspo).topology;}
  if(code<0||code>=int(simple.size()))throw std::invalid_argument("unknown integer ordinary hypothesis");
  auto model=simple[code];if(t::topology_detail::Builder::try_bound(input,model))return model;try{if(input_checks::enabled()){if(source&&!source->bound_to(input))throw std::invalid_argument("ordinary motif source belongs to another input");}auto builder=[&]{return source?*source:t::topology_detail::Builder(input);}();builder.apply(model);builder.finish(model);}
  catch(const t::topology_detail::Failure& error){model.status=error.status;model.reason=error.what();model.nodes.clear();model.topological_order.clear();model.original_to_canonical.clear();model.focal=model.candidate=-1;}t::topology_detail::Builder::remember_bound(input,model);return model;
 }

 struct Context {
  t::Life life;t::Input input;id::CatalogIdentities identities;
  size_t ready_models=0;std::unordered_set<id::AtomKey,id::AtomHash> invariant;
  std::map<int,t::TopologyModel> built;std::unique_ptr<t::topology_detail::Builder> source;
  Context(Visit& visit,int candidate,int side,int opposite):life(visit.lifetime_view()),identities(visit.parents.nrow()){
   input=visit.input(candidate,side,opposite,life);if(visit.reuse_source){try{source=std::make_unique<t::topology_detail::Builder>(input);++visit.source_builds;}catch(const t::topology_detail::Failure&){}catch(const std::bad_alloc&){}}auto unrelated=visit.build(input,1,source.get());if(source)++visit.source_copies;
   if(unrelated.ready()){auto ids=identities.nodes(unrelated);for(const auto& a:pm::structural_atoms(unrelated))invariant.insert(id::atom(a,ids));}
   built.emplace(1,std::move(unrelated));++visit.contexts;
  }
 };
 Model relationship(Context& context,int code,int opposite,bool reduce_invariant_core=false,bool consume=false,const t::TopologyModel* finished=nullptr){

  Model out;
  if(parentless_hspo&&code>=28&&code<=31){
   try{t::full_catalog::detail::require_hspo_parentless(context.input);}
   catch(const t::topology_detail::Failure&){
    // A transferred/cached topology cannot bypass current admission.
    out.topology=build(context.input,code,context.source.get());return out;
   }
  }
  if(finished){out.finished=finished;

  }else{auto at=context.built.find(code);if(at==context.built.end())at=context.built.emplace(code,build(context.input,code,context.source.get())).first,source_copies+=bool(context.source);
   if(consume){out.topology=std::move(at->second);context.built.erase(at);}else out.topology=at->second;}
  const auto& m=out.source();
  if(m.ready()){out.shared_context=context.ready_models++>0;out.core.insert(m.focal);if(!reduce_invariant_core)out.core.insert(m.candidate);
   if(!reduce_invariant_core&&opposite)for(size_t v=0;v<m.nodes.size();++v)if(m.nodes[v].typed_id==opposite)out.core.insert(v);
   auto ids=context.identities.nodes(m);for(const auto& a:pm::structural_atoms(m)){if(!context.invariant.count(id::atom(a,ids))){auto vars=pm::scope(a);out.core.insert(vars.begin(),vars.end());}else ++out.invariant_factors;}
   for(const auto& birth:m.structural_births)if(!reduce_invariant_core&&birth.required_by_motif){out.core.insert(birth.child);for(int p:birth.parents)if(p>=0)out.core.insert(p);}
  }return out;
 }

 struct Raw {
  bool seeded_null=false,used_null=false;
  t::Life life;t::Input input;::relationship_topology::configuration::Boundary boundary;
  id::CatalogIdentities identities;size_t ready_models=0;std::unordered_set<id::AtomKey,id::AtomHash> invariant;
  std::map<std::array<int,2>,::relationship_topology::configuration::ModelPair> built;
  Raw(Visit& visit):life(visit.lifetime_view()),identities(visit.parents.nrow()){
   input=visit.input(0,0,0,life);input.optional_source_birth_budget=-1;input.retained_source_birth_ids.clear();
   for(auto& n:input.nodes)if(n.typed_id>0&&!n.observed){n.observed=true;n.neutral_phi=true;}
   std::vector<int>d,s;for(int v:visit.dam_candidates)d.push_back(local(v));for(int v:visit.sire_candidates)s.push_back(local(v));
   boundary=::relationship_topology::configuration::freeze(input,d,s);if(!boundary.ready)throw std::runtime_error(boundary.reason);
   auto null=::relationship_topology::configuration::build(boundary,-1,-1,nullptr,false);if(!null.ready())throw std::runtime_error(null.reason);
   t::topology_detail::restore_source_founders(null.denominator,input,false);
   auto ids=identities.nodes(null.denominator);for(const auto& a:pm::structural_atoms(null.denominator))invariant.insert(id::atom(a,ids));
   // The invariant construction already produced this exact null pair from
   // this immutable boundary. Reuse it if the scorer subsequently requests it.
   if(visit.reuse_source)try{t::topology_detail::restore_source_founders(null.numerator,input,false);built.emplace(std::array<int,2>{{0,0}},std::move(null));seeded_null=true;}catch(const std::exception&){}

  }
  int local(int v)const {if(!v)return -1;for(size_t i=0;i<input.nodes.size();++i)if(input.nodes[i].typed_id==v)return i;throw std::runtime_error("parent missing from ledger");}
 };
 const ::relationship_topology::configuration::ModelPair& raw_pair(Raw& context,int dam,int sire){
  namespace cfg=::relationship_topology::configuration;
  std::array<int,2> key{{dam,sire}};auto at=context.built.find(key);
  if(at==context.built.end()){auto pair=cfg::build(context.boundary,context.local(dam),context.local(sire),nullptr,false);
   if(pair.ready()){t::topology_detail::restore_source_founders(pair.numerator,context.input,dam||sire);t::topology_detail::restore_source_founders(pair.denominator,context.input,dam||sire);}
   else {pair.numerator.status=pair.denominator.status=pair.status==cfg::BuildStatus::HardInvalid?t::TopologyStatus::Inapplicable:t::TopologyStatus::Unsupported;pair.numerator.reason=pair.denominator.reason=pair.reason;}
   at=context.built.emplace(key,std::move(pair)).first;models_built+=2;
  }
  return at->second;
 }
 Model raw(Raw& context,int dam,int sire,int term,bool reduce_invariant_core=false,const t::TopologyModel* finished=nullptr){

  Model out;std::array<int,2> key{{dam,sire}};
  if(finished){out.finished=finished;

  }else{const auto& pair=raw_pair(context,dam,sire);out.topology=term?pair.denominator:pair.numerator;}

  if(context.seeded_null&&!context.used_null&&key==std::array<int,2>{{0,0}}){context.used_null=true;++raw_null_reuses;}
  out.name_override="raw_"+std::to_string(dam)+"_"+std::to_string(sire)+(term?"_D":"_N");if(!finished)out.topology.name=out.name_override;const auto& m=out.source();
  if(m.ready()){out.shared_context=context.ready_models++>0;out.core.insert(m.focal);if(!reduce_invariant_core)for(size_t v=0;v<m.nodes.size();++v)if((dam&&m.nodes[v].typed_id==dam)||(sire&&m.nodes[v].typed_id==sire))out.core.insert(v);
   auto ids=context.identities.nodes(m);for(const auto& a:pm::structural_atoms(m)){if(!context.invariant.count(id::atom(a,ids))){auto vars=pm::scope(a);out.core.insert(vars.begin(),vars.end());}else ++out.invariant_factors;}}
  return out;
 }
};
}
#endif
