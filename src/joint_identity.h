#ifndef PAP_JOINT_IDENTITY_H
#define PAP_JOINT_IDENTITY_H

// Symbolic evidence scope: selected endpoints and their assigned ancestors.
// A visit-wide union remains available only as a numerical reference.
// This is a symbolic adapter to pedigree likelihood, not another likelihood engine.
#include "candidate_catalog.h"
#include "hypothesis_suite.h"

namespace pap_internal::joint_identity {
namespace t = ::relationship_topology;
struct Form { int code=0; const t::full_catalog::Definition* definition=nullptr; unsigned bits=0; };
struct Model { t::TopologyModel topology; int dam_form=0,sire_form=0; };

inline void apply(const Form& form,t::topology_detail::Builder& builder,
                  const t::Input& input) {
  if(!form.definition)return;
  const auto& d=*form.definition;
  if(!d.parent_numerator||!d.registered_class)
    throw std::logic_error("Identity envelope contains an unreviewed/rival form");
  if(const auto* p=t::full_catalog::detail::compiled_motif(d,form.bits)) {
    if(!p->pre_failure.empty())throw t::topology_detail::Failure(t::TopologyStatus::Inapplicable,p->pre_failure);
    t::full_catalog::detail::apply_motif(*p,builder,input,false);return;
  }
  std::map<std::string,int> roles,bound{{"A",input.focal},{"B",input.candidate}};
  for(size_t r=0;r<d.births.size();++r)for(int s=0;s<2;++s)
    roles[d.births[r].parents[s]]|=1<<(s^((form.bits>>r)&1U));
  for(const auto& [name,mask]:roles)if(mask==3)
    throw t::topology_detail::Failure(t::TopologyStatus::Inapplicable,"Joint parent form has incompatible parental roles");
  std::vector<unsigned char> done(d.births.size());size_t remaining=done.size();
  while(remaining){bool progress=false;
    for(size_t r=0;r<d.births.size();++r)if(!done[r]&&bound.count(d.births[r].child)){
      const auto& birth=d.births[r];
      for(int s=0;s<2;++s){int slot=builder.parent_slot(bound.at(birth.child),1+(s^((form.bits>>r)&1U)));
        auto [at,fresh]=bound.emplace(birth.parents[s],slot);
        if(!fresh)at->second=builder.unify_nodes(slot,at->second);
      }done[r]=1;--remaining;progress=true;
    }
    if(!progress)throw t::topology_detail::Failure(t::TopologyStatus::Unsupported,"Unanchored joint parent form");
  }
  const auto names=t::full_catalog::detail::roles(d);
  for(size_t a=0;a<names.size();++a)for(size_t b=a+1;b<names.size();++b)
    builder.require_distinct_nodes(bound.at(names[a]),bound.at(names[b]),"Joint parent form distinct roles coincide");
}

// Immutable revision/focal facts. Endpoint sets and their ancestry remain
// configuration-local; no other candidate's evidence enters a score.
struct VisitFacts {
  std::set<int> excluded_candidates;
  std::vector<int> parental_roles;
  std::vector<Form> forms{{}};
  VisitFacts(candidate_catalog::Visit& visit,hypothesis_suite::Suite suite,
      bool ancestral,bool focal_role_cavity):parental_roles(visit.parents.nrow()+1) {
    // Freeze before looking at any row, score or winning relationship side.
    if(ancestral){
      // Static CSR candidates may have become descendants after earlier moves.
      // They are no longer legal parent configurations and must not turn the
      // focal's offspring into ancestral identity evidence.
      std::vector<std::vector<int>> children(visit.parents.nrow()+1);
      for(int v=0;v<visit.parents.nrow();++v)for(int s=0;s<2;++s)
        if(int p=visit.parents(v,s))children[p].push_back(v+1);
      std::vector<int> todo{visit.focal};excluded_candidates.insert(visit.focal);
      for(size_t j=0;j<todo.size();++j)for(int child:children[todo[j]])
        if(excluded_candidates.insert(child).second)todo.push_back(child);
    }
    for(int v=0;v<visit.parents.nrow();++v)if(!focal_role_cavity||v+1!=visit.focal)for(int s=0;s<2;++s)
      if(int p=visit.parents(v,s))parental_roles[p]|=1<<s;
    for(const auto& [code,v]:visit.variants)if(hypothesis_suite::contains(suite,code)&&!v.reverse){
      const auto& d=visit.full.at(v.definition);if(d.parent_numerator&&d.registered_class)forms.push_back({code,&d,v.bits});
    }
  }
};

struct Context {
  bool parentless_hspo=false;
  t::Life life;t::Input input;std::vector<Form> forms{{}};
  std::vector<int> roots;
  std::shared_ptr<const VisitFacts> facts;
  std::map<std::array<int,2>,std::vector<Model>> rows;
  std::optional<Model> unrelated_;
  const Model& unrelated(){if(!unrelated_)unrelated_=build(0,0,{},{});return *unrelated_;}
  explicit Context(candidate_catalog::Visit& visit,hypothesis_suite::Suite suite,bool ancestral=false,const std::vector<int>* selected_roots=nullptr,bool ancestral_reference=false,bool focal_role_cavity=true,
      std::shared_ptr<const VisitFacts> shared_facts={}):life(visit.lifetime_view()),
      facts(shared_facts?std::move(shared_facts):std::make_shared<const VisitFacts>(visit,suite,ancestral,focal_role_cavity)) {
    parentless_hspo=suite==hypothesis_suite::Suite::S37;
    const auto& excluded_candidates=facts->excluded_candidates;
    std::set<int> ids{visit.focal};
    for(int candidate:selected_roots?*selected_roots:visit.candidates)if(candidate&&!excluded_candidates.count(candidate))ids.insert(candidate);
    std::vector<int> pending(ids.begin(),ids.end());
    std::map<int,int> depth;for(int v:pending)depth[v]=0;
    for(size_t j=0;j<pending.size();++j){int v=pending[j];if(v==visit.focal|| (ancestral&&!ancestral_reference&&depth[v]>=2))continue;
      for(int s=0;s<2;++s)if(int p=visit.parents(v-1,s))if(ids.insert(p).second){pending.push_back(p);depth[p]=depth[v]+1;}
    }
    input.max_typed_nodes=4096;input.max_nodes=16384;input.optional_source_birth_budget=-1;input.life=&life;
    std::map<int,int> index;
    for(int id:ids){index[id]=int(input.nodes.size());t::TopologyNode n;n.typed_id=id;n.observed=true;n.label="T"+std::to_string(id);input.nodes.push_back(n);roots.push_back(id);}
    input.focal=index.at(visit.focal);input.candidate=-1;
    const auto& parental_roles=facts->parental_roles;
    for(int id:ids)if(t::topology_detail::life_sex_value(life,id)!=1&&t::topology_detail::life_sex_value(life,id)!=2&&
      (parental_roles[id]==1||parental_roles[id]==2))life.sex_overrides.emplace_back(id,parental_roles[id]);
    for(int id:ids)if(id!=visit.focal){auto& n=input.nodes[index.at(id)];bool complete=true;
      for(int s=0;s<2;++s)if(int p=visit.parents(id-1,s))if(!index.count(p))complete=false;
      if(complete)for(int s=0;s<2;++s)if(int p=visit.parents(id-1,s)){n.parents[s]=index.at(p);n.has_birth=true;}}
    forms=facts->forms;
  }
  int local(int id)const{if(!id)return -1;for(size_t v=0;v<input.nodes.size();++v)if(input.nodes[v].typed_id==id)return int(v);throw std::out_of_range("Identity parent outside frozen evidence");}
  Model build(int dam,int sire,const Form& df,const Form& sf)const {
    Model out;out.dam_form=df.code;out.sire_form=sf.code;
    try {
      if(facts->excluded_candidates.count(dam)||facts->excluded_candidates.count(sire))
        throw t::topology_detail::Failure(t::TopologyStatus::Inapplicable,"Focal descendant is not an eligible parent");
      auto source=input;source.candidate=local(dam?dam:sire);
      const auto admit=[&](const Form& form,int parent,int other){
        if(!parentless_hspo||!form.definition||form.definition->relationship!="HS_PO")return;
        auto endpoint=source;endpoint.candidate=local(parent);
        t::full_catalog::detail::require_hspo_parentless(endpoint,other>0);
      };
      // Check the candidate's retained ancestry and the proposed opposite
      // parent before applying either form. The tested PO edge itself is new.
      if(dam)admit(df,dam,sire);
      if(sire)admit(sf,sire,dam);
      if(source.candidate<0)for(size_t v=0;v<source.nodes.size();++v)if(int(v)!=source.focal){source.candidate=int(v);break;}
      if(source.candidate<0){auto frozen=t::configuration::freeze(source,{},{});out.topology=t::configuration::build(frozen,-1,-1,nullptr,false).numerator;return out;}
      t::topology_detail::Builder builder(source);
      if(dam)builder.unify_nodes(builder.parent_slot(source.focal,1),local(dam));
      if(sire)builder.unify_nodes(builder.parent_slot(source.focal,2),local(sire));
      auto endpoint=source;
      if(dam){endpoint.candidate=local(dam);endpoint.proposed_parent_side=1;apply(df,builder,endpoint);}
      if(sire){endpoint.candidate=local(sire);endpoint.proposed_parent_side=2;apply(sf,builder,endpoint);}
      out.topology.name="identity_joint_"+std::to_string(df.code)+"_"+std::to_string(sf.code);
      out.topology.relationship="joint_parent_configuration";
      builder.finish(out.topology);
      t::topology_detail::restore_source_founders(out.topology,input,dam||sire);
    }catch(const t::topology_detail::Failure& e){out.topology.status=e.status;out.topology.reason=e.what();}
    return out;
  }
  const std::vector<Model>& models(int dam,int sire) {
    auto [at,fresh]=rows.try_emplace({{dam,sire}});if(!fresh)return at->second;
    auto& out=at->second;
    std::vector<Form> dams{{}},sires{{}};
    // Biological prechecks prune individual forms; combinations still go
    // through one Builder so both parents' constraints apply jointly.
    if(dam)for(size_t j=1;j<forms.size();++j){auto m=build(dam,sire,forms[j],{});if(m.topology.ready())dams.push_back(forms[j]);else if(m.topology.status!=t::TopologyStatus::Inapplicable)out.push_back(std::move(m));}
    if(sire)for(size_t j=1;j<forms.size();++j){auto m=build(dam,sire,{},forms[j]);if(m.topology.ready())sires.push_back(forms[j]);else if(m.topology.status!=t::TopologyStatus::Inapplicable)out.push_back(std::move(m));}
    for(const auto& d:dams)for(const auto& s:sires)out.push_back(build(dam,sire,d,s));
    return out;
  }
};
} // namespace pap_internal::joint_identity
#endif
