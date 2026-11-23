#include "pap_removed_path.h"
#ifndef PAP_RESOURCE_PLAN_H
#define PAP_RESOURCE_PLAN_H
#include <climits>
#include <type_traits>
namespace pap_internal::pedigree_likelihood::resource_plan {
struct Factor {int kind=0,id=0;std::vector<int> scope;};
struct Step {
  int op=0,variable=-1,auxiliary=-1,parent=-1,select=-1,role=-1,out=-1;
  std::vector<int> inputs,axes,scatter;
  std::size_t cells=0,offset=0;
  double work=0;
  std::vector<std::size_t> axis0,axis1,strides,out_strides;
};
struct Execution {
  int condition=-1,G=0,k=0;
  std::vector<int> domains,original_domains,restricted;
  std::vector<Factor> factors;
  std::vector<Step> steps;
  std::vector<int> scalars;
  std::size_t arena=0,peak_live=0,states_searched=0,workspace_limit=0;
  double work=0,max_step_work=0,predicted_seconds=0;
};
inline bool has(const std::vector<int>& a,int v){return std::find(a.begin(),a.end(),v)!=a.end();}
inline std::size_t cells(const std::vector<int>& axes,const std::vector<int>& dom){
  std::size_t n=1;for(int v:axes){if(n>std::numeric_limits<std::size_t>::max()/std::size_t(dom[v]))return std::numeric_limits<std::size_t>::max();n*=dom[v];}return n;
}
struct State {std::vector<Factor> live;std::vector<Step> steps;double cost=0,work=0;std::size_t peak=0;};
// One deterministic preparation budget for the WHOLE key, including every
// conditioned retry. Charge rejected proposals too, before building a state.
struct SearchBudget {
 static constexpr std::size_t default_limit=8192;
 std::size_t limit=default_limit,attempts=0;bool capped=false;
 explicit SearchBudget(std::size_t n=default_limit):limit(n){}
};
struct SearchAudit {std::size_t keys=0,capped_keys=0,capped_refused=0,attempts=0,admitted=0;};
inline SearchAudit& search_audit(){thread_local SearchAudit x;return x;}

inline std::vector<int> identity(const State& s){
  std::vector<std::vector<int>> fs;for(const auto& f:s.live){std::vector<int> row{f.kind};row.insert(row.end(),f.scope.begin(),f.scope.end());fs.push_back(std::move(row));}
  std::sort(fs.begin(),fs.end());std::vector<int> key;for(const auto& f:fs){key.push_back(int(f.size()));key.insert(key.end(),f.begin(),f.end());}return key;
}
inline bool cheaper(const State& a,const State& b){return std::tie(a.cost,a.peak,a.work)<std::tie(b.cost,b.peak,b.work);}
inline std::vector<State> transitions(const Execution& e,const State& s,const ResourceLimits& limit,SearchBudget& budget){ pap_internal::removed_path("transitions"); }
inline std::shared_ptr<Execution> search(const std::vector<Factor>& factors,const std::vector<int>& domains,int G,int k,const ResourceLimits& limits,int condition=-1,std::size_t max_states=150000,std::size_t width=256,SearchBudget* shared_budget=nullptr){
  // Standalone callers retain their existing accepted-state bound (S37).
  // The broad planner supplies its stricter whole-key proposal budget.
  SearchBudget local_budget(std::numeric_limits<std::size_t>::max());auto& budget=shared_budget?*shared_budget:local_budget;
  // An uncovered variable is already integrated by Program::log_uncovered;
  // conditioning it again would count its domain twice. Only covered genotype
  // variables are eligible for the explicit outer sum.
  if(condition>=0&&(condition>=int(domains.size())||domains[condition]!=G||!std::any_of(factors.begin(),factors.end(),[&](const Factor& f){return has(f.scope,condition);})))return nullptr;
  auto e=std::make_shared<Execution>();e->factors=factors;e->original_domains=e->domains=domains;e->G=G;e->k=k;e->condition=condition;
  if(condition>=0){e->domains[condition]=1;for(const auto& f:factors)if(f.kind==Program::GameteTable&&f.scope[0]==condition){e->domains[f.scope[1]]=2;e->restricted.push_back(f.scope[1]);}}
  std::vector<State> front(1);front[0].live=factors;std::unique_ptr<State> best;
  for(std::size_t depth=0;depth<=domains.size()&&!front.empty();++depth){
    std::map<std::vector<int>,State> next;
    for(const auto& state:front){
      if(budget.capped)break;
      if(state.live.empty()){if(!best||cheaper(state,*best))best=std::make_unique<State>(state);continue;}
      for(auto& candidate:transitions(*e,state,limits,budget)){
        if(++e->states_searched>max_states)return nullptr;
        if(condition>=0&&candidate.work*G>limits.work())continue;
        auto key=identity(candidate);auto it=next.find(key);if(it==next.end())next.emplace(std::move(key),std::move(candidate));else if(cheaper(candidate,it->second))it->second=std::move(candidate);
      }
    }
    if(budget.capped)break;
    front.clear();for(auto& row:next)front.push_back(std::move(row.second));std::stable_sort(front.begin(),front.end(),cheaper);if(front.size()>width)front.resize(width);
  }
  if(!best)return nullptr;e->steps=std::move(best->steps);e->work=best->work*(condition>=0?G:1);e->predicted_seconds=best->cost*(condition>=0?G:1);
  // Compile interval allocation with the same non-overlap certificate used by
  // the established evaluator. Original inputs are external, immutable here.
  Program layout;layout.G=G;
  for(std::size_t i=0;i<factors.size();++i)layout.slots.push_back({Program::PairTable,0,0});
  for(const auto& s:e->steps){Program::Step t;t.inputs=s.inputs;t.cells=s.cells;t.out=s.out;layout.slots.push_back({Program::Produced,0,0});layout.steps.push_back(t);if(s.axes.empty())layout.scalars.push_back(s.out);e->max_step_work=std::max(e->max_step_work,s.work);}
  auto stats=live_arena::pack_program(layout,limits.cells());if(!stats.within_limit)return nullptr;
  e->arena=layout.arena;e->peak_live=stats.concurrent_live_peak_cells;e->scalars=layout.scalars;e->workspace_limit=limits.cells()-e->arena;
  for(const auto& f:factors)if(f.kind==Program::UnaryTable){if(e->workspace_limit<std::size_t(G))return nullptr;e->workspace_limit-=G;}
  for(std::size_t i=0;i<e->steps.size();++i)e->steps[i].offset=layout.steps[i].offset;
  auto all=e->factors;
  for(auto& s:e->steps){
    s.axis0.assign(s.inputs.size(),0);s.axis1=s.axis0;s.strides.assign(s.inputs.size()*s.axes.size(),0);
    std::size_t stride=1;for(int v:s.axes){s.out_strides.push_back(stride);stride*=e->domains[v];}
    for(std::size_t i=0;i<s.inputs.size();++i){const auto& f=all[s.inputs[i]];stride=1;
      for(int v:f.scope){if(v==s.variable)s.axis0[i]=stride;else if(v==s.auxiliary)s.axis1[i]=stride;else s.strides[i*s.axes.size()+std::size_t(std::find(s.axes.begin(),s.axes.end(),v)-s.axes.begin())]=stride;stride*=f.kind==Program::Produced?e->domains[v]:e->original_domains[v];}
    }
    all.push_back({Program::Produced,s.out,s.axes});
  }
  return e;
}
inline std::shared_ptr<Execution> plan(const Program& input,const std::vector<std::vector<int>>& scopes,const std::vector<int>& dom,SearchBudget* supplied=nullptr){
  SearchBudget local;auto& effort=supplied?*supplied:local;
  std::vector<Factor> fs;for(std::size_t i=0;i<scopes.size();++i)fs.push_back({input.slots[i].kind,int(i),scopes[i]});
  auto best=search(fs,dom,input.G,input.k,input.limits,-1,150000,256,&effort);
  if(!best)for(std::size_t g=0;g<dom.size()&&!effort.capped;++g)if(dom[g]==input.G){
    auto candidate=search(fs,dom,input.G,input.k,input.limits,int(g),150000,256,&effort);
    if(candidate&&(!best||candidate->predicted_seconds<best->predicted_seconds))best=std::move(candidate);
  }
  auto& audit=search_audit();++audit.keys;audit.attempts+=effort.attempts;
  audit.capped_keys+=effort.capped;audit.capped_refused+=effort.capped&&!best;audit.admitted+=bool(best);
  return best;
}

struct RunStats {std::uint64_t products=0,multiplies=0,states=0,matrix_calls=0;};
// Exact scalar execution of certified scopes. Fixed variables are substitutions;
// only gamete states proven zero by their original transmission factor vanish.
// Original input tables remain immutable through all conditioning states.
template<class Real>
inline bool execute_precision(const Execution& e,const LocusView& L,const Transmission& T,
                    const std::vector<const double*>& source,Real* arena,double& log_value,RunStats* stats=nullptr, const std::vector<const double*>* source_logs=nullptr){ pap_internal::removed_path("execute_precision"); }
inline bool execute(const Execution& e,const LocusView& L,const Transmission& T,
                    const std::vector<const double*>& source,double* arena,double& log_value,RunStats* stats=nullptr){
  return execute_precision<double>(e,L,T,source,arena,log_value,stats);
}

}
#endif
