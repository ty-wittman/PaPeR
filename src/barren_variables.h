#ifndef PAP_BARREN_VARIABLES_H
#define PAP_BARREN_VARIABLES_H
namespace pap_internal::pedigree_likelihood {
// Remove normalized generative leaves before any numerical ordering/table
// admission. Retained outputs and opaque evidence factors are never removed.
// Normalizers of supplied floating-point priors/population gametes are kept
// explicitly, including exact zeros; no assumption of bitwise unit sums.
inline void prune_barren(Plan& p,std::vector<long long>& names,int& count,
                         size_t& model_unaries,size_t& model_births,int& retained) {
 if(!p.ancestral||!p.limits.barren_pruning||!reversible_keys()||!p.pieces.empty())return;
 for(const auto& b:p.births)if(b.excluded_brood)return;
 p.genotype_variables=count;require_acyclic_plan(p);
 std::vector<int> children(count),owners(count),birth_at(count,-1),prior_at(count,-1);
 std::vector<unsigned char> evidence(count),dead(count),removed_birth(p.births.size());
 for(size_t i=0;i<p.unaries.size();++i){auto u=p.unaries[i];if(u.source==Plan::Prior){++owners[u.var];prior_at[u.var]=int(i);}else evidence[u.var]=1;}
 for(size_t i=0;i<p.births.size();++i){const auto& b=p.births[i];++owners[b.vars[0]];birth_at[b.vars[0]]=int(i);for(int s=1;s<3;++s)if(b.vars[s]>=0)++children[b.vars[s]];}
 for(const auto& f:p.pairs)for(int v:f.vars)if(v>=0)evidence[v]=1;
 if(retained>=0)evidence[retained]=1;
 std::vector<int> todo;
 for(int v=0;v<count;++v)if(!evidence[v]&&!children[v]&&owners[v]==1)todo.push_back(v);
 for(size_t j=0;j<todo.size();++j){int v=todo[j];if(dead[v]||evidence[v]||children[v]||owners[v]!=1)continue;dead[v]=1;++p.barren_removed;
  if(birth_at[v]>=0){int i=birth_at[v];removed_birth[i]=1;const auto b=p.births[i];
   for(int s=1;s<3;++s)if(b.vars[s]<0)++p.barren_population_normalizers;
   for(int s=1;s<3;++s)if(int a=b.vars[s];a>=0){--children[a];if(!children[a]&&!evidence[a]&&owners[a]==1)todo.push_back(a);}
  }else if(prior_at[v]>=0)++p.barren_prior_normalizers;
 }
 if(!p.barren_removed)return;
 std::vector<int> renumber(count,-1);std::vector<long long> kept_names;
 for(int v=0;v<count;++v)if(!dead[v]){renumber[v]=int(kept_names.size());kept_names.push_back(names[v]);}
 std::vector<Plan::Unary> unaries;size_t kept_unaries=0;
 for(size_t i=0;i<p.unaries.size();++i)if(!dead[p.unaries[i].var]){auto u=p.unaries[i];u.var=renumber[u.var];unaries.push_back(u);kept_unaries+=i<model_unaries;}
 std::vector<Plan::Birth> births;size_t kept_births=0;
 for(size_t i=0;i<p.births.size();++i)if(!removed_birth[i]){auto b=p.births[i];for(int& v:b.vars)if(v>=0)v=renumber[v];births.push_back(b);kept_births+=i<model_births;}
 for(auto& f:p.pairs)for(int& v:f.vars)if(v>=0)v=renumber[v];
 if(retained>=0)retained=renumber[retained];
 names=std::move(kept_names);count=int(names.size());p.genotype_variables=count;
 p.unaries=std::move(unaries);p.births=std::move(births);model_unaries=kept_unaries;model_births=kept_births;
 p.family_terms=int(p.births.size()+p.pairs.size());
 p.steiner_individuals=int(std::count_if(names.begin(),names.end(),[](long long x){return x>=2000000000LL&&x<3000000000LL;}));
}
}
#endif
