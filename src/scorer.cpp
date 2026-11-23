#include "pap_removed_path.h"
#include "input_checks.h"
// [[Rcpp::plugins(cpp20)]]
// Production pedigree likelihood scorer session for live Stage-1 and frozen decision audits.
#include <Rcpp.h>
#include <bit>
#include "message_graph.h" // symbolic atoms/compile only; never Engine
#include "four_gene_moments.h"
#include "candidate_catalog.h"
#include "joint_identity.h"
#include "optional_work.h"
#include "scorer.h"
#include "execution_config.h"
#include "score_types.h"
#include "pedigree_likelihood.h"
#include "hypothesis_suite.h"
#include "locus_batches.h"
#include <chrono>
#include <optional>
using namespace Rcpp;
namespace b=pap_internal::side_likelihood;
namespace mc=pap_internal::candidate_catalog;

namespace {
// Undefined equal-infinite contrasts stay unavailable.
inline double support_contrast(double a,double b) {
  if(std::isnan(a)||std::isnan(b)||(a==b&&!std::isfinite(a)))
    return std::numeric_limits<double>::quiet_NaN();
  return a-b;
}

using Clock=std::chrono::steady_clock;
struct InputCheckScope {
 pap_internal::input_checks::State state;pap_internal::input_checks::Scope scope;
 explicit InputCheckScope(const List& input):scope(state){auto flag=[&](const char* name){return input.containsElementNamed(name)&&as<bool>(input[name]);};state.reference=!flag("pedigree_trust_prepared");state.share=flag("pedigree_share_finished_topology");}

};

double elapsed(Clock::time_point a){return std::chrono::duration<double>(Clock::now()-a).count();}
NumericVector program_registry_receipt(){
 const auto& r=pap_internal::pedigree_likelihood::program_registry();
 return NumericVector::create(_["entries"]=double(r.size()),
  _["reused_entries_at_start"]=double(pap_internal::pedigree_likelihood::program_registry_reused_entries()),
  _["charged_bytes"]=double(r.bytes()),_["peak_charged_bytes"]=double(r.peak_bytes),
  _["bounded"]=r.bounded?1.:0.,_["limit_bytes"]=double(r.limit),
  _["hits"]=double(r.hits),_["misses"]=double(r.misses),_["evictions"]=double(r.evictions),
  _["declines"]=double(r.declines),_["insertions"]=double(r.insertions),_["clearings"]=double(r.clearings));
}
struct Locus {
 int k,G;std::vector<int> observed;std::vector<std::array<int,2>> alleles;std::vector<double> prior,logp,trans;
 std::vector<double> prior_probability,linear_p,linear_t;
 std::vector<std::vector<double>> phi;std::vector<int> lookup;

 bool canonical_prior=true,canonical_transmission=true,four_gene_range=true;
};
// Keep the primitive range refusal without constructing unused message tables.
void check_primitive_range(const std::vector<double>& logs) {
  double high=-INFINITY;
  for(double x:logs) {
    if(std::isnan(x)||x==INFINITY)throw std::runtime_error("invalid primitive");
    high=std::max(high,x);
  }
  if(high==-INFINITY)return;
  for(double x:logs)if(x!=-INFINITY&&std::exp(x-high)==0)
    throw std::runtime_error("linear message underflow; explicit refusal");
}
std::vector<Locus> read_loci(List compact){
 List prior=compact["prior_ln"],logp=compact["log_p"],trans=compact["transmission_log_weight"],obs=compact["observations_ln"],indices=compact["child_index"];
 std::vector<Locus> out;
 for(int l=0;l<prior.size();++l){Locus z;z.prior=as<std::vector<double>>(prior[l]);z.logp=as<std::vector<double>>(logp[l]);z.trans=as<std::vector<double>>(trans[l]);z.G=z.prior.size();z.k=z.logp.size();IntegerMatrix ix=indices[l];NumericMatrix o=obs[l];
 if(compact.containsElementNamed("observed_genotypes")){
  IntegerMatrix called=compact["observed_genotypes"];if(called.nrow()!=o.ncol()||called.ncol()!=prior.size())throw std::invalid_argument("observed genotype metadata dimensions");
  z.observed.assign(o.ncol()+1,-1);for(int id=1;id<=o.ncol();++id){int g=called(id-1,l);if(g!=NA_INTEGER&&g>=0){if(g>=z.G)throw std::out_of_range("observed genotype metadata index");z.observed[id]=g;}}
 }
 z.lookup.resize(z.k*z.k);std::vector<std::array<int,2>> allele(z.G);
 for(int a=0;a<z.k;++a)for(int v=0;v<z.k;++v){int g=ix(std::min(a,v),std::max(a,v))-1;z.lookup[a+z.k*v]=g;if(a<=v)allele[g]={{a,v}};}
 z.alleles=allele;for(int g=0;g<z.G;++g)for(int a=0;a<z.k;++a){double expected=.5*((allele[g][0]==a)+(allele[g][1]==a));double actual=z.trans[g+z.G*a];double log_expected=expected==0?-INFINITY:std::log(expected);if(!(actual==log_expected||(std::isfinite(actual)&&std::isfinite(log_expected)&&std::fabs(actual-log_expected)<=16*std::numeric_limits<double>::epsilon())))z.canonical_transmission=false;}
 z.phi.resize(o.ncol()+1,std::vector<double>(z.G,0));for(int id=1;id<=o.ncol();++id)for(int g=0;g<z.G;++g)z.phi[id][g]=o(g,id-1);for(const auto& row:z.phi)check_primitive_range(row);
 for(int a=0;a<z.k;++a)for(int b=a;b<z.k;++b){int g=z.lookup[a+z.k*b];
  double expected=z.logp[a]+z.logp[b]+(a==b?0:std::log(2.));
  double tolerance=16*std::numeric_limits<double>::epsilon()*(1+std::fabs(expected)+std::fabs(z.prior[g]));
  if(!(expected==z.prior[g]||(std::isfinite(expected)&&std::isfinite(z.prior[g])&&std::fabs(expected-z.prior[g])<=tolerance)))z.canonical_prior=false;
 }
 check_primitive_range(z.prior);for(double value:z.prior)z.prior_probability.push_back(value==-INFINITY?0:std::exp(value));
 for(double value:z.logp)z.linear_p.push_back(value==-INFINITY?0:std::exp(value));for(double value:z.trans)z.linear_t.push_back(value==-INFINITY?0:std::exp(value));z.four_gene_range=pap_internal::four_gene::representable(z.logp,z.linear_p)&&pap_internal::four_gene::representable(z.prior,z.prior_probability);out.push_back(std::move(z));}
 return out;
}
struct ScorerInput {
 const std::shared_ptr<const pap_internal::ExecutionConfig> config;
 const std::vector<Locus> loci;
 ScorerInput(List compact,const List& options)
   :config(std::make_shared<const pap_internal::ExecutionConfig>(options)),
    loci(read_loci(compact)) {}
};

// One setup for frozen and live scoring. Native callers may temporarily
// override thread-local planning: bind this scope at every scorer entry.
struct Store {
 const pap_internal::ExecutionConfig& config;
 explicit Store(const pap_internal::ExecutionConfig& options):config(options) {}
 std::array<double,3> likelihood_brood_checks{};
 std::string likelihood_brood_failure;
 int likelihood_round=1;
 struct SharedLikelihoodState {
  std::unique_ptr<pap_internal::pedigree_likelihood::Forest> forest;
  std::vector<std::unique_ptr<pap_internal::pedigree_likelihood::LocusMessages>> messages;
  std::vector<std::unique_ptr<pap_internal::pedigree_likelihood::Transmission>> tables;
  std::vector<int> dam,sire;int round=-1;size_t builds=0;
  pap_internal::pedigree_likelihood::LaneTableCache lane_tables;
  std::vector<std::vector<int>> gamete_allele;
  std::vector<std::vector<double>> gamete_weight;
  std::vector<signed char> gamete_state;
 } likelihood_data;
 pap_internal::topology_life_reuse::Cache life_certificates;
 bool topology_reuse=false,setup_reuse=false;
 int shared_revision=0;
};


struct LiveSession {
 std::unique_ptr<pap_internal::pedigree_likelihood::LocusBatches> locus_batches;bool parallel_messages=false;
 std::shared_ptr<const pap_internal::candidate_catalog::ImmutableInput> immutable;size_t parse_builds=0,parse_hits=0;
 std::shared_ptr<const pap_internal::candidate_catalog::ImmutableInput> catalog(List input){
  if(store.setup_reuse&&immutable){if(pap_internal::input_checks::enabled()){if(!immutable->matches(input)){++parse_builds;return std::make_shared<const mc::ImmutableInput>(input);}}++parse_hits;return immutable;}
  auto result=pap_internal::candidate_catalog::ImmutableInput::parse(input);++parse_builds;if(store.setup_reuse)immutable=result;return result;
 }

 const ScorerInput input;
 Store store;
 LiveSession(List compact,List inputs):input(compact,inputs),store(*input.config) {auto flag=[&](const char* n){return inputs.containsElementNamed(n)&&as<bool>(inputs[n]);};pap_internal::pedigree_likelihood::configure_resource_audit(inputs.containsElementNamed("likelihood_resource_audit_directory")?as<std::string>(inputs["likelihood_resource_audit_directory"]):"");if(flag("likelihood_resource_audit_extended_precision")){
 if(!flag("configuration_diagnostics")||pap_internal::pedigree_likelihood::resource_audit().directory.empty())
  throw std::invalid_argument("unvalidated extended precision requires configuration diagnostics and a W1 resource audit");
 pap_internal::pedigree_likelihood::resource_audit().extended_precision=true;
 }pap_internal::pedigree_likelihood::resource_audit().capture_unique=inputs.containsElementNamed("likelihood_resource_audit_unique")&&as<bool>(inputs["likelihood_resource_audit_unique"]);pap_internal::pedigree_likelihood::resource_audit().capture_three_genotype=inputs.containsElementNamed("likelihood_resource_audit_three_genotype")&&as<bool>(inputs["likelihood_resource_audit_three_genotype"]);
 store.likelihood_data.lane_tables.bounded=flag("likelihood_bounded_lane_tables");
 // Private accepted-move replay, without rescoring a search prefix.
 // Canonical current-state cuts must agree with a direct final-state build.
 if(inputs.containsElementNamed("likelihood_replay_initial_parents")) {
  if(!store.config.pedigree_likelihood||!inputs.containsElementNamed("likelihood_refresh")||
     as<int>(inputs["likelihood_refresh"])!=3||
     (!flag("configuration_diagnostics")&&!flag("private_frozen_decision_audit")))
   throw std::invalid_argument("M_ped2 forest replay requires refresh 3 and diagnostics or frozen audit");
  IntegerMatrix initial=inputs["likelihood_replay_initial_parents"],target=inputs["parents"],moves=inputs["likelihood_replay_moves"];
  if(initial.ncol()!=2||target.ncol()!=2||initial.nrow()!=target.nrow()||moves.ncol()!=7)
   throw std::invalid_argument("invalid M_ped2 forest replay shape");
  const int n=initial.nrow();std::vector<int>d(n+1),s(n+1);
  for(int v=1;v<=n;++v)for(int side=0;side<2;++side) {
   int p=initial(v-1,side);if(p<0||p>n||p==v)throw std::invalid_argument("invalid M_ped2 forest replay parent");
   (side?s:d)[v]=p;
  }
  auto& shared=store.likelihood_data;shared.forest=std::make_unique<pap_internal::pedigree_likelihood::Forest>();
  shared.forest->build(d,s,{});shared.forest->ancestral=store.config.likelihood_ancestral;
  for(int row=0;row<moves.nrow();++row) {
   int c=moves(row,0),nd=moves(row,3),ns=moves(row,4);
   if(c<1||c>n||nd<0||nd>n||ns<0||ns>n||nd==c||ns==c||
      d[c]!=moves(row,1)||s[c]!=moves(row,2)||moves(row,6)!=moves(row,5)+1||
      (row&&moves(row,5)!=moves(row-1,6)))
    throw std::invalid_argument("invalid M_ped2 forest replay move chain");
   d[c]=nd;s[c]=ns;shared.forest->update(d,s);
  }
  if(moves.nrow()&&inputs.containsElementNamed("start_graph_revision")&&
     moves(moves.nrow()-1,6)!=as<int>(inputs["start_graph_revision"]))
   throw std::invalid_argument("M_ped2 forest replay revision does not reach the saved graph");
  for(int v=1;v<=n;++v)if(d[v]!=target(v-1,0)||s[v]!=target(v-1,1))
   throw std::invalid_argument("M_ped2 forest replay does not reach the saved graph");
  shared.dam=std::move(d);shared.sire=std::move(s);shared.messages.resize(input.loci.size());
 }
pap_internal::pedigree_likelihood::program_registry_enabled()=!flag("likelihood_no_program_registry");{const bool bounded=!inputs.containsElementNamed("likelihood_bounded_program_registry")||flag("likelihood_bounded_program_registry");double mib=inputs.containsElementNamed("likelihood_program_registry_mib")?as<double>(inputs["likelihood_program_registry_mib"]):16;if(!std::isfinite(mib)||mib!=std::floor(mib)||mib<0||mib>64)throw std::invalid_argument("mped2_program_registry_mib must be 0..64");const bool retain=flag("likelihood_registry_across_passes")&&!flag("likelihood_no_program_registry");
 const double generation=inputs.containsElementNamed("likelihood_registry_generation")?as<double>(inputs["likelihood_registry_generation"]):0;
 pap_internal::pedigree_likelihood::configure_program_registry(bounded,size_t(mib)*1024*1024,retain,generation);}pap_internal::pedigree_likelihood::program_cache_census()={};pap_internal::pedigree_likelihood::resource_plan::search_audit()={};pap_internal::pedigree_likelihood::program_cache_census_enabled()=flag("likelihood_program_census");pap_internal::pedigree_likelihood::program_cache_trace()=inputs.containsElementNamed("likelihood_program_census_path")?as<std::string>(inputs["likelihood_program_census_path"]):"";{bool reversible=store.config.pedigree_likelihood&&!flag("likelihood_no_reversible_keys");for(const auto& z:input.loci){if(!reversible)break;pap_internal::pedigree_likelihood::LocusView v;v.G=z.G;v.k=z.k;v.prior=&z.prior_probability;v.p=&z.linear_p;v.t=&z.linear_t;v.lookup=&z.lookup;v.phi=&z.phi;reversible=pap_internal::pedigree_likelihood::reversible_locus(v);}pap_internal::pedigree_likelihood::reversible_keys()=reversible;}store.shared_revision=inputs.containsElementNamed("start_graph_revision")?as<int>(inputs["start_graph_revision"]):0;
   store.setup_reuse=flag("pedigree_scorer_setup_reuse");
 store.topology_reuse=flag("pedigree_topology_reuse");
 if(inputs.containsElementNamed("pedigree_reuse_byte_limit")){
  double limit=as<double>(inputs["pedigree_reuse_byte_limit"]);
  if(!std::isfinite(limit)||limit<0||limit>64U*1024U*1024U||limit!=std::floor(limit))
   throw std::invalid_argument("M_ped reuse budget must be an integer from zero through 64 MiB");
  store.life_certificates.limit=size_t(limit);
 }
 if(flag("likelihood_batch_loci")){
  const int workers=inputs.containsElementNamed("likelihood_locus_workers")?as<int>(inputs["likelihood_locus_workers"]):1;
  if(workers<1||workers>64||!store.config.pedigree_likelihood||store.config.likelihood_plan_profile)
   throw std::invalid_argument("mped2_batch_loci requires M_ped2 program kernels, no plan profiler, and workers 1..64");
  parallel_messages=flag("likelihood_parallel_messages");
  double threshold=inputs.containsElementNamed("likelihood_parallel_min_operations")?as<double>(inputs["likelihood_parallel_min_operations"]):100000.;
  if(!std::isfinite(threshold)||threshold<0)throw std::invalid_argument("mped2_parallel_min_operations must be finite and nonnegative");
  locus_batches=std::make_unique<pap_internal::pedigree_likelihood::LocusBatches>(workers,int(input.loci.size()),threshold,flag("pedigree_phase_profile"),pap_internal::pedigree_likelihood::LocusBatches::parse_scheduler(inputs.containsElementNamed("likelihood_scheduler")?as<std::string>(inputs["likelihood_scheduler"]):"static"));
  locus_batches->adaptive_lanes=flag("likelihood_adaptive_lanes");
  locus_batches->resource_diagnostics=flag("likelihood_resource_diagnostics");
  if(!pap_internal::pedigree_likelihood::resource_audit().directory.empty()&&workers!=1)throw std::invalid_argument("M_ped2 resource capture requires one numerical worker");
 }
}
 int synchronize_forest(std::vector<int> d,std::vector<int> s) {
   auto& shared=store.likelihood_data;const bool changed=!shared.forest||shared.dam!=d||shared.sire!=s;
   if(shared.forest&&changed&&(store.config.likelihood_refresh==3||(store.config.likelihood_ancestral&&store.config.likelihood_ancestral_message_reuse&&store.config.likelihood_refresh==1))){
    shared.forest->update(d,s);shared.dam=std::move(d);shared.sire=std::move(s);return 1;
   }else if(!shared.forest||(changed&&(store.config.likelihood_refresh==1||shared.round!=store.likelihood_round))){
    shared.forest=std::make_unique<pap_internal::pedigree_likelihood::Forest>();shared.forest->build(d,s,std::vector<int>{});shared.forest->ancestral=store.config.likelihood_ancestral;
    shared.messages.clear();shared.messages.resize(input.loci.size());shared.dam=std::move(d);shared.sire=std::move(s);shared.round=store.likelihood_round;++shared.builds;return 2;
   }
  return 0;
 }

};
struct ScorerVisit:pap_internal::native_demand::LiveScorer {
 bool locus_batches_ready=false;
 List input;LiveSession& session;const std::vector<Locus>& data;std::optional<pap_internal::candidate_catalog::Visit> visit;Store& store;
 std::map<std::array<int,3>,std::unique_ptr<pap_internal::candidate_catalog::Visit::Context>> contexts;
 std::unique_ptr<pap_internal::candidate_catalog::Visit::Raw> raw;
 // pedigree likelihood visit state. Messages live on a forest that is either built for this
 // visit (refresh 0: the focal's links removed, cuts far from the question
 // individuals) or stored in the session and rebuilt per graph revision (1)
 // or per round (2); a stored forest is read through the visit's context.
 struct LikelihoodEntry {bool focal_table=false;std::vector<double> denominators;bool ready=false;b::PartitionState state=b::PartitionState::unavailable;std::string reason;pap_internal::pedigree_likelihood::Plan plan;std::vector<double> values;std::vector<unsigned char> done;};
 std::unique_ptr<pap_internal::pedigree_likelihood::Forest> likelihood_forest;
 std::vector<std::unique_ptr<pap_internal::pedigree_likelihood::LocusMessages>> likelihood_messages;
 const pap_internal::pedigree_likelihood::Forest* likelihood_active=nullptr;std::vector<std::unique_ptr<pap_internal::pedigree_likelihood::LocusMessages>>* likelihood_active_messages=nullptr;pap_internal::pedigree_likelihood::Context likelihood_ctx;
 pap_internal::optional_demand::Ledger likelihood_optional;
 std::map<std::array<int,5>,std::shared_ptr<LikelihoodEntry>> likelihood_entries;std::map<std::vector<long long>,std::shared_ptr<LikelihoodEntry>> likelihood_plan_entries;size_t likelihood_plan_shared=0,likelihood_cyclic_plans=0,likelihood_rejected_cyclic_plans=0;pap_internal::pedigree_likelihood::PairCache likelihood_pairs;pap_internal::pedigree_likelihood::PieceCache likelihood_pieces;std::vector<double> likelihood_focal_table;pap_internal::pedigree_likelihood::Scaled likelihood_exact_focal_table;std::vector<double> likelihood_arena;size_t likelihood_updates=0,likelihood_lane_batches=0,likelihood_lane_fallbacks=0;std::vector<double> likelihood_lane_arena;std::vector<std::vector<int>> likelihood_groups;double likelihood_class_seconds[3]={0,0,0};size_t likelihood_class_plans[3]={0,0,0};static int likelihood_class(const std::array<int,5>& key){return key[0]==3?(key[1]>0&&key[2]>0?1:0):2;}
 double likelihood_operations=0,likelihood_build_seconds=0,likelihood_topology_seconds=0,likelihood_eval_seconds=0,likelihood_forest_seconds=0;struct PlanRecord{int kind,candidate,side,opposite,hypothesis,variables,genotype_variables,width,steiner,births,pairs;double operations=0,seconds=0,value0=NA_REAL;};std::vector<PlanRecord> likelihood_profile;std::map<std::array<int,5>,size_t> likelihood_profile_index;std::array<double,5> likelihood_scope{};size_t likelihood_questions=0,likelihood_questions_cut=0;std::vector<unsigned char> likelihood_ancestry_cuts;const std::vector<unsigned char>* likelihood_ancestry_cut_view=nullptr;size_t likelihood_rebuilds=0,likelihood_plans=0,likelihood_evaluations=0,likelihood_variables=0,likelihood_family_terms=0,likelihood_steiner=0,likelihood_width_max=0,likelihood_variables_max=0,likelihood_redundant=0;double likelihood_seconds=0,likelihood_message_seconds=0;
 pap_internal::pedigree_likelihood::LocusView likelihood_view(size_t l)const{const auto& z=data.at(l);pap_internal::pedigree_likelihood::LocusView v;v.G=z.G;v.k=z.k;v.prior=&z.prior_probability;v.p=&z.linear_p;v.t=&z.linear_t;v.lookup=&z.lookup;v.phi=&z.phi;v.observed=&z.observed;
  // Sparse gametes once per locus when every genotype transmits at most two alleles.
  auto& m=store.likelihood_data;if(m.gamete_state.size()!=data.size()){m.gamete_state.assign(data.size(),0);m.gamete_allele.assign(data.size(),{});m.gamete_weight.assign(data.size(),{});}
  if(!store.config.likelihood_dense_gametes){if(!m.gamete_state[l])m.gamete_state[l]=pap_internal::pedigree_likelihood::sparse_gametes(v,m.gamete_allele[l],m.gamete_weight[l])?1:-1;if(m.gamete_state[l]>0){v.gamete_allele=&m.gamete_allele[l];v.gamete_weight=&m.gamete_weight[l];}}
  return v;}
 void likelihood_prepare(){
  if(likelihood_active)return;

  auto t=Clock::now();int n=visit->parents.nrow();std::vector<int> d(n+1,0),s(n+1,0);
  for(int v=1;v<=n;++v){d[v]=visit->parents(v-1,0);s[v]=visit->parents(v-1,1);}
  if(store.config.likelihood_refresh==0){
   d[visit->focal]=s[visit->focal]=0;std::vector<int> seeds=visit->candidates;seeds.push_back(visit->focal);
   likelihood_forest=std::make_unique<pap_internal::pedigree_likelihood::Forest>();likelihood_forest->build(d,s,seeds);likelihood_forest->ancestral=store.config.likelihood_ancestral;likelihood_messages.clear();likelihood_messages.resize(data.size());
   likelihood_active=likelihood_forest.get();likelihood_active_messages=&likelihood_messages;++likelihood_rebuilds;
  }else{
   const int action=session.synchronize_forest(std::move(d),std::move(s));
   likelihood_updates+=action==1;likelihood_rebuilds+=action==2;
   likelihood_active=store.likelihood_data.forest.get();
   likelihood_active_messages=&store.likelihood_data.messages;
  }
  if(store.config.likelihood_light_bookkeeping)likelihood_ancestry_cut_view=&likelihood_active->ancestry_cut_flags_cached();
  else{likelihood_ancestry_cuts=likelihood_active->ancestry_cut_flags();likelihood_ancestry_cut_view=&likelihood_ancestry_cuts;}
  likelihood_ctx.reset(*likelihood_active,visit->focal,int(data.size()),store.config.likelihood_latent_brood_exclusion,store.config.likelihood_ancestral_message_reuse);likelihood_forest_seconds+=elapsed(t);
 }
 const pap_internal::pedigree_likelihood::LocusMessages& likelihood_locus(size_t l){
  likelihood_prepare();auto& slot=likelihood_active_messages->at(l);
  if(!slot){auto t=Clock::now();slot=std::make_unique<pap_internal::pedigree_likelihood::LocusMessages>();slot->latent_brood_exclusion=store.config.likelihood_latent_brood_exclusion;slot->ancestral=store.config.likelihood_ancestral;slot->ancestral_reuse=store.config.likelihood_ancestral_message_reuse;pap_internal::pedigree_likelihood::compute_messages(*likelihood_active,likelihood_view(l),*slot);likelihood_message_seconds+=elapsed(t);}
  return *slot;
 }
 const pap_internal::pedigree_likelihood::Transmission& likelihood_table(size_t l){
  auto& tables=store.likelihood_data.tables;if(tables.size()!=data.size())tables.resize(data.size());auto& table=tables[l];
  if(!table){table=std::make_unique<pap_internal::pedigree_likelihood::Transmission>();pap_internal::pedigree_likelihood::build_transmission(likelihood_view(l),*table,store.config.likelihood_implicit_transmission);}
  return *table;
 }
 std::array<int,5> likelihood_entry_key(int hypothesis,pap_internal::native_demand::Tags tags)const{
  int opposite=tags.side==1?tags.sire:tags.dam;
  if(tags.kind==5)return {{5,tags.dam,tags.sire,hypothesis,tags.term}};
  return tags.kind==3?std::array<int,5>{{3,tags.dam,tags.sire,!store.config.likelihood_raw_separate?-1:tags.term,0}}:std::array<int,5>{{1,tags.candidate,tags.side,opposite,hypothesis}};
 }
 bool certify_brood_plan(const ::relationship_topology::TopologyModel& topology,
                         const pap_internal::pedigree_likelihood::Plan& plan,bool raw_table,bool raw) {
  if(!store.config.likelihood_latent_brood_exclusion)return true;
  namespace m=pap_internal::pedigree_likelihood;
  // Structural reference only: never evaluate or retain its numerical tables.
  m::PairCache pairs;m::PieceCache pieces;
  const auto reference=m::plan_model(topology,*likelihood_active,likelihood_ctx,pairs,pieces,raw,
    !store.config.likelihood_no_pieces,raw_table,store.config.likelihood_max_touch,store.config.likelihood_factorized_one,
    store.config.likelihood_high_allele_programs,store.config.likelihood_limits,store.config.likelihood_local_family_pairs,false);
  std::set<std::pair<int,int>> domains;
  for(const auto& z:data)domains.emplace(z.G,z.k);
  for(auto [G,k]:domains) {
    auto inspect=[&](const m::Plan& p,const m::PieceCache& pc) {
      const auto& program=m::program_for(p,G,k);
      int width=program.width;bool ready=program.ready;
      for(const auto& ref:p.pieces){const auto& q=m::program_for(pc.slots.at(ref.slot)->plan,G,k);
        width=std::max(width,q.width);ready=ready&&q.ready;}
      return std::pair<int,bool>{width,ready};
    };
    const auto before=inspect(reference,pieces),after=inspect(plan,likelihood_pieces);
    ++store.likelihood_brood_checks[0];
    if(after.first>before.first)++store.likelihood_brood_checks[1];
    if(before.second&&!after.second)++store.likelihood_brood_checks[2];
    if(after.first>before.first || (before.second&&!after.second)) {
      store.likelihood_brood_failure="M_ped2 brood exclusion certificate: focal "+std::to_string(visit->focal)+
        " k "+std::to_string(k)+" width "+std::to_string(before.first)+" -> "+std::to_string(after.first)+
        " ready "+std::to_string(before.second)+" -> "+std::to_string(after.second);
      return false;
    }
  }
  return true;
 }
 std::map<std::pair<int,int>,std::size_t> likelihood_domain_counts;
 double program_operations(pap_internal::pedigree_likelihood::Plan& plan){
   double work=0;
   if(!store.config.likelihood_light_bookkeeping){for(const auto& z:data)work+=pap_internal::pedigree_likelihood::program_for(plan,z.G,z.k).operations;return work;}
   if(likelihood_domain_counts.empty())for(const auto& z:data)++likelihood_domain_counts[{z.G,z.k}];
   for(const auto& [domain,count]:likelihood_domain_counts)
     work+=count*pap_internal::pedigree_likelihood::program_for(plan,domain.first,domain.second).operations;
   return work;
 }
 void initialize_locus_batches(){
   if(locus_batches_ready)return;
   auto& batches=*session.locus_batches;auto started=Clock::now();
   std::vector<pap_internal::pedigree_likelihood::LocusView> views;std::vector<const pap_internal::pedigree_likelihood::Transmission*> tables;
   for(size_t j=0;j<data.size();++j){views.push_back(likelihood_view(j));tables.push_back(&likelihood_table(j));}
   batches.initialize(std::move(views),std::move(tables),!store.config.likelihood_no_lanes,!store.config.likelihood_narrow_lanes);
   batches.begin_visit(*likelihood_active,visit->focal,*likelihood_active_messages,session.parallel_messages,store.config.likelihood_latent_brood_exclusion,store.config.likelihood_ancestral,store.config.likelihood_ancestral_message_reuse);
   locus_batches_ready=true;likelihood_message_seconds+=elapsed(started);
 }
 void prepare_scores(const std::vector<std::pair<int,pap_internal::native_demand::Tags>>& requests)override{
   if(!prepared_score_batches||!session.locus_batches)return;
   std::vector<pap_internal::pedigree_likelihood::LocusBatches::Request> pending;std::set<LikelihoodEntry*> seen;
   for(const auto& [hypothesis,tags]:requests){
    likelihood_score(hypothesis,tags,0,true);auto& entry=*likelihood_entries.at(likelihood_entry_key(hypothesis,tags));
    if(!entry.ready||entry.done[0]||!seen.insert(&entry).second)continue;
    pending.push_back({&entry.plan,entry.focal_table,&entry.values,&entry.denominators});
   }
   if(pending.empty())return;
   initialize_locus_batches();auto started=Clock::now();auto& batches=*session.locus_batches;
   batches.evaluate_many(std::move(pending),likelihood_pairs,likelihood_pieces);
   for(auto* entry:seen){std::fill(entry->done.begin(),entry->done.end(),1);likelihood_evaluations+=data.size();likelihood_operations+=program_operations(entry->plan);}
   likelihood_pairs.requests=batches.pair_requests();likelihood_pairs.computed=batches.pair_computed();likelihood_ctx.repaired=batches.repairs();likelihood_eval_seconds+=elapsed(started);
   // A wave spans classes; do not mislabel its elapsed time as a per-class time.
   for(double& seconds:likelihood_class_seconds)seconds=std::numeric_limits<double>::quiet_NaN();
 }
 b::PartitionValue likelihood_score(int hypothesis,pap_internal::native_demand::Tags tags,int l,bool prepare_only=false){
  const bool timed=store.config.likelihood_time_evaluations;auto started=timed?Clock::now():Clock::time_point{};int opposite=tags.side==1?tags.sire:tags.dam;
  // Raw numerator and denominator share one plan with the focal retained.
  const bool raw_table=tags.kind==3&&!store.config.likelihood_raw_separate;
  std::array<int,5> key=likelihood_entry_key(hypothesis,tags);
  if(tags.kind==5){++identity_work[2];++identity_work[support_scoring?4:5];}
  auto& entry=likelihood_entries[key];
  likelihood_optional.request(key,bool(entry));
  if(!entry){auto build_started=Clock::now();entry=std::make_shared<LikelihoodEntry>();pap_internal::candidate_catalog::Model built;
   if(tags.kind==5){built.topology=identity_model(tags.dam,tags.sire,hypothesis,tags.term).topology;}
   else if(tags.kind==3){if(!raw)raw=std::make_unique<pap_internal::candidate_catalog::Visit::Raw>(*visit);built=visit->raw(*raw,tags.dam,tags.sire,raw_table?1:tags.term);}
   else{auto& context=contexts[{{tags.candidate,tags.side,opposite}}];if(!context)context=std::make_unique<pap_internal::candidate_catalog::Visit::Context>(*visit,tags.candidate,tags.side,opposite);built=visit->relationship(*context,hypothesis,opposite);}
   const auto& topology=built.source();likelihood_topology_seconds+=elapsed(build_started);
   if(!topology.ready()){entry->state=topology.status==::relationship_topology::TopologyStatus::Inapplicable?b::PartitionState::inapplicable:b::PartitionState::unavailable;entry->reason=topology.reason;}
   else{likelihood_prepare();++likelihood_questions;bool cut=false;for(const auto& node:topology.nodes)if(node.typed_id>0)cut|=likelihood_ancestry_cut_view->at(node.typed_id);likelihood_questions_cut+=cut;try{entry->plan=[&]{return pap_internal::pedigree_likelihood::plan_model(topology,*likelihood_active,likelihood_ctx,likelihood_pairs,likelihood_pieces,tags.kind==3||tags.kind==5,!store.config.likelihood_no_pieces,raw_table,store.config.likelihood_max_touch,store.config.likelihood_factorized_one,store.config.likelihood_high_allele_programs,store.config.likelihood_limits,store.config.likelihood_local_family_pairs,store.config.likelihood_latent_brood_exclusion,store.config.likelihood_ancestral,
      tags.kind==5?identity_context(tags.dam,tags.sire).roots:tags.kind==3?std::vector<int>{tags.dam,tags.sire}:std::vector<int>{tags.candidate,opposite},store.config.likelihood_ancestral_reference);}();
      likelihood_scope[0]+=entry->plan.off_scope_observations_before;likelihood_scope[1]+=entry->plan.descendant_observations_before;likelihood_scope[2]+=entry->plan.barren_removed;likelihood_scope[3]+=entry->plan.off_scope_messages;likelihood_scope[4]+=entry->plan.off_scope_observations_before>0;
      entry->focal_table=raw_table;entry->ready=true;
    // Hypotheses whose plans coincide are one calculation.
    if(!store.config.likelihood_no_dedup){auto plan_key=[&]{return pap_internal::pedigree_likelihood::plan_key(entry->plan,&likelihood_cyclic_plans);}();plan_key.push_back(raw_table);auto at=likelihood_plan_entries.emplace(std::move(plan_key),entry);if(!at.second){entry=at.first->second;++likelihood_plan_shared;likelihood_build_seconds+=elapsed(build_started);goto likelihood_built;}}
    }
    catch(const ::relationship_topology::topology_detail::Failure& failure){
      if(failure.status!=::relationship_topology::TopologyStatus::Inapplicable)throw;if(std::string(failure.what()).find("directed cycle")!=std::string::npos)++likelihood_rejected_cyclic_plans;entry->state=b::PartitionState::inapplicable;entry->reason=failure.what();
      entry->ready=false;entry->values.assign(data.size(),0.0);entry->denominators.assign(data.size(),0.0);entry->done.assign(data.size(),0);goto likelihood_built;
    }
    if(store.config.likelihood_latent_brood_exclusion) {
      bool certified=false;
      try{certified=certify_brood_plan(topology,entry->plan,raw_table,tags.kind==3);}
      catch(...){entry->ready=false;throw;}
      if(!certified){entry->ready=false;entry->reason=store.likelihood_brood_failure;}
    }
    ++likelihood_plans;++likelihood_class_plans[likelihood_class(key)];auto& P=entry->plan;P.width=[&]{const bool actual=!input.containsElementNamed("likelihood_actual_width_domain")||as<bool>(input["likelihood_actual_width_domain"]);return pap_internal::pedigree_likelihood::program_for(P,actual?data.front().G:10,actual?data.front().k:4).width;}();likelihood_variables+=P.variables;likelihood_family_terms+=P.family_terms;likelihood_steiner+=P.steiner_individuals;likelihood_width_max=std::max<size_t>(likelihood_width_max,P.width);likelihood_variables_max=std::max<size_t>(likelihood_variables_max,P.variables);likelihood_redundant+=P.redundant_births;if(store.config.likelihood_plan_profile){likelihood_profile_index[key]=likelihood_profile.size();likelihood_profile.push_back({key[0],key[1],key[2],key[3],key[4],P.variables,P.genotype_variables,P.width,P.steiner_individuals,int(P.births.size()),int(P.pairs.size())});}}
   entry->values.assign(data.size(),0.0);entry->denominators.assign(data.size(),0.0);entry->done.assign(data.size(),0);likelihood_build_seconds+=elapsed(build_started);}
  likelihood_built:
  b::PartitionValue out;
  if(!entry->ready){out.state=entry->state;out.reason=entry->reason;likelihood_seconds+=(timed?elapsed(started):0.);return out;}
  if(prepare_only)return out;
  if(session.locus_batches&&!entry->done[l]){
   auto& batches=*session.locus_batches;
   initialize_locus_batches();
   auto started=timed?Clock::now():Clock::time_point{};batches.evaluate_all(entry->plan,likelihood_pairs,likelihood_pieces,entry->focal_table,entry->values,entry->denominators);
   std::fill(entry->done.begin(),entry->done.end(),1);likelihood_evaluations+=data.size();
   likelihood_pairs.requests=batches.pair_requests();likelihood_pairs.computed=batches.pair_computed();likelihood_ctx.repaired=batches.repairs();
   likelihood_operations+=program_operations(entry->plan);
   double spent=(timed?elapsed(started):0.);likelihood_eval_seconds+=spent;likelihood_class_seconds[likelihood_class(key)]+=spent;
  }
  // Pending loci with this locus's allele count run together in vector lanes.
  if(!entry->done[l]&&!store.config.likelihood_no_lanes&&!store.config.likelihood_plan_profile){
   namespace likelihood_data=pap_internal::pedigree_likelihood;constexpr int MaxW=likelihood_data::kWideLanes;
   if(likelihood_groups.size()!=data.size()){likelihood_groups.assign(data.size(),{});std::map<std::pair<int,int>,std::vector<int>> by;for(size_t x=0;x<data.size();++x)by[{data[x].G,data[x].k}].push_back(int(x));for(auto& [shape,members]:by)for(int x:members)likelihood_groups[x]=members;}
   // Eight AVX-512 lanes when that many loci wait and the CPU has them; else four.
   std::array<int,MaxW> lanes;int count=0;lanes[count++]=l;
   const int cap=!store.config.likelihood_narrow_lanes&&likelihood_data::wide_lanes_available()?MaxW:likelihood_data::kLanes;
   for(int x:likelihood_groups[l]){if(count==cap)break;if(x!=l&&!entry->done[x])lanes[count++]=x;}
   const int W=count>likelihood_data::kLanes?MaxW:likelihood_data::kLanes;
   if(count>=2){
    auto eval_started=timed?Clock::now():Clock::time_point{};for(int i=count;i<W;++i)lanes[i]=lanes[0];
    std::array<likelihood_data::LocusView,MaxW> views;likelihood_data::LaneInput in[MaxW];
    for(int i=0;i<W;++i){views[i]=likelihood_view(lanes[i]);in[i].L=&views[i];in[i].M=&likelihood_locus(lanes[i]);in[i].T=&likelihood_table(lanes[i]);in[i].locus=lanes[i];}
    const auto& fixed=likelihood_data::lane_tables(store.likelihood_data.lane_tables,in,W);double a[MaxW],b[MaxW];
    const bool ok=entry->focal_table?likelihood_data::evaluate_focal_lanes(W,entry->plan,*likelihood_active,in,fixed,likelihood_ctx,likelihood_pairs,likelihood_pieces,int(data.size()),likelihood_lane_arena,visit->focal,a,b):likelihood_data::evaluate_lanes(W,entry->plan,*likelihood_active,in,fixed,likelihood_ctx,likelihood_pairs,likelihood_pieces,int(data.size()),likelihood_lane_arena,a);
    if(ok){const double ops=likelihood_data::program_for(entry->plan,views[0].G,views[0].k).operations;for(int i=0;i<count;++i){entry->values[lanes[i]]=a[i];if(entry->focal_table)entry->denominators[lanes[i]]=b[i];entry->done[lanes[i]]=1;++likelihood_evaluations;likelihood_operations+=ops;}++likelihood_lane_batches;}else ++likelihood_lane_fallbacks;
    const double spent=(timed?elapsed(eval_started):0.);likelihood_eval_seconds+=spent;likelihood_class_seconds[likelihood_class(key)]+=spent;
   }
  }
  if(!entry->done[l]){auto eval_started=timed?Clock::now():Clock::time_point{};const auto view=likelihood_view(l);const auto& messages=likelihood_locus(l);const auto& table=likelihood_table(l);double value=std::numeric_limits<double>::quiet_NaN();if(entry->focal_table){
    // One elimination with the focal kept: the denominator sums its table,
    // the numerator weights it by the focal's emission.
    likelihood_operations+=pap_internal::pedigree_likelihood::program_for(entry->plan,view.G,view.k).operations;double log_z=0;
    if(!pap_internal::pedigree_likelihood::evaluate_focal_table(entry->plan,*likelihood_active,view,messages,table,likelihood_ctx,likelihood_pairs,likelihood_pieces,l,int(data.size()),likelihood_arena,likelihood_focal_table,log_z,&likelihood_exact_focal_table)){value=log_z;entry->denominators[l]=log_z;}
    else if(likelihood_exact_focal_table.wide()){
     double num=-INFINITY,den=-INFINITY;const auto& row=(*view.phi)[visit->focal];
     for(int g=0;g<view.G;++g){double x=likelihood_exact_focal_table.log_at(g);den=pap_internal::pedigree_likelihood::log_add(den,x);num=pap_internal::pedigree_likelihood::log_add(num,x+row[g]);}
     entry->denominators[l]=den+log_z;value=num+log_z;
    }
    else{const auto& row=(*view.phi)[visit->focal];double top=-INFINITY;for(int g=0;g<view.G;++g)top=std::max(top,row[g]);double num=0,den=0;for(int g=0;g<view.G;++g){den+=likelihood_focal_table[g];num+=likelihood_focal_table[g]*std::exp(row[g]-top);}
     entry->denominators[l]=den>0?std::log(den)+log_z:-INFINITY;value=num>0&&std::isfinite(top)?std::log(num)+top+log_z:-INFINITY;}
   }else if(std::isnan(value)){likelihood_operations+=pap_internal::pedigree_likelihood::program_for(entry->plan,view.G,view.k).operations;value=pap_internal::pedigree_likelihood::evaluate(entry->plan,*likelihood_active,view,messages,table,likelihood_ctx,likelihood_pairs,likelihood_pieces,l,int(data.size()),likelihood_arena);}entry->values[l]=value;entry->done[l]=1;++likelihood_evaluations;const double spent=(timed?elapsed(eval_started):0.);likelihood_eval_seconds+=spent;likelihood_class_seconds[likelihood_class(key)]+=spent;if(store.config.likelihood_plan_profile){auto at=likelihood_profile_index.find(key);if(at!=likelihood_profile_index.end()){auto& row=likelihood_profile[at->second];row.seconds+=spent;row.operations+=pap_internal::pedigree_likelihood::program_for(entry->plan,view.G,view.k).operations;if(l==0)row.value0=value;}}}
  double value=entry->focal_table&&tags.term?entry->denominators[l]:entry->values[l];
  if(std::isfinite(value)){out.state=b::PartitionState::positive;out.log_partition=value;}
  else if(value==-INFINITY){out.state=b::PartitionState::zero;out.log_partition=value;}
  else {out.state=b::PartitionState::unavailable;out.reason="M_ped2 evaluation is not finite";}
  likelihood_seconds+=(timed?elapsed(started):0.);return out;
 }
 // Immutable topology ownership used by the Stage-1 catalog handoff.
 std::map<std::array<int,5>,const ::relationship_topology::TopologyModel*> engine_topologies;
 std::vector<std::shared_ptr<const b::Boundary>> engine_boundary_owners;
 std::vector<std::shared_ptr<const pap_internal::predictive_likelihood::RawBoundary>> engine_raw_owners;
 size_t calls=0;
 std::array<double,8> likelihood_availability_counts(bool optional_only=false)const {
  std::array<double,8> out{};
  if(optional_only&&likelihood_optional.extra.empty())return out;
  for(const auto& [key,entry]:likelihood_entries){
   if(likelihood_optional.optional(key)!=optional_only)continue;
   if(!entry->ready){++out[entry->state==b::PartitionState::inapplicable?1:0];continue;}
   for(size_t l=0;l<entry->done.size();++l)if(entry->done[l]){
    ++out[2];if(key[0]==1&&key[4]>=28)++out[7];const double value=entry->values[l];
    if(std::isnan(value)){++out[3];const auto& z=data[l];if(!pap_internal::pedigree_likelihood::program_for(entry->plan,z.G,z.k).ready)++out[4];}
    else if(value==-INFINITY)++out[5];else if(std::isfinite(value))++out[6];
   }
  }return out;
 }
 std::unique_ptr<pap_internal::joint_identity::Context> identity_context_;
 std::map<std::array<int,2>,double> identity_values;
 std::map<std::array<int,2>,std::array<double,2>> identity_partitions;
 std::map<std::array<int,2>,std::unique_ptr<pap_internal::joint_identity::Context>> identity_contexts;
 std::shared_ptr<const pap_internal::joint_identity::VisitFacts> identity_visit_facts;
 int identity_context_builds=0,identity_facts_builds=0;
 std::shared_ptr<const pap_internal::joint_identity::VisitFacts> identity_facts(){
  ++identity_context_builds;
  if(!input.containsElementNamed("likelihood_identity_visit_facts")||!as<bool>(input["likelihood_identity_visit_facts"])){
    ++identity_facts_builds;return {};
  }
  if(!identity_visit_facts){++identity_facts_builds;
    identity_visit_facts=std::make_shared<const pap_internal::joint_identity::VisitFacts>(*visit,identity_suite(),store.config.likelihood_ancestral,focal_role_cavity());
  }
  return identity_visit_facts;
 }
 bool configuration_identity()const{return store.config.likelihood_ancestral&&
   !(input.containsElementNamed("likelihood_identity_union")&&as<bool>(input["likelihood_identity_union"]));}
 std::array<double,6> identity_work{}; // all requests, hits, locus demands, raw reuse, support/chooser locus demands
 bool support_scoring=false;
 bool prefer_joint_support_first()const override {
  const bool s7=input.containsElementNamed("pedigree_hypothesis_suite")?
   as<std::string>(input["pedigree_hypothesis_suite"])=="S7":
   input.containsElementNamed("pedigree_hypothesis_s7")&&as<bool>(input["pedigree_hypothesis_s7"]);
  return s7&&!data.empty();
 }
 double joint_support(int dam,int sire)override {
  const bool previous=support_scoring;support_scoring=true;
  try{double value=joint_identity(dam,sire);support_scoring=previous;return value;}
  catch(...){support_scoring=previous;throw;}
 }
 std::optional<std::vector<int>> selected_identity_roots;int identity_roots_before=-1;
 bool focal_role_cavity()const{return !input.containsElementNamed("likelihood_focal_role_cavity")||as<bool>(input["likelihood_focal_role_cavity"]);}
 pap_internal::joint_identity::Context& identity_context(int dam=-1,int sire=-1){
  if(configuration_identity()&&dam>=0&&sire>=0){
   auto& context=identity_contexts[{{dam,sire}}];
   if(!context){std::vector<int> endpoints;if(dam)endpoints.push_back(dam);if(sire)endpoints.push_back(sire);
    context=std::make_unique<pap_internal::joint_identity::Context>(*visit,identity_suite(),true,&endpoints,store.config.likelihood_ancestral_reference,focal_role_cavity(),identity_facts());}
   return *context;
  }
  if(!identity_context_){
   identity_context_=std::make_unique<pap_internal::joint_identity::Context>(*visit,identity_suite(),store.config.likelihood_ancestral,selected_identity_roots?&*selected_identity_roots:nullptr,store.config.likelihood_ancestral_reference,focal_role_cavity(),identity_facts());
  }return *identity_context_;
 }
 pap_internal::hypothesis_suite::Suite identity_suite()const{
  using namespace pap_internal::hypothesis_suite;
  if(input.containsElementNamed("pedigree_hypothesis_suite"))return parse(as<std::string>(input["pedigree_hypothesis_suite"]));
  return input.containsElementNamed("pedigree_hypothesis_s7")&&as<bool>(input["pedigree_hypothesis_s7"])?Suite::S7:Suite::Legacy;
 }
 const pap_internal::joint_identity::Model& identity_model(int dam,int sire,int member,int term){
  auto& c=identity_context(dam,sire);return term==1?c.unrelated():c.models(dam,sire).at(member);
 }
 std::array<int,2> freeze_identity_roots(const std::vector<int>& roots)override{
  if(!store.config.likelihood_ancestral)throw std::logic_error("survivor roots require ancestral evidence");
  auto ordered=roots;std::sort(ordered.begin(),ordered.end());ordered.erase(std::unique(ordered.begin(),ordered.end()),ordered.end());
  // Configuration-local contexts depend only on their own endpoints. The
  // survivor union has no numerical consumer on this path. In particular,
  // pre-screened identity singles remain valid when survivors are recorded.
  const auto enabled=[&](const char* name){return input.containsElementNamed(name)&&as<bool>(input[name]);};
  if(configuration_identity()&&(enabled("likelihood_skip_union_roots")||enabled("prescreen_identity"))){
    selected_identity_roots=std::move(ordered);return {{0,0}};
  }
  if(selected_identity_roots&&*selected_identity_roots==ordered)
    return {{identity_roots_before,int(identity_context().roots.size())}};
  if(!identity_values.empty())throw std::logic_error("identity roots changed after numerical evaluation");
  if(identity_roots_before<0)identity_roots_before=int(identity_context().roots.size());
  identity_context_.reset();selected_identity_roots=std::move(ordered);
  return {{identity_roots_before,int(identity_context().roots.size())}};
 }
 bool cached_joint_identity(int dam,int sire,double& value)override {
  const auto at=identity_values.find({{dam,sire}});if(at==identity_values.end())return false;
  pap_internal::pair_rescue::identity_call(dam,sire,true);
  ++identity_work[0];++identity_work[1];likelihood_optional.identity(dam,sire);value=at->second;return true;
 }
 double joint_identity(int dam,int sire)override {
  ++identity_work[0];
  likelihood_optional.identity(dam,sire); // promote even when identity_values supplies a cache hit
  const std::array<int,2> key{{dam,sire}};
  const auto cached=identity_values.find(key);
  pap_internal::pair_rescue::identity_call(dam,sire,cached!=identity_values.end());
  if(cached!=identity_values.end()){++identity_work[1];return cached->second;}
  if(configuration_identity()&&!dam&&!sire){identity_partitions[key]={{0.,0.}};return identity_values.emplace(key,0.).first->second;}
  const auto& models=identity_context(dam,sire).models(dam,sire);double maximum=-INFINITY;bool complete=true;
  for(size_t m=0;m<models.size();++m){
   if(models[m].topology.status==::relationship_topology::TopologyStatus::Inapplicable)continue;
   if(!models[m].topology.ready()){complete=false;continue;}
   double value=0;bool applicable=true;
   for(size_t l=0;l<data.size();++l){auto* caller=pap_internal::native_demand::live;
    auto got=(caller?caller:this)->score(int(m),{5,0,visit->focal,0,0,dam,sire,0},int(l),data[l].k);
    if(got.state==b::PartitionState::inapplicable){applicable=false;break;}
    if(!got.available()){complete=false;applicable=false;break;}
    value+=got.log_partition;
   }
   if(applicable)maximum=std::max(maximum,value);
  }
  double denominator=0.;
  if(configuration_identity()&&complete&&std::any_of(models.begin(),models.end(),[](const auto& m){return m.topology.ready();})){
   for(size_t l=0;l<data.size();++l){auto* caller=pap_internal::native_demand::live;
    auto got=(caller?caller:this)->score(0,{5,0,visit->focal,0,0,dam,sire,1},int(l),data[l].k);
    if(!got.available()){complete=false;break;}denominator+=got.log_partition;
   }
  }
  // A zero numerator stays an exact zero. A zero denominator cannot certify a
  // likelihood ratio (including 0/0); retain numerical unavailability.
  if(configuration_identity()&&complete&&denominator==-INFINITY)complete=false;
  identity_partitions[key]={{maximum,denominator}};
  return identity_values.emplace(key,complete?maximum-denominator:double(NAN)).first->second;
 }
 ScorerVisit(List inputs,int focal,LiveSession& live_session)
   :input(inputs),session(live_session),data(session.input.loci),store(session.store) {
  prepared_score_batches=session.locus_batches&&inputs.containsElementNamed("likelihood_hypothesis_batches")&&as<bool>(inputs["likelihood_hypothesis_batches"]);
  visit.emplace(input,focal,session.catalog(input));
  visit->local_inputs=store.config.pedigree_likelihood;
  visit->reuse_source=store.topology_reuse;
  store.shared_revision=input.containsElementNamed("start_graph_revision")?as<int>(input["start_graph_revision"]):0;
 }

 b::PartitionValue score(int hypothesis,pap_internal::native_demand::Tags tags,int l,int)override {
  ++calls;return likelihood_score(hypothesis,tags,l);
 }


};
}
// The coordinator follows the existing native mutation/recheck lifecycle.
// One session owns shared caches; each visit owns fresh plans and scores.
namespace {
struct LiveCoordinator:pap_internal::native_demand::LiveScorer {

 std::array<double,8> pair_extra_availability{};


 void prepare_scores(const std::vector<std::pair<int,pap_internal::native_demand::Tags>>& requests)override{
  ensure_current();

  const auto started=std::chrono::steady_clock::now();
  current->prepare_scores(requests);
  score_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
 }


 double joint_identity(int dam,int sire)override{ensure_current();return current->joint_identity(dam,sire);}
 bool cached_joint_identity(int dam,int sire,double& value)override{ensure_current();return current->cached_joint_identity(dam,sire,value);}
 double joint_support(int dam,int sire)override{ensure_current();return current->joint_support(dam,sire);}
 bool prefer_joint_support_first()const override{return current&&current->prefer_joint_support_first();}
 std::array<int,2> freeze_identity_roots(const std::vector<int>& roots)override{ensure_current();return current->freeze_identity_roots(roots);}
 int single_gate_revision=-1;
 double single_gate_floor=NAN;bool single_gate_exhaustive=false,single_gate_support_only=false;
 std::map<std::array<int,4>,pap_internal::pair_rescue::Score> single_gate_scores;
 pap_internal::pair_rescue::Score single_candidate(int candidate,int side,double support_floor,bool exhaustive,bool support_only)override{
  namespace gate=pap_internal::pair_rescue;
  const int revision=as<int>(next_input["start_graph_revision"]);
  if(single_gate_revision!=revision||single_gate_floor!=support_floor||single_gate_exhaustive!=exhaustive||single_gate_support_only!=support_only){
   single_gate_scores.clear();single_gate_revision=revision;single_gate_floor=support_floor;single_gate_exhaustive=exhaustive;single_gate_support_only=support_only;
  }
  const std::array<int,4> key{{next_focal,candidate,side,int(as<bool>(next_input["recheck_only"]))}};
  if(auto at=single_gate_scores.find(key);at!=single_gate_scores.end()){auto result=at->second;result.cached=true;return result;}
  ensure_current();
  pap_internal::optional_demand::Scope optional_scope(current->likelihood_optional.prescore);
  const int dam=side==1?candidate:0,sire=side==2?candidate:0;
  gate::Score result;result.support_only=support_only;
  const bool identity_rank=current->configuration_identity()&&
    current->input.containsElementNamed("prescreen_identity")&&as<bool>(current->input["prescreen_identity"]);
  if(support_only&&current->store.config.likelihood_ancestral&&!identity_rank){
   // Retained reference: conditional focal Q, not the identity objective.
   // Pre-screen each singleton on its own ancestral question. No full suite,
   // and no numerical identity/root union before survivors have been chosen.
   double numerator=0,denominator=0,null_numerator=0,null_denominator=0;bool complete=true;
   for(size_t l=0;l<current->data.size();++l){
    auto n=score(0,{3,0,next_focal,0,0,dam,sire,0},int(l),current->data[l].k);
    auto d=score(0,{3,0,next_focal,0,0,dam,sire,1},int(l),current->data[l].k);
    auto un=score(0,{3,0,next_focal,0,0,0,0,0},int(l),current->data[l].k);
    auto ud=score(0,{3,0,next_focal,0,0,0,0,1},int(l),current->data[l].k);
    if(!n.available()||!d.available()||!un.available()||!ud.available()){complete=false;break;}
    numerator+=n.log_partition;denominator+=d.log_partition;
    null_numerator+=un.log_partition;null_denominator+=ud.log_partition;
   }
   const auto contrast=support_contrast;
   result.support=complete?contrast(contrast(numerator,denominator),contrast(null_numerator,null_denominator)):NAN;
  }else result.support=support_contrast(joint_identity(dam,sire),joint_identity(0,0));
  result.support_rejected=!std::isnan(result.support)&&result.support<support_floor;
  if(support_only){result.complete=!std::isnan(result.support);result.relationship_skipped=true;
   single_gate_scores.emplace(key,result);return result;}
  if(!exhaustive&&result.support_rejected){
   result.relationship_skipped=true;single_gate_scores.emplace(key,result);return result;
  }
  pap_internal::native_demand::Tags tags{1,0,next_focal,candidate,side,dam,sire};
  const auto boundary=catalog(tags,true);
  double parent=-INFINITY,rival=-INFINITY,parent_error=0,rival_error=0;
  bool complete=true;int parents=0,rivals=0;
  auto consume=[&](int code,bool numerator,::relationship_topology::TopologyStatus status){
    if(status==::relationship_topology::TopologyStatus::Inapplicable)return;
    if(status!=::relationship_topology::TopologyStatus::Ready){complete=false;return;}
    double value=0,error=0;bool applicable=true;
    auto request=tags;request.kind=code>=25?2:1;
    for(std::size_t l=0;l<current->data.size();++l){
      const auto cell=score(code,request,int(l),current->data[l].k);
      if(cell.state==b::PartitionState::inapplicable){applicable=false;break;}
      if(!cell.available()){complete=false;applicable=false;break;}
      if((!std::isfinite(cell.log_partition)&&cell.log_partition!=-INFINITY)||
         !std::isfinite(cell.reference_error)||cell.reference_error<0){complete=false;applicable=false;break;}
      value+=cell.log_partition;error+=cell.reference_error;
    }
    if(!applicable)return;
    if(std::isnan(value)||value==INFINITY||!std::isfinite(error)){complete=false;return;}
    if(numerator){++parents;parent=std::max(parent,value);parent_error=std::max(parent_error,error);}
    else{++rivals;rival=std::max(rival,value);rival_error=std::max(rival_error,error);}
  };
  for(std::size_t code=0;code<boundary->models.size();++code)
    consume(int(code),code==0,boundary->models[code].topology.status);
  if(boundary->full_boundary)for(const auto& member:boundary->full_boundary->members)
    consume(member.model_order,member.parent_numerator,member.topology_status);
  result.relationship=support_contrast(parent,rival);
  result.error=parent_error+rival_error;
  result.complete=complete&&parents&&rivals&&!std::isnan(result.support)&&!std::isnan(result.relationship);
  single_gate_scores.emplace(key,result);return result;
 }



 int diagnostic_round=0;

 bool engine_input_views=false;size_t input_views=0;
 List input,next_input;int next_focal=0;std::unique_ptr<LiveSession> session;std::unique_ptr<ScorerVisit> current;size_t owner_builds=0,closed_visits=0;::relationship_topology::FinishedModelTransfer model_transfer;bool pipeline_reuse=false;size_t pedigree_updates=0;double pedigree_maintenance_seconds=0;double cleanup_seconds=0,owner_prepare_seconds=0,score_seconds=0,current_prepare=0;size_t likelihood_cyclic_plans=0,likelihood_rejected_cyclic_plans=0;std::array<double,5> likelihood_scope{};std::array<double,10> likelihood_totals{};std::array<double,8> likelihood_availability{};std::array<double,6> identity_work{};
 // Passive trace: observe completed native scores after workers have joined.
 // No topology construction, scoring request, or gate is triggered here.


 LiveCoordinator(List c,List i):input(i){auto option=[&](const char* n){return i.containsElementNamed(n)&&as<bool>(i[n]);};prepared_score_batches=option("likelihood_batch_loci")&&option("likelihood_hypothesis_batches");hypothesis_s7=option("pedigree_hypothesis_s7");hypothesis_suite=hypothesis_s7?pap_internal::hypothesis_suite::Suite::S7:pap_internal::hypothesis_suite::Suite::Legacy;if(i.containsElementNamed("pedigree_hypothesis_suite")){hypothesis_suite=pap_internal::hypothesis_suite::parse(as<std::string>(i["pedigree_hypothesis_suite"]));if(i.containsElementNamed("pedigree_hypothesis_exclude")&&Rf_length(i["pedigree_hypothesis_exclude"]))throw std::invalid_argument("named hypothesis suites have fixed membership; do not combine with mped_hypothesis_exclude");}if(i.containsElementNamed("pedigree_hypothesis_exclude"))for(const auto& name:as<std::vector<std::string>>(i["pedigree_hypothesis_exclude"]))hypothesis_exclude.insert(name);exact_step_scaling=option("pedigree_exact_step_scaling");exterior_work_planning=option("pedigree_exterior_work_planning");exterior_state_products=option("pedigree_exterior_state_products");native_catalog=option("pedigree_native_catalog");native_support=native_catalog||option("pedigree_native_support");authority_audit=option("pedigree_authority_audit");skip_c_admissions=native_support||option("pedigree_skip_engine_admissions");pipeline_reuse=option("pedigree_pipeline_reuse");engine_descriptor=i.containsElementNamed("pedigree_engine_descriptor")&&as<bool>(i["pedigree_engine_descriptor"]);engine_topology_reuse=i.containsElementNamed("pedigree_engine_topology_reuse")&&as<bool>(i["pedigree_engine_topology_reuse"]);engine_life_cache.reuse_builder=true;engine_life_cache.flat_topology=i.containsElementNamed("pedigree_engine_flat_topology")&&as<bool>(i["pedigree_engine_flat_topology"]);engine_life_cache.singleton_life=i.containsElementNamed("pedigree_engine_life_certificates")&&as<bool>(i["pedigree_engine_life_certificates"]);construction_budget(i);motif_template_reuse=i.containsElementNamed("pedigree_motif_templates")&&as<bool>(i["pedigree_motif_templates"]);
 engine_input_views=i.containsElementNamed("pedigree_engine_input_views")&&as<bool>(i["pedigree_engine_input_views"]);engine_move_results=i.containsElementNamed("pedigree_engine_move_results")&&as<bool>(i["pedigree_engine_move_results"]);engine_reduction_reuse=i.containsElementNamed("pedigree_engine_reduction_reuse")&&as<bool>(i["pedigree_engine_reduction_reuse"]);engine_compact_results=i.containsElementNamed("pedigree_engine_compact_results")&&as<bool>(i["pedigree_engine_compact_results"]);engine_result_reset=i.containsElementNamed("pedigree_engine_result_reset")&&as<bool>(i["pedigree_engine_result_reset"]);
 session=std::make_unique<LiveSession>(c,i);}

 std::shared_ptr<const pap_internal::ExecutionConfig> execution_config()const override{return session->input.config;}

 void close_visit(){


  native_catalogs.clear();observed_catalogs.clear();if(!current)return;++closed_visits;
  // Count authorized hypotheses only, including aliases as separate scientific
  // questions. Unavailable numerical evidence is distinct from a structural zero.
  auto availability=current->likelihood_availability_counts();
  const auto optional=current->likelihood_availability_counts(true);
  for(size_t j=0;j<optional.size();++j)pair_extra_availability[j]+=optional[j];

  for(size_t j=0;j<availability.size();++j)likelihood_availability[j]+=availability[j];
  for(std::size_t j=0;j<identity_work.size();++j)identity_work[j]+=current->identity_work[j];
  likelihood_cyclic_plans+=current->likelihood_cyclic_plans;likelihood_rejected_cyclic_plans+=current->likelihood_rejected_cyclic_plans;
  for(size_t j=0;j<likelihood_scope.size();++j)likelihood_scope[j]+=current->likelihood_scope[j];
  likelihood_totals[0]+=current->likelihood_build_seconds;likelihood_totals[1]+=current->likelihood_topology_seconds;likelihood_totals[2]+=current->likelihood_forest_seconds;likelihood_totals[3]+=current->likelihood_message_seconds;likelihood_totals[4]+=current->likelihood_eval_seconds;likelihood_totals[5]+=current->likelihood_operations;likelihood_totals[6]+=current->likelihood_lane_batches;likelihood_totals[7]+=current->likelihood_lane_fallbacks;likelihood_totals[8]+=current->likelihood_questions;likelihood_totals[9]+=current->likelihood_questions_cut;

  auto t=Clock::now();{current.reset();}
  double cleanup=elapsed(t);cleanup_seconds+=cleanup;
 }


 bool accepted_move(int focal,int,int,int nd,int ns,const std::vector<int>&,const std::vector<int>&,int revision)override{
  close_visit();if(session)session->store.shared_revision=revision;
  if(auto* cache=pap_internal::relative_life::current)
   if(focal>0&&size_t(focal)<=cache->source_parents.size())cache->source_parents[focal-1]={{nd,ns}};
  if(auto* cache=::relationship_topology::TopologyTemplateCache::active)cache->invalidate(focal);
  return false;
 }




 void begin_visit(int focal,const std::vector<int>& dam,const std::vector<int>& sire,int revision,bool recheck,const List& active_inputs)override{

  close_visit();if(pipeline_reuse)model_transfer.clear();diagnostic_round=active_inputs.containsElementNamed("progress_parentage_round")?as<int>(active_inputs["progress_parentage_round"]):1;if(session)session->store.likelihood_round=diagnostic_round;next_input=[&]{return engine_input_views?List(Rf_shallow_duplicate(active_inputs)):clone(active_inputs);}();input_views+=engine_input_views;IntegerMatrix parents(dam.size()-1,2);for(int v=1;v<int(dam.size());++v){parents(v-1,0)=dam[v];parents(v-1,1)=sire[v];}next_input["parents"]=parents;next_input["start_graph_revision"]=revision;next_input["recheck_only"]=recheck;next_focal=focal;


 }


 // Catalogs have visit ownership. They share the existing pedigree likelihood builder cache,
 // not a C-scope source or a cross-revision topology template.
  using CatalogKey=std::array<int,4>;
 std::map<CatalogKey,std::shared_ptr<const b::Boundary>> native_catalogs;
 // S7 adds reverse PO to ordinary25; compound full variants are not built.
 bool hypothesis_s7=false;pap_internal::hypothesis_suite::Suite hypothesis_suite=pap_internal::hypothesis_suite::Suite::Legacy;
 // Opt-in: full-suite relationships left out by name ("HS_FA"), or only their
 // reverse orientation ("HS_PO:reverse").
 std::set<std::string> hypothesis_exclude;
 std::set<std::array<int,7>> observed_catalogs;
 struct CatalogEvent{int round,revision;pap_internal::native_demand::AuthorityRow row;};
 struct GateEvent{int round,revision;pap_internal::native_demand::GateRow row;};
 std::vector<CatalogEvent> catalog_events;std::vector<GateEvent> gate_events;
 size_t catalog_builds=0,catalog_hits=0,raw_catalog_builds=0;
 void ensure_current(){
  if(current)return;

  auto started=Clock::now();
  current=std::make_unique<ScorerVisit>(next_input,next_focal,*session);++owner_builds;
  current_prepare=elapsed(started);owner_prepare_seconds+=current_prepare;
 }
 std::shared_ptr<const b::Boundary> catalog(pap_internal::native_demand::Tags tags,bool full)override{

  ensure_current();int opposite=tags.side==1?tags.sire:tags.dam;
  CatalogKey key{{tags.candidate,tags.side,opposite,int(full)}};
  if(auto at=native_catalogs.find(key);at!=native_catalogs.end()){++catalog_hits;return at->second;}
  ++catalog_builds;auto& owner=*current;
  pap_internal::topology_life_reuse::Scope life_scope(owner.store.topology_reuse?&owner.store.life_certificates:nullptr);
  auto& visit=*owner.visit;auto& context=owner.contexts[{{tags.candidate,tags.side,opposite}}];
  if(!context)context=std::make_unique<pap_internal::candidate_catalog::Visit::Context>(visit,tags.candidate,tags.side,opposite);
  auto out=std::make_shared<b::Boundary>();out->native_support=true;
  for(const auto& node:context->input.nodes)if(node.observed&&node.typed_id>0)out->evidence_ids.push_back(node.typed_id);
  std::sort(out->evidence_ids.begin(),out->evidence_ids.end());
  auto metadata=[&](int code){
   auto at=context->built.find(code);if(at==context->built.end()){at=context->built.emplace(code,visit.build(context->input,code,context->source.get())).first;visit.source_copies+=bool(context->source);}
   b::Model m;m.topology=at->second;m.supported=m.topology.ready();m.reason=m.topology.reason;return m;
  };
  for(int code=0;code<25;++code)out->models.push_back(metadata(code));
  if(full){auto extra=std::make_shared<b::FullBoundary>();extra->compiled.native_support=true;extra->compiled.evidence_ids=out->evidence_ids;
   extra->tested_side=tags.side;extra->focal_typed_id=tags.focal;extra->candidate_typed_id=tags.candidate;
   auto append=[&](int code,const ::relationship_topology::full_catalog::Definition& def,bool reverse,unsigned bits){
    if(hypothesis_exclude.count(def.relationship)||(reverse&&hypothesis_exclude.count(def.relationship+":reverse")))return;auto m=metadata(code);b::FullMemberDescriptor d;
    d.model_order=code;d.name=m.topology.name;d.relationship=def.relationship;d.reverse=reverse;d.parent_slot_bits=bits;
    d.parent_numerator=def.parent_numerator&&!reverse;d.ownership_review_required=def.parent_numerator&&reverse;
    d.orientation_review_required=def.reverse_review_required&&reverse;d.registered_endpoint_orientation=def.registered_class&&!d.orientation_review_required;
    d.topology_status=m.topology.status;d.reason=m.reason;
    if(m.topology.ready()){d.compiled_index=extra->compiled.models.size();extra->compiled_to_member.push_back(extra->members.size());extra->compiled.models.push_back(std::move(m));}
    extra->members.push_back(std::move(d));
   };
   if(hypothesis_suite==pap_internal::hypothesis_suite::Suite::S7)
    for(int code:{25,26})append(code,pap_internal::candidate_catalog::reverse_po_definition(),true,unsigned(code-25));
   for(const auto& [code,v]:visit.variants)if(pap_internal::hypothesis_suite::contains(hypothesis_suite,code))
    append(code,visit.full.at(v.definition),v.reverse,v.bits);
   out->full_boundary=std::move(extra);
  }
  auto result=native_catalogs.emplace(key,std::move(out)).first->second;
  if(pap_internal::input_checks::share()){
   owner.engine_boundary_owners.push_back(result);
   auto bind=[&](const b::Model& m){owner.engine_topologies[{{1,tags.candidate,tags.side,opposite,m.topology.model_order}}]=&m.topology;};
   for(const auto& m:result->models)bind(m);if(result->full_boundary)for(const auto& m:result->full_boundary->compiled.models)bind(m);
  }return result;
 }
 std::shared_ptr<const pap_internal::predictive_likelihood::RawBoundary> raw_catalog(int focal,const std::vector<std::array<int,2>>& proposals,const std::vector<unsigned char>& needed)override{

  ensure_current();++raw_catalog_builds;auto& owner=*current;
  pap_internal::topology_life_reuse::Scope life_scope(owner.store.topology_reuse?&owner.store.life_certificates:nullptr);
  auto& visit=*owner.visit;if(!owner.raw)owner.raw=std::make_unique<pap_internal::candidate_catalog::Visit::Raw>(visit);
  auto out=std::make_shared<pap_internal::predictive_likelihood::RawBoundary>();out->compiled.native_support=true;
  // Engine mapping metadata only. All evidence remains in pedigree likelihood's actual raw
  // boundary and its original numerator/denominator, consumed by score().
  out->common.input.nodes.resize(visit.parents.nrow());for(int v=0;v<visit.parents.nrow();++v)out->common.input.nodes[v].typed_id=v+1;
  out->common.input.focal=focal-1;out->common.ready=true;
  for(size_t i=0;i<proposals.size();++i){auto proposal=proposals[i];
   pap_internal::predictive_likelihood::Request r;r.dam=proposal[0]-1;r.sire=proposal[1]-1;
   if(!needed.at(i)){r.reason="Not requested: independent pair gate proves ineligibility; no raw score dependency";out->requests.push_back(std::move(r));continue;}
   const auto& pair=visit.raw_pair(*owner.raw,proposal[0],proposal[1]);const auto& n=pair.numerator;const auto& d=pair.denominator;
   r.build_status=n.ready()&&d.ready()?pap_internal::predictive_likelihood::Request::BuildStatus::ready:
    n.status==::relationship_topology::TopologyStatus::Inapplicable||d.status==::relationship_topology::TopologyStatus::Inapplicable?pap_internal::predictive_likelihood::Request::BuildStatus::hard_invalid : pap_internal::predictive_likelihood::Request::BuildStatus::unavailable;
   r.reason=!n.ready()?n.reason:d.reason;
   r.numerator=out->compiled.models.size();b::Model mn;mn.topology=n;mn.supported=n.ready();out->compiled.models.push_back(std::move(mn));
   r.denominator=out->compiled.models.size();b::Model md;md.topology=d;md.supported=d.ready();out->compiled.models.push_back(std::move(md));
   out->requests.push_back(std::move(r));
  }
  if(pap_internal::input_checks::share()){owner.engine_raw_owners.push_back(out);for(size_t i=0;i<proposals.size();++i){const auto& r=out->requests[i];if(!r.ready())continue;for(int term=0;term<2;++term)owner.engine_topologies[{{3,proposals[i][0],proposals[i][1],term,0}}]=&out->compiled.models[term?r.denominator:r.numerator].topology;}}
  return out;
 }
 void observe_raw_catalog(int focal,const pap_internal::predictive_likelihood::RawBoundary& boundary)override{
  if(!authority_audit)return;int revision=as<int>(next_input["start_graph_revision"]);
  auto typed=[&](int v){return v>=0&&size_t(v)<boundary.common.input.nodes.size()?boundary.common.input.nodes[v].typed_id:0;};
  for(const auto& r:boundary.requests)for(int term=0;term<2;++term){size_t index=term?r.denominator:r.numerator;
   const b::Model* m=index<boundary.compiled.models.size()?&boundary.compiled.models[index]:nullptr;
   int status=m?int(m->topology.status):int(r.build_status==pap_internal::predictive_likelihood::Request::BuildStatus::hard_invalid?::relationship_topology::TopologyStatus::Inapplicable : ::relationship_topology::TopologyStatus::Unsupported);
   catalog_events.push_back({diagnostic_round,revision,{{3,0,focal,0,0,typed(r.dam),typed(r.sire),term},term,status,m&&m->supported,term?"raw_D":"raw_N","raw",r.reason}});
  }
 }
 void observe_catalog(pap_internal::native_demand::Tags tags,const b::Boundary& boundary)override{
  if(!authority_audit)return;int revision=as<int>(next_input["start_graph_revision"]);
  if(!observed_catalogs.insert(std::array<int,7>{{diagnostic_round,tags.focal,revision,tags.candidate,tags.side,tags.dam,tags.sire}}).second)return;
  auto add=[&](int code,const std::string& name,const std::string& relationship,int status,bool supported,const std::string& reason){catalog_events.push_back({diagnostic_round,revision,{tags,code,status,supported,name,relationship,reason}});};
  for(const auto& m:boundary.models)add(m.topology.model_order,m.topology.name,m.topology.relationship,int(m.topology.status),m.supported,m.reason);
  if(boundary.full_boundary)for(const auto& d:boundary.full_boundary->members){bool supported=d.compiled_index<boundary.full_boundary->compiled.models.size()&&boundary.full_boundary->compiled.models[d.compiled_index].supported;add(d.model_order,d.name,d.relationship,int(d.topology_status),supported,d.reason);}
 }
 void observe_gate(pap_internal::native_demand::Tags tags,int code,int locus,int k,const b::Model& m,const b::PartitionValue& value)override{
  if(!authority_audit)return;
  bool old_applicable=m.topology.status!=::relationship_topology::TopologyStatus::Inapplicable,new_applicable=value.state!=b::PartitionState::inapplicable;
  bool old_available=m.topology.ready()&&m.supported,new_available=new_applicable&&value.available();
  if(old_applicable==new_applicable&&old_available==new_available)return;
  gate_events.push_back({diagnostic_round,as<int>(next_input["start_graph_revision"]),{tags,code,locus,k,int(m.topology.status),int(value.state),m.supported,value.log_partition,value.reference_error,m.reason,value.reason}});
 }
 List authority_receipt(){
  std::vector<int> round,revision,kind,focal,candidate,side,dam,sire,code,status;std::vector<bool> supported;std::vector<std::string> name,relationship,reason;
  for(const auto& e:catalog_events){const auto& r=e.row;round.push_back(e.round);revision.push_back(e.revision);kind.push_back(r.tags.kind);focal.push_back(r.tags.focal);candidate.push_back(r.tags.candidate);side.push_back(r.tags.side);dam.push_back(r.tags.dam);sire.push_back(r.tags.sire);code.push_back(r.code);status.push_back(r.status);supported.push_back(r.supported);name.push_back(r.name);relationship.push_back(r.relationship);reason.push_back(r.reason);}
  List out=List::create(_["catalog_builds"]=double(catalog_builds),_["catalog_hits"]=double(catalog_hits),_["raw_catalog_builds"]=double(raw_catalog_builds));
  if(authority_audit){out["catalog"]=DataFrame::create(_["round"]=round,_["revision"]=revision,_["kind"]=kind,_["focal"]=focal,_["candidate"]=candidate,_["side"]=side,_["dam"]=dam,_["sire"]=sire,_["code"]=code,_["status"]=status,_["supported"]=supported,_["name"]=name,_["relationship"]=relationship,_["reason"]=reason);
   std::vector<int> gr,gv,gf,gc,gs,gd,gt,gcode,gl,gk,old_status,state;std::vector<bool> old_supported;std::vector<double> value,error;std::vector<std::string> old_reason,new_reason;
   for(const auto& e:gate_events){const auto& r=e.row;gr.push_back(e.round);gv.push_back(e.revision);gf.push_back(r.tags.focal);gc.push_back(r.tags.candidate);gs.push_back(r.tags.side);gd.push_back(r.tags.dam);gt.push_back(r.tags.sire);gcode.push_back(r.code);gl.push_back(r.locus);gk.push_back(r.k);old_status.push_back(r.status);state.push_back(r.state);old_supported.push_back(r.supported);value.push_back(r.value);error.push_back(r.error);old_reason.push_back(r.old_reason);new_reason.push_back(r.reason);}
   out["gates"]=DataFrame::create(_["round"]=gr,_["revision"]=gv,_["focal"]=gf,_["candidate"]=gc,_["side"]=gs,_["dam"]=gd,_["sire"]=gt,_["code"]=gcode,_["locus"]=gl,_["k"]=gk,_["C_status"]=old_status,_["C_supported"]=old_supported,_["pedigree_state"]=state,_["log_partition"]=value,_["reference_error"]=error,_["C_reason"]=old_reason,_["pedigree_reason"]=new_reason);
  }return out;
 }
 b::PartitionValue score(int h,pap_internal::native_demand::Tags tags,int l,int k)override{

  ensure_current();
  auto t=Clock::now();auto out=current->score(h,tags,l,k);score_seconds+=elapsed(t);

  return out;
 }

};
}
List search_rounds_cpp(SEXP compiled,List inputs);
List visit_candidates_cpp(SEXP compiled,List inputs);


// [[Rcpp::export]]
List run_parentage_cpp(SEXP compiled,List compact,List inputs){
 InputCheckScope internal_checks(inputs);
 bool relative_enabled=inputs.containsElementNamed("pedigree_relative_life_cache")&&as<bool>(inputs["pedigree_relative_life_cache"]);
 pap_internal::relative_life::Cache relative_cache;
 if(relative_enabled){IntegerMatrix parents=inputs["parents"];relative_cache.source_parents.resize(parents.nrow());for(int v=0;v<parents.nrow();++v)relative_cache.source_parents[v]={{parents(v,0),parents(v,1)}};}
 pap_internal::relative_life::Scope relative_scope(relative_enabled?&relative_cache:nullptr);
 bool templates_enabled=inputs.containsElementNamed("pedigree_topology_templates")&&as<bool>(inputs["pedigree_topology_templates"]),bound_enabled=inputs.containsElementNamed("pedigree_bound_topologies")&&as<bool>(inputs["pedigree_bound_topologies"]);
 ::relationship_topology::TopologyTemplateCache topology_templates;topology_templates.shape_enabled=templates_enabled;topology_templates.bound_enabled=bound_enabled;
 ::relationship_topology::TopologyTemplateCache::Scope topology_scope(templates_enabled||bound_enabled?&topology_templates:nullptr);


 auto started=Clock::now();LiveCoordinator owner(compact,inputs);::relationship_topology::FinishedModelTransfer::Scope transfer_scope(owner.pipeline_reuse?&owner.model_transfer:nullptr);pap_internal::native_demand::LiveScope scope(&owner);
 List result;{const bool frozen=inputs.containsElementNamed("private_frozen_decision_audit")&&as<bool>(inputs["private_frozen_decision_audit"]);
 result=frozen?visit_candidates_cpp(compiled,inputs):search_rounds_cpp(compiled,inputs);owner.close_visit();
 if(frozen){result["frozen_state_audit"]=true;
  IntegerMatrix before=inputs["parents"],after=result["final_parents"];
  if(before.nrow()!=after.nrow()||before.ncol()!=after.ncol()||!std::equal(before.begin(),before.end(),after.begin()))
    throw std::logic_error("frozen decision audit changed the parent matrix");
  std::vector<std::array<int,4>> keys;
  if(owner.session&&owner.session->store.likelihood_data.forest){const auto& f=*owner.session->store.likelihood_data.forest;
   for(const auto& c:f.cuts)keys.push_back({c[0],c[2],f.families[c[1]].dam,f.families[c[1]].sire});}
  std::sort(keys.begin(),keys.end());IntegerMatrix cuts(keys.size(),4);
  for(std::size_t i=0;i<keys.size();++i)for(int j=0;j<4;++j)cuts(i,j)=keys[i][j];result["likelihood_forest_cuts"]=cuts;
 }}
 if(owner.session)result["scorer_setup_reuse"]=List::create(_["parse_builds"]=double(owner.session->parse_builds),_["parse_hits"]=double(owner.session->parse_hits));


 result["pedigree_visits"]=List();
 result["pedigree_hypothesis_suite"]=pap_internal::hypothesis_suite::name(owner.hypothesis_suite);
 result["pedigree_hypothesis_codes"]=wrap(pap_internal::hypothesis_suite::codes(owner.hypothesis_suite));
 result["likelihood_availability"]=NumericVector::create(_["topology_unavailable_models"]=owner.likelihood_availability[0],_["inapplicable_models"]=owner.likelihood_availability[1],_["evaluated_hypothesis_locus_cells"]=owner.likelihood_availability[2],_["unavailable_hypothesis_locus_cells"]=owner.likelihood_availability[3],_["program_refused_hypothesis_locus_cells"]=owner.likelihood_availability[4],_["structural_zero_hypothesis_locus_cells"]=owner.likelihood_availability[5],_["finite_hypothesis_locus_cells"]=owner.likelihood_availability[6],_["compound_hypothesis_locus_cells"]=owner.likelihood_availability[7]);
 if((inputs.containsElementNamed("pair_rescue")&&as<bool>(inputs["pair_rescue"]))||(inputs.containsElementNamed("pair_rescue_shadow")&&as<bool>(inputs["pair_rescue_shadow"]))) {
  NumericVector extra(owner.pair_extra_availability.begin(),owner.pair_extra_availability.end());
  extra.attr("names")=as<NumericVector>(result["likelihood_availability"]).attr("names");result["likelihood_optional_prescore_availability"]=extra;
 }
 result["likelihood_program_registry"]=program_registry_receipt();
 result["likelihood_ancestral_scope"]=NumericVector::create(_["off_scope_observation_ports_before"]=owner.likelihood_scope[0],_["descendant_observation_ports_before"]=owner.likelihood_scope[1],_["off_scope_observation_ports_after"]=0.,_["descendant_observation_ports_after"]=0.,_["barren_variables_removed"]=owner.likelihood_scope[2],_["off_scope_messages_replaced"]=owner.likelihood_scope[3],_["questions_with_off_scope_observations"]=owner.likelihood_scope[4]);
 result["likelihood_identity_work"]=NumericVector::create(_["requests"]=owner.identity_work[0],_["cache_hits"]=owner.identity_work[1],_["kind5_locus_requests"]=owner.identity_work[2],_["raw_reuses"]=owner.identity_work[3],_["support_locus_requests"]=owner.identity_work[4],_["chooser_locus_requests"]=owner.identity_work[5]);
 {const auto& a=pap_internal::pedigree_likelihood::resource_plan::search_audit();result["likelihood_plan_search"]=NumericVector::create(_["keys"]=double(a.keys),_["capped_keys"]=double(a.capped_keys),_["capped_refused"]=double(a.capped_refused),_["attempts"]=double(a.attempts),_["admitted"]=double(a.admitted),_["attempt_limit"]=double(pap_internal::pedigree_likelihood::resource_plan::SearchBudget::default_limit));}
 result["likelihood_plan_keys"]=NumericVector::create(_["cyclic_plans"]=double(owner.likelihood_cyclic_plans),_["rejected_cyclic_plans"]=double(owner.likelihood_rejected_cyclic_plans));
 result["likelihood_work"]=NumericVector::create(_["model_build_inclusive_seconds"]=owner.likelihood_totals[0],_["topology_seconds"]=owner.likelihood_totals[1],_["forest_seconds"]=owner.likelihood_totals[2],_["message_initialization_seconds"]=owner.likelihood_totals[3],_["evaluation_inclusive_seconds"]=owner.likelihood_totals[4],_["operations"]=owner.likelihood_totals[5],_["lane_batches"]=owner.likelihood_totals[6],_["lane_fallbacks"]=owner.likelihood_totals[7],_["questions"]=owner.likelihood_totals[8],_["questions_with_ancestry_cut"]=owner.likelihood_totals[9]);
 if(!owner.session->store.config.likelihood_time_evaluations){NumericVector times=result["likelihood_work"];times["evaluation_inclusive_seconds"]=NA_REAL;result["likelihood_work"]=times;}
 result["pedigree_owner_prepare_seconds"]=owner.owner_prepare_seconds;result["pedigree_owner_cleanup_seconds"]=owner.cleanup_seconds;result["pedigree_score_seconds"]=owner.score_seconds;result["pedigree_seconds"]=elapsed(started);


 result["topology_templates"]=List::create(_["shape_requests"]=double(topology_templates.shape_requests),_["shape_hits"]=double(topology_templates.shape_hits),_["shape_builds"]=double(topology_templates.shape_builds),_["bound_requests"]=double(topology_templates.bound_requests),_["bound_hits"]=double(topology_templates.bound_hits),_["invalidations"]=double(topology_templates.invalidations),_["closure_misses"]=double(topology_templates.closure_misses),_["refusals"]=double(topology_templates.refusals),_["bytes"]=double(topology_templates.bytes),_["lookup_seconds"]=topology_templates.lookup_seconds);
 if(owner.session&&owner.session->locus_batches){const auto& p=*owner.session->locus_batches;result["likelihood_narrowed_batches"]=double(p.narrowed_batches());{IntegerVector alleles;CharacterVector reason;NumericVector cells;for(const auto& [key,n]:p.refusal_counts()){alleles.push_back(key.first);reason.push_back(key.second);cells.push_back(double(n));}{const auto& a=pap_internal::pedigree_likelihood::resource_audit();result["likelihood_resource_audit"]=List::create(_["captured"]=double(a.captured),_["capture_seconds"]=a.seconds);}result["likelihood_parallel_refusals"]=DataFrame::create(_["alleles"]=alleles,_["reason"]=reason,_["evaluated_cells"]=cells);if(p.resource_diagnostics){
 IntegerVector kk,step,dim,inputs,gg,aa,select;CharacterVector source,reason,shape,axes;NumericVector cells,arena,work,total,count;
 for(const auto& [key,d]:p.resource_refusals()){const auto& q=*d.program;
 kk.push_back(q.k);step.push_back(q.resource_step);dim.push_back(q.eliminated_dim);inputs.push_back(q.input_count);gg.push_back(q.genotype_axes);aa.push_back(q.allele_axes);select.push_back(q.selection);
 source.push_back(key.first);reason.push_back(q.refusal?q.refusal:"nonfinite value");
 std::string text;for(int v:key.second){if(!text.empty())text+=",";text+=std::to_string(v);}shape.push_back(text);
 text.clear();for(int v:q.output_dims){if(!text.empty())text+="x";text+=std::to_string(v);}axes.push_back(text);
 cells.push_back(q.required_cells);arena.push_back(q.required_arena);work.push_back(q.step_work);total.push_back(q.required_operations);count.push_back(double(d.count));
 }
 result["likelihood_resource_refusals"]=DataFrame::create(_["alleles"]=kk,_["source"]=source,_["reason"]=reason,_["shape"]=shape,_["step"]=step,_["eliminated_dim"]=dim,_["input_count"]=inputs,_["genotype_axes"]=gg,_["allele_axes"]=aa,_["selection"]=select,_["output_dims"]=axes,_["table_cells"]=cells,_["arena_cells"]=arena,_["step_operations"]=work,_["program_operations_through_step"]=total,_["evaluated_cells"]=count);
}}{const auto c=p.cache_counts();const auto s=p.scheduling_counts();result["likelihood_scheduler_work"]=NumericVector::create(_["cache_owners"]=double(s[0]),_["owner_claims"]=double(s[1]),_["owner_migrations"]=double(s[2]),_["lane_batches"]=double(s[3]),_["lane_fallbacks"]=double(s[4]),_["active_lane_cells"]=double(s[5]),_["allocated_lane_cells"]=double(s[6]),_["pair_slot_copies"]=double(s[7]),_["piece_slot_copies"]=double(s[8]),_["pair_requests"]=double(c[0]),_["pair_computed"]=double(c[1]),_["piece_requests"]=double(c[2]),_["piece_computed"]=double(c[3]),_["context_repairs"]=double(c[4]));}result["likelihood_parallel"]=List::create(_["scheduler"]=p.scheduler_name(),_["workers"]=p.worker_count(),_["jobs"]=double(p.jobs),_["parallel_jobs"]=double(p.parallel_jobs),_["serial_jobs"]=double(p.serial_jobs),_["prepared_plans"]=double(p.prepared_plans),_["prepared_pieces"]=double(p.prepared_pieces),_["preparation_seconds"]=p.preparation_seconds,_["execution_seconds"]=p.execution_seconds,_["message_seconds"]=p.message_seconds,_["worker_profile"]=p.profile,_["worker_busy_seconds"]=p.profile?p.worker_busy_seconds:NA_REAL,_["worker_critical_seconds"]=p.profile?p.worker_critical_seconds:NA_REAL,_["dispatch_gap_seconds"]=p.profile?p.dispatch_gap_seconds:NA_REAL,_["imbalance_worker_seconds"]=p.profile?p.imbalance_seconds:NA_REAL);}if(owner.session&&owner.session->store.config.likelihood_latent_brood_exclusion){const auto& c=owner.session->store.likelihood_brood_checks;
 if(c[1]||c[2])throw std::logic_error(owner.session->store.likelihood_brood_failure);
 result["likelihood_brood_exclusion"]=NumericVector::create(_["paired_program_domains"]=c[0],_["width_increases"]=c[1],_["new_program_refusals"]=c[2]);}
result["pedigree_visits"]=List();result["pedigree_hypothesis_suite"]=pap_internal::hypothesis_suite::name(owner.hypothesis_suite);result["pedigree_hypothesis_codes"]=wrap(pap_internal::hypothesis_suite::codes(owner.hypothesis_suite));result["likelihood_availability"]=NumericVector::create(_["topology_unavailable_models"]=owner.likelihood_availability[0],_["inapplicable_models"]=owner.likelihood_availability[1],_["evaluated_hypothesis_locus_cells"]=owner.likelihood_availability[2],_["unavailable_hypothesis_locus_cells"]=owner.likelihood_availability[3],_["program_refused_hypothesis_locus_cells"]=owner.likelihood_availability[4],_["structural_zero_hypothesis_locus_cells"]=owner.likelihood_availability[5],_["finite_hypothesis_locus_cells"]=owner.likelihood_availability[6],_["compound_hypothesis_locus_cells"]=owner.likelihood_availability[7]);if((inputs.containsElementNamed("pair_rescue")&&as<bool>(inputs["pair_rescue"]))||
   (inputs.containsElementNamed("pair_rescue_shadow")&&as<bool>(inputs["pair_rescue_shadow"]))) {
 NumericVector extra(owner.pair_extra_availability.begin(),owner.pair_extra_availability.end());
 extra.attr("names")=as<NumericVector>(result["likelihood_availability"]).attr("names");
 result["likelihood_optional_prescore_availability"]=extra;
}

if(owner.session)result["likelihood_program_registry"]=program_registry_receipt();if(owner.session){const auto& c=owner.session->store.likelihood_data.lane_tables;
result["likelihood_lane_table_cache"]=NumericVector::create(_["entries"]=double(c.values.size()),_["bounded"]=c.bounded?1.:0.,_["charged_bytes"]=double(c.bytes),_["peak_charged_bytes"]=double(c.peak_bytes),_["limit_bytes"]=double(c.limit),_["hits"]=double(c.hits),_["misses"]=double(c.misses),_["evictions"]=double(c.evictions),_["declines"]=double(c.declines));}

result["likelihood_ancestral_scope"]=NumericVector::create(_["off_scope_observation_ports_before"]=owner.likelihood_scope[0],_["descendant_observation_ports_before"]=owner.likelihood_scope[1],_["off_scope_observation_ports_after"]=0.,_["descendant_observation_ports_after"]=0.,_["barren_variables_removed"]=owner.likelihood_scope[2],_["off_scope_messages_replaced"]=owner.likelihood_scope[3],_["questions_with_off_scope_observations"]=owner.likelihood_scope[4]);result["likelihood_identity_work"]=NumericVector::create(_["requests"]=owner.identity_work[0],_["cache_hits"]=owner.identity_work[1],_["kind5_locus_requests"]=owner.identity_work[2],_["raw_reuses"]=owner.identity_work[3],_["support_locus_requests"]=owner.identity_work[4],_["chooser_locus_requests"]=owner.identity_work[5]);{const auto& a=pap_internal::pedigree_likelihood::resource_plan::search_audit();result["likelihood_plan_search"]=NumericVector::create(_["keys"]=double(a.keys),_["capped_keys"]=double(a.capped_keys),_["capped_refused"]=double(a.capped_refused),_["attempts"]=double(a.attempts),_["admitted"]=double(a.admitted),_["attempt_limit"]=double(pap_internal::pedigree_likelihood::resource_plan::SearchBudget::default_limit));}result["likelihood_plan_keys"]=NumericVector::create(_["cyclic_plans"]=double(owner.likelihood_cyclic_plans),_["rejected_cyclic_plans"]=double(owner.likelihood_rejected_cyclic_plans));result["likelihood_work"]=NumericVector::create(_["model_build_inclusive_seconds"]=owner.likelihood_totals[0],_["topology_seconds"]=owner.likelihood_totals[1],_["forest_seconds"]=owner.likelihood_totals[2],_["message_initialization_seconds"]=owner.likelihood_totals[3],_["evaluation_inclusive_seconds"]=owner.likelihood_totals[4],_["operations"]=owner.likelihood_totals[5],_["lane_batches"]=owner.likelihood_totals[6],_["lane_fallbacks"]=owner.likelihood_totals[7],_["questions"]=owner.likelihood_totals[8],_["questions_with_ancestry_cut"]=owner.likelihood_totals[9]);if(owner.session&&!owner.session->store.config.likelihood_time_evaluations){NumericVector times=result["likelihood_work"];times["evaluation_inclusive_seconds"]=NA_REAL;result["likelihood_work"]=times;}
return result;
}
