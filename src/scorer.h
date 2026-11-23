#ifndef PAP_SCORER_H
#define PAP_SCORER_H
// Live scorer protocol shared by the public Stage-1 engine and pedigree likelihood scorer.
#include "score_types.h"
#include "topology_cache.h"
#include "motif_templates.h"
#include <vector>
#include "pair_rescue.h"
namespace pap_internal { struct ExecutionConfig; }
namespace pap_internal::native_demand {
struct Tags {int kind=0,phase=0,focal=0,candidate=0,side=0,dam=0,sire=0,term=-1,refinement=0;};
struct AuthorityRow {Tags tags;int code=0,status=0;bool supported=false;std::string name,relationship,reason;};
struct GateRow {Tags tags;int code=0,locus=0,k=0,status=0,state=0;bool supported=false;double value=0,error=0;std::string old_reason,reason;};
struct LiveScorer {
 virtual std::shared_ptr<const ExecutionConfig> execution_config()const{return {};}
 bool prepared_score_batches=false;
 virtual void prepare_scores(const std::vector<std::pair<int,Tags>>&){};
 // Absolute whole-panel parent envelope on one immutable focal evidence union.
 virtual double joint_identity(int,int){return std::numeric_limits<double>::quiet_NaN();}
 virtual bool cached_joint_identity(int,int,double&){return false;}
 virtual double joint_support(int dam,int sire){return joint_identity(dam,sire);}
 virtual bool prefer_joint_support_first()const{return false;}
 virtual pair_rescue::Score single_candidate(int,int,double,bool,bool){return {};}
 virtual std::array<int,2> freeze_identity_roots(const std::vector<int>&){return {{-1,-1}};}

 bool native_support=false,native_catalog=false,authority_audit=false,exact_step_scaling=false,exterior_work_planning=false,exterior_state_products=false,shared_exact_cavities=false,shared_exact_reference=false;
 virtual std::shared_ptr<const side_likelihood::Boundary> catalog(Tags,bool){return {};}
 virtual std::shared_ptr<const predictive_likelihood::RawBoundary> raw_catalog(int,const std::vector<std::array<int,2>>&,const std::vector<unsigned char>&){return {};}
 virtual void observe_raw_catalog(int,const predictive_likelihood::RawBoundary&){}
 virtual void observe_catalog(Tags,const side_likelihood::Boundary&){}
 virtual void observe_gate(Tags,int,int,int,const side_likelihood::Model&,const side_likelihood::PartitionValue&){}
 bool skip_c_admissions=false;bool engine_result_reset=false;std::size_t invalidated_result_rows=0;bool engine_move_results=false,engine_reduction_reuse=false,engine_compact_results=false;std::size_t compact_packets_avoided=0;void construction_budget(const Rcpp::List& flags){if(!flags.containsElementNamed("pedigree_reuse_byte_limit"))return;double n=Rcpp::as<double>(flags["pedigree_reuse_byte_limit"]);if(!std::isfinite(n)||n<0||n>64U*1024U*1024U||n!=std::floor(n))throw std::invalid_argument("M_ped reuse budget must be an integer from zero through 64 MiB");engine_life_cache.limit=static_cast<size_t>(n);motif_cache.limit=std::min(size_t(n),size_t(1024U*1024U));}bool motif_template_reuse=false;motif_templates::Cache motif_cache;bool engine_topology_reuse=false;topology_life_reuse::Cache engine_life_cache;bool engine_descriptor = false;  virtual ~LiveScorer()=default;virtual void begin_visit(int,const std::vector<int>&,const std::vector<int>&,int,bool,const Rcpp::List&){};virtual bool accepted_move(int,int,int,int,int,const std::vector<int>&,const std::vector<int>&,int){return false;};virtual side_likelihood::PartitionValue score(int hypothesis,Tags tags,int locus,int k)=0;};
inline thread_local LiveScorer* live=nullptr;
inline bool engine_descriptor_live(){return live && live->engine_descriptor;}
struct LiveScope {LiveScorer* previous;topology_life_reuse::Scope topology;motif_templates::Scope motifs;explicit LiveScope(LiveScorer* next):previous(live),topology(next&&next->engine_topology_reuse?&next->engine_life_cache:nullptr),motifs(next&&next->motif_template_reuse?&next->motif_cache:nullptr){live=next;}~LiveScope(){live=previous;}};

}
#endif
