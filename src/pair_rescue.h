#ifndef PAP_PAIR_RESCUE_H
#define PAP_PAIR_RESCUE_H
#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <set>
#include <utility>
#include <vector>
namespace pap_internal::pair_rescue {
struct Score {
  double support=std::numeric_limits<double>::quiet_NaN();
  double relationship=std::numeric_limits<double>::quiet_NaN();
  double error=0;
  bool complete=false,cached=false,support_rejected=false,relationship_skipped=false,support_only=false;
};
struct Candidate {
  int id=0;Score score;
  bool feasible=false,incumbent=false,passing=false,retained=false;
};
// Unknown numerical evidence is never a negative screening certificate.
inline bool unknown(const Score& s) {
  if(s.support_only)return std::isnan(s.support);
  if(s.support_rejected)return false; // complete negative screen, unrequested relationship panel
  return !s.complete||std::isnan(s.support)||std::isnan(s.relationship)||!std::isfinite(s.error);
}
inline void select(std::vector<Candidate>& candidates,double threshold,double margin,int k) {
  std::vector<std::size_t> order;
  for(std::size_t i=0;i<candidates.size();++i) {
    auto& c=candidates[i];c.passing=c.retained=false;
    if(!c.feasible)continue;
    const bool pending=unknown(c.score);
    c.passing=!c.score.support_rejected&&(pending||(c.score.support>=threshold-margin&&
      (c.score.support_only||c.score.relationship+c.score.error>=threshold-margin)));
    // Incumbents remain in the exact retention comparison; unknown rows are
    // retained outside the cap rather than being rejected for numerical failure.
    c.retained=pending||c.incumbent;
    if(!pending)order.push_back(i);
  }
  // Score.support is the singleton screen value supplied by the caller.
  // The default identity screen uses I(single)-I(empty), with I(empty)=0.
  // A failed support gate can skip its relationship panel while retaining its
  // exact rank and potential role as a rescued mate.
  auto rank=[&](std::size_t i){return candidates[i].score.support;};
  std::sort(order.begin(),order.end(),[&](auto a,auto b){
    const double x=rank(a),y=rank(b);return x!=y?x>y:candidates[a].id<candidates[b].id;
  });
  for(int j=0;j<std::min(k,int(order.size()));++j)candidates[order[j]].retained=true;
}
inline bool pair(const Candidate& dam,const Candidate& sire,bool current) {
  return current||(dam.retained&&sire.retained&&(dam.passing||sire.passing));
}
// Restore only the rank heuristic: the independent support floor still holds.
inline bool rank_omitted_pair(const Candidate& dam,const Candidate& sire) {
  return dam.id!=sire.id&&dam.feasible&&sire.feasible&&
    (!dam.retained||!sire.retained)&&(dam.passing||sire.passing);
}
// A rank cap is a proposal heuristic, not an identity rejection certificate.
// Close a provisional winner under all one-parent replacements/completions.
inline std::vector<std::pair<int,int>> completions(int dam,int sire,
    std::span<const int> dams,std::span<const int> sires) {
  std::vector<std::pair<int,int>> out;
  if(dam){out.emplace_back(dam,0);for(int s:sires)if(dam!=s)out.emplace_back(dam,s);}
  if(sire){out.emplace_back(0,sire);for(int d:dams)if(d!=sire)out.emplace_back(d,sire);}
  return out;
}
struct Visit {
  int round=0,focal=0,revision=0,recheck=0;
  int singles=0,cache_hits=0,unknown=0,passing=0,retained=0,support_shortcuts=0;
  int pairs=0,kept=0,failed_both=0,rank_pruned=0,trio_pruned=0;
  int configurations_before=0,configurations_after=0,roots_before=-1,roots_after=-1;
  int completion_round=0,completion_rows=0;
  int rank_capped_dams=0,rank_capped_sires=0;
  int completion_identity_requests=0,completion_identity_hits=0,completion_identity_computations=0;
  int dam_floor_pass=0,sire_floor_pass=0;
  double seconds=0;
};
// Passive preparation accounting. Cache probes that miss are not counted twice when
// followed by joint_identity(). The scope also covers support counterfactuals.
inline thread_local const std::set<std::pair<int,int>>* completion_ids=nullptr;
inline thread_local Visit* completion_work=nullptr;
struct CompletionScope {
  const std::set<std::pair<int,int>>* old_ids=completion_ids;
  Visit* old_work=completion_work;
  explicit CompletionScope(const std::set<std::pair<int,int>>& ids) {
    completion_ids=&ids;completion_work=nullptr;
  }
  ~CompletionScope(){completion_ids=old_ids;completion_work=old_work;}
};
inline void identity_call(int dam,int sire,bool hit) {
  if(!completion_work||!completion_ids||!completion_ids->count({dam,sire}))return;
  ++completion_work->completion_identity_requests;
  if(hit)++completion_work->completion_identity_hits;
  else ++completion_work->completion_identity_computations;
}
struct Diagnostic {int round=0,focal=0,revision=0,recheck=0,side=0;Candidate candidate;};
struct PairDiagnostic {
  int round=0,focal=0,revision=0,recheck=0,dam=0,sire=0;
  bool incumbent=false,hard_feasible=false,trio_me_pass=false,candidate_keep=false;
  double trio_ln=std::numeric_limits<double>::quiet_NaN();
};
// Compact sufficient conditions for retaining an exact-grid identity
// competitor under every support-only margin/K choice. No truth is used.
inline int required_rank(const Candidate& c,const std::vector<Candidate>& side){
  if(c.incumbent||unknown(c.score))return 0;
  int rank=1;for(const auto& x:side)if(x.feasible&&!unknown(x.score)&&
    (x.score.support>c.score.support||(x.score.support==c.score.support&&x.id<c.id)))++rank;
  return rank;
}
inline double required_margin(const Candidate& d,const Candidate& s,double threshold){
  if(unknown(d.score)||unknown(s.score))return 0;
  return std::max(0.,std::min(threshold-d.score.support,threshold-s.score.support));
}
struct Requirement {
  int round=0,focal=0,revision=0,recheck=0,dam=0,sire=0,dam_rank=0,sire_rank=0;
  double margin=0;bool current=false,selected=false,trio_pass=false;
};
}
#endif
