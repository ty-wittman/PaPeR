#include "pap_removed_path.h"
#ifndef PAP_RESOURCE_AUDIT_H
#define PAP_RESOURCE_AUDIT_H
#include <chrono>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
namespace pap_internal::pedigree_likelihood {
// Resource diagnostic only. Capture original factors, not compiled strides/order,
// for independently planned offline exact contractions. Disabled by default.
struct ResourceAudit {
  std::string directory;
  std::map<std::vector<int>,std::string> baseline_refusals;
  bool capture_all=false; // native reference-fixture tests only
  bool extended_precision=false; // unvalidated, W1 factor-audit probes only
  std::size_t extended_attempts=0,extended_successes=0,extended_peak_bytes=0;
  bool capture_unique=false;
  bool capture_three_genotype=false; // original factors for corrected admissions
  std::set<std::vector<int>> captured_shapes;
  std::size_t captured=0;
  double seconds=0;
};
inline ResourceAudit& resource_audit(){thread_local ResourceAudit x;return x;}
inline void configure_resource_audit(std::string directory){resource_audit()=ResourceAudit{};resource_audit().directory=std::move(directory);}
template<class PlanT> std::string baseline_resource_refusal(const PlanT& plan,int G,int k){ pap_internal::removed_path("baseline_resource_refusal"); }
template<class PlanT,class ProgramT,class LocusT,class TransmissionT,class PiecesT>
inline void capture_resource_program(const PlanT& plan,const ProgramT& p,const LocusT& l,const TransmissionT& t,
                                    PiecesT* pieces,const std::vector<const double*>& ptr,
                                    double input_scale,double output_scale,const double* result){
  auto& state=resource_audit();if(state.directory.empty()||(!state.capture_all&&p.k<=9&&!p.exact&&!p.coalesced_inputs))return;
  const auto start=std::chrono::steady_clock::now();
  const auto& generic=p;
  bool three=false;
  if(state.capture_three_genotype)for(const auto& step:generic.steps)
    three|=std::count(step.dims.begin(),step.dims.end(),p.G)>=3;
  if(state.capture_three_genotype&&!three)return;
  std::string reason=state.capture_three_genotype?"budget-admitted three-genotype program":state.capture_all?"native reference fixture":p.exact&&!plan.limits.resource_planning?"automatic exact resource retry":p.coalesced_inputs?"coalesced repeated-factor inputs":baseline_resource_refusal(plan,p.G,p.k);
  if(reason.empty()&&pieces)for(const auto& ref:plan.pieces){reason=baseline_resource_refusal(pieces->slots[ref.slot]->plan,p.G,p.k);if(!reason.empty()){reason="dependency: "+reason;break;}}
  if(reason.empty()){state.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return;}
  if((state.capture_unique||state.capture_three_genotype)&&!state.captured_shapes.insert(plan_shape(plan,p.G,p.k)).second)return;
  if(p.audit_scopes.empty()||p.audit_domains.empty())throw std::logic_error("resource audit lacks original factor scopes");
  const auto id=state.captured+1;
  const std::string filename="case-"+std::to_string(id)+".bin",path=state.directory+"/"+filename;
  std::ofstream f(path+".tmp",std::ios::binary);
  auto u64=[&](std::uint64_t x){f.write(reinterpret_cast<const char*>(&x),8);};
  auto real=[&](double x){f.write(reinterpret_cast<const char*>(&x),8);};
  auto text=[&](const std::string& x){u64(x.size());f.write(x.data(),x.size());};
  f.write("MP2AUD1\n",8);u64(p.k);u64(p.G);u64(plan.retained);u64(p.audit_domains.size());
  for(int d:p.audit_domains)u64(d);
  text(pieces?"main":"piece");text(reason);real(input_scale);real(output_scale);
  u64(p.audit_scopes.size());
  for(std::size_t i=0;i<p.audit_scopes.size();++i){
    const auto& scope=p.audit_scopes[i];u64(scope.size());std::size_t cells=1;
    for(int v:scope){u64(v);cells*=p.audit_domains[v];}u64(cells);
    const auto kind=p.slots[i].kind;
    auto generated=[&](auto read){double buffer[8192];for(std::size_t at=0;at<cells;at+=8192){const auto n=std::min(std::size_t(8192),cells-at);for(std::size_t j=0;j<n;++j)buffer[j]=read(at+j);f.write(reinterpret_cast<const char*>(buffer),n*sizeof(double));}};
    if(kind==ProgramT::PairingTable&&t.implicit){generated([&](std::size_t j){return double(int(j%l.G)==l.lookup->at(j/l.G));});}
    else if(kind==ProgramT::OneTable&&t.implicit){generated([&](std::size_t j){return t.one_value(l,j);});}
    else if(kind==ProgramT::OneGameteTable){generated([&](std::size_t j){int b=t.partner[j];return b<0?0.:l.p->at(b);});}
    else {if(!ptr[i])throw std::logic_error("missing resource-audit input");f.write(reinterpret_cast<const char*>(ptr[i]),cells*sizeof(double));}
  }
  std::size_t cells=1;for(int v=0;v<plan.retained;++v)cells*=p.audit_domains[v];u64(cells);
  if(result)f.write(reinterpret_cast<const char*>(result),cells*sizeof(double));else real(1.);
  f.close();if(!f||std::rename((path+".tmp").c_str(),path.c_str()))throw std::runtime_error("cannot save resource-audit case");
  ++state.captured;
  double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();state.seconds+=elapsed;
  std::ofstream index(state.directory+"/index.tsv",std::ios::app);
  index<<filename<<'\t'<<(pieces?"main":"piece")<<'\t'<<p.k<<'\t'<<reason<<'\t'<<elapsed<<'\n';
  if(!index)throw std::runtime_error("cannot save resource-audit index");
}
template<class ProgramT> inline void capture_resource_refusal(const ProgramT& p,const std::string& source){ pap_internal::removed_path("capture_resource_refusal"); }

}
#endif
