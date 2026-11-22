#ifndef PAP_LOCUS_BATCHES_H
#define PAP_LOCUS_BATCHES_H
#include "native_locus_executor.h"
#include "pedigree_likelihood.h"
#include <chrono>
#include <atomic>

namespace pap_internal::pedigree_likelihood {
// Run-owned threads and numeric scratch. Catalogs and compiled programs are
// supplied by the caller; no R values or biological model builders enter jobs.
class LocusBatches {
public:
  enum class Scheduler { Static, Dynamic, Sticky };
  static Scheduler parse_scheduler(const std::string& name) {
    if(name=="static")return Scheduler::Static;
    if(name=="dynamic")return Scheduler::Dynamic;
    if(name=="sticky")return Scheduler::Sticky;
    throw std::invalid_argument("mped2_scheduler must be static, dynamic or sticky");
  }
  struct Request {const Plan* plan;bool raw;std::vector<double>* values;std::vector<double>* denominators;};
  struct RefusalDetail {std::shared_ptr<const Program> program;std::size_t count=0;};
  using RefusalDetails=std::map<std::pair<std::string,std::vector<int>>,RefusalDetail>;
private:
  using Clock = std::chrono::steady_clock;
  struct Batch { std::vector<int> loci; int width; const LaneTables* tables = nullptr;std::array<const LaneTables*,2> narrow{}; };
  struct Worker {
    std::vector<Batch> batches;
    Context context;
    PairCache pairs;
    PieceCache pieces;
    std::vector<double> arena, lane_arena, focal_table;
    std::size_t messages = 0, lanes = 0, fallbacks = 0, narrowed = 0;
    // Passive run-lifetime counts; read/merge only after workers have joined.
    std::map<std::pair<int, std::string>, std::size_t> refusals;
    RefusalDetails refusal_details;
    int preferred=0,last_worker=-1;
    std::size_t claims=0,migrations=0,active_lanes=0,lane_slots=0;
  };
  NativeLocusWorkerExecutor<NativeNoopWorkerGuard> executor;
  // Static owners retain the original per-worker grouping. Dynamic/sticky
  // owners hold one complete SIMD bundle, including all numeric caches.
  std::vector<Worker> owners;
  std::vector<double> worker_seconds;
  Scheduler scheduler;
  std::atomic<std::size_t> cursor{0};
  std::unique_ptr<std::atomic<bool>[]> claimed;
  std::array<std::size_t,5> past_cache{};
  std::size_t pair_slot_copies=0,piece_slot_copies=0;
  std::vector<LocusView> views;
  std::vector<const Transmission*> transmission;
  std::vector<std::unique_ptr<LocusMessages>>* messages = nullptr;
  std::map<std::vector<int>, LaneTables> fixed_tables;
  std::vector<std::pair<int,int>> domains;
  const Forest* forest = nullptr;
  std::vector<Request> requests;
  int focal = 0;
  double minimum_operations=0;
  bool message_job = false;
  static double seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
  struct Scope {
    std::size_t* old_count; bool old_prepared;
    explicit Scope(std::size_t& count):old_count(worker_message_count),old_prepared(prepared_programs_only) {
      worker_message_count=&count; prepared_programs_only=true;
    }
    ~Scope(){worker_message_count=old_count;prepared_programs_only=old_prepared;}
  };
  static void execute(void* payload,int id) { static_cast<LocusBatches*>(payload)->run_worker(id); }
  void record_refusal(Worker& w, const Plan& plan, int locus) {
    const auto& L=views[locus];
    const Program& main=program_for(plan,L.G,L.k);
    const Plan* refused_plan=nullptr;
    std::string source;
    std::string reason;
    if(!main.ready) {source="main";refused_plan=&plan;reason=std::string("main: ")+(main.refusal?main.refusal:"unspecified");}
    else if(main.exact&&!transmission[locus]->two_support){source="transmission";reason="exact resource plan requires at most two supported gametes per genotype";}
    else for(const auto& ref:plan.pieces) {
      const auto& slot=*w.pieces.slots[ref.slot];
      if(!slot.ready[locus] || !std::isnan(slot.values[locus].scale))continue;
      const Program& piece=program_for(slot.plan,L.G,L.k);
      source="piece";refused_plan=&slot.plan;
      reason=std::string("piece: ")+(piece.refusal?piece.refusal:"nonfinite value");
      break; // Same first unavailable dependency as the scalar evaluator.
    }
    if(reason.empty())reason="nonfinite result without program refusal";
    ++w.refusals[{L.k,reason}];
    if(refused_plan)capture_resource_refusal(*refused_plan->programs[L.k],source);
    if(resource_diagnostics&&refused_plan){
      auto& d=w.refusal_details[{source,plan_shape(*refused_plan,L.G,L.k)}];
      d.program=refused_plan->programs[L.k];++d.count;
    }
  }
  void dispatch(bool threaded=true) {
    cursor.store(0,std::memory_order_relaxed);
    if(scheduler==Scheduler::Sticky)for(std::size_t j=0;j<owners.size();++j)claimed[j].store(false,std::memory_order_relaxed);
    if(executor.effective_workers==1||!threaded){for(int id=0;id<executor.effective_workers;++id)run_worker(id);}
    else { executor.dispatch({this,&execute,"M_ped2 batch"});executor.wait(); }
    for(auto& w:owners){forest->messages_computed+=w.messages;w.messages=0;
      if(scheduler==Scheduler::Sticky&&w.last_worker>=0)w.preferred=w.last_worker;
    }
  }
  void run_worker(int id) {
    struct Measure {double& value;bool enabled;Clock::time_point start;
      Measure(double& v,bool on):value(v),enabled(on){if(on)start=Clock::now();}
      ~Measure(){if(enabled)value=std::chrono::duration<double>(Clock::now()-start).count();}
    } measure(worker_seconds[id],profile);
    auto run=[&](std::size_t j){auto& w=owners[j];
      if(!message_job){++w.claims;if(w.last_worker>=0&&w.last_worker!=id)++w.migrations;}
      w.last_worker=id;run_owner(w);
    };
    if(scheduler==Scheduler::Static){run(id);return;}
    // Each claim covers a complete SIMD bundle across all authorized requests.
    // Claims are exclusive; executor publication/join orders cache ownership
    // between waves. No thread-local numerical state is retained by a bundle.
    if(scheduler==Scheduler::Sticky)for(std::size_t j=0;j<owners.size();++j)
      if(owners[j].preferred==id&&!claimed[j].exchange(true,std::memory_order_relaxed))run(j);
    for(;;){const auto j=cursor.fetch_add(1,std::memory_order_relaxed);if(j>=owners.size())break;
      if(scheduler!=Scheduler::Sticky||!claimed[j].exchange(true,std::memory_order_relaxed))run(j);
    }
  }
  void run_owner(Worker& w) {
    Scope scope(w.messages);
    if(message_job){for(const auto& batch:w.batches)for(int l:batch.loci){auto& m=*messages->at(l);if(!m.F)compute_messages(*forest,views[l],m);else m.totals(views[l]);}return;}
    for(const auto& request:requests){
      const auto* plan=request.plan;auto* values=request.values;auto* denominators=request.denominators;const bool focal_form=request.raw;
    for(const auto& batch:w.batches) {
      const int count=int(batch.loci.size()), W=batch.width;
      bool ok=false; double a[kWideLanes],b[kWideLanes];
      if(count>1&&W>1) {
        LaneInput in[kWideLanes];
        for(int i=0;i<W;++i){int l=batch.loci[i<count?i:0];in[i]={&views[l],messages->at(l).get(),transmission[l],l};}
        const Program* compiled=adaptive_lanes&&W==kWideLanes?&program_for(*plan,in[0].L->G,in[0].L->k):nullptr;
        if(compiled&&!lane_memory_fits(*plan,*compiled,W)&&lane_memory_fits(*plan,*compiled,kLanes)){
          ok=true;++w.narrowed;
          for(int half=0;half<2&&ok;++half){LaneInput small[kLanes];std::vector<int> ids;
            for(int j=0;j<kLanes;++j){int i=half*kLanes+j;small[j]=in[i<count?i:half*kLanes];if(i<count)ids.push_back(small[j].locus);}
            ok=focal_form?evaluate_focal_lanes(kLanes,*plan,*forest,small,*batch.narrow[half],w.context,w.pairs,w.pieces,int(views.size()),w.lane_arena,focal,a+half*kLanes,b+half*kLanes):
              evaluate_lanes(kLanes,*plan,*forest,small,*batch.narrow[half],w.context,w.pairs,w.pieces,int(views.size()),w.lane_arena,a+half*kLanes);
          }
        }else ok=focal_form?evaluate_focal_lanes(W,*plan,*forest,in,*batch.tables,w.context,w.pairs,w.pieces,int(views.size()),w.lane_arena,focal,a,b):
          evaluate_lanes(W,*plan,*forest,in,*batch.tables,w.context,w.pairs,w.pieces,int(views.size()),w.lane_arena,a);
        if(ok){++w.lanes;w.active_lanes+=count;w.lane_slots+=W;for(int i=0;i<count;++i){(*values)[batch.loci[i]]=a[i];if(focal_form)(*denominators)[batch.loci[i]]=b[i];}}
        else ++w.fallbacks;
      }
      if(ok)continue;
      for(int l:batch.loci) {
        const auto& L=views[l];const auto& M=*messages->at(l);const auto& T=*transmission[l];
        if(!focal_form){(*values)[l]=evaluate(*plan,*forest,L,M,T,w.context,w.pairs,w.pieces,l,int(views.size()),w.arena);continue;}
        double z=0;
        if(!evaluate_focal_table(*plan,*forest,L,M,T,w.context,w.pairs,w.pieces,l,int(views.size()),w.arena,w.focal_table,z)) {
          (*values)[l]=(*denominators)[l]=z;continue;
        }
        const auto& row=(*L.phi)[focal];double top=-INFINITY,num=0,den=0;
        for(int g=0;g<L.G;++g)top=std::max(top,row[g]);
        for(int g=0;g<L.G;++g){den+=w.focal_table[g];num+=w.focal_table[g]*std::exp(row[g]-top);}
        (*denominators)[l]=den>0?std::log(den)+z:-INFINITY;
        (*values)[l]=num>0&&std::isfinite(top)?std::log(num)+top+z:-INFINITY;
      }
    }
    // Inspect completed values only: no new programs, cache fills or questions.
    for(const auto& batch:w.batches)for(int l:batch.loci)
      if(std::isnan((*values)[l]))record_refusal(w,*plan,l);
    }
  }
public:
  // Diagnostic lookup after the executor barrier. Inspect only already-built
  // programs and already-evaluated piece dependencies; never request evidence.
  struct Failure {std::shared_ptr<const Program> program;std::string source,reason;int piece=-1;};
  Failure failure(const Plan& plan,int locus)const {
    Failure out;const int k=views.at(locus).k;
    if(k<int(plan.programs.size())&&plan.programs[k]&&!plan.programs[k]->ready){out.program=plan.programs[k];out.source="main";out.reason=out.program->refusal?out.program->refusal:"unspecified";return out;}
    if(k<int(plan.programs.size())&&plan.programs[k]&&plan.programs[k]->exact&&!transmission[locus]->two_support){out.source="transmission";out.reason="exact resource plan requires at most two supported gametes per genotype";return out;}
    for(const auto& w:owners)for(const auto& batch:w.batches)
      if(std::find(batch.loci.begin(),batch.loci.end(),locus)!=batch.loci.end()){
        for(const auto& ref:plan.pieces){const auto& slot=*w.pieces.slots.at(ref.slot);
          if(slot.ready.empty()||!slot.ready[locus]||!std::isnan(slot.values[locus].scale))continue;
          out.piece=ref.slot;out.source="piece";
          if(k<int(slot.plan.programs.size()))out.program=slot.plan.programs[k];
          out.reason=out.program&&out.program->refusal?out.program->refusal:"nonfinite piece value";return out;
        }
      }
    out.source="numerical";out.reason="nonfinite result without program refusal";return out;
  }
  std::size_t jobs=0,parallel_jobs=0,serial_jobs=0,prepared_plans=0,prepared_pieces=0;
  double preparation_seconds=0,execution_seconds=0,message_seconds=0;
  bool profile=false;
  bool resource_diagnostics=false,adaptive_lanes=false;
  double worker_busy_seconds=0,worker_critical_seconds=0,dispatch_gap_seconds=0,imbalance_seconds=0;
  explicit LocusBatches(int count,int loci,double threshold=0,bool profiling=false,Scheduler scheduling=Scheduler::Static):executor(count,loci),worker_seconds(executor.effective_workers),scheduler(scheduling),minimum_operations(threshold),profile(profiling) {
    if(scheduler==Scheduler::Static)owners.resize(executor.effective_workers);
  }
  const char* scheduler_name()const{return scheduler==Scheduler::Static?"static":scheduler==Scheduler::Dynamic?"dynamic":"sticky";}
  int worker_count()const{return executor.effective_workers;}
  // Called once per run. Fixed input views and interleaved transmission tables
  // remain alive across visits; each locus has exactly one numeric cache owner.
  void initialize(std::vector<LocusView> input,std::vector<const Transmission*> tables,bool lanes,bool wide) {
    if(!views.empty())return;
    views=std::move(input);transmission=std::move(tables);
    std::map<std::pair<int,int>,std::vector<int>> groups;
    for(int l=0;l<int(views.size());++l)groups[{views[l].G,views[l].k}].push_back(l);
    std::vector<Batch> batches;
    const int cap=lanes?(wide&&wide_lanes_available()?kWideLanes:kLanes):1;
    for(const auto& [domain,ls]:groups){domains.push_back(domain);for(std::size_t at=0;at<ls.size();at+=cap){Batch b;const auto end=std::min(ls.size(),at+cap);b.loci.assign(ls.begin()+at,ls.begin()+end);b.width=b.loci.size()>kLanes?kWideLanes:(b.loci.size()>1?kLanes:1);batches.push_back(std::move(b));}}
    auto weight=[&](const Batch& b){double G=views[b.loci[0]].G;return G*G*G*b.loci.size();};
    std::stable_sort(batches.begin(),batches.end(),[&](const Batch& a,const Batch& b){return weight(a)>weight(b);});
    std::vector<double> load(executor.effective_workers);
    if(scheduler!=Scheduler::Static){owners.resize(batches.size());claimed=std::make_unique<std::atomic<bool>[]>(batches.size());}
    std::size_t owner=0;
    for(auto& b:batches){
      if(b.width>1){LaneInput in[kWideLanes];for(int j=0;j<b.width;++j){int l=b.loci[j<int(b.loci.size())?j:0];in[j]={&views[l],nullptr,transmission[l],l};}b.tables=&lane_tables(fixed_tables,in,b.width);}
      if(adaptive_lanes&&b.width==kWideLanes)for(int half=0;half<2;++half){LaneInput in[kLanes];for(int j=0;j<kLanes;++j){int i=half*kLanes+j;int l=b.loci[i<int(b.loci.size())?i:half*kLanes];in[j]={&views[l],nullptr,transmission[l],l};}b.narrow[half]=&lane_tables(fixed_tables,in,kLanes);}
      auto id=std::min_element(load.begin(),load.end())-load.begin();load[id]+=weight(b);
      auto& w=owners[scheduler==Scheduler::Static?id:owner++];w.preferred=int(id);w.batches.push_back(std::move(b));
    }
  }
  void begin_visit(const Forest& F,int focal_id,std::vector<std::unique_ptr<LocusMessages>>& slots,bool parallel_messages,bool latent_brood_exclusion=false,bool ancestral=false,bool ancestral_reuse=false) {
    forest=&F;focal=focal_id;messages=&slots;
    for(auto& w:owners){
      past_cache[0]+=w.pairs.requests;past_cache[1]+=w.pairs.computed;
      past_cache[2]+=w.pieces.requests;past_cache[3]+=w.pieces.computed;past_cache[4]+=w.context.repaired;
      w.context.reset(F,focal_id,int(views.size()),latent_brood_exclusion,ancestral_reuse);w.pairs.clear();w.pieces.clear();
    }
    // Allocate owners on caller. Message contents are independently initialized
    // by locus; stale values retain the original Forest revision semantics.
    for(auto& m:slots){if(!m)m=std::make_unique<LocusMessages>();
      if(m->F&&m->latent_brood_exclusion!=latent_brood_exclusion)throw std::logic_error("M_ped2 brood policy changed within a message store");
      m->latent_brood_exclusion=latent_brood_exclusion;
      if(m->F&&m->ancestral!=ancestral)throw std::logic_error("M_ped2 evidence policy changed within message store");
      m->ancestral=ancestral;m->ancestral_reuse=ancestral_reuse;}
    const auto t=Clock::now();message_job=true;
    if(parallel_messages)dispatch();
    else for(std::size_t l=0;l<views.size();++l)if(!slots[l]->F)compute_messages(F,views[l],*slots[l]);
    message_job=false;message_seconds+=seconds(t);
  }
  void evaluate_all(const Plan& P,PairCache& source_pairs,PieceCache& source_pieces,bool raw,
                    std::vector<double>& out,std::vector<double>& den) {
    evaluate_many({Request{&P,raw,&out,&den}},source_pairs,source_pieces);
  }
  void evaluate_many(std::vector<Request> next,PairCache& source_pairs,PieceCache& source_pieces){
    if(next.empty())return;
    auto t=Clock::now();requests=std::move(next);
    // Compile once on caller for every domain, including newly interned pieces.
    // Existing programs are shared const objects; worker slots contain numeric
    // values and arenas only. No per-job topology/catalog reconstruction.
    const std::size_t old=owners.front().pieces.slots.size();
    for(const auto& r:requests)for(auto [G,k]:domains)(void)program_for(*r.plan,G,k);
    for(std::size_t j=old;j<source_pieces.slots.size();++j)
      for(auto [G,k]:domains)(void)program_for(source_pieces.slots[j]->plan,G,k);
    for(auto& w:owners){
      while(w.pairs.slots.size()<source_pairs.slots.size()){
        const auto& src=source_pairs.slots[w.pairs.slots.size()];PairCache::Slot slot;
        slot.family=src.family;slot.dam=src.dam;slot.sire=src.sire;slot.edges=src.edges;
        slot.values.resize(views.size());slot.ready.assign(views.size(),0);w.pairs.slots.push_back(std::move(slot));++pair_slot_copies;
      }
      while(w.pieces.slots.size()<source_pieces.slots.size()){
        auto slot=std::make_unique<PieceCache::Slot>();slot->plan=source_pieces.slots[w.pieces.slots.size()]->plan;
        slot->values.resize(views.size());slot->ready.assign(views.size(),0);w.pieces.slots.push_back(std::move(slot));++piece_slot_copies;
      }
    }
    prepared_plans+=requests.size();prepared_pieces+=source_pieces.slots.size()-old;
    preparation_seconds+=seconds(t);
    double operations=0;for(const auto& r:requests)for(const auto& L:views)operations+=program_for(*r.plan,L.G,L.k).operations;
    const bool threaded=executor.effective_workers>1&&operations>=minimum_operations;
    if(threaded)++parallel_jobs;else ++serial_jobs;
    t=Clock::now();dispatch(threaded);const double wall=seconds(t);execution_seconds+=wall;++jobs;
    if(profile){double sum=0,critical=0;for(double t:worker_seconds){sum+=t;critical=std::max(critical,t);}
      // Inline shards execute sequentially; only concurrent jobs have imbalance.
      if(!threaded)critical=sum;
      worker_busy_seconds+=sum;worker_critical_seconds+=critical;
      dispatch_gap_seconds+=std::max(0.,wall-critical);
      if(threaded)imbalance_seconds+=executor.effective_workers*critical-sum;
    }
  }
  std::size_t narrowed_batches()const{std::size_t n=0;for(const auto& w:owners)n+=w.narrowed;return n;}
  std::size_t pair_requests()const{std::size_t n=0;for(const auto& w:owners)n+=w.pairs.requests;return n;}
  std::size_t pair_computed()const{std::size_t n=0;for(const auto& w:owners)n+=w.pairs.computed;return n;}
  std::size_t repairs()const{std::size_t n=0;for(const auto& w:owners)n+=w.context.repaired;return n;}
  // Run totals; the original per-visit accessors above retain their contract.
  std::array<std::size_t,5> cache_counts()const {
    auto n=past_cache;for(const auto& w:owners){n[0]+=w.pairs.requests;n[1]+=w.pairs.computed;
      n[2]+=w.pieces.requests;n[3]+=w.pieces.computed;n[4]+=w.context.repaired;}return n;
  }
  std::array<std::size_t,9> scheduling_counts()const {
    std::array<std::size_t,9> n{};n[0]=owners.size();n[7]=pair_slot_copies;n[8]=piece_slot_copies;
    for(const auto& w:owners){n[1]+=w.claims;n[2]+=w.migrations;n[3]+=w.lanes;n[4]+=w.fallbacks;n[5]+=w.active_lanes;n[6]+=w.lane_slots;}return n;
  }
  std::map<std::pair<int,std::string>,std::size_t> refusal_counts()const {
    std::map<std::pair<int,std::string>,std::size_t> out;
    for(const auto& w:owners)for(const auto& [key,n]:w.refusals)out[key]+=n;
    return out;
  }
  RefusalDetails resource_refusals()const {
    RefusalDetails out;
    for(const auto& w:owners)for(const auto& [key,d]:w.refusal_details){auto& x=out[key];x.program=d.program;x.count+=d.count;}
    return out;
  }
};
}
#endif
