#pragma once
#include "plan_evaluator.h"
namespace pap_internal::pedigree_likelihood {
template <int W> struct Lanes { typedef double Vector __attribute__((vector_size(W * sizeof(double)))); };

struct LaneInput {
  const LocusView* L = nullptr;
  const LocusMessages* M = nullptr;
  const Transmission* T = nullptr;
  int locus = 0;
};

// Interleaved transmission tables of one group of loci; fixed for a run.
struct LaneTables {
  std::vector<int> loci;
  std::vector<double> t, one, pairing;
};

LaneTables make_lane_tables(const LaneInput* in,int W,std::vector<int> key);
inline const LaneTables& lane_tables(std::map<std::vector<int>, LaneTables>& cache, const LaneInput* in, int W) {
  std::vector<int> key(W);for(int l=0;l<W;++l)key[l]=in[l].locus;
  auto at=cache.find(key);if(at!=cache.end())return at->second;
  return cache.emplace(key,make_lane_tables(in,W,key)).first->second;
}

// Scalar dispatch consumes a table synchronously before the next lookup. It
// may release all cached tables safely; LocusBatches holds long-lived table
// pointers and deliberately keeps its separate fixed grouping cache.
struct LaneTableCache {
  std::map<std::vector<int>,LaneTables> values;
  LaneTables uncached;
  bool bounded=false;std::size_t limit=8U*1024U*1024U;
  std::size_t bytes=0,peak_bytes=0,hits=0,misses=0,evictions=0,declines=0;
  static std::size_t charge(const std::vector<int>& key,const LaneTables& v);
  const LaneTables& get(const LaneInput* in,int W){
    std::vector<int> key(W);for(int l=0;l<W;++l)key[l]=in[l].locus;
    auto at=values.find(key);if(at!=values.end()){++hits;return at->second;}
    ++misses;auto built=make_lane_tables(in,W,key);auto n=charge(key,built);
    if(bounded&&n>limit){++declines;uncached=std::move(built);return uncached;}
    if(bounded&&bytes>limit-n){evictions+=values.size();values.clear();bytes=0;}
    auto inserted=values.emplace(std::move(key),std::move(built));bytes+=n;
    peak_bytes=std::max(peak_bytes,bytes);return inserted.first->second;
  }
};
inline const LaneTables& lane_tables(LaneTableCache& cache,const LaneInput* in,int W){return cache.get(in,W);}

// Returns false when any lane's value is zero or a lane cannot be run; the
// caller then evaluates the loci one at a time. log_z receives each lane's
// scale; for piece plans *result points at the interleaved retained table.
// On x86-64 Linux with GCC an AVX2 clone is chosen at run time when the CPU
// has it; elsewhere the portable build is used.
#if defined(__GNUC__) && !defined(__clang__) && defined(__x86_64__) && defined(__linux__) && !defined(__SANITIZE_THREAD__)
#define PAP_LANE_CLONES __attribute__((target_clones("avx2", "default")))
#else
#define PAP_LANE_CLONES
#endif
template <int W>
__attribute__((always_inline)) inline bool run_lanes_body(const Plan& plan, const Program& P, const Forest& F, const LaneInput* in,
                      const LaneTables& fixed_tables, Context& ctx, PairCache& cache, PieceCache* pieces, int loci,
                      std::vector<double>& arena, double* log_z, const double** result) {
  typedef typename Lanes<W>::Vector LaneVector;
  const int G = in[0].L->G;
  const double none = -std::numeric_limits<double>::infinity();
  thread_local std::vector<const Scaled*> piece_values_tls; auto& piece_values = piece_values_tls;
  if (!plan.pieces.empty()) {
    piece_values.assign(plan.pieces.size() * W, nullptr);
    for (std::size_t i = 0; i < plan.pieces.size(); ++i) {
      pieces->get_lanes(plan.pieces[i].slot, W, in, fixed_tables, loci, F, ctx, cache);
      for (int l = 0; l < W; ++l) {
        const Scaled* w = &pieces->get(plan.pieces[i].slot, in[l].locus, loci, F, *in[l].L, *in[l].M, *in[l].T, ctx, cache);
        if (!std::isfinite(w->scale)||w->wide()) return false;
        piece_values[i * W + l] = w;
      }
    }
  }
  // Interleaved copies of per-locus input tables follow the program arena.
  std::size_t extra = 0;
  for (const auto& slot : P.slots) {
    if (slot.kind == Program::PairTable) { const auto& pr = plan.pairs[slot.index]; extra += static_cast<std::size_t>(pr.vars[0] >= 0 ? G : 1) * (pr.vars[1] >= 0 ? G : 1); }
    else if (slot.kind == Program::PieceTable) { std::size_t c = 1; for (int v : plan.pieces[slot.index].vars) c *= plan.allele[v] ? in[0].L->k : G; extra += c; }
  }
  if (P.arena>P.limits.cells()/W || extra>P.limits.cells()/W-P.arena) return false;
  if (arena.size() < (P.arena + extra) * W) arena.resize((P.arena + extra) * W);
  double* A = arena.data();
  double* spare = A + P.arena * W;
  for (int l = 0; l < W; ++l) log_z[l] += P.log_uncovered;
  auto settle = [&](double* x, std::size_t cells, const double* top) {
    for (int l = 0; l < W; ++l) {
      if (top[l] > 1e-30 && top[l] < 1e30) continue;
      const double inverse = 1.0 / top[l];
      for (std::size_t c = 0; c < cells; ++c) x[c * W + l] *= inverse;
      log_z[l] += std::log(top[l]);
    }
  };
  const std::size_t npos = std::numeric_limits<std::size_t>::max();
  for (int v = 0; v < plan.variables; ++v) if (P.unary_offset[v] != npos) std::fill(A + P.unary_offset[v] * W, A + (P.unary_offset[v] + G) * W, 1.0);
  for (const auto& u : plan.unaries) {
    double* x = A + P.unary_offset[u.var] * W;
    for (int l = 0; l < W; ++l) {
      const LocusView& L = *in[l].L;
      switch (u.source) {
        case Plan::Emission: {
          const Scaled& e = in[l].M->emission_of(u.index, L);
          if (e.wide()||!std::isfinite(e.scale)) return false;
          for (int g = 0; g < G; ++g) x[g * W + l] *= e.cells[g];
          log_z[l] += e.scale; break;
        }
        case Plan::Prior: for (int g = 0; g < G; ++g) x[g * W + l] *= (*L.prior)[g]; break;
        case Plan::Incoming: case Plan::Outgoing: {
          const Scaled& s = u.source == Plan::Incoming ? ctx.to_individual(u.index, in[l].locus, L, *in[l].M) : ctx.to_family(u.index, in[l].locus, L, *in[l].M);
          if (s.wide()||!std::isfinite(s.scale)) return false;
          for (int g = 0; g < G; ++g) x[g * W + l] *= s.cells[g];
          log_z[l] += s.scale; break;
        }
      }
    }
  }
  thread_local std::vector<const double*> ptr_tls; auto& ptr = ptr_tls;
  ptr.assign(P.slots.size(), nullptr);
  for (std::size_t i = 0; i < P.slots.size(); ++i) {
    const auto& slot = P.slots[i];
    switch (slot.kind) {
      case Program::UnaryTable: {
        double* x = A + slot.offset * W; double top[W] = {};
        for (int g = 0; g < G; ++g) for (int l = 0; l < W; ++l) top[l] = std::max(top[l], x[g * W + l]);
        for (int l = 0; l < W; ++l) if (!(top[l] > 0.0)) return false;
        settle(x, G, top); ptr[i] = x; break;
      }
      case Program::GameteTable: ptr[i] = fixed_tables.t.data(); break;
      case Program::PairingTable: ptr[i] = fixed_tables.pairing.data(); break;
      case Program::OneTable: ptr[i] = fixed_tables.one.data(); break;
      case Program::OneGameteTable: break;
      case Program::PairTable: case Program::PieceTable: {
        double* x = spare;
        std::size_t cells = 0;
        for (int l = 0; l < W; ++l) {
          const Scaled& w = slot.kind == Program::PairTable ?
            cache.get(plan.pairs[slot.index].slot, in[l].locus, loci, F, *in[l].L, *in[l].M, ctx) : *piece_values[slot.index * W + l];
          if (!std::isfinite(w.scale)) return false;
          if(w.wide())return false;log_z[l] += w.scale; cells = w.cells.size();
          for (std::size_t c = 0; c < cells; ++c) x[c * W + l] = w.cells[c];
        }
        spare += cells * W; ptr[i] = x; break;
      }
      case Program::Produced: break;
    }
  }
  const int* lookup = in[0].L->lookup->data();
  const int* partner = in[0].T->partner.data();
  const std::size_t Gs = static_cast<std::size_t>(G), GK = Gs * static_cast<std::size_t>(in[0].L->k);
  thread_local std::vector<std::size_t> assignment_tls; auto& assignment = assignment_tls;
  thread_local std::vector<const double*> base_tls; auto& base = base_tls;
  thread_local std::vector<std::size_t> idx_tls, a0_tls; auto& idx = idx_tls; auto& a0 = a0_tls;
  for (const auto& step : P.steps) {
    const std::size_t n = step.inputs.size(), m = step.dims.size();
    double* out = A + step.offset * W;
    LaneVector top = {};
    base.resize(n); idx.assign(n, 0); a0.resize(n);
    for (std::size_t i = 0; i < n; ++i) { base[i] = ptr[step.inputs[i]]; a0[i] = step.axis0[i] * W; }
    assignment.assign(m, 0);
    const std::size_t* strides = step.strides.data();
    const int sel = step.select, role = step.select_role;
    auto read = [&](std::size_t i, int state, LaneVector& value) {
      const std::size_t at=idx[i]+state*step.axis0[i];
      const auto kind=P.slots[step.inputs[i]].kind;
      if(kind==Program::OneGameteTable){for(int l=0;l<W;++l){const int b=in[l].T->partner[at];value[l]=b<0?0.:(*in[l].L->p)[b];}}
      else if(in[0].T->implicit && (kind==Program::OneTable || kind==Program::PairingTable)) {
        for(int l=0;l<W;++l) value[l]=kind==Program::OneTable?in[l].T->one_value(*in[l].L,at):double(int(at%Gs)==(*in[l].L->lookup)[at/Gs]);
      } else __builtin_memcpy(&value,base[i]+at*W,sizeof(value));
    };
    bool structured=false;
    if(in[0].T->implicit || plan.factorized_one) for(std::size_t i=0;i<n;++i){const auto kind=P.slots[step.inputs[i]].kind;structured|=kind==Program::OneGameteTable||(in[0].T->implicit&&(kind==Program::OneTable||(kind==Program::PairingTable&&int(i)!=sel)));}
    for (std::size_t cell = 0; cell < step.cells; ++cell) {
      LaneVector sum = {};
      if (n == 0) {
        // A retained boundary can become independent after an evidence
        // factor is excluded. The empty product broadcasts one to all lanes.
        sum = LaneVector{} + double(step.dim);
      } else if (sel >= 0) {
        const std::size_t at = idx[sel];
        const int g = role == 0 ? lookup[at / Gs] : role == 1 ? partner[at % Gs + Gs * (at / GK)] : partner[at % Gs + Gs * (at / Gs)];
        if (g >= 0) {
          LaneVector value = LaneVector{} + 1.0;
          for (std::size_t i = 0; i < n; ++i) if (static_cast<int>(i) != sel) {
            LaneVector x; if(structured)read(i,g,x);else __builtin_memcpy(&x,base[i]+idx[i]*W+g*a0[i],sizeof(x)); value *= x;
          }
          sum = value;
        }
      } else if (structured) {
        for(int g=0;g<step.dim;++g){LaneVector value=LaneVector{}+1.;for(std::size_t i=0;i<n;++i){LaneVector x;read(i,g,x);value*=x;}sum+=value;}
      } else {
        const double* q[16];
        for (std::size_t i = 0; i < n; ++i) q[i] = base[i] + idx[i] * W;
        for (int g = 0; g < step.dim; ++g) {
          LaneVector value; __builtin_memcpy(&value, q[0], sizeof(value)); q[0] += a0[0];
          for (std::size_t i = 1; i < n; ++i) { LaneVector x; __builtin_memcpy(&x, q[i], sizeof(x)); value *= x; q[i] += a0[i]; }
          sum += value;
        }
      }
      __builtin_memcpy(out + cell * W, &sum, sizeof(sum));
      if(P.coalesced_inputs&&step.variable<0&&step.dim==1)for(int l=0;l<W;++l){
        if(!std::isfinite(sum[l])||sum[l]<0)return false;
        if(sum[l]>0)continue;bool positive=true;
        for(std::size_t i=0;i<n&&positive;++i){LaneVector value;read(i,0,value);positive=value[l]>0;}
        if(positive)return false; // scalar retry reports numerical unavailability
      }
      for (int l = 0; l < W; ++l) top[l] = std::max(top[l], sum[l]);
      for (std::size_t j = 0; j < m; ++j) {
        if (++assignment[j] < static_cast<std::size_t>(step.dims[j])) {
          for (std::size_t i = 0; i < n; ++i) idx[i] += strides[i * m + j];
          break;
        }
        assignment[j] = 0;
        for (std::size_t i = 0; i < n; ++i) idx[i] -= (step.dims[j] - 1) * strides[i * m + j];
      }
    }
    double tops[W];
    for (int l = 0; l < W; ++l) { tops[l] = top[l]; if (!(tops[l] > 0.0)) return false; }
    settle(out, step.cells, tops);
    ptr[step.out] = out;
  }
  for (int sc : P.scalars)
    for (int l = 0; l < W; ++l) {
      const double v = ptr[sc][l];
      if (!(v > 0.0)) return false;
      log_z[l] += std::log(v);
    }
  if (result) *result = P.result >= 0 ? ptr[P.result] : nullptr;
  return true;
}

PAP_LANE_CLONES inline bool run_lanes4(const Plan& plan, const Program& P, const Forest& F, const LaneInput* in, const LaneTables& fixed_tables,
                       Context& ctx, PairCache& cache, PieceCache* pieces, int loci, std::vector<double>& arena,
                       double* log_z, const double** result) {
  return run_lanes_body<kLanes>(plan, P, F, in, fixed_tables, ctx, cache, pieces, loci, arena, log_z, result);
}
#if defined(__GNUC__) && !defined(__clang__) && defined(__x86_64__)
#define PAP_WIDE_LANES 1
__attribute__((target("avx512f"))) inline bool run_lanes8(const Plan& plan, const Program& P, const Forest& F, const LaneInput* in,
                       const LaneTables& fixed_tables, Context& ctx, PairCache& cache, PieceCache* pieces, int loci,
                       std::vector<double>& arena, double* log_z, const double** result) {
  return run_lanes_body<kWideLanes>(plan, P, F, in, fixed_tables, ctx, cache, pieces, loci, arena, log_z, result);
}
#endif
// Eight lanes only where the CPU runs them natively.
inline bool wide_lanes_available() {
#ifdef PAP_WIDE_LANES
  static const bool available = __builtin_cpu_supports("avx512f");
  return available;
#else
  return false;
#endif
}
inline bool lane_memory_fits(const Plan& plan,const Program& P,int W) {
  if(!P.ready||P.exact)return false;
  std::size_t available=P.limits.cells()/W;if(P.arena>available)return false;available-=P.arena;
  for(const auto& slot:P.slots){std::size_t cells=0;
    if(slot.kind==Program::PairTable){const auto& r=plan.pairs[slot.index];cells=std::size_t(r.vars[0]>=0?P.G:1)*(r.vars[1]>=0?P.G:1);}
    else if(slot.kind==Program::PieceTable){cells=1;for(int v:plan.pieces[slot.index].vars){const std::size_t d=plan.allele[v]?P.k:P.G;if(cells>available/d)return false;cells*=d;}}
    if(cells>available)return false;available-=cells;
  }
  return true;
}
inline bool run_lanes(int W, const Plan& plan, const Program& P, const Forest& F, const LaneInput* in, const LaneTables& fixed_tables,
                      Context& ctx, PairCache& cache, PieceCache* pieces, int loci, std::vector<double>& arena,
                      double* log_z, const double** result) {
  if(P.exact)return false; // exact resource programs use the scalar route
#ifdef PAP_WIDE_LANES
  if (W == kWideLanes) return run_lanes8(plan, P, F, in, fixed_tables, ctx, cache, pieces, loci, arena, log_z, result);
#endif
  if (W != kLanes) return false;
  return run_lanes4(plan, P, F, in, fixed_tables, ctx, cache, pieces, loci, arena, log_z, result);
}

// Piece tables for a lane group: lanes not yet computed run as one batch; a
// lane batch that hits a zero leaves its lanes to the per-locus path.
template <class LaneInputT, class LaneTablesT>
inline void PieceCache::get_lanes(int id, int W, const LaneInputT* in, const LaneTablesT& fixed, int loci, const Forest& F,
                                  Context& ctx, PairCache& pairs) {
  Slot& s = *slots[id];
  if (s.ready.empty()) { s.values.resize(loci); s.ready.assign(loci, 0); }
  int missing = 0;
  for (int l = 0; l < W; ++l) missing += !s.ready[in[l].locus];
  if (missing < 2) return;
  const Program& P = program_for(s.plan, in[0].L->G, in[0].L->k);
  if (!P.ready) return;
  for (const auto& step : P.steps) if (step.inputs.size() > 16) return;
  double log_z[kWideLanes] = {};
  const double* table = nullptr;
  Workspace workspace(*this,s);
  if (!run_lanes(W, s.plan, P, F, in, fixed, ctx, pairs, nullptr, loci, workspace.arena, log_z, &table)) return;
  std::size_t cells = 1;
  for (int v = 0; v < s.plan.retained; ++v) cells *= s.plan.allele[v] ? in[0].L->k : in[0].L->G;
  for (int l = 0; l < W; ++l) {
    const int locus = in[l].locus;
    if (s.ready[locus]) continue;
    Scaled& out = s.values[locus];
    if (table) { out.cells.resize(cells); for (std::size_t c = 0; c < cells; ++c) out.cells[c] = table[c * W + l]; }
    else out.cells.assign(1, 1.0);
    out.scale = log_z[l];
    s.ready[locus] = 1; ++computed; operations += P.operations;
  }
}

// Value of a plan on W loci of one allele count; false sends the caller
// to evaluate each locus alone.
bool evaluate_lanes(int W, const Plan& plan, const Forest& F, const LaneInput* in, const LaneTables& fixed_tables, Context& ctx,
                           PairCache& cache, PieceCache& pieces, int loci, std::vector<double>& arena, double* value);

// The focal-table form for raw questions: each lane's table over the focal
// genotype, read into numerator (weighted by the focal emission) and
// denominator.
bool evaluate_focal_lanes(int W, const Plan& plan, const Forest& F, const LaneInput* in, const LaneTables& fixed_tables, Context& ctx,
                                 PairCache& cache, PieceCache& pieces, int loci, std::vector<double>& arena, int focal,
                                 double* numerator, double* denominator);

// A plan with the focal retained: its table over the focal's genotype (unit
// maximum) and the log scale, outside components included.
bool evaluate_focal_table(const Plan& plan, const Forest& F, const LocusView& L, const LocusMessages& M,
                                 const Transmission& T, Context& ctx, PairCache& cache, PieceCache& pieces, int locus, int loci,
                                 std::vector<double>& arena, std::vector<double>& table, double& log_z,Scaled* exact_table=nullptr);

double evaluate(const Plan& plan, const Forest& F, const LocusView& L, const LocusMessages& M,
                       const Transmission& T, Context& ctx, PairCache& cache, PieceCache& pieces, int locus, int loci,
                       std::vector<double>& arena);


}
