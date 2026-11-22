#pragma once
#include "plan_compiler.h"
namespace pap_internal::pedigree_likelihood {
double outside_log_z(const Plan& plan, const LocusView& L, const LocusMessages& M, Context& ctx, int locus);

// Range-safe execution of the SAME compiled instructions. No replanning,
// widened scope, work admission, or extended-precision retry is involved.
inline bool run_log_program(const Plan& plan,const Program& compiled,const Forest& F,const LocusView& L,
 const LocusMessages& M,const Transmission& T,Context& ctx,PairCache& cache,PieceCache* pieces,
 int locus,int loci,std::vector<double>& arena,double& log_z,const double** result,Scaled* retained){
  const Program& P=compiled;const int G=L.G;
  std::vector<double>().swap(arena);
  // Unaries use the admitted workspace prefix, just as in double execution.
  // Resource programs own only that prefix here; their intermediates are
  // allocated once in the log value type below.
  std::vector<double> exact_prefix(P.exact?P.arena-P.exact->arena:0);
  if(!P.exact)arena.resize(P.arena);
  double* storage=P.exact?exact_prefix.data():arena.data();
  std::vector<double*> unary(plan.variables,nullptr);
  for(int v=0;v<plan.variables;++v)if(P.unary_offset[v]!=std::numeric_limits<std::size_t>::max()){
    unary[v]=storage+P.unary_offset[v];std::fill_n(unary[v],G,0.);
  }
  log_z+=P.log_uncovered;
  for(const auto& u:plan.unaries){auto* row=unary[u.var];
    if(u.source==Plan::Prior){for(int g=0;g<G;++g)row[g]+=std::log((*L.prior)[g]);}
    else {const Scaled& x=u.source==Plan::Emission?M.emission_of(u.index,L):u.source==Plan::Incoming?ctx.to_individual(u.index,locus,L,M):ctx.to_family(u.index,locus,L,M);
      log_z+=x.scale;for(int g=0;g<G;++g)row[g]+=x.log_at(g);}
  }
  std::vector<const double*> inputs(P.slots.size(),nullptr),logs(P.slots.size(),nullptr);
  for(std::size_t i=0;i<P.slots.size();++i){const auto& slot=P.slots[i];
    if(slot.kind==Program::UnaryTable)logs[i]=unary[slot.index];
    else if(slot.kind==Program::GameteTable)inputs[i]=L.t->data();
    else if(slot.kind==Program::PairingTable)inputs[i]=T.implicit?nullptr:T.pairing.data();
    else if(slot.kind==Program::OneTable)inputs[i]=T.implicit?nullptr:T.one.data();
    else if(slot.kind==Program::PairTable||slot.kind==Program::PieceTable){
      const Scaled& x=slot.kind==Program::PairTable?cache.get(plan.pairs[slot.index].slot,locus,loci,F,L,M,ctx):pieces->get(plan.pieces[slot.index].slot,locus,loci,F,L,M,T,ctx,cache);
      log_z+=x.scale;if(x.wide())logs[i]=x.log_cells.data();else inputs[i]=x.cells.data();
    }
  }
  // Convert ordinary immutable sources once, rather than calling log at
  // every contraction multiply. These copies are charged to the existing
  // live workspace cap; skip a copy if it would exceed that cap.
  std::size_t free=P.arena<plan.limits.cells()?plan.limits.cells()-P.arena:0;
  std::vector<std::vector<double>> source_logs(P.audit_scopes.size());
  for(std::size_t i=0;i<P.audit_scopes.size();++i)if(!logs[i]){
    std::size_t count=1;for(int v:P.audit_scopes[i])count*=P.audit_domains[v];
    if(count>free)continue;
    const auto kind=P.slots[i].kind;auto& row=source_logs[i];row.resize(count);
    for(std::size_t j=0;j<count;++j){double value;
      if(kind==Program::PairingTable&&T.implicit)value=int(j%G)==(*L.lookup)[j/G];
      else if(kind==Program::OneTable&&T.implicit)value=T.one_value(L,j);
      else if(kind==Program::OneGameteTable){int b=T.partner[j];value=b<0?0.:(*L.p)[b];}
      else value=inputs[i][j];
      row[j]=std::log(value);
    }
    logs[i]=row.data();free-=count;
  }
  if(!std::isfinite(log_z))return false;
  if(P.exact){
    // LogPositive contains one double, so this has exactly the admitted live
    // workspace size. Release the old workspace before allocating its log view.
    std::vector<double>().swap(arena);std::vector<LogPositive> work(P.exact->arena);
    double value=0;bool ok=resource_plan::execute_precision<LogPositive>(*P.exact,L,T,inputs,work.data(),value,nullptr,&logs);
    log_z+=value;return ok;
  }
  // Generic instructions already encode indicator selections, axis strides
  // and retained outputs. Only multiplication/summation representation changes.
  auto read=[&](int id,std::size_t at){if(logs[id])return logs[id][at];
    const auto kind=P.slots[id].kind;
    if(kind==Program::Produced)return arena[P.slots[id].offset+at];
    if(kind==Program::PairingTable&&T.implicit)return int(at%G)==(*L.lookup)[at/G]?0.:-INFINITY;
    if(kind==Program::OneTable&&T.implicit)return std::log(T.one_value(L,at));
    if(kind==Program::OneGameteTable){int b=T.partner[at];return b<0?-INFINITY:std::log((*L.p)[b]);}
    return std::log(inputs[id][at]);};
  for(const auto& step:P.steps){
    const std::size_t n=step.inputs.size(),m=step.dims.size();
    std::vector<std::size_t> index(n,0),assignment(m,0);double top=-INFINITY;
    for(std::size_t c=0;c<step.cells;++c){double sum=-INFINITY;
      int first=0,last=step.dim;
      if(step.select>=0){const auto at=index[step.select];const int role=step.select_role;
        first=role==0?(*L.lookup)[at/G]:role==1?T.partner[at%G+G*(at/(G*L.k))]:T.partner[at%G+G*(at/G)];last=first<0?first:first+1;}
      for(int g=first;g<last;++g){double v=0.;for(std::size_t i=0;i<n;++i)if(int(i)!=step.select)v+=read(step.inputs[i],index[i]+g*step.axis0[i]);sum=log_add(sum,v);}
      arena[step.offset+c]=sum;top=std::max(top,sum);
      for(std::size_t j=0;j<m;++j){if(++assignment[j]<std::size_t(step.dims[j])){for(std::size_t i=0;i<n;++i)index[i]+=step.strides[i*m+j];break;}
        assignment[j]=0;for(std::size_t i=0;i<n;++i)index[i]-=(step.dims[j]-1)*step.strides[i*m+j];}
    }
    if(top==-INFINITY){log_z=-INFINITY;return false;}
    for(std::size_t c=0;c<step.cells;++c)arena[step.offset+c]-=top;
    log_z+=top;
  }
  for(int id:P.scalars)log_z+=read(id,0);
  if(P.result>=0){std::size_t count=1;for(int v=0;v<plan.retained;++v)count*=plan.allele[v]?L.k:G;
    std::vector<double> row(count);for(std::size_t i=0;i<count;++i)row[i]=read(P.result,i);
    Scaled local;store_logs(local,std::move(row),0.);log_z+=local.scale;local.scale=log_z;
    if(retained){*retained=std::move(local);if(result)*result=retained->wide()?nullptr:retained->cells.data();}
    else {if(local.wide()){log_z=NAN;return false;}std::copy(local.cells.begin(),local.cells.end(),arena.begin());if(result)*result=arena.data();}
  }
  return std::isfinite(log_z);
}

// Runs a compiled program; adds every scale to log_z. Returns false when the
// value is zero. For piece plans *result points at the retained table.
inline bool run_program(const Plan& plan, const Program& P, const Forest& F, const LocusView& L, const LocusMessages& M,
                        const Transmission& T, Context& ctx, PairCache& cache, PieceCache* pieces, int locus, int loci,
                        std::vector<double>& arena, double& log_z, const double** result,Scaled* retained_log=nullptr) {
  const double starting_log_z=log_z;
  const int G = L.G;
  const double none = -std::numeric_limits<double>::infinity();
  // Piece tables first: they run their own programs (which hold no pieces,
  // so they never touch this scratch).
  thread_local std::vector<const Scaled*> piece_values;
  if (!plan.pieces.empty()) {
    piece_values.assign(plan.pieces.size(), nullptr);
    for (std::size_t i = 0; i < plan.pieces.size(); ++i) {
      piece_values[i] = &pieces->get(plan.pieces[i].slot, locus, loci, F, L, M, T, ctx, cache);
      if (piece_values[i]->scale == none) return false;
    }
  }
  auto log_route=[&](){log_z=starting_log_z;return run_log_program(plan,P,F,L,M,T,ctx,cache,pieces,locus,loci,arena,log_z,result,retained_log);};
  if(!plan.pieces.empty())for(const auto* w:piece_values)if(w&&w->wide())return log_route();
  for(const auto& pair:plan.pairs)if(cache.get(pair.slot,locus,loci,F,L,M,ctx).wide())return log_route();
  for(const auto& u:plan.unaries)if(u.source!=Plan::Prior){const Scaled& x=u.source==Plan::Emission?M.emission_of(u.index,L):u.source==Plan::Incoming?ctx.to_individual(u.index,locus,L,M):ctx.to_family(u.index,locus,L,M);if(x.wide())return log_route();}
  // Tables are rescaled only when their maximum leaves [1e-30, 1e30]; the
  // scale moves into log_z. Products of up to ten tables stay in range.
  auto settle = [&](double* x, std::size_t cells, double top) {
    if (top > 1e-30 && top < 1e30) return;
    const double inverse = 1.0 / top;
    for (std::size_t c = 0; c < cells; ++c) x[c] *= inverse;
    log_z += std::log(top);
  };
  if (arena.size() < P.arena) arena.resize(P.arena);
  double* A = arena.data();
  log_z += P.log_uncovered;
  const std::size_t npos = std::numeric_limits<std::size_t>::max();
  for (int v = 0; v < plan.variables; ++v) if (P.unary_offset[v] != npos) std::fill(A + P.unary_offset[v], A + P.unary_offset[v] + G, 1.0);
  for (const auto& u : plan.unaries) {
    double* x = A + P.unary_offset[u.var];
    switch (u.source) {
      case Plan::Emission: {
        const Scaled& e = M.emission_of(u.index, L);
        if (!std::isfinite(e.scale)) return false;
        for (int g = 0; g < G; ++g) {double v=x[g]*e.cells[g];if(v<std::numeric_limits<double>::min()&&x[g]>0&&e.cells[g]>0)return log_route();x[g]=v;}
        log_z += e.scale; break;
      }
      case Plan::Prior: for (int g = 0; g < G; ++g) {double v=x[g]*(*L.prior)[g];if(v<std::numeric_limits<double>::min()&&x[g]>0&&(*L.prior)[g]>0)return log_route();x[g]=v;} break;
      case Plan::Incoming: case Plan::Outgoing: {
        const Scaled& s = u.source == Plan::Incoming ? ctx.to_individual(u.index, locus, L, M) : ctx.to_family(u.index, locus, L, M);
        if (s.scale == none) return false;
        for (int g = 0; g < G; ++g) {double v=x[g]*s.cells[g];if(v<std::numeric_limits<double>::min()&&x[g]>0&&s.cells[g]>0)return log_route();x[g]=v;}
        log_z += s.scale; break;
      }
    }
  }
  thread_local std::vector<const double*> ptr;
  ptr.assign(P.slots.size(), nullptr);
  for (std::size_t i = 0; i < P.slots.size(); ++i) {
    const auto& slot = P.slots[i];
    switch (slot.kind) {
      case Program::UnaryTable: {
        double* x = A + slot.offset; double top = 0.0;
        for (int g = 0; g < G; ++g) top = std::max(top, x[g]);
        if (!(top > 0.0)) return plan.limits.contraction_log_retry?log_route():false;
        settle(x, G, top); ptr[i] = x; break;
      }
      case Program::GameteTable: ptr[i] = L.t->data(); break;
      case Program::PairingTable: ptr[i] = T.implicit ? nullptr : T.pairing.data(); break;
      case Program::OneTable: ptr[i] = T.implicit ? nullptr : T.one.data(); break;
      case Program::OneGameteTable: break;
      case Program::PairTable: {
        const Scaled& w = cache.get(plan.pairs[slot.index].slot, locus, loci, F, L, M, ctx);
        if (w.scale == none) return false;
        log_z += w.scale; ptr[i] = w.cells.data(); break;
      }
      case Program::PieceTable: {
        const Scaled& w = *piece_values[slot.index];
        log_z += w.scale; ptr[i] = w.cells.data(); break;
      }
      case Program::Produced: break;
    }
  }
  const double audit_input_scale=log_z;
  // Packed offsets may overwrite original unary factors. Preserve the audit
  // inputs before the first contraction; execution keeps its original pointers.
  std::vector<double> audit_unaries;
  std::vector<const double*> audit_ptr;
  if(P.limits.live_arena && !resource_audit().directory.empty()) {
    audit_ptr=ptr;
    std::size_t count=0;for(std::size_t i=0;i<P.audit_scopes.size();++i)count+=P.slots[i].kind==Program::UnaryTable;
    audit_unaries.resize(count*G);std::size_t at=0;
    for(std::size_t i=0;i<P.audit_scopes.size();++i)if(P.slots[i].kind==Program::UnaryTable){
      std::copy_n(ptr[i],G,audit_unaries.data()+at);audit_ptr[i]=audit_unaries.data()+at;at+=G;
    }
  }
  thread_local std::vector<std::size_t> index, assignment;
  if(P.exact){
    const auto start=std::chrono::steady_clock::now();resource_plan::RunStats stats;double value=0;
    const std::size_t prefix=P.arena-P.exact->arena;
    bool ok=resource_plan::execute(*P.exact,L,T,ptr,A+prefix,value,resource_audit().directory.empty()?nullptr:&stats);
    // Retry the same certified instructions at extended precision when the
    // double arithmetic cannot preserve positive entries. Release the reusable
    // double arena before allocating the wider intermediates; immutable input
    // factors remain external, with only its unary prefix copied. Charge both
    // actual live allocation phases to the unchanged numerical storage limit.
    std::vector<double> retry_unaries;
    if(resource_audit().extended_precision&&!ok&&std::isnan(value)&&std::numeric_limits<long double>::max_exponent>std::numeric_limits<double>::max_exponent){
      ++resource_audit().extended_attempts;
      const auto cells=P.limits.cells();
      const std::size_t wide_ratio=(sizeof(long double)+sizeof(double)-1)/sizeof(double);
      if(prefix<=cells&&arena.capacity()<=cells-prefix&&P.exact->arena<=(cells-prefix)/wide_ratio){
        resource_audit().extended_peak_bytes=std::max(resource_audit().extended_peak_bytes,
          std::max((arena.capacity()+prefix)*sizeof(double),prefix*sizeof(double)+P.exact->arena*sizeof(long double)));
        retry_unaries.assign(A,A+prefix);
        for(std::size_t j=0;j<P.audit_scopes.size();++j)if(P.slots[j].kind==Program::UnaryTable)
          ptr[j]=retry_unaries.data()+P.slots[j].offset;
        std::vector<double>().swap(arena);
        std::vector<long double> wide(P.exact->arena);
        ok=resource_plan::execute_precision<long double>(*P.exact,L,T,ptr,wide.data(),value,resource_audit().directory.empty()?nullptr:&stats);
        resource_audit().extended_successes+=ok;
      }
    }
    log_z+=value;
    if(!resource_audit().directory.empty()){
      capture_resource_program(plan,P,L,T,pieces,audit_ptr.empty()?ptr:audit_ptr,audit_input_scale,log_z,nullptr);
      std::ofstream receipt(resource_audit().directory+"/execution.tsv",std::ios::app);
      receipt<<P.profile_identity<<'\t'<<P.k<<'\t'<<P.exact->condition<<'\t'<<std::setprecision(17)<<P.operations<<'\t'<<stats.multiplies<<'\t'<<stats.states<<'\t'<<P.arena*8<<'\t'<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\t'<<ok<<'\t'<<stats.matrix_calls<<'\n';
    }
    return ok;
  }
  for (const auto& step : P.steps) {
    const std::size_t n = step.inputs.size(), m = step.dims.size();
    double* out = A + step.offset;
    double top = 0.0;
    // Fixed input counts let the compiler unroll the product; the arithmetic
    // and its order are the same as the general loop below.
    const int* lookup = L.lookup->data();
    const int* partner = T.partner.data();
    const std::size_t Gs = static_cast<std::size_t>(G), GK = Gs * static_cast<std::size_t>(L.k);
    auto fixed = [&](auto count) {
      constexpr int N = decltype(count)::value;
      const double* base[N]; std::size_t a0[N], idx[N];
      for (int i = 0; i < N; ++i) { base[i] = ptr[step.inputs[i]]; a0[i] = step.axis0[i]; idx[i] = 0; }
      assignment.assign(m, 0);
      const std::size_t* strides = step.strides.data();
      const int sel = step.select, role = step.select_role;
      for (std::size_t cell = 0; cell < step.cells; ++cell) {
        const double* q[N];
        for (int i = 0; i < N; ++i) q[i] = base[i] ? base[i] + idx[i] : nullptr;
        double sum = 0.0;
        if (sel >= 0) {
          // The indicator is one at a single state of the eliminated variable.
          const std::size_t at = idx[sel];
          const int g = role == 0 ? lookup[at / Gs] : role == 1 ? partner[at % Gs + Gs * (at / GK)] : partner[at % Gs + Gs * (at / Gs)];
          if (g >= 0) {
            double value = 1.0;
            for (int i = 0; i < N; ++i) if (i != sel) value *= q[i][g * a0[i]];
            sum = value;
          }
        } else {
          // Two interleaved partial sums shorten the dependency chain.
          double s0 = 0.0, s1 = 0.0;
          int g = 0;
          for (; g + 1 < step.dim; g += 2) {
            double v0 = q[0][0], v1 = q[0][a0[0]];
            for (int i = 1; i < N; ++i) { v0 *= q[i][0]; v1 *= q[i][a0[i]]; }
            s0 += v0; s1 += v1;
            for (int i = 0; i < N; ++i) q[i] += 2 * a0[i];
          }
          if (g < step.dim) {
            double v = q[0][0];
            for (int i = 1; i < N; ++i) v *= q[i][0];
            s0 += v;
          }
          sum = s0 + s1;
        }
        out[cell] = sum; top = std::max(top, sum);
        for (std::size_t j = 0; j < m; ++j) {
          if (++assignment[j] < static_cast<std::size_t>(step.dims[j])) {
            for (int i = 0; i < N; ++i) idx[i] += strides[i * m + j];
            break;
          }
          assignment[j] = 0;
          for (int i = 0; i < N; ++i) idx[i] -= (step.dims[j] - 1) * strides[i * m + j];
        }
      }
    };
    bool structured = false;
    if(T.implicit || plan.factorized_one) for (std::size_t i=0;i<n;++i) {
      const auto kind=P.slots[step.inputs[i]].kind;
      structured |= kind==Program::OneGameteTable || (T.implicit && (kind==Program::OneTable || (kind==Program::PairingTable && int(i)!=step.select)));
    }
    if (structured) {
      index.assign(n,0); assignment.assign(m,0);
      auto read = [&](std::size_t i,int state) {
        const std::size_t at=index[i]+state*step.axis0[i];
        const auto kind=P.slots[step.inputs[i]].kind;
        if(kind==Program::OneGameteTable){const int b=partner[at];return b<0?0.:(*L.p)[b];}
        if(kind==Program::OneTable) return T.one_value(L,at);
        if(kind==Program::PairingTable) return double(int(at%Gs)==lookup[at/Gs]);
        return ptr[step.inputs[i]][at];
      };
      for(std::size_t cell=0;cell<step.cells;++cell) {
        double sum=0;
        if(step.select>=0) {
          const auto at=index[step.select];const int role=step.select_role;
          const int g=role==0?lookup[at/Gs]:role==1?partner[at%Gs+Gs*(at/GK)]:partner[at%Gs+Gs*(at/Gs)];
          if(g>=0){sum=1;for(std::size_t i=0;i<n;++i)if(int(i)!=step.select)sum*=read(i,g);}
        } else {
          double s0=0,s1=0;int g=0;
          for(;g+1<step.dim;g+=2){double v0=1,v1=1;for(std::size_t i=0;i<n;++i){v0*=read(i,g);v1*=read(i,g+1);}s0+=v0;s1+=v1;}
          if(g<step.dim){double v=1;for(std::size_t i=0;i<n;++i)v*=read(i,g);s0+=v;}
          sum=s0+s1;
        }
        out[cell]=sum;top=std::max(top,sum);
        for(std::size_t j=0;j<m;++j){if(++assignment[j]<std::size_t(step.dims[j])){for(std::size_t i=0;i<n;++i)index[i]+=step.strides[i*m+j];break;}
          assignment[j]=0;for(std::size_t i=0;i<n;++i)index[i]-=(step.dims[j]-1)*step.strides[i*m+j];}
      }
    } else switch (n) {
      case 1: fixed(std::integral_constant<int, 1>{}); break;
      case 2: fixed(std::integral_constant<int, 2>{}); break;
      case 3: fixed(std::integral_constant<int, 3>{}); break;
      case 4: fixed(std::integral_constant<int, 4>{}); break;
      case 5: fixed(std::integral_constant<int, 5>{}); break;
      case 6: fixed(std::integral_constant<int, 6>{}); break;
      case 7: fixed(std::integral_constant<int, 7>{}); break;
      case 8: fixed(std::integral_constant<int, 8>{}); break;
      case 9: fixed(std::integral_constant<int, 9>{}); break;
      case 10: fixed(std::integral_constant<int, 10>{}); break;
      default: {
        index.assign(n, 0); assignment.assign(m, 0);
        std::vector<const double*> inp(n); std::vector<std::size_t> ax(n);
        for (std::size_t i = 0; i < n; ++i) { inp[i] = ptr[step.inputs[i]]; ax[i] = step.axis0[i]; }
        for (std::size_t cell = 0; cell < step.cells; ++cell) {
          double sum = 0.0;
          for (int g = 0; g < step.dim; ++g) {
            double value = 1.0;
            for (std::size_t i = 0; i < n; ++i) value *= inp[i][index[i] + g * ax[i]];
            sum += value;
          }
          out[cell] = sum; top = std::max(top, sum);
          for (std::size_t j = 0; j < m; ++j) {
            if (++assignment[j] < static_cast<std::size_t>(step.dims[j])) {
              for (std::size_t i = 0; i < n; ++i) index[i] += step.strides[i * m + j];
              break;
            }
            assignment[j] = 0;
            for (std::size_t i = 0; i < n; ++i) index[i] -= (step.dims[j] - 1) * step.strides[i * m + j];
          }
        }
      }
    }
    if(P.coalesced_inputs&&step.variable<0&&step.dim==1){
      // Newly introduced products must not turn positive operands into a
      // structural zero. Refuse numerical range loss; the scalar caller
      // distinguishes NaN (unavailable) from an exact zero return.
      for(std::size_t cell=0;cell<step.cells;++cell){
        if(!std::isfinite(out[cell])||out[cell]<0){log_z=std::numeric_limits<double>::quiet_NaN();return false;}
        if(out[cell]>0)continue;bool positive=true;
        for(std::size_t i=0;i<n&&positive;++i){std::size_t rest=cell,at=0;
          for(std::size_t j=0;j<m;++j){at+=(rest%step.dims[j])*step.strides[i*m+j];rest/=step.dims[j];}
          const int slot=step.inputs[i];const auto kind=P.slots[slot].kind;double value;
          if(kind==Program::OneGameteTable){const int b=T.partner[at];value=b<0?0.:(*L.p)[b];}
          else if(kind==Program::OneTable&&T.implicit)value=T.one_value(L,at);
          else if(kind==Program::PairingTable&&T.implicit)value=double(int(at%G)==(*L.lookup)[at/G]);
          else value=ptr[slot][at];
          positive=value>0;
        }
        if(positive){if(plan.limits.contraction_log_retry)return log_route();log_z=std::numeric_limits<double>::quiet_NaN();return false;}
      }
    }
    if (!(top > 0.0)) return plan.limits.contraction_log_retry?log_route():false;
    settle(out, step.cells, top);
    ptr[step.out] = out;
  }
  for (int sc : P.scalars) {
    const double v = ptr[sc][0];
    if (!(v > 0.0)) return plan.limits.contraction_log_retry?log_route():false;
    log_z += std::log(v);
  }
  if (result) *result = P.result >= 0 ? ptr[P.result] : nullptr;
  if(!resource_audit().directory.empty())capture_resource_program(plan,P,L,T,pieces,audit_ptr.empty()?ptr:audit_ptr,audit_input_scale,log_z,P.result>=0?ptr[P.result]:nullptr);
  return true;
}



// ---------------------------------------------------------------------------
// Several loci with one allele count run one compiled program together. Every
// table holds the loci as interleaved lanes (cell * W + lane) and each lane
// keeps its own scale, so the product and sum vectorize across loci. The
// arithmetic per lane is the same as run_program's.
constexpr int kLanes = 4;      // portable / AVX2 width
constexpr int kWideLanes = 8;  // AVX-512 width

}
