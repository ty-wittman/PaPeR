#include "plan_evaluator.h"
namespace pap_internal::pedigree_likelihood {
double outside_log_z(const Plan& plan, const LocusView& L, const LocusMessages& M, Context& ctx, int locus) {
  if(plan.ancestral){double z=0.;
    if(plan.barren_prior_normalizers)z+=plan.barren_prior_normalizers*std::log(std::accumulate(L.prior->begin(),L.prior->end(),0.));
    if(plan.barren_population_normalizers)z+=plan.barren_population_normalizers*std::log(std::accumulate(L.p->begin(),L.p->end(),0.));
    return z;}
  double z = 0.0;
  M.totals(L);
  if (std::isfinite(M.total_log_z)) {
    z = M.total_log_z;
    for (int c : plan.touched_components) z -= M.component_log_z[c];
  } else {
    std::size_t at = 0;
    for (int c = 0; c < static_cast<int>(M.component_log_z.size()); ++c) {
      if (at < plan.touched_components.size() && plan.touched_components[at] == c) { ++at; continue; }
      z += M.component_log_z[c];
    }
  }
  if (plan.family_side_untouched) z += ctx.family_side_log_z(locus, L, M);
  return z;
}

const Scaled& PieceCache::get(int id, int locus, int loci, const Forest& F, const LocusView& L, const LocusMessages& M,
                                     const Transmission& T, Context& ctx, PairCache& pairs) {
  ++requests;
  Slot& s = *slots[id];
  if (s.ready.empty()) { s.values.resize(loci); s.ready.assign(loci, 0); }
  Scaled& out = s.values[locus];
  if (s.ready[locus]) return out;
  const Program& P = program_for(s.plan, L.G, L.k);
  double log_z = 0.0;
  const double* table = nullptr;
  out.cells.clear();out.log_cells.clear();
  Workspace workspace(*this,s);
  if (!P.ready) out.scale = std::numeric_limits<double>::quiet_NaN();
  else if (!run_program(s.plan, P, F, L, M, T, ctx, pairs, nullptr, locus, loci, workspace.arena, log_z, &table, &out)) out.scale = std::isnan(log_z)?log_z:-std::numeric_limits<double>::infinity();
  else {
    std::size_t cells = 1;
    for (int v = 0; v < s.plan.retained; ++v) cells *= s.plan.allele[v] ? L.k : L.G;
    if (!out.wide() && table) {if(table!=out.cells.data())out.cells.assign(table, table + cells);}
    else if(!out.wide())out.cells.assign(1, 1.0);   // a piece touching no model individual is a scalar
    out.scale = log_z;
  }
  s.ready[locus] = 1; ++computed; operations += P.operations;
  return out;
}
}
