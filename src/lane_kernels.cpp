#include "lane_kernels.h"
namespace pap_internal::pedigree_likelihood {
LaneTables make_lane_tables(const LaneInput* in,int W,std::vector<int> key){
  LaneTables x;x.loci=std::move(key);
  auto weave = [&](std::vector<double>& out, auto&& source) {
    const std::size_t size = source(0).size();
    out.assign(size * W, 0.0);
    for (int l = 0; l < W; ++l) {
      const auto& s = source(l);
      for (std::size_t c = 0; c < size; ++c) out[c * W + l] = s[c];
    }
  };
  weave(x.t, [&](int l) -> const std::vector<double>& { return *in[l].L->t; });
  weave(x.one, [&](int l) -> const std::vector<double>& { return in[l].T->one; });
  weave(x.pairing, [&](int l) -> const std::vector<double>& { return in[l].T->pairing; });
  return x;
}



bool evaluate_lanes(int W, const Plan& plan, const Forest& F, const LaneInput* in, const LaneTables& fixed_tables, Context& ctx,
                           PairCache& cache, PieceCache& pieces, int loci, std::vector<double>& arena, double* value) {
  const Program& P = program_for(plan, in[0].L->G, in[0].L->k);
  if (!P.ready) return false;
  for (const auto& step : P.steps) if (step.inputs.size() > 16) return false;
  double log_z[kWideLanes] = {};
  if (!run_lanes(W, plan, P, F, in, fixed_tables, ctx, cache, &pieces, loci, arena, log_z, nullptr)) return false;
  for (int l = 0; l < W; ++l) value[l] = log_z[l] + outside_log_z(plan, *in[l].L, *in[l].M, ctx, in[l].locus);
  return true;
}

bool evaluate_focal_lanes(int W, const Plan& plan, const Forest& F, const LaneInput* in, const LaneTables& fixed_tables, Context& ctx,
                                 PairCache& cache, PieceCache& pieces, int loci, std::vector<double>& arena, int focal,
                                 double* numerator, double* denominator) {
  const Program& P = program_for(plan, in[0].L->G, in[0].L->k);
  if (!P.ready || plan.retained != 1) return false;
  for (const auto& step : P.steps) if (step.inputs.size() > 16) return false;
  double log_z[kWideLanes] = {};
  const double* table = nullptr;
  if (!run_lanes(W, plan, P, F, in, fixed_tables, ctx, cache, &pieces, loci, arena, log_z, &table) || !table) return false;
  const int G = in[0].L->G;
  for (int l = 0; l < W; ++l) {
    const double z = log_z[l] + outside_log_z(plan, *in[l].L, *in[l].M, ctx, in[l].locus);
    const auto& row = (*in[l].L->phi)[focal];
    double top = -std::numeric_limits<double>::infinity();
    for (int g = 0; g < G; ++g) top = std::max(top, row[g]);
    double num = 0.0, den = 0.0;
    for (int g = 0; g < G; ++g) { const double x = table[g * W + l]; den += x; num += x * std::exp(row[g] - top); }
    denominator[l] = den > 0.0 ? std::log(den) + z : -std::numeric_limits<double>::infinity();
    numerator[l] = num > 0.0 && std::isfinite(top) ? std::log(num) + top + z : -std::numeric_limits<double>::infinity();
  }
  return true;
}

bool evaluate_focal_table(const Plan& plan, const Forest& F, const LocusView& L, const LocusMessages& M,
                                 const Transmission& T, Context& ctx, PairCache& cache, PieceCache& pieces, int locus, int loci,
                                 std::vector<double>& arena, std::vector<double>& table, double& log_z,Scaled* exact_table) {
  const Program& P = program_for(plan, L.G, L.k);
  if (!P.ready || plan.retained != 1) { log_z=std::numeric_limits<double>::quiet_NaN(); return false; }
  for (const auto& ref : plan.pieces) {
    const Scaled& w = pieces.get(ref.slot, locus, loci, F, L, M, T, ctx, cache);
    if (std::isnan(w.scale)) { log_z=std::numeric_limits<double>::quiet_NaN(); return false; }
  }
  log_z = 0.0;
  const double* result = nullptr;
  Scaled retained;
  if(exact_table){exact_table->cells.clear();exact_table->log_cells.clear();}
  if (!run_program(plan, P, F, L, M, T, ctx, cache, &pieces, locus, loci, arena, log_z, &result,&retained) || (!result&&!retained.wide())) {
    if(!std::isnan(log_z))log_z = -std::numeric_limits<double>::infinity(); return false;
  }
  if(retained.wide()){if(!exact_table){log_z=NAN;return false;}*exact_table=std::move(retained);table.clear();}
  else {table.assign(result,result+L.G);if(exact_table)exact_table->cells=table;}
  log_z += outside_log_z(plan, L, M, ctx, locus);
  return true;
}

double evaluate(const Plan& plan, const Forest& F, const LocusView& L, const LocusMessages& M,
                       const Transmission& T, Context& ctx, PairCache& cache, PieceCache& pieces, int locus, int loci,
                       std::vector<double>& arena) {
  const Program& P = program_for(plan, L.G, L.k);
  if (!P.ready) return std::numeric_limits<double>::quiet_NaN();
  for (const auto& ref : plan.pieces) {
    const Scaled& w = pieces.get(ref.slot, locus, loci, F, L, M, T, ctx, cache);
    if (std::isnan(w.scale)) return w.scale;
  }
  double log_z = 0.0;
  if (!run_program(plan, P, F, L, M, T, ctx, cache, &pieces, locus, loci, arena, log_z, nullptr))
    return std::isnan(log_z)?log_z:-std::numeric_limits<double>::infinity();
  const double value=log_z + outside_log_z(plan, L, M, ctx, locus);
  return value;
}
auto LaneTableCache::charge(const std::vector<int>& key,const LaneTables& v) -> std::size_t {
    return sizeof(LaneTables)+sizeof(std::vector<int>)+128+
      retained_vector(key)+retained_vector(v.loci)+retained_vector(v.t)+
      retained_vector(v.one)+retained_vector(v.pairing);
  }



}
