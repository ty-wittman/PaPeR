#pragma once
#include "pedigree_forest.h"
namespace pap_internal::pedigree_likelihood {
struct Scaled {
  // Exactly one payload is live. Log cells preserve positive entries whose
  // normalized range cannot be represented by double; -Inf is a true zero.
  std::vector<double> cells, log_cells;
  double scale = 0.0;
  bool wide() const { return !log_cells.empty(); }
  std::size_t size() const { return wide()?log_cells.size():cells.size(); }
  double log_at(std::size_t i) const { return wide()?log_cells[i]:std::log(cells[i]); }
  template<class Real> Real at(std::size_t i) const {
    if constexpr(std::is_same_v<Real,LogPositive>)return LogPositive::from_log(log_at(i));
    else { if(wide())throw StorageRange();return cells[i]; }
  }
};
void store_logs(Scaled& out, std::vector<double> logs, double scale);
inline void normalize(Scaled& x);
// Source-state missingness, never a forest cut: two assigned parents retain
// their brood even when one of their forest ports was cut to break a loop.
inline bool latent_brood_upward(const Forest& F, int e) {
  const auto& edge=F.edges[e];
  if(edge.role>=2)return false;
  const auto& f=F.families[edge.family];
  return f.dam==0 || f.sire==0;
}

inline void normalize(Scaled& x) {
  if(x.wide()){store_logs(x,x.log_cells,x.scale);return;}
  double top = 0.0;
  for (double v : x.cells) top = std::max(top, v);
  if (!(top > 0.0)) { x.scale = -std::numeric_limits<double>::infinity(); return; }
  if (!std::isfinite(top)) throw std::runtime_error("M_ped2 message is not finite");
  for(double v:x.cells)if(v>0 && v/top<std::numeric_limits<double>::min()){std::vector<double> logs;for(double w:x.cells)logs.push_back(std::log(w));store_logs(x,std::move(logs),x.scale);return;}
  for (double& v : x.cells) v /= top;
  x.scale += std::log(top);
}

// Stored messages of one locus, computed on demand and kept until the
// forest marks them stale. Each message records the revision it was computed
// at; it is current while that is no older than the edge's stale mark.
struct LocusMessages {
  const Forest* F = nullptr;
  bool latent_brood_exclusion = false;
  bool ancestral = false;
  bool ancestral_reuse = false;
  mutable std::vector<Scaled> family_cells;      // individual -> family, per edge
  mutable std::vector<Scaled> individual_cells;  // family -> individual, per edge
  mutable std::vector<int> family_stamp, individual_stamp;
  mutable std::vector<double> component_log_z;
  mutable double total_log_z = 0.0;   // sum over components
  mutable int totals_revision = -1;
  mutable std::vector<Scaled> emissions;   // per individual, exp(phi - top) and top
  // In ancestral mode the outgoing value does not depend on the mate family,
  // and a family's downward value does not depend on the receiving sibling.
  // Keep exactly one of each, retaining only current families after an update.
  mutable std::vector<Scaled> ancestor_cells, birth_cells;
  mutable std::vector<int> ancestor_stamp, birth_stamp;   // birth: per stable family id
  mutable Scaled unit;
  mutable int shared_revision = -1;
  bool shared() const { return ancestral && ancestral_reuse; }
  void attach(const Forest& forest);
  void sync_shared() const;
  bool suppressed(int e,bool side) const { return side ? F->edges[e].role==2 : F->edges[e].role<2; }
  Scaled& cells_for(int e,bool side) const {
    if(!shared())return (side?family_cells:individual_cells)[e];
    if(suppressed(e,side))return unit;
    return side?ancestor_cells[F->edges[e].individual]:birth_cells[F->family_stable[F->edges[e].family]];
  }
  int& stamp_for(int e,bool side) const {
    return shared() ? (side?ancestor_stamp[F->edges[e].individual]:birth_stamp[F->family_stable[F->edges[e].family]]) :
      (side?family_stamp:individual_stamp)[e];
  }
  const Scaled& emission_of(int id, const LocusView& L) const;
  bool current(int e, bool to_family_side) const;
  const Scaled& to_family(int e, const LocusView& L) const { ensure(e, true, L); return cells_for(e,true); }
  const Scaled& to_individual(int e, const LocusView& L) const { ensure(e, false, L); return cells_for(e,false); }
  void ensure(int e, bool to_family_side, const LocusView& L) const;
  void totals(const LocusView& L) const;
};

inline void emission(const LocusView& L, int v, Scaled& out) {
  const auto& row = (*L.phi)[v];
  out.log_cells.clear();out.cells.resize(L.G);
  double top = -std::numeric_limits<double>::infinity();
  for (int g = 0; g < L.G; ++g) top = std::max(top, row[g]);
  if(top==-INFINITY){out.cells.assign(L.G,0.);out.scale=-INFINITY;return;}
  if (!std::isfinite(top)) throw std::runtime_error("M_ped2 invalid emission");
  for (int g = 0; g < L.G; ++g) {out.cells[g] = std::exp(row[g] - top);if(std::isfinite(row[g])&&out.cells[g]<std::numeric_limits<double>::min()){store_logs(out,row,0.);return;}}
  out.scale = top;
}
inline void multiply_into(Scaled& x, const Scaled& y) {
  bool wide=x.wide()||y.wide();
  if(!wide)for(std::size_t i=0;i<x.cells.size();++i)
    if(x.cells[i]>0&&y.cells[i]>0&&x.cells[i]*y.cells[i]<std::numeric_limits<double>::min()){wide=true;break;}
  if(wide){std::vector<double> logs(x.size());for(std::size_t i=0;i<logs.size();++i)logs[i]=x.log_at(i)+y.log_at(i);store_logs(x,std::move(logs),x.scale+y.scale);}
  else {for(std::size_t i=0;i<x.cells.size();++i)x.cells[i]*=y.cells[i];x.scale+=y.scale;}
}
inline void multiply_prior(Scaled& out,const LocusView& L){Scaled prior;prior.cells=*L.prior;multiply_into(out,prior);}
inline double log_sum(const Scaled& x){double sum=-INFINITY;for(std::size_t i=0;i<x.size();++i)sum=log_add(sum,x.log_at(i));return sum;}
// Individual v along edge e: emission, prior when founder/cut, other incoming.
// `incoming(e)` returns the family->individual message or null to skip a link.
template <class Incoming>
inline void individual_message_with(const Forest& F, const LocusView& L, int v, int exclude, bool founder,
                                    Incoming&& incoming, Scaled& out) {
  emission(L, v, out);
  if (founder) multiply_prior(out,L);
  for (int e : F.individual_edges[v]) if (e != exclude) if (const Scaled* m = incoming(e)) multiply_into(out, *m);
  normalize(out);
}
inline void individual_message(const Forest& F, const LocusView& L, const LocusMessages& M,
                               int v, int exclude, Scaled& out) {
  individual_message_with(F, L, v, exclude, F.child_edge[v] < 0, [&](int e) { return &M.to_individual(e, L); }, out);
}

// Per-child likelihood over (dam genotype, sire genotype) given the child's
// message; an absent side uses population gametes. Table is G_d x G_s with
// absent sides of extent 1.
template<class Real>
inline void child_table_values(const LocusView& L, const Scaled& child, bool dam, bool sire,
                        std::vector<Real>& out) {
  const int G = L.G, k = L.k;
  if constexpr(std::is_same_v<Real,double>){if(child.wide())throw StorageRange();
    double low=1.,weight=1.;for(double x:child.cells)if(x>0)low=std::min(low,x);
    const auto& weights=*L.p;for(double x:weights)if(x>0)weight=std::min(weight,x);
    for(double x:*L.t)if(x>0)weight=std::min(weight,x);
    if((long double)low*weight*weight<std::numeric_limits<double>::min())throw StorageRange();}
  const auto& t = *L.t; const auto& lookup = *L.lookup; const auto& p = *L.p;
  const int nd = dam ? G : 1, ns = sire ? G : 1;
  if (L.gamete_allele) {
    // A parent genotype transmits at most two alleles: sum only those.
    const int* al = L.gamete_allele->data(); const double* wt = L.gamete_weight->data();
    std::vector<Real> Q(static_cast<std::size_t>(k) * ns, 0.0);
    for (int y = 0; y < ns; ++y)
      for (int a = 0; a < k; ++a) {
        Real acc = 0.0;
        if (sire) {
          for (int r = 0; r < 2; ++r) { const int b = al[2 * y + r]; if (b >= 0) acc += wt[2 * y + r] * child.at<Real>(lookup[a + k * b]); }
        } else for (int b = 0; b < k; ++b) acc += p[b] * child.at<Real>(lookup[a + k * b]);
        Q[a + static_cast<std::size_t>(k) * y] = acc;
      }
    out.assign(static_cast<std::size_t>(nd) * ns, 0.0);
    for (int y = 0; y < ns; ++y)
      for (int x = 0; x < nd; ++x) {
        Real acc = 0.0;
        if (dam) {
          for (int r = 0; r < 2; ++r) { const int a = al[2 * x + r]; if (a >= 0) acc += wt[2 * x + r] * Q[a + static_cast<std::size_t>(k) * y]; }
        } else for (int a = 0; a < k; ++a) acc += p[a] * Q[a + static_cast<std::size_t>(k) * y];
        out[x + static_cast<std::size_t>(nd) * y] = acc;
      }
    return;
  }
  // Q(a, s) = sum_b gamete_s(b) child(lookup(a,b))
  std::vector<Real> Q(static_cast<std::size_t>(k) * ns, 0.0);
  for (int y = 0; y < ns; ++y)
    for (int a = 0; a < k; ++a) {
      Real acc = 0.0;
      for (int b = 0; b < k; ++b) acc += (sire ? t[y + G * b] : p[b]) * child.at<Real>(lookup[a + k * b]);
      Q[a + static_cast<std::size_t>(k) * y] = acc;
    }
  out.assign(static_cast<std::size_t>(nd) * ns, 0.0);
  for (int y = 0; y < ns; ++y)
    for (int x = 0; x < nd; ++x) {
      Real acc = 0.0;
      for (int a = 0; a < k; ++a) acc += (dam ? t[x + G * a] : p[a]) * Q[a + static_cast<std::size_t>(k) * y];
      out[x + static_cast<std::size_t>(nd) * y] = acc;
    }
}

inline void child_table_scaled(const LocusView& L,const Scaled& child,bool dam,bool sire,Scaled& out){
  out.log_cells.clear();out.scale=0.;
  try{child_table_values<double>(L,child,dam,sire,out.cells);}
  catch(const StorageRange&){std::vector<LogPositive> values;child_table_values<LogPositive>(L,child,dam,sire,values);std::vector<double> logs;for(auto x:values)logs.push_back(x.log);store_logs(out,std::move(logs),0.);}
}
// Message from family f along edge `target`.
// `outgoing(e)` returns the individual->family message or null to skip a link.
template <class Real,class Outgoing>
inline void family_message_values(const Forest& F, const LocusView& L, int f, int target, Outgoing&& outgoing, Scaled& out) {
  const int G = L.G, k = L.k;
  const int de = F.family_parent_edge[f][0], se = F.family_parent_edge[f][1];
  const bool dam = de >= 0, sire = se >= 0;
  const int nd = dam ? G : 1, ns = sire ? G : 1;
  std::vector<Real> W(static_cast<std::size_t>(nd) * ns, 1.0), table;
  double scale = 0.0;
  for (int e : F.family_child_edges[f]) {
    if (e == target) continue;
    const Scaled* in = outgoing(e);
    if (!in) continue;
    const Scaled& c = *in;
    child_table_values<Real>(L, c, dam, sire, table);
    Real top = 0.0;
    for (std::size_t i = 0; i < W.size(); ++i) { W[i] = storage_product(W[i],table[i]); top = std::max(top, W[i]); }
    if (!(top > 0.0)) { out.log_cells.clear();out.cells.assign(G, 0.0); out.scale = -std::numeric_limits<double>::infinity(); return; }
    for (Real& w : W) {Real scaled=w/top;if constexpr(std::is_same_v<Real,double>){if(w>0. && scaled<std::numeric_limits<double>::min())throw StorageRange();}w=scaled;}
    scale += c.scale + positive_log(top);
  }
  const Scaled* md = dam && de != target ? outgoing(de) : nullptr;
  const Scaled* ms = sire && se != target ? outgoing(se) : nullptr;
  auto weight_dam = [&](int x) { return md ? md->at<Real>(x) : 1.0; };
  auto weight_sire = [&](int y) { return ms ? ms->at<Real>(y) : 1.0; };
  if constexpr(std::is_same_v<Real,double>){
    long double lo=1.;for(double w:W)if(w>0)lo=std::min(lo,(long double)w);
    for(const Scaled* m:{md,ms})if(m){if(m->wide())throw StorageRange();double low=1.;for(double x:m->cells)if(x>0)low=std::min(low,x);lo*=low;}
    if(lo>0 && lo<std::numeric_limits<double>::min())throw StorageRange();
  }
  if (md) scale += md->scale;
  if (ms) scale += ms->scale;
  std::vector<Real> result(G,0.);
  if (target == de) {
    for (int x = 0; x < G; ++x) { Real acc = 0.0; for (int y = 0; y < ns; ++y) acc += W[x + static_cast<std::size_t>(nd) * y] * weight_sire(y); result[x] = acc; }
  } else if (target == se) {
    for (int y = 0; y < G; ++y) { Real acc = 0.0; for (int x = 0; x < nd; ++x) acc += W[x + static_cast<std::size_t>(nd) * y] * weight_dam(x); result[y] = acc; }
  } else if (L.gamete_allele) {  // target is a child; parents transmit at most two alleles
    const auto& lookup = *L.lookup; const auto& p = *L.p;
    const int* al = L.gamete_allele->data(); const double* wt = L.gamete_weight->data();
    std::vector<Real> X(static_cast<std::size_t>(k) * ns, 0.0);  // X(a, y)
    for (int y = 0; y < ns; ++y)
      for (int x = 0; x < nd; ++x) {
        const Real w = W[x + static_cast<std::size_t>(nd) * y] * weight_dam(x);
        if (dam) { for (int r = 0; r < 2; ++r) { const int a = al[2 * x + r]; if (a >= 0) X[a + static_cast<std::size_t>(k) * y] += wt[2 * x + r] * w; } }
        else for (int a = 0; a < k; ++a) X[a + static_cast<std::size_t>(k) * y] += p[a] * w;
      }
    std::vector<Real> acc(static_cast<std::size_t>(k) * k, 0.0);   // acc(a, b)
    for (int y = 0; y < ns; ++y) {
      const Real ws = weight_sire(y);
      auto add = [&](int b, double tb) { for (int a = 0; a < k; ++a) acc[a + static_cast<std::size_t>(k) * b] += tb * ws * X[a + static_cast<std::size_t>(k) * y]; };
      if (sire) { for (int r = 0; r < 2; ++r) { const int b = al[2 * y + r]; if (b >= 0) add(b, wt[2 * y + r]); } }
      else for (int b = 0; b < k; ++b) add(b, p[b]);
    }
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) result[lookup[a + k * b]] += acc[a + static_cast<std::size_t>(k) * b];
  } else {  // target is a child: distribute parental pair weight through gametes
    const auto& t = *L.t; const auto& lookup = *L.lookup; const auto& p = *L.p;
    std::vector<Real> X(static_cast<std::size_t>(k) * ns, 0.0);  // X(a, y)
    for (int y = 0; y < ns; ++y)
      for (int a = 0; a < k; ++a) {
        Real acc = 0.0;
        for (int x = 0; x < nd; ++x) acc += (dam ? t[x + G * a] : p[a]) * W[x + static_cast<std::size_t>(nd) * y] * weight_dam(x);
        X[a + static_cast<std::size_t>(k) * y] = acc;
      }
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) {
        Real acc = 0.0;
        for (int y = 0; y < ns; ++y) acc += (sire ? t[y + G * b] : p[b]) * weight_sire(y) * X[a + static_cast<std::size_t>(k) * y];
        result[lookup[a + k * b]] += acc;
      }
  }
  out.log_cells.clear();
  if constexpr(std::is_same_v<Real,double>){out.cells=std::move(result);out.scale=scale;normalize(out);}
  else {std::vector<double> logs;for(auto v:result)logs.push_back(v.log);store_logs(out,std::move(logs),scale);}

}

template<class Outgoing>
inline void family_message_with(const Forest& F,const LocusView& L,int f,int target,Outgoing&& outgoing,Scaled& out){
  try{family_message_values<double>(F,L,f,target,outgoing,out);}
  catch(const StorageRange&){family_message_values<LogPositive>(F,L,f,target,outgoing,out);}
}
inline void family_message(const Forest& F, const LocusView& L, const LocusMessages& M,
                           int f, int target, Scaled& out) {
  family_message_with(F, L, f, target, [&](int e) { return &M.to_family(e, L); }, out);
}

// Computes a stale message after its stale inputs, without recursion.






void compute_messages(const Forest& F, const LocusView& L, LocusMessages& M);

// ---------------------------------------------------------------------------
// One visit's view of a stored forest. The visit removes the focal's own birth
// link (its parents are what the questions decide); every stored message whose
// source side contains that link is repaired on demand along the path to the
// question that reads it, and cached for the visit. With no link to remove the
// stored messages are used as they are.
struct Context {
  const Forest* F = nullptr;
  bool latent_brood_exclusion = false;
  bool ancestral_reuse = false;
  mutable std::vector<signed char> focal_dependency;
  int focal = 0, removed = -1, lower = -1, focal_component = -1;
  bool focal_below = false;   // focal side is the subtree under `lower`
  struct Repair {
    std::unordered_map<int, Scaled> to_family, to_individual;
    bool side_ready = false; double side_log_z = 0.0;
  };
  std::vector<Repair> repairs;
  std::size_t repaired = 0;
  void reset(const Forest& forest, int focal_id, int loci, bool exclude_latent_brood=false,bool reuse_ancestral=false);
  int child_edge(int id) const { const int e = F->child_edge[id]; return e == removed ? -1 : e; }
  // Side of the removed link holding the focal's former birth family.
  bool family_side(int x) const;
  // Component in the visit's context; the family side gets id F->components.
  int component(int x) const { return family_side(x) ? F->components : F->node_component[x]; }
  bool depends_on_focal(int v) const;
  bool dirty(int e, bool to_family) const;
  const Scaled& to_individual(int e, int locus, const LocusView& L, const LocusMessages& M) {
    if(M.ancestral&&F->edges[e].role<2)return M.to_individual(e,L);
    if(latent_brood_exclusion && latent_brood_upward(*F,e))return M.to_individual(e,L);
    if (!dirty(e, false)) return M.to_individual(e, L);
    auto& memo = repairs[locus].to_individual;
    const int key=ancestral_reuse?F->edges[e].family:e;
    if (auto at = memo.find(key); at != memo.end()) return at->second;
    Scaled out;
    family_message_with(*F, L, F->edges[e].family, e,
      [&](int x) -> const Scaled* { return x == removed || (M.ancestral&&F->edges[x].role==2) ? nullptr : &to_family(x, locus, L, M); }, out);
    ++repaired;
    return memo.emplace(key, std::move(out)).first->second;
  }
  const Scaled& to_family(int e, int locus, const LocusView& L, const LocusMessages& M) {
    if(M.ancestral&&F->edges[e].role==2)return M.to_family(e,L);
    if (!dirty(e, true)) return M.to_family(e, L);
    auto& memo = repairs[locus].to_family;
    const int key=ancestral_reuse?F->edges[e].individual:e;
    if (auto at = memo.find(key); at != memo.end()) return at->second;
    Scaled out;
    const int v = F->edges[e].individual;
    individual_message_with(*F, L, v, e, child_edge(v) < 0,
      [&](int x) -> const Scaled* { return x == removed || (M.ancestral&&F->edges[x].role<2) ? nullptr : &to_individual(x, locus, L, M); }, out);
    ++repaired;
    return memo.emplace(key, std::move(out)).first->second;
  }
  // log Z of the family side when no question individual lies there.
  double family_side_log_z(int locus, const LocusView& L, const LocusMessages& M);
};

// ---------------------------------------------------------------------------
// Small exact variable elimination. Genotype variables have G states and
// gamete (allele) variables k states; the order greedily minimizes the table
// size of each elimination and is computed once per hypothesis and allele count.

// Greedy elimination orders scored by the work they cause: a step costs
// cells x states x tables, or cells x (tables - 1) when a pairing indicator
// selects the eliminated state. Three greedy rules (smallest clique, least
// fill, least step work) are simulated and the cheapest order is kept.

}
