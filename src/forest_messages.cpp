#include "forest_messages.h"
namespace pap_internal::pedigree_likelihood {
void store_logs(Scaled& out, std::vector<double> logs, double scale) {
  double top=-INFINITY;for(double x:logs)top=std::max(top,x);
  out.scale=scale;out.log_cells.clear();out.cells.clear();
  if(top==-INFINITY){out.cells.assign(logs.size(),0.);out.scale=-INFINITY;return;}
  if(!std::isfinite(top))throw std::runtime_error("M_ped2 invalid log factor");
  out.scale+=top;
  bool wide=false;for(double& x:logs){x-=top;if(std::isfinite(x)&&x<std::log(std::numeric_limits<double>::min()))wide=true;}
  if(wide){std::vector<double>().swap(out.cells);out.log_cells=std::move(logs);}
  else {std::vector<double>().swap(out.log_cells);out.cells.reserve(logs.size());for(double x:logs)out.cells.push_back(std::exp(x));}
}

void LocusMessages::ensure(int e0, bool side0, const LocusView& L) const {
  if(shared()) {
    sync_shared();if(unit.size()!=std::size_t(L.G)){unit.cells.assign(L.G,1.);unit.log_cells.clear();unit.scale=0.;}
  }
  if (current(e0, side0)) return;
  const std::size_t E = F->edges.size();
  if (!shared() && family_cells.size() < E) {
    family_cells.resize(E); individual_cells.resize(E); family_stamp.resize(E, -1); individual_stamp.resize(E, -1);
  }
  std::vector<std::pair<int, bool>> stack{{e0, side0}};
  while (!stack.empty()) {
    const auto [e, side] = stack.back();
    if (current(e, side)) { stack.pop_back(); continue; }
    if(ancestral && (side ? F->edges[e].role==2 : F->edges[e].role<2)) {
      (side?family_cells:individual_cells)[e].log_cells.clear();
      (side?family_cells:individual_cells)[e].cells.assign(L.G,1.0);
      (side?family_cells:individual_cells)[e].scale=0.0;
      (side?family_stamp:individual_stamp)[e]=F->revision;stack.pop_back();continue;
    }
    if(!side && latent_brood_exclusion && latent_brood_upward(*F,e)) {
      individual_cells[e].log_cells.clear();individual_cells[e].cells.assign(L.G,1.0);individual_cells[e].scale=0.0;
      individual_stamp[e]=F->revision;stack.pop_back();continue;
    }
    bool ready = true;
    if (side) {
      for (int x : F->individual_edges[F->edges[e].individual])
        if (x != e && (!ancestral || F->edges[x].role==2) && !current(x, false)) { stack.push_back({x, false}); ready = false; }
    } else {
      const int f = F->edges[e].family;
      for (int x : F->family_parent_edge[f]) if (x >= 0 && x != e && !current(x, true)) { stack.push_back({x, true}); ready = false; }
      for (int x : F->family_child_edges[f]) if (!ancestral && x != e && !current(x, true)) { stack.push_back({x, true}); ready = false; }
    }
    if (!ready) continue;
    if (side) individual_message(*F, L, *this, F->edges[e].individual, e, cells_for(e,side));
    else if(ancestral)family_message_with(*F,L,F->edges[e].family,e,
      [&](int x)->const Scaled*{return F->edges[x].role<2?&to_family(x,L):nullptr;},cells_for(e,side));
    else family_message(*F, L, *this, F->edges[e].family, e, cells_for(e,side));
    stamp_for(e,side) = F->revision;
    if (worker_message_count) ++*worker_message_count; else ++F->messages_computed;
    stack.pop_back();
  }
}

void LocusMessages::totals(const LocusView& L) const {
  if (totals_revision == F->revision) return;
  component_log_z.assign(F->components, 0.0);
  for (int c = 0; c < F->components; ++c) {
    int r = F->root[c];
    Scaled belief; emission(L, r, belief);
    if (F->child_edge[r] < 0) multiply_prior(belief,L);
    for (int e : F->individual_edges[r]) multiply_into(belief, to_individual(e, L));
    if(belief.wide())component_log_z[c]=log_sum(belief)+belief.scale;
    else {double sum=0.;for(double v:belief.cells)sum+=v;component_log_z[c]=sum>0?std::log(sum)+belief.scale:-INFINITY;}
  }
  total_log_z = 0.0;
  for (double z : component_log_z) total_log_z += z;
  totals_revision = F->revision;
}

void compute_messages(const Forest& F, const LocusView& L, LocusMessages& M) {
  M.attach(F);
  if(!M.ancestral)M.totals(L);
}
auto LocusMessages::attach(const Forest& forest) -> void  {
    if (F != &forest) {
      ancestor_cells.clear(); birth_cells.clear(); ancestor_stamp.clear();
      birth_stamp.clear(); shared_revision = -1;
    }
    F = &forest; totals_revision = -1;
  }

auto LocusMessages::sync_shared() const -> void  {
    if (shared_revision == F->revision) return;
    if (ancestor_cells.size() != std::size_t(F->n+1)) {
      ancestor_cells.resize(F->n+1); ancestor_stamp.assign(F->n+1,-1);
    }
    // Birth messages are indexed by stable family id; new ids start unset.
    // Release the messages of families that no longer exist.
    if (shared_revision >= 0 && shared_revision + 1 == F->revision && F->vanished_at == F->revision) {
      for (int id : F->vanished) if (std::size_t(id) < birth_cells.size()) {
        birth_cells[id] = Scaled(); birth_stamp[id] = -1;
      }
    } else {
      for (std::size_t id = 0; id < birth_cells.size() && id < F->stable_epoch.size(); ++id)
        if (F->stable_epoch[id] < 0 && birth_stamp[id] >= 0) { birth_cells[id] = Scaled(); birth_stamp[id] = -1; }
    }
    if (birth_cells.size() < F->stable_epoch.size()) {
      birth_cells.resize(F->stable_epoch.size()); birth_stamp.resize(F->stable_epoch.size(), -1);
    }
    shared_revision=F->revision;
  }

auto LocusMessages::current(int e, bool to_family_side) const -> bool  {
    if(shared()) {
      sync_shared();if(suppressed(e,to_family_side))return unit.size()>0;
      const auto& edge=F->edges[e];const auto& f=F->families[edge.family];
      const int needed=to_family_side?F->ancestry_revision[edge.individual]:
        std::max(F->ancestry_revision[f.dam],F->ancestry_revision[f.sire]);
      const int stamp=stamp_for(e,to_family_side);
      // A birth message is kept only while its family has existed since it
      // was computed, as the former per-revision remapping kept it.
      if(!to_family_side&&stamp<F->stable_epoch[F->family_stable[edge.family]])return false;
      return stamp>=0&&stamp>=needed;
    }
    const auto& stamp = to_family_side ? family_stamp : individual_stamp;
    if (e >= static_cast<int>(stamp.size()) || stamp[e] < 0) return false;
    return stamp[e] >= (to_family_side ? F->stale_to_family[e] : F->stale_to_individual[e]);
  }

auto Context::reset(const Forest& forest, int focal_id, int loci, bool exclude_latent_brood,bool reuse_ancestral) -> void  {
    F = &forest; focal = focal_id; removed = lower = focal_component = -1; focal_below = false;
    latent_brood_exclusion=exclude_latent_brood;
    ancestral_reuse=forest.ancestral&&reuse_ancestral;
    focal_dependency.clear();
    repairs.assign(loci, Repair{}); repaired = 0;
    if (focal_id <= 0 || focal_id > forest.n) return;
    focal_component = forest.component[focal_id];
    removed = forest.child_edge[focal_id];
    if (removed < 0) return;
    lower = forest.lower(removed);
    focal_below = lower == Forest::ind(focal_id);
  }

auto Context::family_side(int x) const -> bool  {
    if (removed < 0 || F->node_component[x] != focal_component) return false;
    return focal_below ? !F->within(x, lower) : F->within(x, lower);
  }

auto Context::depends_on_focal(int v) const -> bool  {
    if(!v)return false;
    if(focal_dependency.empty())focal_dependency.assign(F->n+1,-1);
    std::vector<int> stack{v};
    while(!stack.empty()) {
      const int id=stack.back();
      if(focal_dependency[id]>=0){stack.pop_back();continue;}
      if(id==focal){focal_dependency[id]=1;stack.pop_back();continue;}
      const int birth=F->child_edge[id];
      if(birth<0){focal_dependency[id]=0;stack.pop_back();continue;}
      bool ready=true,affected=false;
      for(int e:F->family_parent_edge[F->edges[birth].family])if(e>=0){
        const int parent=F->edges[e].individual;
        if(focal_dependency[parent]<0){stack.push_back(parent);ready=false;}
        else affected|=focal_dependency[parent]>0;
      }
      if(ready){focal_dependency[id]=affected;stack.pop_back();}
    }
    return focal_dependency[v]>0;
  }

auto Context::dirty(int e, bool to_family) const -> bool  {
    if (removed < 0 || e == removed || F->node_component[Forest::ind(F->edges[e].individual)] != focal_component) return false;
    if(ancestral_reuse) {
      if(to_family)return F->edges[e].role<2&&depends_on_focal(F->edges[e].individual);
      if(F->edges[e].role<2)return false;
      for(int p:F->family_parent_edge[F->edges[e].family])if(p>=0&&depends_on_focal(F->edges[p].individual))return true;
      return false;
    }
    const int source = to_family ? Forest::ind(F->edges[e].individual) : Forest::fam(F->edges[e].family);
    const int below = F->lower(e);
    return source == below ? F->within(lower, below) : !F->within(lower, below);
  }

auto Context::family_side_log_z(int locus, const LocusView& L, const LocusMessages& M) -> double  {
    if (removed < 0) return 0.0;
    auto& r = repairs[locus];
    if (r.side_ready) return r.side_log_z;
    const int f = F->edges[removed].family;
    int e = -1;
    for (int x : F->family_parent_edge[f]) if (x >= 0) { e = x; break; }
    if (e < 0) for (int x : F->family_child_edges[f]) if (x != removed) { e = x; break; }
    double z = 0.0;
    if (e >= 0) {
      const Scaled& a = to_individual(e, locus, L, M);
      const Scaled& b = to_family(e, locus, L, M);
      double sum = -INFINITY;
      if(a.wide()||b.wide()){for(int g=0;g<L.G;++g)sum=log_add(sum,a.log_at(g)+b.log_at(g));}
      else {double acc=0.;for(int g=0;g<L.G;++g)acc+=a.cells[g]*b.cells[g];sum=std::log(acc);}
      z = std::isfinite(sum) && a.scale != -std::numeric_limits<double>::infinity() && b.scale != -std::numeric_limits<double>::infinity()
        ? sum + a.scale + b.scale : -std::numeric_limits<double>::infinity();
    }
    r.side_ready = true; r.side_log_z = z;
    return z;
  }


const Scaled& LocusMessages::emission_of(int id, const LocusView& L) const {
  if (emissions.size() <= static_cast<std::size_t>(id)) emissions.resize(L.phi->size());
  Scaled& out = emissions[id];
  if (out.size()==0) {
    emission(L,id,out);
  }
  return out;
}
}
