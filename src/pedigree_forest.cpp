#include "pedigree_forest.h"
namespace pap_internal::pedigree_likelihood {

auto Forest::build(const std::vector<int>& dam, const std::vector<int>& sire, const std::vector<int>& seeds) -> void  {
    stable_index.clear(); stable_epoch.clear();
    build_families(dam, sire, seeds);
    for (int& e : stable_epoch) e = 0;
  }

auto Forest::build_families(const std::vector<int>& dam, const std::vector<int>& sire,
                      const std::vector<int>& seeds) -> void  {
    cached_cut_flags.clear();
    n = static_cast<int>(dam.size()) - 1;
    families.clear(); family_index.clear();
    family_of.assign(n + 1, -1);
    current_dam.assign(n + 1, 0); current_sire.assign(n + 1, 0);
    for (int c = 1; c <= n; ++c) {
      int d = dam[c] > 0 ? dam[c] : 0, s = sire[c] > 0 ? sire[c] : 0;
      current_dam[c] = d; current_sire[c] = s;
      if (!d && !s) continue;
      auto at = family_index.emplace(std::make_pair(d, s), static_cast<int>(families.size()));
      if (at.second) { Family f; f.dam = d; f.sire = s; families.push_back(f); }
      families[at.first->second].children.push_back(c);
      family_of[c] = at.first->second;
    }
    // Current-pedigree-only ordering: families are numbered by the first
    // child ID in build's ascending scan. Seeds never affect cut selection.
    (void)seeds;
    const int F = static_cast<int>(families.size());
    family_stable.assign(F, -1);
    for (int f = 0; f < F; ++f) {
      auto at = stable_index.emplace(std::make_pair(families[f].dam, families[f].sire),
                                     static_cast<int>(stable_epoch.size()));
      if (at.second) stable_epoch.push_back(-1);
      family_stable[f] = at.first->second;
    }
    std::vector<int> order_f(F); std::iota(order_f.begin(), order_f.end(), 0);
    // Union-find over individuals plus family nodes (index n+1+f).
    std::vector<int> link(n + 2 + F);
    std::iota(link.begin(), link.end(), 0);
    auto find = [&](int x) { while (link[x] != x) { link[x] = link[link[x]]; x = link[x]; } return x; };
    edges.clear(); cuts.clear(); cut_incidences = 0;
    family_parent_edge.assign(F, {{-1, -1}});
    family_child_edges.assign(F, {});
    individual_edges.assign(n + 1, {});
    child_edge.assign(n + 1, -1);
    for (int f : order_f) {
      for (int role = 0; role < 2; ++role) {
        const int p = role == 0 ? families[f].dam : families[f].sire;
        if (p > 0) attach(f, p, role, link, find);
      }
      for (int c : families[f].children) attach(f, c, 2, link, find);
    }
    revision = 0;ancestry_revision.assign(n+1,0);
    stale_to_family.assign(edges.size(), 0); stale_to_individual.assign(edges.size(), 0);
    reorient();
  }

auto Forest::update(const std::vector<int>& dam, const std::vector<int>& sire) -> bool  {
    std::vector<int> moved;
    for(int c=1;c<=n;++c)if(std::max(0,dam[c])!=current_dam[c]||std::max(0,sire[c])!=current_sire[c])moved.push_back(c);
    if(moved.empty())return false;
    auto cut_keys=[&](){std::set<std::array<int,4>> out;for(const auto& c:cuts){const auto& f=families[c[1]];out.insert({c[0],f.dam,f.sire,c[2]});}return out;};
    auto old_cuts=cut_keys();auto prior_ancestry=ancestry_revision;const int next=revision+1;
    const size_t old_edges=edges.size();
    const std::vector<int> old_stable=family_stable;
    build_families(dam,sire,std::vector<int>{});revision=next;++updates;
    {std::vector<unsigned char> now(stable_epoch.size());for(int id:family_stable)now[id]=1;
     vanished.clear();vanished_at=next;
     for(int id:old_stable)if(!now[id]){stable_epoch[id]=-1;vanished.push_back(id);}
     for(int id:family_stable)if(stable_epoch[id]<0)stable_epoch[id]=next;}
    removed_edges+=old_edges;added_edges+=edges.size();
    stale_to_family.assign(edges.size(),next);stale_to_individual.assign(edges.size(),next);
    stale_marks+=2*edges.size();ancestry_revision=std::move(prior_ancestry);
    if(old_cuts!=cut_keys())std::fill(ancestry_revision.begin(),ancestry_revision.end(),next);
    else {
      std::vector<std::vector<int>> children(n+1);for(int c=1;c<=n;++c)for(int p:{current_dam[c],current_sire[c]})if(p>0)children[p].push_back(c);
      std::vector<unsigned char> seen(n+1);for(int c:moved)seen[c]=1;
      for(size_t j=0;j<moved.size();++j){int v=moved[j];ancestry_revision[v]=next;for(int c:children[v])if(!seen[c]){seen[c]=1;moved.push_back(c);}}
    }
    return true;
  }

auto Forest::ancestry_cut_flags() const -> std::vector<unsigned char>  {
    std::vector<unsigned char> flags(n+1),state(n+1);
    for(const auto& c:cuts){if(c[2]==2)flags[c[0]]=1;else for(int child:families[c[1]].children)flags[child]=1;}
    std::function<void(int)> visit=[&](int v){if(state[v]==2)return;if(state[v]==1)throw std::logic_error("assigned ancestry cycle");state[v]=1;
      for(int p:{current_dam[v],current_sire[v]})if(p>0){visit(p);flags[v]|=flags[p];}state[v]=2;};
    for(int v=1;v<=n;++v)visit(v);return flags;
  }

auto Forest::reorient() -> void  {
    const int F = static_cast<int>(families.size());
    const int nodes = 2 * std::max(n, F) + 2;
    component.assign(n + 1, -1);
    order.clear(); up_node.assign(nodes, -1); up_edge.assign(nodes, -1); depth.assign(nodes, 0); root.clear(); components = 0;
    std::vector<unsigned char> seen(nodes, 0);
    for (int start = 1; start <= n; ++start) {
      if (seen[ind(start)]) continue;
      root.push_back(start);
      std::deque<int> q{ind(start)}; seen[ind(start)] = 1;
      while (!q.empty()) {
        int x = q.front(); q.pop_front(); order.push_back(x);
        if (x % 2 == 0) component[x / 2] = components;
        visit_neighbors(x, [&](int y, int e) {
          if (seen[y]) return;
          seen[y] = 1; up_node[y] = x; up_edge[y] = e; depth[y] = depth[x] + 1; q.push_back(y);
        });
      }
      ++components;
    }
    node_component.assign(nodes, -1);
    for (int x : order) node_component[x] = x % 2 == 0 ? component[x / 2] : node_component[up_node[x]];
    std::vector<std::vector<int>> below(nodes);
    for (int x : order) if (up_node[x] >= 0) below[up_node[x]].push_back(x);
    tin.assign(nodes, -1); tout.assign(nodes, -1);
    int clock = 0;
    std::vector<std::pair<int, std::size_t>> stack;
    for (int r : root) {
      stack.push_back({ind(r), 0}); tin[ind(r)] = clock++;
      while (!stack.empty()) {
        auto& [x, next] = stack.back();
        if (next < below[x].size()) { const int y = below[x][next++]; tin[y] = clock++; stack.push_back({y, 0}); }
        else { tout[x] = clock - 1; stack.pop_back(); }
      }
    }
  }

auto Forest::path(int x, int y) const -> std::vector<int>  {
    std::vector<int> a, b;
    while (depth[x] > depth[y]) { a.push_back(x); x = up_node[x]; }
    while (depth[y] > depth[x]) { b.push_back(y); y = up_node[y]; }
    while (x != y) { a.push_back(x); b.push_back(y); x = up_node[x]; y = up_node[y]; }
    a.push_back(x);
    a.insert(a.end(), b.rbegin(), b.rend());
    return a;
  }



}
