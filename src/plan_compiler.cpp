#include "plan_compiler.h"
namespace pap_internal::pedigree_likelihood {
std::vector<int> choose_order(const std::vector<std::vector<int>>& scopes, const std::vector<unsigned char>& pairing,
                                     const std::vector<int>& dom, int retained, int* width) {
  const int variables = static_cast<int>(dom.size());
  std::vector<int> best_order;
  double best_cost = std::numeric_limits<double>::infinity();
  int best_width = 0;
  for (int rule = 0; rule < 3; ++rule) {
    // Tables containing each variable, their pairing flags, and adjacency.
    std::vector<std::vector<int>> tables(scopes);
    std::vector<unsigned char> selecting(pairing), alive(tables.size(), 1), done(variables, 0);
    std::vector<std::vector<int>> holding(variables);
    std::vector<std::vector<unsigned char>> adjacent(variables, std::vector<unsigned char>(variables, 0));
    for (std::size_t t = 0; t < tables.size(); ++t)
      for (int u : tables[t]) {
        holding[u].push_back(static_cast<int>(t));
        for (int x : tables[t]) if (x != u) adjacent[u][x] = 1;
      }
    std::vector<int> order;
    double cost = 0.0;
    int w = 0;
    std::vector<int> joint;
    auto gather = [&](int v, int& n, bool& select) {
      joint.clear(); n = 0; select = false;
      for (int t : holding[v]) if (alive[t]) { ++n; select = select || selecting[t]; }
      for (int u = 0; u < variables; ++u) if (!done[u] && adjacent[v][u]) joint.push_back(u);
    };
    for (int round = retained; round < variables; ++round) {
      int pick = -1; double pick_score = std::numeric_limits<double>::infinity();
      for (int v = retained; v < variables; ++v) {
        if (done[v]) continue;
        int n; bool select;
        gather(v, n, select);
        double cells = 1.0;
        for (int u : joint) cells *= dom[u];
        double score;
        if (rule == 0) score = cells * dom[v];
        else if (rule == 1) {
          double fill = 0.0;
          for (std::size_t x = 0; x < joint.size(); ++x)
            for (std::size_t y = x + 1; y < joint.size(); ++y)
              if (!adjacent[joint[x]][joint[y]]) fill += static_cast<double>(dom[joint[x]]) * dom[joint[y]];
          score = fill * 1e6 + cells * dom[v];
        } else score = (n == 0 ? 0.0 : select && n >= 2 ? cells * (n - 1) : cells * dom[v] * n) + cells;
        if (score < pick_score - 1e-9) { pick_score = score; pick = v; }
      }
      int n; bool select;
      gather(pick, n, select);
      done[pick] = 1; order.push_back(pick);
      if (n == 0) continue;
      double cells = 1.0;
      for (int u : joint) cells *= dom[u];
      cost += select && n >= 2 ? cells * (n - 1) : cells * dom[pick] * n;
      w = std::max(w, static_cast<int>(joint.size()) + 1);
      for (int t : holding[pick]) alive[t] = 0;
      const int t = static_cast<int>(tables.size());
      tables.push_back(joint); selecting.push_back(0); alive.push_back(1);
      for (int u : joint) {
        holding[u].push_back(t);
        for (int x : joint) if (x != u) adjacent[u][x] = 1;
      }
    }
    if (cost < best_cost) { best_cost = cost; best_order = order; best_width = w; }
  }
  if (width) *width = std::max(*width, best_width);
  return best_order;
}

void build_transmission(const LocusView& L, Transmission& T, bool implicit ) {
  const int G = L.G, k = L.k;
  const auto& t = *L.t; const auto& lookup = *L.lookup; const auto& p = *L.p;
  const std::size_t g = static_cast<std::size_t>(G);
  T.implicit = implicit;
  T.two_support=sparse_gametes(L,T.sparse_allele,T.sparse_weight);
  if(T.two_support)for(double w:T.sparse_weight)if(!std::isfinite(w)||w<0)T.two_support=false;
  T.alleles.assign(G, {{-1,-1}});
  for (int a=0;a<k;++a) for(int b=a;b<k;++b) T.alleles[lookup[a+k*b]]={{a,b}};
  T.pairing.clear(); T.one.clear();
  if (!implicit) T.pairing.assign(g * k * k, 0.0);
  for (int a = 0; a < k; ++a)
    for (int b = 0; !implicit && b < k; ++b) T.pairing[lookup[a + k * b] + g * (a + static_cast<std::size_t>(k) * b)] = 1.0;
  T.partner.assign(g * k, -1);
  for (int a = 0; a < k; ++a)
    for (int b = 0; b < k; ++b) T.partner[lookup[a + k * b] + g * b] = a;
  if (implicit) return;
  T.one.assign(g * g, 0.0);
  for (int d = 0; d < G; ++d)
    for (int a = 0; a < k; ++a) {
      const double x = t[d + G * a];
      if (x == 0.0) continue;
      for (int b = 0; b < k; ++b) T.one[lookup[a + k * b] + g * d] += x * p[b];
    }
}

bool cyclic_plan(const Plan& p) {
  const int V=p.genotype_variables?p.genotype_variables:p.variables;
  std::vector<int> degree(V);std::vector<std::vector<int>> children(V);
  for(const auto& b:p.births)if(!b.excluded_brood)for(int s=1;s<3;++s)if(b.vars[s]>=0){
    if(b.vars[0]<0||b.vars[0]>=V||b.vars[s]>=V)
      throw std::out_of_range("M_ped2 birth outside genotype-variable domain");
    ++degree[b.vars[0]];children[b.vars[s]].push_back(b.vars[0]);
  }
  std::vector<int> ready;for(int v=0;v<V;++v)if(!degree[v])ready.push_back(v);
  size_t visited=0;for(size_t j=0;j<ready.size();++j){++visited;for(int c:children[ready[j]])if(!--degree[c])ready.push_back(c);}
  return visited!=size_t(V);
}

void require_acyclic_plan(const Plan& p) {
  if(cyclic_plan(p))throw ::relationship_topology::topology_detail::Failure(
    ::relationship_topology::TopologyStatus::Inapplicable,"M_ped2 inheritance plan contains a directed cycle");
}

std::uint64_t mix_hash(std::uint64_t h, std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  x ^= x >> 31; x *= 0xbf58476d1ce4e5b9ULL; x ^= x >> 27; x *= 0x94d049bb133111ebULL; x ^= x >> 31;
  return h ^ x;
}

std::vector<long long> plan_key(const Plan& p, std::size_t* cyclic_plans ) {
  if(cyclic_plan(p)){if(cyclic_plans)++*cyclic_plans;require_acyclic_plan(p);}
  const int V = p.genotype_variables;
  // Reversible trees: births with one parent (and a population mate) form
  // trees; when a tree's top is a founder carrying only its prior, the
  // parent-offspring joint is symmetric under Hardy-Weinberg priors, so the
  // prior may sit at any node. Those links are keyed undirected and the
  // prior as a mark on the tree.
  std::vector<unsigned char> undirected(p.births.size(), 0), dropped(p.unaries.size(), 0);
  std::vector<int> parent_of(V, -1), root(V, -1), tops;
  if (reversible_keys()) {
    std::vector<int> birth_of(V, -1), priors(V, 0), prior_at(V, -1);
    std::vector<unsigned char> two_parent(V, 0);
    for (std::size_t i = 0; i < p.births.size(); ++i) {
      const auto& b = p.births[i];
      if ((b.vars[1] >= 0) != (b.vars[2] >= 0)) { parent_of[b.vars[0]] = b.vars[1] >= 0 ? b.vars[1] : b.vars[2]; birth_of[b.vars[0]] = static_cast<int>(i); }
      else two_parent[b.vars[0]] = 1;
    }
    for (std::size_t i = 0; i < p.unaries.size(); ++i)
      if (p.unaries[i].source == Plan::Prior && p.unaries[i].var < V) { priors[p.unaries[i].var]++; prior_at[p.unaries[i].var] = static_cast<int>(i); }
    auto top_of = [&](int v) {
      // A terminating chain has at most V-1 edges. A V-edge walk must
      // revisit a genotype variable; reversible-tree keying cannot apply.
      for (int step = 0; step < V; ++step) {
        if (parent_of[v] < 0) return v;
        v = parent_of[v];
      }
      return -1;
    };
    bool cyclic = false;
    for (int v = 0; v < V; ++v) {
      root[v] = top_of(v);
      if (root[v] < 0) { cyclic = true; break; }
    }
    if (cyclic) {
      if (cyclic_plans) ++*cyclic_plans;
      throw ::relationship_topology::topology_detail::Failure(::relationship_topology::TopologyStatus::Inapplicable,"M_ped2 one-parent walk contains a directed cycle");
    } else {
      // Commit reversible marks only after checking the complete plan.
      // Cycles are rejected before any prior can be dropped or link marked.
      std::vector<unsigned char> used(V, 0);
      for (int c = 0; c < V; ++c) {
        if (birth_of[c] < 0) continue;
        const int t = root[c];
        if (birth_of[t] >= 0 || two_parent[t] || priors[t] != 1 || t < p.retained) continue;
        undirected[birth_of[c]] = 1;
        if (!used[t]) { used[t] = 1; tops.push_back(t); dropped[prior_at[t]] = 1; }
      }
    }
  }
  std::vector<std::uint64_t> colour(V, 0), next(V);
  for (int v = 0; v < V; ++v) colour[v] = mix_hash(1, v < p.retained ? static_cast<std::uint64_t>(v) + 1 : 0);
  std::vector<std::vector<std::uint64_t>> unary(V);
  for (std::size_t i = 0; i < p.unaries.size(); ++i) {
    const auto& u = p.unaries[i];
    if (!dropped[i] && u.var >= 0 && u.var < V) unary[u.var].push_back(mix_hash(mix_hash(2, u.source), static_cast<std::uint64_t>(u.index)));
  }
  for (int v = 0; v < V; ++v) {
    std::sort(unary[v].begin(), unary[v].end());
    for (auto x : unary[v]) colour[v] = mix_hash(colour[v], x);
  }
  auto c = [&](int v) -> std::uint64_t { return v < 0 ? 0x51ULL : colour[v]; };
  std::vector<std::vector<std::uint64_t>> seen(V);
  int classes = -1;
  for (int round = 0; round < V; ++round) {
    for (auto& x : seen) x.clear();
    for (std::size_t i = 0; i < p.births.size(); ++i) {
      const auto& b = p.births[i];
      const int ch = b.vars[0], d = b.vars[1], s = b.vars[2];
      if (undirected[i]) {
        const int par = d >= 0 ? d : s;
        seen[ch].push_back(mix_hash(6, c(par))); seen[par].push_back(mix_hash(6, c(ch)));
        continue;
      }
      const std::uint64_t lo = std::min(c(d), c(s)), hi = std::max(c(d), c(s));
      seen[ch].push_back(mix_hash(mix_hash(mix_hash(3, 0), lo), hi));
      if (d >= 0) seen[d].push_back(mix_hash(mix_hash(mix_hash(3, 1), c(ch)), c(s)));
      if (s >= 0) seen[s].push_back(mix_hash(mix_hash(mix_hash(3, 1), c(ch)), c(d)));
    }
    for (const auto& x : p.pairs)
      for (int r = 0; r < 2; ++r)
        if (x.vars[r] >= 0) seen[x.vars[r]].push_back(mix_hash(mix_hash(mix_hash(mix_hash(4, x.slot), r), c(x.vars[0])), c(x.vars[1])));
    for (const auto& x : p.pieces) {
      std::uint64_t all = mix_hash(5, x.slot);
      for (int v : x.vars) all = mix_hash(all, c(v));
      for (std::size_t r = 0; r < x.vars.size(); ++r) if (x.vars[r] >= 0) seen[x.vars[r]].push_back(mix_hash(all, r));
    }
    for (int v = 0; v < V; ++v) {
      std::sort(seen[v].begin(), seen[v].end());
      std::uint64_t h = colour[v];
      for (auto x : seen[v]) h = mix_hash(h, x);
      next[v] = h;
    }
    colour.swap(next);
    std::vector<std::uint64_t> distinct(colour);
    std::sort(distinct.begin(), distinct.end());
    const int now = static_cast<int>(std::unique(distinct.begin(), distinct.end()) - distinct.begin());
    if (now == classes) break;
    classes = now;
  }
  std::vector<int> order(V), label(V);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) { return colour[a] != colour[b] ? colour[a] < colour[b] : a < b; });
  for (int i = 0; i < V; ++i) label[order[i]] = i;
  auto L = [&](int v) -> long long { return v < 0 ? -1 : v < V ? label[v] : 1000000LL + v; };
  std::vector<std::array<long long, 3>> units, births;
  for (std::size_t i = 0; i < p.unaries.size(); ++i)
    if (!dropped[i]) units.push_back({{L(p.unaries[i].var), static_cast<long long>(p.unaries[i].source), p.unaries[i].index}});
  for (std::size_t i = 0; i < p.births.size(); ++i) {
    const auto& b = p.births[i];
    long long d = L(b.vars[1]), s = L(b.vars[2]);
    if (undirected[i]) {
      long long x = L(b.vars[0]), y = d >= 0 ? d : s;
      if (x > y) std::swap(x, y);
      births.push_back({{-2, x, y}});
      continue;
    }
    if (d > s) std::swap(d, s);
    births.push_back({{L(b.vars[0]), d, s}});
  }
  // Each reversible tree's prior mark sits at its smallest label.
  for (int t : tops) {
    long long smallest = L(t);
    for (int v = 0; v < V; ++v) {
      if (root[v] == t) smallest = std::min(smallest, L(v));
    }
    births.push_back({{-3, smallest, 0}});
  }
  std::vector<std::array<long long, 3>> pairs;
  for (const auto& x : p.pairs) pairs.push_back({{x.slot, L(x.vars[0]), L(x.vars[1])}});
  std::vector<std::vector<long long>> pieces;
  for (const auto& x : p.pieces) {
    std::vector<long long> row{x.slot};
    for (int v : x.vars) row.push_back(L(v));
    pieces.push_back(std::move(row));
  }
  std::sort(units.begin(), units.end()); std::sort(births.begin(), births.end());
  std::sort(pairs.begin(), pairs.end()); std::sort(pieces.begin(), pieces.end());
  std::vector<long long> k{p.limits.mib,p.limits.work_millions,p.limits.live_arena,p.limits.shared_piece_scratch,p.retained, V, p.family_side_untouched ? 1 : 0,
                           static_cast<long long>(units.size()), static_cast<long long>(births.size()),
                           static_cast<long long>(pairs.size()), static_cast<long long>(pieces.size()),
                           static_cast<long long>(p.touched_components.size())};
  for (const auto& x : units) k.insert(k.end(), x.begin(), x.end());
  for (const auto& x : births) k.insert(k.end(), x.begin(), x.end());
  for (const auto& x : pairs) k.insert(k.end(), x.begin(), x.end());
  for (const auto& x : pieces) { k.push_back(static_cast<long long>(x.size())); k.insert(k.end(), x.begin(), x.end()); }
  k.insert(k.end(), p.touched_components.begin(), p.touched_components.end());
  for (int v = 0; v < p.retained && v < V; ++v) k.push_back(label[v]);
  if(p.limits.total_work_limit)k.insert(k.end(),{-70015,1});
  if(p.limits.s37_refusal_repair)k.insert(k.end(),{-70010,1});
  if(p.limits.resource_planning){k.push_back(-70001);k.push_back(1);}
  if(p.limits.resource_min_alleles!=20){k.push_back(-70003);k.push_back(p.limits.resource_min_alleles);}
  if(p.ancestral)k.insert(k.end(),{-70009,1});
  if(p.barren_prior_normalizers||p.barren_population_normalizers)k.insert(k.end(),{-70014,p.barren_prior_normalizers,p.barren_population_normalizers});
  if(p.local_family_pairs)k.insert(k.end(),{-70006,1});
  if(p.latent_brood_exclusion)k.insert(k.end(),{-70007,1});
  for(const auto& b:p.births)if(b.excluded_brood)k.insert(k.end(),{-70008,L(b.vars[0]),L(b.vars[1]),L(b.vars[2])});
  return k;
}

void assign_gametes(Plan& plan, int count) {
  plan.genotype_variables = count;
  for (auto& b : plan.births) if (b.vars[1] >= 0 && b.vars[2] >= 0) b.gametes = {{count++, count++}};
  plan.variables = count;
  plan.allele.assign(count, 0);
  for (int v = plan.genotype_variables; v < count; ++v) plan.allele[v] = 1;
}

Plan finish_plan(Plan plan,const std::vector<long long>& name,int count,
                       std::size_t model_unaries,std::size_t model_births,PieceCache& pieces,
                       bool split,int max_touch,int focal_var) {
  plan.genotype_variables=count;require_acyclic_plan(plan);
  constexpr long long Steiner=2000000000LL;
  const auto limits=plan.limits;const bool factorized_one=plan.factorized_one;
  const bool high_allele_programs=plan.high_allele_programs;
  // Split forest-only variables into connected pieces; each becomes a cached
  // table over the model individuals it touches.
  auto internal = [&](int v) { return split && v >= 0 && name[v] >= Steiner; };
  std::vector<int> link(count);
  std::iota(link.begin(), link.end(), 0);
  auto find = [&](int x) { while (link[x] != x) { link[x] = link[link[x]]; x = link[x]; } return x; };
  auto join = [&](const std::vector<int>& vs) {
    int root = -1;
    for (int v : vs) if (internal(v)) { if (root < 0) root = find(v); else link[find(v)] = root; }
  };
  for (std::size_t i = model_births; i < plan.births.size(); ++i) join({plan.births[i].vars[0], plan.births[i].vars[1], plan.births[i].vars[2]});
  for (const auto& p : plan.pairs) join({p.vars[0], p.vars[1]});
  // A piece is cached only when its table is small: at most two touched model
  // individuals (G^2 cells). Wider pieces stay in the main elimination.
  std::map<int, std::set<int>> touch_of;
  auto touch = [&](const std::vector<int>& vs) {
    int root = -1;
    for (int v : vs) if (internal(v)) { root = find(v); break; }
    if (root < 0) return;
    for (int v : vs) if (v >= 0 && !internal(v)) touch_of[root].insert(v);
    touch_of[root];
  };
  for (std::size_t i = model_unaries; i < plan.unaries.size(); ++i) touch({plan.unaries[i].var});
  for (std::size_t i = model_births; i < plan.births.size(); ++i) touch({plan.births[i].vars[0], plan.births[i].vars[1], plan.births[i].vars[2]});
  for (const auto& p : plan.pairs) touch({p.vars[0], p.vars[1]});
  std::vector<unsigned char> pieced(count, 0);
  for (int v = 0; v < count; ++v) if (internal(v)) pieced[v] = touch_of[find(v)].size() <= static_cast<std::size_t>(max_touch);
  struct Piece { std::vector<Plan::Unary> unaries; std::vector<Plan::Birth> births; std::vector<Plan::Pair> pairs; std::set<long long> touch; std::set<long long> inner; };
  std::map<int, Piece> by_root;
  Plan main;
  main.ancestral=plan.ancestral;
  main.barren_removed=plan.barren_removed;main.barren_prior_normalizers=plan.barren_prior_normalizers;main.barren_population_normalizers=plan.barren_population_normalizers;
  main.off_scope_observations_before=plan.off_scope_observations_before;main.descendant_observations_before=plan.descendant_observations_before;main.off_scope_messages=plan.off_scope_messages;
  main.local_family_pairs=plan.local_family_pairs;
  main.latent_brood_exclusion=plan.latent_brood_exclusion;
  main.restored_family_pairs=plan.restored_family_pairs;
  main.redundant_births = plan.redundant_births; main.steiner_individuals = plan.steiner_individuals;
  main.family_terms = plan.family_terms; main.touched_components = plan.touched_components;
  main.family_side_untouched = plan.family_side_untouched;
  auto piece_of = [&](const std::vector<int>& vs) -> Piece* {
    for (int v : vs) if (v >= 0 && pieced[v]) return &by_root[find(v)];
    return nullptr;
  };
  auto note = [&](Piece& q, const std::vector<int>& vs) { for (int v : vs) if (v >= 0) (pieced[v] ? q.inner : q.touch).insert(name[v]); };
  for (std::size_t i = 0; i < plan.unaries.size(); ++i) {
    const auto& u = plan.unaries[i];
    if (Piece* q = i >= model_unaries ? piece_of({u.var}) : nullptr) { q->unaries.push_back(u); note(*q, {u.var}); }
    else main.unaries.push_back(u);
  }
  for (std::size_t i = 0; i < plan.births.size(); ++i) {
    const auto& b = plan.births[i];
    const std::vector<int> vs{b.vars[0], b.vars[1], b.vars[2]};
    if (Piece* q = i >= model_births ? piece_of(vs) : nullptr) { q->births.push_back(b); note(*q, vs); }
    else main.births.push_back(b);
  }
  for (const auto& p : plan.pairs) {
    const std::vector<int> vs{p.vars[0], p.vars[1]};
    if (Piece* q = piece_of(vs)) { q->pairs.push_back(p); note(*q, vs); }
    else main.pairs.push_back(p);
  }
  // Renumber the main plan without forest-only variables.
  std::vector<int> renumber(count, -1);
  int kept = 0;
  // retain_focal: the focal is variable 0 and is kept, so one elimination
  // serves every focal emission (raw numerator and denominator).
  if (focal_var >= 0) renumber[focal_var] = kept++;
  for (int v = 0; v < count; ++v) if (!pieced[v] && v != focal_var) renumber[v] = kept++;
  std::map<long long, int> by_name;
  for (int v = 0; v < count; ++v) by_name[name[v]] = v;
  auto main_var = [&](int v) { return v < 0 ? -1 : renumber[v]; };
  for (auto& u : main.unaries) u.var = main_var(u.var);
  for (auto& b : main.births) for (int& v : b.vars) v = main_var(v);
  for (auto& p : main.pairs) for (int& v : p.vars) v = main_var(v);
  for (auto& [root, q] : by_root) {
    // Piece plan: touched model variables first (sorted by name), then its own.
    Plan sub;sub.ancestral=plan.ancestral;
    std::map<long long, int> local_var;
    for (long long n : q.touch) local_var.emplace(n, static_cast<int>(local_var.size()));
    sub.limits=limits;sub.factorized_one = factorized_one;sub.high_allele_programs=high_allele_programs;
    sub.retained = static_cast<int>(q.touch.size());
    for (long long n : q.inner) local_var.emplace(n, static_cast<int>(local_var.size()));
    auto sub_var = [&](int v) { return v < 0 ? -1 : local_var.at(name[v]); };
    auto nm = [&](int v) { return v < 0 ? -1LL : name[v]; };
    std::vector<std::vector<long long>> terms;
    for (auto u : q.unaries) { terms.push_back({1, nm(u.var), u.source, u.index}); u.var = sub_var(u.var); sub.unaries.push_back(u); }
    for (auto b : q.births) { terms.push_back({2, nm(b.vars[0]), nm(b.vars[1]), nm(b.vars[2])});
      if(b.excluded_brood)terms.back().push_back(-70008);
      for (int& v : b.vars) v = sub_var(v); sub.births.push_back(b); }
    for (auto p : q.pairs) { terms.push_back({3, p.slot, nm(p.vars[0]), nm(p.vars[1])}); for (int& v : p.vars) v = sub_var(v); sub.pairs.push_back(p); }
    std::sort(terms.begin(), terms.end());
    std::vector<long long> key(q.touch.begin(), q.touch.end());
    key.insert(key.begin(), {limits.mib,limits.work_millions,limits.live_arena,limits.shared_piece_scratch,static_cast<long long>(q.touch.size())});
    for (const auto& t : terms) { key.push_back(-1); key.insert(key.end(), t.begin(), t.end()); }
    if(limits.total_work_limit)key.insert(key.end(),{-70015,1});
    if(limits.s37_refusal_repair)key.insert(key.end(),{-70010,1});
    if(limits.resource_planning){key.push_back(-70001);key.push_back(1);}
    if(limits.resource_min_alleles!=20){key.push_back(-70003);key.push_back(limits.resource_min_alleles);}
    assign_gametes(sub, static_cast<int>(local_var.size()));
    sub.ready = true;
    Plan::PieceRef ref;
    ref.slot = pieces.slot(std::move(key), std::move(sub));
    for (long long n : q.touch) ref.vars.push_back(renumber[by_name.at(n)]);
    main.pieces.push_back(std::move(ref));
  }
  // Two-parent births factor through gametes: t(d, a) t(s, b) [c = ab].
  assign_gametes(main, kept);
  if (focal_var >= 0) main.retained = 1;
  main.limits=limits;main.factorized_one = factorized_one;main.high_allele_programs=high_allele_programs;
  main.ready = true;
  return main;
}

AncestralClosure ancestral_closure(const ::relationship_topology::TopologyModel& source,
                                         const Forest& F,const std::vector<int>& roots) { pap_internal::removed_path("ancestral_closure"); }

Plan ancestral_plan(const ::relationship_topology::TopologyModel& source,const Forest& F,
                           PieceCache& pieces,const std::vector<int>& roots,bool split,
                           int max_touch,bool retain_focal,bool factorized_one,
                           bool high_allele_programs,ResourceLimits limits) {
  auto closure=ancestral_closure(source,F,roots);
  auto& nodes=closure.nodes;auto& active=closure.active;
  Plan plan;plan.ancestral=true;plan.limits=limits;plan.factorized_one=factorized_one;
  plan.high_allele_programs=high_allele_programs;
  std::vector<int> var(nodes.size(),-1);std::vector<long long> name;
  const int original=int(source.nodes.size());
  for(int i=0;i<int(nodes.size());++i)if(active[i]){
    var[i]=int(name.size());name.push_back(nodes[i].typed_id>0?
      (i<original?1000000000LL:2000000000LL)+nodes[i].typed_id:i);
  }
  // Original model terms stay in the main question. Newly exposed ancestor
  // terms use the same bounded cached-piece splitter as full evidence.
  std::size_t model_unaries=0,model_births=0;
  for(int pass=0;pass<2;++pass){
    for(int i=0;i<int(nodes.size());++i)if(active[i]&&((i<original)==(pass==0))){
      const auto& n=nodes[i];
      if(n.typed_id>0&&n.observed&&!n.neutral_phi)plan.unaries.push_back({var[i],Plan::Emission,n.typed_id});
      const int d=n.has_birth&&n.parents[0]>=0?var[n.parents[0]]:-1;
      const int s=n.has_birth&&n.parents[1]>=0?var[n.parents[1]]:-1;
      if(d<0&&s<0)plan.unaries.push_back({var[i],Plan::Prior,0});
      else plan.births.push_back({{{var[i],d,s}}});
    }
    if(pass==0){model_unaries=plan.unaries.size();model_births=plan.births.size();}
  }
  plan.steiner_individuals=int(nodes.size())-original;
  plan.family_terms=int(plan.births.size());
  return finish_plan(std::move(plan),name,int(name.size()),model_unaries,0,
    pieces,split,max_touch,retain_focal?var[source.focal]:-1);
}

Plan plan_model(const ::relationship_topology::TopologyModel& source, const Forest& F, const Context& ctx, PairCache& cache,
                       PieceCache& pieces, bool defer_cut_births , bool split , bool retain_focal ,
                       int max_touch , bool factorized_one , bool high_allele_programs , ResourceLimits limits ,
                       bool local_family_pairs , bool latent_brood_exclusion ,
                       bool ancestral ,const std::vector<int>& ancestral_roots ,bool ancestral_reference ) {
  limits.validate();if(ancestral&&ancestral_reference)return ancestral_plan(source,F,pieces,ancestral_roots,split,max_touch,retain_focal,factorized_one,high_allele_programs,limits);
  Plan plan;plan.ancestral=ancestral;plan.limits=limits;plan.factorized_one=factorized_one;plan.high_allele_programs=high_allele_programs;
  plan.local_family_pairs=local_family_pairs;
  plan.latent_brood_exclusion=latent_brood_exclusion;
  // Under Hardy-Weinberg priors (reversible_keys) a latent founder whose only
  // role is to parent one model individual transmits population gametes:
  // sum_z prior(z) t(b | z) = p_b. Such founders are dropped and the birth
  // takes a population mate.
  ::relationship_topology::TopologyModel reduced;
  const ::relationship_topology::TopologyModel* model = &source;
  if (reversible_keys()) {
    std::vector<int> children(source.nodes.size(), 0);
    for (const auto& node : source.nodes) if (node.has_birth) for (int q : node.parents) if (q >= 0) ++children[q];
    bool any = false;
    for (std::size_t i = 0; i < source.nodes.size(); ++i) {
      const auto& node = source.nodes[i];
      if (node.typed_id <= 0 && !node.has_birth && children[i] == 1 && static_cast<int>(i) != source.focal &&
          static_cast<int>(i) != source.candidate) { any = true; break; }
    }
    if (any) {
      reduced = source;
      for (std::size_t i = 0; i < reduced.nodes.size(); ++i) {
        const auto& node = reduced.nodes[i];
        if (node.typed_id > 0 || node.has_birth || children[i] != 1 || static_cast<int>(i) == reduced.focal ||
            static_cast<int>(i) == reduced.candidate) continue;
        for (auto& other : reduced.nodes)
          if (other.has_birth) for (int& q : other.parents) if (q == static_cast<int>(i)) q = -1;
      }
      model = &reduced;
    }
  }
  const ::relationship_topology::TopologyModel& m = *model;
  const int nodes = static_cast<int>(m.nodes.size());
  // Ancestor-only observation scope. Incidental descendants/collaterals on
  // a forest Steiner path are not evidence for this parentage question.
  std::unordered_set<int> named,allowed;
  const int focal_id=m.nodes[m.focal].typed_id;
  if(ancestral){
    named.insert(focal_id);for(int id:ancestral_roots)if(id>0)named.insert(id);
    if(ancestral_roots.empty())for(const auto& node:m.nodes)if(node.typed_id>0)named.insert(node.typed_id);
    std::vector<int> todo(named.begin(),named.end());allowed=named;
    for(size_t j=0;j<todo.size();++j){int v=todo[j];if(v==focal_id)continue;
      for(int p:{F.current_dam.at(v),F.current_sire.at(v)})if(p>0&&allowed.insert(p).second)todo.push_back(p);}
  }
  auto observation=[&](int id){if(!ancestral||allowed.count(id))return true;
    ++plan.off_scope_observations_before;
    std::unordered_set<int> seen;std::vector<int> todo{id};bool descendant=false;
    for(size_t j=0;j<todo.size()&&!descendant;++j)for(int p:{F.current_dam.at(todo[j]),F.current_sire.at(todo[j])})if(p>0){
      if(named.count(p)){descendant=true;break;}if(seen.insert(p).second)todo.push_back(p);}
    plan.descendant_observations_before+=descendant;return false;};
  auto typed = [&](int i) { return i >= 0 ? m.nodes[i].typed_id : 0; };
  // A typed birth that repeats an uncut current birth is read from the stored
  // layer. With defer_cut_births (raw questions, which only change the
  // focal's parents) a current birth is deferred even when the forest cut it,
  // so every configuration sees the same cuts. Relationship hypotheses keep
  // cut current births exact: sibling and avuncular structures depend on them.
  std::vector<unsigned char> birth(nodes, 0);
  for (int i = 0; i < nodes; ++i) {
    const auto& node = m.nodes[i];
    if (!node.has_birth) continue;
    birth[i] = 1;
    const int c = node.typed_id;
    if (c <= 0 || c > F.n || i == m.focal) continue;
    const int f = F.family_of[c];
    if (f < 0) continue;
    bool same = true;
    for (int s = 0; s < 2; ++s) {
      const int p = node.parents[s], id = typed(p);
      if (p >= 0 && id <= 0) { same = false; break; }      // latent parent
      if ((s == 0 ? F.families[f].dam : F.families[f].sire) != id) { same = false; break; }
    }
    if (same && !defer_cut_births) {
      const int e = ctx.child_edge(c);
      same = e >= 0 && F.edges[e].family == f;
      for (int s = 0; s < 2 && same; ++s)
        if ((s == 0 ? F.families[f].dam : F.families[f].sire) > 0 && F.family_parent_edge[f][s] < 0) same = false;
    }
    if (same) { birth[i] = 0; ++plan.redundant_births; }
  }
  std::vector<unsigned char> active(nodes, 0);
  std::unordered_set<int> declared_births;
  if(latent_brood_exclusion)for(const auto& node:m.nodes)if(node.typed_id>0)declared_births.insert(node.typed_id);
  for (int i = 0; i < nodes; ++i) {
    if (!birth[i]) continue;
    active[i] = 1;
    for (int p : m.nodes[i].parents) if (p >= 0) active[p] = 1;
  }
  if (m.focal >= 0) active[m.focal] = 1;
  if (m.candidate >= 0) active[m.candidate] = 1;
  // Fixed declared evidence ports also survive an unrelated/null topology.
  if(ancestral){std::set<int> roots(ancestral_roots.begin(),ancestral_roots.end());
    for(int i=0;i<nodes;++i)if(m.nodes[i].typed_id>0&&roots.count(m.nodes[i].typed_id))active[i]=1;}
  // A cut source-parent port must not turn an exterior brood into a
  // population-mate message when both actual parents are declared typed
  // individuals of this question. Expose those existing variables in every
  // hypothesis (including U/raw null), then keep the brood as a bounded
  // joint factor. No new observations or pedigree individuals are selected.
  std::set<int> paired_families;
  if (local_family_pairs) {
    std::unordered_map<int,int> declared;
    for (int i=0;i<nodes;++i) if(m.nodes[i].typed_id>0) declared.emplace(m.nodes[i].typed_id,i);
    for (const auto& cut:F.cuts) {
      if(cut[2]>=2)continue;
      const auto& family=F.families[cut[1]];
      auto d=declared.find(family.dam),s=declared.find(family.sire);
      if(d==declared.end()||s==declared.end())continue;
      active[d->second]=active[s->second]=1;
      paired_families.insert(cut[1]);
    }
    plan.restored_family_pairs=static_cast<int>(paired_families.size());
  }
  std::vector<int> var(nodes, -1);
  int count = 0;
  // Canonical names: 1e9 + typed id (model individual), latent model nodes
  // below 1e9, 2e9 + id (forest individual), 3e9 + edge (forest parent).
  constexpr long long Model = 1000000000LL, Steiner = 2000000000LL, Fresh = 3000000000LL;
  std::vector<long long> name;
  for (int i = 0; i < nodes; ++i) if (active[i]) { var[i] = count++; name.push_back(m.nodes[i].typed_id > 0 ? Model + m.nodes[i].typed_id : i); }
  // Local individuals are kept sparse: a question never sweeps the pedigree.
  std::unordered_map<int, int> local_of;
  std::unordered_set<int> modeled_birth;
  std::vector<int> model_ids;
  for (int i = 0; i < nodes; ++i) {
    if (!active[i]) continue;
    const int id = m.nodes[i].typed_id;
    if (id <= 0) continue;
    if (id > F.n) throw std::out_of_range("M_ped2 typed node outside the pedigree");
    if (!local_of.emplace(id, var[i]).second) throw std::logic_error("M_ped2 typed id owned by two model nodes");
    model_ids.push_back(id);
    if (birth[i]) modeled_birth.insert(id);
  }
  auto local = [&](int id) { const auto at = local_of.find(id); return at == local_of.end() ? -1 : at->second; };
  // Steiner subtree of the model's typed individuals in each context component:
  // the union of forest paths between neighbours in Euler order.
  std::map<int, std::vector<int>> groups;
  for (int id : model_ids) groups[ctx.component(Forest::ind(id))].push_back(Forest::ind(id));
  std::unordered_set<int> region;
  for (auto& [component, members] : groups) {
    std::sort(members.begin(), members.end(), [&](int a, int b) { return F.tin[a] < F.tin[b]; });
    region.insert(members.begin(), members.end());
    for (std::size_t i = 1; i < members.size(); ++i) for (int x : F.path(members[i - 1], members[i])) region.insert(x);
  }
  // A re-parented typed individual leaves its forest family; the rest of that
  // family is still evidence, so the family is evaluated locally without it.
  for (int id : model_ids) if (modeled_birth.count(id) && ctx.child_edge(id) >= 0) region.insert(Forest::fam(F.edges[ctx.child_edge(id)].family));
  for (int f:paired_families) region.insert(Forest::fam(f));
  std::vector<int> individuals, families;
  for (int x : region) (x % 2 == 0 ? individuals : families).push_back(x / 2);
  std::sort(individuals.begin(), individuals.end());
  std::sort(families.begin(), families.end());
  const int model_count = count;
  for (int id : individuals) if (local(id) < 0) { local_of.emplace(id, count++); name.push_back(Steiner + id); }
  plan.steiner_individuals = count - model_count;
  std::set<int> touched;
  for (int x : region) {
    const int c = ctx.component(x);
    if (c < F.components) touched.insert(c);
  }
  plan.touched_components.assign(touched.begin(), touched.end());
  if (ctx.removed >= 0) {
    bool side = false;
    for (int x : region) if (ctx.family_side(x)) { side = true; break; }
    plan.family_side_untouched = !side;
  }

  for (int i = 0; i < nodes; ++i) {
    if (!active[i]) continue;
    const auto& node = m.nodes[i];
    if (node.typed_id > 0 && node.observed && !node.neutral_phi && observation(node.typed_id)) plan.unaries.push_back({var[i], Plan::Emission, node.typed_id});
    if (birth[i]) {
      const int d = node.parents[0] >= 0 ? var[node.parents[0]] : -1, s = node.parents[1] >= 0 ? var[node.parents[1]] : -1;
      if (d < 0 && s < 0) plan.unaries.push_back({var[i], Plan::Prior, 0});
      else plan.births.push_back({{{var[i], d, s}}});
    } else if (node.typed_id <= 0) {
      plan.unaries.push_back({var[i], Plan::Prior, 0});
    }
  }
  auto local_family = [&](int f) { return region.count(Forest::fam(f)) != 0; };
  std::size_t model_unaries = plan.unaries.size(), model_births = plan.births.size();
  for (int id : individuals) {
    const int v = local(id);
    if (v >= model_count && observation(id)) plan.unaries.push_back({v, Plan::Emission, id});
    if (!modeled_birth.count(id) && ctx.child_edge(id) < 0) plan.unaries.push_back({v, Plan::Prior, 0});
    for (int e : F.individual_edges[id]) {
      if (e == ctx.removed || local_family(F.edges[e].family)) continue;
      if (F.edges[e].role == 2 && modeled_birth.count(id)) throw std::logic_error("M_ped2 re-parented family was not localized");
      if(ancestral&&F.edges[e].role<2)continue; // certified upward unit message
      if(ancestral&&!allowed.count(id)){plan.unaries.push_back({v,Plan::Prior,0});++plan.off_scope_messages;}
      else plan.unaries.push_back({v, Plan::Incoming, e});
    }
  }
  for (int f : families) {
    const bool exclude_brood=latent_brood_exclusion &&
      (F.families[f].dam==0 || F.families[f].sire==0);
    std::array<int, 2> parent{{-1, -1}};
    for (int s = 0; s < 2; ++s) {
      const int e = F.family_parent_edge[f][s];
      if (e < 0) {
        if(paired_families.count(f)) {
          parent[s]=local(s==0?F.families[f].dam:F.families[f].sire);
          if(parent[s]<0)throw std::logic_error("M_ped2 restored family lost its declared parent");
        }
        continue;
      }
      int v = local(F.edges[e].individual);
      if (v < 0) { v = count++; name.push_back(Fresh + e);
        if(ancestral&&!allowed.count(F.edges[e].individual)){plan.unaries.push_back({v,Plan::Prior,0});++plan.off_scope_messages;}
        else plan.unaries.push_back({v, Plan::Outgoing, e}); }
      parent[s] = v;
    }
    std::vector<int> outside;
    for (int e : F.family_child_edges[f]) {
      if (e == ctx.removed) continue;
      const int c = F.edges[e].individual;
      if (modeled_birth.count(c)) continue;
      const int v = local(c);
      if (v >= 0) {
        // Literal question births retain their inheritance factors. Only the
        // exterior source-family contribution is excluded. A localized
        // exterior child keeps its own evidence with a population prior.
        if (parent[0] < 0 && parent[1] < 0) plan.unaries.push_back({v, Plan::Prior, 0});
        else {
          const bool excluded=exclude_brood&&!declared_births.count(c);
          plan.births.push_back({{{v,parent[0],parent[1]}},{{-1,-1}},excluded});
          if(excluded)plan.unaries.push_back({v,Plan::Prior,0});
        }
      } else outside.push_back(e);
    }
    // Exterior child->family vectors are 1 in descendant-free messages.
    // Eliminate the resulting constant factor before planning; preserve the
    // supplied population-gamete normalizer for a missing parent. Custom
    // transmissions without the existing Mendelian certificate stay generic.
    if (!outside.empty() && (parent[0] >= 0 || parent[1] >= 0)) {
      if(ancestral&&limits.ancestral_unit_pairs&&reversible_keys())
        plan.barren_population_normalizers+=int(outside.size())*((parent[0]<0)+(parent[1]<0));
      else plan.pairs.push_back({cache.slot(f,parent[0]>=0,parent[1]>=0,outside),parent});
    }
  }
  plan.family_terms = static_cast<int>(plan.births.size() + plan.pairs.size());
  if(ancestral)for(const auto& u:plan.unaries){
    if(u.source==Plan::Emission&&!allowed.count(u.index))throw std::logic_error("off-scope ancestral observation");
    if((u.source==Plan::Incoming||u.source==Plan::Outgoing)&&!allowed.count(F.edges[u.index].individual))throw std::logic_error("off-scope ancestral message");
  }
  int focal_var=retain_focal&&m.focal>=0?var[m.focal]:-1;
  prune_barren(plan,name,count,model_unaries,model_births,focal_var);
  return finish_plan(std::move(plan),name,count,model_unaries,model_births,pieces,
    split,max_touch,focal_var);
}

std::shared_ptr<const Program> compile_scalar_program(const Plan& plan, int G, int k) {
  auto out = std::make_shared<Program>();
  Program& P = *out;
  plan.limits.validate();P.limits=plan.limits;P.G = G; P.k = k;P.model_variables=plan.variables;P.model_births=int(plan.births.size());P.model_pieces=int(plan.pieces.size());P.retained_axes=plan.retained;
  int single_births=0;
  if(plan.factorized_one)for(const auto& b:plan.births)single_births+=b.gametes[0]<0;
  const int V = plan.variables+single_births;
  int next_gamete=plan.variables;
  const std::size_t npos = std::numeric_limits<std::size_t>::max();
  std::vector<int> dom(V);
  for (int v = 0; v < V; ++v) dom[v] = v>=plan.variables || plan.allele[v] ? k : G;
  std::vector<std::vector<int>> scope;   // per slot
  std::vector<std::vector<int>> excluded_scopes;
  std::vector<int> excluded_axes;
  P.unary_offset.assign(V, npos);
  for (const auto& u : plan.unaries)
    if (P.unary_offset[u.var] == npos) {
      P.unary_offset[u.var] = P.arena; P.arena += G;
      P.slots.push_back({Program::UnaryTable, u.var, P.unary_offset[u.var]}); scope.push_back({u.var});
    }
  for (const auto& b : plan.births) {
    if(b.excluded_brood) {
      // Keep the old ordering envelope and partition, while evaluating only
      // the child's added prior. Deleted couplings must not make a heuristic
      // choose a wider contraction. No numerical table is allocated here.
      const int parent=b.vars[1]>=0?b.vars[1]:b.vars[2];
      if(plan.factorized_one){const int a=next_gamete++;
        excluded_scopes.push_back({parent,a});excluded_scopes.push_back({b.vars[0],a});excluded_axes.push_back(a);}
      else excluded_scopes.push_back({b.vars[0],parent});
      continue;
    }
    if (b.gametes[0] >= 0) {
      P.slots.push_back({Program::GameteTable, 0, 0}); scope.push_back({b.vars[1], b.gametes[0]});
      P.slots.push_back({Program::GameteTable, 0, 0}); scope.push_back({b.vars[2], b.gametes[1]});
      P.slots.push_back({Program::PairingTable, 0, 0}); scope.push_back({b.vars[0], b.gametes[0], b.gametes[1]});
    } else if(plan.factorized_one) {
      // Sum the known parent's transmitted allele, then pair it with a
      // population gamete. Both factors have G*k domains, never G*G.
      const int a=next_gamete++;
      P.slots.push_back({Program::GameteTable,0,0});scope.push_back({b.vars[1]>=0?b.vars[1]:b.vars[2],a});
      P.slots.push_back({Program::OneGameteTable,0,0});scope.push_back({b.vars[0],a});
    } else {
      P.slots.push_back({Program::OneTable, 0, 0}); scope.push_back({b.vars[0], b.vars[1] >= 0 ? b.vars[1] : b.vars[2]});
    }
  }
  for (std::size_t i = 0; i < plan.pairs.size(); ++i) {
    std::vector<int> sc;
    for (int v : plan.pairs[i].vars) if (v >= 0) sc.push_back(v);
    P.slots.push_back({Program::PairTable, static_cast<int>(i), 0}); scope.push_back(sc);
  }
  for (std::size_t i = 0; i < plan.pieces.size(); ++i) {
    P.slots.push_back({Program::PieceTable, static_cast<int>(i), 0}); scope.push_back(plan.pieces[i].vars);
  }
  std::vector<unsigned char> covered(V, 0);
  for (const auto& sc : scope) for (int v : sc) covered[v] = 1;
  for (int v = plan.retained; v < V; ++v) if (!covered[v]&&std::find(excluded_axes.begin(),excluded_axes.end(),v)==excluded_axes.end()) P.log_uncovered += std::log(static_cast<double>(dom[v]));
  const std::size_t original_scope_count=scope.size();
  auto finish=[&]() -> std::shared_ptr<const Program> {
    // Preserve original factors only when a refused program needs a bounded
    // exact retry. Already admitted plans retain their existing metadata cost.
    if(!P.ready&&P.audit_scopes.empty()) {
      P.audit_scopes.assign(scope.begin(),scope.begin()+original_scope_count);P.audit_domains=dom;
    }
    return out;
  };
  if(!resource_audit().directory.empty()||plan.limits.resource_planning){P.audit_scopes=scope;P.audit_domains=dom;}
  const int R = plan.retained;
  auto resource_step = [&](const Program::Step& step, const std::vector<int>& axes, double cells, int select) {
    P.resource_step=int(P.steps.size());P.eliminated_dim=step.dim;P.input_count=int(step.inputs.size());
    P.output_dims=step.dims;P.required_cells=cells;P.required_arena=double(P.arena)+cells;
    P.genotype_axes=P.allele_axes=0;P.selection=select;
    for(int v:axes)if(dom[v]==G)++P.genotype_axes;else ++P.allele_axes;
    P.step_work=select>=0?cells*(P.input_count-1):cells*step.dim*P.input_count;
    P.required_operations=P.operations+P.step_work;
  };
  auto memory_refused = [&]() {
    // Domain count is not a resource bound: 55^3 is only 1.27 MiB.
    // Admit any exact shape that fits the configured table, arena and work
    // limits. Large G^3 factors (including k30/k40) still exceed these limits.
    if(!std::isfinite(P.required_cells)||P.required_cells>plan.limits.cells()){P.refusal=plan.limits.table_reason();return true;}
    return false;
  };
  std::vector<unsigned char> pairing(scope.size(), 0);
  for (std::size_t i = 0; i < P.slots.size(); ++i) pairing[i] = P.slots[i].kind == Program::PairingTable;
  std::vector<std::vector<int>> ordering_storage;
  std::vector<unsigned char> ordering_marks;
  if(!excluded_scopes.empty()){
    ordering_storage=scope;ordering_storage.insert(ordering_storage.end(),excluded_scopes.begin(),excluded_scopes.end());
    ordering_marks=pairing;ordering_marks.resize(ordering_storage.size(),0);
  }
  const auto& ordering_scope=excluded_scopes.empty()?scope:ordering_storage;
  const auto& ordering_pairing=excluded_scopes.empty()?pairing:ordering_marks;
  if (!plan.limits.domain_order && !single_births && plan.order.empty() && V > R) {
    plan.order_width = 0; plan.order = choose_order(ordering_scope, ordering_pairing, dom, R, &plan.order_width);
#ifdef PAP_LIKELIHOOD_DUMP
    if (const char* path = std::getenv("MPED2_DUMP_SCOPES")) {
      if (FILE* f = std::fopen(path, "a")) {
        std::fprintf(f, "plan %d %d %d", V, R, static_cast<int>(scope.size()));
        for (int v = 0; v < V; ++v) std::fprintf(f, " %d", plan.allele[v] ? 1 : 0);
        for (std::size_t i = 0; i < scope.size(); ++i) {
          std::fprintf(f, " | %d", pairing[i] ? 1 : 0);
          for (int v : scope[i]) std::fprintf(f, " %d", v);
        }
        std::fprintf(f, "\n"); std::fclose(f);
      }
    }
#endif
  }
  const bool domain_order=plan.limits.domain_order || single_births || (plan.high_allele_programs && k>9);
  int factorized_width=0;
  const std::vector<int> factorized_order=domain_order?choose_order(ordering_scope,ordering_pairing,dom,R,&factorized_width):std::vector<int>{};
  const std::vector<int>& order = domain_order?factorized_order:plan.order;
  P.width = std::max(P.width, domain_order?factorized_width:plan.order_width);
  const int E = static_cast<int>(order.size());
  std::vector<int> position(V, E);
  for (int i = 0; i < E; ++i) position[order[i]] = i;
  std::vector<std::vector<int>> bucket(E);
  std::vector<int> kept;   // tables over retained variables only
  auto place = [&](int slot) {
    if (scope[slot].empty()) { P.scalars.push_back(slot); return; }
    int first = E;
    for (int v : scope[slot]) first = std::min(first, position[v]);
    if (first == E) kept.push_back(slot); else bucket[first].push_back(slot);
  };
  for (int i = 0; i < static_cast<int>(P.slots.size()); ++i) place(i);
  // A refused sum can repeat identical-scope messages inside every eliminated
  // state. Multiply those messages once first, using the existing dim=1 product
  // operation. This retains every factor and every axis, and only changes a
  // program that the original construction would refuse. Pairing selectors are
  // left intact so their one-compatible-term execution remains available.
  auto coalesce_inputs = [&](std::vector<int>& inputs) {
    std::map<std::vector<int>,std::vector<int>> groups;
    for(int slot:inputs)if(P.slots[slot].kind!=Program::PairingTable){
      auto axes=scope[slot];std::sort(axes.begin(),axes.end());groups[axes].push_back(slot);
    }
    std::vector<Program::Step> products;std::vector<std::vector<int>> axes;
    double extra=0;
    for(const auto& [key,members]:groups)for(std::size_t begin=0;begin+1<members.size();begin+=10){
      const auto end=std::min(begin+10,members.size());
      auto joint=key;std::sort(joint.begin(),joint.end(),[&](int a,int b){return position[a]<position[b];});
      Program::Step step;step.dim=1;step.variable=-1;step.inputs.assign(members.begin()+begin,members.begin()+end);
      double cells=1;for(int v:joint){step.dims.push_back(dom[v]);cells*=dom[v];}
      if(!std::isfinite(cells)||cells>plan.limits.cells()||cells*step.inputs.size()>plan.limits.work())return;
      extra+=cells;if(extra+double(P.arena)>plan.limits.cells())return;
      step.cells=std::size_t(cells);step.axis0.assign(step.inputs.size(),0);step.strides.assign(step.inputs.size()*joint.size(),0);
      for(std::size_t i=0;i<step.inputs.size();++i){std::size_t stride=1;
        for(int v:scope[step.inputs[i]]){auto j=std::find(joint.begin(),joint.end(),v)-joint.begin();step.strides[i*joint.size()+j]=stride;stride*=dom[v];}}
      products.push_back(std::move(step));axes.push_back(std::move(joint));
    }
    if(products.empty())return;
    std::map<int,int> replacement;std::set<int> absorbed;
    for(std::size_t j=0;j<products.size();++j){auto& step=products[j];
      step.offset=P.arena;P.arena+=step.cells;step.out=int(P.slots.size());
      replacement[step.inputs.front()]=step.out;
      for(std::size_t i=1;i<step.inputs.size();++i)absorbed.insert(step.inputs[i]);
      P.operations+=double(step.cells)*step.inputs.size();
      P.slots.push_back({Program::Produced,0,step.offset});scope.push_back(std::move(axes[j]));P.steps.push_back(std::move(step));
    }
    std::vector<int> next;for(int slot:inputs)if(!absorbed.count(slot)){
      auto at=replacement.find(slot);next.push_back(at==replacement.end()?slot:at->second);}
    inputs=std::move(next);P.coalesced_inputs=true;
  };
  for (int at = 0; at < E; ++at) {
    if (bucket[at].empty()) continue;
    const int best = order[at];
    std::vector<int> joint;
    for (int slot : bucket[at]) for (int v : scope[slot]) if (v != best && std::find(joint.begin(), joint.end(), v) == joint.end()) joint.push_back(v);
    // Produced tables are laid out in elimination order, so the next
    // variable summed out of them is their contiguous axis.
    std::sort(joint.begin(), joint.end(), [&](int a, int b) { return position[a] != position[b] ? position[a] < position[b] : a < b; });
    Program::Step step;step.variable=best;
    step.inputs = bucket[at];
    step.dim = dom[best];
    const std::size_t m = joint.size(), n = step.inputs.size();
    double cells = 1.0;
    for (int v : joint) { step.dims.push_back(dom[v]); cells *= dom[v]; }
    for(std::size_t i=0;i<n;++i)if(P.slots[step.inputs[i]].kind==Program::PairingTable&&n>=2&&n<=10){
      step.select=int(i);const auto& sc=scope[step.inputs[i]];
      step.select_role=int(std::find(sc.begin(),sc.end(),best)-sc.begin());
    }
    resource_step(step,joint,cells,step.select);
    const bool modern=plan.high_allele_programs&&k>9;
    const double work_used=plan.limits.total_work_limit?P.required_operations:
      modern?P.step_work:cells*step.dim*n;
    const double work_cap=(modern||plan.limits.total_work_limit)?plan.limits.work():4e8;
    if(work_used>work_cap){
      const auto before=P.steps.size();if(plan.limits.s37_refusal_repair)coalesce_inputs(bucket[at]);
      if(P.steps.size()!=before){--at;continue;}
      P.refusal=(modern||plan.limits.total_work_limit)?plan.limits.work_reason():"elimination work exceeds 400M operations";return finish();
    }
    if(memory_refused())return finish();
    step.cells = static_cast<std::size_t>(cells);
    step.axis0.assign(n, 0); step.strides.assign(n * m, 0);
    for (std::size_t i = 0; i < n; ++i) {
      std::size_t stride = 1;
      for (int v : scope[step.inputs[i]]) {
        if (v == best) step.axis0[i] = stride;
        else step.strides[i * m + (std::find(joint.begin(), joint.end(), v) - joint.begin())] = stride;
        stride *= dom[v];
      }
    }
    const double work=step.select >= 0 ? cells * (n - 1) : cells * step.dim * n;
    // Pairing reductions select one Mendelian-compatible term. Charge the
    // compiled work, rather than a dense sum that will never be executed.
    resource_step(step,joint,cells,step.select);
    if(plan.high_allele_programs && k>9 && work>plan.limits.work()){P.refusal=plan.limits.work_reason();return finish();}
    P.operations += work;
    if(!plan.limits.live_arena && (step.cells>plan.limits.cells() || P.arena>plan.limits.cells()-step.cells)){P.refusal=plan.limits.arena_reason();return finish();}
    if(step.cells>std::numeric_limits<std::size_t>::max()-P.arena)throw std::overflow_error("program arena offsets overflow");
    step.offset = P.arena; P.arena += step.cells;
    step.out = static_cast<int>(P.slots.size());
    P.slots.push_back({Program::Produced, 0, step.offset});
    scope.push_back(joint);
    P.steps.push_back(std::move(step));
    place(P.steps.back().out);
  }
  if (R > 0) {
    // Product of the remaining tables over the retained variables, in order.
    Program::Step step;
    step.inputs = kept;
    step.dim = 1;
    const std::size_t m = R, n = kept.size();
    double cells = 1.0;
    for (int v = 0; v < R; ++v) { step.dims.push_back(dom[v]); cells *= dom[v]; }
    std::vector<int> retained_axes(R);std::iota(retained_axes.begin(),retained_axes.end(),0);
    resource_step(step,retained_axes,cells,-1);
    if(memory_refused())return finish();
    if((plan.limits.total_work_limit?P.required_operations:P.step_work)>plan.limits.work()){P.refusal=plan.limits.work_reason();return finish();}
    step.cells = static_cast<std::size_t>(cells);
    P.operations += cells * n;
    step.axis0.assign(n, 0); step.strides.assign(n * m, 0);
    for (std::size_t i = 0; i < n; ++i) {
      std::size_t stride = 1;
      for (int v : scope[kept[i]]) { step.strides[i * m + v] = stride; stride *= dom[v]; }
    }
    if(!plan.limits.live_arena && (step.cells>plan.limits.cells() || P.arena>plan.limits.cells()-step.cells)){P.refusal=plan.limits.arena_reason();return finish();}
    if(step.cells>std::numeric_limits<std::size_t>::max()-P.arena)throw std::overflow_error("program arena offsets overflow");
    step.offset = P.arena; P.arena += step.cells;
    step.out = static_cast<int>(P.slots.size());
    P.slots.push_back({Program::Produced, 0, step.offset});
    P.result = step.out;
    P.steps.push_back(std::move(step));
  }
  if(plan.limits.live_arena) {
    P.arena_stats=live_arena::pack_program(P,plan.limits.cells());
    if(!P.arena_stats.within_limit) {
      P.refusal=plan.limits.arena_reason();
      const int at=P.arena_stats.first_exceeding_step;
      P.resource_step=at;
      if(at>=0) {
        const auto& step=P.steps[at];P.eliminated_dim=step.dim;P.input_count=int(step.inputs.size());
        P.output_dims=step.dims;P.required_cells=step.cells;P.selection=step.select;
        P.genotype_axes=P.allele_axes=0;
        for(int d:step.dims)if(d==G)++P.genotype_axes;else ++P.allele_axes;
        P.step_work=step.select>=0?double(step.cells)*(P.input_count-1):double(step.cells)*step.dim*P.input_count;
        P.required_operations=0;
        for(int i=0;i<=at;++i){const auto& q=P.steps[i];P.required_operations+=q.select>=0?double(q.cells)*(q.inputs.size()-1):double(q.cells)*q.dim*q.inputs.size();}
        // The high-water extent through this allocation, not monotonic sum.
        P.required_arena=0;
        for(int i=0;i<=at;++i)P.required_arena=std::max(P.required_arena,double(P.steps[i].offset)+P.steps[i].cells);
        for(auto offset:P.unary_offset)if(offset!=npos)P.required_arena=std::max(P.required_arena,double(offset)+G);
      } else P.required_arena=P.arena;
      return finish();
    }
  }
  if(plan.limits.total_work_limit&&P.operations>plan.limits.work()){P.refusal=plan.limits.work_reason();return finish();}
  P.ready = true;P.profile_identity=program_profile_identity(P);
  return finish();
}

std::shared_ptr<const Program> compile_program(const Plan& plan,int G,int k){
  auto scalar=compile_scalar_program(plan,G,k);
  if(scalar->ready||plan.retained)return scalar;
  std::shared_ptr<resource_plan::Execution> selected;resource_plan::SearchBudget search_effort;
  if(plan.limits.resource_planning&&k>=plan.limits.resource_min_alleles) {
    const auto started=std::chrono::steady_clock::now();
    selected=resource_plan::plan(*scalar,scalar->audit_scopes,scalar->audit_domains,&search_effort);
    if(program_cache_census_enabled()){program_cache_census()[9]+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();++program_cache_census()[10];}
  } else if(plan.limits.s37_refusal_repair) {
    // A previously refused scalar question may use the established exact
    // backend at any allele count when its actual work/storage fit. Bound
    // automatic preparation by 32,768 states and a 64-state beam, and do not
    // launch a per-genotype conditioned search from the default path.
    std::vector<resource_plan::Factor> factors;
    for(std::size_t j=0;j<scalar->audit_scopes.size();++j)
      factors.push_back({scalar->slots[j].kind,int(j),scalar->audit_scopes[j]});
    selected=resource_plan::search(factors,scalar->audit_domains,G,k,plan.limits,-1,32768,64);
  }
  if(!selected){
    if(search_effort.capped){auto refused=std::make_shared<Program>(*scalar);
      refused->search_capped=true;refused->search_attempts=search_effort.attempts;
      refused->refusal="exact resource-plan search effort cap reached";return refused;}
    return scalar;
  }
  auto result=std::make_shared<Program>(*scalar);auto& p=*result;
  p.slots.resize(p.audit_scopes.size());p.steps.clear();p.scalars.clear();
  p.arena=0;for(auto& slot:p.slots)if(slot.kind==Program::UnaryTable){slot.offset=p.arena;p.unary_offset[slot.index]=p.arena;p.arena+=G;}
  if(p.arena>plan.limits.cells()||selected->arena>plan.limits.cells()-p.arena)return scalar;
  p.search_capped=search_effort.capped;p.search_attempts=search_effort.attempts;
  p.arena+=selected->arena;p.operations=selected->work;p.exact=std::move(selected);p.ready=true;p.refusal=nullptr;
  p.resource_step=-1;p.output_dims.clear();p.required_cells=p.required_arena=p.step_work=p.required_operations=0;
  p.profile_identity=program_profile_identity(p);p.profile_identity=mix_hash(p.profile_identity,0x52504c414e303031ULL);
  p.profile_identity=mix_hash(p.profile_identity,p.exact->condition);
  for(const auto& f:p.exact->factors){p.profile_identity=mix_hash(p.profile_identity,f.kind);for(int v:f.scope)p.profile_identity=mix_hash(p.profile_identity,v);p.profile_identity=mix_hash(p.profile_identity,-1);}
  for(const auto& s:p.exact->steps){for(int v:{s.op,s.variable,s.auxiliary,s.parent,s.select})p.profile_identity=mix_hash(p.profile_identity,v);for(int id:s.inputs)p.profile_identity=mix_hash(p.profile_identity,id);p.profile_identity=mix_hash(p.profile_identity,-1);}
  return result;
}

std::vector<int> plan_shape(const Plan& plan, int G, int k) {
  std::vector<int> key{G, k, plan.limits.mib,plan.limits.work_millions,int(plan.limits.live_arena),int(plan.limits.shared_piece_scratch),int(plan.factorized_one), int(plan.high_allele_programs), plan.variables, plan.retained, static_cast<int>(plan.unaries.size()),
                       static_cast<int>(plan.births.size()), static_cast<int>(plan.pairs.size()), static_cast<int>(plan.pieces.size())};
  for (unsigned char a : plan.allele) key.push_back(a);
  for (const auto& u : plan.unaries) key.push_back(u.var);
  for (const auto& b : plan.births) { key.insert(key.end(), b.vars.begin(), b.vars.end()); key.insert(key.end(), b.gametes.begin(), b.gametes.end()); }
  for (const auto& x : plan.pairs) key.insert(key.end(), x.vars.begin(), x.vars.end());
  for (const auto& x : plan.pieces) { key.push_back(static_cast<int>(x.vars.size())); key.insert(key.end(), x.vars.begin(), x.vars.end()); }
  // Preserve old keys when OFF; complete numerical-admission identity when ON.
  if(plan.limits.total_work_limit)key.insert(key.end(),{-70015,1});
  if(plan.limits.domain_order)key.insert(key.end(),{-70012,1});
  if(plan.limits.contraction_log_retry)key.insert(key.end(),{-70013,1});
  if(plan.limits.s37_refusal_repair)key.insert(key.end(),{-70010,1});
  if(plan.limits.resource_planning){key.push_back(-70001);key.push_back(1);}
  if(plan.limits.resource_min_alleles!=20){key.push_back(-70003);key.push_back(plan.limits.resource_min_alleles);}
  for(std::size_t n=0;n<plan.births.size();++n)if(plan.births[n].excluded_brood)key.insert(key.end(),{-70008,int(n)});
  return key;
}
auto PairCache::slot(int family, bool dam, bool sire, const std::vector<int>& edges) -> int  {
    std::vector<int> key{family, dam ? 1 : 0, sire ? 1 : 0};
    key.insert(key.end(), edges.begin(), edges.end());
    auto at = index.emplace(std::move(key), static_cast<int>(slots.size()));
    if (at.second) { Slot s; s.family = family; s.dam = dam; s.sire = sire; s.edges = edges; slots.push_back(std::move(s)); }
    return at.first->second;
  }

auto PieceCache::slot(std::vector<long long> key, Plan&& plan) -> int  {
    auto at = index.emplace(std::move(key), static_cast<int>(slots.size()));
    if (at.second) { auto x = std::make_unique<Slot>(); x->plan = std::move(plan); slots.push_back(std::move(x)); }
    return at.first->second;
  }

auto ShapeHash::operator()(const std::vector<int>& v) const -> std::size_t  {
    std::uint64_t h = 1469598103934665603ULL;
    for (int x : v) h = mix_hash(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(x)));
    return static_cast<std::size_t>(h);
  }


bool& reversible_keys() { thread_local bool enabled = false; return enabled; }

bool& program_registry_enabled() { thread_local bool enabled = true; return enabled; }

std::array<double,12>& program_cache_census(){thread_local std::array<double,12> x{};return x;}

bool& program_cache_census_enabled(){thread_local bool on=false;return on;}

std::string& program_cache_trace(){thread_local std::string path;return path;}

const Program& program_for(const Plan& plan, int G, int k) {
  using Clock=std::chrono::steady_clock;
  const bool audit=program_cache_census_enabled();const auto started=audit?Clock::now():Clock::time_point{};
  auto spent=[&](){return std::chrono::duration<double>(Clock::now()-started).count();};
  if(audit)++program_cache_census()[0];
  if (prepared_programs_only) {
    if (static_cast<int>(plan.programs.size()) <= k || !plan.programs[k] || plan.programs[k]->G != G || !(plan.programs[k]->limits==plan.limits))
      throw std::logic_error("locus job attempted unprepared program compilation");
    return *plan.programs[k];
  }
  if (static_cast<int>(plan.programs.size()) <= k) plan.programs.resize(k + 1);
  auto& p = plan.programs[k];
  if (p && p->G == G && p->limits==plan.limits) {if(audit){++program_cache_census()[1];program_cache_census()[6]+=spent();}return *p;}
  auto& registry = program_registry();
  auto key = program_registry_enabled()?plan_shape(plan,G,k):std::vector<int>{};
  if(program_registry_enabled()){
    auto hit = registry.lookup(key);
    if (hit) { p = std::move(hit);if(audit){++program_cache_census()[2];program_cache_census()[7]+=spent();}return *p; }
  }
  const double search_before=program_cache_census()[9];
  p = compile_program(plan, G, k);
  if(audit){++program_cache_census()[3];program_cache_census()[4]+=!p->ready;const double seconds=spent();program_cache_census()[8]+=seconds;
    if(!program_cache_trace().empty()){std::ofstream out(program_cache_trace(),std::ios::app);out<<G<<'\t'<<k<<'\t'<<plan.variables<<'\t'<<plan.retained<<'\t'<<p->ready<<'\t'<<std::setprecision(17)<<seconds<<'\t'<<program_cache_census()[9]-search_before<<'\t';for(int x:key)out<<x<<',';out<<'\n';}}
  if(!program_registry_enabled())return *p;
  const auto clears=registry.clearings;
  registry.insert(std::move(key), p);
  if(audit)program_cache_census()[5]+=registry.clearings-clears;
  return *p;
}
}
