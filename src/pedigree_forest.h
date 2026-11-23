#include "pap_removed_path.h"
#pragma once
// pedigree likelihood: local hypothesis models read stored link messages.
//
// The current typed pedigree is cut into a forest in typed-ID order. A parent
// incidence that would close a loop is replaced by an
// absent (population-gamete) parent. The visit context removes the focal's
// incoming birth without changing that cut choice. Every directed family message
// is computed once per locus by two-pass sum-product, with absolute log scales.
//
// A question evaluates a TopologyModel exactly over its structurally active
// typed individuals, its latent stand-ins and the forest paths (Steiner
// subtree) joining those individuals. Everything else enters through the
// stored messages. Loops inside the model are exact; loops through the rest of
// the pedigree are cut once, identically for every question of the visit.
#include "resource_limits.h"
#include "resource_audit.h"
#include "live_arena.h"
#include "relationship_topology.h"
#include "positive_logs.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <list>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace pap_internal::pedigree_likelihood {
// Worker-only receipts: no mutation of shared Forest telemetry or plan compilation.
inline thread_local std::size_t* worker_message_count = nullptr;
inline thread_local bool prepared_programs_only = false;

// Locus tables shared with the pedigree likelihood scorer (see read_loci).
struct LocusView {
  int G = 0, k = 0;
  const std::vector<double>* prior = nullptr;  // genotype prior, linear
  const std::vector<double>* p = nullptr;      // allele frequencies, linear
  const std::vector<double>* t = nullptr;      // gamete weights t[g + G a]
  const std::vector<int>* lookup = nullptr;    // child genotype of alleles a + k b
  const std::vector<std::vector<double>>* phi = nullptr;  // log emission, [id][g]
  // Optional sparse gametes: genotype g transmits allele gamete_allele[2g]
  // with weight gamete_weight[2g] and gamete_allele[2g+1] (-1 if none) with
  // gamete_weight[2g+1]; alleles ascending. Absent: dense t is used.
  const std::vector<int>* gamete_allele = nullptr;
  const std::vector<double>* gamete_weight = nullptr;
  const std::vector<int>* observed = nullptr; // observed genotype by typed ID; -1 missing
};

// Hardy-Weinberg prior consistent with the allele frequencies and Mendelian
// gametes: the condition under which single-parent links are reversible.
inline bool reversible_locus(const LocusView& L) {
  const int k = L.k, G = L.G;
  for (int a = 0; a < k; ++a)
    for (int b = a; b < k; ++b) {
      const int g = (*L.lookup)[a + k * b];
      const double expected = (*L.p)[a] * (*L.p)[b] * (a == b ? 1.0 : 2.0);
      if (std::abs((*L.prior)[g] - expected) > 1e-12 * std::max(1e-300, expected)) return false;
      for (int x = 0; x < k; ++x) {
        const double want = a == b ? (x == a ? 1.0 : 0.0) : (x == a || x == b ? 0.5 : 0.0);
        if ((*L.t)[g + G * x] != want) return false;
      }
    }
  return true;
}

// Sparse form of t when every genotype transmits at most two alleles.
inline bool sparse_gametes(const LocusView& L, std::vector<int>& allele, std::vector<double>& weight) {
  allele.assign(2 * L.G, -1); weight.assign(2 * L.G, 0.0);
  for (int g = 0; g < L.G; ++g) {
    int at = 0;
    for (int a = 0; a < L.k; ++a) {
      const double w = (*L.t)[g + L.G * a];
      if (w == 0.0) continue;
      if (at == 2) return false;
      allele[2 * g + at] = a; weight[2 * g + at] = w; ++at;
    }
  }
  return true;
}

// Bipartite forest over individuals and nuclear families. A family is one mate
// pair (or one known parent with an unknown mate) and all their children, so
// sibships are not loops. Edges are dam/sire/child incidences; an incidence
// that would close a genuine loop is cut (absent parent or founder child).
// `update` rebuilds canonical incidences from the current pedigree and
// invalidates numerical messages. Cut selection is independent of move history.
struct Forest {
  struct Family { int dam = 0, sire = 0; std::vector<int> children; };
  struct Edge { int family = -1, individual = 0, role = 0; };  // role 0 dam, 1 sire, 2 child
  int n = 0;
  bool ancestral = false;
  std::vector<int> ancestry_revision;
  std::vector<Family> families;
  std::vector<Edge> edges;
  std::vector<std::array<int, 2>> family_parent_edge;   // per family: dam/sire edge or -1
  std::vector<std::vector<int>> family_child_edges;     // per family
  std::vector<std::vector<int>> individual_edges;       // per individual
  std::vector<int> child_edge;                          // per individual, -1 founder/cut
  std::vector<int> family_of;                           // current birth family, cut or not; -1 none
  std::vector<int> component;                           // per individual
  int components = 0;
  // Traversal over nodes: 2v = individual v, 2f+1 = family f.
  std::vector<int> order, up_node, up_edge, root, depth;
  std::vector<int> tin, tout, node_component;    // Euler intervals and component of every node
  std::vector<std::array<int, 3>> cuts;   // individual, family, role of each cut incidence
  std::size_t cut_incidences = 0;
  // Persistent state for in-place updates.
  std::map<std::pair<int, int>, int> family_index;      // live (dam, sire) -> family
  // Stable family ids survive update(): one per (dam, sire) key seen since
  // the last cold build, so per-locus birth messages need no remapping after
  // a move. stable_epoch is the revision since which a stable family has
  // existed continuously (-1 while absent).
  std::map<std::pair<int, int>, int> stable_index;
  std::vector<int> family_stable;   // per current family
  std::vector<int> stable_epoch;    // per stable id
  std::vector<int> vanished;        // stable ids that disappeared at revision vanished_at
  int vanished_at = -1;
  std::vector<int> current_dam, current_sire;
  int revision = 0;
  std::vector<int> stale_to_family, stale_to_individual;   // per edge: revision its message last changed
  std::size_t updates = 0, removed_edges = 0, added_edges = 0, stale_marks = 0;
  mutable std::size_t messages_computed = 0;
  mutable std::vector<unsigned char> cached_cut_flags;
  mutable std::size_t cut_flag_builds = 0;

  static int ind(int v) { return 2 * v; }
  static int fam(int f) { return 2 * f + 1; }

  void build(const std::vector<int>& dam, const std::vector<int>& sire, int focal) {
    build(dam, sire, std::vector<int>{focal});
  }
  // Seeds are accepted for source compatibility; cuts use current IDs only.
  void build(const std::vector<int>& dam, const std::vector<int>& sire, const std::vector<int>& seeds);
  void build_families(const std::vector<int>& dam, const std::vector<int>& sire,
                      const std::vector<int>& seeds);

  // Canonical rebuild after a move. Family/edge IDs can change, so every
  // numerical message is invalidated; a warm forest equals a cold build.
  bool update(const std::vector<int>& dam, const std::vector<int>& sire);

  std::vector<unsigned char> ancestry_cut_flags() const;

  // Host-side telemetry only. build/update invalidates this with the graph.
  const std::vector<unsigned char>& ancestry_cut_flags_cached() const {
    if(cached_cut_flags.empty()){cached_cut_flags=ancestry_cut_flags();++cut_flag_builds;}
    return cached_cut_flags;
  }

  // Components, BFS order and Euler intervals of the current links.
  void reorient();
  // The rooted-tree child end of edge e, and subtree membership.
  int lower(int e) const { const int x = ind(edges[e].individual); return up_edge[x] == e ? x : fam(edges[e].family); }
  bool within(int x, int sub) const { return tin[sub] <= tin[x] && tin[x] <= tout[sub]; }
  // Forest nodes on the path between x and y (same component), inclusive.
  std::vector<int> path(int x, int y) const;
  template <class Fn> void visit_neighbors(int x, Fn&& fn) const {
    if (x % 2 == 0) { for (int e : individual_edges[x / 2]) fn(fam(edges[e].family), e); }
    else {
      int f = x / 2;
      for (int e : family_parent_edge[f]) if (e >= 0) fn(ind(edges[e].individual), e);
      for (int e : family_child_edges[f]) fn(ind(edges[e].individual), e);
    }
  }

 private:
  // Links individual v to family f in `role` unless that closes a loop.
  template <class Find> void attach(int f, int v, int role, std::vector<int>& link, Find& find) {
    int a = find(v), b = find(n + 1 + f);
    if (a == b) { ++cut_incidences; cuts.push_back({{v, f, role}}); return; }
    link[a] = b;
    int e = static_cast<int>(edges.size());
    edges.push_back(Edge{f, v, role});
    individual_edges[v].push_back(e);
    if (role < 2) family_parent_edge[f][role] = e;
    else { family_child_edges[f].push_back(e); child_edge[v] = e; }
    ++added_edges;
  }
};

// A scaled positive vector: value = cells * exp(scale).

}
