#ifndef PAP_PEDIGREE_GRAPH_SNAPSHOT_H
#define PAP_PEDIGREE_GRAPH_SNAPSHOT_H

#include <algorithm>
#include <functional>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pap_internal {

constexpr int kNoEvidenceRow = -1;

// Immutable, non-owning graph view. Node IDs are one-based native pedigree
// IDs. Evidence rows are zero-based native observation rows; -1 means that a
// graph node is latent/untyped. The referenced vectors must outlive the view.
class PedigreeGraphSnapshot {
 public:
  PedigreeGraphSnapshot(
      int node_count, const std::vector<int>& dam,
      const std::vector<int>& sire,
      const std::vector<std::vector<int>>& children,
      const std::vector<int>& node_to_evidence_row, int graph_revision)
      : node_count_(node_count), dam_(&dam), sire_(&sire), children_(&children),
        node_to_evidence_row_(&node_to_evidence_row),
        graph_revision_(graph_revision) {
    if (node_count_ < 0 || graph_revision_ < 0 ||
        dam_->size() != static_cast<std::size_t>(node_count_ + 1) ||
        sire_->size() != static_cast<std::size_t>(node_count_ + 1) ||
        children_->size() != static_cast<std::size_t>(node_count_ + 1) ||
        node_to_evidence_row_->size() !=
          static_cast<std::size_t>(node_count_ + 1))
      throw std::invalid_argument("pedigree graph snapshot dimensions are invalid");
  }

  int node_count() const { return node_count_; }
  int graph_revision() const { return graph_revision_; }

  int parent(int node, int side) const {
    require_node(node);
    if (side != 1 && side != 2)
      throw std::out_of_range("pedigree parent side must be one or two");
    return side == 1 ? (*dam_)[static_cast<std::size_t>(node)] :
      (*sire_)[static_cast<std::size_t>(node)];
  }

  std::pair<int, int> parents(int node) const {
    require_node(node);
    return std::make_pair((*dam_)[static_cast<std::size_t>(node)],
                          (*sire_)[static_cast<std::size_t>(node)]);
  }

  const std::vector<int>& children(int node) const {
    require_node(node);
    return (*children_)[static_cast<std::size_t>(node)];
  }


  std::vector<int> topological_order() const {
    std::vector<int> indegree(static_cast<std::size_t>(node_count_ + 1), 0);
    for (int node = 1; node <= node_count_; ++node) {
      const std::pair<int, int> direct = parents(node);
      indegree[static_cast<std::size_t>(node)] =
        (direct.first > 0 ? 1 : 0) + (direct.second > 0 ? 1 : 0);
    }
    std::priority_queue<int, std::vector<int>, std::greater<int>> ready;
    for (int node = 1; node <= node_count_; ++node)
      if (indegree[static_cast<std::size_t>(node)] == 0) ready.push(node);
    std::vector<int> order;
    order.reserve(static_cast<std::size_t>(node_count_));
    while (!ready.empty()) {
      const int node = ready.top();
      ready.pop();
      order.push_back(node);
      for (int child : children(node)) {
        int& degree = indegree[static_cast<std::size_t>(child)];
        --degree;
        if (degree == 0) ready.push(child);
      }
    }
    if (order.size() != static_cast<std::size_t>(node_count_))
      throw std::invalid_argument("pedigree graph snapshot contains an ancestry cycle");
    return order;
  }

  std::vector<int> descendants(const std::vector<int>& roots) const {
    std::vector<char> seen(static_cast<std::size_t>(node_count_ + 1), 0);
    std::queue<int> pending;
    for (int root : roots) if (valid_node(root)) pending.push(root);
    while (!pending.empty()) {
      const int node = pending.front();
      pending.pop();
      if (seen[static_cast<std::size_t>(node)]) continue;
      seen[static_cast<std::size_t>(node)] = 1;
      for (int child : children(node)) pending.push(child);
    }
    return sorted_seen(seen);
  }

  std::vector<int> affected_component(const std::vector<int>& roots) const {
    std::vector<char> seen(static_cast<std::size_t>(node_count_ + 1), 0);
    std::queue<int> pending;
    for (int root : roots) if (valid_node(root)) pending.push(root);
    while (!pending.empty()) {
      const int node = pending.front();
      pending.pop();
      if (seen[static_cast<std::size_t>(node)]) continue;
      seen[static_cast<std::size_t>(node)] = 1;
      const std::pair<int, int> direct = parents(node);
      if (direct.first > 0) pending.push(direct.first);
      if (direct.second > 0) pending.push(direct.second);
      for (int child : children(node)) pending.push(child);
    }
    return sorted_seen(seen);
  }

 private:
  int node_count_;
  const std::vector<int>* dam_;
  const std::vector<int>* sire_;
  const std::vector<std::vector<int>>* children_;
  const std::vector<int>* node_to_evidence_row_;
  int graph_revision_;

  bool valid_node(int node) const { return node >= 1 && node <= node_count_; }

  void require_node(int node) const {
    if (!valid_node(node))
      throw std::out_of_range("pedigree graph node is out of range");
  }

  std::vector<int> sorted_seen(const std::vector<char>& seen) const {
    std::vector<int> out;
    for (int node = 1; node <= node_count_; ++node)
      if (seen[static_cast<std::size_t>(node)]) out.push_back(node);
    return out;
  }
};

struct DirtyRegion {
  std::vector<int> changed_nodes;
  std::vector<int> affected_nodes;
  std::vector<int> publication_nodes;
  int graph_revision_before = 0;
  int graph_revision_after = 0;

  void normalize() {
    normalize_nodes(changed_nodes);
    normalize_nodes(affected_nodes);
    normalize_nodes(publication_nodes);
  }

  bool revision_advanced_once() const {
    return graph_revision_after == graph_revision_before + 1;
  }

 private:
  static void normalize_nodes(std::vector<int>& nodes) {
    nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
                               [](int node) { return node <= 0; }),
                nodes.end());
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
  }
};

inline std::vector<int> union_nodes(std::vector<int> first,
                                    const std::vector<int>& second) {
  first.insert(first.end(), second.begin(), second.end());
  std::sort(first.begin(), first.end());
  first.erase(std::unique(first.begin(), first.end()), first.end());
  return first;
}

}  // namespace pap_internal

#endif
