#include "visit_graph.h"
namespace pap_internal::search {

VisitGraph::VisitGraph(EngineState& compiled_, int focal_,
           int graph_revision_, const IntegerMatrix& parents,
           EngineCounters* telemetry_override )
      : compiled(compiled_),
        telemetry(telemetry_override == nullptr ? compiled_.telemetry :
                  *telemetry_override),
        focal(focal_),
        dam(static_cast<std::size_t>(compiled_.n + 1), 0),
        sire(static_cast<std::size_t>(compiled_.n + 1), 0),
        reverse(static_cast<std::size_t>(compiled_.n + 1)),
        node_to_evidence_row(
          static_cast<std::size_t>(compiled_.n + 1), kNoEvidenceRow) {
    if (parents.nrow() != compiled.n || parents.ncol() != 2)
      engine_stop("simp batch parents must be an n by 2 integer matrix");
    initialize_parents(parents.begin(), false);
  }

VisitGraph::VisitGraph(EngineState& compiled_, int focal_,
           int graph_revision_, const std::vector<int>& parents_column_major,
           EngineCounters& telemetry_override )
      : compiled(compiled_),
        telemetry(telemetry_override), focal(focal_), graph_revision(graph_revision_),
        dam(static_cast<std::size_t>(compiled_.n + 1), 0),
        sire(static_cast<std::size_t>(compiled_.n + 1), 0),
        reverse(static_cast<std::size_t>(compiled_.n + 1)),
        node_to_evidence_row(
          static_cast<std::size_t>(compiled_.n + 1), kNoEvidenceRow) {
    const std::size_t expected = static_cast<std::size_t>(compiled_.n) * 2U;
    if (parents_column_major.size() != expected)
      engine_stop("native scheduled worker parent matrix has incompatible dimensions");
    initialize_parents(parents_column_major.data(), true);
  }

void VisitGraph::initialize_parents(const int* parents, bool scheduled_worker) {
  for (int child = 1; child <= compiled.n; ++child) {
    const int d = parents[child - 1];
    const int s = parents[compiled.n + child - 1];
    if (d < 0 || d > compiled.n || s < 0 || s > compiled.n)
      engine_stop(scheduled_worker ?
        "native scheduled worker parent IDs are out of range" :
        "simp batch parent IDs are out of range");
    if (d == child || s == child || (d > 0 && d == s))
      engine_stop(scheduled_worker ?
        "native scheduled worker parent matrix has an invalid edge" :
        "simp batch parents contain a self, duplicate, or invalid edge");
    dam[child] = d;
    sire[child] = s;
    node_to_evidence_row[static_cast<std::size_t>(child)] = child - 1;
    if (d > 0) reverse[d].push_back(child);
    if (s > 0) reverse[s].push_back(child);
  }
}

auto VisitGraph::parent(int child, int side) const -> int  {
    return side == 1 ? dam[child] : sire[child];
  }

auto VisitGraph::set_focal(int child, int revision) -> void  {
    focal = child;
    graph_revision = revision;
    ibd_cache.clear();
  }

auto VisitGraph::remove_reverse_edge(int parent_id, int child) -> void  {
    if (parent_id <= 0) return;
    std::vector<int>& children = reverse[parent_id];
    children.erase(std::remove(children.begin(), children.end(), child),
                   children.end());
  }

auto VisitGraph::add_reverse_edge(int parent_id, int child) -> void  {
    if (parent_id <= 0) return;
    std::vector<int>& children = reverse[parent_id];
    if (std::find(children.begin(), children.end(), child) == children.end()) {
      children.push_back(child);
      std::sort(children.begin(), children.end());
    }
  }

auto VisitGraph::apply_parent_move(int child, int new_dam, int new_sire) -> void  {
    remove_reverse_edge(dam[child], child);
    remove_reverse_edge(sire[child], child);
    dam[child] = new_dam;
    sire[child] = new_sire;
    add_reverse_edge(new_dam, child);
    add_reverse_edge(new_sire, child);
    ibd_cache.clear();
  }

auto VisitGraph::ancestors(int node) const -> std::vector<int>  {
    std::vector<int> distance(static_cast<std::size_t>(compiled.n + 1), -1);
    std::queue<int> queue;
    distance[node] = 0;
    queue.push(node);
    while (!queue.empty()) {
      const int current = queue.front();
      queue.pop();
      const int next_distance = distance[current] + 1;
      for (int parent_id : {dam[current], sire[current]}) {
        if (parent_id <= 0) continue;
        if (distance[parent_id] < 0 || distance[parent_id] > next_distance) {
          distance[parent_id] = next_distance;
          queue.push(parent_id);
        }
      }
    }
    return distance;
  }

auto VisitGraph::ibd(int first, int second) -> IbdClass  {

    const std::pair<int, int> key(std::min(first, second),
                                  std::max(first, second));
    const auto found = ibd_cache.find(key);
    if (found != ibd_cache.end()) {

      return found->second;
    }

    IbdClass out;
    if (first == second) {
      out.k0 = 0.0; out.k1 = 0.0; out.k2 = 1.0;
      out.relation = ROUTE_DUPLICATE;
      ibd_cache[key] = out;

      return out;
    }
    std::vector<int> p1, p2;
    if (dam[first] > 0) p1.push_back(dam[first]);
    if (sire[first] > 0) p1.push_back(sire[first]);
    if (dam[second] > 0) p2.push_back(dam[second]);
    if (sire[second] > 0) p2.push_back(sire[second]);
    int shared_count = 0;
    for (int x : p1) for (int y : p2) if (x == y) ++shared_count;
    if (p1.size() == 2 && p2.size() == 2 && shared_count == 2) {
      out.k0 = .25; out.k1 = .5; out.k2 = .25;
      out.relation = ROUTE_FULL_SIBLINGS;
      ibd_cache[key] = out;

      return out;
    }
    const std::vector<int> d1 = ancestors(first);
    const std::vector<int> d2 = ancestors(second);
    const int direct12 = d1[second];
    const int direct21 = d2[first];
    int direct = std::numeric_limits<int>::max();
    if (direct12 >= 0) direct = std::min(direct, direct12);
    if (direct21 >= 0) direct = std::min(direct, direct21);
    if (direct != std::numeric_limits<int>::max()) {
      out.k1 = std::min(1.0, std::pow(2.0, 1.0 - direct));
      out.k0 = 1.0 - out.k1;
      out.k2 = 0.0;
      out.relation = direct == 1 ? ROUTE_PARENT_OFFSPRING :
        ROUTE_DIRECT_ANCESTOR;
      ibd_cache[key] = out;

      return out;
    }
    double r = 0.0;
    for (int node = 1; node <= compiled.n; ++node) {
      if (node == first || node == second) continue;
      if (d1[node] >= 0 && d2[node] >= 0)
        r += std::pow(2.0, -static_cast<double>(d1[node] + d2[node]));
    }
    if (r == 0.0) {
      out.k0 = 1.0; out.k1 = 0.0; out.k2 = 0.0;
      out.relation = ROUTE_UNCONNECTED;
    } else {
      out.k1 = std::min(1.0, 2.0 * r);
      out.k0 = 1.0 - out.k1;
      out.k2 = 0.0;
      out.relation = shared_count == 1 ? ROUTE_HALF_SIBLINGS :
        ROUTE_SHARED_ANCESTRY;
    }
    ibd_cache[key] = out;

    return out;
  }
}
