#ifndef PAP_TOPOLOGY_TYPES_H
#define PAP_TOPOLOGY_TYPES_H

// Value types for the Stage-1 relationship-topology builder and pedigree likelihood.
// This header has no likelihood factors, numerical plans, admissions, or
// resource limits.
#include "topology_templates.h"
#include <array>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace relationship_topology {

enum class TopologyStatus { Ready, Inapplicable, Unsupported };

struct TopologyNode {
  int typed_id;
  bool observed;
  bool neutral_phi;
  bool has_birth;
  std::array<int, 2> parents;
  std::string label;
  int source_input_index;

  TopologyNode()
      : typed_id(0), observed(false), neutral_phi(false), has_birth(false),
        parents{{-1, -1}}, source_input_index(-1) {}
};

// Life-history rows use original typed_id - 1 indexing. Borrowed vectors are
// synchronous inputs only; retaining boundaries must copy them into owned
// vectors. A null Input::life means life-history checks are disabled.
struct Life {
  std::vector<int> sex;
  std::vector<int> birth_lo;
  std::vector<int> birth_hi;
  std::vector<int> year_last;
  const std::vector<int>* borrowed_sex = nullptr;
  const std::vector<int>* borrowed_birth_lo = nullptr;
  const std::vector<int>* borrowed_birth_hi = nullptr;
  const std::vector<int>* borrowed_year_last = nullptr;
  std::vector<std::pair<int, int>> sex_overrides;
  std::array<double, 2> min_age;
  std::array<double, 2> max_age;

  Life() : min_age{{0.0, 0.0}},
           max_age{{std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()}} {}
};

struct Input {
  std::vector<TopologyNode> nodes;
  int focal;
  int candidate;
  int proposed_parent_side;
  const Life* life;
  std::size_t max_typed_nodes;
  std::size_t max_nodes;
  int optional_source_birth_budget = -1;
  std::vector<int> retained_source_birth_ids;
  std::vector<std::array<int, 2>> structural_source_parents;

  Input() : focal(-1), candidate(-1), proposed_parent_side(0), life(nullptr),
            max_typed_nodes(128U), max_nodes(256U) {}
};

struct StructuralBirth {
  int child = -1;
  std::array<int, 2> parents{{-1, -1}};
  bool required_by_motif = false;
  bool protected_source = false;
  bool cut = false;
};

struct TopologyModel {
  std::string name;
  std::string relationship;
  int model_order;
  int orientation_index;
  bool focal_is_descendant;
  int bridge_side;
  int shared_side;
  TopologyStatus status;
  std::string reason;
  std::vector<TopologyNode> nodes;
  int focal;
  int candidate;
  std::vector<int> original_to_canonical;
  std::vector<int> topological_order;
  std::vector<std::string> provenance;
  std::vector<StructuralBirth> structural_births;

  TopologyModel() : model_order(-1), orientation_index(-1),
                    focal_is_descendant(true), bridge_side(0), shared_side(0),
                    status(TopologyStatus::Unsupported), focal(-1), candidate(-1) {}
  bool ready() const { return status == TopologyStatus::Ready; }
};

using TopologyTemplateCache = pap_internal::topology_templates::Cache<TopologyModel>;

// Visit-scoped handoff of finished, life-validated topology models. This is
// exact source-input reuse, not a topology-template or cross-mutation cache.
struct FinishedModelTransfer {
  struct Entry { std::string input; TopologyModel model; };
  std::map<std::string, Entry> entries;
  bool reader = false;
  std::size_t offers = 0, requests = 0, hits = 0, mismatches = 0, missing = 0;
  double key_seconds = 0, copy_seconds = 0;
  std::size_t bytes = 0, peak_bytes = 0, refusals = 0;
  static inline thread_local FinishedModelTransfer* active = nullptr;
  struct Scope {
    FinishedModelTransfer* old;
    explicit Scope(FinishedModelTransfer* value) : old(active) { active = value; }
    ~Scope() { active = old; }
  };
  struct Read {
    FinishedModelTransfer* owner;
    bool old = false;
    Read() : owner(active) { if (owner) { old = owner->reader; owner->reader = true; } }
    ~Read() { if (owner) owner->reader = old; }
  };
  void clear() { entries.clear(); bytes = 0; }
};

}  // namespace relationship_topology

// Topology-only catalog published by the serial Stage-1 owner. The numerical
// pedigree likelihood program and its result remain owned by the live scorer.
#endif
