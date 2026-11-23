#include "pap_removed_path.h"
#ifndef PAP_PREDICTIVE_MODELS_H
#define PAP_PREDICTIVE_MODELS_H

// Symbolic raw Stage-1 configuration likelihoods on ONE caller-declared local
// observation union.  No production hook, kernel, scorer or state mutation.
// The caller must mask BOTH incoming focal edges before freezing this union;
// every other typed birth, observation owner and supplied latent identity stays
// fixed.  Parent proposals change only the focal birth.  The raw estimand is
// log Z(Phi_A,C | H) - log Z(C | H), not a relationship partition or a published
// endpoint cavity.  A is retained as a latent genotype when C uses its children.

#include "relationship_topology.h"
#include "predictive_error_bounds.h"
#include <memory>
#include <utility>

namespace relationship_topology {
namespace configuration {

enum class BuildStatus { Ready, HardInvalid, ConditioningZero, Unavailable };
enum class Support { Positive, Zero, Unavailable };

// Tokens own no data. Their stable distinct addresses prevent accidentally
// reusing a with-Phi_A payload/partition as a without-Phi_A input.
struct BoundaryToken {};
// The caller supplies the immutable primitive-table snapshot identity and
// current locus. A new token is required when either changes. Numeric support
// is not a property of the graph alone (in particular at exact-zero error).
struct PrimitiveToken {
  const void* snapshot;
  int locus;
  PrimitiveToken(const void* immutable_snapshot, int locus_index)
      : snapshot(immutable_snapshot), locus(locus_index) {}
  bool valid() const { return snapshot != NULL && locus >= 0; }
};
struct ObservationToken {
  std::shared_ptr<const BoundaryToken> boundary;
  std::shared_ptr<const PrimitiveToken> primitive;
  bool includes_focal_observation;
  ObservationToken(std::shared_ptr<const BoundaryToken> b,
                   std::shared_ptr<const PrimitiveToken> p, bool includes)
      : boundary(std::move(b)), primitive(std::move(p)),
        includes_focal_observation(includes) {}
};

struct Boundary {
  Input input;
  std::shared_ptr<const Life> owned_life;
  std::shared_ptr<const BoundaryToken> token;
  std::vector<int> dam_candidates;
  std::vector<int> sire_candidates;
  std::vector<int> evidence_ids;
  // Original input indices; components are independent under the common
  // focal-masked DAG. Attaching a focal parent joins only affected components.
  std::vector<std::vector<int>> components;
  bool ready = false;
  std::string reason;
};

struct ComponentSupport {
  std::shared_ptr<const BoundaryToken> boundary;
  std::shared_ptr<const PrimitiveToken> primitive;
  std::vector<Support> components;
};

struct ModelPair {
  BuildStatus status = BuildStatus::Unavailable;
  std::string reason;
  TopologyModel numerator;
  TopologyModel denominator;
  std::shared_ptr<const ObservationToken> numerator_token;
  std::shared_ptr<const ObservationToken> denominator_token;
  std::vector<int> included_input_indices;
  std::vector<std::size_t> canceled_components;
  std::vector<int> common_evidence_ids;
  // Local indices in the frozen input. -1 is a missing parent, never a typed ID.
  int dam = -1;
  int sire = -1;
  bool ready() const { return status == BuildStatus::Ready; }
};

namespace detail {

inline std::vector<int> sorted_unique(std::vector<int> values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
  return values;
}

inline bool contains(const std::vector<int>& values, int value) {
  return std::binary_search(values.begin(), values.end(), value);
}

inline bool directed_acyclic(const std::vector<TopologyNode>& nodes) {
  std::vector<int> degree(nodes.size(), 0);
  std::vector<std::vector<int>> children(nodes.size());
  std::vector<int> ready;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    for (int parent : nodes[i].parents)
      if (parent >= 0) {
        ++degree[i];
        children[static_cast<std::size_t>(parent)].push_back(static_cast<int>(i));
      }
    if (degree[i] == 0) ready.push_back(static_cast<int>(i));
  }
  std::size_t visited = 0U;
  while (!ready.empty()) {
    const int parent = ready.back();
    ready.pop_back();
    ++visited;
    for (int child : children[static_cast<std::size_t>(parent)])
      if (--degree[static_cast<std::size_t>(child)] == 0) ready.push_back(child);
  }
  return visited == nodes.size();
}

inline void validate_one_focal_life(const Input& input) {
  if (input.life == NULL) return;
  const int id = input.nodes[static_cast<std::size_t>(input.focal)].typed_id;
  const auto& life = *input.life;
  // A sole population founder has no parent roles or age edges. Still validate
  // its declared interval and field coverage without inventing a second ID.
  (void)topology_detail::life_sex_value(life, id);
  (void)topology_detail::life_value(life.year_last, life.borrowed_year_last, id);
  const int lo = topology_detail::life_value(life.birth_lo, life.borrowed_birth_lo, id);
  const int hi = topology_detail::life_value(life.birth_hi, life.borrowed_birth_hi, id);
  if (lo != INT_MIN && hi != INT_MIN && lo > hi)
    throw topology_detail::Failure(TopologyStatus::Inapplicable,
                                  "Sole focal birth interval is infeasible");
}

// Dependency-cut raw inputs retain the ORIGINAL bounded structural rows as a
// separate validity view. Numerical component cancellation cannot discard a
// parent role, age constraint or descendant path needed to reject a proposal.
// This is metadata-only validation using the existing native Builder. It does
// not restore any cut factor to the numerical model or extract extra ancestry.
inline void validate_structural_source_proposal(const Input& input,
                                               int dam = -1, int sire = -1) {
  if (input.structural_source_parents.empty()) {
    if (input.optional_source_birth_budget >= 0)
      throw topology_detail::Failure(TopologyStatus::Unsupported,
                                    "Budgeted raw input lacks its structural source validity snapshot");
    return;  // Original exact path.
  }
  if (input.structural_source_parents.size() != input.nodes.size() ||
      input.focal < 0 || static_cast<std::size_t>(input.focal) >= input.nodes.size())
    throw topology_detail::Failure(TopologyStatus::Unsupported,
                                  "Raw structural source snapshot has invalid ownership/dimensions");
  const auto focal_parents = input.structural_source_parents[static_cast<std::size_t>(input.focal)];
  if (focal_parents[0] >= 0 || focal_parents[1] >= 0)
    throw topology_detail::Failure(TopologyStatus::Unsupported,
                                  "Raw structural source focal birth must be masked on both sides");
  if (dam < -1 || sire < -1 ||
      (dam >= 0 && static_cast<std::size_t>(dam) >= input.nodes.size()) ||
      (sire >= 0 && static_cast<std::size_t>(sire) >= input.nodes.size()))
    throw topology_detail::Failure(TopologyStatus::Unsupported,
                                  "Raw structural parent proposal is outside bounded source");
  if (dam == input.focal || sire == input.focal || (dam >= 0 && dam == sire))
    throw topology_detail::Failure(TopologyStatus::Inapplicable,
                                  "Self or identical typed structural parent proposal");
  Input structural = input;
  // Disable dependency cuts while proving structural validity. Every snapshot
  // row is owned by this same local input; no global pedigree is consulted.
  structural.optional_source_birth_budget = -1;
  structural.retained_source_birth_ids.clear();
  structural.structural_source_parents.clear();
  for (std::size_t child = 0U; child < structural.nodes.size(); ++child) {
    auto& node = structural.nodes[child];
    node.parents = input.structural_source_parents[child];
    node.has_birth = node.parents[0] >= 0 || node.parents[1] >= 0;
    for (int parent : node.parents)
      if (parent < -1 || (parent >= 0 && static_cast<std::size_t>(parent) >= structural.nodes.size()))
        throw topology_detail::Failure(TopologyStatus::Unsupported,
                                      "Raw structural source parent is outside bounded source");
  }
  if (!directed_acyclic(structural.nodes))
    throw topology_detail::Failure(TopologyStatus::Unsupported,
                                  "Raw masked structural source contains a directed cycle");
  structural.candidate = -1;
  for (std::size_t node = 0U; node < structural.nodes.size(); ++node)
    if (static_cast<int>(node) != structural.focal && structural.nodes[node].typed_id > 0) {
      structural.candidate = static_cast<int>(node); break;
    }
  if (structural.candidate < 0) {
    if (structural.nodes.size() != 1U || dam >= 0 || sire >= 0)
      throw topology_detail::Failure(TopologyStatus::Unsupported,
                                    "Nontrivial raw structural validity requires two real typed anchors");
    validate_one_focal_life(structural);
    return;
  }
  topology_detail::Builder builder(structural);
  if (dam >= 0) builder.unify_nodes(builder.parent_slot(structural.focal, 1), dam);
  if (sire >= 0) builder.unify_nodes(builder.parent_slot(structural.focal, 2), sire);
  TopologyModel verified;
  verified.name = "raw_structural_validity_only";
  verified.relationship = "raw_configuration";
  builder.finish(verified);
}

}  // namespace detail

// Input.nodes is the immutable bounded union, not a candidate-specific extract.
// The supplied numeric topology caps apply to EACH affected contraction. The
// symbolic union may be larger and then cancel independent positive components.
// max_union_nodes bounds only this topology/ownership metadata; it admits no
// numeric contraction and never changes Input's original per-model caps.
inline Boundary freeze(const Input& input, std::vector<int> dam_candidates,
                       std::vector<int> sire_candidates,
                       std::size_t max_union_nodes = 4096U) {
  Boundary result;
  try {
    if (input.nodes.empty() || input.max_nodes == 0U || input.max_typed_nodes == 0U ||
        max_union_nodes == 0U || input.nodes.size() > max_union_nodes)
      throw topology_detail::Failure(TopologyStatus::Unsupported,
                                    "Empty or unsupported symbolic configuration boundary");
    if (input.focal < 0 || static_cast<std::size_t>(input.focal) >= input.nodes.size() ||
        input.nodes[static_cast<std::size_t>(input.focal)].typed_id <= 0)
      throw topology_detail::Failure(TopologyStatus::Unsupported,
                                    "Raw configuration focal is outside typed boundary");
    if (input.nodes[static_cast<std::size_t>(input.focal)].parents[0] >= 0 ||
        input.nodes[static_cast<std::size_t>(input.focal)].parents[1] >= 0)
      throw topology_detail::Failure(TopologyStatus::Unsupported,
                                    "Both incoming focal edges must be masked before freezing C");
    result.input = input;
    result.input.nodes[static_cast<std::size_t>(input.focal)].has_birth = false;
    if (input.life != NULL) {
      result.owned_life = std::make_shared<const Life>(
        topology_detail::materialize_life(*input.life));
      result.input.life = result.owned_life.get();
    }
    result.dam_candidates = detail::sorted_unique(std::move(dam_candidates));
    result.sire_candidates = detail::sorted_unique(std::move(sire_candidates));
    std::set<int> typed;
    std::vector<std::vector<int>> adjacency(input.nodes.size());
    for (std::size_t i = 0; i < input.nodes.size(); ++i) {
      auto& node = result.input.nodes[i];
      if (node.typed_id < 0 || (node.typed_id == 0 && node.observed) ||
          (node.typed_id > 0 && (!node.observed || !typed.insert(node.typed_id).second)))
        throw topology_detail::Failure(TopologyStatus::Unsupported,
                                      "Common typed IDs require exactly one observation owner");
      if (node.typed_id > 0) result.evidence_ids.push_back(node.typed_id);
      for (int parent : node.parents) {
        if (parent < -1 || (parent >= 0 && static_cast<std::size_t>(parent) >= input.nodes.size()))
          throw topology_detail::Failure(TopologyStatus::Unsupported,
                                        "Retained parent is outside common configuration union");
        if (parent >= 0) {
          node.has_birth = true;
          adjacency[i].push_back(parent);
          adjacency[static_cast<std::size_t>(parent)].push_back(static_cast<int>(i));
        }
      }
    }
    for (const auto* candidates : {&result.dam_candidates, &result.sire_candidates})
      for (int index : *candidates)
        if (index < 0 || static_cast<std::size_t>(index) >= input.nodes.size() ||
            input.nodes[static_cast<std::size_t>(index)].typed_id <= 0)
          throw topology_detail::Failure(TopologyStatus::Unsupported,
                                        "Candidate is outside frozen typed configuration union");
    if (!detail::directed_acyclic(result.input.nodes))
      throw topology_detail::Failure(TopologyStatus::Unsupported,
                                    "Masked common background contains a directed cycle");
    detail::validate_structural_source_proposal(result.input);
    std::vector<unsigned char> visited(input.nodes.size(), 0U);
    for (std::size_t i = 0; i < input.nodes.size(); ++i) {
      if (visited[i]) continue;
      std::vector<int> component, pending{static_cast<int>(i)};
      while (!pending.empty()) {
        const int node = pending.back();
        pending.pop_back();
        if (visited[static_cast<std::size_t>(node)]) continue;
        visited[static_cast<std::size_t>(node)] = 1U;
        component.push_back(node);
        for (int neighbor : adjacency[static_cast<std::size_t>(node)])
          if (!visited[static_cast<std::size_t>(neighbor)]) pending.push_back(neighbor);
      }
      result.components.push_back(detail::sorted_unique(std::move(component)));
    }
    std::sort(result.evidence_ids.begin(), result.evidence_ids.end());
    result.token = std::make_shared<const BoundaryToken>();
    result.ready = true;
  } catch (const topology_detail::Failure& error) {
    result.reason = error.what();
  } catch (const std::bad_alloc&) {
    result.reason = "Configuration boundary allocation failed";
  } catch (const std::length_error&) {
    result.reason = "Configuration boundary allocation length exceeded";
  }
  return result;
}

// A support certificate MUST come from the same frozen C and current locus:
// exact component support, or a sufficient positive-primitive proof. A failed
// sufficient proof is Unavailable, never Zero. This helper validates ownership;
// it does not manufacture evidence that supplied supports are correct.
inline ComponentSupport support_certificate(const Boundary& boundary,
                                            std::shared_ptr<const PrimitiveToken> primitive,
                                            std::vector<Support> supports) { pap_internal::removed_path("support_certificate"); }

inline ModelPair build(const Boundary& boundary, int dam = -1, int sire = -1,
                       const ComponentSupport* support = NULL,
                       bool cancel_unaffected = true,
                       std::shared_ptr<const PrimitiveToken> primitive = nullptr) {
  ModelPair result;
  result.dam = dam;
  result.sire = sire;
  result.common_evidence_ids = boundary.evidence_ids;
  try {
    if (!boundary.ready)
      throw topology_detail::Failure(TopologyStatus::Unsupported, boundary.reason);
    if (dam < -1 || sire < -1 ||
        (dam >= 0 && !detail::contains(boundary.dam_candidates, dam)) ||
        (sire >= 0 && !detail::contains(boundary.sire_candidates, sire)))
      throw topology_detail::Failure(TopologyStatus::Unsupported,
                                    "Proposal is outside frozen active candidate union");
    if (dam == boundary.input.focal || sire == boundary.input.focal ||
        (dam >= 0 && dam == sire))
      throw topology_detail::Failure(TopologyStatus::Inapplicable,
                                    "Self or identical typed parent proposal");
    detail::validate_structural_source_proposal(boundary.input, dam, sire);
    for (std::size_t c = 0; c < boundary.components.size(); ++c) {
      const auto& component = boundary.components[c];
      const bool affected = detail::contains(component, boundary.input.focal) ||
        (dam >= 0 && detail::contains(component, dam)) ||
        (sire >= 0 && detail::contains(component, sire));
      if (!cancel_unaffected || affected) {
        result.included_input_indices.insert(result.included_input_indices.end(),
                                             component.begin(), component.end());
      } else {
        result.canceled_components.push_back(c);
        if (primitive == nullptr || !primitive->valid() || support == NULL ||
            support->boundary != boundary.token || support->primitive != primitive ||
            support->components.size() != boundary.components.size() ||
            support->components[c] == Support::Unavailable)
          throw topology_detail::Failure(TopologyStatus::Unsupported,
                                        "Canceled common-C component has no positive support certificate");
        if (support->components[c] == Support::Zero) {
          result.status = BuildStatus::ConditioningZero;
          result.reason = "Common-C component has zero normalizer; predictive 0/0 cannot cancel";
          return result;
        }
      }
    }
    result.included_input_indices = detail::sorted_unique(std::move(result.included_input_indices));
    Input input;
    input.life = boundary.input.life;
    input.max_typed_nodes = boundary.input.max_typed_nodes;
    input.max_nodes = boundary.input.max_nodes;
    std::vector<int> remap(boundary.input.nodes.size(), -1);
    for (int original : result.included_input_indices) {
      remap[static_cast<std::size_t>(original)] = static_cast<int>(input.nodes.size());
      input.nodes.push_back(boundary.input.nodes[static_cast<std::size_t>(original)]);
    }
    for (auto& node : input.nodes)
      for (int& parent : node.parents)
        if (parent >= 0) {
          parent = remap[static_cast<std::size_t>(parent)];
          if (parent < 0)
            throw topology_detail::Failure(TopologyStatus::Unsupported,
                                          "Cancellation split a retained birth factor");
        }
    input.focal = remap[static_cast<std::size_t>(boundary.input.focal)];
    input.candidate = -1;
    for (std::size_t i = 0; i < input.nodes.size(); ++i)
      if (static_cast<int>(i) != input.focal && input.nodes[i].typed_id > 0) {
        input.candidate = static_cast<int>(i);
        break;
      }
    TopologyModel model;
    model.name = dam >= 0 ? (sire >= 0 ? "raw_pair" : "raw_dam_only") :
                           (sire >= 0 ? "raw_sire_only" : "raw_null");
    model.relationship = "raw_configuration";
    model.model_order = dam >= 0 ? (sire >= 0 ? 3 : 1) : (sire >= 0 ? 2 : 0);
    if (input.candidate < 0) {
      if (input.nodes.size() != 1U || dam >= 0 || sire >= 0 ||
          input.max_nodes < 1U || input.max_typed_nodes < 1U)
        throw topology_detail::Failure(TopologyStatus::Unsupported,
                                      "A nontrivial configuration needs two real typed anchors");
      detail::validate_one_focal_life(input);
      model.nodes = input.nodes;
      model.focal = 0;
      model.candidate = -1;
      model.original_to_canonical = {0};
      model.topological_order = {0};
      model.status = TopologyStatus::Ready;
      model.reason = "Sole focal normalized population founder; no invented observation anchor";
    } else {
      topology_detail::Builder builder(input);
      if (dam >= 0)
        builder.unify_nodes(builder.parent_slot(input.focal, 1), remap[static_cast<std::size_t>(dam)]);
      if (sire >= 0)
        builder.unify_nodes(builder.parent_slot(input.focal, 2), remap[static_cast<std::size_t>(sire)]);
      builder.finish(model);
    }
    model.provenance.push_back("Same frozen focal-masked C for null, dam, sire and pair");
    model.provenance.push_back("Only proposed focal birth changes; every other retained typed birth remains");
    if (!boundary.input.structural_source_parents.empty())
      model.provenance.push_back("Original bounded structural source validated before numerical dependency cuts");
    if (!result.canceled_components.empty())
      model.provenance.push_back("Independent common-C components cancel with certified positive normalizers");
    result.numerator = model;
    result.numerator.name += "_with_focal_Phi";
    result.denominator = std::move(model);
    result.denominator.name += "_without_focal_Phi";
    result.denominator.nodes[static_cast<std::size_t>(result.denominator.focal)].observed = false;
    result.denominator.provenance.push_back("Exclude only Phi_A; retain A genotype and all descendant birth factors");
    result.numerator_token = std::make_shared<const ObservationToken>(boundary.token, primitive, true);
    result.denominator_token = std::make_shared<const ObservationToken>(boundary.token, primitive, false);
    result.status = BuildStatus::Ready;
  } catch (const topology_detail::Failure& error) {
    result.status = error.status == TopologyStatus::Inapplicable ?
                    BuildStatus::HardInvalid : BuildStatus::Unavailable;
    result.reason = error.what();
  } catch (const std::bad_alloc&) {
    result.status = BuildStatus::Unavailable;
    result.reason = "Configuration topology allocation failed";
  } catch (const std::length_error&) {
    result.status = BuildStatus::Unavailable;
    result.reason = "Configuration topology allocation length exceeded";
  }
  return result;
}

enum class PredictiveStatus { DefinedFinite, DefinedGeneticZero, ConditioningZero, Unavailable };
struct Partition {
  Support support = Support::Unavailable;
  double ln = std::numeric_limits<double>::quiet_NaN();
  double reference_error = 0.0;
  bool message_approximate = false;
};
struct PredictiveResult {
  PredictiveStatus status = PredictiveStatus::Unavailable;
  double ln = std::numeric_limits<double>::quiet_NaN();
  std::string reason;
  double reference_error = 0.0;
  bool message_approximate = false;
};

inline PredictiveResult conditional(const Partition& numerator, const Partition& denominator) {
  PredictiveResult result;
  result.message_approximate = numerator.message_approximate || denominator.message_approximate;
  const auto valid = [](const Partition& p) {
    return p.support == Support::Positive ? std::isfinite(p.ln) :
      p.support == Support::Zero && p.ln == -std::numeric_limits<double>::infinity();
  };
  if (!valid(numerator) || !valid(denominator)) {
    result.reason = "Unscored, unavailable or nonfinite applicable partition";
  } else if (denominator.support == Support::Zero) {
    if (numerator.support == Support::Positive) {
      result.reason = "Positive numerator with zero own conditioning partition violates support inclusion";
    } else {
      result.status = PredictiveStatus::ConditioningZero;
      result.reason = "Undefined predictive 0/0";
    }
  } else if (numerator.support == Support::Zero) {
    result.status = PredictiveStatus::DefinedGeneticZero;
    result.ln = -std::numeric_limits<double>::infinity();
  } else {
    result.ln = numerator.ln - denominator.ln;
    if (std::isfinite(result.ln)) {
      result.status = PredictiveStatus::DefinedFinite;
      result.reference_error = pap_internal::predictive_error_bounds::predictive_error(
        result.ln, numerator.reference_error, denominator.reference_error);
    }
    else result.reason = "Finite partition subtraction exceeds numeric range";
  }
  return result;
}

}  // namespace configuration
}  // namespace relationship_topology

#endif
