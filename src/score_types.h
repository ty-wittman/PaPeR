#ifndef PAP_SCORE_TYPES_H
#define PAP_SCORE_TYPES_H
#include "topology_types.h"
namespace pap_internal::side_likelihood {

struct FullBoundary;

struct Model {
  ::relationship_topology::TopologyModel topology;
  bool supported = false;
  std::string reason;
};

struct Boundary {
  mutable bool native_support = false;
  std::vector<int> evidence_ids;
  std::vector<Model> models;
  std::string failure;
  std::shared_ptr<const FullBoundary> full_boundary;
};

enum class PartitionState : unsigned char {
  inapplicable = 0, positive = 1, zero = 2, unavailable = 3
};

struct PartitionValue {
  PartitionState state = PartitionState::unavailable;
  double log_partition = std::numeric_limits<double>::quiet_NaN();
  double reference_error = 0.0;
  bool message_approximate = false;
  std::string reason;
  bool available() const { return state != PartitionState::unavailable; }
};

struct PanelCell {
  PartitionState state = PartitionState::unavailable;
  double log_partition = std::numeric_limits<double>::quiet_NaN();
  double reference_error = 0.0;
  std::size_t loci = 0U;
  std::string reason;
};

struct FullMemberDescriptor {
  int model_order = -1;
  std::string name, relationship;
  bool parent_numerator = false;
  bool reverse = false;
  bool ownership_review_required = false;
  bool registered_endpoint_orientation = true;
  bool orientation_review_required = false;
  unsigned int parent_slot_bits = 0U;
  ::relationship_topology::TopologyStatus topology_status =
    ::relationship_topology::TopologyStatus::Unsupported;
  std::string reason;
  std::size_t compiled_index = std::numeric_limits<std::size_t>::max();
};

struct FullBoundary {
  Boundary compiled;
  std::vector<FullMemberDescriptor> members;
  std::vector<std::size_t> compiled_to_member;
  int tested_side = 0;
  int focal_typed_id = 0;
  int candidate_typed_id = 0;
  std::string failure;
};

struct FullPanelTotals {
  const FullBoundary* boundary = nullptr;
  std::size_t expected_loci = 0U;
  std::size_t loci = 0U;
  std::vector<PanelCell> cells;
  std::string failure;

  FullPanelTotals() = default;
  FullPanelTotals(const FullBoundary& source, std::size_t count)
      : boundary(&source), expected_loci(count), cells(source.compiled.models.size()) {}
};

}  // namespace pap_internal::side_likelihood

namespace pap_internal::predictive_likelihood {

struct CommonBoundary {
  ::relationship_topology::Input input;
  bool ready = false;
};

struct Request {
  enum class BuildStatus : unsigned char { ready, hard_invalid, unavailable };
  int dam = -1, sire = -1;
  bool phased_refused = false;
  BuildStatus build_status = BuildStatus::unavailable;
  std::string reason;
  std::size_t numerator = std::numeric_limits<std::size_t>::max();
  std::size_t denominator = std::numeric_limits<std::size_t>::max();
  bool ready() const {
    return build_status == BuildStatus::ready &&
      numerator != std::numeric_limits<std::size_t>::max() &&
      denominator != std::numeric_limits<std::size_t>::max();
  }
};

struct RawBoundary {
  CommonBoundary common;
  side_likelihood::Boundary compiled;
  std::vector<Request> requests;
  std::string failure;
  bool ready() const { return common.ready && failure.empty(); }
};

}  // namespace pap_internal::predictive_likelihood
#endif
