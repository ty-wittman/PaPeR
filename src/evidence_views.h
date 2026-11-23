#pragma once
#include "search_common.h"
namespace pap_internal::search {
extern thread_local bool native_error_scope;

template <typename... Args>
[[noreturn]] static void engine_stop(const char* message, Args&&... args) {
  if (native_error_scope) {
    std::ostringstream rendered;
    rendered << message;
    using unpack = int[];
    (void)unpack{0, ((void)(rendered << " " << args), 0)...};
    throw std::runtime_error(rendered.str());
  }
  Rcpp::stop(message, std::forward<Args>(args)...);
}
struct NativeErrorScope {
  const bool previous;
  NativeErrorScope();
  ~NativeErrorScope();
  NativeErrorScope(const NativeErrorScope&) = delete;
  NativeErrorScope& operator=(const NativeErrorScope&) = delete;
};













bool finite_or_neg_inf(double value);

double log_add(double first, double second);

double log_sum(const std::vector<double>& values, bool sorted,
                      const char* context);

int genotype_index_zero(int first, int second, int k);

void normalize_in_place(std::vector<double>& value);

enum RelationCode {
  ROUTE_NOT_APPLICABLE = 0,
  ROUTE_UNCONNECTED = 1,
  ROUTE_PARENT_OFFSPRING = 2,
  ROUTE_FULL_SIBLINGS = 3,
  ROUTE_HALF_SIBLINGS = 4,
  ROUTE_DIRECT_ANCESTOR = 5,
  ROUTE_SHARED_ANCESTRY = 6,
  ROUTE_DUPLICATE = 7
};

struct IbdClass {
  double k0 = 1.0;
  double k1 = 0.0;
  double k2 = 0.0;
  int relation = ROUTE_UNCONNECTED;
};

struct SideParentSlots {
  int dam = 0;
  int sire = 0;
  bool has_dam() const;
  bool has_sire() const;
  bool has_any() const;
};

struct RelationshipContext {
  int focal = 0;
  int candidate = 0;
  // The proposal is deliberately carried with the ordinary descriptor rather
  // than reconstructed by each motif.  The tested candidate edge is absent
  // from the alternative topology; the other fixed proposed parent remains a
  // typed focal-parent slot in U/FS/HS/GP/FA/HA.
  int proposed_parent_side = 0;
  int proposed_dam = 0;
  int proposed_sire = 0;
  int focal_endpoint_excluded_child = 0;
  int candidate_endpoint_excluded_child = 0;
  SideParentSlots focal_parents;
  SideParentSlots candidate_parents;
  SideParentSlots focal_dam_parents;
  SideParentSlots focal_sire_parents;
  SideParentSlots candidate_dam_parents;
  SideParentSlots candidate_sire_parents;
  unsigned char fs_orientation_mask = 0;
  unsigned char hs_orientation_mask = 0;
  unsigned char gp_orientation_mask = 0;
  unsigned char fa_orientation_mask = 0;
  unsigned char ha_orientation_mask = 0;
  bool has_relevant_context = false;
  std::uint64_t signature = 0;
  // The immutable boundary identity owns its exact graph and primitive panel.
  // This is proof provenance, not a global cache key or evidence exclusion.
  const void* common_panel_owner = nullptr;
  std::uint64_t common_graph_revision = 0U;
  std::uint64_t common_message_revision = 0U;
  // Immutable, proposal-aware original-observation boundary. No published
  // marginal or per-class exclusion can redefine its evidence ownership.
  std::shared_ptr<const pap_internal::side_likelihood::Boundary>
    common_boundary;
};

struct LifeHistoryView {
  int n = 0;
  // Optional one-based counts for the current immutable pass graph. Production
  // supplies SearchState ownership; private standalone callers use null.
  const std::vector<std::array<int, 2>>* role_counts = nullptr;
  const std::vector<int>* sex = nullptr;
  const std::vector<int>* birth_lo = nullptr;
  const std::vector<int>* birth_hi = nullptr;
  const std::vector<int>* year_last = nullptr;
  const std::vector<double>* min_age = nullptr;
  const std::vector<double>* max_age = nullptr;
};

enum class SideEvidencePolicy : unsigned char {
  connector_direct = 0,
  endpoint_retain_external_descendants = 1
};







}
