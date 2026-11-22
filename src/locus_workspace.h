#pragma once
#include "visit_graph.h"
namespace pap_internal::search {
struct LocusWorkspace {

  VisitGraph& pass;
  const LocusTables& locus;
  int locus_index;

  std::shared_ptr<const pap_internal::side_likelihood::Boundary>
    common_last_boundary;
  std::array<double, 25> common_last_partition;
  std::array<double, 25> common_last_error_bound;
  bool common_last_scaled_fast = false;
  bool common_last_pair_fast = false;
  pap_internal::side_staging::Phase common_last_phase =
    pap_internal::side_staging::Phase::complete;
  const pap_internal::side_staging::Prefix* common_last_prefix = nullptr;
  std::array<unsigned char, 25> common_last_status;
  std::array<std::string, 25> common_last_reason;
  bool common_last_supported = false;
  std::shared_ptr<const pap_internal::predictive_likelihood::RawBoundary>
    raw_last_boundary;
  bool raw_last_scaled_fast = false;
  // Scheduling-owned immutable Phi/P/T payloads, bounded per worker.
  std::vector<::relationship_topology::configuration::PredictiveResult> raw_last_values;
  bool common_primitive_certificate_checked = false;
  bool common_primitive_certificate_valid = false;
  double common_population_log_bound = 0.0;
  bool common_genotype_primitives_checked = false;

  LocusWorkspace(VisitGraph& pass_, int locus_index_);


  IbdClass ibd(int first, int second);

  // Per-locus partitions on the immutable shared relationship boundary.
  const std::array<double, 25>& common_side_partitions(
      const RelationshipContext& context,
      pap_internal::side_staging::Phase phase =
        pap_internal::side_staging::Phase::complete,
      const pap_internal::side_staging::Prefix* prefix = nullptr) {
    namespace stage = pap_internal::side_staging;
    if (common_last_boundary == context.common_boundary &&
        common_last_phase == phase && common_last_prefix == prefix &&
        common_last_scaled_fast == pass.options.relationship_common_scaled_fast &&
        common_last_pair_fast == pass.options.relationship_common_pair_fast)
      return common_last_partition;
    common_last_boundary = context.common_boundary;
    common_last_phase = phase;
    common_last_prefix = prefix;
    common_last_scaled_fast = pass.options.relationship_common_scaled_fast;
    common_last_pair_fast = pass.options.relationship_common_pair_fast;
    common_last_partition.fill(NA_REAL);
    common_last_error_bound.fill(0.0);
    // 0=inapplicable, 1=positive, 2=zero, 3=unavailable.
    common_last_status.fill(3);
    common_last_supported = common_last_boundary != nullptr &&
      common_last_boundary->failure.empty() &&
      common_last_boundary->models.size() == 25U;

    if (!common_last_supported) {

      return common_last_partition;
    }
    const auto& boundary = *common_last_boundary;
    const bool saved_prefix = phase == stage::Phase::remainder && prefix &&
      prefix->reusable(context.common_boundary.get(), locus_index);
    const stage::Phase evaluation_phase = phase == stage::Phase::remainder && !saved_prefix ?
      stage::Phase::complete : phase;
    if (pass.message_scoring && pass.options.private_test_message_comparison_decline) {
      // Keep the private comparison seam's deliberate one-model refusal.
      common_last_supported = false;
      common_last_reason[1] = "Private one-model message comparison decline";
      pass.telemetry.message_partition_unavailable += 1.0;

      return common_last_partition;
    }
    if (saved_prefix) {
      for (std::size_t model = 0U; model != 3U; ++model) {
        common_last_partition[model] = prefix->partition[model];
        common_last_error_bound[model] = prefix->error[model];
        common_last_status[model] = prefix->status[model];
        common_last_reason[model] = boundary.models[model].reason;
      }
      common_last_supported = prefix->supported;
    }
    const pap_internal::native_demand::Tags full_tags{
      1, int(evaluation_phase), context.focal, context.candidate,
      context.proposed_parent_side, context.proposed_dam, context.proposed_sire};
    auto* live = pap_internal::native_demand::live;
    if (!live)
      throw std::logic_error("Stage-1 ordinary scoring requires its live M_ped2 coordinator");
    boundary.native_support = live->native_support;
    if (locus_index == 0 && live->prepared_score_batches) {
      std::vector<std::pair<int, pap_internal::native_demand::Tags>> requests;
      for (std::size_t m = 0; m < boundary.models.size(); ++m)
        if (stage::selected(m, evaluation_phase))
          requests.push_back({boundary.models[m].topology.model_order, full_tags});
      live->prepare_scores(requests);
    }
    for (std::size_t m = 0; m < boundary.models.size(); ++m) {
      if (!stage::selected(m, evaluation_phase)) continue;
      pap_internal::native_demand::Tags tags{
        1, int(evaluation_phase), context.focal, context.candidate,
        context.proposed_parent_side, context.proposed_dam, context.proposed_sire};
      auto value = live->score(boundary.models[m].topology.model_order,
                               tags, locus_index, locus.k);
      live->observe_gate(tags, boundary.models[m].topology.model_order,
                         locus_index, locus.k, boundary.models[m], value);
      common_last_partition[m] = value.log_partition;
      common_last_status[m] = int(value.state);
      common_last_error_bound[m] = value.reference_error;
      common_last_reason[m] = live->engine_move_results ?
        std::move(value.reason) : value.reason;
      if (!value.available()) common_last_supported = false;
    }
    if (common_last_status[1] == 0 ||
        (common_last_status[0] == 1 && common_last_status[1] == 2))
      common_last_supported = false;
    return common_last_partition;
  }

  ::relationship_topology::configuration::PredictiveResult raw_common_predictive(
      const std::shared_ptr<const pap_internal::predictive_likelihood::RawBoundary>& request,
      std::size_t configuration_index, bool force_reference = false);

  std::vector<pap_internal::side_likelihood::PartitionValue>
  common_full_partitions(const RelationshipContext& context,
      pap_internal::panel_staging::Phase phase =
        pap_internal::panel_staging::Phase::complete);

  // Reduce biological orientations to the six relationship totals.
  std::array<double, 6> bounded_contextual_side_members(
      const RelationshipContext& context,
      std::array<double, 6>* partition_out,
      std::array<double, 22>* orientation_partition_out,
      pap_internal::side_staging::Phase phase,
      const pap_internal::side_staging::Prefix* prefix);

  double child_cell(int child, int dam_genotype, int sire_genotype) const;

  double child_genotype_cell(int child, int dam_genotype, int sire_genotype,
                             int child_genotype) const;



};
}
