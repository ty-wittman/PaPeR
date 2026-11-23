#include "locus_workspace.h"
namespace pap_internal::search {

LocusWorkspace::LocusWorkspace(VisitGraph& pass_, int locus_index_)
    : pass(pass_), locus(pass_.compiled.locus[locus_index_]), locus_index(locus_index_) {}

auto LocusWorkspace::ibd(int first, int second) -> IbdClass  {

    const IbdClass out = pass.ibd(first, second);

    return out;
  }

auto LocusWorkspace::raw_common_predictive(
      const std::shared_ptr<const pap_internal::predictive_likelihood::RawBoundary>& request,
      std::size_t configuration_index, bool force_reference ) -> ::relationship_topology::configuration::PredictiveResult  {

    namespace cfg = ::relationship_topology::configuration;
    if (!request || !request->ready() || configuration_index >= request->requests.size()) {
      cfg::PredictiveResult value;
      value.reason = request ? request->failure : "Missing raw common background";
      return value;
    }
    const bool fast = pass.raw_scaled_fast && !force_reference;
    if (!force_reference && raw_last_boundary == request &&
        raw_last_scaled_fast == fast)
      return raw_last_values[configuration_index];
    if (!force_reference) {
      raw_last_boundary = request;
      raw_last_scaled_fast = fast;
      raw_last_values.assign(request->requests.size(), cfg::PredictiveResult());
    }
    auto* live = pap_internal::native_demand::live;
    if (!live)
      throw std::logic_error("Stage-1 raw scoring requires its live M_ped2 coordinator");
    if (locus_index == 0 && live->prepared_score_batches) {
      std::vector<std::pair<int, pap_internal::native_demand::Tags>> jobs;
      for (std::size_t i = 0; i < request->requests.size(); ++i) {
        if (force_reference && i != configuration_index) continue;
        const auto& model = request->requests[i];
        if (!model.ready() && !(live->native_support && model.phased_refused)) continue;
        const int dam = model.dam < 0 ? 0 :
          request->common.input.nodes[model.dam].typed_id;
        const int sire = model.sire < 0 ? 0 :
          request->common.input.nodes[model.sire].typed_id;
        const int focal = request->common.input.nodes[
          request->common.input.focal].typed_id;
        jobs.push_back({-1, {3, 0, focal, 0, 0, dam, sire, 0}});
        jobs.push_back({-1, {3, 0, focal, 0, 0, dam, sire, 1}});
      }
      live->prepare_scores(jobs);
    }
    for (std::size_t i = 0; i < request->requests.size(); ++i) {
      if (force_reference && i != configuration_index) continue;
      const auto& model = request->requests[i];
      cfg::PredictiveResult value;
        if (!model.ready() && !(live->native_support && model.phased_refused)) {
        value.reason = model.reason;
      } else {
        const int dam = model.dam < 0 ? 0 :
          request->common.input.nodes[model.dam].typed_id;
        const int sire = model.sire < 0 ? 0 :
          request->common.input.nodes[model.sire].typed_id;
        const int focal = request->common.input.nodes[
          request->common.input.focal].typed_id;
        const auto numerator = live->score(
          -1, {3, 0, focal, 0, 0, dam, sire, 0}, locus_index, locus.k);
        const auto denominator = live->score(
          -1, {3, 0, focal, 0, 0, dam, sire, 1}, locus_index, locus.k);
        const auto convert = [](const pap_internal::side_likelihood::PartitionValue& scored) {
          return cfg::Partition{
            scored.state == pap_internal::side_likelihood::PartitionState::positive ? cfg::Support::Positive :
              scored.state == pap_internal::side_likelihood::PartitionState::zero ? cfg::Support::Zero :
                cfg::Support::Unavailable,
            scored.log_partition, scored.reference_error, scored.message_approximate};
        };
        value = cfg::conditional(convert(numerator), convert(denominator));
      }
      if (force_reference) {
        if (i == configuration_index) return value;
      } else {
        raw_last_values[i] = value;
      }
    }
    return raw_last_values[configuration_index];
  }

auto LocusWorkspace::common_full_partitions(const RelationshipContext& context,
      pap_internal::panel_staging::Phase phase ) -> std::vector<pap_internal::side_likelihood::PartitionValue>  {
    namespace b = pap_internal::side_likelihood;
    if (!context.common_boundary || !context.common_boundary->full_boundary)
      throw std::runtime_error("Full common request has no owned hypothesis boundary");
    const auto& full = *context.common_boundary->full_boundary;
    const auto& boundary = full.compiled;
    std::vector<b::PartitionValue> out(boundary.models.size());
    if (!full.failure.empty()) {
      for (auto& value : out) value.reason = full.failure;
      return out;
    }
    if (phase != pap_internal::panel_staging::Phase::complete &&
        !pap_internal::panel_staging::valid_map(full))
      throw std::runtime_error("Staged full request has an invalid immutable model map");
    auto* live = pap_internal::native_demand::live;
    if (!live)
      throw std::logic_error("Stage-1 full scoring requires its live M_ped2 coordinator");
    const pap_internal::native_demand::Tags full_tags{
      2, int(phase), context.focal, context.candidate,
      context.proposed_parent_side, context.proposed_dam, context.proposed_sire};
    boundary.native_support = live->native_support;
    if (locus_index == 0 && live->prepared_score_batches) {
      std::vector<std::pair<int, pap_internal::native_demand::Tags>> requests;
      for (std::size_t m = 0; m < boundary.models.size(); ++m)
        if (pap_internal::panel_staging::selected(full, m, phase))
          requests.push_back({boundary.models[m].topology.model_order, full_tags});
      live->prepare_scores(requests);
    }
    for (std::size_t m = 0; m < boundary.models.size(); ++m) {
      if (!pap_internal::panel_staging::selected(full, m, phase)) continue;
      const auto tags = pap_internal::native_demand::Tags{
        2, int(phase), context.focal, context.candidate,
        context.proposed_parent_side, context.proposed_dam, context.proposed_sire};
      auto value = live->score(boundary.models[m].topology.model_order,
                               tags, locus_index, locus.k);
      live->observe_gate(tags, boundary.models[m].topology.model_order,
                         locus_index, locus.k, boundary.models[m], value);
      out[m].state = static_cast<b::PartitionState>(value.state);
      out[m].log_partition = value.log_partition;
      out[m].reference_error = value.reference_error;
      out[m].reason = live->engine_move_results ?
        std::move(value.reason) : value.reason;
    }
    return out;
  }

auto LocusWorkspace::bounded_contextual_side_members(
      const RelationshipContext& context,
      std::array<double, 6>* partition_out,
      std::array<double, 22>* orientation_partition_out,
      pap_internal::side_staging::Phase phase,
      const pap_internal::side_staging::Prefix* prefix) -> std::array<double, 6>  {
    if (!context.common_boundary)
      throw std::logic_error("Stage-1 ordinary members require their live topology catalog");
    const auto& raw = common_side_partitions(context, phase, prefix);
    std::array<double, 6> partition;
    std::array<double, 22> orientations;
    partition.fill(R_NegInf);
    orientations.fill(R_NegInf);
    partition[0] = raw[1];
    partition[1] = raw[2];
    for (int orientation = 0; orientation < 22; ++orientation)
      orientations[orientation] = raw[static_cast<std::size_t>(orientation + 3)];
    const int begin[] = {0, 0, 0, 2, 10, 14};
    const int end[] = {0, 0, 2, 10, 14, 22};
    for (int relationship = 2; relationship < 6; ++relationship)
      for (int o = begin[relationship]; o < end[relationship]; ++o)
        if (!std::isnan(orientations[o]))
          partition[relationship] = std::max(partition[relationship], orientations[o]);
    if (partition_out != nullptr) *partition_out = partition;
    if (orientation_partition_out != nullptr) *orientation_partition_out = orientations;
    std::array<double, 6> out;
    out.fill(R_NegInf);
    // A zero U cannot define even the U/U display representation. Keep
    // that undefined boundary distinct from a positive U's neutral gauge.
    out[0] = R_finite(partition[0]) ? 0.0 : R_NegInf;
    if (common_last_supported && R_finite(partition[0]))
      for (int relationship = 1; relationship < 6; ++relationship)
        if (R_finite(partition[relationship]))
          out[relationship] = partition[relationship] - partition[0];
    return out;
    }

auto LocusWorkspace::child_cell(int child, int dam_genotype, int sire_genotype) const -> double  {
    const LatentMendelianChildSupport support = locus_mendelian_child_support(
      locus, dam_genotype, sire_genotype);

    LogSumExpAccumulator accumulator;
    for (int i = 0; i < support.count; ++i) {
      const int child_genotype = support.genotype[i];
      accumulator.add_log(
        child_genotype_cell(
          child, dam_genotype, sire_genotype, child_genotype),
        "simp compiled child-observation contraction");
    }
    return accumulator.value();
  }

auto LocusWorkspace::child_genotype_cell(int child, int dam_genotype, int sire_genotype,
                             int child_genotype) const -> double  { pap_internal::removed_path("child_genotype_cell"); }
}
