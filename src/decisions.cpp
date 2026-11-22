#include "decisions.h"
#include "scoring.h"
namespace pap_internal::search {


void compute_diagnostics(
    const SearchInputs& input, int current_index,
    std::vector<ConfigurationScore>& configurations,
    const std::vector<SideReduction>* precomputed,
    std::vector<SideReduction>* numerical_out) {

  const SingleConfigurationIndices single_indices =
    single_configuration_indices(configurations);
  if (current_index < 0 || current_index >= static_cast<int>(configurations.size()))
    engine_stop("whole-pass diagnostics received an invalid current configuration");
  const ConfigurationScore& incumbent = configurations[current_index];
  for (auto& row : configurations) set_incumbent_policy(input, incumbent, row);
  // The descriptor-only chooser seam intentionally supplies no relationship
  // registry or member storage. Production simp and full calls both supply
  // the shared payload and therefore take this reduction.
  if (!input.full_relationship_delta.empty()) {
    int null_index = -1;
    for (std::size_t i = 0; i < configurations.size(); ++i)
      if (configurations[i].nparents == 0) { null_index = static_cast<int>(i); break; }
    if (null_index < 0) engine_stop("full relationship reduction lost null row");
    const double null_objective = input.elide_raw() ? 0.0 : configurations[null_index].objective_ln;
    std::vector<SideReduction> reduction(configurations.size());
    for (std::size_t i = 0; i < configurations.size(); ++i) {
      if (configurations[i].raw_support_rejected || configurations[i].gate_deferred) {
        // An independently proven support negative needs no relationship answer,
        // and a row deferred by the gate window has none. Do not publish an
        // empty/unrequested member panel as an evaluated one.
        configurations[i].full_operational_ln = NA_REAL;
        configurations[i].full_alternative_ln = NA_REAL;
        configurations[i].full_relationship_ln = NA_REAL;
        configurations[i].relationship_comparison_available = false;
        configurations[i].full_proven_rejected = false;
        continue;
      }
      if(precomputed) {
        if(precomputed->size()!=configurations.size())engine_stop("focal memo reduction shape mismatch");
        reduction[i]=(*precomputed)[i];
      } else if(configurations[i].native_full_reduction) {
        reduction[i]=*configurations[i].native_full_reduction;
      } else {
      finalize_contextual_side_members(
        configurations[i].full_relationship_members);
      align_contextual_operational_display_reference(
        configurations[i].full_relationship_members, null_objective);
      reduction[i] = full_reduce_members(
        configurations[i].full_relationship_members, configurations[i].nparents,
        null_objective);
      }
      configurations[i].full_operational_ln = reduction[i].operational_ln;
      configurations[i].full_alternative_ln = reduction[i].alternative_ln;
      configurations[i].full_relationship_ln = reduction[i].relationship_ln;
      configurations[i].relationship_comparison_available =
        reduction[i].comparison_available;
      configurations[i].full_proven_rejected =
        configurations[i].nparents > 0 &&
        full_relationship_known_rejection(reduction[i], input.Tassign_ln);
    }
    if(numerical_out)*numerical_out=reduction;
    for (std::size_t i = 0; i < configurations.size(); ++i)
      configurations[i].full_relationship_gate_pass =
        full_relationship_gate_for_configuration(
          configurations, reduction, static_cast<int>(i),
          input.Tassign_ln);
  }
  int null_index = -1;
  for (std::size_t i = 0; i < configurations.size(); ++i)
    if (configurations[i].nparents == 0) { null_index = static_cast<int>(i); break; }
  if (null_index < 0) engine_stop("whole-pass configuration set lost its null row");
  const double null_objective = configurations[null_index].objective_ln;
  const auto comparable = [](double value) { return !std::isnan(value); };
  const auto row_contrast = [&](const ConfigurationScore& numerator,
                            const ConfigurationScore& denominator) {
    return numerator.raw_comparison_available && denominator.raw_comparison_available ?
      configuration_contrast(numerator, denominator) : NA_REAL;
  };
  for (ConfigurationScore& configuration : configurations) {
    if(input.joint_identity_objective&&input.null_joint_identity&&configuration.nparents==0){
      configuration.joint_identity=true;configuration.joint_identity_ln=0.;
      configuration.full_operational_ln=configuration.full_alternative_ln=configuration.full_relationship_ln=0.;
      configuration.relationship_comparison_available=true;configuration.full_relationship_gate_pass=true;
    }
    const bool relationship_available = configuration.relationship_comparison_available;
    const bool relationship_rejected = configuration.nparents > 0 &&
      (configuration.full_proven_rejected ||
       (input.joint_identity_objective &&
        (configuration.full_stage_proven_rejected||configuration.relationship_stage_proven_rejected)) ||
       (relationship_available && !configuration.full_relationship_gate_pass));
    bool support_rejected = configuration.raw_support_rejected;
    const bool incumbent_only = (configuration.dam == 0 || configuration.incumbent_dam) &&
      (configuration.sire == 0 || configuration.incumbent_sire) && configuration.nparents > 0;
    if (incumbent_only) support_rejected = false;
    bool support_pass = true;
    configuration.support_comparison_available = (input.joint_identity_objective&&input.null_joint_identity&&configuration.nparents==0) ||
      incumbent_only || configuration.raw_comparison_available;
    const int dam_index = configuration.nparents == 2 ?
      single_indices.dam_only(configuration.dam) : -1;
    const int sire_index = configuration.nparents == 2 ?
      single_indices.sire_only(configuration.sire) : -1;
    const auto consume_support = [&](double support) {
      if (!comparable(support)) configuration.support_comparison_available = false;
      else if (support <= input.Tassign_ln) support_rejected = true;
      support_pass = support_pass && comparable(support) && support > input.Tassign_ln;
    };
    auto* joint_live=pap_internal::native_demand::live;
    if(joint_live&&input.joint_support_gates&&configuration.hard_feasible&&
       configuration.status==0&&configuration.trio_me_pass){
      if(configuration.gate_deferred||
         (!(input.configuration_diagnostics||input.configuration_locus_diagnostics)&&
          (relationship_rejected||!relationship_available))){
        configuration.dam_support_ln=configuration.sire_support_ln=configuration.pair_support_ln=NA_REAL;
        configuration.support_comparison_available=false;support_pass=false;
      }else{
        const auto joint=pap_internal::joint_support_gate::evaluate(
          configuration.dam,configuration.sire,configuration.incumbent_dam,
          configuration.incumbent_sire,input.Tassign_ln,
          [&](int d,int s){return joint_live->joint_support(d,s);},input.support_available_baselines);
        configuration.dam_support_ln=joint.support[0];
        configuration.sire_support_ln=joint.support[1];
        configuration.pair_support_ln=joint.support[2];
        configuration.support_comparison_available=joint.available;
        support_pass=joint.pass;support_rejected=joint.rejected;
      }
    } else if (configuration.nparents == 1) {
      const double support = row_contrast(configuration, configurations[null_index]);
      configuration.dam_support_ln = configuration.dam > 0 ? support : R_NegInf;
      configuration.sire_support_ln = configuration.sire > 0 ? support : R_NegInf;
      if (!(configuration.incumbent_dam || configuration.incumbent_sire)) consume_support(support);
    } else if (configuration.nparents == 2) {
      configuration.dam_support_ln = sire_index < 0 ? NA_REAL :
        row_contrast(configuration, configurations[sire_index]);
      configuration.sire_support_ln = dam_index < 0 ? NA_REAL :
        row_contrast(configuration, configurations[dam_index]);
      bool baseline_available = configurations[null_index].raw_comparison_available &&
        dam_index >= 0 && sire_index >= 0 &&
        configurations[dam_index].raw_comparison_available &&
        configurations[sire_index].raw_comparison_available;
      configuration.pair_support_ln = baseline_available && configuration.raw_comparison_available ?
        log_contrast(configuration.objective_ln,
          std::max(null_objective, std::max(configurations[dam_index].objective_ln,
                                          configurations[sire_index].objective_ln))) : NA_REAL;
      if (!configuration.incumbent_dam) consume_support(configuration.dam_support_ln);
      if (!configuration.incumbent_sire) consume_support(configuration.sire_support_ln);
      if (!(configuration.incumbent_dam || configuration.incumbent_sire))
        consume_support(configuration.pair_support_ln);
    }
    // An AND gate can be disproved by an independently scored negative, even
    // when another question remains unavailable. That proof is distinct from
    // complete evaluation, and never inferred from an identity margin.
    configuration.independently_rejected = relationship_rejected || support_rejected;
    configuration.relationship_comparison_available = relationship_available &&
      configuration.support_comparison_available;
    configuration.relationship_pass = configuration.full_relationship_gate_pass &&
      configuration.support_comparison_available && support_pass;
  }
  // Identity competes against a different independently admissible row.
  // A pair already rejected by its relationship, support or hard gates
  // cannot block an assignable singleton with its raw objective. Keep
  // unknown relationship comparisons as competitors and exclude no row
  // merely because its own identity margin will be insufficient. Same-dam
  // and same-sire alternatives remain distinct rows; equal surviving top
  // values still give a zero margin. The top-two calculation stays linear.
  if(input.joint_identity_objective){
    auto* live=pap_internal::native_demand::live;
    for(auto& row:configurations){
      if(input.null_joint_identity&&row.nparents==0){row.joint_identity=true;row.joint_identity_ln=0.;continue;}
      row.joint_identity_ln=NA_REAL;
      if(!live||!identity_competitor(row))continue;
      // An already available exact score is a free bound, including a
      // counterfactual requested by another row's support gate. Unanswered
      // gates/window rows without that bound remain uncertainty competitors;
      // they cannot authorize a positive margin or request a contraction.
      if(live->cached_joint_identity(row.dam,row.sire,row.joint_identity_ln))continue;
      if(row.gate_deferred||(row.nparents>0&&!row.relationship_pass))continue;
      row.joint_identity_ln=live->joint_identity(row.dam,row.sire);
    }
  }
  const std::vector<double> identity_margins =
    identity_margin_values(configurations, true);
  for (std::size_t i = 0; i < configurations.size(); ++i) {
    ConfigurationScore& configuration = configurations[i];
    configuration.identity_margin_ln = identity_margins[i];
    // Keep hard/state validity separate from assignment confidence.  In
    // particular, a finite incumbent can remain a valid current state even
    // when a newly proposed configuration is not sufficiently separated from
    // its best competing identity.
    const bool base_valid = configuration.hard_feasible &&
      configuration.status == 0 &&
      configuration.trio_me_pass;
    configuration.hard_valid = base_valid;
    configuration.assignment_confident = base_valid &&
      ((configuration.joint_identity||configuration.full_identity) ? R_finite(decision_objective(configuration)) :
        (configuration.raw_comparison_available && R_finite(configuration.objective_ln))) &&
      configuration.relationship_comparison_available &&
      configuration.relationship_pass &&
      (configuration.nparents == 0 ||
       configuration.identity_margin_ln > input.identity_margin_ln);
    // Missing/undefined likelihoods are not biological invalidity or evidence
    // for removal. A proven independent negative still disqualifies retention.
    configuration.eligible_for_retention = configuration.nparents > 0 &&
      configuration.hard_valid && !configuration.independently_rejected;
    configuration.eligible = configuration.assignment_confident;
  }
  if (current_index < 0 ||
      current_index >= static_cast<int>(configurations.size()))
    engine_stop("whole-pass diagnostics received an invalid current configuration");
  const ConfigurationScore& current = configurations[current_index];
  const auto retains_current_sides = [&](const ConfigurationScore& proposal) {
    return (current.dam == 0 || proposal.dam == current.dam) &&
      (current.sire == 0 || proposal.sire == current.sire);
  };
  const bool current_hard_invalid = !current.hard_valid;
  const bool current_relationship_rejected = current.independently_rejected;
  for (ConfigurationScore& configuration : configurations) {
    const bool changed = configuration.dam != current.dam ||
      configuration.sire != current.sire;
    configuration.eligible_for_addition =
      configuration.assignment_confident && changed &&
      configuration.nparents == current.nparents + 1 &&
      retains_current_sides(configuration);
    configuration.eligible_for_replacement =
      configuration.assignment_confident && changed && current.nparents > 0 &&
      configuration.nparents == current.nparents;
    configuration.eligible_for_removal =
      configuration.assignment_confident && changed &&
      configuration.nparents < current.nparents &&
      (current_hard_invalid || current_relationship_rejected);
  }
}

int choose_configuration(
    const SearchInputs& input, int current_index,
    const std::vector<ConfigurationScore>& configurations,
    bool& near_tie, int& rejection_reason) {

  std::vector<int> eligible;
  for (std::size_t i = 0; i < configurations.size(); ++i)
    if (configurations[i].eligible &&
        (configurations[i].nparents >= configurations[current_index].nparents ||
         configurations[i].eligible_for_removal))
      eligible.push_back(static_cast<int>(i));
  // A current, supported parent is a valid state even when its identity
  // margin is insufficient for a fresh assignment.  Keep it in the chooser's
  // candidate set so ambiguity on the opposite side cannot fall through to
  // null.
  const bool current_retention_eligible =
    configurations[current_index].eligible_for_retention;
  if (current_retention_eligible &&
      std::find(eligible.begin(), eligible.end(), current_index) == eligible.end())
    eligible.push_back(current_index);
  bool has_relationship = false;
  for (int index : eligible)
    if (configurations[index].nparents > 0) { has_relationship = true; break; }
  if (has_relationship) {
    std::vector<int> filtered;
    for (int index : eligible)
      if (configurations[index].nparents > 0) filtered.push_back(index);
    eligible.swap(filtered);
  }
  if (eligible.empty()) {
    rejection_reason = 3;
    return current_index;
  }
  auto better = [&](int first, int second) {
    const double left = decision_objective(configurations[first]);
    const double right = decision_objective(configurations[second]);
    if (!objective_tie(left, right, input.tie_tolerance))
      return left > right;
    if (first == current_index && second != current_index) return true;
    if (second == current_index && first != current_index) return false;
    if (configurations[first].nparents != configurations[second].nparents)
      return configurations[first].nparents < configurations[second].nparents;
    const int first_dam = id_rank(input, configurations[first].dam);
    const int second_dam = id_rank(input, configurations[second].dam);
    if (first_dam != second_dam) return first_dam < second_dam;
    const int first_sire = id_rank(input, configurations[first].sire);
    const int second_sire = id_rank(input, configurations[second].sire);
    return first_sire < second_sire;
  };
  int chosen = eligible.front();
  for (std::size_t i = 1; i < eligible.size(); ++i)
    if (better(eligible[i], chosen)) chosen = eligible[i];
  near_tie = false;
  for (int index : eligible) if (index != chosen) {
    const bool tied=objective_tie(decision_objective(configurations[index]),
                                        decision_objective(configurations[chosen]),
                                        input.tie_tolerance);
    if(!tied)continue;
    near_tie = true;
    break;
  }
  const bool current_eligible = current_retention_eligible &&
    !(configurations[current_index].nparents == 0 && has_relationship);
  if (chosen == current_index) {
    rejection_reason = 1;
  } else if (current_eligible && configuration_contrast(configurations[chosen],
      configurations[current_index],true) <= input.tol) {
    rejection_reason = 2;
  } else {
    rejection_reason = 0;
  }
  return chosen;
}

}
