#include "execution_config.h"
#include "hypothesis_suite.h"
using namespace Rcpp;
namespace pap_internal {
ExecutionConfig::ExecutionConfig(const List& inputs) {
auto flag=[&](const char* n){return inputs.containsElementNamed(n)&&Rcpp::as<bool>(inputs[n]);};
 pedigree_likelihood=flag("pedigree_likelihood");
 likelihood_implicit_transmission=flag("likelihood_implicit_transmission");
 likelihood_factorized_one=flag("likelihood_factorized_one");
 likelihood_high_allele_programs=flag("likelihood_high_allele_programs");
{auto read_limit=[&](const char* name,int& target){if(inputs.containsElementNamed(name)){double x=as<double>(inputs[name]);if(!std::isfinite(x)||x!=std::floor(x)||x<0||x>5000)throw std::invalid_argument(std::string(name)+" must be a supported integer limit");target=int(x);}};
 read_limit("likelihood_memory_mib",likelihood_limits.mib);read_limit("likelihood_work_millions",likelihood_limits.work_millions);
 {auto read_order=[&](const char* name,int& target){if(inputs.containsElementNamed(name)){double x=as<double>(inputs[name]);if(!std::isfinite(x)||x!=std::floor(x)||x<1||x>1000000)throw std::invalid_argument(std::string(name)+" must be a bounded integer");target=int(x);}};
 read_order("likelihood_resource_min_alleles",likelihood_limits.resource_min_alleles);}
 likelihood_limits.s37_refusal_repair=inputs.containsElementNamed("pedigree_hypothesis_suite")&&
   pap_internal::hypothesis_suite::parse(as<std::string>(inputs["pedigree_hypothesis_suite"]))==pap_internal::hypothesis_suite::Suite::S37;
 likelihood_local_family_pairs=flag("likelihood_local_family_pairs");
 likelihood_latent_brood_exclusion=flag("likelihood_latent_brood_exclusion");
 const std::string evidence=inputs.containsElementNamed("likelihood_evidence")?as<std::string>(inputs["likelihood_evidence"]):(pedigree_likelihood?"ancestral":"full");
 if(evidence!="full"&&evidence!="ancestral")throw std::invalid_argument("mped2_evidence must be full or ancestral");
 likelihood_time_evaluations=!flag("likelihood_profile_on_demand")||flag("pedigree_phase_profile")||flag("likelihood_plan_profile")||flag("pedigree_engine_descriptor_profile");likelihood_light_bookkeeping=flag("likelihood_light_bookkeeping");likelihood_ancestral=evidence=="ancestral";likelihood_ancestral_message_reuse=flag("likelihood_ancestral_message_reuse");likelihood_ancestral_reference=flag("likelihood_ancestral_reference");
 likelihood_limits.ancestral_unit_pairs=!flag("likelihood_keep_ancestral_unit_pairs");
 likelihood_limits.barren_pruning=!inputs.containsElementNamed("likelihood_barren_pruning")||flag("likelihood_barren_pruning");
 if(likelihood_ancestral&&(!pedigree_likelihood||likelihood_local_family_pairs||likelihood_latent_brood_exclusion))
   throw std::invalid_argument("ancestral evidence requires M_ped2 with both earlier family-evidence experiments OFF");
 likelihood_limits.total_work_limit=flag("likelihood_total_work_limit");likelihood_limits.domain_order=!inputs.containsElementNamed("likelihood_domain_order")||flag("likelihood_domain_order");likelihood_limits.contraction_log_retry=!inputs.containsElementNamed("likelihood_contraction_log_retry")||flag("likelihood_contraction_log_retry");
 likelihood_limits.resource_planning=flag("likelihood_resource_planning");likelihood_limits.live_arena=flag("likelihood_live_arena");likelihood_limits.shared_piece_scratch=flag("likelihood_shared_piece_scratch");likelihood_limits.validate();
 if((likelihood_limits.mib!=64||likelihood_limits.work_millions!=400)&&!likelihood_high_allele_programs)throw std::invalid_argument("expanded M_ped2 limits require mped2_high_allele_programs");}
likelihood_plan_profile=flag("likelihood_plan_profile");likelihood_no_pieces=flag("likelihood_no_pieces");likelihood_raw_separate=flag("likelihood_raw_separate");likelihood_no_dedup=flag("likelihood_no_dedup");likelihood_no_lanes=flag("likelihood_no_lanes");likelihood_narrow_lanes=flag("likelihood_narrow_lanes");if(inputs.containsElementNamed("likelihood_max_touch")){likelihood_max_touch=as<int>(inputs["likelihood_max_touch"]);if(likelihood_max_touch<0||likelihood_max_touch>4)throw std::invalid_argument("mped2_max_touch must be 0..4");}likelihood_dense_gametes=flag("likelihood_dense_gametes");if(inputs.containsElementNamed("likelihood_refresh")){likelihood_refresh=as<int>(inputs["likelihood_refresh"]);if(likelihood_refresh<0||likelihood_refresh>3)throw std::invalid_argument("mped2_refresh must be 0 (visit), 1 (graph revision), 2 (round) or 3 (in-place update)");}
 if(inputs.containsElementNamed("relationship_common_pair_fast")) evaluation.relationship_common_pair_fast=as<bool>(inputs["relationship_common_pair_fast"]);
 if(inputs.containsElementNamed("relationship_common_scaled_fast")) evaluation.relationship_common_scaled_fast=as<bool>(inputs["relationship_common_scaled_fast"]);
 if(inputs.containsElementNamed("private_test_message_comparison_decline")) evaluation.private_test_message_comparison_decline=as<bool>(inputs["private_test_message_comparison_decline"]);
}
}
