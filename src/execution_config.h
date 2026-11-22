#pragma once
#include <Rcpp.h>
#include <memory>
#include "resource_limits.h"
namespace pap_internal {
// Fixed run options are parsed once. Visit-specific arithmetic overrides are separate.
struct EvaluationOptions {
 bool relationship_common_pair_fast=true;
 bool relationship_common_scaled_fast=true;
 bool private_test_message_comparison_decline=false;
};
struct ExecutionConfig {
 pap_internal::pedigree_likelihood::ResourceLimits likelihood_limits;
 bool likelihood_ancestral=false,likelihood_ancestral_reference=false;
 bool likelihood_local_family_pairs=false,likelihood_latent_brood_exclusion=false;
 bool likelihood_ancestral_message_reuse=false,likelihood_light_bookkeeping=false,likelihood_time_evaluations=true;
 bool likelihood_high_allele_programs=false,likelihood_factorized_one=false;
 bool likelihood_implicit_transmission=false,pedigree_likelihood=false,likelihood_no_dedup=false;
 bool likelihood_no_lanes=false,likelihood_narrow_lanes=false,likelihood_dense_gametes=false;
 bool likelihood_plan_profile=false,likelihood_no_pieces=false,likelihood_raw_separate=false;
 int likelihood_refresh=1,likelihood_max_touch=2;
 EvaluationOptions evaluation;
 explicit ExecutionConfig(const Rcpp::List& inputs);
};
}
