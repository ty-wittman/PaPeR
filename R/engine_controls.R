# pedigree likelihood execution defaults. Explicit internal inputs take precedence.
.engine_profile <- function() list(
  pedigree_likelihood = TRUE, pedigree_hypothesis_suite = "S7", likelihood_evidence = "ancestral",
  likelihood_native_setup = TRUE,
  likelihood_total_work_limit = TRUE, likelihood_direct_suite = TRUE,
  likelihood_catalog_limit_key = TRUE,
  likelihood_registry_across_passes = TRUE, likelihood_actual_width_domain = TRUE,
  likelihood_light_bookkeeping = TRUE, likelihood_bounded_lane_tables = TRUE,
  likelihood_profile_on_demand = TRUE,

  complete_identity_competitors = TRUE,
  completion_rank_only = TRUE,
  prescreen_identity = TRUE,
  likelihood_skip_unused_raw = TRUE, likelihood_confidence_raw_elision = TRUE,
  likelihood_skip_union_roots = TRUE,
  likelihood_identity_visit_facts = TRUE,
  likelihood_ancestral_message_reuse = TRUE,
  likelihood_focal_role_cavity = TRUE,
  likelihood_check_pedigree_joint = TRUE,
  likelihood_domain_order = TRUE,
  likelihood_contraction_log_retry = TRUE,
  support_available_baselines = TRUE,
  null_joint_identity = TRUE,
  likelihood_current_confidence_records = TRUE,
  likelihood_confidence_dependency_records = TRUE,
  likelihood_batch_loci = TRUE, likelihood_hypothesis_batches = TRUE,
  likelihood_parallel_messages = FALSE, likelihood_refresh = 1L,
  complete_recheck_scope = "local",
  local_conditioning = TRUE,
  local_core_only = TRUE, local_expansion = 0L,
  message_scoring = FALSE, stored_message_scoring = FALSE)

# Accepted pedigree likelihood execution defaults. Explicit controls always take precedence.
# Native controls require the matching installed DLL.
.engine_defaults <- function(inputs) {
  ancestral <- isTRUE(inputs$pedigree_likelihood) &&
    (is.null(inputs$likelihood_evidence) || identical(inputs$likelihood_evidence, "ancestral"))
  prescreen <- ancestral && !isTRUE(inputs$pair_rescue_shadow) &&
    !isTRUE(inputs$pair_rescue_calibration)
  defaults <- list(

    pedigree_trust_prepared = TRUE,
    pedigree_share_finished_topology = TRUE,
    pedigree_metrics_direct = TRUE,
    likelihood_evidence = if (isTRUE(inputs$pedigree_likelihood)) "ancestral" else "full",
    likelihood_ancestral_reference = FALSE,
    likelihood_actual_width_domain = FALSE,
    likelihood_bounded_program_registry = TRUE,
    likelihood_program_registry_mib = 16L,
    likelihood_keep_ancestral_unit_pairs = FALSE,
    likelihood_barren_pruning = TRUE,
    likelihood_local_family_pairs = FALSE,
    likelihood_latent_brood_exclusion = FALSE,
    incumbent_full_scoring = TRUE,
    identity_objective = "joint",
    joint_support_gates = TRUE,
    likelihood_identity_union = FALSE,
    pair_rescue = prescreen,
    pair_rescue_support_only = ancestral,
    pair_rescue_margin_ln = if (ancestral) log(100) else 0,
    pair_rescue_top_k = if (ancestral) 10L else 50L,
    pedigree_relative_life_cache = TRUE,
    pedigree_scorer_setup_reuse = TRUE,
    pedigree_engine_descriptor = TRUE,
    pedigree_skip_engine_admissions = TRUE,
    pedigree_native_support = TRUE,
    pedigree_native_catalog = TRUE,
    pedigree_topology_reuse = TRUE,
    pedigree_local_dirty_k = 0L)
  if (isTRUE(inputs$pedigree_likelihood)) {
    profile <- .engine_profile()
    if (!is.null(inputs$pedigree_hypothesis_s7) && is.null(inputs$pedigree_hypothesis_suite))
      profile$pedigree_hypothesis_suite <- NULL
    defaults <- utils::modifyList(defaults, profile)
  }
  missing <- setdiff(names(defaults), names(inputs))
  inputs[missing] <- defaults[missing]
  inputs$likelihood_evidence <- match.arg(inputs$likelihood_evidence, c("full", "ancestral"))
  inputs$identity_objective <- match.arg(inputs$identity_objective,
    c("joint", "legacy_side_u"))
  if (inputs$likelihood_evidence == "ancestral" && !isTRUE(inputs$pedigree_likelihood))
    stop("ancestral evidence requires mped2 = TRUE", call. = FALSE)
  inputs
}

# Deterministic cache ownership only; this never consumes the scientific RNG.
.next_registry_generation <- local({
  counter <- 0
  function() {
    counter <<- counter + 1
    if (counter >= 2^52) counter <<- 1
    counter
  }
})
