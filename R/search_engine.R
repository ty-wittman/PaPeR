
.engine_inputs <- function(
    snapshot, state, pool, ids, universe, Tassign_ln, margin_ln, Tfilter_ln,
    tol, max_rounds, shortlist_settings, run_check_pedigree,
    relationship_mode = c("simp", "full"), diagnostics = list(),
    likelihood_decrease_tolerance = 1e-8) {
  relationship_mode <- match.arg(relationship_mode)
  if (!identical(shortlist_settings$mode, "sequoia50"))
    stop("native Stage-1 requires shortlist_mode='sequoia50'", call. = FALSE)
  if (is.null(diagnostics)) diagnostics <- list()
  if (!is.list(diagnostics) || (length(diagnostics) &&
      (is.null(names(diagnostics)) || any(!nzchar(names(diagnostics))))))
    stop("native Stage-1 diagnostics must be a named list", call. = FALSE)
  supported <- c(
    "configuration_diagnostics", "complete_recheck", "complete_recheck_scope",
    "native_input_callback", "native_state_callback", "configuration_locus_diagnostics",
    "configuration_locus_focal_position", "adaptive_selective_publisher", "adaptive_selective_threshold",
    "search_loop_engine", "benchmark_opt5_diag_mode", "phase_timing",
    "lazy_independent_gates", "incumbent_full_scoring", "message_scoring",
    "stored_message_scoring", "raw_scaled_fast", "defer_intermediate_publication",
    "bounded_kernel_profile", "focal_start_position", "focal_end_position",
    "progress_locus_interval", "progress_heartbeat_seconds", "progress_locus_heartbeats",
    "progress_callback", "visit_progress_callback")
  unknown <- setdiff(names(diagnostics), supported)
  if (length(unknown))
    stop("unsupported native Stage-1 diagnostic field(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  input <- .visit_inputs(
    snapshot, state, pool, ids, universe, Tassign_ln, margin_ln, tol,
    shortlist_settings$mode)
  logical_flag <- function(name, default = FALSE) {
    value <- diagnostics[[name]] %||% default
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(name, " must be one non-missing logical value", call. = FALSE)
    isTRUE(value)
  }
  engine <- diagnostics$search_loop_engine %||% "bounded_messages"
  if (!identical(engine, "bounded_messages"))
    stop("native Stage-1 accepts bounded_messages only", call. = FALSE)
  input$search_loop_engine <- engine
  # This is deliberately a relationship extension flag, not an engine selector.
  # Both modes retain the identical fixed-50 candidate arrays, frozen message
  # publication, native whole pass, chooser, moves, and convergence lifecycle.
  input$relationship_mode <- if (identical(relationship_mode, "full")) 1L else 0L
  # The native payload has one shared ordinary relationship universe.  Full
  # adds its reviewed compound rows to that same reduction; it does not own a
  # separate ordinary-parentage scorer or candidate path.
  input$full_relationship_registry <- .full_native_registry_payload()
  input$configuration_diagnostics <- logical_flag("configuration_diagnostics")
  # Complete-pair rechecks are an additional native lifecycle pass, not a
  # second scorer. The ordinary pass retains Sequoia's complete-pair skip;
  # the engine alone opens selected stale complete focals with recheck_only.
  complete_recheck <- as.character(diagnostics$complete_recheck %||%
    "first_dirty_wave_convergence")
  complete_recheck_scope <- as.character(
    diagnostics$complete_recheck_scope %||% "all_complete")
  if (length(complete_recheck) != 1L || !complete_recheck %in%
      c("none", "end_round", "convergence", "first_dirty_wave_convergence"))
    stop("complete_recheck must be one of none, end_round, convergence, or first_dirty_wave_convergence",
         call. = FALSE)
  if (length(complete_recheck_scope) != 1L ||
      !complete_recheck_scope %in% c("all_complete", "local"))
    stop("complete_recheck_scope must be all_complete or local",
         call. = FALSE)
  input$complete_recheck <- complete_recheck
  input$complete_recheck_scope <- complete_recheck_scope
  input$recheck_only <- FALSE
  input$recheck_focals <- integer()
  input$configuration_locus_diagnostics <- logical_flag(
    "configuration_locus_diagnostics")
  input$configuration_locus_focal_position <- as.integer(
    diagnostics$configuration_locus_focal_position %||% 0L)
  # Choose full or dirty-region message publication from the changed evidence.
  input$adaptive_selective_publisher <- logical_flag(
    "adaptive_selective_publisher", TRUE)
  threshold <- diagnostics$adaptive_selective_threshold %||% 0.5
  if (length(threshold) != 1L || is.na(threshold) ||
      !is.finite(threshold) || threshold < 0 || threshold > 1)
    stop("adaptive_selective_threshold must be one finite number in [0, 1]",
         call. = FALSE)
  input$adaptive_selective_threshold <- as.numeric(threshold)
  input$benchmark_opt5_diag_mode <- logical_flag("benchmark_opt5_diag_mode")
  input$phase_timing <- logical_flag("phase_timing")
  input$lazy_independent_gates <- logical_flag("lazy_independent_gates", TRUE)
  # Leave an unspecified policy absent so only pedigree likelihood supplies its accepted default.
  if ("incumbent_full_scoring" %in% names(diagnostics))
    input$incumbent_full_scoring <- logical_flag("incumbent_full_scoring")
  input$message_scoring <- logical_flag("message_scoring", FALSE)
  input$stored_message_scoring <- logical_flag("stored_message_scoring", FALSE)
  input$message_scoring <- input$message_scoring || input$stored_message_scoring
  input$raw_scaled_fast <- logical_flag("raw_scaled_fast", FALSE)
  input$defer_intermediate_publication <- logical_flag(
    "defer_intermediate_publication", FALSE)
  input$bounded_kernel_profile <- logical_flag("bounded_kernel_profile")
  input$focal_start_position <- as.integer(
    diagnostics$focal_start_position %||% 1L)
  input$focal_end_position <- as.integer(diagnostics$focal_end_position %||% 0L)
  input$progress_locus_interval <- as.integer(
    diagnostics$progress_locus_interval %||% 5L)
  input$progress_heartbeat_seconds <- as.numeric(
    diagnostics$progress_heartbeat_seconds %||% 15)
  input$progress_locus_heartbeats <- logical_flag(
    "progress_locus_heartbeats", TRUE)
  input$progress_callback <- diagnostics$progress_callback %||% NULL
  input$visit_progress_callback <- diagnostics$visit_progress_callback %||% NULL
  expected_year <- vapply(state$birth_year_masses, function(x)
    if (is.null(x)) Inf else sum(as.numeric(names(x)) * x), numeric(1))
  input$expected_birth_year <- as.numeric(expected_year[ids])
  input$id_labels <- as.character(ids)
  input$Tfilter_ln <- as.numeric(Tfilter_ln)
  input$max_rounds <- as.integer(max_rounds)
  input$run_check_pedigree <- isTRUE(run_check_pedigree)
  input$likelihood_decrease_tolerance <- as.numeric(
    likelihood_decrease_tolerance)
  input
}

.hydrate_messages <- function(template, state, result) {
  ids <- .parentage_state_ids(state)
  values <- as.list.environment(template, all.names = TRUE)
  matrices <- result$final_messages_ln
  if (!is.list(matrices) || length(matrices) != length(values$messages))
    stop("native Stage-1 final message payload is malformed", call. = FALSE)
  values$messages <- lapply(matrices, function(x) {
    if (!is.matrix(x) || ncol(x) != length(ids))
      stop("native Stage-1 final message matrix is malformed", call. = FALSE)
    lapply(seq_len(ncol(x)), function(i) x[, i])
  })
  names(values$messages) <- names(template$messages)
  likelihood <- result$final_locus_likelihood_ln
  if (!is.matrix(likelihood) ||
      !identical(dim(likelihood), c(length(values$loci), length(ids))))
    stop("native Stage-1 locus-likelihood payload is malformed", call. = FALSE)
  values$individual_likelihood_ln <- lapply(seq_len(nrow(likelihood)),
    function(i) setNames(as.numeric(likelihood[i, ]), ids))
  names(values$individual_likelihood_ln) <- names(template$messages)
  values$message_revision <- as.integer(state$message_revision)
  values$graph_revision <- as.integer(state$graph_revision)
  values$message_signature <- paste0(
    "message-revision=", state$message_revision)
  values$context_token <- paste0(values$context_identity$schema_version, ":",
    length(ids), ":", length(values$loci), ":", state$message_revision)
  values$parents_at_construction <- state$parents
  values$topological_order <- .message_topological_order(state$parents)
  routing <- .message_snapshot_families(ids, state$parents)
  values$families <- routing$families
  values$families_by_parent <- routing$by_parent
  values$graph_provenance <- list(
    message_revision = as.integer(state$message_revision),
    graph_revision = as.integer(state$graph_revision),
    parents_at_construction = state$parents,
    purpose = values$purpose)
  values$diagnostic_metadata$publication <-
    "native Stage-1 final result hydration"
  snapshot <- list2env(values, parent = emptyenv())
  class(snapshot) <- class(template)
  .message_snapshot_lock(snapshot)
  snapshot
}

.accepted_history <- function(result, ids, initial_revision) {
  out <- list(); at <- 0L
  checks <- result$check_pedigree
  revision <- as.integer(initial_revision)
  if (nrow(checks)) for (i in seq_len(nrow(checks))) {
    if (!identical(as.character(checks$action[i]), "drop_complete_pair")) next
    revision <- revision + 1L
    at <- at + 1L
    out[[at]] <- list(
      graph_revision = revision,
      message_revision = as.integer(result$message_revision_history$message_revision[1L]),
      child = ids[checks$child[i]], move_type = "replace_pair",
      old = c(dam = if (checks$old_dam[i] > 0L) ids[checks$old_dam[i]] else NA,
              sire = if (checks$old_sire[i] > 0L) ids[checks$old_sire[i]] else NA),
      proposed = c(dam = NA_character_, sire = NA_character_),
      affected_nodes = ids, edges_added = character(),
      edges_removed = character(),
      score = list(total = -checks$pair_support_ln[i]),
      label = "native-CheckPedigree")
  }
  for (round in seq_along(result$accepted_move_blocks)) {
    block <- result$accepted_move_blocks[[round]]
    if (!nrow(block)) next
    for (i in seq_len(nrow(block))) {
      at <- at + 1L
      decode <- function(index) if (index > 0L) ids[index] else NA_character_
      out[[at]] <- list(
        graph_revision = as.integer(block$graph_revision_after[i]),
        message_revision = as.integer(
          result$round_summary$message_revision_start[round]),
        child = ids[block$child[i]], move_type = "native_parentage_move",
        old = c(dam = decode(block$old_dam[i]),
                sire = decode(block$old_sire[i])),
        proposed = c(dam = decode(block$new_dam[i]),
                     sire = decode(block$new_sire[i])),
        # Hydration starts a conservative final cache epoch.  All nodes are
        # dependencies, so no stale pre-engine cache can survive.
        affected_nodes = ids, edges_added = character(),
        edges_removed = character(),
        score = list(genetic = block$objective_delta_ln[i], age = 0,
                     total = block$objective_delta_ln[i]),
        label = "native-stage1-engine")
    }
  }
  out
}

# Convert native move blocks into ID-labeled audit rows.
.named_moves <- function(moves, ids, round) {
  nonempty <- is.data.frame(moves) && nrow(moves) > 0L
  if (!nonempty) moves <- data.frame()
  decode <- function(x) {
    out <- rep(NA_character_, length(x)); keep <- x > 0L
    out[keep] <- ids[x[keep]]; out
  }
  column <- function(value, empty) list(value = value, empty = empty)
  # One projection definition also preserves the established empty-row schema.
  fields <- list(
    round = column(as.integer(round), integer()),
    child = column(ids[moves$child], character()),
    previous_dam = column(decode(moves$old_dam), character()),
    previous_sire = column(decode(moves$old_sire), character()),
    proposed_dam = column(decode(moves$new_dam), character()),
    proposed_sire = column(decode(moves$new_sire), character()),
    objective_delta_ln = column(moves$objective_delta_ln, numeric()),
    dam_support_ln = column(moves$dam_support_ln, NULL),
    sire_support_ln = column(moves$sire_support_ln, NULL),
    pair_support_ln = column(moves$pair_support_ln, NULL),
    identity_margin_ln = column(moves$identity_margin_ln, NULL),
    graph_revision_before = column(moves$graph_revision_before, integer()),
    graph_revision_after = column(moves$graph_revision_after, integer()))
  values <- lapply(fields, `[[`, if (nonempty) "value" else "empty")
  if (!nonempty) values <- Filter(Negate(is.null), values)
  do.call(data.frame, c(values, list(stringsAsFactors = FALSE)))
}

.hydrate_state <- function(
    state, result, pool, snapshot) {
  ids <- .parentage_state_ids(state)
  initial_revision <- state$graph_revision
  parents <- result$final_parents
  if (!is.matrix(parents) || !identical(dim(parents), c(length(ids), 2L)))
    stop("native Stage-1 final parents are malformed", call. = FALSE)
  storage.mode(parents) <- "integer"
  dimnames(parents) <- list(ids, c("dam", "sire"))
  state$parents <- parents
  indexes <- .derive_parentage_indexes(ids, parents, state$sex)
  expected_role <- unname(as.integer(indexes$inferred_role))
  if (!identical(expected_role, as.integer(result$final_inferred_role)))
    stop("native Stage-1 inferred roles disagree with final parents",
         call. = FALSE)
  state$reverse_offspring <- indexes$reverse_offspring
  state$inferred_role <- indexes$inferred_role
  state$role_counts <- indexes$role_counts
  state$graph_revision <- as.integer(result$final_graph_revision)
  state$message_revision <- as.integer(result$final_message_revision)
  # Native search revisions are authoritative.  Final hydration initializes
  # conservative coherent revision stamps rather than reconstructing mutation
  # side effects through the accepted-move sequence.
  state$node_revision[] <- state$graph_revision
  state$component_revision[] <- state$graph_revision
  state$node_likelihood_ln <- setNames(
    as.numeric(result$node_likelihood_ln), ids)
  state$total_likelihood_ln <- as.numeric(result$total_likelihood_ln)
  state$candidate_pool <- pool
  state$accepted_move_history <- .accepted_history(
    result, ids, initial_revision)
  state$round_history <- lapply(seq_len(nrow(result$round_summary)), function(i)
    result$round_summary[i, , drop = FALSE])
  state$focal_order <- ids[as.integer(result$final_focal_order)]
  state$state_overlay_cache <- new.env(parent = emptyenv())
  state$configuration_score_cache <- new.env(parent = emptyenv())
  state$message_score_cache <- new.env(parent = emptyenv())
  state$message_snapshot <- .hydrate_messages(
    snapshot, state, result)
  state$current_messages <- state$message_snapshot$messages
  .validate_state_graph(state)
  state
}

.materialize_candidates <- function(
    pool, state, input, last_result, ids, shortlist_settings) {
  diagnostics <- .message_whole_pass_diagnostics(
    last_result, ids, input, pool, state, shortlist_settings)
  for (o in seq_along(ids)) for (side in c("dam", "sire")) {
    key <- .candidate_pool_key(ids[o], side)
    flat <- if (side == "dam") input$dam_ids else input$sire_ids
    offsets <- if (side == "dam") input$dam_offsets else input$sire_offsets
    active <- if (offsets[o + 1L] > offsets[o])
      flat[seq.int(offsets[o] + 1L, offsets[o + 1L])] else integer()
    pool$active_rows[[key]] <- which(
      pool$rows$offspring == ids[o] & pool$rows$side == side &
        pool$rows$candidate %in% ids[active])
  }
  pool$shortlist_diagnostics <- list(
    schema_version = "candidate-diagnostics-v2", mode = "sequoia50",
    settings = shortlist_settings, per_offspring_side = diagnostics)
  screen_key <- .candidate_pool_key(pool$diagnostics$offspring,
                                    pool$diagnostics$side)
  shortlist_key <- .candidate_pool_key(diagnostics$offspring,
                                       diagnostics$side)
  extra <- diagnostics[setdiff(names(diagnostics), c("offspring", "side"))]
  extra <- extra[match(screen_key, shortlist_key), , drop = FALSE]
  pool$diagnostics <- cbind(pool$diagnostics, extra)
  list(pool = pool, diagnostics = diagnostics)
}

.prepare_search <- function(
    state, pool, universe, snapshot, Tassign_ln, margin_ln, Tfilter_ln,
    tol, max_rounds, shortlist_settings, likelihood_controls,
    run_check_pedigree = TRUE, relationship_mode = c("simp", "full"),
    pass_diagnostics = list()) {
  relationship_mode <- match.arg(relationship_mode)
  ids <- .parentage_state_ids(state)
  input_started <- proc.time()[["elapsed"]]
  input <- .engine_inputs(
    snapshot, state, pool, ids, universe, Tassign_ln, margin_ln, Tfilter_ln,
    tol, max_rounds, shortlist_settings, run_check_pedigree,
    relationship_mode = relationship_mode, diagnostics = pass_diagnostics)
  # The generic input builder supplies older defaults. Replace its implicit
  # values with the pedigree likelihood profile; explicit diagnostics and controls win.
  profile <- .engine_profile()
  fields <- setdiff(names(profile), union(names(pass_diagnostics), names(likelihood_controls)))
  input[fields] <- profile[fields]
  input[names(likelihood_controls)] <- likelihood_controls
  input$pedigree_likelihood <- TRUE
  input <- .engine_defaults(input)
  input_seconds <- proc.time()[["elapsed"]] - input_started
  compact <- .message_compact_snapshot_inputs(snapshot)
  list(input = input, compact = compact, input_seconds = input_seconds)
}

.run_search <- function(
    ctx, state, pool, universe, snapshot, Tassign_ln, margin_ln, Tfilter_ln,
    tol, max_rounds, shortlist_settings, likelihood_controls,
    run_check_pedigree = TRUE, relationship_mode = c("simp", "full"),
    pass_diagnostics = list()) {
  relationship_mode <- match.arg(relationship_mode)
  ids <- .parentage_state_ids(state)
  initial <- .parentage_state_data_frame(state)[, c("id", "dam", "sire")]
  setup <- .prepare_search(
    state, pool, universe, snapshot, Tassign_ln, margin_ln, Tfilter_ln,
    tol, max_rounds, shortlist_settings, likelihood_controls, run_check_pedigree,
    relationship_mode, pass_diagnostics)
  input <- setup$input; compact <- setup$compact
  input_seconds <- setup$input_seconds
  started <- proc.time()[["elapsed"]]
  state_callback <- pass_diagnostics$native_state_callback
  if (!is.null(state_callback))
    state_callback(list(event = "start", compact = compact, input = input))
  result <- run_parentage_cpp(snapshot$compiled, compact, input)
  native_seconds <- proc.time()[["elapsed"]] - started
  final_active <- result$final_active_input
  for (name in names(final_active)) input[[name]] <- final_active[[name]]
  state <- .hydrate_state(
    state, result, pool, snapshot)
  materialized <- .materialize_candidates(
    pool, state, input, result$last_pass, ids, shortlist_settings)
  pool <- materialized$pool
  state$candidate_pool <- pool

  move_blocks <- result$accepted_move_blocks
  accepted_raw <- if (length(move_blocks)) do.call(rbind, move_blocks) else
    data.frame()
  recheck_move_blocks <- result$accepted_recheck_move_blocks %||% list()
  accepted_recheck_raw <- if (length(recheck_move_blocks))
    do.call(rbind, recheck_move_blocks) else data.frame()
  changes <- if (length(move_blocks)) do.call(rbind, lapply(
    seq_along(move_blocks), function(round)
      .named_moves(move_blocks[[round]], ids, round))) else
        .named_moves(data.frame(), ids, integer())
  recheck_changes <- if (length(recheck_move_blocks)) do.call(rbind, lapply(
    seq_along(recheck_move_blocks), function(wave) {
      phase <- result$recheck_phases
      round <- if (is.data.frame(phase) && nrow(phase) >= wave)
        as.integer(phase$round[wave]) else as.integer(wave)
      .named_moves(recheck_move_blocks[[wave]], ids, round)
    })) else .named_moves(data.frame(), ids, integer())
  if (nrow(recheck_changes)) changes <- rbind(changes, recheck_changes)
  check <- result$check_pedigree
  decode <- function(x) {
    out <- rep(NA_character_, length(x)); keep <- x > 0L
    out[keep] <- ids[x[keep]]; out
  }
  if (nrow(check)) {
    check$child <- ids[check$child]
    check$old_dam <- decode(check$old_dam)
    check$old_sire <- decode(check$old_sire)
    check$active_pair_check_pass <- as.logical(check$active_pair_check_pass)
    check$individual_side_check_active <-
      as.logical(check$individual_side_check_active)
  }
  rounds <- result$round_summary
  rounds$publication_performed <- as.logical(rounds$publication_performed)
  rounds$unexplained_likelihood_decrease <-
    as.logical(rounds$unexplained_likelihood_decrease)
  for (field in c("total_likelihood_current", "publication_deferred"))
    if (field %in% names(rounds)) rounds[[field]] <- as.logical(rounds[[field]])
  rounds$source_converged <- as.logical(rounds$source_converged)
  if ("recheck_converged" %in% names(rounds))
    rounds$recheck_converged <- as.logical(rounds$recheck_converged)
  focal_blocks <- result$focal_summary_blocks
  if (nrow(rounds)) for (round in seq_len(nrow(rounds))) {
    focal <- focal_blocks[[round]]
    scored <- focal$complete_pair_skipped == 0L
    rounds$complete_pairs_at_start[round] <- sum(!scored)
    rounds$complete_pairs_skipped[round] <- sum(!scored)
    rounds$incomplete_individuals_considered[round] <- sum(scored)
    rounds$focals_scored[round] <- sum(scored)
    rounds$active_dam_candidates[round] <-
      sum(focal$dam_candidate_count[scored])
    rounds$active_sire_candidates[round] <-
      sum(focal$sire_candidate_count[scored])
    rounds$candidate_sides_active[round] <-
      rounds$active_dam_candidates[round] +
        rounds$active_sire_candidates[round]
    rounds$configurations_scored[round] <-
      sum(focal$configurations_total[scored])
  }
  likelihood_current <- result$likelihood_history_current
  if (is.null(likelihood_current))
    likelihood_current <- rep(NA, length(result$likelihood_history))
  if (length(likelihood_current) != length(result$likelihood_history))
    stop("native likelihood-current flags must cover the likelihood history", call. = FALSE)
  likelihood <- data.frame(
    round = seq.int(0L, length(result$likelihood_history) - 1L),
    total_likelihood_ln = as.numeric(result$likelihood_history),
    total_likelihood_current = as.logical(likelihood_current),
    stringsAsFactors = FALSE)
  list(
    schema_version = "simp-native-stage1-rounds-v1",
    source_contract = "parents:CheckPedigree-Parentage-UpdateAllProbs-TotLL",
    status = as.character(result$status), converged = isTRUE(result$converged),
    rounds = as.integer(result$rounds), max_rounds = as.integer(max_rounds),
    benchmark_cap = FALSE, initial_pedigree = initial,
    pedigree = .parentage_state_data_frame(state)[, c("id", "dam", "sire")],
    state = state, candidate_pool = pool, snapshot = NULL,
    initial_total_likelihood_ln = as.numeric(result$initial_total_likelihood_ln),
    check_pedigree = check, round_summary = rounds, changes = changes,
    message_revision_history = result$message_revision_history,
    likelihood_history = likelihood, accepted_moves_raw = accepted_raw,
    accepted_recheck_moves_raw = accepted_recheck_raw,
    pass_results = lapply(seq_along(move_blocks), function(round) list(
      accepted_moves = move_blocks[[round]],
      focal_summary = focal_blocks[[round]],
      configuration_diagnostics = if (length(result$configuration_diagnostics_blocks))
        result$configuration_diagnostics_blocks[[round]] else NULL,
      diagnostic_graph_states = if (length(result$graph_state_diagnostic_blocks))
        result$graph_state_diagnostic_blocks[[round]] else NULL,
      relationship_provenance = if (length(result$relationship_provenance_blocks))
        result$relationship_provenance_blocks[[round]] else NULL,
      configuration_locus_diagnostics = if (
        length(result$configuration_locus_diagnostics_blocks))
        result$configuration_locus_diagnostics_blocks[[round]] else NULL,
      shortlist_diagnostics = materialized$diagnostics)),
    recheck_pass_results = lapply(seq_along(recheck_move_blocks), function(wave)
      list(
        accepted_moves = recheck_move_blocks[[wave]],
        focal_summary = if (length(result$recheck_focal_summary_blocks))
          result$recheck_focal_summary_blocks[[wave]] else NULL,
        configuration_diagnostics = if (
          length(result$recheck_configuration_diagnostics_blocks))
          result$recheck_configuration_diagnostics_blocks[[wave]] else NULL,
        diagnostic_graph_states = if (
          length(result$recheck_graph_state_diagnostic_blocks))
          result$recheck_graph_state_diagnostic_blocks[[wave]] else NULL,
        relationship_provenance = if (
          length(result$recheck_relationship_provenance_blocks))
          result$recheck_relationship_provenance_blocks[[wave]] else NULL,
        recheck_phase = if (is.data.frame(result$recheck_phases) &&
                            nrow(result$recheck_phases) >= wave)
          result$recheck_phases[wave, , drop = FALSE] else NULL)),
    recheck_phases = result$recheck_phases,
    native_result = result, static_input = input,
    native_engine_seconds = native_seconds,
    static_input_seconds = input_seconds)
}
