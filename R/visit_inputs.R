# One compiled, fixed-50 Complex = "simp" parentage pass.
#
# This file owns only the search-pass boundary.  Message publication,
# convergence rounds, final-output LLRs, and Complex = "full" remain outside it.

.candidate_side_inputs <- function(pool, state, child, side, ids) {
  key <- .candidate_pool_key(child, side)
  side_rows <- if (!is.null(pool$row_indices_by_key) &&
                   !is.null(pool$row_indices_by_key[[key]])) {
    pool$row_indices_by_key[[key]]
  } else {
    which(pool$rows$offspring == child & pool$rows$side == side)
  }
  rows <- pool$shortlist_rows[[key]]
  if (is.null(rows)) rows <- integer()
  active <- if (length(rows)) pool$rows$candidate[rows] else character()
  current <- .current_parent_configuration(state, child)[[side]]
  current <- current[!is.na(current) & nzchar(current)]
  # The immutable sparse graph contains discovery survivors only.  An
  # incumbent is added locally for a later-round comparison and is never
  # written back into that graph.
  active <- unique(c(active, current))
  indices <- match(active, ids)
  if (anyNA(indices)) stop("whole-pass candidate ID is absent from state",
                           call. = FALSE)
  matched <- match(active, pool$rows$candidate[side_rows])
  absent <- is.na(matched)
  # A current supplied/inferred parent must remain scorable for temporary
  # removal or replacement even when it failed a permanent screen and was not
  # retained in the compact survivor table.  It receives no screening bonus;
  # the compiled scorer evaluates its actual local configuration as before.
  if (any(absent & !active %in% current))
    stop("whole-pass candidate metadata is absent from candidate pool",
         call. = FALSE)
  row_index <- rep.int(NA_integer_, length(active))
  row_index[!absent] <- side_rows[matched[!absent]]
  quick <- rep.int(-Inf, length(active))
  if (any(!absent)) {
    quick[!absent] <- pool$rows$po_u_llr_ln[row_index[!absent]]
    quick[is.na(quick)] <- -Inf
  }
  list(ids = as.integer(indices), quick = as.numeric(quick),
       forced = current,
       rows = as.integer(row_index),
       permanent = as.logical(active %in% pool$rows$candidate[rows]))
}

.visit_inputs <- function(snapshot, state, pool, ids, universe,
                                    Tassign_ln, margin_ln, tol,
                                    shortlist_mode) {
  n <- length(ids)
  focal_order <- unname(state$id_to_index[state$focal_order])
  focal_order <- as.integer(focal_order)
  dam_offsets <- integer(n + 1L)
  sire_offsets <- integer(n + 1L)
  dam_ids <- integer(); sire_ids <- integer()
  dam_quick <- numeric(); sire_quick <- numeric()
  dam_permanent <- logical(); sire_permanent <- logical()
  forced <- list(dam = logical(n), sire = logical(n))
  for (o in seq_len(n)) {
    child <- ids[o]
    dam <- .candidate_side_inputs(pool, state, child, "dam", ids)
    sire <- .candidate_side_inputs(pool, state, child, "sire", ids)
    retained_dam <- pool$rows$candidate[
      pool$shortlist_rows[[.candidate_pool_key(child, "dam")]]]
    retained_sire <- pool$rows$candidate[
      pool$shortlist_rows[[.candidate_pool_key(child, "sire")]]]
    forced$dam[o] <- length(dam$forced) > 0L &&
      !dam$forced %in% retained_dam
    forced$sire[o] <- length(sire$forced) > 0L &&
      !sire$forced %in% retained_sire
    dam_offsets[o] <- length(dam_ids)
    sire_offsets[o] <- length(sire_ids)
    dam_ids <- c(dam_ids, dam$ids)
    sire_ids <- c(sire_ids, sire$ids)
    dam_quick <- c(dam_quick, dam$quick)
    sire_quick <- c(sire_quick, sire$quick)
    dam_permanent <- c(dam_permanent, dam$permanent)
    sire_permanent <- c(sire_permanent, sire$permanent)

  }
  dam_offsets[n + 1L] <- length(dam_ids)
  sire_offsets[n + 1L] <- length(sire_ids)
  rank <- integer(n)
  rank[order(ids)] <- seq_len(n)
  lh <- state$life_history
  list(
    parents = state$parents,
    focal_order = focal_order,
    dam_offsets = as.integer(dam_offsets),
    sire_offsets = as.integer(sire_offsets),
    dam_ids = as.integer(dam_ids), sire_ids = as.integer(sire_ids),
    dam_quick = as.numeric(dam_quick), sire_quick = as.numeric(sire_quick),
    dam_permanent = as.logical(dam_permanent),
    sire_permanent = as.logical(sire_permanent),
    id_order_rank = as.integer(rank),
    sex = as.integer(state$sex[ids]),
    birth_lo = as.integer(lh$BY.min[match(ids, lh$ID)]),
    birth_hi = as.integer(lh$BY.max[match(ids, lh$ID)]),
    year_last = as.integer(state$year_last[ids]),
    min_age = as.numeric(unname(state$min_age_parent[c("M", "P")])),
    max_age = as.numeric(unname(state$max_age_parent[c("M", "P")])),
    Tassign_ln = as.numeric(Tassign_ln),
    identity_margin_ln = as.numeric(margin_ln), tol = as.numeric(tol),
    tie_tolerance = 1e-6,
    max_mendelian_errors = as.integer(pool$MaxMendelE),
    start_graph_revision = as.integer(state$graph_revision),
    message_revision = as.integer(state$message_revision),
    shortlist_mode = as.character(shortlist_mode),
    search_loop_engine = "bounded_messages",
    benchmark_opt5a_mode = TRUE,
    adaptive_selective_publisher = FALSE,
    adaptive_selective_threshold = 0,
    forced = forced)
}

.message_whole_pass_diagnostics <- function(result, ids, input, pool, state,
                                         shortlist_settings) {
  summary <- result$focal_summary
  if (!is.data.frame(summary) || !nrow(summary))
    stop("compiled whole-pass result did not return focal summaries",
         call. = FALSE)
  focal_ids <- ids[summary$focal]
  interleave <- function(dam, sire)
    as.integer(t(cbind(dam = dam, sire = sire)))
  permanent_count <- function(focal, side)
    sum(pool$rows$offspring == ids[focal] & pool$rows$side == side &
        pool$rows$permanent_pass)
  permanent_dam <- vapply(summary$focal, permanent_count, integer(1),
                           side = "dam")
  permanent_sire <- vapply(summary$focal, permanent_count, integer(1),
                           side = "sire")
  active_size <- data.frame(
    offspring = rep(focal_ids, each = 2L),
    side = rep(c("dam", "sire"), nrow(summary)),
    permanent_pool_size = interleave(permanent_dam, permanent_sire),
    initial_active_size = interleave(
      summary$dam_candidate_count, summary$sire_candidate_count),
    assigned_parent_forced_inclusion = as.logical(interleave(
      input$forced$dam[summary$focal], input$forced$sire[summary$focal])),
    expansion_count = 0L,
    final_active_size = interleave(
      summary$dam_candidate_count, summary$sire_candidate_count),
    best_excluded_quick_score_ln = NA_real_,
    quick_score_boundary_gap_ln = NA_real_,
    co_parent_conditioned_rank_change = NA_integer_,
    winning_identity_margin_ln = rep(summary$identity_margin_ln, each = 2L),
    winning_relationship_margin_ln = NA_real_,
    safety_cap = 50L, cap_hit = FALSE, resolved = TRUE,
    stop_reason = "sequoia50",
    stop_detail = "fixed-50 whole-pass compiled search",
    strict_compatibility_truncation = as.logical(interleave(
      permanent_dam > 50L,
      permanent_sire > 50L)),
    single_scores_evaluated = NA_real_, pair_cells_evaluated = NA_real_,
    pair_cells_reused = NA_real_, pair_cells_evaluated_this_call = NA_real_,
    pair_cells_reused_this_call = NA_real_, pair_cells_possible = NA_real_,
    pair_grid_cells_possible = NA_real_, pair_cells_avoided = NA_real_,
    pair_configurations_possible = NA_real_,
    pair_configurations_evaluated = NA_real_,
    pair_configurations_avoided = NA_real_,
    pair_configurations_evaluated_this_call = NA_real_,
    stringsAsFactors = FALSE)
  active_size$active_sizes_evaluated <- I(lapply(seq_len(nrow(summary)),
    function(i) list(dam = summary$dam_candidate_count[i],
                     sire = summary$sire_candidate_count[i])))
  active_size$expansion_reasons <- I(rep(list(character()), nrow(active_size)))
  active_size
}

