# Admission policy is independent of the numerical estimator. In particular,
# a point estimate alone never authorizes a second search.
.parentage_error_bound_decision <- function(fit, supplied,
                                            direction = c("raise_only", "lower_only", "two_sided")) {
  direction <- match.arg(direction)
  if (length(supplied) > 1L) {
    if (is.null(names(supplied)) || anyDuplicated(names(supplied)))
      stop("multiple supplied error rates must have unique names")
    decisions <- lapply(names(supplied), function(name) {
      interval <- if (identical(fit$status, "fitted")) fit$interval[name, ] else NULL
      decision <- .parentage_error_bound_decision(list(status = fit$status, interval = interval), supplied[[name]], direction)
      # Foreign remains in the joint fit/profile likelihood, but is never a
      # calibration target. Preserve its supplied value under every policy.
      if (name == "Foreign") {
        decision$adjust <- decision$raise <- decision$lower <- FALSE
        decision$used_rate <- supplied[[name]]
        decision$status <- "supplied_foreign_retained"
      }
      decision
    })
    names(decisions) <- names(supplied)
    collect <- function(field) setNames(vapply(decisions, `[[`, 0.0, field), names(supplied))
    adjust <- any(vapply(decisions, `[[`, FALSE, "adjust"))
    return(list(adjust = adjust, raise = any(vapply(decisions, `[[`, FALSE, "raise")),
      lower = any(vapply(decisions, `[[`, FALSE, "lower")), direction = direction,
      used_rate = collect("used_rate"), lower_95 = collect("lower_95"), upper_95 = collect("upper_95"),
      status = if (adjust) "parameters_adjusted_to_profile_bounds" else if (identical(fit$status, "fitted"))
        "supplied_rate_retained" else fit$status, parameters = decisions))
  }
  if (!is.numeric(supplied) || length(supplied) != 1L ||
      !is.finite(supplied) || supplied <= 0 || supplied > 1)
    stop("supplied error must be one finite, strictly positive probability; supply a small positive value", call. = FALSE)
  interval <- if (identical(fit$status, "fitted")) fit$interval else NULL
  get_bound <- function(name) {
    x <- unname(interval[name])
    if (length(x) == 1L && is.finite(x) && x >= 0 && x <= 1) x else NA_real_
  }
  lower <- get_bound("lower"); upper <- get_bound("upper")
  effective_direction <- direction
  if (direction == "two_sided") effective_direction <- if (is.finite(lower) && lower > supplied)
    "raise_only" else "lower_only"
  bound <- if (effective_direction == "raise_only") lower else upper
  if (is.finite(bound)) bound <- max(.parentage_positive_rate_floor, bound)
  adjust <- is.finite(bound) &&
    if (effective_direction == "raise_only") bound > supplied else bound < supplied
  list(adjust = adjust, raise = adjust && effective_direction == "raise_only",
       lower = adjust && effective_direction == "lower_only", direction = direction,
       used_rate = if (adjust) bound else supplied, lower_95 = lower, upper_95 = upper,
       status = if (adjust) if (effective_direction == "raise_only") "raised_to_lower_bound" else
         "lowered_to_upper_bound" else if (identical(fit$status, "fitted"))
         "supplied_rate_retained" else fit$status)
}

# Read the converged question records. Complete pairs are closed in ordinary
# rounds, so their final-round row carries no scores: recover the most recent
# scored row carried by that pair through the recheck lifecycle. Keep its actual
# evaluation revision for audit; do not relabel it as a new evaluation.
.parentage_confidence_from_records <- function(ids, parents, input, blocks,
                                               converged, extra_ln = log(10),
                                               final_revision = NULL, current_records = TRUE,
                                               dependency_state = NULL, allow_partial = FALSE) {
  if (!isTRUE(converged) || !identical(input$likelihood_evidence %||% "ancestral", "ancestral"))
    return(NULL)
  if (is.data.frame(parents)) {
    at <- match(ids, parents$id)
    p <- cbind(match(parents$dam[at], ids), match(parents$sire[at], ids))
    p[is.na(p)] <- 0L
  } else p <- as.matrix(parents)
  blocks <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, blocks)
  if (!length(blocks)) return(NULL)
  d <- do.call(rbind, blocks)
  needed <- c("focal", "selected_dam", "selected_sire", "graph_revision_after",
    "identity_margin_ln", "full_relationship_ln", "dam_support_ln", "sire_support_ln", "pair_support_ln")
  if (!all(needed %in% names(d))) return(NULL)
  if (current_records) {
    if (is.null(final_revision)) final_revision <- max(d$graph_revision_after)
    # Skip is an explicit lifecycle status, not a numerical finiteness test.
    if (!"complete_pair_skipped" %in% names(d)) return(NULL)
    d <- d[!is.na(d$complete_pair_skipped) & d$complete_pair_skipped == 0L, , drop = FALSE]
  } else d <- d[is.finite(d$identity_margin_ln) & is.finite(d$full_relationship_ln), , drop = FALSE]
  d <- d[order(d$graph_revision_after, seq_len(nrow(d))), , drop = FALSE]
  d <- d[!duplicated(d$focal, fromLast = TRUE), , drop = FALSE]
  focals <- which(rowSums(p > 0L) == 2L)
  r <- d[match(focals, d$focal), , drop = FALSE]
  usable <- !is.na(r$focal) & r$selected_dam == p[focals, 1L] & r$selected_sire == p[focals, 2L]
  if (current_records && !is.null(dependency_state)) {
    state <- dependency_state
    valid <- isTRUE(state$enabled) && isTRUE(state$converged) &&
      identical(state$schema_version, "complete-trio-dependency-v1") &&
      isTRUE(all.equal(state$final_graph_revision, final_revision)) &&
      identical(unname(state$final_parents), unname(p)) &&
      length(state$last_evaluated_revision) == length(ids) + 1L &&
      length(state$recheck_pending) == length(ids) + 1L
    if (valid) usable <- usable & state$recheck_pending[focals + 1L] == 0L &
      state$last_evaluated_revision[focals + 1L] == r$graph_revision_after
    else usable[] <- FALSE
    # Incomplete numerical records have no reusable confidence certificate.
    usable <- usable & rowSums(is.finite(as.matrix(r[c("identity_margin_ln", "full_relationship_ln",
      "dam_support_ln", "sire_support_ln", "pair_support_ln")]))) == 5L
  } else if (current_records) usable <- usable & !is.na(r$graph_revision_after) &
    r$graph_revision_after == final_revision & is.finite(r$identity_margin_ln) & is.finite(r$full_relationship_ln)
  usable[is.na(usable)] <- FALSE
  if (!allow_partial && any(!usable)) return(NULL)
  needs_audit <- focals[!usable]
  r <- r[usable, , drop = FALSE]; reused <- focals[usable]
  identity_threshold <- input$identity_margin_ln + extra_ln
  threshold <- input$Tassign_ln + extra_ln
  supports <- cbind(r$dam_support_ln, r$sire_support_ln, r$pair_support_ln)
  ok <- r$identity_margin_ln > identity_threshold & r$full_relationship_ln > threshold &
    rowSums(is.finite(supports) & supports > input$Tassign_ln) == 3L
  confident <- rep(FALSE, length(ids)); confident[reused] <- ok
  list(confident = confident, diagnostics = data.frame(id = ids[reused], confident = ok,
    identity_ln = r$identity_margin_ln, relationship_ln = r$full_relationship_ln,
    dam_support_ln = r$dam_support_ln, sire_support_ln = r$sire_support_ln,
    pair_support_ln = r$pair_support_ln, evaluated_revision = r$graph_revision_after),
    identity_threshold_ln = identity_threshold, relationship_threshold_ln = threshold,
    extra_ln = extra_ln, frozen_state_nonmutation = TRUE, native_audit_calls = 0L,
    native_audit_focals = 0L, needs_audit = as.integer(needs_audit), reused_focals = as.integer(reused),
    source = "converged_question_records")
}

# Certify confidence on the fixed final pedigree, never on an earlier accepted
# visit whose ancestry might have changed. Each audit is explicitly nonmutating.
.parentage_confident_final_trios <- function(compact, input, parents,
                                             extra_ln = log(10), batch_size = 16L,
                                             focals = NULL) {
  n <- compact$n_individuals
  ids <- names(compact$id_index)
  if (is.data.frame(parents)) {
    at <- match(ids, parents$id)
    p <- cbind(match(parents$dam[at], ids), match(parents$sire[at], ids))
    p[is.na(p)] <- 0L
  } else p <- as.matrix(parents)
  stopifnot(identical(dim(p), c(n, 2L)), !anyNA(p))
  input$parents <- p
  input$start_graph_revision <- input$message_revision <- 1L
  input$private_frozen_decision_audit <- TRUE
  # Confidence reads only joint identity, support and relationship contrasts.
  # Preserve configuration rows without requesting unused conditional raw Q.
  input$private_decision_only_diagnostics <- isTRUE(input$likelihood_confidence_raw_elision %||% TRUE)
  input$configuration_diagnostics <- TRUE
  input$configuration_locus_diagnostics <- FALSE
  # A captured search checkpoint may record only one diagnostic position.
  # Every explicitly requested audit focal needs its own configuration rows.
  input$diagnostic_start_position <- 1L
  input$diagnostic_end_position <- 0L

  input$recheck_only <- TRUE
  input$diagnostic_visit_recheck_order <- NULL
  input$progress_callback <- input$phase_checkpoint_callback <- input$checkpoint_callback <- NULL
  input$visit_progress_callback <- NULL
  input$likelihood_replay_initial_parents <- input$likelihood_replay_moves <- NULL
  input$pedigree_probe_requests <- input$pedigree_identity_configurations <- NULL
  input$focal_start_position <- input$focal_end_position <- 1L
  confident <- rep(FALSE, n)
  threshold <- input$Tassign_ln + extra_ln
  identity_threshold <- input$identity_margin_ln + extra_ln
  rows <- list();audit_work <- list()
  if (length(batch_size) != 1L || !is.numeric(batch_size) ||
      !is.finite(batch_size) || batch_size < 1L || batch_size > 512L ||
      batch_size != floor(batch_size)) stop("batch_size must be an integer in 1..512")
  complete <- which(rowSums(p > 0L) == 2L)
  if (is.null(focals)) focals <- complete
  if (!is.numeric(focals) || anyNA(focals) || any(!is.finite(focals)) ||
      any(focals != floor(focals)) || any(!focals %in% complete))
    stop("audit focals must identify complete typed trios", call. = FALSE)
  focals <- unique(as.integer(focals))
  # Reuse one frozen engine across a small bounded block, avoiding O(N)
  # graph preparation per trio without retaining a whole-panel diagnostic grid.
  batches <- split(focals, ceiling(seq_along(focals) / batch_size))
  if (length(focals)) {
    compiled <- build_snapshot_cpp(compact)
    invisible(publish_messages_cpp(compiled, p, 1L, 1L, 1L))
  }
  for (batch in batches) {
    input$diagnostic_visit_order <- input$recheck_focals <- as.integer(batch)
    z <- run_parentage_cpp(compiled, compact, input)
    stopifnot(isTRUE(z$frozen_state_audit),
              identical(unname(z$final_parents), unname(p)), nrow(z$accepted_moves) == 0L)
    if (is.numeric(z$likelihood_work))
      audit_work[[length(audit_work)+1L]] <- z$likelihood_work[c("operations", "questions", "questions_with_ancestry_cut")]
    d <- as.data.frame(z$configuration_diagnostics)
    for (f in batch) {
      r <- d[d$focal == f & d$dam == p[f, 1L] & d$sire == p[f, 2L], , drop = FALSE]
      if (nrow(r) != 1L) next
      ok <- isTRUE(r$hard_valid == 1L) &&
        is.finite(r$identity_margin_ln) && r$identity_margin_ln > identity_threshold &&
        is.finite(r$full_relationship_ln) && r$full_relationship_ln > threshold &&
        all(is.finite(c(r$dam_support_ln, r$sire_support_ln, r$pair_support_ln))) &&
        min(r$dam_support_ln, r$sire_support_ln, r$pair_support_ln) > input$Tassign_ln
      confident[f] <- ok
      rows[[length(rows) + 1L]] <- data.frame(id = ids[f], confident = ok,
        identity_ln = r$identity_margin_ln, relationship_ln = r$full_relationship_ln,
        dam_support_ln = r$dam_support_ln, sire_support_ln = r$sire_support_ln,
        pair_support_ln = r$pair_support_ln)
    }
  }
  list(confident = confident, diagnostics = if (length(rows)) do.call(rbind, rows) else data.frame(),
       identity_threshold_ln = identity_threshold, relationship_threshold_ln = threshold,
       extra_ln = extra_ln, frozen_state_nonmutation = TRUE,
       native_audit_calls = length(batches), native_audit_focals = length(focals),
       numerical_work = if (length(audit_work)) Reduce(`+`, audit_work) else NULL,
       batch_size = as.integer(batch_size))
}

# The native convergence lifecycle owns freshness. Reuse a trio's actual
# evaluation revision while unrelated graph revisions advance; audit only the
# missing, incomplete or pending records. No earlier score is relabelled.
.parentage_final_confidence <- function(compact, input, parents, blocks,
                                        converged, final_revision, dependency_state = NULL,
                                        extra_ln = log(10)) {
  started <- proc.time()[["elapsed"]]
  ids <- names(compact$id_index)
  dependency <- isTRUE(input$likelihood_confidence_dependency_records %||% TRUE)
  current <- isTRUE(input$likelihood_current_confidence_records %||% TRUE)
  records <- .parentage_confidence_from_records(ids, parents, input, blocks, converged,
    extra_ln, final_revision, current_records = current,
    dependency_state = if (dependency && current) dependency_state else NULL,
    allow_partial = dependency && current)
  record_seconds <- proc.time()[["elapsed"]] - started
  if (is.null(records)) {
    audit <- .parentage_confident_final_trios(compact, input, parents, extra_ln)
    audit$record_seconds <- record_seconds
    audit$audit_seconds <- proc.time()[["elapsed"]] - started - record_seconds
    audit$reused_focals <- integer()
    return(audit)
  }
  needs <- records$needs_audit
  records$record_seconds <- record_seconds; records$audit_seconds <- 0
  if (!length(needs)) return(records)
  started <- proc.time()[["elapsed"]]
  audit <- .parentage_confident_final_trios(compact, input, parents, extra_ln, focals = needs)
  records$audit_seconds <- proc.time()[["elapsed"]] - started
  records$confident[needs] <- audit$confident[needs]
  a <- audit$diagnostics; d <- records$diagnostics
  if (nrow(a)) {
    a$evaluated_revision <- NA_integer_
    d <- rbind(d, a[names(d)])
    d <- d[order(match(d$id, ids)), , drop = FALSE]; rownames(d) <- NULL
  }
  records$diagnostics <- d
  records$native_audit_calls <- audit$native_audit_calls
  records$native_audit_focals <- audit$native_audit_focals
  records$numerical_work <- audit$numerical_work
  records$source <- "dependency_current_records_plus_missing_or_stale_audits"
  records
}

