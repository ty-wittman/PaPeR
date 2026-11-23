# Structured Stage-1 candidate pools.

.candidate_pool_key <- function(offspring, side)
  paste(as.character(offspring), as.character(side), sep = "\r")

.candidate_pool_empty_rows <- function() {
  data.frame(
    offspring = character(), side = character(), candidate = character(),
    initial_order = integer(), self_pass = logical(), sex_pass = logical(),
    age_pass = logical(), year_last_pass = logical(), age_year_pass = logical(),
    oh_count = integer(), oh_pass = logical(), po_u_llr_ln = numeric(),
    po_u_llr_log10 = numeric(), tfilter_pass = logical(),
    permanent_pass = logical(), permanent_failure_reason = character(),
    relationship_score_ln = numeric(),
    shortlist_rank = integer(),
    permanent_rank = integer(), stringsAsFactors = FALSE)
}

.candidate_diagnostic_frame <- function(rows) {
  if (!length(rows)) return(do.call(rbind, rows))
  columns <- lapply(names(rows[[1L]]), function(column) {
    unlist(lapply(rows, `[[`, column), use.names = FALSE)
  })
  names(columns) <- names(rows[[1L]])
  out <- do.call(data.frame, c(columns, list(stringsAsFactors = FALSE)))
  attributes(out) <- list(names = names(out), row.names = .row_names_info(out, 0L),
                          class = class(out))
  out
}

# Production sparse Sequoia-style candidate frontend.  Genetic work is done
# exactly once for each hard-feasible unordered pair.  The resulting CSR is
# immutable; rows are only the oriented, retained parent arcs, never an
# offspring-side universe or rejected-pair trace.
.build_candidate_pool <- function(ctx, lh, ids, minp, maxp, T,
                                  Tfilter_log10, universe,
                                  oh_calibration, me_calibration) {
  candidate_group_index <- getOption("PaPeR.candidate_group_index", TRUE)
  if (!is.logical(candidate_group_index) || length(candidate_group_index) != 1L ||
      is.na(candidate_group_index))
    stop("PaPeR.candidate_group_index must be one non-missing logical value",
         call. = FALSE)
  candidate_diagnostic_batch <- getOption(
    "PaPeR.candidate_diagnostic_batch", TRUE)
  if (!is.logical(candidate_diagnostic_batch) ||
      length(candidate_diagnostic_batch) != 1L ||
      is.na(candidate_diagnostic_batch))
    stop(paste0("PaPeR.candidate_diagnostic_batch must be one ",
                "non-missing logical value"), call. = FALSE)
  candidate_group_ordinal <- getOption(
    "PaPeR.candidate_group_ordinal", TRUE)
  if (!is.logical(candidate_group_ordinal) ||
      length(candidate_group_ordinal) != 1L ||
      is.na(candidate_group_ordinal))
    stop(paste0("PaPeR.candidate_group_ordinal must be one ",
                "non-missing logical value"), call. = FALSE)
  if (!inherits(ctx, "PaPeR_inference_context"))
    stop("ctx must be an inference context", call. = FALSE)
  if (!inherits(universe, "PaPeR_relationship_universe"))
    stop("universe must be a relationship universe", call. = FALSE)
  if (!is.numeric(Tfilter_log10) || length(Tfilter_log10) != 1L ||
      is.na(Tfilter_log10) || is.nan(Tfilter_log10) ||
      (is.infinite(Tfilter_log10) && Tfilter_log10 > 0))
    stop("Tfilter must be one numeric base-10 threshold", call. = FALSE)
  ids <- as.character(ids)
  if (!identical(ids, as.character(lh$ID)))
    stop("ids must match normalized life-history order", call. = FALSE)
  n <- length(ids)
  started <- proc.time()[["elapsed"]]
  tfilter_ln <- Tfilter_log10 * log(10)

  # Direction masks are biological only: sex is intentionally deferred until
  # after the symmetric genotype screen.  The native frontend fuses exact
  # contiguous-support/Year.last feasibility with OH early-exit screening and
  # returns survivors in the original i<j order.
  birth_lo <- as.integer(lh$BY.min)
  birth_hi <- as.integer(lh$BY.max)
  exact <- !is.na(lh$BirthYear)
  birth_lo[exact] <- as.integer(lh$BirthYear[exact])
  birth_hi[exact] <- as.integer(lh$BirthYear[exact])
  year_last <- as.integer(lh$Year.last)
  min_age <- .normalize_age_vector(minp, "MinAgeParent")
  max_age <- .normalize_age_vector(maxp, "MaxAgeParent")
  fused_started <- proc.time()[["elapsed"]]
  passed <- screen_candidate_pairs_sparse_cpp(
    ctx$G, birth_lo, birth_hi, year_last, min_age, max_age, as.integer(T))
  fused_seconds <- proc.time()[["elapsed"]] - fused_started
  hard_seconds <- as.numeric(attr(
    passed, "candidate_hard_feasibility_seconds"))
  oh_seconds <- as.numeric(attr(passed, "candidate_oh_seconds"))

  # Unconditional PO/U is symmetric.  This is one compiled batch for each OH
  # survivor, deliberately excluding MZ, ordinary registries, pedigrees, and
  # frozen-message likelihoods.
  qpo_started <- proc.time()[["elapsed"]]
  if (nrow(passed)) {
    delta_rows <- match(c("PO", "U"), rownames(ctx$Delta))
    if (anyNA(delta_rows)) stop("compact registry lacks PO/U rows", call. = FALSE)
    qpo <- score_basis_pairs_cpp(ctx$markers, ctx$obsidx, passed$left,
                                 passed$right,
                                 ctx$Delta[delta_rows, , drop = FALSE])
    passed$po_u_llr_ln <- .log_likelihood_contrast_ln(
      qpo[, 1L], qpo[, 2L], "sparse PO-versus-unrelated screening contrast")
  } else {
    passed$po_u_llr_ln <- numeric()
  }
  qpo_seconds <- proc.time()[["elapsed"]] - qpo_started

  orient_started <- proc.time()[["elapsed"]]
  sex <- as.integer(lh$Sex)
  add_arc <- function(bit, side, parent_from_left) {
    keep <- bitwAnd(passed$direction_mask, bit) != 0L
    if (!any(keep)) return(NULL)
    parent <- if (parent_from_left) passed$left[keep] + 1L else
      passed$right[keep] + 1L
    child <- if (parent_from_left) passed$right[keep] + 1L else
      passed$left[keep] + 1L
    role_ok <- if (side == "dam") is.na(sex[parent]) | sex[parent] == 1L else
      is.na(sex[parent]) | sex[parent] == 2L
    if (!any(role_ok)) return(NULL)
    hit <- which(keep)[role_ok]
    data.frame(
      offspring = ids[child[role_ok]], side = side,
      candidate = ids[parent[role_ok]], initial_order = NA_integer_,
      self_pass = TRUE, sex_pass = TRUE, age_pass = TRUE,
      year_last_pass = TRUE, age_year_pass = TRUE,
      oh_count = passed$oh_count[hit], oh_pass = TRUE,
      po_u_llr_ln = passed$po_u_llr_ln[hit],
      po_u_llr_log10 = passed$po_u_llr_ln[hit] / log(10),
      tfilter_pass = FALSE, permanent_pass = FALSE,
      permanent_failure_reason = "none",
      relationship_score_ln = passed$po_u_llr_ln[hit],
      shortlist_rank = NA_integer_, permanent_rank = NA_integer_,
      stringsAsFactors = FALSE)
  }
  rows <- do.call(rbind, Filter(Negate(is.null), list(
    add_arc(1L, "dam", TRUE), add_arc(2L, "sire", TRUE),
    add_arc(4L, "dam", FALSE), add_arc(8L, "sire", FALSE))))
  if (is.null(rows)) rows <- .candidate_pool_empty_rows()

  all_keys <- unlist(lapply(ids, function(id)
    c(.candidate_pool_key(id, "dam"), .candidate_pool_key(id, "sire"))),
    use.names = FALSE)
  diagnostics <- vector("list", length(all_keys))
  plain_numeric_scalar <- function(value) {
    if (!(typeof(value) %in% c("integer", "double")) ||
        length(value) != 1L || !is.null(attributes(value)))
      return(FALSE)
    TRUE
  }
  diagnostic_operands <- list(
    oh_calibration$target_retention,
    oh_calibration$achieved_retention,
    me_calibration$T,
    me_calibration$target_retention,
    me_calibration$achieved_retention,
    me_calibration$expected_false_me_count,
    Tfilter_log10, hard_seconds, oh_seconds)
  batch_diagnostics <- isTRUE(candidate_diagnostic_batch) &&
    length(all_keys) > 0L && !anyNA(ids) && all(nzchar(ids)) &&
    !anyDuplicated(ids) && all(!grepl("\r", ids, fixed = TRUE)) &&
    all(vapply(diagnostic_operands, plain_numeric_scalar, logical(1)))
  diagnostic_row <- if (batch_diagnostics) {
    function(..., stringsAsFactors = FALSE) list(...)
  } else data.frame
  keys <- paste(rows$offspring, rows$side, sep = "\r")
  fast_group_index <- isTRUE(candidate_group_index) &&
    !anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids) &&
    all(!grepl("\r", ids, fixed = TRUE)) &&
    !anyNA(all_keys) && all(nzchar(all_keys)) && !anyDuplicated(all_keys) &&
    !anyNA(rows$offspring) && !anyNA(rows$side) &&
    all(nzchar(rows$offspring)) && all(nzchar(rows$side)) &&
    all(!grepl("\r", rows$offspring, fixed = TRUE)) &&
    all(rows$side %in% c("dam", "sire")) && all(keys %in% all_keys)
  prefilter_indices <- if (fast_group_index) {
    split(seq_along(keys), factor(keys, levels = all_keys), drop = FALSE)
  } else NULL
  tfilter_examined <- 0L; tfilter_rejected <- 0L
  for (at in seq_along(all_keys)) {
    key <- all_keys[[at]]
    idx <- if (fast_group_index) prefilter_indices[[at]] else which(keys == key)
    # Screen every side before ranking and applying the fixed-50 cap.
    apply_tfilter <- !(is.infinite(tfilter_ln) && tfilter_ln < 0)
    keep <- if (apply_tfilter) rows$po_u_llr_ln[idx] >= tfilter_ln else
      ifelse(is.na(rows$po_u_llr_ln[idx]), NA, TRUE)
    known_rejection <- !is.na(keep) & !keep
    if (apply_tfilter) {
      tfilter_examined <- tfilter_examined + length(idx)
      tfilter_rejected <- tfilter_rejected + sum(known_rejection)
    }
    retained <- idx[is.na(keep) | keep]
    if (length(retained)) {
      rows$tfilter_pass[retained] <- keep[is.na(keep) | keep]
      rows$permanent_pass[retained] <- TRUE
    }
    if (length(idx) && any(known_rejection))
      rows$permanent_failure_reason[idx[known_rejection]] <- "tfilter"
    split_key <- strsplit(key, "\r", fixed = TRUE)[[1L]]
    diagnostics[[at]] <- diagnostic_row(
      offspring = split_key[[1L]], side = split_key[[2L]],
      n_initial = NA_integer_, n_after_self = NA_integer_,
      n_after_sex = length(idx), n_after_age = length(idx),
      n_after_year_last = length(idx), n_after_age_year_support = length(idx),
      n_after_oh = length(idx), n_after_tfilter = length(retained),
      n_permanent = length(retained), OH_threshold = as.integer(T),
      OH_target_retention = oh_calibration$target_retention,
      OH_achieved_retention = oh_calibration$achieved_retention,
      MaxMendelE = me_calibration$T,
      ME_target_retention = me_calibration$target_retention,
      ME_achieved_retention = me_calibration$achieved_retention,
      expected_false_ME_count = me_calibration$expected_false_me_count,
      Tfilter_log10 = Tfilter_log10,
      PO_U_min_before_Tfilter = if (length(idx)) min(rows$po_u_llr_ln[idx]) else NA_real_,
      PO_U_median_before_Tfilter = if (length(idx)) stats::median(rows$po_u_llr_ln[idx]) else NA_real_,
      PO_U_max_before_Tfilter = if (length(idx)) max(rows$po_u_llr_ln[idx]) else NA_real_,
      n_current_excluded_role_consistency = 0L, n_current_excluded_cycles = 0L,
      n_evaluated_co_parent_me_tuples = 0L, n_failed_co_parent_me_tuples = 0L,
      n_active_excluded_role_consistency = 0L, n_active_excluded_cycles = 0L,
      runtime_co_parent_me_seconds = 0, runtime_self_seconds = 0,
      runtime_sex_seconds = 0, runtime_age_seconds = hard_seconds,
      runtime_year_last_seconds = 0, runtime_age_year_support_seconds = 0,
      runtime_oh_seconds = oh_seconds, runtime_tfilter_seconds = 0,
      stringsAsFactors = FALSE)
  }
  rows <- rows[rows$permanent_pass, , drop = FALSE]
  rownames(rows) <- NULL
  postfilter_indices <- if (fast_group_index) {
    postfilter_keys <- paste(rows$offspring, rows$side, sep = "\r")
    split(seq_along(postfilter_keys),
          factor(postfilter_keys, levels = all_keys), drop = FALSE)
  } else NULL
  row_indices_by_key <- setNames(vector("list", length(all_keys)), all_keys)
  shortlist_rows <- setNames(vector("list", length(all_keys)), all_keys)
  for (at in seq_along(all_keys)) {
    key <- all_keys[[at]]
    slot <- if (fast_group_index && isTRUE(candidate_group_ordinal)) at else key
    idx <- if (fast_group_index) postfilter_indices[[slot]] else
      which(rows$offspring == sub("\r.*$", "", key) &
            rows$side == sub("^.*\r", "", key))
    row_indices_by_key[[slot]] <- idx
    ord <- if (length(idx)) idx[order(-rows$po_u_llr_ln[idx],
                                      rows$candidate[idx])] else integer()
    if (length(ord)) rows$permanent_rank[ord] <- seq_along(ord)
    chosen <- head(ord, 50L)
    if (length(chosen)) rows$shortlist_rank[chosen] <- seq_along(chosen)
    shortlist_rows[[slot]] <- chosen
  }
  orient_seconds <- proc.time()[["elapsed"]] - orient_started
  # The symmetric screen also applies to the native candidate stream. Keep
  # unknown contrasts conservatively, and retain the unfiltered `passed` for
  # OH/QPO work accounting. CSR remains independent of sex and pedigree state.
  edge_keep <- is.na(passed$po_u_llr_ln) |
    passed$po_u_llr_ln >= tfilter_ln
  if (is.infinite(tfilter_ln) && tfilter_ln < 0)
    edge_keep[] <- TRUE
  edges <- passed[edge_keep, , drop = FALSE]
  owner <- c(edges$left + 1L, edges$right + 1L)
  neighbor <- c(edges$right + 1L, edges$left + 1L)
  score <- rep(edges$po_u_llr_ln, 2L)
  mask <- c(edges$direction_mask, bitwShiftR(edges$direction_mask, 2L))
  ord <- order(owner, neighbor)
  degree <- tabulate(owner, n)
  csr <- list(offsets = as.integer(c(0L, cumsum(degree))),
              candidate_ids = as.integer(neighbor[ord]),
              po_u_llr_ln = as.numeric(score[ord]),
              direction_mask = as.integer(mask[ord]))
  total_pairs <- as.double(n) * as.double(n - 1L) / 2
  hard_pairs <- as.numeric(attr(passed, "candidate_hard_feasible_pairs"))
  oh_survivors <- as.numeric(attr(passed, "candidate_oh_survivors"))
  rejected_pairs <- as.numeric(attr(passed, "oh_rejected_pairs"))
  rejected_loci_total <- as.numeric(attr(
    passed, "oh_rejected_loci_examined_total"))
  total_oh_loci <- as.numeric(attr(passed, "oh_total_loci_examined"))
  mean_rejected_loci <- if (rejected_pairs > 0) {
    rejected_loci_total / rejected_pairs
  } else NA_real_
  me_stats <- new.env(parent = emptyenv())
  me_stats$evaluated <- 0L; me_stats$failed <- 0L; me_stats$runtime_seconds <- 0
  me_stats$by_side <- list()
  telemetry <- list(
    candidate_unordered_pairs_total = total_pairs,
    candidate_hard_feasible_pairs = hard_pairs,
    candidate_hard_rejected_pairs = total_pairs - hard_pairs,
    oh_pairs_examined = hard_pairs,
    oh_pairs_early_rejected = rejected_pairs,
    oh_pairs_passed = oh_survivors, oh_loci_examined = total_oh_loci,
    mean_oh_loci_before_rejection = mean_rejected_loci,
    qpo_pairs_scored = nrow(passed), sparse_candidate_edges = nrow(edges),
    maximum_sparse_degree = max(c(0L, degree)), mean_sparse_degree = mean(degree),
    tfilter_examined = tfilter_examined, tfilter_rejected = tfilter_rejected,
    dam_shortlist_edges = sum(rows$side == "dam" & !is.na(rows$shortlist_rank)),
    sire_shortlist_edges = sum(rows$side == "sire" & !is.na(rows$shortlist_rank)),
    maximum_dam_shortlist = max(c(0L, vapply(shortlist_rows[grepl("\rdam$", names(shortlist_rows))], length, integer(1)))),
    maximum_sire_shortlist = max(c(0L, vapply(shortlist_rows[grepl("\rsire$", names(shortlist_rows))], length, integer(1)))),
    incumbent_forced_inclusions_later_rounds = 0L,
    candidate_fused_frontend_wall_seconds = fused_seconds,
    candidate_hard_feasibility_seconds = hard_seconds,
    candidate_oh_screening_wall_seconds = oh_seconds,
    candidate_po_u_wall_seconds = qpo_seconds,
    candidate_discovery_wall_seconds = hard_seconds + oh_seconds + qpo_seconds,
    candidate_orientation_shortlist_wall_seconds = orient_seconds,
    oh_rejected_loci_examined = rejected_loci_total,
    oh_rejected_loci_examined_total = rejected_loci_total)
  pool <- list(
    schema_version = "candidate-sparse-csr-v2", rows = rows,
    diagnostics = if (batch_diagnostics)
      .candidate_diagnostic_frame(diagnostics) else do.call(rbind, diagnostics),
    shortlist_rows = shortlist_rows,
    row_indices_by_key = row_indices_by_key, active_rows = list(), csr = csr,
    pair_edges = edges, telemetry = telemetry,
    shortlist_policy = "sequoia50-sparse-csr", calibration = list(
      OH = oh_calibration, ME = me_calibration), relationship_universe = universe,
    OH_threshold = as.integer(T), Tfilter_log10 = Tfilter_log10,
    Tfilter_ln = tfilter_ln, MaxMendelE = as.integer(me_calibration$T),
    # This field is retained only for downstream diagnostics: production never
    # constructs logical candidate rows for the 2*N^2 universe.
    initial_candidate_rows = 0L, permanent_candidate_rows = as.integer(nrow(rows)),
    me_cache = new.env(parent = emptyenv()), me_stats = me_stats,
    overlay_diagnostics = list(), state_mask_cache = new.env(parent = emptyenv()))
  class(pool) <- c("PaPeR_candidate_pool", "list")
  pool
}
