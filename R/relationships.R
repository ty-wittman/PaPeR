# Typed-parent relationship screening and diagnostics.

#' Screen typed individuals for a relationship configuration
#'
#' Scores only the canonical typed-parent relationship universe. This is a
#' diagnostic and candidate-screening layer for typed individuals.
#' @param GenoM public one-column-per-locus character genotypes.
#' @param tables per-locus likelihood tables from \code{build_tables}.
#' @param focal focal relationship name, or \code{"RELATED"} for all retained
#'   non-unrelated typed relationships.
#' @param alts optional alternative relationship names. Unrelated alternatives are always retained. Duplicate detection is separate.
#' @param thresh minimum public base-10 contrast to report.
#' @param LifeHistData optional life-history data for an age-difference screen.
#' @param max_agediff maximum absolute birth-year gap; \code{-1} disables it.
#' @param min_agediff minimum absolute birth-year gap.
#' @param oh_max maximum tolerated opposing-homozygote loci; \code{-1} disables.
#' @param Complex \code{"simp"} for the ordinary universe or \code{"full"}
#'   for static ordinary screening. Contextual full rows require the internal
#'   Section-2 conditional scorer and are not unconditional screening rows.
#' @param operation internal resolver operation used by duplicate and diagnostic
#'   wrappers; ordinary callers should use the default \code{"pair"}.
#' @return A data frame of typed pair diagnostics.
#' The public \code{llr} and \code{duplicate_guard_llr} columns are base 10;
#' \code{thresh} is converted once from base 10 to the internal natural-log
#' comparison threshold. The result carries \code{log_base = 10}.
#' @keywords internal
#' @noRd
rel_edges <- function(GenoM, tables, focal, alts = NULL, thresh = 0,
                      LifeHistData = NULL, max_agediff = -1L, min_agediff = 0L,
                      oh_max = NULL, Complex = c("simp", "full"),
                      operation = c("pair", "duplicate", "diagnostic")) {
  Complex <- match.arg(Complex)
  operation <- match.arg(operation)
  GenoM <- .prepare_for_tables(GenoM, tables)
  universe <- .relationship_universe(Complex, focal, alternatives = alts,
                                     operation = operation)
  ctx <- .new_inference_context(GenoM, tables, Complex)
  by <- if (is.null(LifeHistData)) rep(NA_integer_, nrow(GenoM)) else
    as.integer(LifeHistData$BirthYear[match(rownames(GenoM), LifeHistData$ID)])
  if (is.null(oh_max))
    oh_max <- as.integer(ceiling(0.8 * sum(po_mismatch_prob(tables)$p_um)))
  e <- screen_rel_pairs_cpp(
    tables, ctx$G, universe$score_delta,
    universe$focal_indices - 1L, universe$alternative_indices - 1L,
    .threshold_log10_to_ln(thresh), by, as.integer(max_agediff), as.integer(min_agediff),
    as.integer(oh_max))
  focal_labels <- .relationship_winner_labels(e$focal_win, universe, "focal")
  alternative_labels <- .relationship_winner_labels(
    e$alternative_win, universe, "alternative")
  out <- data.frame(
    id1 = rownames(ctx$G)[e$i], id2 = rownames(ctx$G)[e$j],
    llr = .llr_log10(e$llr),
    focal_winner = focal_labels$winner,
    focal_winner_set = focal_labels$set,
    alternative_winner = alternative_labels$winner,
    alternative_winner_set = alternative_labels$set,
    duplicate_guard_llr = rep(NA_real_, length(e$llr)),
    duplicate_guard_status = rep("not_evaluated", length(e$llr)),
    stringsAsFactors = FALSE)

  attr(out, "relationship_universe") <- universe
  attr(out, "log_base") <- .public_log_base
  attr(out, "threshold_log10") <- thresh
  attr(out, "threshold_ln") <- .threshold_log10_to_ln(thresh)
  out
}
