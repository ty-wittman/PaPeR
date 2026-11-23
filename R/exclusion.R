#' Per-locus false-mismatch probability for a true PO pair
#' @param tables per-locus tables from \code{build_tables}.
#' @return data.frame(locus, p_fm, p_um, p_me_true_trio): P(OH | true PO),
#'   P(OH | unrelated), and P(trio ME | true parent-parent-offspring trio),
#'   conditional on the relevant observations being nonmissing.
#' @keywords internal
#' @noRd
po_mismatch_prob <- function(tables) {
  .require_positive_error_rates(tables)
  po_mismatch_prob_cpp(tables)
}

#' Poisson-binomial pmf over 0..length(probs) successes
#' @keywords internal
#' @noRd
.pois_binom_pmf <- function(probs) {
  if (!is.numeric(probs) || is.object(probs) || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs < 0 | probs > 1))
    stop("Poisson-binomial probabilities must be finite values in [0, 1]",
         call. = FALSE)
  d <- 1
  for (p in probs) d <- c(d * (1 - p), 0) + c(0, d * p)
  d
}

.validate_retention_target <- function(retain) {
  if (!is.numeric(retain) || length(retain) != 1L || is.na(retain) ||
      !is.finite(retain) || retain <= 0 || retain > 1)
    stop("retain must be one finite value in (0, 1]", call. = FALSE)
  as.numeric(retain)
}

.calibrate_screen_threshold <- function(probs, retain, name, threshold = NULL) {
  retain <- .validate_retention_target(retain)
  probs <- as.numeric(probs)
  pmf <- .pois_binom_pmf(probs)
  if (is.null(threshold)) {
    cdf <- cumsum(pmf)
    threshold <- which(cdf >= retain)[1L] - 1L
  } else {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        is.na(threshold) || !is.finite(threshold) || threshold != floor(threshold) ||
        threshold < 0)
      stop(name, " threshold must be one nonnegative integer",
           call. = FALSE)
    threshold <- as.integer(threshold)
  }
  keep_n <- if (threshold >= length(pmf) - 1L) length(pmf) else
    threshold + 1L
  retained <- sum(pmf[seq_len(keep_n)])
  list(T = as.integer(threshold), target_retention = retain,
       retained = retained, achieved_retention = retained,
       expected_fm = sum(probs), expected_me = sum(probs),
       expected_false_me_count = sum(probs), pmf = pmf,
       probabilities = probs)
}

#' Calibrate the OH exclusion threshold T (generalizes Sequoia's CalcMaxMismatch)
#'
#' Chooses the smallest integer T such that a true PO pair keeps >= \code{retain}
#' of its probability mass at OH count <= T, using the exact Poisson-binomial
#' distribution of per-locus false mismatches. Because per-pair missingness only
#' lowers the expected mismatch count, the full-panel T is conservative (never
#' over-excludes true PO) for pairs with missing loci.
#'
#' @param p_fm per-locus false-mismatch probabilities (from \code{po_mismatch_prob}).
#' @param retain target retention of true PO pairs (default 0.999).
#' @return list(T, retained, expected_fm, pmf): threshold and diagnostics.
#' @keywords internal
#' @noRd
calibrate_maxmismatch <- function(p_fm, retain = 0.999) {
  if (is.data.frame(p_fm)) p_fm <- p_fm$p_fm
  .calibrate_screen_threshold(p_fm, retain, "OH")
}

#' Calibrate the maximum tolerated trio Mendelian-error count
#'
#' Chooses the smallest integer threshold retaining at least \code{retain} of
#' the exact true-trio probability mass, using the heterogeneous
#' Poisson-binomial distribution of \code{p_me_true_trio}. The probabilities are
#' conditional on all three calls being nonmissing. A full-panel threshold is
#' conservative because missing calls reduce the number of counted errors.
#'
#' @param p_me_true_trio per-locus true-trio Mendelian-error probabilities, or a
#'   data frame returned by \code{po_mismatch_prob}.
#' @param retain target retention (default 0.999).
#' @return list with threshold, target/achieved retention, expected false ME
#'   count, pmf, and the probabilities used for calibration.
#' @keywords internal
#' @noRd
calibrate_max_mendelian_error <- function(p_me_true_trio, retain = 0.999) {
  if (is.data.frame(p_me_true_trio)) {
    if (!"p_me_true_trio" %in% names(p_me_true_trio))
      stop("data frame must contain p_me_true_trio", call. = FALSE)
    p_me_true_trio <- p_me_true_trio$p_me_true_trio
  }
  out <- .calibrate_screen_threshold(
    p_me_true_trio, retain, "MaxMendelE")
  out$expected_me <- out$expected_false_me_count
  out$p_me_true_trio <- out$probabilities
  out
}

