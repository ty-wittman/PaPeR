# Independent called-locus Mendelian-error observations. No truth pedigree is
# accepted here; parents are the complete, confidently assigned pass-1 trios.
.parentage_error_observations <- function(GenoM, tables, parents,
                                           confident = NULL) {
  g <- .prepare_for_tables(GenoM, tables)
  if (is.data.frame(parents)) {
    if (!all(c("id", "dam", "sire") %in% names(parents)))
      stop("parents must have id, dam and sire columns", call. = FALSE)
    at <- match(rownames(g), parents$id)
    if (anyNA(at)) stop("parents must include every genotyped individual", call. = FALSE)
    p <- cbind(match(parents$dam[at], rownames(g)),
               match(parents$sire[at], rownames(g)))
    p[is.na(p)] <- 0L
  } else {
    p <- as.matrix(parents)
    if (!identical(dim(p), c(nrow(g), 2L)) || anyNA(p) ||
        any(p != floor(p) | p < 0L | p > nrow(g)))
      stop("parents must be an n-individual by two parent-index matrix", call. = FALSE)
  }
  if (is.null(confident)) confident <- rep(TRUE, nrow(g))
  if (!is.logical(confident) || length(confident) != nrow(g) || anyNA(confident))
    stop("confident must be a nonmissing logical value per individual", call. = FALSE)
  child <- which(rowSums(p > 0L) == 2L & confident)
  L <- length(tables)
  called <- incompatible <- matrix(FALSE, length(child), L)
  dam_oh <- sire_oh <- integer(length(child))
  if (length(child)) for (l in seq_len(L)) {
    a <- g[, 2L * l - 1L]; b <- g[, 2L * l]
    d <- p[child, 1L]; s <- p[child, 2L]
    called[, l] <- a[child] > 0L & a[d] > 0L & a[s] > 0L
    in_parent <- function(x, q) x == a[q] | x == b[q]
    compatible <- (in_parent(a[child], d) & in_parent(b[child], s)) |
      (in_parent(b[child], d) & in_parent(a[child], s))
    incompatible[, l] <- called[, l] & !compatible
    opposite <- function(q) a[child] > 0L & a[q] > 0L &
      a[child] == b[child] & a[q] == b[q] & a[child] != a[q]
    dam_oh <- dam_oh + opposite(d); sire_oh <- sire_oh + opposite(s)
  }
  list(child = child, ids = rownames(g)[child], parents = p[child, , drop = FALSE],
       called = called, incompatible = incompatible,
       count = rowSums(incompatible), n_called = rowSums(called),
       dam_oh = dam_oh, sire_oh = sire_oh)
}

# The estimator corrects ME selection only. Audit the additional OH hard
# ceiling without claiming that it is part of the fitted likelihood.
.parentage_error_ceiling_audit <- function(fit, T) {
  obs <- fit$observations
  available <- length(T) == 1L && is.finite(T) &&
    !is.null(obs$dam_oh) && !is.null(obs$sire_oh)
  fit$opposite_homozygote_ceiling <- if (!available)
    list(available = FALSE, modeled = FALSE, reason = 'missing_ceiling_or_observations') else
    list(available = TRUE, modeled = FALSE, T = T, trios = length(obs$dam_oh),
      dam_at_T = sum(obs$dam_oh == T), sire_at_T = sum(obs$sire_oh == T),
      any_at_T = sum(obs$dam_oh == T | obs$sire_oh == T),
      any_above_T = sum(obs$dam_oh > T | obs$sire_oh > T))
  fit
}

# Rows are independent marginal count distributions for the composite fit.
# Keeping probabilities only through the admission ceiling is intentional:
# rowSums is the selection probability used by the truncated likelihood.
.parentage_error_count_pmf <- function(probabilities, called, ceiling) {
  storage.mode(called) <- "logical"
  error_count_pmf_cpp(as.double(probabilities), called, as.integer(ceiling))
}

.parentage_error_scale_spec <- function(tables) {
  .require_positive_error_rates(tables)
  model <- tables[[1L]]$ErrorModel %||% "legacy_sequoia"
  if (model == "legacy_sequoia") {
    initial <- tables[[1L]]$Err
    rebuild <- function(rate) {
      ans <- lapply(tables, function(x) {
        profile <- if (isTRUE(x$ErrProfileDefault)) NULL else x$ErrProfile
        build_locus_table(x$p, Err = rate, ErrProfile = profile)
      })
      names(ans) <- names(tables); ans
    }
    return(list(model = unname(model[1L]), parameter = "Err", initial = unname(initial[1L]),
                upper = .25, rebuild = rebuild,
                rates = function(theta) c(Err = theta)))
  }
  # All public MH callers use the two-parameter category/profile estimator.
  # The former common-scale MH branch had no callers after that dispatch.
  stop("scalar calibration specification is SNP-only", call. = FALSE)
}

#' Estimate a panel-wide called-genotype error level from assigned trios
#'
#' Fits the distribution of Mendelian-inconsistency counts conditional on the
#' pass-1 admission ceiling. Missing loci are marginalized separately for each
#' trio. Only complete, confidently assigned trios may be supplied.
#'
#' @param GenoM public or prepared microhaplotype observations.
#' @param tables pass-1 likelihood tables.
#' @param parents assigned parent data frame or integer parent-index matrix.
#' @param MaxMendelE pass-1 Mendelian-inconsistency admission ceiling.
#' @param confident optional nonmissing logical vector selecting confidently
#'   assigned individuals. No extra threshold is introduced by this function.
#' @return Calibration, observed/predicted count histogram and rebuilt tables.
#' @details This count estimator is used by explicit pedigree likelihood calibration policies. It estimates one scalar
#'   Err from total ME counts. Dropout/foreign panels instead use per-locus
#'   homozygous/heterozygous offspring ME categories to profile the two rates
#'   separately, retaining the ceiling correction. Rates stay
#'   on the low-error branch (each at most 0.25). The likelihood corrects the
#'   explicit ceiling but cannot remove all assignment-score ascertainment or
#'   dependence between trios sharing parents. Profile intervals are composite
#'   likelihood diagnostics, not calibrated coverage guarantees.
#'   Assignment calibration also reports confident trios at the opposite-
#'   homozygote ceiling T; that additional selection is audited, not modeled.
#'   Supplied rates must be positive; fitted rates and profile endpoints are
#'   bounded below by `1e-12` and report when the fit reaches that bound.
#' @keywords internal
#' @noRd
estimate_parentage_error <- function(GenoM, tables, parents, MaxMendelE,
                                      confident = NULL) {
  .require_positive_error_rates(tables)
  if (!is.numeric(MaxMendelE) || length(MaxMendelE) != 1L || is.na(MaxMendelE) ||
      !is.finite(MaxMendelE) || MaxMendelE < 0 || MaxMendelE != floor(MaxMendelE))
    stop("MaxMendelE must be one nonnegative integer", call. = FALSE)
  if (length(tables) && identical(tables[[1L]]$ErrorModel, "dropout_foreign"))
    return(.estimate_parentage_mh_error(GenoM, tables, parents, MaxMendelE, confident))
  ceiling <- as.integer(min(MaxMendelE, length(tables)))
  spec <- .parentage_error_scale_spec(tables)
  obs <- .parentage_error_observations(GenoM, tables, parents, confident = confident)
  nonempty <- obs$n_called > 0L
  obs <- lapply(obs, function(x) if (is.matrix(x)) x[nonempty, , drop = FALSE] else x[nonempty])
  base <- list(model = spec$model, parameter = spec$parameter, initial = spec$initial,
               initial_rates = spec$rates(spec$initial), ceiling = ceiling,
               observations = obs, trios = length(obs$count), tables = tables)
  if (length(obs$count) < 10L || ceiling == 0L) {
    base$status <- if (ceiling == 0L) "unidentifiable_zero_ceiling" else "too_few_trios"
    return(base)
  }
  if (any(obs$count > ceiling))
    stop("a supplied confident trio exceeds the pass-1 MaxMendelE ceiling", call. = FALSE)
  count_probability <- function(theta) {
    t <- spec$rebuild(theta)
    p <- vapply(t, `[[`, 0.0, "p_me_true_trio")
    .parentage_error_count_pmf(p, obs$called, ceiling)
  }
  objective <- function(theta, truncated = TRUE) {
    pmf <- count_probability(theta)
    selected <- pmf[cbind(seq_along(obs$count), obs$count + 1L)]
    if (truncated) selected <- selected / rowSums(pmf)
    if (any(!is.finite(selected) | selected <= 0)) return(Inf)
    -sum(log(selected))
  }
  floor <- .parentage_positive_rate_floor
  fit <- stats::optimize(objective, c(floor, spec$upper), tol = 1e-8)
  endpoints <- c(floor, fit$minimum, spec$upper)
  values <- vapply(endpoints, objective, 0.0)
  theta <- endpoints[which.min(values)]; minimum <- min(values)
  # Diagnostic only: this fit deliberately omits the admission denominator.
  naive <- stats::optimize(function(x) objective(x, FALSE), c(floor, spec$upper), tol = 1e-8)$minimum
  cutoff <- minimum + stats::qchisq(.95, 1L) / 2
  cross <- function(x) objective(x) - cutoff
  lower <- if (theta == floor || cross(floor) <= 0) floor else
    stats::uniroot(cross, c(floor, theta), tol = 1e-7)$root
  upper <- if (theta == spec$upper || cross(spec$upper) <= 0) spec$upper else
    stats::uniroot(cross, c(theta, spec$upper), tol = 1e-7)$root
  assumed <- count_probability(spec$initial); fitted <- count_probability(theta)
  histogram <- data.frame(count = 0:ceiling,
    observed = tabulate(obs$count + 1L, nbins = ceiling + 1L),
    assumed_unconditional = colSums(assumed),
    assumed_given_ceiling = colSums(assumed / rowSums(assumed)),
    fitted_given_ceiling = colSums(fitted / rowSums(fitted)))
  c(base[setdiff(names(base), "tables")], list(status = "fitted", estimate = theta,
    estimated_rates = spec$rates(theta), interval = c(lower = lower, upper = upper),
    positive_rate_floor = floor, at_positive_rate_floor = theta == floor,
    profile_excludes_initial = spec$initial < lower || spec$initial > upper,
    loglik_gain = objective(spec$initial) - minimum,
    naive_estimate = naive, naive_rates = spec$rates(naive),
    histogram = histogram, mean_count = mean(obs$count),
    ceiling_count = sum(obs$count == ceiling),
    boundary_fit = theta %in% c(floor, spec$upper), tables = spec$rebuild(theta)))
}

