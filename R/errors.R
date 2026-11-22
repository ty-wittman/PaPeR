# One error model and strictly positive, panel-wide rates on every public path.
.require_positive_error_rates <- function(tables) {
  for (i in seq_along(tables)) {
    t <- tables[[i]]
    fields <- if (identical(t$ErrorModel, 'dropout_foreign')) c('Dropout','Foreign') else 'Err'
    for (field in fields) {
      value <- t[[field]]
      if (length(value) != 1L || !is.numeric(value) || is.complex(value) ||
          !is.finite(value) || value <= 0 || value > 1)
        stop(field, ' must be strictly positive at locus ', i,
          '; zero error/dropout/foreign rates are unsupported. Supply a small positive value.', call. = FALSE)
    }
  }
  models <- vapply(tables, function(t) t$ErrorModel %||% 'legacy_sequoia', '')
  if (length(unique(models)) > 1L)
    stop('A panel must use one error model and panel-wide rates; mixed error models are unsupported.', call. = FALSE)
  if (length(tables)) for (field in if (models[1L] == 'dropout_foreign') c('Dropout','Foreign') else 'Err') {
    values <- vapply(tables, function(t) as.numeric(t[[field]]), 0.0)
    if (length(unique(values)) > 1L)
      stop('Panel-wide ', field, ' required; per-locus error rates are unsupported.', call. = FALSE)
  }
  invisible(tables)
}

.parentage_positive_rate_floor <- 1e-12

.validate_Err <- function(Err, n_loci) {
  if (!is.numeric(Err) || is.object(Err) || length(Err) == 0L)
    stop("Err must be one numeric panel-wide scalar")
  if (length(Err) != 1L)
    stop("Err must be one panel-wide scalar; per-locus Err is unsupported", call. = FALSE)
  if (anyNA(Err) || any(!is.finite(Err)) || any(Err <= 0) || any(Err > 1))
    stop("Err must be strictly positive and finite, at most 1; supply a small positive value")
  rep(as.numeric(Err), n_loci)
}

.profile_from_matrix <- function(x) {
  if (!is.numeric(x) || !is.matrix(x) || !identical(dim(x), c(3L, 3L)))
    stop("the ErrProfile matrix adapter must be a numeric 3 x 3 matrix")
  if (anyNA(x) || any(!is.finite(x)) || any(x < 0) || any(x > 1))
    stop("ErrProfile matrix entries must be finite, non-missing, and between 0 and 1")
  tol <- 1e-12
  paired <- c(x[1, 2] - x[3, 2], x[2, 1] - x[2, 3],
              x[1, 3] - x[3, 1], x[1, 1] - x[3, 3])
  if (any(abs(paired) > tol))
    stop("ErrProfile matrix must be allele-label symmetric")
  stochastic <- all(abs(rowSums(x) - 1) <= tol)
  implied <- NULL
  if (stochastic) {
    implied <- 1 - diag(x)
    if (max(implied) - min(implied) > tol)
      stop("a row-stochastic ErrProfile matrix must imply one common Err")
    implied <- mean(implied)
  }
  list(
    weights = c(w_Het = x[1, 2], w_Hom = x[2, 1], w_Far = x[1, 3]),
    implied = implied
  )
}

.normalize_profile <- function(ErrProfile, Err) {
  if (is.null(ErrProfile)) {
    E <- Err[1L]
    return(list(
      weights = c(w_Het = E - (E / 2)^2,
                  w_Hom = E / 2,
                  w_Far = (E / 2)^2),
      is_default = TRUE,
      implied = NULL
    ))
  }
  if (is.matrix(ErrProfile)) {
    ans <- .profile_from_matrix(ErrProfile)
    ans$is_default <- FALSE
    return(ans)
  }
  if (is.list(ErrProfile)) ErrProfile <- unlist(ErrProfile, use.names = TRUE)
  wanted <- c("w_Het", "w_Hom", "w_Far")
  if (!is.numeric(ErrProfile) || length(ErrProfile) != 3L ||
      is.null(names(ErrProfile)) || !setequal(names(ErrProfile), wanted))
    stop("ErrProfile must contain exactly w_Het, w_Hom, and w_Far")
  ErrProfile <- ErrProfile[wanted]
  if (anyNA(ErrProfile) || any(!is.finite(ErrProfile)) ||
      any(ErrProfile < 0) || any(ErrProfile > 1))
    stop("ErrProfile weights must be finite, non-missing, and between 0 and 1")
  list(weights = structure(as.numeric(ErrProfile), names = wanted),
       is_default = FALSE, implied = NULL)
}

#' Normalize a scalar microhaplotype genotype-error specification
#'
#' \code{Err} is the probability that an observed, called unordered diploid
#' genotype differs from its true genotype. One positive panel-wide scalar
#' is replicated across loci; rate vectors are unsupported.
#'
#' @param Err one positive panel-wide scalar. Default \code{1e-4}.
#' @param n_loci number of loci.
#' @param ErrProfile optional panel-wide named weights \code{w_Het}, \code{w_Hom},
#'   and \code{w_Far}, or an allele-label-symmetric Sequoia-style 3 x 3 matrix.
#' @return list with per-locus \code{Err}, canonical \code{ErrProfile}, and
#'   \code{profile_default}.
#' @details With total error \code{E}, the default raw weights are
#'   \code{w_Het = E - (E/2)^2}, \code{w_Hom = E/2}, and
#'   \code{w_Far = (E/2)^2}. Each wrong target genotype receives its class
#'   weight, then each true-genotype row is normalized so its off-diagonal sum
#'   is exactly \code{E}. In the 3 x 3 adapter, rows are true AA, AB, BB and
#'   columns are observed AA, AB, BB. A row-stochastic adapter implies
#'   \code{Err}; an explicitly supplied \code{Err} must agree.
#' @keywords internal
#' @noRd
expand_errors <- function(Err = 1e-4, n_loci, ErrProfile = NULL) {
  Err_was_missing <- missing(Err)
  if (!is.numeric(n_loci) || length(n_loci) != 1L || is.na(n_loci) ||
      !is.finite(n_loci) || n_loci < 1 || n_loci != as.integer(n_loci))
    stop("n_loci must be one positive integer")
  n_loci <- as.integer(n_loci)
  totals <- .validate_Err(Err, n_loci)
  profile <- .normalize_profile(ErrProfile, totals)
  if (!is.null(profile$implied)) {
    if (!Err_was_missing && any(abs(totals - profile$implied) > 1e-12))
      stop("Err conflicts with the error rate implied by the row-stochastic ErrProfile matrix")
    if (any(profile$implied <= 0)) stop("Err must be strictly positive; supply a small positive value", call. = FALSE)
    if (Err_was_missing) totals[] <- profile$implied
  }
  list(Err = totals, ErrProfile = profile$weights,
       profile_default = profile$is_default)
}

.validate_dropout_foreign_rates <- function(x, n_loci, name) {
  if (!is.numeric(x) || is.object(x) || !is.null(dim(x)) || length(x) == 0L)
    stop(name, " must be one numeric panel-wide scalar",
         call. = FALSE)
  if (length(x) != 1L)
    stop(name, " must be one panel-wide scalar; per-locus rates are unsupported", call. = FALSE)
  if (anyNA(x) || any(!is.finite(x)) || any(x <= 0) || any(x > 1))
    stop(name, " must be strictly positive and finite, at most 1; supply a small positive value",
         call. = FALSE)
  rep(as.numeric(x), n_loci)
}

.validate_foreign_q <- function(q, p, locus = "locus") {
  if (!is.numeric(q) || is.object(q) || !is.null(dim(q)) || length(q) != length(p) ||
      anyNA(q) || any(!is.finite(q)) || any(q < 0))
    stop("ForeignQ for ", locus,
         " must be a finite, nonnegative vector matching the allele universe",
         call. = FALSE)
  if (!is.null(names(q))) {
    expected <- if (is.null(names(p))) as.character(seq_along(p)) else names(p)
    if (!identical(as.character(names(q)), as.character(expected)))
      stop("ForeignQ names for ", locus,
           " must exactly match the declared allele universe in order",
           call. = FALSE)
  }
  total <- sum(q)
  if (!is.finite(total) || total <= 0 || abs(total - 1) > 1e-12)
    stop("ForeignQ for ", locus,
         " must be normalized to sum to one (within 1e-12)", call. = FALSE)
  q_names <- names(q)
  q <- as.numeric(q)
  if (!is.null(names(p))) names(q) <- names(p)
  else if (!is.null(q_names)) names(q) <- q_names
  q
}

.validate_dropout_foreign_spec <- function(Dropout, Foreign, ForeignQ,
                                            freqs, n_loci = length(freqs)) {
  dropout <- .validate_dropout_foreign_rates(Dropout, n_loci, "Dropout")
  foreign <- .validate_dropout_foreign_rates(Foreign, n_loci, "Foreign")
  if (is.null(ForeignQ)) {
    q <- lapply(freqs, function(p) {
      p <- as.numeric(p)
      p / sum(p)
    })
    q_mode <- rep("p", n_loci)
  } else {
    if (n_loci == 1L && is.numeric(ForeignQ)) {
      q <- list(ForeignQ)
    } else if (is.list(ForeignQ) && length(ForeignQ) == n_loci) {
      q <- ForeignQ
    } else {
      stop("ForeignQ must be one vector for one locus or an exact per-locus list",
           call. = FALSE)
    }
    q <- Map(.validate_foreign_q, q, freqs,
             if (is.null(names(freqs))) seq_len(n_loci) else names(freqs))
    q_mode <- rep("supplied", n_loci)
  }
  list(Dropout = dropout, Foreign = foreign, ForeignQ = q,
       ForeignQMode = q_mode)
}
