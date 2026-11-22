# Two observable Mendelian-inconsistency categories identify the two error
# mechanisms without treating an observed incompatible trio as an error label.
# Integrate latent true genotypes exactly using the established emission kernel.
.parentage_mh_me_categories <- function(table) {
  k <- table$k; g <- table$geno; G <- nrow(g); p <- table$p
  index <- matrix(0L, k, k)
  for (j in seq_len(G)) index[g[j, 1L], g[j, 2L]] <- index[g[j, 2L], g[j, 1L]] <- j
  E <- matrix(emission_probabilities_cpp(table, rep(seq_len(G) - 1L, G),
                       rep(seq_len(G) - 1L, each = G)), G, G)
  r <- t(vapply(seq_len(k), function(a) colSums(E[index[a, ], , drop = FALSE] * p), numeric(G)))
  contains <- vapply(seq_len(k), function(a) as.numeric(g[, 1L] == a | g[, 2L] == a), numeric(G))
  containment <- r %*% contains
  prior <- p[g[, 1L]] * p[g[, 2L]] * ifelse(g[, 1L] == g[, 2L], 1, 2)
  marginal <- colSums(E * prior); compatible <- numeric(G); pp <- outer(p, p)
  for (j in seq_len(G)) {
    u <- g[j, 1L]; v <- g[j, 2L]
    K <- outer(containment[, u], containment[, v])
    if (u != v) K <- K + t(K) - outer(r[, j], r[, j])
    compatible[j] <- sum(pp * matrix(E[as.vector(index), j], k, k) * K)
  }
  me <- marginal - compatible
  if (any(me < -1e-11)) stop('invalid multiallelic ME category probability')
  me <- pmax(me, 0)
  c(homozygous = sum(me[g[, 1L] == g[, 2L]]),
    heterozygous = sum(me[g[, 1L] != g[, 2L]]))
}

.parentage_mh_rebuild <- function(tables, rates) {
  ans <- lapply(tables, function(t) {
    x <- build_locus_table(t$p, ErrorModel = 'dropout_foreign',
      Dropout = unname(rates['Dropout']), Foreign = unname(rates['Foreign']), ForeignQ = t$ForeignQ)
    x$ForeignQMode <- t$ForeignQMode; x
  })
  names(ans) <- names(tables); ans
}

.parentage_bernstein <- function(x, degree) stats::dbinom(0:degree, degree, x)

# Emission is degree (1,2) in (D,F); the three-emission probability is degree
# (3,6). Interpolation in a bounded Bernstein basis compiles that exact small
# polynomial once, rather than reconstructing tables at every profile point.
.parentage_mh_category_model <- function(tables, upper = .25, cache = NULL) {
  # These are the complete inputs read by rebuild(). Error rates are the
  # polynomial variables, so changing only them cannot invalidate this basis.
  tables <- lapply(tables,function(t)t[c('p','ForeignQ','ForeignQMode')])
  key <- list(upper=upper,tables=tables)
  if (!is.null(cache) && identical(cache$key,key)) {
    cache$hits <- cache$hits + 1L
    return(cache$probability)
  }
  started <- proc.time()[['elapsed']]

  d <- seq(.parentage_positive_rate_floor / upper, 1, length.out = 4L)
  f <- seq(.parentage_positive_rate_floor / upper, 1, length.out = 7L)
  grid <- expand.grid(Dropout = d, Foreign = f)
  basis <- function(rate) as.vector(outer(.parentage_bernstein(rate[1L], 3L),
                                         .parentage_bernstein(rate[2L], 6L)))
  B <- t(apply(grid, 1L, basis)); values <- matrix(0, nrow(grid), 2L * length(tables))
  for (i in seq_len(nrow(grid))) {
    t <- .parentage_mh_rebuild(tables, setNames(upper * unlist(grid[i, ]), c('Dropout', 'Foreign')))
    v <- vapply(t, .parentage_mh_me_categories, numeric(2L))
    stopifnot(max(abs(colSums(v) - vapply(t, `[[`, 0.0, 'p_me_true_trio'))) < 1e-10)
    values[i, ] <- as.vector(v)
  }
  coefficients <- solve(B, values)
  probability <- .parentage_mh_category_evaluator(coefficients,upper,tables)
  if (!is.null(cache)) {
    cache$builds <- cache$builds + 1L
    cache$seconds <- cache$seconds + proc.time()[['elapsed']] - started
    bytes <- as.numeric(object.size(list(key,coefficients,tables))) + 8192
    if (bytes <= cache$limit) {
      cache$key <- key;cache$probability <- probability;cache$bytes <- bytes
    } else {
      cache$key <- cache$probability <- NULL;cache$bytes <- 0;cache$declines <- cache$declines + 1L
    }
  }
  probability
}
.parentage_mh_category_evaluator <- function(coefficients,upper,tables) {
  force(coefficients);force(upper);force(tables)
  basis <- function(rate) as.vector(outer(.parentage_bernstein(rate[1L],3L),.parentage_bernstein(rate[2L],6L)))
  function(rates) {
    if (length(rates) != 2L || any(!is.finite(rates)) || any(rates <= 0))
      stop('multiallelic rates must be two finite probabilities greater than zero; supply small positive values')
    # The estimator searches the existing low-error branch, but the supplied
    # rates used for diagnostic likelihood gain may be outside that branch.
    # Evaluate those valid physical rates directly; never feed x>1 to dbinom's
    # Bernstein basis. This does not expand the optimization/profile domain.
    if (any(rates < 0 | rates > upper))
      return(vapply(.parentage_mh_rebuild(tables,
        setNames(rates, c('Dropout', 'Foreign'))), .parentage_mh_me_categories, numeric(2L)))
    p <- matrix(as.vector(basis(rates / upper) %*% coefficients), 2L, length(tables))
    if (any(!is.finite(p)) || any(p < -1e-10) || any(colSums(p) > 1 + 1e-10))
      stop('invalid compiled multiallelic category probability')
    pmax(p, 0)
  }
}

.estimate_parentage_mh_error <- function(GenoM, tables, parents, MaxMendelE, confident = NULL) {
  if (!all(vapply(tables, function(t) identical(t$ErrorModel, 'dropout_foreign'), logical(1))))
    stop('multiallelic calibration requires dropout_foreign tables')
  .require_positive_error_rates(tables)
  floor <- .parentage_positive_rate_floor
  supplied <- c(Dropout = tables[[1L]]$Dropout, Foreign = tables[[1L]]$Foreign)
  if (length(MaxMendelE) != 1L || !is.finite(MaxMendelE) || MaxMendelE < 0 || MaxMendelE != floor(MaxMendelE))
    stop('MaxMendelE must be one nonnegative integer')
  ceiling <- as.integer(min(MaxMendelE, length(tables)))
  obs <- .parentage_error_observations(GenoM, tables, parents, confident)
  use <- obs$n_called > 0L
  obs <- lapply(obs, function(x) if (is.matrix(x)) x[use, , drop = FALSE] else x[use])
  base <- list(model = 'dropout_foreign', parameter = c('Dropout', 'Foreign'),
    initial_rates = supplied, ceiling = ceiling, trios = length(obs$count), observations = obs,
    likelihood = 'ceiling_truncated_ME_category_composite',
    interval_method = '95% two-sided nuisance-profile; shared-trio coverage approximate')
  if (length(obs$count) < 10L || ceiling == 0L) {
    base$status <- if (ceiling == 0L) 'unidentifiable_zero_ceiling' else 'too_few_trios'
    return(base)
  }
  if (any(obs$count > ceiling)) stop('confident trio exceeds pass-1 MaxMendelE ceiling')
  g <- .prepare_for_tables(GenoM, tables)
  homo <- g[obs$child, seq(1L, ncol(g), 2L), drop = FALSE] ==
    g[obs$child, seq(2L, ncol(g), 2L), drop = FALSE]
  categories <- obs$incompatible * ifelse(homo, 1L, 2L)
  # Retain per-locus counts and called masks: no genotype error is inferred
  # directly from its Mendelian incompatibility category.
  counts <- rbind(colSums(obs$called & categories == 0L),
                  colSums(obs$called & categories == 1L),
                  colSums(obs$called & categories == 2L))
  probability <- .parentage_mh_category_model(tables,cache=attr(tables,'pap_category_cache'))
  objective <- function(rates) {
    if (length(rates) != 2L || any(!is.finite(rates)) || any(rates <= 0 | rates > 1)) return(Inf)
    p <- probability(rates); q <- rbind(1 - colSums(p), p)
    if (any(!is.finite(q))) return(Inf)
    if (any(q[counts > 0L] <= 0)) return(Inf)
    numerator <- sum(counts[counts > 0L] * log(q[counts > 0L]))
    selection <- rowSums(.parentage_error_count_pmf(colSums(p), obs$called, ceiling))
    if (any(!is.finite(selection)) || any(selection <= 0)) return(Inf)
    -numerator + sum(log(selection))
  }
  starts <- list(pmax(floor, pmin(supplied, .25)), c(.01, .01), c(.1, .001))
  fits <- lapply(starts, function(s) tryCatch({
    if (!is.finite(objective(s))) return(NULL)
    stats::nlminb(s, objective, lower = rep(floor, 2L), upper = c(.25, .25),
      control = list(rel.tol = 1e-10, x.tol = 1e-9, iter.max = 300L))
  }, error = function(e) NULL))
  fits <- Filter(function(x) !is.null(x) && is.finite(x$objective) &&
    all(is.finite(x$par)), fits)
  if (!length(fits)) { base$status <- 'fit_not_converged';return(base) }
  fit <- fits[[which.min(vapply(fits, `[[`, 0.0, 'objective'))]]
  rates <- setNames(fit$par, names(supplied)); minimum <- fit$objective
  if (!is.finite(minimum) || fit$convergence != 0L) {
    base$status <- 'fit_not_converged';return(base)
  }
  cutoff <- minimum + stats::qchisq(.95, 1L) / 2
  interval <- matrix(NA_real_, 2L, 2L, dimnames = list(names(supplied), c('lower', 'upper')))
  for (j in 1:2) {
    profile <- function(x) {
      eval <- function(y) {r <- rates;r[j] <- x;r[3L-j] <- y;objective(r)}
      opt <- stats::optimize(eval, c(floor, .25), tol = 1e-9)
      min(opt$objective, eval(floor), eval(.25)) - cutoff
    }
    interval[j, 1L] <- if (rates[j] == floor || profile(floor) <= 0) floor else
      stats::uniroot(profile, c(floor, rates[j]), tol = 1e-7)$root
    interval[j, 2L] <- if (rates[j] == .25 || profile(.25) <= 0) .25 else
      stats::uniroot(profile, c(rates[j], .25), tol = 1e-7)$root
  }
  c(base, list(status = 'fitted', estimated_rates = rates, interval = interval,
    loglik_gain = objective(supplied) - minimum, category_counts = counts,
    optimizer_convergence = fit$convergence, positive_rate_floor = floor,
    at_positive_rate_floor = rates == floor))
}
