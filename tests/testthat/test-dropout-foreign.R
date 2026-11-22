# Independent dense test-only oracles for the explicit dropout/foreign model.
# These helpers intentionally do not call the production emission kernel.

.df_genotypes <- function(k) {
  do.call(rbind, lapply(seq_len(k), function(a) cbind(a, b = a:k)))
}

.df_oracle <- function(p, Dropout, Foreign, q = p) {
  k <- length(p)
  geno <- .df_genotypes(k)
  G <- nrow(geno)
  key <- paste(geno[, 1L], geno[, 2L], sep = ":")
  D <- matrix(0, G, G)
  for (t in seq_len(G)) {
    a <- geno[t, 1L]
    b <- geno[t, 2L]
    if (a == b) {
      D[t, t] <- 1
    } else {
      D[t, match(paste(a, a, sep = ":"), key)] <- Dropout / 2
      D[t, match(paste(a, b, sep = ":"), key)] <- 1 - Dropout
      D[t, match(paste(b, b, sep = ":"), key)] <- Dropout / 2
    }
  }
  Q <- outer(seq_len(k), seq_len(k), function(x, y)
    (1 - Foreign) * (x == y) + Foreign * q[y])
  R <- matrix(0, G, G)
  for (z in seq_len(G)) for (o in seq_len(G)) {
    x <- geno[z, ]; y <- geno[o, ]
    R[z, o] <- if (y[1L] == y[2L]) {
      Q[x[1L], y[1L]] * Q[x[2L], y[1L]]
    } else {
      Q[x[1L], y[1L]] * Q[x[2L], y[2L]] +
        Q[x[1L], y[2L]] * Q[x[2L], y[1L]]
    }
  }
  list(geno = geno, D = D, R = R, K = D %*% R)
}

.df_index <- function(a, b, k) {
  if (a > b) { z <- a; a <- b; b <- z }
  (a - 1L) * k - (a - 1L) * (a - 2L) / 2L + b - a + 1L
}

.df_founder_prior <- function(geno, p) {
  vapply(seq_len(nrow(geno)), function(i) {
    if (geno[i, 1L] == geno[i, 2L]) p[geno[i, 1L]]^2 else
      2 * p[geno[i, 1L]] * p[geno[i, 2L]]
  }, numeric(1))
}

.df_recursive_dense <- function(x, G) {
  if (is.matrix(x) && is.numeric(x) && identical(dim(x), rep(as.integer(G), 2L)))
    return(TRUE)
  if (!is.list(x)) return(FALSE)
  any(vapply(x, .df_recursive_dense, logical(1), G = G))
}

.df_descriptor <- function(oracle, p, observed) {
  k <- length(p)
  d <- vapply(seq_len(k), function(a)
    oracle$K[.df_index(a, a, k), observed], numeric(1))
  r <- vapply(seq_len(k), function(a)
    sum(vapply(seq_len(k), function(b)
      p[b] * oracle$K[.df_index(a, b, k), observed], numeric(1))), numeric(1))
  list(d = d, r = r, h = sum(p * d), m = sum(p * r))
}

.df_pair_states <- function(oracle, p, left, right) {
  tab <- list(k = length(p), G = nrow(oracle$geno), geno = oracle$geno, p = p)
  # The existing dense Jacquard construction is an independent test oracle;
  # only its p/genotype fields are used here.
  basis <- dense_jacquard_basis_oracle(tab)
  vapply(basis, function(x)
    sum(x * outer(oracle$K[, left], oracle$K[, right])), numeric(1))
}

.df_transmission <- function(oracle, p, observed) {
  k <- length(p)
  prior <- .df_founder_prior(oracle$geno, p)
  u <- numeric(k)
  for (g in seq_len(nrow(oracle$geno))) {
    copies <- oracle$geno[g, ]
    for (a in copies)
      u[a] <- u[a] + prior[g] * oracle$K[g, observed] / 2
  }
  u
}

.df_trio <- function(oracle, p, observed_child, observed_dam, observed_sire) {
  k <- length(p)
  prior <- .df_founder_prior(oracle$geno, p)
  m <- as.numeric(crossprod(prior, oracle$K))
  child <- .df_descriptor(oracle, p, observed_child)
  dam <- .df_transmission(oracle, p, observed_dam)
  sire <- .df_transmission(oracle, p, observed_sire)
  h_pp <- 0
  for (a in seq_len(k)) for (b in seq_len(k))
    h_pp <- h_pp + dam[a] * sire[b] * oracle$K[.df_index(a, b, k), observed_child]
  c(H_PP = h_pp, H_D0 = sum(dam * child$r) * m[observed_sire],
    H_0S = sum(sire * child$r) * m[observed_dam],
    H_U = m[observed_child] * m[observed_dam] * m[observed_sire])
}

.df_disjoint <- function(geno, x, y)
  !any(geno[x, ] %in% geno[y, ])

.df_oh <- function(oracle, p) {
  geno <- oracle$geno
  prior <- .df_founder_prior(geno, p)
  m <- as.numeric(crossprod(prior, oracle$K))
  unrelated <- sum(vapply(seq_len(nrow(geno)), function(x)
    sum(vapply(seq_len(nrow(geno)), function(y)
      m[x] * m[y] * .df_disjoint(geno, x, y), numeric(1))), numeric(1)))
  true_po <- 0
  for (g in seq_len(nrow(geno))) {
    transmitted <- if (geno[g, 1L] == geno[g, 2L]) geno[g, 1L] else geno[g, ]
    tx_prob <- if (length(transmitted) == 1L) 1 else .5
    for (a in transmitted) for (b in seq_along(p)) {
      child <- .df_index(a, b, length(p))
      true_po <- true_po + prior[g] * tx_prob * p[b] * sum(vapply(
        seq_len(nrow(geno)), function(x)
          oracle$K[g, x] * sum(vapply(seq_len(nrow(geno)), function(y)
            oracle$K[child, y] * .df_disjoint(geno, x, y), numeric(1))),
        numeric(1)))
    }
  }
  c(true_po = true_po, unrelated = unrelated)
}

.df_me <- function(oracle, p) {
  geno <- oracle$geno
  prior <- .df_founder_prior(geno, p)
  compatible <- function(child, dam, sire) {
    any(vapply(geno[dam, ], function(a)
      any(vapply(geno[sire, ], function(b)
        identical(sort(c(a, b)), unname(geno[child, ])), logical(1))), logical(1)))
  }
  out <- 0
  for (dam in seq_len(nrow(geno))) for (sire in seq_len(nrow(geno)))
    for (od in seq_len(nrow(geno))) for (os in seq_len(nrow(geno)))
      for (oc in seq_len(nrow(geno))) {
        child_mass <- 0
        for (a in geno[dam, ]) for (b in geno[sire, ]) {
          child <- .df_index(a, b, length(p))
          child_mass <- child_mass + .25 * oracle$K[child, oc]
        }
        if (!compatible(oc, od, os))
          out <- out + prior[dam] * prior[sire] * oracle$K[dam, od] *
            oracle$K[sire, os] * child_mass
      }
  out
}

.df_sample_observations <- function(truth, K) {
  vapply(truth, function(x)
    sample.int(nrow(K), 1L, prob = K[x, ]), integer(1))
}

.df_observed_trio_compatible <- function(geno, child, dam, sire) {
  any(vapply(geno[dam, ], function(a)
    any(vapply(geno[sire, ], function(b)
      identical(sort(c(a, b)), unname(geno[child, ])), logical(1))), logical(1)))
}

.df_descriptor_missing <- function(oracle, p, observed) {
  if (observed < 1L)
    return(list(d = rep(1, length(p)), r = rep(1, length(p)), h = 1, m = 1))
  .df_descriptor(oracle, p, observed)
}

.df_transmission_missing <- function(oracle, p, observed) {
  if (observed < 1L) return(p)
  .df_transmission(oracle, p, observed)
}

.df_pair_states_missing <- function(oracle, p, left, right) {
  tab <- list(k = length(p), G = nrow(oracle$geno), geno = oracle$geno, p = p)
  basis <- dense_jacquard_basis_oracle(tab)
  e_left <- if (left < 1L) rep(1, nrow(oracle$K)) else oracle$K[, left]
  e_right <- if (right < 1L) rep(1, nrow(oracle$K)) else oracle$K[, right]
  vapply(basis, function(x)
    sum(x * outer(e_left, e_right)), numeric(1))
}

.df_trio_missing <- function(oracle, p, observed_child, observed_dam,
                             observed_sire) {
  k <- length(p)
  prior <- .df_founder_prior(oracle$geno, p)
  m <- as.numeric(crossprod(prior, oracle$K))
  child <- .df_descriptor_missing(oracle, p, observed_child)
  dam <- .df_transmission_missing(oracle, p, observed_dam)
  sire <- .df_transmission_missing(oracle, p, observed_sire)
  child_emission <- if (observed_child < 1L) rep(1, nrow(oracle$K)) else
    oracle$K[, observed_child]
  h_pp <- sum(vapply(seq_len(k), function(a)
    sum(vapply(seq_len(k), function(b)
      dam[a] * sire[b] * child_emission[.df_index(a, b, k)], numeric(1))),
    numeric(1)))
  dam_m <- if (observed_dam < 1L) 1 else m[observed_dam]
  sire_m <- if (observed_sire < 1L) 1 else m[observed_sire]
  c(H_PP = h_pp, H_D0 = sum(dam * child$r) * sire_m,
    H_0S = sum(sire * child$r) * dam_m,
    H_U = (if (observed_child < 1L) 1 else m[observed_child]) *
      dam_m * sire_m)
}

test_that("dropout_foreign dense composition is stochastic through k=30", {
  for (k in c(1L, 2L, 3L, 4L, 6L, 10L)) {
    p <- seq_len(k) / sum(seq_len(k))
    oracle <- .df_oracle(p, .23, .17)
    tab <- PaPeR:::build_locus_table(
      p, ErrorModel = "dropout_foreign", Dropout = .23, Foreign = .17)
    idx <- seq_len(tab$G) - 1L
    got <- PaPeR:::emission_probabilities_cpp(
      tab, rep(idx, each = tab$G), rep.int(idx, times = tab$G))
    got <- matrix(got, tab$G, tab$G, byrow = TRUE)
    expect_equal(got, oracle$K, tolerance = 1e-13)
    expect_true(all(is.finite(got)) && all(got >= 0))
    expect_equal(rowSums(got), rep(1, tab$G), tolerance = 1e-13)
    expect_equal(tab$mObs, as.numeric(crossprod(.df_founder_prior(oracle$geno, p),
                                                 oracle$K)), tolerance = 1e-13)
  }
  k <- 30L
  tab <- PaPeR:::build_locus_table(
    rep(1 / k, k), ErrorModel = "dropout_foreign", Dropout = .23, Foreign = .17)
  set.seed(81726)
  cells <- cbind(sample.int(tab$G, 250L, replace = TRUE),
                 sample.int(tab$G, 250L, replace = TRUE))
  got <- PaPeR:::emission_probabilities_cpp(
    tab, cells[, 1L] - 1L, cells[, 2L] - 1L)
  expect_true(all(is.finite(got)) && all(got >= 0))
  row_sums <- numeric(tab$G)
  for (truth in seq_len(tab$G) - 1L) {
    row <- PaPeR:::emission_probabilities_cpp(
      tab, rep(truth, tab$G), seq_len(tab$G) - 1L)
    expect_true(all(is.finite(row)) && all(row >= 0))
    row_sums[truth + 1L] <- sum(row)
  }
  expect_equal(row_sums, rep(1, tab$G), tolerance = 1e-12)
  expect_false(.df_recursive_dense(tab, tab$G))
})

test_that("k=2 dropout_foreign has the stated closed-form matrix", {
  delta <- .4
  rho <- .2
  a <- 1 - rho / 2
  b <- rho / 2
  ab <- a * b
  aa <- a^2
  bb <- b^2
  middle <- a^2 + b^2
  expected <- rbind(
    c(aa, 2 * ab, bb),
    c((1 - delta) * ab + delta * middle / 2,
      (1 - delta) * middle + 2 * delta * ab,
      (1 - delta) * ab + delta * middle / 2),
    c(bb, 2 * ab, aa))
  tab <- PaPeR:::build_locus_table(
    c(.3, .7), ErrorModel = "dropout_foreign", Dropout = delta,
    Foreign = rho, ForeignQ = c(.5, .5))
  idx <- seq_len(tab$G) - 1L
  got <- matrix(PaPeR:::emission_probabilities_cpp(
    tab, rep(idx, each = tab$G), rep.int(idx, times = tab$G)),
    nrow = tab$G, byrow = TRUE)
  expect_equal(got, expected, tolerance = 1e-14)
  legacy <- PaPeR:::build_locus_table(c(.3, .7), Err = .4)
  legacy_idx <- seq_len(legacy$G) - 1L
  legacy_matrix <- matrix(PaPeR:::emission_probabilities_cpp(
    legacy, rep(legacy_idx, each = legacy$G),
    rep.int(legacy_idx, times = legacy$G)), nrow = legacy$G, byrow = TRUE)
  expect_gt(max(abs(got - legacy_matrix)), 1e-8)
})

test_that("q=p gives the independent transmitted-state identity s=p", {
  for (p in list(c(.97, .02, .01), c(.01, .09, .2, .7),
                c(.001, .009, .09, .2, .7))) {
    rho <- .37
    Q <- outer(seq_along(p), seq_along(p), function(x, y)
      (1 - rho) * (x == y) + rho * p[y])
    s <- vapply(seq_along(p), function(y) sum(p * Q[, y]), numeric(1))
    expect_equal(s, p, tolerance = 1e-15)
  }
})

test_that("positive-rate component limits and far-error scaling match enumeration", {
  p <- c(.2, .3, .5)
  identity <- PaPeR:::build_locus_table(
    p, ErrorModel = "dropout_foreign", Dropout = 1e-12, Foreign = 1e-12)
  positive_oracle <- .df_oracle(p, 1e-12, 1e-12)
  expect_equal(identity$mObs, as.numeric(crossprod(
    .df_founder_prior(positive_oracle$geno, p), positive_oracle$K)), tolerance = 1e-14)
  expect_equal(PaPeR:::emission_probabilities_cpp(
    identity, 0:5, 0:5), diag(.df_oracle(p, 1e-12, 1e-12)$K), tolerance = 1e-14)
  for (pars in list(c(1, 1e-12), c(1e-12, .2), c(.3, .2))) {
    oracle <- .df_oracle(p, pars[1L], pars[2L])
    tab <- PaPeR:::build_locus_table(
      p, ErrorModel = "dropout_foreign", Dropout = pars[1L], Foreign = pars[2L])
    cells <- data.frame(true = rep(0:5, each = 6L),
                        observed = rep.int(0:5, times = 6L))
    got <- PaPeR:::emission_probabilities_cpp(
      tab, cells$true, cells$observed)
    expect_equal(got, as.vector(t(oracle$K)), tolerance = 1e-14)
  }
  q <- c(.6, .1, .3)
  oracle_q <- .df_oracle(p, .3, .2, q)
  tab_q <- PaPeR:::build_locus_table(
    p, ErrorModel = "dropout_foreign", Dropout = .3, Foreign = .2,
    ForeignQ = q)
  idx <- seq_len(tab_q$G) - 1L
  got_q <- PaPeR:::emission_probabilities_cpp(
    tab_q, rep(idx, each = tab_q$G), rep.int(idx, times = tab_q$G))
  expect_equal(matrix(got_q, tab_q$G, tab_q$G, byrow = TRUE),
               oracle_q$K, tolerance = 1e-14)
  low <- PaPeR:::build_locus_table(
    c(.7, .3), ErrorModel = "dropout_foreign", Dropout = 1e-12, Foreign = .01)
  high <- PaPeR:::build_locus_table(
    c(.7, .3), ErrorModel = "dropout_foreign", Dropout = 1e-12, Foreign = .02)
  one_low <- PaPeR:::emission_probabilities_cpp(low, 0L, 1L)
  one_high <- PaPeR:::emission_probabilities_cpp(high, 0L, 1L)
  two_low <- PaPeR:::emission_probabilities_cpp(low, 0L, 2L)
  two_high <- PaPeR:::emission_probabilities_cpp(high, 0L, 2L)
  expect_true(all(c(one_low, one_high, two_low, two_high) > 0))
  expect_equal(one_high / one_low, 2, tolerance = .02)
  expect_equal(two_high / two_low, 4, tolerance = 1e-12)
})

test_that("OH and true-trio ME calibration match dense enumeration", {
  p <- c(.2, .3, .5)
  oracle <- .df_oracle(p, .2, .1)
  tab <- PaPeR:::build_locus_table(
    p, ErrorModel = "dropout_foreign", Dropout = .2, Foreign = .1)
  oh <- .df_oh(oracle, p)
  expect_equal(tab$p_oh_true_po, oh[["true_po"]], tolerance = 1e-12)
  expect_equal(tab$p_oh_unrelated, oh[["unrelated"]], tolerance = 1e-12)
  expect_equal(tab$p_me_true_trio, .df_me(oracle, p), tolerance = 1e-12)
})

test_that("latent dense Monte Carlo calibrates OH and true-trio ME scalars", {
  p <- c(.2, .3, .5)
  oracle <- .df_oracle(p, .2, .1)
  tab <- PaPeR:::build_locus_table(
    p, ErrorModel = "dropout_foreign", Dropout = .2, Foreign = .1)
  geno <- oracle$geno
  prior <- .df_founder_prior(geno, p)
  N <- 6000L
  set.seed(8172601)

  po_parent <- sample.int(nrow(geno), N, replace = TRUE, prob = prior)
  po_parent_geno <- geno[po_parent, , drop = FALSE]
  po_tx <- ifelse(po_parent_geno[, 1L] == po_parent_geno[, 2L],
                  po_parent_geno[, 1L],
                  ifelse(runif(N) < .5, po_parent_geno[, 1L],
                         po_parent_geno[, 2L]))
  po_other <- sample.int(length(p), N, replace = TRUE, prob = p)
  po_child <- vapply(seq_len(N), function(i)
    as.integer(.df_index(po_tx[i], po_other[i], length(p))), integer(1))
  po_parent_obs <- .df_sample_observations(po_parent, oracle$K)
  po_child_obs <- .df_sample_observations(po_child, oracle$K)
  po_event <- vapply(seq_len(N), function(i)
    .df_disjoint(geno, po_parent_obs[i], po_child_obs[i]), logical(1))

  unrelated_a <- sample.int(nrow(geno), N, replace = TRUE, prob = prior)
  unrelated_b <- sample.int(nrow(geno), N, replace = TRUE, prob = prior)
  unrelated_a_obs <- .df_sample_observations(unrelated_a, oracle$K)
  unrelated_b_obs <- .df_sample_observations(unrelated_b, oracle$K)
  unrelated_event <- vapply(seq_len(N), function(i)
    .df_disjoint(geno, unrelated_a_obs[i], unrelated_b_obs[i]), logical(1))

  me_dam <- sample.int(nrow(geno), N, replace = TRUE, prob = prior)
  me_sire <- sample.int(nrow(geno), N, replace = TRUE, prob = prior)
  me_dam_geno <- geno[me_dam, , drop = FALSE]
  me_sire_geno <- geno[me_sire, , drop = FALSE]
  me_dam_tx <- ifelse(runif(N) < .5, me_dam_geno[, 1L], me_dam_geno[, 2L])
  me_sire_tx <- ifelse(runif(N) < .5, me_sire_geno[, 1L], me_sire_geno[, 2L])
  me_child <- vapply(seq_len(N), function(i)
    as.integer(.df_index(me_dam_tx[i], me_sire_tx[i], length(p))), integer(1))
  me_dam_obs <- .df_sample_observations(me_dam, oracle$K)
  me_sire_obs <- .df_sample_observations(me_sire, oracle$K)
  me_child_obs <- .df_sample_observations(me_child, oracle$K)
  me_event <- !vapply(seq_len(N), function(i)
    .df_observed_trio_compatible(geno, me_child_obs[i], me_dam_obs[i],
                                 me_sire_obs[i]), logical(1))

  expect_gt(sum(po_event), 0L)
  expect_gt(sum(unrelated_event), 0L)
  expect_gt(sum(me_event), 0L)
  expect_lt(sum(po_event), N)
  expect_lt(sum(unrelated_event), N)
  expect_lt(sum(me_event), N)
  tol <- 6 * sqrt(.25 / N) # six worst-case binomial standard errors
  expect_lt(abs(mean(po_event) - tab$p_oh_true_po), tol)
  expect_lt(abs(mean(unrelated_event) - tab$p_oh_unrelated), tol)
  expect_lt(abs(mean(me_event) - tab$p_me_true_trio), tol)
})

test_that("public routing, missingness, provenance, and compact schema hold", {
  p <- list(L1 = c(`1` = .2, `2` = .8), L2 = c(`1` = .4, `2` = .3, `3` = .3))
  tables <- build_tables(p, ErrorModel = "dropout_foreign",
                         Dropout = .2, Foreign = .1)
  expect_true(all(vapply(tables, function(x)
    identical(x$schema_version, "compact-locus-v2"), logical(1))))
  expect_true(all(vapply(tables, function(x)
    !any(c("Err", "ErrProfile", "Z_hom", "Z_het") %in% names(x)), logical(1))))
  expect_false(any(vapply(tables, .df_recursive_dense,
                          logical(1), G = tables[[1L]]$G)))
  GenoM <- matrix(c("1/1", "1/2", "2/2"), ncol = 1L,
                  dimnames = list(c("a", "b", "c"), "L1"))
  pair <- .calc_pair_LL(GenoM, "a", "b", tables[1L])
  trio <- .calc_trio_LL(GenoM, "a", "b", "c", tables[1L])
  full <- .calc_pair_LL(GenoM, "a", "b", tables[1L], Complex = "full")
  expect_true(all(is.finite(pair)) && all(is.finite(full)) && all(is.finite(trio)))
  missing <- matrix(c("0", "1/1"), ncol = 1L,
                    dimnames = list(c("missing", "called"), "L1"))
  expect_true(all(is.finite(.calc_pair_LL(missing, "missing", "called",
                                         tables[1L]))))
  expect_equal(tables[[1L]]$ForeignQMode, "p")
  expect_error(build_tables(p, ErrorModel = "dropout_foreign", Dropout = 0),
               "Dropout and Foreign")
  expect_error(build_tables(p, ErrorModel = "dropout_foreign", Dropout = 0,
                            Foreign = 0, Err = 0), "Err and ErrProfile")
  expect_error(build_tables(p, Dropout = 0), "require ErrorModel")
  expect_error(PaPeR:::build_locus_table(c(.5, .5),
    ErrorModel = "dropout_foreign", Dropout = 1e-12, Foreign = 1e-12,
    ForeignQ = c(.4, .4)), "normalized")
  expect_error(PaPeR:::build_locus_table(c(`A` = .5, `B` = .5),
    ErrorModel = "dropout_foreign", Dropout = 1e-12, Foreign = 1e-12,
    ForeignQ = c(B = .5, A = .5)), "names")
  expect_error(PaPeR:::build_locus_table(c(.5, .5),
    ErrorModel = "legacy_sequoia", Dropout = 0), "require ErrorModel")
  q <- c(.7, .3)
  supplied <- PaPeR:::build_locus_table(
    c(.7, .3), ErrorModel = "dropout_foreign", Dropout = 1e-12,
    Foreign = 1e-12, ForeignQ = q)
  expect_identical(supplied$ForeignQMode, "supplied")
  expect_equal(supplied$mObs, PaPeR:::build_locus_table(
    c(.7, .3), ErrorModel = "dropout_foreign", Dropout = 1e-12,
    Foreign = 1e-12)$mObs, tolerance = 0)
})

test_that("allele permutation, monomorphic loci, and bounded Monte Carlo agree", {
  p <- c(.15, .35, .5)
  q <- c(.6, .1, .3)
  tab <- PaPeR:::build_locus_table(
    p, ErrorModel = "dropout_foreign", Dropout = .13, Foreign = .21,
    ForeignQ = q)
  perm <- c(3L, 1L, 2L)
  perm_tab <- PaPeR:::build_locus_table(
    p[perm], ErrorModel = "dropout_foreign", Dropout = .13, Foreign = .21,
    ForeignQ = q[perm])
  map <- vapply(seq_len(tab$G), function(i) {
    z <- sort(match(tab$geno[i, ], perm))
    as.integer(.df_index(z[1L], z[2L], length(p)) - 1L)
  }, integer(1))
  idx <- seq_len(tab$G) - 1L
  map_r <- seq_len(tab$G)
  expect_equal(PaPeR:::emission_probabilities_cpp(perm_tab,
    map[rep(map_r, each = tab$G)], map[rep.int(map_r, times = tab$G)]),
    PaPeR:::emission_probabilities_cpp(tab,
      rep(idx, each = tab$G), rep.int(idx, times = tab$G)), tolerance = 1e-14)
  mono <- PaPeR:::build_locus_table(1, ErrorModel = "dropout_foreign",
                                        Dropout = 1, Foreign = 1)
  expect_equal(mono$G, 1L)
  expect_equal(mono$mObs, 1)

  set.seed(81726)
  oracle <- .df_oracle(c(.2, .3, .5), .2, .1)
  true <- sample.int(nrow(oracle$geno), 2000L, replace = TRUE,
                     prob = .df_founder_prior(oracle$geno, c(.2, .3, .5)))
  observed <- vapply(true, function(x)
    sample.int(nrow(oracle$geno), 1L, prob = oracle$K[x, ]), integer(1))
  empirical <- tabulate(observed, nbins = nrow(oracle$geno)) / length(observed)
  expect_true(sum(empirical > 0) > 3L)
  expect_equal(empirical, as.numeric(crossprod(
    .df_founder_prior(oracle$geno, c(.2, .3, .5)), oracle$K)), tolerance = .06)
})
