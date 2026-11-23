gkey <- function(x) paste(sort(x), collapse = ":")
gidx_error <- function(a, b, k) {
  if (a > b) { z <- a; a <- b; b <- z }
  (a - 1) * k - (a - 1) * (a - 2) / 2 + (b - a) + 1
}

class_of_transition <- function(true, observed) {
  overlap <- sum(pmin(tabulate(true, nbins = max(c(true, observed))),
                      tabulate(observed, nbins = max(c(true, observed)))))
  if (2L - overlap == 2L) return("w_Far")
  if (observed[1] == observed[2]) "w_Hom" else "w_Het"
}

# Independent dense oracle retained only in tests.  Production tables must not
# contain this G x G object; the compact C++ kernel is queried cellwise below.
dense_error_oracle <- function(tab) {
  E <- matrix(0, tab$G, tab$G)
  for (true_g in seq_len(tab$G)) {
    for (observed_g in seq_len(tab$G)) {
      if (true_g == observed_g) {
        E[true_g, observed_g] <- 1 - tab$Err
      } else if (tab$Err > 0) {
        cls <- class_of_transition(tab$geno[true_g, ], tab$geno[observed_g, ])
        Z <- if (tab$geno[true_g, 1L] == tab$geno[true_g, 2L]) {
          tab$Z_hom
        } else {
          tab$Z_het
        }
        E[true_g, observed_g] <- tab$Err * tab$ErrProfile[[cls]] / Z
      }
    }
  }
  E
}

implicit_emission_matrix <- function(tab) {
  idx <- seq_len(tab$G) - 1L
  matrix(
    PaPeR:::emission_probabilities_cpp(
      tab, rep(idx, each = tab$G), rep.int(idx, times = tab$G)),
    nrow = tab$G, ncol = tab$G, byrow = TRUE
  )
}

test_that("omitted Err is the documented Sequoia default", {
  freqs <- list(c(0.6, 0.4), c(0.5, 0.3, 0.2))
  expect_equal(build_tables(freqs), build_tables(freqs, Err = 1e-4))
})

test_that("Err is panel-wide and rate vectors are unsupported", {
  freqs <- list(c(0.6, 0.4), c(0.5, 0.3, 0.2))
  expect_error(build_tables(freqs, Err = rep(0.01, length(freqs))), "panel-wide")
  expect_error(build_tables(freqs, Err = c(.01,.02)), "panel-wide")
})

test_that("biallelic default is exactly the Sequoia matrix", {
  E <- 0.02
  tab <- PaPeR:::build_locus_table(c(0.55, 0.45), Err = E)
  got <- implicit_emission_matrix(tab)
  expected <- matrix(c(
    1-E, E-(E/2)^2, (E/2)^2,
    E/2, 1-E, E/2,
    (E/2)^2, E-(E/2)^2, 1-E
  ), 3, 3, byrow = TRUE)
  expect_equal(got, expected, tolerance = 1e-15)
  expect_equal(dense_error_oracle(tab), expected, tolerance = 1e-15)
})

test_that("the implicit production kernel matches a dense independent oracle", {
  for (k in c(2L, 3L, 4L, 6L)) {
    E <- 0.013
    tab <- PaPeR:::build_locus_table(rep(1 / k, k), Err = E)
    G <- as.integer(k * (k + 1L) / 2L)
    got <- implicit_emission_matrix(tab)
    expect_equal(got, dense_error_oracle(tab), tolerance = 1e-14)
    expect_true(all(is.finite(got)))
    expect_true(all(got >= 0))
    expect_equal(diag(got), rep(1 - E, G), tolerance = 1e-14)
    expect_equal(rowSums(got), rep(1, G), tolerance = 1e-14)
    expect_equal(rowSums(got - diag(diag(got))), rep(E, G), tolerance = 1e-14)
  }

  k <- 30L
  tab <- PaPeR:::build_locus_table(rep(1 / k, k), Err = 0.013)
  set.seed(2401)
  cells <- cbind(sample.int(tab$G, 200L, replace = TRUE),
                 sample.int(tab$G, 200L, replace = TRUE))
  emission_cell <- function(true_index, observed_index) {
    true <- tab$geno[true_index, ]; observed <- tab$geno[observed_index, ]
    if (true_index == observed_index) return(1 - tab$Err)
    overlap <- sum(pmin(tabulate(true, nbins = tab$k),
                        tabulate(observed, nbins = tab$k)))
    cls <- if (overlap == 0L) "w_Far" else if (observed[1L] == observed[2L]) "w_Hom" else "w_Het"
    tab$Err * tab$ErrProfile[[cls]] / if (true[1L] == true[2L]) tab$Z_hom else tab$Z_het
  }
  got <- PaPeR:::emission_probabilities_cpp(tab, cells[, 1L] - 1L,
                                                 cells[, 2L] - 1L)
  expect_equal(got, vapply(seq_len(nrow(cells)), function(i)
    emission_cell(cells[i, 1L], cells[i, 2L]), numeric(1)), tolerance = 1e-14)
})

test_that("transition classes are equal and different homozygotes are far", {
  tab <- PaPeR:::build_locus_table(rep(1/4, 4), Err = 0.02,
                           ErrProfile = c(w_Het = .04, w_Hom = .02, w_Far = .01))
  E <- implicit_emission_matrix(tab)
  for (t in seq_len(tab$G)) {
    cls <- vapply(seq_len(tab$G), function(o) {
      if (o == t) return("diagonal")
      class_of_transition(tab$geno[t, ], tab$geno[o, ])
    }, character(1))
    for (z in c("w_Het", "w_Hom", "w_Far")) {
      vals <- E[t, cls == z]
      if (length(vals) > 1L) expect_equal(vals, rep(vals[1], length(vals)))
    }
  }
  i11 <- match("1:1", apply(tab$geno, 1, gkey))
  i22 <- match("2:2", apply(tab$geno, 1, gkey))
  i34 <- match("3:4", apply(tab$geno, 1, gkey))
  expect_equal(E[i11, i22], E[i11, i34])
})

test_that("allele-label permutation preserves corresponding probabilities", {
  k <- 4L
  tab <- PaPeR:::build_locus_table(rep(1/k, k), Err = 0.03,
                           ErrProfile = c(w_Het = .05, w_Hom = .03, w_Far = .01))
  perm <- c(3L, 1L, 4L, 2L)
  mapped <- apply(tab$geno, 1, function(g) gkey(perm[g]))
  idx <- match(mapped, apply(tab$geno, 1, gkey))
  E <- implicit_emission_matrix(tab)
  expect_equal(E, E[idx, idx], tolerance = 1e-15)
})

test_that("custom profiles and 3 x 3 adapters are equivalent", {
  profile <- c(w_Het = .04, w_Hom = .02, w_Far = .01)
  adapter <- matrix(c(
    1, .04, .01,
    .02, 1, .02,
    .01, .04, 1
  ), 3, 3, byrow = TRUE)
  direct <- PaPeR:::build_locus_table(rep(1/3, 3), Err = 0.02, ErrProfile = profile)
  adapted <- PaPeR:::build_locus_table(rep(1/3, 3), Err = 0.02, ErrProfile = adapter)
  expect_equal(implicit_emission_matrix(direct), implicit_emission_matrix(adapted))

  E <- 0.02
  stochastic <- matrix(c(
    1-E, E-(E/2)^2, (E/2)^2,
    E/2, 1-E, E/2,
    (E/2)^2, E-(E/2)^2, 1-E
  ), 3, 3, byrow = TRUE)
  expect_equal(PaPeR:::build_locus_table(c(0.5, 0.5), ErrProfile = stochastic)$Err, E)
  expect_error(PaPeR:::build_locus_table(c(0.5, 0.5), Err = 0.03,
                                 ErrProfile = stochastic), "conflicts")
})

test_that("missing calls are marginalized while homozygotes remain called", {
  tab <- PaPeR:::build_locus_table(c(0.5, 0.5), Err = 0.01)
  G <- prepare_genotypes(matrix(c("0", "0", "1/2", "1/1"), ncol = 1,
                                dimnames = list(c("missing_first", "missing_second",
                                                  "heterozygote", "homozygote"), "MH1")),
                         n_alleles = 2L)
  expect_equal(unname(PaPeR:::calc_pair_LL_cpp(list(tab), G, 0L, 2L)),
               rep(log(tab$mObs[gidx_error(1, 2, 2)]), nrow(compact_legacy_delta())))
  expect_equal(unname(PaPeR:::calc_pair_LL_cpp(list(tab), G, 1L, 3L)),
               rep(log(tab$mObs[gidx_error(1, 1, 2)]), nrow(compact_legacy_delta())))
  called <- PaPeR:::calc_pair_LL_cpp(list(tab), G, 2L, 3L)
  expect_true(any(called != 0))
  expect_gt(PaPeR:::emission_probabilities_cpp(
    tab, gidx_error(1, 2, 2) - 1L, gidx_error(1, 1, 2) - 1L), 0)
})

test_that("compact locus storage has no dense error matrix", {
  tab <- PaPeR:::build_locus_table(rep(1 / 30, 30), Err = 0.01)
  expect_false("E" %in% names(tab))
  expect_true(all(c("Err", "ErrProfile", "Z_hom", "Z_het", "k", "geno") %in%
                  names(tab)))
  compact <- tab[c("Err", "ErrProfile", "Z_hom", "Z_het", "k", "geno")]
  expect_lt(as.numeric(object.size(compact)),
            as.numeric(object.size(dense_error_oracle(tab))) / 10)
})

test_that("invalid and legacy specifications fail clearly", {
  bad <- list(NA_real_, Inf, -0.1, 1.1, c(0.1, 0.2, 0.3))
  for (x in bad)
    expect_error(build_tables(list(c(0.5, 0.5), c(0.5, 0.5)), Err = x))
  expect_error(build_tables(list(c(0.5, 0.5)),
                            Err = list(dropout = 0.01, other_err = 0.01)),
               "Err must")
  expect_error(build_tables(list(c(0.5, 0.5)), legacy = 0.01),
               "unused argument")
  expect_error(PaPeR:::build_locus_table(c(0.5, 0.5), Err = 0.01,
                                 ErrProfile = c(w_Het = 0, w_Hom = 0, w_Far = 0)),
               "non-normalizable")
  expect_error(PaPeR:::build_locus_table(c(0.5, 0.5), Err = 0.01,
                                 ErrProfile = c(w_Het = 1.1, w_Hom = .1, w_Far = .1)),
               "between 0 and 1")
  expect_error(PaPeR:::build_locus_table(c(0.5, 0.5), ErrProfile = matrix(0, 4, 4)),
               "3 x 3")
  asymmetric <- diag(3); asymmetric[1, 2] <- 0.01
  expect_error(PaPeR:::build_locus_table(c(0.5, 0.5), ErrProfile = asymmetric),
               "symmetric")
  expect_error(PaPeR:::build_locus_table(1, Err = 0.01), "monomorphic")
})

test_that("production likelihoods reject the retired error floor", {
  expect_error(
    build_tables(list(c(0.5, 0.5)), Err = 0.001, err_floor = 0.01),
    "unused argument"
  )
})
