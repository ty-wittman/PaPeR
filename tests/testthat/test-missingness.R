# Independent exhaustive reference calculations for the Stage-1 §2.5
# missing-observation contract.  These helpers deliberately do not call the
# production pair or trio likelihood entry points.

missing_gidx <- function(a, b, k) {
  if (a > b) { z <- a; a <- b; b <- z }
  (a - 1L) * k - (a - 1L) * (a - 2L) / 2L + (b - a) + 1L
}

missing_dense_emission <- function(tab) {
  out <- matrix(0, tab$G, tab$G)
  for (truth in seq_len(tab$G)) {
    gt <- tab$geno[truth, ]
    for (observed in seq_len(tab$G)) {
      go <- tab$geno[observed, ]
      if (truth == observed) {
        out[truth, observed] <- 1 - tab$Err
      } else if (tab$Err > 0) {
        overlap <- sum(pmin(tabulate(gt, nbins = tab$k),
                            tabulate(go, nbins = tab$k)))
        cls <- if (overlap == 0L) "w_Far" else if (go[1L] == go[2L]) {
          "w_Hom"
        } else {
          "w_Het"
        }
        z <- if (gt[1L] == gt[2L]) tab$Z_hom else tab$Z_het
        out[truth, observed] <- tab$Err * tab$ErrProfile[[cls]] / z
      }
    }
  }
  out
}

missing_gamete_probability <- function(allele, parent) {
  sum(parent == allele) / 2
}

missing_two_parent_probability <- function(child, dam, sire) {
  child <- sort(as.integer(child))
  a <- child[1L]; b <- child[2L]
  if (a == b) {
    missing_gamete_probability(a, dam) * missing_gamete_probability(a, sire)
  } else {
    missing_gamete_probability(a, dam) * missing_gamete_probability(b, sire) +
      missing_gamete_probability(b, dam) * missing_gamete_probability(a, sire)
  }
}

missing_one_parent_probability <- function(child, parent, p) {
  child <- sort(as.integer(child))
  a <- child[1L]; b <- child[2L]
  if (a == b) {
    missing_gamete_probability(a, parent) * p[a]
  } else {
    missing_gamete_probability(a, parent) * p[b] +
      missing_gamete_probability(b, parent) * p[a]
  }
}

missing_log_probability <- function(x) if (x == 0) -Inf else log(x)

missing_emission_vector <- function(E, index) {
  if (is.na(index)) rep(1, nrow(E)) else E[, index]
}

missing_trio_oracle <- function(tab, o_index, d_index, s_index) {
  G <- tab$G
  E <- missing_dense_emission(tab)
  e_o <- missing_emission_vector(E, o_index)
  e_d <- missing_emission_vector(E, d_index)
  e_s <- missing_emission_vector(E, s_index)
  h_pp <- h_d0 <- h_0s <- 0

  for (g_d in seq_len(G)) {
    for (g_s in seq_len(G)) {
      for (g_o in seq_len(G)) {
        transmission <- missing_two_parent_probability(
          tab$geno[g_o, ], tab$geno[g_d, ], tab$geno[g_s, ])
        h_pp <- h_pp + tab$pi[g_d] * e_d[g_d] * tab$pi[g_s] * e_s[g_s] *
          e_o[g_o] * transmission
        h_d0 <- h_d0 + tab$pi[g_d] * e_d[g_d] * tab$pi[g_s] *
          e_o[g_o] * transmission
        h_0s <- h_0s + tab$pi[g_d] * tab$pi[g_s] * e_s[g_s] *
          e_o[g_o] * transmission
      }
    }
  }
  h_d0 <- h_d0 * sum(tab$pi * e_s)
  h_0s <- h_0s * sum(tab$pi * e_d)
  h_u <- sum(tab$pi * e_o) * sum(tab$pi * e_d) * sum(tab$pi * e_s)
  vapply(c(H_PP = h_pp, H_D0 = h_d0, H_0S = h_0s, H_U = h_u),
         missing_log_probability, numeric(1))
}

missing_one_parent_oracle <- function(tab, o_index, parent_index) {
  E <- missing_dense_emission(tab)
  e_o <- E[, o_index]
  e_parent <- E[, parent_index]
  value <- 0
  for (g_parent in seq_len(tab$G)) {
    for (g_o in seq_len(tab$G)) {
      value <- value + tab$pi[g_parent] * e_parent[g_parent] * e_o[g_o] *
        missing_one_parent_probability(tab$geno[g_o, ], tab$geno[g_parent, ], tab$p)
    }
  }
  missing_log_probability(value)
}

missing_pair_oracle <- function(tab, first_index, second_index, delta) {
  E <- missing_dense_emission(tab)
  e_first <- missing_emission_vector(E, first_index)
  e_second <- missing_emission_vector(E, second_index)
  psi <- dense_relationship_oracle(tab, delta)
  missing_log_probability(sum(psi * outer(e_first, e_second)))
}

missing_public_genotypes <- function(tab, indices, ids = c("o", "d", "s")) {
  calls <- vapply(indices, function(index) {
    if (is.na(index)) "0" else paste(tab$geno[index, ], collapse = "/")
  }, character(1))
  matrix(calls, ncol = 1L, dimnames = list(ids, "MH1"))
}

test_that("pair likelihoods marginalize every missing-observation pattern", {
  tab <- PaPeR:::build_locus_table(c(0.55, 0.3, 0.15), Err = 0.02)
  observed_a <- missing_gidx(1L, 3L, tab$k)
  observed_b <- missing_gidx(2L, 2L, tab$k)
  patterns <- list(
    called = c(observed_a, observed_b),
    first_missing = c(NA_integer_, observed_b),
    second_missing = c(observed_a, NA_integer_),
    both_missing = c(NA_integer_, NA_integer_)
  )

  for (pattern in patterns) {
    G <- prepare_genotypes(
      missing_public_genotypes(tab, pattern, ids = c("first", "second")),
      n_alleles = tab$k)
    got <- PaPeR:::calc_pair_LL_cpp(list(tab), G, 0L, 1L)
    delta <- compact_legacy_delta()
    expected <- vapply(seq_len(nrow(delta)), function(r)
      missing_pair_oracle(tab, pattern[1L], pattern[2L], delta[r, ]), numeric(1))
    names(expected) <- rownames(delta)
    expect_equal(got, expected, tolerance = 1e-12)
  }

  # The active ordinary registry is outbred.  The one-sided exact calculations
  # therefore agree across its rows, but that equality is established here
  # rather than used as a hard-coded missing-locus shortcut.
  one_sided <- patterns$first_missing
  G <- prepare_genotypes(
    missing_public_genotypes(tab, one_sided, ids = c("first", "second")),
    n_alleles = tab$k)
  got <- PaPeR:::calc_pair_LL_cpp(list(tab), G, 0L, 1L)
  expect_equal(unname(got), rep(log(tab$mObs[observed_b]), length(got)),
               tolerance = 1e-12)

  ctx <- PaPeR:::.new_inference_context(G, list(MH1 = tab), Complex = "simp")
  score <- PaPeR:::.score_unconditional_pairs(ctx, "first", "second")
  active<-match(colnames(score),ctx$registry$config_id)
  expect_false(anyNA(active))
  expected_config <- vapply(ctx$registry$Delta[active], function(delta)
    missing_pair_oracle(tab, one_sided[1L], one_sided[2L], delta), numeric(1))
  expect_equal(unname(score[1L, ]), unname(expected_config), tolerance = 1e-12)

  # The public ordinary screen may cancel a one-sided term, but an asymmetric
  # Jacquard state need not. D3 has an inbred first individual and D5 is its
  # focal reversal; both directions must be marginalized rather than skipped.
  custom_delta <- rbind(
    D3 = c(0, 0, 1, rep(0, 6L)),
    D5 = c(rep(0, 4L), 1, rep(0, 4L)),
    U = c(rep(0, 8L), 1)
  )
  asymmetric_observed <- missing_gidx(1L, 1L, tab$k)
  first_observed <- c(asymmetric_observed, NA_integer_)
  second_observed <- rev(first_observed)
  screen_missing <- function(pattern, focal) {
    G_missing <- prepare_genotypes(
      missing_public_genotypes(tab, pattern, ids = c("first", "second")),
      n_alleles = tab$k)
    PaPeR:::screen_rel_pairs_cpp(
      list(tab), G_missing, custom_delta, focal, 2L, -Inf,
      rep(NA_integer_, 2L), -1L, 0L, -1L)
  }
  oracle_contrast <- function(pattern, focal) {
    missing_pair_oracle(tab, pattern[1L], pattern[2L],
                         custom_delta[focal + 1L, ]) -
      missing_pair_oracle(tab, pattern[1L], pattern[2L], custom_delta[3L, ])
  }

  d3_first <- screen_missing(first_observed, focal = 0L)
  d3_second <- screen_missing(second_observed, focal = 0L)
  expect_equal(d3_first$llr,
               oracle_contrast(first_observed, focal = 0L), tolerance = 1e-12)
  expect_equal(d3_second$llr,
               oracle_contrast(second_observed, focal = 0L), tolerance = 1e-12)
  expect_gt(abs(d3_first$llr - d3_second$llr), 1e-6)

  # Reversing focal individuals requires D3 -> D5 and swaps the input rows. It
  # preserves the appropriate directional marginal, not an accidental
  # row-order-invariant value.
  d5_reversed_first <- screen_missing(second_observed, focal = 1L)
  d5_reversed_second <- screen_missing(first_observed, focal = 1L)
  expect_equal(d5_reversed_first$llr,
               oracle_contrast(second_observed, focal = 1L),
               tolerance = 1e-12)
  expect_equal(d5_reversed_second$llr,
               oracle_contrast(first_observed, focal = 1L),
               tolerance = 1e-12)
  expect_equal(d3_first$llr, d5_reversed_first$llr, tolerance = 1e-12)
  expect_equal(d3_second$llr, d5_reversed_second$llr, tolerance = 1e-12)
})

test_that("trio likelihoods exactly marginalize all eight observation patterns", {
  patterns <- expand.grid(o = c(TRUE, FALSE), d = c(TRUE, FALSE),
                          s = c(TRUE, FALSE), KEEP.OUT.ATTRS = FALSE)
  for (k in c(2L, 3L, 6L)) {
    p <- seq_len(k)^-0.7
    p <- p / sum(p)
    o_index <- missing_gidx(1L, k, k)
    d_index <- missing_gidx(1L, 1L, k)
    s_index <- missing_gidx(k, k, k)
    for (Err in c(1e-12, 0.02)) {
      tab <- PaPeR:::build_locus_table(p, Err = Err)
      for (row in seq_len(nrow(patterns))) {
        observed <- c(
          if (patterns$o[row]) o_index else NA_integer_,
          if (patterns$d[row]) d_index else NA_integer_,
          if (patterns$s[row]) s_index else NA_integer_
        )
        G <- prepare_genotypes(missing_public_genotypes(tab, observed), n_alleles = k)
        got <- PaPeR:::calc_trio_LL_cpp(list(tab), G, 0L, 1L, 2L)
        expected <- missing_trio_oracle(tab, observed[1L], observed[2L], observed[3L])
        expect_equal(got, expected, tolerance = 1e-11)
        expect_false(anyNA(got))
        expect_false(any(is.nan(got)))

        if (!patterns$o[row])
          expect_equal(unname(got), rep(unname(got["H_PP"]), length(got)),
                       tolerance = 1e-12)
        if (patterns$o[row] && patterns$d[row] && !patterns$s[row]) {
          expect_equal(unname(got["H_PP"]), unname(got["H_D0"]), tolerance = 1e-12)
          expect_equal(unname(got["H_0S"]), unname(got["H_U"]), tolerance = 1e-12)
        }
        if (patterns$o[row] && !patterns$d[row] && patterns$s[row]) {
          expect_equal(unname(got["H_PP"]), unname(got["H_0S"]), tolerance = 1e-12)
          expect_equal(unname(got["H_D0"]), unname(got["H_U"]), tolerance = 1e-12)
        }
        if (patterns$o[row] && !patterns$d[row] && !patterns$s[row])
          expect_equal(unname(got), rep(unname(got["H_PP"]), length(got)),
                       tolerance = 1e-12)
      }
    }
  }
})

test_that("dam and sire missing branches are symmetric and retain structural zeros", {
  tab <- PaPeR:::build_locus_table(c(0.5, 0.3, 0.2), Err = 0.02)
  child <- missing_gidx(1L, 3L, tab$k)
  parent <- missing_gidx(1L, 2L, tab$k)
  dams_missing <- prepare_genotypes(
    missing_public_genotypes(tab, c(child, NA_integer_, parent)), n_alleles = tab$k)
  sires_missing <- prepare_genotypes(
    missing_public_genotypes(tab, c(child, parent, NA_integer_)), n_alleles = tab$k)
  left <- PaPeR:::calc_trio_LL_cpp(list(tab), dams_missing, 0L, 1L, 2L)
  right <- PaPeR:::calc_trio_LL_cpp(list(tab), sires_missing, 0L, 1L, 2L)
  expect_equal(unname(left["H_PP"]), unname(right["H_PP"]), tolerance = 1e-12)
  expect_equal(unname(left["H_0S"]), unname(right["H_D0"]), tolerance = 1e-12)
  expect_equal(unname(left["H_D0"]), unname(right["H_0S"]), tolerance = 1e-12)
  expect_equal(unname(left["H_U"]), unname(right["H_U"]), tolerance = 1e-12)

  zero_tab <- PaPeR:::build_locus_table(c(0.5, 0.3, 0.2), Err = 1e-12)
  impossible <- prepare_genotypes(
    missing_public_genotypes(zero_tab, c(missing_gidx(3L, 3L, 3L),
                                         missing_gidx(1L, 1L, 3L),
                                         missing_gidx(2L, 2L, 3L))),
    n_alleles = zero_tab$k)
  got_zero <- PaPeR:::calc_trio_LL_cpp(list(zero_tab), impossible, 0L, 1L, 2L)
  expect_true(all(is.finite(got_zero[c("H_PP", "H_D0", "H_0S")])))
  expect_true(is.finite(got_zero["H_U"]))
  expect_false(any(is.nan(got_zero)))
})

test_that("the optimized one-parent branch matches an independent k=30 oracle", {
  k <- 30L
  p <- seq_len(k)^-0.8
  p <- p / sum(p)
  tab <- PaPeR:::build_locus_table(p, Err = 0.01)
  observed <- c(missing_gidx(1L, k, k), missing_gidx(1L, 1L, k), NA_integer_)
  G <- prepare_genotypes(missing_public_genotypes(tab, observed), n_alleles = k)
  got <- PaPeR:::calc_trio_LL_cpp(list(tab), G, 0L, 1L, 2L)
  one_parent <- missing_one_parent_oracle(tab, observed[1L], observed[2L])
  unrelated <- log(tab$mObs[observed[1L]] * tab$mObs[observed[2L]])
  expected <- c(H_PP = one_parent, H_D0 = one_parent,
                H_0S = unrelated, H_U = unrelated)
  expect_equal(got, expected, tolerance = 1e-11)
  expect_false(anyNA(got))
  expect_false(any(is.nan(got)))
})
