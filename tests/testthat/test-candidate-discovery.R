# Section 4.1 candidate-pool and calibrated trio-ME gates.

.section41_oracle <- function(tab) {
  k <- tab$k
  G <- tab$G
  p <- tab$p
  genotype_index <- function(a, b) {
    aa <- min(a, b); bb <- max(a, b)
    which(tab$geno[, 1L] == aa & tab$geno[, 2L] == bb)[1L] - 1L
  }
  emission <- function(true_index, observed_index) {
    truth <- tab$geno[true_index + 1L, ]
    observed <- tab$geno[observed_index + 1L, ]
    if (all(truth == observed)) return(1 - tab$Err)
    if (tab$Err == 0) return(0)
    overlap <- if (truth[1L] == truth[2L]) sum(observed == truth[1L]) else
      (observed[1L] == truth[1L] || observed[2L] == truth[1L]) +
      (observed[1L] == truth[2L] || observed[2L] == truth[2L])
    weight <- if (overlap == 0) tab$ErrProfile[["w_Far"]] else
      if (observed[1L] == observed[2L]) tab$ErrProfile[["w_Hom"]] else
        tab$ErrProfile[["w_Het"]]
    tab$Err * weight / if (truth[1L] == truth[2L]) tab$Z_hom else tab$Z_het
  }
  r <- matrix(0, k, G)
  for (a in seq_len(k)) for (o in 0:(G - 1L))
    r[a, o + 1L] <- sum(vapply(seq_len(k), function(c) {
      emission(genotype_index(a, c), o)
    }, numeric(1)) * p)
  compatible <- function(dam, sire, offspring) {
    ox <- tab$geno[offspring + 1L, 1L]
    oy <- tab$geno[offspring + 1L, 2L]
    dam <- tab$geno[dam + 1L, ]; sire <- tab$geno[sire + 1L, ]
    ((ox %in% dam) && (oy %in% sire)) ||
      ((oy %in% dam) && (ox %in% sire))
  }
  p_compatible <- 0
  for (a in seq_len(k)) for (b in seq_len(k))
    for (o in 0:(G - 1L)) {
      e <- emission(genotype_index(a, b), o)
      if (e == 0) next
      for (dam in 0:(G - 1L)) for (sire in 0:(G - 1L))
        if (compatible(dam, sire, o))
          p_compatible <- p_compatible + p[a] * p[b] * e *
            r[a, dam + 1L] * r[b, sire + 1L]
    }
  1 - p_compatible
}

test_that("compact true-trio ME probability matches an independent oracle", {
  for (k in c(2L, 3L, 6L)) {
    p <- seq_len(k)^(-1.3); p <- p / sum(p)
    tab <- PaPeR:::build_locus_table(p, Err = 0.017)
    expect_equal(tab$p_me_true_trio, .section41_oracle(tab), tolerance = 1e-10)
  }
  profile <- c(w_Het = 0.02, w_Hom = 0.01, w_Far = 0.003)
  tab <- PaPeR:::build_locus_table(
    c(.91, .06, .02, .01), Err = 0.02, ErrProfile = profile)
  expect_equal(tab$p_me_true_trio, .section41_oracle(tab), tolerance = 1e-10)
  expect_equal(PaPeR:::build_locus_table(c(.8, .2), Err = 1e-12)$p_me_true_trio, 0)
  k30 <- PaPeR:::build_locus_table(rep(1 / 30, 30), Err = .012)
  expect_true(is.finite(k30$p_me_true_trio))
  expect_true(k30$p_me_true_trio >= 0 && k30$p_me_true_trio <= 1)
})

test_that("ME calibration is exact Poisson-binomial and missingness is conservative", {
  probs <- c(.01, .2, .4, .03)
  pmf <- PaPeR:::.pois_binom_pmf(probs)
  expected <- which(cumsum(pmf) >= .9)[1L] - 1L
  got <- calibrate_max_mendelian_error(probs, retain = .9)
  expect_equal(got$T, expected)
  expect_gte(got$achieved_retention, .9)
  expect_equal(got$expected_false_me_count, sum(probs))
  expect_equal(got$p_me_true_trio, probs)
})

test_that("co-parent ME batches remain tuple-specific", {
  x <- rbind(o = c("1/1", "1/1"), a = c("1/2", "1/2"),
             b = c("2/2", "2/2"), c = c("1/1", "1/1"))
  colnames(x) <- c("L1", "L2")
  G <- prepare_genotypes(x, c(2L, 2L))
  bad <- PaPeR:::count_mendelian_errors_batch_cpp(G, 0L, 1L, 2L)
  good <- PaPeR:::count_mendelian_errors_batch_cpp(G, 0L, 1L, 3L)
  expect_equal(unname(bad), 2L)
  expect_equal(unname(good), 0L)
})

