# Typed-parent OH exclusion filter: direct-comparison correctness and calibration.

test_that("calibration retains at least the target and grows with error", {
  p_lo <- rep(0.005, 157)
  p_hi <- rep(0.03, 157)
  c_lo <- calibrate_maxmismatch(p_lo, retain = 0.999)
  c_hi <- calibrate_maxmismatch(p_hi, retain = 0.999)
  expect_gte(c_lo$retained, 0.999)
  expect_gte(c_hi$retained, 0.999)
  expect_gt(c_hi$T, c_lo$T)
  expect_equal(sum(PaPeR:::.pois_binom_pmf(p_lo)), 1)
})

test_that("po_mismatch_prob retains the expected probability ordering", {
  set.seed(5)
  freqs <- lapply(c(2, 4, 6, 8), function(k) {
    g <- rgamma(k, 0.8, 1)
    (g + 1e-3) / sum(g + 1e-3)
  })
  pm <- po_mismatch_prob(build_tables(freqs, Err = 0.01))
  expect_true(all(pm$p_um >= pm$p_fm))
  expect_true(all(pm$p_fm >= 0 & pm$p_um <= 1))
})

test_that("safe-range probabilities and threshold equal an independent mask reference", {
  set.seed(55)
  freqs <- lapply(c(3L, 6L, 8L), function(k) {
    x <- rgamma(k, 0.8, 1)
    x / sum(x)
  })
  tables <- build_tables(freqs, Err = 0.012)
  old_reference <- function(tab) {
    masks <- apply(tab$geno, 1L, function(g)
      bitwOr(bitwShiftL(1L, g[1L] - 1L), bitwShiftL(1L, g[2L] - 1L)))
    disjoint <- outer(masks, masks, Vectorize(function(a, b) bitwAnd(a, b) == 0L))
    Ypo <- dense_observed_relationship_oracle(tab, compact_legacy_delta()["PO", ])
    Yu <- dense_observed_relationship_oracle(tab, compact_legacy_delta()["U", ])
    c(p_fm = sum(Ypo[disjoint]) / sum(Ypo),
      p_um = sum(Yu[disjoint]) / sum(Yu))
  }
  old <- t(vapply(tables, old_reference, numeric(2)))
  new <- po_mismatch_prob(tables)
  expect_equal(new$p_fm, old[, "p_fm"], tolerance = 1e-12)
  expect_equal(new$p_um, old[, "p_um"], tolerance = 1e-12)
  expect_equal(calibrate_maxmismatch(rep(new$p_fm, 25L), 0.999)$T,
               calibrate_maxmismatch(rep(old[, "p_fm"], 25L), 0.999)$T)
})

test_that("compact OH scalars cover error profiles without observed-pair work", {
  cases <- list(
    list(p = c(.8, .2), Err = 1e-12),
    list(p = c(.91, .06, .03), Err = .017),
    list(p = c(.6, .15, .1, .06, .04, .025, .015, .01), Err = .025,
         ErrProfile = c(w_Het = .02, w_Hom = .01, w_Far = .003))
  )
  for (case in cases) {
    tab <- do.call(PaPeR:::build_locus_table, case)
    expect_true(all(c("p_oh_true_po", "p_oh_unrelated") %in% names(tab)))
    got <- po_mismatch_prob(list(L = tab))
    expect_equal(unname(got$p_fm), unname(tab$p_oh_true_po), tolerance = 0)
    expect_equal(unname(got$p_um), unname(tab$p_oh_unrelated), tolerance = 0)
    expect_true(all(unlist(got[c("p_fm", "p_um")]) >= 0 &
                    unlist(got[c("p_fm", "p_um")]) <= 1))
  }
  error_free <- PaPeR:::build_locus_table(c(.7, .2, .1), Err = 1e-12)
  expect_gt(unname(po_mismatch_prob(list(error_free))$p_fm), 0)
  expect_lt(unname(po_mismatch_prob(list(error_free))$p_fm), 1e-12)
  k30 <- PaPeR:::build_locus_table(rep(1 / 30, 30), Err = .012)
  expect_true(is.finite(k30$p_oh_true_po))
  expect_true(is.finite(k30$p_oh_unrelated))
})

