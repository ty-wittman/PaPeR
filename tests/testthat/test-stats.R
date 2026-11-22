# Exact and stratified multiallelic marker statistics.

pic_direct <- function(p) {
  pair <- 0
  if (length(p) > 1L) for (i in seq_len(length(p) - 1L))
    for (j in (i + 1L):length(p))
      pair <- pair + 2 * p[i]^2 * p[j]^2
  1 - sum(p^2) - pair
}

test_that("MHStats returns exact ALL frequencies and zero-only likelihood adjustment", {
  GenoM <- rbind(
    i1 = c("1/2", "1/2"),
    i2 = c("2/2", "1/2"),
    i3 = c("1/3", "1/2"),
    i4 = c("0",   "1/2")
  )
  colnames(GenoM) <- c("MH1", "MH2")
  out <- MHStats(GenoM, n_alleles = c(4L, 2L),
                 pseudocount = 0.5, warn_call_rate = 0, warn_Ae = 0)

  expect_s3_class(out, "PaPeR_MHStats")
  expect_identical(names(out), "ALL")
  expect_s3_class(out[[1]], "PaPeR_MHStats_slice")
  expect_identical(
    names(out[[1]]),
    c("group", "locus", "n_individuals", "n_typed", "call_rate",
      "n_alleles", "n_alleles_all", "n_alleles_universe", "Ae", "PIC", "allele_frequencies",
      "likelihood_frequencies", "n_zero_alleles", "zero_adjusted",
      "low_call_rate", "low_Ae")
  )

  p <- c(`1` = 2/6, `2` = 3/6, `3` = 1/6, `4` = 0)
  q <- c(`1` = 2/6.5, `2` = 3/6.5, `3` = 1/6.5, `4` = 0.5/6.5)
  expect_equal(out$ALL$allele_frequencies[[1]], p)
  expect_equal(out$ALL$likelihood_frequencies[[1]], q)
  expect_equal(out$ALL$allele_frequencies[[2]],
               out$ALL$likelihood_frequencies[[2]])
  expect_equal(out$ALL$call_rate, c(0.75, 1))
  expect_identical(out$ALL$n_alleles, c(3L, 2L))
  expect_identical(out$ALL$n_alleles_all, c(3L, 2L))
  expect_identical(out$ALL$n_alleles_universe, c(4L, 2L))
  expect_equal(out$ALL$Ae[1], 1 / sum(p^2))
  expect_equal(out$ALL$PIC[1], unname(pic_direct(p)))
  expect_true(out$ALL$zero_adjusted[1])
  expect_false(out$ALL$zero_adjusted[2])

  different_floor <- MHStats(GenoM, n_alleles = c(4L, 2L),
                             pseudocount = 10,
                             warn_call_rate = 0, warn_Ae = 0)
  expect_equal(different_floor$ALL$allele_frequencies,
               out$ALL$allele_frequencies)
  expect_equal(different_floor$ALL[c("call_rate", "n_alleles", "Ae", "PIC")],
               out$ALL[c("call_rate", "n_alleles", "Ae", "PIC")])
  expect_false(isTRUE(all.equal(
    different_floor$ALL$likelihood_frequencies[[1]],
    out$ALL$likelihood_frequencies[[1]]
  )))
})

test_that("Ae and PIC match hand calculations and thresholds are strict", {
  biallelic <- matrix(rep("1/2", 4), ncol = 1,
                      dimnames = list(paste0("i", 1:4), "B"))
  expect_warning(
    b <- MHStats(biallelic, warn_call_rate = 0.20, warn_Ae = 3),
    "Ae below 3"
  )
  expect_equal(b$ALL$Ae, 2)
  expect_equal(b$ALL$PIC, 0.375)

  triallelic <- matrix(c("1/2", "2/3", "1/3"), ncol = 1,
                       dimnames = list(paste0("t", 1:3), "T"))
  expect_no_warning(
    t <- MHStats(triallelic, warn_call_rate = 0.20, warn_Ae = 3)
  )
  expect_equal(t$ALL$Ae, 3)
  expect_equal(t$ALL$PIC, 16/27)
  expect_false(t$ALL$low_Ae)
})

test_that("ALL is first and cohort-by-sex slices follow genotype-row appearance", {
  GenoM <- matrix(c("1/1", "1/2", "2/2", "2/3", "3/3"), ncol = 1,
                  dimnames = list(paste0("i", 1:5), "MH"))
  by <- data.frame(
    ID = c("i3", "i1", "i5", "i2", "i4"),
    BirthYear = c(2021, 2021, 2020, 2020, 2021),
    Sex = c(1, 2, 3, 1, 1)
  )
  expect_warning(
    out <- MHStats(GenoM, n_alleles = 3L, by = by,
                   warn_call_rate = 0, warn_Ae = 0),
    "excluded from grouped slices"
  )
  expect_identical(
    names(out),
    c("ALL", "Cohort=2021;Sex=2", "Cohort=2020;Sex=1",
      "Cohort=2021;Sex=1")
  )
  expect_true(all(c("Cohort", "Sex") %in% names(out[[1]])))
  expect_true(all(is.na(out$ALL$Cohort)))
  expect_true(all(is.na(out$ALL$Sex)))
  expect_identical(out[[2]]$Cohort, "2021")
  expect_identical(out[[2]]$Sex, 2L)
})

test_that("Cohort and BirthYear are interchangeable aliases", {
  GenoM <- matrix(c("1/1", "1/2", "2/2", "2/3"), ncol = 1,
                  dimnames = list(paste0("i", 1:4), "MH"))
  cohort <- data.frame(ID = rownames(GenoM), Cohort = c("B", "A", "B", "A"))
  birth <- data.frame(ID = rownames(GenoM), BirthYear = c("B", "A", "B", "A"))
  a <- MHStats(GenoM, by = cohort, strata = "Cohort",
               warn_call_rate = 0, warn_Ae = 0)
  b <- MHStats(GenoM, by = birth, strata = "BirthYear",
               warn_call_rate = 0, warn_Ae = 0)
  expect_equal(a, b)
  expect_identical(names(a), c("ALL", "Cohort=B", "Cohort=A"))

  both <- data.frame(ID = rownames(GenoM), Cohort = c(1, 1, 2, 2),
                     BirthYear = c(1, 1, 2, 3))
  expect_error(MHStats(GenoM, by = both), "conflicting")
  both$BirthYear <- both$Cohort
  expect_identical(
    names(MHStats(GenoM, by = both, warn_call_rate = 0, warn_Ae = 0)),
    c("ALL", "Cohort=1", "Cohort=2")
  )
})

test_that("stratum allele counts expose dropout relative to ALL", {
  GenoM <- matrix(c("1/1", "1/2", "2/2", "3/3"), ncol = 1,
                  dimnames = list(paste0("i", 1:4), "MH"))
  by <- data.frame(ID = rownames(GenoM), Cohort = c(1L, 1L, 2L, 2L))
  out <- MHStats(GenoM, n_alleles = 4L, by = by,
                 warn_call_rate = 0, warn_Ae = 0)
  expect_identical(out$ALL$n_alleles, 3L)
  expect_identical(out$ALL$n_alleles_all, 3L)
  expect_identical(out[["Cohort=1"]]$n_alleles, 2L)
  expect_identical(out[["Cohort=1"]]$n_alleles_all, 3L)
  expect_length(out[["Cohort=1"]]$allele_frequencies[[1]], 4L)
  expect_equal(unname(out[["Cohort=1"]]$allele_frequencies[[1]]),
               c(0.75, 0.25, 0, 0))
})

test_that("call-rate warning uses individuals and respects the 20 percent boundary", {
  at_boundary <- matrix(c("1/2", rep("0", 4)), ncol = 1,
                        dimnames = list(paste0("i", 1:5), "MH"))
  expect_no_warning(
    a <- MHStats(at_boundary, n_alleles = 2L,
                 warn_call_rate = 0.20, warn_Ae = 0)
  )
  expect_equal(a$ALL$call_rate, 0.20)
  expect_false(a$ALL$low_call_rate)

  below <- matrix(c("1/2", rep("0", 5)), ncol = 1,
                  dimnames = list(paste0("j", 1:6), "MH"))
  expect_warning(
    b <- MHStats(below, n_alleles = 2L,
                 warn_call_rate = 0.20, warn_Ae = 0),
    "call rate below 0.2"
  )
  expect_equal(b$ALL$call_rate, 1/6)
  expect_true(b$ALL$low_call_rate)
})

test_that("untyped and zero-pseudocount slices are rejected by build_tables", {
  untyped <- matrix("0", nrow = 2, ncol = 1,
                    dimnames = list(c("i1", "i2"), "MH"))
  expect_warning(
    u <- MHStats(untyped, n_alleles = 2L, warn_Ae = 0),
    "call rate below 0.2"
  )
  expect_true(all(is.na(u$ALL$allele_frequencies[[1]])))
  expect_true(is.na(u$ALL$Ae))
  expect_true(is.na(u$ALL$PIC))
  expect_error(build_tables(u$ALL), "unavailable or zero")

  with_zero <- matrix(c("1/1", "1/2"), ncol = 1,
                      dimnames = list(c("a", "b"), "MH"))
  z <- MHStats(with_zero, n_alleles = 3L, pseudocount = 0,
               warn_call_rate = 0, warn_Ae = 0)
  expect_equal(z$ALL$allele_frequencies[[1]],
               z$ALL$likelihood_frequencies[[1]])
  expect_false(z$ALL$zero_adjusted)
  expect_error(build_tables(z$ALL), "unavailable or zero")
})

test_that("build_tables consumes an MHStats slice and retains legacy lists", {
  GenoM <- matrix(c("1/1", "1/2", "2/3"), ncol = 1,
                  dimnames = list(c("a", "b", "c"), "MH"))
  out <- MHStats(GenoM, n_alleles = 4L,
                 warn_call_rate = 0, warn_Ae = 0)
  tab <- build_tables(out$ALL, Err = 0.01)
  expect_equal(tab[[1]]$p, out$ALL$likelihood_frequencies[[1]])
  expect_error(build_tables(out), "select one MHStats group slice")
  expect_length(build_tables(list(MH = c(0.5, 0.5)), Err = 1e-12), 1L)
})

test_that("selected panels retain locus identities through likelihood construction", {
  GenoM <- cbind(
    L1 = c(a = "1/1", b = "1/2", c = "2/3"),
    L2 = c(a = "1/2", b = "2/2", c = "1/3")
  )
  stats <- MHStats(GenoM, n_alleles = c(L1 = 3L, L2 = 3L),
                   warn_call_rate = 0, warn_Ae = 0)
  tables <- build_tables(stats$ALL, Err = 0.01)
  expect_identical(names(tables), c("L1", "L2"))
  expect_equal(tables[[1]]$p, stats$ALL$likelihood_frequencies[[1]])
  expect_error(.calc_pair_LL(GenoM[, c("L2", "L1")], "a", "b", tables),
               "locus names and order")
  expect_error(.calc_pair_LL(GenoM, "a", "b", tables[c("L2", "L1")]),
               "locus names and order")

  duplicate_locus <- stats$ALL
  duplicate_locus$locus[2L] <- "L1"
  expect_error(build_tables(duplicate_locus), "unique, non-missing locus names")
  wrong_length <- stats$ALL
  wrong_length$likelihood_frequencies[[1L]] <- c(`1` = 0.5, `2` = 0.5)
  expect_error(build_tables(wrong_length), "allele universe")
})

test_that("frequency summaries are invariant to rows and transform with allele labels", {
  GenoM <- cbind(
    MH = c(a = "1/1", b = "1/2", c = "2/3", d = "3/3")
  )
  baseline <- MHStats(GenoM, n_alleles = 3L,
                      warn_call_rate = 0, warn_Ae = 0)$ALL
  reordered <- MHStats(GenoM[c("d", "b", "a", "c"), , drop = FALSE],
                       n_alleles = 3L, warn_call_rate = 0, warn_Ae = 0)$ALL
  expect_equal(reordered$allele_frequencies, baseline$allele_frequencies)
  expect_equal(reordered$likelihood_frequencies, baseline$likelihood_frequencies)

  relabel <- c(`1` = 3L, `2` = 1L, `3` = 2L)
  relabeled <- apply(GenoM, 2L, function(x) {
    vapply(strsplit(x, "/", fixed = TRUE), function(g) {
      paste(relabel[g], collapse = "/")
    }, character(1))
  })
  rownames(relabeled) <- rownames(GenoM)
  relabeled_stats <- MHStats(relabeled, n_alleles = 3L,
                              warn_call_rate = 0, warn_Ae = 0)$ALL
  expect_equal(unname(relabeled_stats$allele_frequencies[[1L]][relabel]),
               unname(baseline$allele_frequencies[[1L]]))
  expect_equal(unname(relabeled_stats$likelihood_frequencies[[1L]][relabel]),
               unname(baseline$likelihood_frequencies[[1L]]))
})

test_that("likelihood panels reject invalid frequencies before table construction", {
  bad <- list(
    negative = c(0.5, -0.1, 0.6),
    nonfinite = c(0.5, Inf, 0.5),
    zero_total = c(0, 0),
    retained_zero = c(0.5, 0.5, 0)
  )
  for (p in bad)
    expect_error(build_tables(list(MH = p)), "likelihood frequencies")
})

test_that("metadata and controls are validated", {
  GenoM <- matrix(c("1/1", "1/2"), ncol = 1,
                  dimnames = list(c("a", "b"), "MH"))
  expect_error(MHStats(GenoM, by = data.frame(ID = "a", Cohort = 1)),
               "exactly match")
  expect_error(MHStats(GenoM, by = data.frame(ID = c("a", "a"), Cohort = 1:2)),
               "unique")
  expect_error(MHStats(GenoM, by = data.frame(ID = c("a", "b"), Sex = c(1, 4))),
               "Sex must")
  expect_error(MHStats(GenoM, by = data.frame(ID = c("a", "b"), Cohort = 1:2),
                       strata = "Sex"),
               "Sex stratum")
  expect_error(MHStats(GenoM, pseudocount = -1), "pseudocount")
  expect_error(MHStats(GenoM, warn_call_rate = 1.1), "must not exceed")
  expect_error(MHStats(GenoM, warn_Ae = NA_real_), "warn_Ae")

  by_matrix <- cbind(ID = c("b", "a"), BirthYear = c("2", "1"))
  out <- MHStats(GenoM, by = by_matrix, warn_call_rate = 0, warn_Ae = 0)
  expect_identical(names(out), c("ALL", "Cohort=1", "Cohort=2"))
})
