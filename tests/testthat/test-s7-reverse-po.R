test_that("public trio likelihood matches an independent co-parent sum", {
  source(test_path("helper-parentage-fixture.R"), local = TRUE)
  calls <- rbind(D = c("1/1", "1/2"), F = c("1/2", "2/2"),
                 C = c("2/2", "1/2"))
  colnames(calls) <- paste0("L", seq_len(ncol(calls)))
  freq <- c(.37, .63)
  f <- .scientific_fixture(calls, dam = c(F = "D"), sire = c(F = "C"),
    p = freq, focal = "F", Complex = "full",
    birth_year = c(D = 0L, F = 2L, C = 0L))

  alleles <- rbind(c(1L, 1L), c(1L, 2L), c(2L, 2L))
  prior <- apply(alleles, 1L, function(g)
    prod(freq[g]) * if (g[1L] == g[2L]) 1 else 2)
  reference <- vapply(seq_len(ncol(calls)), function(locus) {
    emission <- exp(f$compact_inputs$observations_ln[[locus]])
    total <- 0
    for (child in 1:3) for (dam in 1:3) for (sire in 1:3) {
      transmitted <- 0
      for (a in 1:2) for (b in 1:2)
        transmitted <- transmitted + .25 *
          as.numeric(identical(sort(c(alleles[dam, a], alleles[sire, b])),
                               alleles[child, ]))
      total <- total + prior[dam] * prior[sire] * transmitted *
        emission[child, 2L] * emission[dam, 1L] * emission[sire, 3L]
    }
    log(total)
  }, numeric(1))
  got <- calc_trio_LL(calls, "F", "D", "C", control = f$public_control)
  expect_equal(got[["H_PP"]], sum(reference), tolerance = 1e-10)
})

test_that("public PO likelihood is direction-symmetric and ages orient the assignment", {
  source(test_path("helper-parentage-fixture.R"), local = TRUE)
  pattern <- rbind(
    c("1/1", "1/1", "1/2", "1/2", "1/2", "1/2", "2/2", "2/2"),
    c("1/1", "1/2", "1/1", "1/2", "1/2", "2/2", "1/2", "2/2"))
  calls <- pattern[, rep(1:8, 10)]
  rownames(calls) <- c("F", "C")
  colnames(calls) <- paste0("L", seq_len(ncol(calls)))
  f <- .scientific_fixture(calls, sire_candidates = "C", focal = "F",
    p = c(.5, .5), Complex = "full", birth_year = c(F = NA_integer_, C = NA_integer_))
  alleles <- rbind(c(1L, 1L), c(1L, 2L), c(2L, 2L))
  prior <- c(.25, .5, .25)
  one_parent <- matrix(0, 3L, 3L)
  for (parent in 1:3) for (child in 1:3) for (a in 1:2) for (b in 1:2)
    if (identical(sort(c(alleles[parent, a], b)), alleles[child, ]))
      one_parent[child, parent] <- one_parent[child, parent] +
        mean(alleles[parent, ] == a) * .5
  po_likelihood <- function(parent_row, child_row) {
    sum(vapply(seq_len(ncol(calls)), function(locus) {
      emission <- exp(f$compact_inputs$observations_ln[[locus]])
      total <- 0
      for (parent in 1:3) for (child in 1:3)
        total <- total + prior[parent] * one_parent[child, parent] *
          emission[parent, parent_row] * emission[child, child_row]
      log(total)
    }, numeric(1)))
  }
  expected <- po_likelihood(1L, 2L)
  expect_equal(expected, po_likelihood(2L, 1L), tolerance = 1e-12)
  pair <- calc_pair_LL(calls, "F", "C", control = f$public_control)
  expect_equal(pair[["PO"]], expected, tolerance = 1e-10)

  wide <- f
  wide$life_history$Sex <- c(1L, 2L)
  wide$life_history$BY.min <- c(0L, 0L)
  wide$life_history$BY.max <- c(10L, 10L)
  ambiguous <- .public_assignment(wide, control = list(
    hypotheses = "S7", Tfilter = -Inf, Tassign = .5, margin = .5))
  expect_true(all(is.na(ambiguous$dam)))
  expect_true(all(is.na(ambiguous$sire)))

  oriented <- f
  oriented$life_history$Sex <- c(1L, 2L)
  oriented$life_history$BirthYear <- c(10L, 0L)
  oriented$life_history$BY.min <- oriented$life_history$BY.max <- NA_integer_
  assigned <- .public_assignment(oriented, control = list(
    hypotheses = "S7", Tfilter = -Inf, Tassign = .5, margin = .5))
  expect_identical(assigned$sire[match("F", assigned$id)], "C")
  expect_true(is.na(assigned$sire[match("C", assigned$id)]))
})
