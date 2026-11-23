test_that("public calls prepare canonically and preserve metadata", {
  x <- rbind(i1 = c("1/1", "1/3"), i2 = c("3/1", "0"),
             i3 = c("0/3", "3/0"), i4 = c("0/0", "2/2"))
  colnames(x) <- c("MH_A", "MH_B")
  G <- prepare_genotypes(x, c(MH_A = 3L, MH_B = 4L))
  expect_s3_class(G, "PaPeR_prepared_genotypes")
  expect_identical(prepare_genotypes(G), G)
  expect_equal(unname(G["i1", 3:4]), c(1L, 3L))
  expect_equal(unname(G["i2", 1:2]), c(1L, 3L))
  expect_equal(unname(G["i3", ]), rep(0L, 4L))
  expect_equal(unname(G["i4", 1:2]), c(0L, 0L))
  expect_identical(attr(G, "locus_names"), c("MH_A", "MH_B"))
  expect_identical(attr(G, "n_alleles"), c(3L, 4L))
})

test_that("n_alleles is inferred from the largest observed code", {
  x <- rbind(a = c("1/3", "2/2"), b = c("3/1", "1/4"))
  colnames(x) <- c("A", "B")
  G <- prepare_genotypes(x)
  expect_identical(attr(G, "n_alleles"), c(3L, 4L))
  s <- MHStats(x, warn_call_rate = 0, warn_Ae = 0)$ALL
  expect_identical(s$locus, c("A", "B"))
  expect_equal(lengths(s$allele_frequencies), c(3L, 4L))
})

test_that("the prepared representation records the complete input contract", {
  x <- rbind(alpha = c("33/7", "0/2", "0"),
             beta  = c("7/33", "2/2", "1/1"))
  colnames(x) <- c("L1", "L2", "L3")
  G <- prepare_genotypes(x, c(L1 = 40L, L2 = 2L, L3 = 2L))
  expect_identical(attr(G, "row_ids"), rownames(x))
  expect_identical(attr(G, "locus_columns"),
                   matrix(c(1L, 2L, 3L, 4L, 5L, 6L), ncol = 2L,
                          byrow = TRUE,
                          dimnames = list(colnames(x), c("copy1", "copy2"))))
  expect_equal(unname(G["alpha", 1:2]), c(7L, 33L))
  expect_equal(unname(G["alpha", 3:6]), c(0L, 0L, 0L, 0L))
  expect_identical(prepare_genotypes(G), G)

  broken <- G
  attr(broken, "row_ids") <- rev(attr(G, "row_ids"))
  expect_error(prepare_genotypes(broken), "invalid prepared")
  broken <- G
  attr(broken, "locus_columns")[1L, ] <- c(2L, 1L)
  expect_error(prepare_genotypes(broken), "copy-column mapping")
  expect_error(prepare_genotypes(G, c(40L, 2L, 3L)), "do not match")
})

test_that("partial and completely missing loci follow the explicit missing contract", {
  x <- rbind(i1 = c("0/3", "0/0", "0"),
             i2 = c("1/2", "1/1", "0"))
  colnames(x) <- c("partial", "zero_pair", "untyped")
  G <- prepare_genotypes(x, c(3L, 3L, 5L))
  expect_equal(unname(G["i1", ]), rep(0L, 6L))
  expect_identical(attr(G, "n_alleles"), c(3L, 3L, 5L))
  expect_error(prepare_genotypes(x[, "untyped", drop = FALSE]), "all-missing")
  expect_error(prepare_genotypes(
    matrix("0/4", 1L, 1L, dimnames = list("i", "L")), 3L),
    "exceeds n_alleles")
})

test_that("explicit row and locus permutations preserve the canonical calls", {
  x <- rbind(a = c("1/4", "2/3"), b = c("4/1", "3/2"), c = c("0", "1/1"))
  colnames(x) <- c("L1", "L2")
  G <- prepare_genotypes(x, c(L1 = 4L, L2 = 3L))
  xp <- x[c("c", "a", "b"), c("L2", "L1"), drop = FALSE]
  Gp <- prepare_genotypes(xp, c(L2 = 3L, L1 = 4L))
  expect_equal(as.integer(Gp[c("a", "b", "c"), c(3L, 4L, 1L, 2L)]),
               as.integer(G))
})

test_that("character data frames are accepted when their IDs are explicit", {
  x <- data.frame(L1 = c("1/2", "0"), L2 = c("2/2", "1/3"),
                  row.names = c("i1", "i2"), check.names = FALSE,
                  stringsAsFactors = FALSE)
  G <- prepare_genotypes(x, c(L1 = 2L, L2 = 3L))
  expect_equal(as.integer(G["i1", ]), c(1L, 2L, 2L, 2L))
  expect_identical(attr(G, "row_ids"), c("i1", "i2"))
})

test_that("invalid public genotype inputs fail at the boundary", {
  good <- matrix("1/1", 1, 1, dimnames = list("i1", "L1"))
  no_rows <- good; rownames(no_rows) <- NULL
  no_cols <- good; colnames(no_cols) <- NULL
  dup_rows <- rbind(i = "1/1", i = "1/1"); colnames(dup_rows) <- "L"
  dup_cols <- cbind(L = c(i = "1/1"), L = c(i = "1/1"))
  auto_rows <- data.frame(L = "1/1")
  expect_error(prepare_genotypes(no_rows), "row names")
  expect_error(prepare_genotypes(no_cols), "column names")
  expect_error(prepare_genotypes(dup_rows), "row names")
  expect_error(prepare_genotypes(dup_cols), "column names")
  expect_error(prepare_genotypes(auto_rows), "explicit individual row names")
  invalid <- c(NA, "", "1", "1/", "/1", "1/2/3", "-1/2", "1/-2",
               "1.0/2", "1e2/2", "a/2", " 1/2", "1 /2", "1/ 2",
               "01/2", "+1/2")
  for (z in invalid) expect_error(prepare_genotypes(
    matrix(z, 1, 1, dimnames = list("i", "L")), 3L))
  expect_error(prepare_genotypes(matrix("1/4", 1, 1,
                                        dimnames = list("i", "L")), 3L),
               "exceeds n_alleles")
  expect_error(prepare_genotypes(matrix("0/4", 1, 1,
                                        dimnames = list("i", "L")), 3L),
               "exceeds n_alleles")
  expect_error(prepare_genotypes(matrix("0", 1, 1,
                                        dimnames = list("i", "L"))), "all-missing")
})

