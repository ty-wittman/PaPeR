test_that("public pre-screen settings preserve the exhaustive final score", {
  source(test_path("helper-parentage-fixture.R"), local = TRUE)
  calls <- rbind(
    A = rep(c("1/2", "1/2", "1/1"), 8),
    B = rep(c("1/1", "1/2", "1/1"), 8),
    C = rep(c("2/2", "1/2", "1/2"), 8),
    D = rep(c("2/2", "2/2", "2/2"), 8))
  colnames(calls) <- paste0("L", seq_len(ncol(calls)))
  f <- .scientific_fixture(calls, dam_candidates = c("B", "D"),
    sire_candidates = "C", focal = "A", Complex = "full",
    birth_year = c(A = 3L, B = 1L, C = 1L, D = 1L))
  viable_dams <- sum(f$life_history$Sex == 1L &
    f$life_history$BirthYear < f$life_history$BirthYear[f$life_history$ID == "A"],
    na.rm = TRUE)
  expect_identical(viable_dams, 2L)
  common <- list(evidence = "full", hypotheses = "S7", CalcLLR = TRUE)
  exhaustive <- .public_assignment(f,
    control = c(common, list(prescreen_K = Inf)))
  permissive <- .public_assignment(f,
    control = c(common, list(prescreen_K = 1L, prescreen_margin_ln = 1e6)))
  at <- match("A", permissive$id)
  expect_lt(attr(permissive, "settings")$prescreen_K, viable_dams)
  expect_identical(permissive$dam[at], "B")
  expect_identical(permissive$sire[at], "C")
  expect_true(all(is.finite(unlist(permissive[at,
    c("LLRdam", "LLRsire", "LLRpair")], use.names = FALSE))))
  expect_identical(exhaustive[, c("id", "dam", "sire")],
                   permissive[, c("id", "dam", "sire")])
  expect_equal(exhaustive[, c("LLRdam", "LLRsire", "LLRpair")],
               permissive[, c("LLRdam", "LLRsire", "LLRpair")],
               tolerance = 1e-8)
  expect_identical(attr(exhaustive, "settings")$prescreen_K, Inf)
  expect_identical(attr(permissive, "settings")$prescreen_K, 1L)
  expect_identical(attr(permissive, "settings")$prescreen_margin_ln, 1e6)
})
