test_that("public ancestral pre-screen preserves a clear assignment under a restrictive top-K", {
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

  exhaustive <- .public_assignment(f,
    control = list(evidence = "ancestral", hypotheses = "S7",
      prescreen_K = Inf, CalcLLR = TRUE))
  screened <- .public_assignment(f,
    control = list(evidence = "ancestral", hypotheses = "S7",
      prescreen_K = 1L, prescreen_margin_ln = 1e6, CalcLLR = TRUE))
  expect_lt(attr(screened, "settings")$prescreen_K, viable_dams)
  at <- match("A", screened$id)
  expect_identical(screened$dam[at], "B")
  expect_identical(screened$sire[at], "C")
  expect_true(all(is.finite(unlist(screened[at,
    c("LLRdam", "LLRsire", "LLRpair")], use.names = FALSE))))
  expect_identical(exhaustive[, c("id", "dam", "sire")],
                   screened[, c("id", "dam", "sire")])
  expect_equal(exhaustive[, c("LLRdam", "LLRsire", "LLRpair")],
               screened[, c("LLRdam", "LLRsire", "LLRpair")],
               tolerance = 1e-8)
})
