test_that("public fixed scoring excludes descendants only in ancestral mode", {
  source(test_path("helper-parentage-fixture.R"), local = TRUE)
  calls <- rbind(
    D = c("1/1", "1/2", "2/2", "1/2"),
    S = c("1/2", "1/1", "1/2", "2/2"),
    F = c("1/2", "1/1", "1/2", "1/2"),
    B = c("1/1", "1/2", "1/2", "2/2"),
    C = c("2/2", "2/2", "2/2", "2/2"),
    X = c("2/2", "2/2", "2/2", "2/2"))
  colnames(calls) <- paste0("L", 1:4)
  f <- .scientific_fixture(calls,
    dam = c(F = "D", B = "D"), sire = c(F = "S", B = "S", C = "F", X = "S"),
    dam_candidates = "D", sire_candidates = "S", focal = "F", Complex = "full",
    birth_year = c(D = 0L, S = 0L, F = 2L, B = 2L, C = 4L, X = 4L))
  original <- serialize(list(f$calls, f$pedigree, f$life_history), NULL)
  ancestral <- .public_fixed_score(f,
    control = list(llr_evidence = "ancestral", hypotheses = "S7"))
  full <- .public_fixed_score(f,
    control = list(llr_evidence = "full", hypotheses = "S7"))
  focal <- match("F", ancestral$id)
  expect_identical(serialize(list(f$calls, f$pedigree, f$life_history), NULL), original)
  expect_identical(ancestral$DamLLRStatus[focal], "available")
  expect_identical(ancestral$SireLLRStatus[focal], "available")
  expect_true(is.finite(ancestral$LLRdam[focal]))
  expect_true(is.finite(ancestral$LLRsire[focal]))

  # Alter only descendant C. The public ancestral contrast is fixed to the
  # focal's assigned-parent cavity, while full evidence includes C's calls.
  descendant <- f
  descendant$calls["C", ] <- c("1/1", "1/1", "1/1", "1/1")
  ancestral_after <- .public_fixed_score(descendant,
    control = list(llr_evidence = "ancestral", hypotheses = "S7"))
  full_after <- .public_fixed_score(descendant,
    control = list(llr_evidence = "full", hypotheses = "S7"))
  expect_equal(ancestral_after$LLRdam[focal], ancestral$LLRdam[focal], tolerance = 1e-9)
  expect_gt(abs(full_after$LLRdam[focal] - full$LLRdam[focal]), 1e-6)
})
