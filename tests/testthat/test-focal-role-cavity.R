test_that("public fixed scoring treats an unknown-sex parent symmetrically across roles", {
  calls <- rbind(A = rep("1/2", 6), X = rep("1/1", 6))
  colnames(calls) <- paste0("L", 1:6)
  life <- data.frame(ID = rownames(calls), Sex = c(NA_integer_, NA_integer_),
                     BirthYear = c(3L, 0L))
  control <- list(error_model = "sequoia", Err = .01,
                  frequencies = rep(list(c(.4, .6)), 6), error_calibration = "none")
  pedigree <- data.frame(id = rownames(calls), dam = c("X", NA_character_),
                         sire = NA_character_)
  original <- serialize(list(calls, life, pedigree), NULL)
  dam <- calc_llr(pedigree, calls, life, control = control)
  expect_identical(serialize(list(calls, life, pedigree), NULL), original)
  pedigree$dam[1] <- NA_character_
  pedigree$sire[1] <- "X"
  sire <- calc_llr(pedigree, calls, life, control = control)
  expect_identical(dam$DamLLRStatus[1], "available")
  expect_identical(sire$SireLLRStatus[1], "available")
  expect_true(is.finite(dam$LLRdam[1]))
  expect_equal(dam$LLRdam[1], sire$LLRsire[1], tolerance = 1e-9)
  life$Sex[2] <- 1L
  expect_error(calc_llr(pedigree, calls, life, control = control),
               "known_sex_role_conflict")
})
