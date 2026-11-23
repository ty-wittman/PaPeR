test_that("named hypothesis suites are supported public opt-ins", {
  expect_identical(pap_control()$hypotheses, "S7")
  for (suite in c("S7", "S37", "S45"))
    expect_identical(pap_control(list(hypotheses = suite))$hypotheses, suite)
  expect_error(pap_control(list(hypotheses = "S209")), "S7, S37, S45")
})

test_that("public search accepts each named suite without changing a clear pedigree", {
  calls <- rbind(D = rep("1/1", 30), S = rep("2/2", 30),
                 A = rep("1/2", 30), B = rep("1/2", 30))
  colnames(calls) <- paste0("L", seq_len(ncol(calls)))
  life <- data.frame(ID = rownames(calls), Sex = c(1L, 2L, 1L, 2L),
                     BirthYear = c(0L, 0L, 2L, 2L))
  frequencies <- setNames(rep(list(c(.5, .5)), ncol(calls)), colnames(calls))
  expected <- data.frame(id = rownames(calls),
    dam = c(NA_character_, NA_character_, "D", "D"),
    sire = c(NA_character_, NA_character_, "S", "S"),
    stringsAsFactors = FALSE)

  for (suite in c("S7", "S37", "S45")) {
    out <- assign_parents(calls, life, control = list(
      error_model = "sequoia", Err = .001, frequencies = frequencies,
      error_calibration = "none", hypotheses = suite))
    expect_identical(out[, c("id", "dam", "sire")], expected)
    expect_identical(attr(out, "settings")$hypotheses, suite)
  }
})
