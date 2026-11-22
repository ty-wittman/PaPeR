final_llr_toy <- function() {
  g <- rbind(D = c("1/1", "1/2", "2/2"), S = c("2/2", "1/2", "1/2"),
             F = c("1/2", "1/1", "1/2"), X = c("1/1", "1/2", "1/2"))
  colnames(g) <- paste0("L", 1:3)
  list(g = g, lh = data.frame(ID = rownames(g), Sex = c(1L, 2L, 1L, 2L),
                              BirthYear = c(0L, 0L, 5L, 5L)),
       ped = data.frame(id = rownames(g), dam = c(NA, NA, "D", "D"),
                        sire = c(NA, NA, "S", NA)),
       ctl = list(frequencies = rep(list(c(.6, .4)), 3), error_calibration = "none"))
}

test_that("fixed final LLRs run without search or candidate discovery and preserve the pedigree", {
  x <- final_llr_toy(); original <- serialize(x, NULL)
  local_mocked_bindings(.pap_search = function(...) stop("search called"),
    .build_candidate_pool = function(...) stop("candidate discovery called"))
  for (mode in c("ancestral", "full")) for (suite in c("S7", "S37", "S45")) {
    ctl <- c(x$ctl, list(llr_evidence = mode, hypotheses = suite))
    out <- calc_llr(x$ped, x$g, x$lh, .01, .001, control = ctl)
    expect_identical(out[, c("id", "dam", "sire")], x$ped)
    expect_identical(serialize(x, NULL), original)
    expect_true(all(is.finite(out$LLRdam[3:4])))
    expect_true(is.finite(out$LLRsire[3]))
    expect_equal(out$LLRpair[3], min(out$LLRdam[3], out$LLRsire[3]))
    expect_true(is.na(out$LLRsire[4])); expect_true(is.na(out$LLRpair[4]))
    expect_equal(out$OHdam[3:4], c(0L, 0L)); expect_equal(out$MEpair[3], 0L)
    expect_identical(out$DamLLRStatus, c("not_assigned", "not_assigned", "available", "available"))
    expect_identical(out$SireLLRStatus[4], "not_assigned")
    expect_identical(attr(out, "final_llr_evidence"), mode)
    expect_identical(attr(out, "final_llr_hypotheses"), suite)
    shuffled <- calc_llr(x$ped[c(4, 2, 3, 1), ], x$g, x$lh, .01, .001, control = ctl)
    expect_equal(shuffled$LLRdam[match(out$id, shuffled$id)], out$LLRdam, tolerance = 1e-9)
  }
})


test_that("missing calls, weak supplied assignments and sire-only links remain scorable", {
  x <- final_llr_toy()
  x$g[,] <- "0"
  z <- calc_llr(x$ped, x$g, x$lh, .01, .001, control = x$ctl)
  expect_equal(z$LLRdam[3:4], c(0, 0), tolerance = 1e-10)
  expect_equal(z$LLRsire[3], 0, tolerance = 1e-10)
  x <- final_llr_toy(); x$ped$dam[3] <- NA_character_
  x$g["F", ] <- "1/1"; x$g["S", ] <- "2/2"
  z <- calc_llr(x$ped, x$g, x$lh, .01, .001,
    control = c(x$ctl, list(Tfilter = 100, MaxMendelE = 0L, Tassign = 1000, margin = 1000)))
  expect_true(is.finite(z$LLRsire[3])); expect_lt(z$LLRsire[3], 0)
  expect_equal(z$OHsire[3], 3L); expect_true(is.na(z$LLRdam[3])); expect_true(is.na(z$LLRpair[3]))
  expect_identical(z$SireLLRStatus[3], "available")
  ordinary <- calc_llr(x$ped, x$g, x$lh, .01, .001, control = x$ctl)
  expect_equal(z$LLRsire, ordinary$LLRsire, tolerance = 1e-9)
})

test_that("fixed scoring also accepts the Sequoia SNP error model", {
  x <- final_llr_toy()
  out <- calc_llr(x$ped, x$g, x$lh,
    control = c(x$ctl, list(error_model = "sequoia", Err = .01)))
  expect_true(all(is.finite(out$LLRdam[3:4])))
  expect_true(is.finite(out$LLRsire[3]))
  expect_identical(out[, c("id", "dam", "sire")], x$ped)
  expect_identical(attr(out, "log_base"), 10L)
})

test_that("unavailable side scores propagate to the pair without discarding the other side", {
  x <- final_llr_toy(); native <- run_parentage_cpp
  local_mocked_bindings(run_parentage_cpp = function(compiled, compact, input) {
    expect_equal(input$likelihood_memory_mib, 128)
    expect_equal(input$likelihood_work_millions, 800)
    z <- native(compiled, compact, input)
    d <- as.data.frame(z$configuration_diagnostics)
    d$full_dam_relationship_ln[d$focal == 3L] <- NA_real_
    z$configuration_diagnostics <- d
    z
  })
  out <- calc_llr(x$ped, x$g, x$lh, .01, .001,
    control = c(x$ctl, list(memory_mib = 128, work_millions = 800)))
  expect_true(is.na(out$LLRdam[3])); expect_true(is.na(out$LLRpair[3]))
  expect_true(is.finite(out$LLRsire[3]))
  expect_identical(out$DamLLRStatus[3], "unavailable")
  expect_identical(out$PairLLRStatus[3], "unavailable")
  expect_identical(out$SireLLRStatus[3], "available")
})

test_that("CalcLLR preserves the search and uses its final calibrated rates", {
  g <- rbind(D = rep("1/1", 20), S = rep("2/2", 20), F = rep("1/2", 20), X = rep("1/2", 20))
  colnames(g) <- paste0("L", 1:20)
  lh <- data.frame(ID = rownames(g), Sex = c(1, 2, 1, 2), BirthYear = c(0, 0, 4, 4))
  for (calibration in c("none", "raise_guarded")) {
    ctl <- list(frequencies = rep(list(c(.5, .5)), 20), error_calibration = calibration)
    before <- assign_parents(g, lh, .005, .0001, control = ctl)
    after <- assign_parents(g, lh, .005, .0001, control = c(ctl, list(CalcLLR = TRUE)))
    expect_identical(after[, c("id", "dam", "sire")], before[, c("id", "dam", "sire")])
    for (field in c("calibration", "rounds", "converged", "unavailable"))
      expect_identical(attr(after, field), attr(before, field))
    # Wall-clock measurements vary between independent searches.
    science <- function(x) {
      if (!is.list(x)) return(x)
      x <- x[!grepl("seconds$", names(x))]
      lapply(x, science)
    }
    expect_identical(lapply(attr(after, "calibration_passes"), science),
                     lapply(attr(before, "calibration_passes"), science))
    cal <- tail(attr(after, "calibration"), 1)
    scored <- calc_llr(before, g, lh, cal$Dropout_used, cal$Foreign_used, control = ctl)
    expect_equal(after$LLRdam, scored$LLRdam, tolerance = 1e-9)
    expect_equal(after$LLRsire, scored$LLRsire, tolerance = 1e-9)
    expect_equal(after$LLRpair, pmin(after$LLRdam, after$LLRsire), tolerance = 1e-12)
  }
})

test_that("fixed scoring validates life history and pedigree instead of changing links", {
  x <- final_llr_toy()
  expect_error(calc_llr(x$ped, x$g, x$lh[, c("ID", "Sex")], .01, .001, control = x$ctl), "BirthYear")
  bad <- x$ped; bad$dam[1] <- "F"
  expect_error(calc_llr(bad, x$g, x$lh, .01, .001, control = x$ctl), "structural validation")
  bad <- x$ped; bad$dam[3] <- "untyped"
  expect_error(calc_llr(bad, x$g, x$lh, .01, .001, control = x$ctl), "absent|unknown|typed")
  empty <- calc_llr(x$ped[FALSE, ], x$g, x$lh, .01, .001, control = x$ctl)
  expect_equal(nrow(empty), 0L)
})

test_that("CalcLLR consumes the final tables returned by calibration", {
  g <- rbind(D = rep("1/1", 20), S = rep("2/2", 20), F = rep("1/2", 20), X = rep("1/2", 20))
  colnames(g) <- paste0("L", 1:20)
  lh <- data.frame(ID = rownames(g), Sex = c(1, 2, 1, 2), BirthYear = c(0, 0, 4, 4))
  ctl <- list(frequencies = rep(list(c(.5, .5)), 20))
  calibrate <- .pap_calibrate
  local_mocked_bindings(.pap_calibrate = function(...) {
    z <- calibrate(...)
    # Isolate the lifecycle boundary from the estimator: the final tables
    # deliberately differ from the supplied ones.
    z$final_tables <- .pap_tables(g, .05, .0001, pap_control(ctl))
    z
  })
  out <- assign_parents(g, lh, .005, .0001, control = c(ctl, list(CalcLLR = TRUE)))
  ped <- out[, c("id", "dam", "sire")]
  final <- calc_llr(ped, g, lh, .05, .0001, control = ctl)
  supplied <- calc_llr(ped, g, lh, .005, .0001, control = ctl)
  expect_equal(out$LLRdam, final$LLRdam, tolerance = 1e-9)
  expect_equal(out$LLRsire, final$LLRsire, tolerance = 1e-9)
  expect_gt(max(abs(out$LLRdam - supplied$LLRdam), na.rm = TRUE), .001)
})

test_that("final-LLR evidence is its own setting, independent of the search evidence", {
  x <- final_llr_toy()
  a <- calc_llr(x$ped, x$g, x$lh, .01, .001, control = x$ctl)
  expect_identical(attr(a, "final_llr_evidence"), "ancestral")
  f <- calc_llr(x$ped, x$g, x$lh, .01, .001, control = c(x$ctl, list(llr_evidence = "full")))
  s <- calc_llr(x$ped, x$g, x$lh, .01, .001, control = c(x$ctl, list(evidence = "full")))
  expect_identical(attr(s, "final_llr_evidence"), "ancestral")
  expect_equal(s$LLRdam, a$LLRdam, tolerance = 1e-12)
  expect_identical(attr(f, "final_llr_evidence"), "full")
  expect_error(pap_control(list(llr_evidence = "descendant")), "ancestral, full")
})
