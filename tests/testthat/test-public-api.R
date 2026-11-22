# The public interface: assign_parents(GenoM, LifeHistData, Dropout, Foreign,
# control) and pap_control(). A two-parent, two-offspring toy pedigree with 30
# fully informative loci keeps every call fast.

toy <- function() {
  g <- rbind(D = rep("1/1", 30), S = rep("2/2", 30), A = rep("1/2", 30), B = rep("1/2", 30))
  colnames(g) <- paste0("L", 1:30)
  list(g = g, lh = data.frame(ID = rownames(g), Sex = c(1L, 2L, 1L, 2L),
                              BirthYear = c(0L, 0L, 2L, 2L)),
       freqs = rep(list(c(.5, .5)), 30))
}

test_that("pap_control returns every setting and rejects unknown or invalid ones", {
  s <- pap_control()
  expect_identical(s$error_model, "dropout_foreign")
  expect_identical(s$error_calibration, "raise_guarded")
  expect_identical(s$evidence, "ancestral")
  expect_identical(s$hypotheses, "S7")
  expect_equal(s$Tassign, 0.5); expect_equal(s$Tfilter, -2); expect_equal(s$prescreen_K, 10)
  expect_equal(s$prescreen_margin_ln, log(100)); expect_equal(s$max_rounds, 41)
  expect_false(s$CalcLLR)
  expect_identical(attr(pap_control(list(Tassign = 1)), "supplied"), "Tassign")
  expect_error(pap_control(list(Tasign = 1)), "unknown control setting.*Tasign")
  expect_error(pap_control(list(1)), "named list")
  expect_error(pap_control(list(error_calibration = "two_sided")), "raise_guarded, raise_only, none")
  expect_error(pap_control(list(evidence = "descendant")), "ancestral, full")
  expect_error(pap_control(list(hypotheses = "S99")), "S7, S37, S45")
  expect_error(pap_control(list(memory_mib = 100)), "64, 128, 256 or 512")
  expect_error(pap_control(list(work_millions = 1)), "400, 800")
  expect_error(pap_control(list(prescreen_K = 0)), "prescreen_K")
  expect_error(pap_control(list(Tfilter = Inf)), "Tfilter")
  expect_error(pap_control(list(error_model = "sequoia")), "requires control\\$Err")
  expect_error(pap_control(list(Err = .01)), "only to error_model")
  expect_false(s$verbose)
  expect_error(pap_control(list(verbose = NA)), "verbose must be TRUE or FALSE")
  expect_equal(pap_control(list(prescreen_K = Inf))$prescreen_K, Inf)
})

test_that("a resolved settings list keeps the names the user set", {
  s <- pap_control(list(prescreen_K = 5))
  again <- pap_control(s)
  expect_identical(attr(again, "supplied"), "prescreen_K")
  controls <- .engine_controls(again, 1)
  expect_identical(controls$pair_rescue_top_k, 5L)
  expect_null(controls$pair_rescue_margin_ln)
  default <- .engine_controls(pap_control(), 1)
  expect_identical(default$pedigree_hypothesis_suite, "S7")
  expect_identical(default$likelihood_evidence, "ancestral")
  expect_null(default$pair_rescue)
  expect_null(default$likelihood_memory_mib)
  off <- .engine_controls(pap_control(list(prescreen_K = Inf, memory_mib = 128)), 1)
  expect_false(off$pair_rescue); expect_identical(off$pair_rescue_top_k, 100000L)
  expect_identical(off$likelihood_memory_mib, 128)
  expect_true(off$likelihood_high_allele_programs)
})

test_that("both error models assign the toy parents and report one calibration row per search", {
  x <- toy()
  mh <- assign_parents(x$g, x$lh, 0.005, 1e-4, control = list(frequencies = x$freqs))
  snp <- assign_parents(x$g, x$lh, control = list(frequencies = x$freqs,
                                                  error_model = "sequoia", Err = 0.001))
  for (ped in list(mh, snp)) {
    expect_identical(names(ped), c("id", "dam", "sire"))
    expect_identical(ped$dam, c(NA, NA, "D", "D")); expect_identical(ped$sire, c(NA, NA, "S", "S"))
    expect_true(attr(ped, "converged"))
    cal <- attr(ped, "calibration")
    expect_s3_class(cal, "data.frame")
    expect_identical(nrow(cal), length(attr(ped, "calibration_passes")))
    expect_true(all(c("search", "confident_trios", "decision", "converged", "rounds") %in% names(cal)))
    expect_identical(attr(ped, "settings")$error_calibration, "raise_guarded")
  }
  expect_true(all(c("Dropout_used", "Dropout_estimate", "Dropout_lower", "Dropout_upper",
                    "Foreign_used") %in% names(attr(mh, "calibration"))))
  expect_true("Err_used" %in% names(attr(snp, "calibration")))
  expect_equal(attr(mh, "calibration")$Foreign_used, rep(1e-4, nrow(attr(mh, "calibration"))))
})

test_that("frequencies are estimated from GenoM unless supplied", {
  x <- toy()
  estimated <- assign_parents(x$g, x$lh, 0.005, 1e-4, control = list(error_calibration = "none"))
  supplied <- assign_parents(x$g, x$lh, 0.005, 1e-4,
    control = list(error_calibration = "none", frequencies = x$freqs))
  expect_identical(estimated[, 1:3], supplied[, 1:3])
  expect_error(assign_parents(x$g, x$lh), "Dropout and Foreign are required")
  expect_error(assign_parents(x$g, x$lh, 0.005, 1e-4,
    control = list(error_model = "sequoia", Err = .01)), "apply only")
  expect_error(assign_parents(x$g, x$lh, 0, 1e-4), "strictly positive")
})

test_that("calibration policies, CalcLLR and explicit screens reach the result", {
  x <- toy()
  none <- assign_parents(x$g, x$lh, 0.005, 1e-4, control = list(error_calibration = "none"))
  expect_identical(attr(none, "calibration")$decision, "not_calibrated")
  raise <- assign_parents(x$g, x$lh, 0.005, 1e-4, control = list(error_calibration = "raise_only"))
  expect_identical(raise[, 1:3], none[, 1:3])
  llr <- assign_parents(x$g, x$lh, 0.005, 1e-4,
    control = list(error_calibration = "none", CalcLLR = TRUE))
  expect_true(all(c("LLRdam", "LLRsire", "LLRpair", "OHdam", "OHsire") %in% names(llr)))
  expect_true(all(is.finite(llr$LLRdam[3:4])))
  expect_true(all(llr$LLRpair[3:4] > 0))
  fixed <- assign_parents(x$g, x$lh, 0.005, 1e-4,
    control = list(error_calibration = "none", T = 3L, MaxMendelE = 4L))
  expect_identical(attr(fixed, "MaxMismatch"), 3L)
  expect_identical(as.integer(attr(fixed, "MaxMendelE")), 4L)
})

test_that("a supplied starting pedigree is validated and used", {
  x <- toy()
  start <- data.frame(id = c("A", "B"), dam = c("D", NA), sire = c(NA, "S"))
  ped <- assign_parents(x$g, x$lh, 0.005, 1e-4,
    control = list(error_calibration = "none", pedigree = start))
  expect_identical(ped$dam[3:4], c("D", "D")); expect_identical(ped$sire[3:4], c("S", "S"))
  bad <- data.frame(id = "D", dam = "A", sire = NA)
  expect_error(assign_parents(x$g, x$lh, 0.005, 1e-4,
    control = list(error_calibration = "none", pedigree = bad)), "pedigree")
})

test_that("every calibration search restarts from the original inputs", {
  # Mock the search, the confident-trio audit and the fit: the driver alone is
  # under test. Pass 1 lowers Dropout once (supplied > 2 x upper, 70% confident).
  env <- new.env(parent = environment(.pap_calibrate)); calls <- list()
  pedigree <- data.frame(id = paste0("i", 1:10), dam = rep("d", 10), sire = rep("s", 10))
  env$.pap_search <- function(GenoM, LifeHistData, tables, settings, generation, progress = NULL) {
    calls[[length(calls) + 1L]] <<- list(tables = tables, settings = settings, generation = generation)
    list(pedigree = pedigree, converged = TRUE, rounds = 1L, MaxMismatch = 5L, MaxMendelE = 5L,
         graph_revision = 0L, focal_summaries = list(), complete_trio_record_state = NULL,
         unavailable = NULL, native_seconds = 0, checkpoint = list(compact = list(), input = list()))
  }
  env$.parentage_final_confidence <- function(...) list(confident = c(rep(TRUE, 7), rep(FALSE, 3)))
  env$estimate_parentage_error <- function(GenoM, tables, ...) list(status = "fitted",
    estimated_rates = c(Dropout = .005, Foreign = .001),
    interval = rbind(Dropout = c(lower = .001, upper = .01), Foreign = c(lower = .0005, upper = .002)))
  run <- .pap_calibrate; environment(run) <- env
  tables <- build_tables(list(c(.5, .5)), ErrorModel = "dropout_foreign", Dropout = .1, Foreign = .05)
  settings <- pap_control(list(pedigree = pedigree[1, ]))
  z <- run(NULL, NULL, tables, settings)
  expect_length(calls, 2L)
  for (a in calls) expect_identical(a$settings$pedigree, pedigree[1, ])
  expect_identical(calls[[1]]$generation, calls[[2]]$generation)
  expect_equal(calls[[1]]$tables[[1]]$Dropout, .1)
  expect_equal(calls[[2]]$tables[[1]]$Dropout, .01)
  expect_equal(calls[[2]]$tables[[1]]$Foreign, .05)
  expect_equal(z$passes[[1]]$confident_fraction, .7)
  expect_identical(z$passes[[2]]$stop_reason, "supplied_rate_retained")
  table <- .pap_calibration_table(z$passes)
  expect_equal(table$Dropout_used, c(.1, .01)); expect_equal(table$Foreign_used, c(.05, .05))
  expect_equal(table$Dropout_upper, c(.01, .01))
})

test_that("the pair utilities take the same error arguments as assign_parents", {
  x <- toy()
  g <- rbind(x$g, A2 = x$g["A", ])
  tables <- build_tables(x$freqs, ErrorModel = "dropout_foreign", Dropout = .005, Foreign = 1e-4)
  ctl <- list(frequencies = x$freqs)
  expect_equal(find_duplicates(g, .005, 1e-4, thresh = 1, control = ctl),
               .find_duplicates(g, tables, thresh = 1))
  expect_identical(find_duplicates(g, .005, 1e-4, thresh = 1, control = ctl)[1, c("ID1", "ID2")],
                   data.frame(ID1 = "A", ID2 = "A2"))
  expect_equal(calc_pair_LL(x$g, "A", "D", .005, 1e-4, control = ctl),
               .calc_pair_LL(x$g, "A", "D", tables))
  expect_equal(calc_trio_LL(x$g, "A", "D", "S", .005, 1e-4, control = ctl),
               .calc_trio_LL(x$g, "A", "D", "S", tables))
  ped <- data.frame(id = c("A", "B"), dam = "D", sire = "S")
  expect_equal(get_maybe_rel(x$g, .005, 1e-4, Pedigree = ped, thresh = -Inf, control = ctl),
               .get_maybe_rel(x$g, tables, Pedigree = ped, thresh = -Inf))
  snp <- calc_trio_LL(x$g, "A", "D", "S", control = c(ctl, list(error_model = "sequoia", Err = .001)))
  expect_true(snp[["H_PP"]] > snp[["H_U"]])
  expect_error(calc_pair_LL(x$g, "A", "D"), "Dropout and Foreign are required")
  expect_error(find_duplicates(x$g, 0, 1e-4), "strictly positive")
})

test_that("score_pedigree counts links and projects untyped true parents to none", {
  truth <- data.frame(id = c("D", "S", "A", "B"), dam = c(NA, NA, "D", "X"), sire = c(NA, NA, "S", "S"))
  ped <- data.frame(id = c("D", "S", "A", "B"), dam = c(NA, NA, "D", "D"), sire = c(NA, NA, "S", NA))
  raw <- score_pedigree(ped, truth)
  expect_identical(c(raw$dam$correct, raw$dam$wrong, raw$sire$correct, raw$sire$missed), c(1L, 1L, 1L, 1L))
  typed <- score_pedigree(ped, truth, typed_ids = c("D", "S", "A", "B"))
  expect_identical(typed$dam$false_positive_true_null, 1L)
  expect_identical(typed$dam$recall, 1)
})

test_that("verbose mode reports each pass and search without changing the result", {
  x <- toy()
  quiet <- assign_parents(x$g, x$lh, 0.005, 1e-4, control = list(frequencies = x$freqs))
  run <- function(dynamic) {
    withr::local_options(PaPeR.progress_dynamic = dynamic)
    shown <- testthat::capture_messages(loud <- assign_parents(x$g, x$lh, 0.005, 1e-4,
      control = list(frequencies = x$freqs, verbose = TRUE)))
    expect_identical(loud[, 1:3], quiet[, 1:3])
    paste(shown, collapse = "")
  }
  for (dynamic in c(FALSE, TRUE)) {
    text <- run(dynamic)
    expect_match(text, "PaPeR: 4 individuals, 30 loci; Dropout 0.005, Foreign 1e-04")
    expect_match(text, "search 1: [0-9]+ rounds, 2 dams, 2 sires")
    expect_match(text, "done: 2 dams, 2 sires")
  }
  # Log mode writes one line per round; terminal mode shows the phase status.
  expect_match(run(FALSE), "search 1  round 1: 2 moves")
  expect_match(run(TRUE), "search 1  preparing candidates")
  # On a terminal the status line is redrawn in place.
  p <- .pap_progress(TRUE); p$dynamic <- TRUE; p$search <- 1L
  bar <- testthat::capture_messages(.pap_visit_callback(p)(c(2L, 0L, 5L, 10L, 3L)))
  expect_match(paste(bar, collapse = ""), "search 1  round 2 +\\[#####-----\\]  50%")
  # Redraws wait two seconds.
  again <- testthat::capture_messages(.pap_visit_callback(p)(c(2L, 0L, 6L, 10L, 3L)))
  expect_length(again, 0L)
})
