# Independent arithmetic and legacy-caller witnesses for a scored 0/0 ratio.
# Native common-background Check acceptance is tested separately: these R
# descriptor checks do not certify the older frozen-message evidence model.

.zero_contract_rows <- function(objective, dam = c(NA, "D", NA, "D"),
                                sire = c(NA, NA, "S", "S"),
                                applicable = rep(TRUE, length(objective))) {
  data.frame(dam = as.character(dam), sire = as.character(sire),
    nparents = as.integer(!is.na(dam)) + as.integer(!is.na(sire)),
    objective_ln = as.numeric(objective), comparison_complete = applicable,
    relationship_gate_pass = FALSE, hard_feasible = TRUE,
    state_active = TRUE,
    trio_me_pass = TRUE, eligible = FALSE, stringsAsFactors = FALSE)
}

.zero_contract_finish <- function(z, threshold = 1, margin = .1)
  PaPeR:::.search_finish_configuration_contrasts(z, threshold, margin)

test_that("only a zero-over-zero likelihood ratio is undefined", {
  contrast <- PaPeR:::.log_likelihood_contrast_ln
  result <- contrast(c(-3, -Inf, -3, -Inf), c(-5, -5, -Inf, -Inf))
  expect_identical(result, c(2, -Inf, Inf, NA_real_))
  expect_false(any(is.nan(result)))
  expect_identical(contrast(-Inf, c(-5, -Inf)), c(-Inf, NA_real_))
  expect_identical(contrast(c(-5, -Inf), -Inf), c(Inf, NA_real_))
  expect_identical(contrast(setNames(c(-5, -Inf), c("a", "b")), -Inf),
    setNames(c(Inf, NA_real_), c("a", "b")))
  for (gauge in c(-1000, -2.75, 10))
    expect_equal(contrast(-3 + gauge, -5 + gauge), 2, tolerance = 1e-12)
  for (bad in list(NA_real_, NaN, Inf, "zero", structure(0, class = "score")))
    expect_error(contrast(bad, 0), "finite or -Inf")
  expect_error(contrast(1:2, 1:3), "equal lengths or length one")
})

.zero_contract_pool_fixture <- function() {
  calls <- matrix("1/1", 3L, 1L, dimnames = list(c("O", "D", "S"), "L1"))
  # Strictly positive supplied frequencies; squared rare genotype mass may
  # underflow in this unconditional screen. This does not assert prevalence.
  tables <- build_tables(list(L1 = c(1e-200, 1)), Err = 1e-200)
  ids <- rownames(calls)
  lh <- PaPeR:::.normalize_life_history(data.frame(ID = ids,
    Sex = c(3L, 1L, 2L), BirthYear = c(5L, 0L, 0L)), ids)
  ctx <- PaPeR:::.new_inference_context(calls, tables, "simp")
  universe <- PaPeR:::.relationship_universe("simp", "PO", operation = "parent")
  pm <- po_mismatch_prob(tables)
  oh <- PaPeR:::.calibrate_screen_threshold(pm$p_fm, .999, "OH", threshold = 99L)
  me <- PaPeR:::.calibrate_screen_threshold(pm$p_me_true_trio, .999,
    "MaxMendelE", threshold = 99L)
  me$p_me_true_trio <- pm$p_me_true_trio
  list(calls = calls, tables = tables, ids = ids, lh = lh, ctx = ctx,
    universe = universe, oh = oh, me = me)
}

test_that("an oversized packed window counts only defined cheap rejections", {
  ids <- c("O", sprintf("D%02d", 1:51))
  calls <- matrix("1/1", length(ids), 1L, dimnames = list(ids, "L1"))
  tables <- build_tables(list(L1 = c(1e-200, 1)), Err = 1e-200)
  ctx <- PaPeR:::.new_inference_context(calls, tables, "simp")
  lh <- PaPeR:::.normalize_life_history(data.frame(ID = ids,
    Sex = c(3L, rep(1L, 51L)), BirthYear = c(5L, rep(0L, 51L))), ids)
  universe <- PaPeR:::.relationship_universe("simp", "PO", operation = "parent")
  pm <- po_mismatch_prob(tables)
  oh <- PaPeR:::.calibrate_screen_threshold(pm$p_fm, .999, "OH", threshold = 99L)
  me <- PaPeR:::.calibrate_screen_threshold(pm$p_me_true_trio, .999,
    "MaxMendelE", threshold = 99L)
  me$p_me_true_trio <- pm$p_me_true_trio
  pool <- PaPeR:::.build_candidate_pool(ctx, lh, ids,
    c(M = 1L, P = 1L), c(M = 10L, P = 10L), 99L, -2, universe, oh, me)
  focal <- which(pool$rows$offspring == "O" & pool$rows$side == "dam")
  expect_length(focal, 51L)
  expect_true(all(is.na(pool$rows$tfilter_pass[focal])))
  expect_true(all(pool$rows$permanent_pass[focal]))
  expect_false(anyNA(pool$rows$candidate))
  expect_equal(pool$telemetry$tfilter_examined, 51)
  expect_equal(pool$telemetry$tfilter_rejected, 0)
  shortlist <- pool$shortlist_rows[[PaPeR:::.candidate_pool_key("O", "dam")]]
  expect_equal(pool$rows$candidate[shortlist], ids[2:51])
})
