.duplicate_report_fixture <- function() {
  calls <- rbind(A = c("1/2", "1/1", "0", "2/2"),
                 B = c("2/1", "2/2", "1/1", "2/2"),
                 COPY = c("1/2", "1/1", "0", "2/2"),
                 MISSING = rep("0", 4L))
  colnames(calls) <- paste0("L", 1:4)
  list(calls = calls, tables = build_tables(setNames(rep(list(c(.6, .4)), 4L),
                                                   colnames(calls)), Err = .01))
}

test_that("explicit duplicate reports have likelihood contrasts and missing-aware QC", {
  f <- .duplicate_report_fixture(); before <- serialize(f, NULL)
  for (mode in c("simp", "full")) {
    d <- .find_duplicates(f$calls, f$tables, thresh = -Inf, Complex = mode)
    expect_equal(nrow(d), choose(nrow(f$calls), 2))
    ab <- d[d$ID1 == "A" & d$ID2 == "B", ]
    expect_equal(ab$jointly_called, 3L)
    expect_equal(ab$genotype_mismatches, 1L)
    expect_equal(ab$opposing_homozygotes, 1L)
    expect_equal(ab$mismatch_rate, 1/3)
    expect_equal(ab$missing_id1, 1L)
    expect_equal(ab$missing_id2, 0L)
    copy <- d[d$ID1 == "A" & d$ID2 == "COPY", ]
    expect_equal(copy$genotype_mismatches, 0L)
    missing <- d[d$ID2 == "MISSING", ]
    expect_true(all(missing$jointly_called == 0L))
    expect_true(all(is.na(missing$mismatch_rate)))
    ctx <- PaPeR:::.new_inference_context(f$calls, f$tables, mode,
                                             include_duplicate = TRUE)
    ll <- PaPeR:::.score_unconditional_pairs(ctx, d$ID1, d$ID2)
    for (other in c("PO", "FS", "HS", "GP", "FA", "HA", "U"))
      expect_equal(d[[paste0("LLR_MZ_vs_", other)]],
                   unname((ll[, "MZ"] - ll[, other]) / log(10)))
    expect_equal(d$LLR, d$LLR_MZ_vs_U)
    expect_equal(d$LLR_MZ_vs_best, apply(d[, paste0("LLR_MZ_vs_",
      c("PO", "FS", "HS", "GP", "FA", "HA", "U"))], 1L, min))
    expect_identical(attr(d, "log_base"), 10L)
    empty <- .find_duplicates(f$calls, f$tables, thresh = Inf, Complex = mode)
    expect_equal(nrow(empty), 0L)
    expect_identical(names(empty), names(d))
  }
  expect_identical(serialize(f, NULL), before)
})

test_that("ordinary inference has no implicit duplicate screen or MZ requirement", {
  f <- .duplicate_report_fixture()
  for (mode in c("simp", "full")) {
    ctx <- PaPeR:::.new_inference_context(f$calls, f$tables, mode)
    expect_false("MZ" %in% rownames(ctx$Delta))
    ordinary <- rel_edges(f$calls, f$tables, "PO", alts = "U", thresh = -Inf,
                          oh_max = -1L, Complex = mode)
    expect_true(any(ordinary$id1 == "A" & ordinary$id2 == "COPY"))
    expect_true(all(is.na(ordinary$duplicate_guard_llr)))
    expect_true(all(ordinary$duplicate_guard_status == "not_evaluated"))
    maybe <- .get_maybe_rel(f$calls, f$tables, thresh = -Inf, Complex = mode)
    expect_true(any(maybe$ID1 == "A" & maybe$ID2 == "COPY"))
  }
})

test_that("the active duplicate kernel agrees with positive-error shared-genotype probabilities", {
  calls <- matrix(c("1/1", "1/1", "2/2"), ncol = 1L,
                  dimnames = list(c("A", "COPY", "DISCORDANT"), "L1"))
  tables <- build_tables(list(L1 = rep(1 / 3, 3)), Err = 1e-12)
  for (mode in c("simp", "full")) {
    ctx <- PaPeR:::.new_inference_context(calls, tables, mode,
                                             include_duplicate = TRUE)
    ll <- PaPeR:::.score_unconditional_pairs(
      ctx, rep("A", 2L), c("COPY", "DISCORDANT"))
    # One shared founder genotype; each observed call has its own error draw.
    E<-dense_emission_oracle_compact(tables[[1]])
    expected<-vapply(c(1L,4L),function(o)sum(tables[[1]]$pi*E[,1L]*E[,o]),0.0)
    expect_equal(unname(ll[,"MZ"]),log(expected),tolerance=1e-12)
    expect_true(all(is.finite(ll[,"MZ"])))
    expect_lt(ll[2L,"MZ"],ll[1L,"MZ"])
    expect_true(is.finite(ll[2L, "U"]))
  }
})

