test_that("locus tables use the compact schema and reject dense legacy fields", {
  p <- c(0.5, 0.3, 0.2)
  tab <- PaPeR:::build_locus_table(p, Err = 0.01)
  expect_identical(tab$schema_version, "compact-locus-v1")
  expect_true(all(c("k", "G", "geno", "p", "pi", "mObs", "Err",
                    "ErrProfile", "Z_hom", "Z_het") %in% names(tab)))
  expect_false(any(c("S", "YS", "X", "Delta", "relnames") %in% names(tab)))
  expect_false(any(vapply(Filter(is.numeric, tab), function(x)
    is.matrix(x) && identical(dim(x), c(tab$G, tab$G)), logical(1))))
  legacy <- tab
  legacy$S <- vector("list", 9L)
  expect_error(PaPeR:::.new_inference_context(
    matrix("1/1", nrow = 1L, dimnames = list("a", "L1")), list(L1 = legacy)),
    "obsolete dense")
})

test_that("inference context prepares genotype indices once and locks inputs", {
  x <- matrix(c("1/1", "1/2", "0", "1/2", "2/2", "1/1",
                "2/2", "0", "1/2"), nrow = 3L, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), paste0("L", 1:3)))
  tables <- build_tables(rep(list(c(0.55, 0.45)), 3L), Err = 0.01)
  names(tables) <- colnames(x)
  ctx <- PaPeR:::.new_inference_context(x, tables, "simp")
  expect_s3_class(ctx, "PaPeR_inference_context")
  expect_equal(dim(ctx$obsidx), c(3L, 3L))
  expect_identical(ctx$obsidx[1L, 3L], -1L)
  expect_identical(ctx$obsidx[3L, 2L], -1L)
  expect_error(ctx$Delta <- matrix(0, 1, 9), "locked binding")
  expect_error(ctx$likelihood_cache <- new.env(), "locked binding")
  expect_equal(ctx$unique_observed_states[[1L]], 0:2)
})

test_that("compiled batch scoring preserves duplicate requests and orientations", {
  x <- matrix(c("1/1", "1/2", "2/2", "1/2", "1/1", "2/2"), nrow = 3L,
              byrow = TRUE, dimnames = list(c("a", "b", "c"), c("L1", "L2")))
  tables <- build_tables(list(L1 = c(.6, .4), L2 = c(.55, .45)), Err = .01)
  ctx <- PaPeR:::.new_inference_context(x, tables, "simp")
  batch <- PaPeR:::.score_unconditional_pairs(
    ctx, c("a", "a", "b", "a"), c("b", "b", "c", "b"))
  raw_batch <- PaPeR:::score_basis_pairs_cpp(
    ctx$markers, ctx$obsidx, c(0L, 0L, 1L, 0L), c(1L, 1L, 2L, 1L), ctx$Delta)
  counters <- attr(raw_batch, "compact_batch_counters")
  expect_equal(counters$requested_pairs, 4)
  expect_equal(counters$unique_individual_pairs, 2)
  expect_equal(counters$requested_pair_locus_cells, 8)
  expect_lte(counters$unique_ordered_observed_state_pairs, 4)
  expect_equal(counters$reused_pair_locus_cells,
               8 - counters$unique_ordered_observed_state_pairs)
  expect_equal(batch[1L, ], batch[2L, ], tolerance = 0)
  expect_equal(batch[1L, ], batch[4L, ], tolerance = 0)
  scalar_ab <- PaPeR:::score_basis_pairs_cpp(
    ctx$markers, ctx$obsidx, 0L, 1L, ctx$Delta)
  expect_equal(unname(batch[1L, ]), unname(scalar_ab[1L, ]), tolerance = 1e-12)
  reversed <- PaPeR:::.score_unconditional_pairs(ctx, "b", "a")
  scalar_ba <- PaPeR:::score_basis_pairs_cpp(
    ctx$markers, ctx$obsidx, 1L, 0L, ctx$Delta)
  expect_equal(unname(reversed[1L, ]), unname(scalar_ba[1L, ]), tolerance = 1e-12)
})

