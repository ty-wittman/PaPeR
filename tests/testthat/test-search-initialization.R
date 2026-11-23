.section5_pool_fixture <- function(calls, T = 99L, MaxMendelE = 99L,
                                   Tfilter = -Inf, year_last = NULL) {
  x <- as.matrix(calls)
  if (is.null(colnames(x))) colnames(x) <- paste0("L", seq_len(ncol(x)))
  tables <- build_tables(rep(list(c(.5, .5)), ncol(x)), Err = .01)
  names(tables) <- colnames(x)
  ids <- rownames(x)
  if (is.null(year_last)) year_last <- rep(NA_integer_, length(ids))
  sex <- setNames(rep(3L, length(ids)), ids)
  sex[intersect(ids, c("D", "X"))] <- 1L
  sex[intersect(ids, c("S", "Y"))] <- 2L
  lh <- PaPeR:::.normalize_life_history(data.frame(
    ID = ids, Sex = unname(sex), BirthYear = c(O = 2005L,
      setNames(rep(2000L, length(ids) - 1L), setdiff(ids, "O")))[ids],
    Year.last = year_last, stringsAsFactors = FALSE), ids)
  ctx <- PaPeR:::.new_inference_context(x, tables, "simp")
  universe <- PaPeR:::.relationship_universe(
    "simp", "PO", operation = "parent")
  pm <- po_mismatch_prob(tables)
  oh <- PaPeR:::.calibrate_screen_threshold(
    pm$p_fm, .999, "OH", threshold = T)
  me <- PaPeR:::.calibrate_screen_threshold(
    pm$p_me_true_trio, .999, "MaxMendelE", threshold = MaxMendelE)
  me$p_me_true_trio <- pm$p_me_true_trio
  pool <- PaPeR:::.build_candidate_pool(
    ctx, lh, ids, c(M = 1L, P = 1L), c(M = 10L, P = 10L),
    T, Tfilter, universe, oh, me)
  list(x = x, tables = tables, ids = ids, lh = lh, ctx = ctx,
       universe = universe, pool = pool)
}

.section5_initialize <- function(fixture, pedigree = NULL) {
  PaPeR:::.initialize_parentage_state(
    fixture$ids, genotype_rows = seq_along(fixture$ids),
    life_history = fixture$lh, min_age_parent = c(1L, 1L),
    max_age_parent = c(10L, 10L), pedigree = pedigree,
    ctx = fixture$ctx, candidate_pool = fixture$pool)
}

test_that("empty, one-parent, and two-parent supplied scaffolds initialize", {
  calls <- rbind(O = rep("1/2", 4L), D = rep("1/1", 4L),
                 S = rep("2/2", 4L), X = rep("1/1", 4L))
  f <- .section5_pool_fixture(calls)
  empty <- .section5_initialize(f, data.frame(
    id = character(), dam = character(), sire = character()))
  expect_equal(empty$parents, matrix(0L, 4L, 2L,
    dimnames = list(f$ids, c("dam", "sire"))))
  one <- .section5_initialize(f, data.frame(
    id = "O", dam = "D", sire = NA_character_))
  expect_identical(PaPeR:::.current_parent_configuration(one, "O"),
                   c(dam = "D", sire = NA_character_))
  two <- .section5_initialize(f, data.frame(id = "O", dam = "D", sire = "S"))
  expect_identical(PaPeR:::.current_parent_configuration(two, "O"),
                   c(dam = "D", sire = "S"))
  expect_identical(two$candidate_pool, f$pool)
})

test_that("supplied structural conflicts fail explicitly without edge deletion", {
  calls <- rbind(O = rep("1/2", 2L), D = rep("1/1", 2L),
                 S = rep("2/2", 2L))
  f <- .section5_pool_fixture(calls)
  expect_error(.section5_initialize(f, data.frame(
    id = "O", dam = "S", sire = NA_character_)),
    "known_sex_role_conflict")
  expect_error(.section5_initialize(f, data.frame(
    id = "O", dam = "D", sire = "D")), "same_parent_both_sides")
  expect_error(.section5_initialize(f, data.frame(
    id = "O", dam = "O", sire = NA_character_)), "self_parentage")
  expect_error(.section5_initialize(f, data.frame(
    id = c("O", "D"), dam = c("D", "O"), sire = NA_character_)),
    "ancestry_cycle")
  expect_error(.section5_initialize(f, data.frame(
    id = c("O", "O"), dam = c("D", NA), sire = NA_character_)),
    "unique nonempty")

  late_lh <- f$lh
  late_lh$BirthYear[late_lh$ID == "D"] <- 2004L
  late_lh$BY.min[late_lh$ID == "D"] <- 2004L
  late_lh$BY.max[late_lh$ID == "D"] <- 2004L
  expect_error(PaPeR:::.initialize_parentage_state(
    f$ids, life_history = late_lh, min_age_parent = 2L, max_age_parent = 10L,
    pedigree = data.frame(id = "O", dam = "D", sire = NA_character_)),
    "impossible_age_direction")

  last_lh <- f$lh
  last_lh$Year.last[last_lh$ID == "D"] <- 2004L
  expect_error(PaPeR:::.initialize_parentage_state(
    f$ids, life_history = last_lh, min_age_parent = 1L, max_age_parent = 10L,
    pedigree = data.frame(id = "O", dam = "D", sire = NA_character_)),
    "year_last_conflict")
})

test_that("OH, ME, permanent and retained-universe gates remain without duplicate screening", {
  oh <- .section5_pool_fixture(rbind(
    O = rep("2/2", 3L), D = rep("1/1", 3L), S = rep("1/2", 3L)), T = 0L)
  # The compact directional pool omits OH-failing rows; supplied validation
  # therefore reports the retained-universe boundary first.
  expect_error(.section5_initialize(oh, data.frame(
    id = "O", dam = "D", sire = NA_character_)), "candidate_absent")

  me <- .section5_pool_fixture(rbind(
    O = rep("1/1", 3L), D = rep("1/2", 3L), S = rep("2/2", 3L)),
    MaxMendelE = 0L)
  expect_error(.section5_initialize(me, data.frame(
    id = "O", dam = "D", sire = "S")), "mendelian_error_threshold")

  duplicate <- .section5_pool_fixture(rbind(
    O = rep("1/2", 8L), D = rep("1/2", 8L), S = rep("2/2", 8L)))
  # Duplicate evidence cannot invalidate a supplied, biologically legal parent.
  for (flag in list(TRUE, NA, NULL)) {
    duplicate$pool$rows$duplicate_veto <- flag
    state <- .section5_initialize(duplicate, data.frame(
      id = "O", dam = "D", sire = NA_character_))
    expect_identical(PaPeR:::.current_parent_configuration(state, "O")[["dam"]], "D")
  }

  weak_calls <- rbind(
    O = rep("1/2", 3L), D = rep("1/1", 3L), S = rep("2/2", 3L),
    matrix(rep("1/1", 3L * 52L), nrow = 52L, byrow = TRUE))
  rownames(weak_calls) <- c("O", "D", "S", sprintf("X%02d", 1:52))
  weak <- .section5_pool_fixture(weak_calls,
    Tfilter = 1e6)
  expect_gt(weak$pool$telemetry$tfilter_examined, 50L)
  expect_gt(weak$pool$telemetry$tfilter_rejected, 50L)
  expect_error(.section5_initialize(weak, data.frame(
    id = "O", dam = "D", sire = NA_character_)),
    "candidate_absent")

  absent <- .section5_pool_fixture(rbind(
    O = rep("1/2", 3L), D = rep("1/1", 3L), S = rep("2/2", 3L)))
  absent$pool$rows <- absent$pool$rows[!(
    absent$pool$rows$offspring == "O" & absent$pool$rows$side == "dam" &
      absent$pool$rows$candidate == "D"), , drop = FALSE]
  # The compact directional pool omits OH-failing rows; supplied validation
  # therefore reports the retained-universe boundary first.
  expect_error(.section5_initialize(absent, data.frame(
    id = "O", dam = "D", sire = NA_character_)), "candidate_absent")
})

