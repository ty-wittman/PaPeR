retention_policy_fixture <- function() {
  x <- rbind(
    O = c("1/3", "1/3", "2/2"),
    D1 = c("1/1", "1/1", "2/3"),
    D2 = c("2/2", "2/2", "1/1"),
    S1 = c("1/2", "3/3", "1/3"),
    S2 = c("2/3", "2/2", "3/3"))
  colnames(x) <- paste0("L", 1:3)
  tables <- build_tables(rep(list(c(.5, .3, .2)), 3L), Err = .02)
  lh <- data.frame(ID = rownames(x), Sex = c(3L, 1L, 1L, 2L, 2L),
                   BirthYear = c(2005L, 2000L, 2000L, 2000L, 2000L))
  list(GenoM = x, LifeHistData = lh, tables = tables)
}

retention_policy_run <- function(GenoM, LifeHistData, tables, ...) {
  args <- c(list(
    GenoM = GenoM, LifeHistData = LifeHistData, tables = tables,
    Tfilter = -Inf, Tassign = -Inf, margin = -Inf, quiet = TRUE,
    Complex = "simp", .phase = "candidate_only"), list(...))
  do.call(assign_parents, args)
}

test_that("heterogeneous calibration preserves independently enumerated mass", {
  probs <- c(.01, .2, .4, .03)
  states <- as.matrix(expand.grid(rep(list(0:1), length(probs))))
  mass <- apply(states, 1L, function(s) prod(ifelse(s == 1L, probs, 1-probs)))
  expected <- vapply(0:length(probs), function(k) sum(mass[rowSums(states) == k]), numeric(1))
  target <- 0.9999^(1/600)
  for (cal in list(calibrate_maxmismatch(probs, target),
                  calibrate_max_mendelian_error(probs, target))) {
    expect_equal(cal$pmf, expected, tolerance = 1e-15)
    expect_identical(cal$T, as.integer(which(cumsum(expected) >= target)[1L]-1L))
  }
})

