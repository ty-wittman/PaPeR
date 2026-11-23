test_that("joint support, not separately normalized raw support, controls assignment", {
  source(test_path("helper-parentage-fixture.R"), local = TRUE)
  prior <- c(.09, .42, .49)
  alleles <- rbind(c(0, 0), c(0, 1), c(1, 1))
  transmission <- array(0, c(3, 3, 3))
  for (d in 1:3) for (s in 1:3) for (a in 1:2) for (b in 1:2) {
    child <- 1 + alleles[d, a] + alleles[s, b]
    transmission[child, d, s] <- transmission[child, d, s] + .25
  }
  one_parent <- apply(transmission, c(1, 2), function(x) sum(x * prior))
  states <- expand.grid(A = 1:3, B = 1:3, C = 1:3)
  threshold_ln <- .2 * log(10)

  # A has an observed offspring C. The exact sum distinguishes the joint
  # evidence contrast from separately normalizing both hypotheses first.
  oracle <- function(f, ancestral) {
    totals <- c(joint = 0, raw = 0)
    for (locus in seq_len(ncol(f$calls))) {
      emission <- exp(f$compact_inputs$observations_ln[[locus]])
      w <- prior[states$B] * emission[cbind(states$B, 2L)]
      if (ancestral) w <- w / 3 else
        w <- w * one_parent[cbind(states$C, states$A)] * emission[cbind(states$C, 3L)]
      unrelated <- w * prior[states$A]
      parent <- w * one_parent[cbind(states$A, states$B)]
      numerator <- sum(parent * emission[cbind(states$A, 1L)])
      null <- sum(unrelated * emission[cbind(states$A, 1L)])
      totals <- totals + c(
        joint = log(numerator / null),
        raw = log(numerator / sum(parent)) - log(null / sum(unrelated)))
    }
    totals
  }

  for (case in 1:2) {
    loci <- if (case == 1L) 3L else 4L
    calls <- rbind(A = rep("1/2", loci),
      B = rep(if (case == 1L) "1/1" else "2/2", loci), C = rep("1/1", loci))
    colnames(calls) <- paste0("L", seq_len(loci))
    f <- .scientific_fixture(calls, dam = c(C = "A"),
      dam_candidates = "B", p = c(.3, .7), Complex = "full",
      birth_year = c(A = 3L, B = 1L, C = 5L))
    full <- oracle(f, ancestral = FALSE)
    expect_identical(unname(full[["joint"]] > threshold_ln), case == 1L)
    expect_identical(unname(full[["raw"]] > threshold_ln), case == 2L)

    result <- .public_assignment(f, control = list(
      pedigree = f$pedigree, hypotheses = "S7", evidence = "full",
      Tassign = .2, margin = 0, Tfilter = -Inf))
    at <- match("A", result$id)
    expect_identical(is.na(result$dam[at]), case == 2L)
    if (case == 1L) expect_identical(result$dam[at], "B")
  }
})

test_that("public assignment follows configuration-local support with ancestral or full evidence", {
  source(test_path("helper-parentage-fixture.R"), local = TRUE)
  calls <- rbind(A = c("1/2", "2/2", "1/1"),
    B = c("1/1", "1/2", "0"), C = c("2/2", "1/2", "1/2"),
    G = c("1/2", "1/1", "1/2"), H = c("1/2", "2/2", "1/1"))
  colnames(calls) <- paste0("L", 1:3)
  freq <- c(.31, .69)
  prior <- c(freq[1]^2, 2 * prod(freq), freq[2]^2)
  f <- .scientific_fixture(calls, dam = c(B = "G", C = "G"),
    sire = c(A = "C", C = "H"), dam_candidates = "B", sire_candidates = "C",
    focal = "A", p = freq, Complex = "full",
    birth_year = c(A = 4L, B = 2L, C = 2L, G = 0L, H = 0L))

  alleles <- rbind(c(0, 0), c(0, 1), c(1, 1))
  transmission <- array(0, c(3, 3, 3))
  for (d in 1:3) for (s in 1:3) for (a in 1:2) for (b in 1:2)
    transmission[1 + alleles[d, a] + alleles[s, b], d, s] <-
      transmission[1 + alleles[d, a] + alleles[s, b], d, s] + .25
  one_parent <- apply(transmission, c(1, 2), function(x) sum(x * prior))
  states <- as.matrix(expand.grid(rep(list(1:3), 5)))
  joint <- local <- numeric(4)
  scopes <- list(1L, c(1L, 2L, 4L), c(1L, 3L, 4L, 5L), 1:5)
  for (locus in seq_len(ncol(calls))) {
    emission <- exp(f$compact_inputs$observations_ln[[locus]])
    w <- prior[states[, 4]] * prior[states[, 5]] *
      one_parent[cbind(states[, 2], states[, 4])] *
      transmission[cbind(states[, 3], states[, 4], states[, 5])]
    for (configuration in 1:4) {
      local_weight <- w
      for (id in scopes[[configuration]])
        local_weight <- local_weight * emission[cbind(states[, id], id)]
      focal_factor <- switch(configuration,
        prior[states[, 1]],
        one_parent[cbind(states[, 1], states[, 2])],
        one_parent[cbind(states[, 1], states[, 3])],
        transmission[cbind(states[, 1], states[, 2], states[, 3])])
      local[configuration] <- local[configuration] +
        log(sum(local_weight * focal_factor)) -
        log(sum(local_weight * prior[states[, 1]]))
    }
    for (id in 1:5) w <- w * emission[cbind(states[, id], id)]
    joint <- joint + log(c(
      sum(w * prior[states[, 1]]),
      sum(w * one_parent[cbind(states[, 1], states[, 2])]),
      sum(w * one_parent[cbind(states[, 1], states[, 3])]),
      sum(w * transmission[cbind(states[, 1], states[, 2], states[, 3])])))
  }

  threshold_ln <- .05 * log(10)
  for (mode in c("ancestral", "full")) {
    expected <- if (mode == "ancestral") local else joint
    support <- c(dam = expected[4] - expected[3],
      sire = expected[4] - expected[2],
      pair = expected[4] - max(expected[1:3]))
    expect_true(all(is.finite(support)))
    expect_true(all(support > threshold_ln))
    result <- .public_assignment(f, control = list(
      pedigree = f$pedigree, hypotheses = "S7", evidence = mode,
      Tassign = .05, margin = 0, Tfilter = -Inf))
    at <- match("A", result$id)
    expect_identical(result$dam[at], "B")
    expect_identical(result$sire[at], "C")
  }
  expect_gt(max(abs((local - local[1]) - (joint - joint[1]))), 1e-4)
})
