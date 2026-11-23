test_that("public final LLR agrees with tiny pedigree enumeration", {
  calls <- rbind(D = c("1/1", "1/2"), A = c("1/2", "2/2"))
  colnames(calls) <- c("L1", "L2")
  freq <- c(.6, .4); dropout <- .01; foreign <- .001
  geno <- rbind(c(1L, 1L), c(1L, 2L), c(2L, 2L))
  prior <- apply(geno, 1L, function(g) prod(freq[g]) * if (g[1] == g[2]) 1 else 2)
  transmission <- array(0, c(3L, 3L, 3L))
  for (c in 1:3) for (d in 1:3) for (s in 1:3)
    for (a in 1:2) for (b in 1:2)
      transmission[c, d, s] <- transmission[c, d, s] + .25 *
        identical(sort(c(geno[d, a], geno[s, b])), geno[c, ])
  emission <- function(g, observed) {
    precursor <- if (g[1] == g[2]) list(g) else list(g, rep(g[1], 2), rep(g[2], 2))
    weight <- if (g[1] == g[2]) 1 else c(1 - dropout, dropout / 2, dropout / 2)
    r <- function(a, b) (1 - foreign) * (a == b) + foreign * freq[b]
    sum(vapply(seq_along(precursor), function(i) {
      p <- precursor[[i]]; value <- r(p[1], observed[1]) * r(p[2], observed[2])
      if (observed[1] != observed[2]) value <- value + r(p[1], observed[2]) * r(p[2], observed[1])
      weight[i] * value
    }, numeric(1)))
  }
  founder <- c(0L, 0L)
  graphs <- list(
    U = list(p = list(founder, founder), typed = c(1L, 2L)),
    PO = list(p = list(founder, founder, c(1L, 2L)), typed = c(1L, 3L)),
    FS = list(p = list(founder, founder, c(1L, 2L), c(1L, 2L)), typed = c(3L, 4L)),
    HS = list(p = list(founder, founder, founder, c(1L, 2L), c(1L, 3L)), typed = c(4L, 5L)),
    GP = list(p = list(founder, founder, founder, c(1L, 2L), c(4L, 3L)), typed = c(1L, 5L)),
    FA = list(p = list(founder, founder, founder, c(1L, 2L), c(1L, 2L), c(5L, 3L)), typed = c(4L, 6L)),
    HA = list(p = list(founder, founder, founder, founder, c(1L, 2L), c(1L, 3L), c(6L, 4L)), typed = c(5L, 7L)))
  score <- function(graph) {
    states <- expand.grid(rep(list(1:3), length(graph$p)))
    weight <- rep(1, nrow(states))
    for (v in seq_along(graph$p)) {
      parents <- graph$p[[v]]
      weight <- weight * if (parents[1] == 0L) prior[states[[v]]] else
        transmission[cbind(states[[v]], states[[parents[1]]], states[[parents[2]]])]
    }
    sum(vapply(seq_len(ncol(calls)), function(locus) {
      likelihood <- weight
      for (typed in 1:2) {
        observed <- as.integer(strsplit(calls[typed, locus], "/", fixed = TRUE)[[1]])
        e <- apply(geno, 1L, emission, observed = observed)
        likelihood <- likelihood * e[states[[graph$typed[typed]]]]
      }
      log(sum(likelihood))
    }, numeric(1)))
  }
  values <- vapply(graphs, score, numeric(1))
  expected <- (values["PO"] - max(values[names(values) != "PO"])) / log(10)
  ped <- data.frame(id = c("D", "A"), dam = c(NA, "D"), sire = NA_character_)
  life <- data.frame(ID = c("D", "A"), Sex = c(1L, 2L), BirthYear = c(0L, 10L))
  got <- calc_llr(ped, calls, life, dropout, foreign,
    control = list(frequencies = rep(list(freq), 2), llr_evidence = "full", hypotheses = "S7"))
  expect_equal(got$LLRdam[2], unname(expected), tolerance = 1e-10)
  expect_identical(got[, c("id", "dam", "sire")], ped)
})
