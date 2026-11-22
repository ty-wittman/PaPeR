# Stage-1 likelihood core: invariants, scalar error model, and
# cross-checks against brute force and (if available) CKMRsim.

gidx <- function(a, b, k) { if (a > b) { t <- a; a <- b; b <- t }
  off <- (a - 1) * k - (a - 1) * (a - 2) / 2; off + (b - a) + 1 }

# Test-only dense observed-given-true oracle.  The production locus table uses
# only its compact scalar/profile kernel and never retains this matrix.
dense_emission_oracle_core <- function(tab) {
  out <- matrix(0, tab$G, tab$G)
  for (true_g in seq_len(tab$G)) {
    true <- tab$geno[true_g, ]
    for (observed_g in seq_len(tab$G)) {
      observed <- tab$geno[observed_g, ]
      if (true_g == observed_g) {
        out[true_g, observed_g] <- 1 - tab$Err
      } else if (tab$Err > 0) {
        overlap <- sum(pmin(tabulate(true, nbins = tab$k),
                            tabulate(observed, nbins = tab$k)))
        cls <- if (overlap == 0L) "w_Far" else if (observed[1L] == observed[2L]) {
          "w_Hom"
        } else {
          "w_Het"
        }
        Z <- if (true[1L] == true[2L]) tab$Z_hom else tab$Z_het
        out[true_g, observed_g] <- tab$Err * tab$ErrProfile[[cls]] / Z
      }
    }
  }
  out
}

# Independent small-state founder and transmission oracle for §2.3.  These
# functions deliberately do not use the package's genotype tables or C++ code.
unordered_genotypes_ref <- function(k) {
  do.call(rbind, lapply(seq_len(k), function(a) cbind(a, a:k)))
}

founder_probability_ref <- function(g, p) {
  g <- sort(as.integer(g))
  if (g[1L] == g[2L]) p[g[1L]]^2 else 2 * p[g[1L]] * p[g[2L]]
}

gamete_probability_ref <- function(x, parent) {
  0.5 * (x == parent[1L]) + 0.5 * (x == parent[2L])
}

two_parent_probability_ref <- function(offspring, dam, sire) {
  offspring <- sort(as.integer(offspring))
  a <- offspring[1L]; b <- offspring[2L]
  if (a == b) {
    gamete_probability_ref(a, dam) * gamete_probability_ref(a, sire)
  } else {
    gamete_probability_ref(a, dam) * gamete_probability_ref(b, sire) +
      gamete_probability_ref(b, dam) * gamete_probability_ref(a, sire)
  }
}

one_parent_probability_ref <- function(offspring, parent, p) {
  offspring <- sort(as.integer(offspring))
  a <- offspring[1L]; b <- offspring[2L]
  if (a == b) {
    gamete_probability_ref(a, parent) * p[a]
  } else {
    gamete_probability_ref(a, parent) * p[b] +
      gamete_probability_ref(b, parent) * p[a]
  }
}

two_parent_probability_bruteforce <- function(offspring, dam, sire) {
  offspring <- sort(as.integer(offspring))
  mean(vapply(dam, function(a) {
    mean(vapply(sire, function(b) identical(sort(c(a, b)), offspring), logical(1)))
  }, numeric(1)))
}

one_parent_probability_bruteforce <- function(offspring, parent, p) {
  offspring <- sort(as.integer(offspring))
  sum(vapply(parent, function(a) {
    0.5 * sum(p * vapply(seq_along(p), function(b)
      identical(sort(c(a, b)), offspring), logical(1)))
  }, numeric(1)))
}

test_that("per-locus tables satisfy structural invariants", {
  p <- c(0.5, 0.3, 0.2); Err <- 0.02
  tab <- PaPeR:::build_locus_table(p, Err)
  expect_equal(sum(tab$pi), 1)                              # HWE prior normalized
  expect_true(all(abs(rowSums(dense_emission_oracle_core(tab)) - 1) < 1e-12))
  expect_equal(sum(tab$mObs), 1)
  for (nm in rownames(compact_legacy_delta()))       # outbred: marginal is pi
    expect_true(max(abs(rowSums(dense_relationship_oracle(
      tab, compact_legacy_delta()[nm, ])) - tab$pi)) < 1e-12)
})

test_that("HWE founder priors and one- and two-parent transmission normalize", {
  set.seed(2301)
  spectra <- c(
    list(c(1), c(0.5, 0.5), c(0.7, 0.2, 0.1),
         c(0.35, 0.25, 0.15, 0.1, 0.08, 0.07)),
    lapply(c(2L, 3L, 6L), function(k) {
      x <- rgamma(k, shape = 0.7, rate = 1)
      x / sum(x)
    })
  )

  for (p in spectra) {
    k <- length(p)
    genotypes <- unordered_genotypes_ref(k)
    pi_ref <- apply(genotypes, 1L, founder_probability_ref, p = p)
    expect_equal(sum(pi_ref), 1, tolerance = 1e-14)

    tab <- if (k == 1L) PaPeR:::build_locus_table(p,
      ErrorModel = "dropout_foreign", Dropout = .01, Foreign = .001) else
      PaPeR:::build_locus_table(p, Err = 1e-12)
    expect_equal(unname(tab$geno), unname(genotypes))
    expect_equal(tab$pi, pi_ref, tolerance = 1e-14)

    for (dam_i in seq_len(nrow(genotypes))) {
      dam <- genotypes[dam_i, ]
      one_parent <- apply(genotypes, 1L, one_parent_probability_ref,
                          parent = dam, p = p)
      expect_equal(sum(one_parent), 1, tolerance = 1e-14)
      expect_equal(one_parent,
                   apply(genotypes, 1L, one_parent_probability_bruteforce,
                         parent = dam, p = p), tolerance = 1e-14)

      for (sire_i in seq_len(nrow(genotypes))) {
        sire <- genotypes[sire_i, ]
        two_parent <- apply(genotypes, 1L, two_parent_probability_ref,
                            dam = dam, sire = sire)
        expect_equal(sum(two_parent), 1, tolerance = 1e-14)
        expect_equal(two_parent,
                     apply(genotypes, 1L, two_parent_probability_bruteforce,
                           dam = dam, sire = sire), tolerance = 1e-14)
      }
    }
  }
})

test_that("Mendelian kernels preserve allele order and impossible states", {
  p2 <- c(0.75, 0.25)
  genotypes2 <- unordered_genotypes_ref(2L)
  expect_equal(apply(genotypes2, 1L, founder_probability_ref, p = p2),
               c(0.5625, 0.375, 0.0625))
  expect_equal(apply(genotypes2, 1L, two_parent_probability_ref,
                     dam = c(1L, 2L), sire = c(1L, 2L)),
               c(0.25, 0.5, 0.25))
  expect_equal(apply(genotypes2, 1L, one_parent_probability_ref,
                     parent = c(1L, 2L), p = p2),
               c(0.375, 0.5, 0.125))

  p <- c(0.6, 0.3, 0.1)
  dam <- c(1L, 2L); sire <- c(3L, 3L)
  expect_equal(two_parent_probability_ref(c(1L, 3L), dam, sire), 0.5)
  expect_equal(one_parent_probability_ref(c(1L, 3L), dam, p), 0.05)
  expect_equal(two_parent_probability_ref(c(3L, 3L), c(1L, 1L), c(2L, 2L)), 0)
  expect_equal(one_parent_probability_ref(c(3L, 3L), c(1L, 1L), p), 0)
  expect_equal(two_parent_probability_ref(c(3L, 1L), rev(dam), rev(sire)),
               two_parent_probability_ref(c(1L, 3L), dam, sire))
  expect_equal(one_parent_probability_ref(c(3L, 1L), rev(dam), p),
               one_parent_probability_ref(c(1L, 3L), dam, p))
})

test_that("compiled trio terms equal an independent positive-error transmission oracle", {
  p <- c(0.6, 0.3, 0.1); k <- length(p)
  tab <- PaPeR:::build_locus_table(p, Err = 1e-12)
  genotypes <- tab$geno
  emissions <- dense_emission_oracle_core(tab)
  fixtures <- list(
    c(o = gidx(1, 3, k), d = gidx(1, 2, k), s = gidx(3, 3, k)),
    c(o = gidx(2, 2, k), d = gidx(1, 2, k), s = gidx(2, 3, k)),
    c(o = gidx(3, 3, k), d = gidx(1, 1, k), s = gidx(2, 2, k))
  )

  for (fixture in fixtures) {
    o <- genotypes[fixture[["o"]], ]
    d <- genotypes[fixture[["d"]], ]
    s <- genotypes[fixture[["s"]], ]
    public <- matrix(vapply(list(o, d, s), function(g)
      paste(sort(g), collapse = "/"), character(1)), ncol = 1L,
      dimnames = list(c("o", "d", "s"), "MH1"))
    prepared <- prepare_genotypes(public, n_alleles = k)
    got <- PaPeR:::calc_trio_LL_cpp(list(tab), prepared, 0L, 1L, 2L)
    expected <- setNames(rep(0,4),c('H_PP','H_D0','H_0S','H_U'))
    for (oi in seq_len(tab$G)) for (di in seq_len(tab$G)) for (si in seq_len(tab$G)) {
      ot<-genotypes[oi,];dt<-genotypes[di,];st<-genotypes[si,]
      mass<-founder_probability_ref(dt,p)*founder_probability_ref(st,p)*
        emissions[oi,fixture[['o']]]*emissions[di,fixture[['d']]]*emissions[si,fixture[['s']]]
      expected<-expected+mass*c(two_parent_probability_ref(ot,dt,st),
        one_parent_probability_ref(ot,dt,p),one_parent_probability_ref(ot,st,p),
        founder_probability_ref(ot,p))
    }
    expect_equal(got, log(expected), tolerance = 1e-14)
  }
})

test_that("non-inbred Jacquard relationships match CKMRsim make_matrix_X_l", {
  skip_if_not_installed("CKMRsim")
  mmX <- get("make_matrix_X_l", getNamespace("CKMRsim"))
  p <- c(0.5, 0.30, 0.20)
  tab <- PaPeR:::build_locus_table(p, Err = 1e-12)             # X is error-free
  kaps <- list(U = c(1,0,0), PO = c(0,1,0), FS = c(.25,.5,.25),
               DEG2 = c(.5,.5,0), HA = c(.75,.25,0), MZ = c(0,0,1))
  for (nm in names(kaps))
    expect_true(max(abs(dense_relationship_oracle(tab, compact_legacy_delta()[nm, ]) -
                          mmX(p, kaps[[nm]]))) < 1e-12)
})

test_that("inbreeding-aware relationships: proper joint and inbred marginal", {
  p <- c(0.5, 0.30, 0.20)
  tab <- PaPeR:::build_locus_table(p, Err = 1e-12)
  # every retained relationship X is a proper joint with HWE marginal.
  for (nm in rownames(compact_legacy_delta()))
    expect_lt(abs(sum(dense_relationship_oracle(tab, compact_legacy_delta()[nm, ])) - 1), 1e-12)
  expect_lt(max(abs(rowSums(dense_relationship_oracle(
    tab, compact_legacy_delta()["PO", ])) - tab$pi)), 1e-12)
})

test_that("pairwise Y lookup equals brute-force sum over true genotypes", {
  p <- c(0.4,0.3,0.2,0.1); k <- length(p)
  tab <- PaPeR:::build_locus_table(p, Err = 0.02)
  E <- dense_emission_oracle_core(tab)
  brute <- function(o1, o2) sapply(rownames(compact_legacy_delta()),
    function(nm) sum(outer(E[,o1], E[,o2]) *
                       dense_relationship_oracle(tab, compact_legacy_delta()[nm, ])))
  for (rep in 1:20) {
    o1 <- sample(tab$G, 1); o2 <- sample(tab$G, 1)
    y <- sapply(rownames(compact_legacy_delta()), function(nm)
      dense_observed_relationship_oracle(tab, compact_legacy_delta()[nm, ])[o1, o2])
    expect_lt(max(abs(y - brute(o1, o2))), 1e-12)
  }
})

test_that("trio LL matches an independent brute-force implementation", {
  set.seed(9); k <- 6; p <- rgamma(k, 0.8, 1); p <- (p + 1e-3) / sum(p + 1e-3)
  tab <- PaPeR:::build_locus_table(p, Err = 0.02)
  brute <- function(oO, oD, oS) {
    G <- tab$G; pi <- tab$pi; E <- dense_emission_oracle_core(tab)
    mObs <- tab$mObs; pp <- tab$p; geno <- tab$geno
    L_PP <- 0; L_D0 <- 0
    for (gD in 1:G) { fD <- pi[gD]*E[gD,oD]; if (fD == 0) next; cD <- geno[gD,]
      inner <- 0
      for (gS in 1:G) { fS <- pi[gS]*E[gS,oS]; if (fS == 0) next; cS <- geno[gS,]
        W <- 0.25*(E[gidx(cD[1],cS[1],k),oO]+E[gidx(cD[1],cS[2],k),oO]+
                   E[gidx(cD[2],cS[1],k),oO]+E[gidx(cD[2],cS[2],k),oO])
        inner <- inner + fS*W }
      L_PP <- L_PP + fD*inner
      W1 <- 0.5*sum(sapply(cD, function(c) sum(pp*sapply(1:k, function(x) E[gidx(c,x,k),oO]))))
      L_D0 <- L_D0 + fD*W1 }
    c(H_PP = log(L_PP), H_D0 = log(L_D0*mObs[oS]), H_U = log(mObs[oO]*mObs[oD]*mObs[oS]))
  }
  for (rep in 1:30) {
    oO <- sample(tab$G,1); oD <- sample(tab$G,1); oS <- sample(tab$G,1)
    raw <- matrix(c(tab$geno[oO,], tab$geno[oD,], tab$geno[oS,]), 3, 2, byrow = TRUE)
    pub <- matrix(apply(raw, 1, function(x) paste(sort(x), collapse = "/")),
                  ncol = 1, dimnames = list(c("o", "d", "s"), "MH1"))
    Gm <- prepare_genotypes(pub, k)
    cpp <- PaPeR:::calc_trio_LL_cpp(list(tab), Gm, 0, 1, 2)
    br <- brute(oO, oD, oS)
    expect_lt(max(abs(cpp[c("H_PP","H_D0","H_U")] - br)), 1e-10)
  }
})
