# Literal log-domain enumeration, with positive observation error. No native
# likelihood evaluator is used to construct the expected values.
.numeric_positive_trio <- function(tab, obs) {
  geno<-tab$geno;p<-tab$p;G<-tab$G
  E<-dense_emission_oracle_compact(tab)
  pi<-apply(geno,1,function(g)sum(log(p[g]))+if(g[1]!=g[2])log(2) else 0)
  emit<-function(g,o)if(is.na(o))0 else log(E[g,o])
  gamete<-function(a,g)sum(g==a)/2
  po<-function(c,d)if(c[1]==c[2])gamete(c[1],d)*p[c[1]] else
    gamete(c[1],d)*p[c[2]]+gamete(c[2],d)*p[c[1]]
  trio<-function(c,d,s)if(c[1]==c[2])gamete(c[1],d)*gamete(c[1],s) else
    gamete(c[1],d)*gamete(c[2],s)+gamete(c[2],d)*gamete(c[1],s)
  terms<-matrix(-Inf,G^3,4);row<-0L
  for(o in seq_len(G))for(d in seq_len(G))for(s in seq_len(G)){
    row<-row+1L
    common<-pi[d]+pi[s]+emit(o,obs[1])+emit(d,obs[2])+emit(s,obs[3])
    terms[row,]<-common+c(log(trio(geno[o,],geno[d,],geno[s,])),
      log(po(geno[o,],geno[d,])),log(po(geno[o,],geno[s,])),pi[o])
  }
  lse<-function(x){m<-max(x);if(is.infinite(m))m else m+log(sum(exp(x-m)))}
  setNames(apply(terms,2,lse),c('H_PP','H_D0','H_0S','H_U'))
}

test_that("structural zeros are explicit while positive probabilities retain finite logs", {
  expect_identical(PaPeR:::.log_probability_ln(0), -Inf)
  expect_equal(PaPeR:::.log_probability_ln(1e-300), log(1e-300))
  expect_error(PaPeR:::.log_probability_ln(-1), "nonnegative finite")
  expect_error(PaPeR:::.log_probability_ln(Inf), "nonnegative finite")

  expect_identical(PaPeR:::.log_likelihood_contrast_ln(-Inf, -Inf), NA_real_)
  expect_identical(PaPeR:::.log_likelihood_contrast_ln(-Inf, -3), -Inf)
  expect_identical(PaPeR:::.log_likelihood_contrast_ln(-3, -Inf), Inf)
  expect_identical(PaPeR:::.log_likelihood_contrast_ln(-Inf, c(-Inf, -3)),
                   c(NA_real_, -Inf))
  expect_equal(PaPeR:::.log_likelihood_contrast_ln(-3, -4), 1)
})

test_that("native relationship screening preserves defined contrasts and membership", {
  tab <- PaPeR:::build_locus_table(c(.5, .5), Err = 1e-12)
  G <- matrix(c(1L, 1L, 2L, 2L), nrow = 2L, byrow = TRUE)
  registry <- PaPeR:::.catalog_relationship_registry("simp")
  Delta <- do.call(rbind, registry$Delta[match(c("PO", "U"), registry$config_id)])
  screen <- function(focal, alternative, threshold)
    PaPeR:::screen_rel_pairs_cpp(list(tab), G, Delta,
      as.integer(focal), as.integer(alternative), threshold,
      rep(NA_integer_, 2L), -1L, 0L, -1L)
  positive_over_zero <- screen(1L, 0L, 0)
  expect_equal(nrow(positive_over_zero), 1L)
  expected<-PaPeR:::calc_pair_LL_cpp(list(tab),G,0L,1L)
  expect_equal(positive_over_zero$llr,unname(expected["U"]-expected["PO"]),tolerance=1e-12)
  finite <- screen(1L, 1L, -.1)
  expect_equal(nrow(finite), 1L)
  expect_identical(finite$llr, 0)
  # Positive error makes both comparisons defined, so -Inf admits both.
  expect_equal(nrow(screen(0L, 0L, -Inf)), 1L)
  expect_equal(nrow(screen(0L, 1L, -Inf)), 1L)
})

test_that("likelihood kernels preserve tiny positive evidence and structural zeros", {
  p <- c(1e-150, 1 - 1e-150)
  tab <- PaPeR:::build_locus_table(p, Err = 1e-200)
  G <- matrix(c(1L, 1L, 2L, 2L), nrow = 2L, byrow = TRUE,
              dimnames = list(c("rare", "common"), c("locus.1", "locus.2")))
  pair <- PaPeR:::calc_pair_LL_cpp(list(tab), G, 0L, 1L)
  exact_pair<-.numeric_positive_trio(tab,c(1L,3L,NA_integer_))
  expect_equal(unname(pair[c("PO","U")]),unname(exact_pair[c("H_PP","H_U")]),tolerance=1e-11)
  expect_false(any(is.nan(pair)))

  trio <- PaPeR:::calc_trio_LL_cpp(list(tab), G, 0L, 0L, 1L)
  expected <- .numeric_positive_trio(tab,c(1L,1L,3L))
  expect_equal(trio, expected, tolerance = 1e-11)
  expect_false(any(is.nan(trio)))
})

test_that("thousands of independent loci accumulate in natural-log space", {
  n_loci <- 2000L
  p <- c(1e-150, 1 - 1e-150)
  tables <- lapply(seq_len(n_loci), function(i) PaPeR:::build_locus_table(p, Err = 1e-200))
  names(tables) <- paste0("L", seq_len(n_loci))
  calls <- rbind(rep("1/1", n_loci), rep("2/2", n_loci))
  rownames(calls) <- c("rare", "common")
  colnames(calls) <- names(tables)
  pair <- .calc_pair_LL(calls, "rare", "common", tables, Complex = "simp")
  exact<-.numeric_positive_trio(tables[[1]],c(1L,3L,NA_integer_))
  expected <- n_loci * exact[["H_U"]]
  expect_true(is.finite(pair["U"]))
  expect_equal(unname(pair["U"]), expected, tolerance = 1e-8)
  expect_equal(unname(pair["PO"]),n_loci*exact[["H_PP"]],tolerance=1e-8)
  expect_false(any(is.nan(pair)))
})

test_that("invalid probability cells fail instead of being logged as NaN", {
  tab <- PaPeR:::build_locus_table(c(0.5, 0.5), Err = 1e-12)
  tab$schema_version <- "obsolete-dense-schema"
  G <- matrix(c(1L, 1L, 1L, 1L), nrow = 2L, byrow = TRUE)
  expect_error(PaPeR:::calc_pair_LL_cpp(list(tab), G, 0L, 1L),
               "compact-locus-v1")
})
