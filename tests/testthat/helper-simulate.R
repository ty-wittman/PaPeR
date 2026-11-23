# Simulated pedigrees and genotypes for the tests and the example data.
# Base R only. data-raw/make_example_data.R sources this file to build
# `pap_example`.
#
#   sim_pedigree()        overlapping-generation pedigree with birth years and sex
#   sim_allele_freqs()    per-locus allele frequencies (k alleles each)
#   drop_genotypes()      true genotypes: Hardy-Weinberg founders, Mendelian transmission
#   add_errors()          observed genotypes: dropout, miscalls and missing calls
#   as_public_genotypes() PaPeR's "a/b" character format
# -----------------------------------------------------------------------------
# 1. Pedigree with overlapping generations, birth years, and sex
# -----------------------------------------------------------------------------
# Cohort-based: found a base cohort, then each year sample mothers & fathers from
# the pool of living, mature, correctly-sexed adults. This yields realistic age
# gaps (so the age prior has signal) and half-sibs through both sexes.
sim_pedigree <- function(n_founders   = 200,   # size of base (year-0) cohort
                         n_years       = 8,     # number of reproductive years after founding
                         n_per_year    = 150,   # offspring produced each year
                         mat_age       = 2,     # minimum age at first reproduction
                         max_age       = 10,    # maximum age an adult can reproduce
                         sex_ratio     = 0.5,   # P(female)
                         poly_dam      = FALSE, # a dam can have multiple mates within a year?
                         poly_sire     = TRUE,  # a sire can have multiple mates within a year?
                         seed          = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # founders
  ped <- data.frame(
    id        = sprintf("F%04d", seq_len(n_founders)),
    dam       = NA_character_,
    sire      = NA_character_,
    sex       = ifelse(runif(n_founders) < sex_ratio, 1L, 2L),  # 1=F, 2=M
    birthyear = 0L,
    stringsAsFactors = FALSE
  )

  next_id <- 1L
  for (yr in seq_len(n_years)) {
    # adults alive & mature this year
    age   <- yr - ped$birthyear
    adult <- age >= mat_age & age <= max_age
    dams  <- ped$id[adult & ped$sex == 1L]
    sires <- ped$id[adult & ped$sex == 2L]
    if (length(dams) == 0 || length(sires) == 0) next

    n_off <- n_per_year
    # sample a mother for each offspring
    md <- sample(dams,  n_off, replace = poly_dam  || n_off > length(dams))
    ms <- sample(sires, n_off, replace = poly_sire || n_off > length(sires))

    kids <- data.frame(
      id        = sprintf("Y%d_%04d", yr, seq_len(n_off)),
      dam       = md,
      sire      = ms,
      sex       = ifelse(runif(n_off) < sex_ratio, 1L, 2L),
      birthyear = yr,
      stringsAsFactors = FALSE
    )
    ped <- rbind(ped, kids)
    next_id <- next_id + n_off
  }
  rownames(ped) <- NULL
  ped
}
# LifeHistData in Sequoia's format
as_lifehist <- function(ped) {
  data.frame(ID = ped$id, Sex = ped$sex, BirthYear = ped$birthyear,
             stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# 2. Per-locus allele frequencies
# -----------------------------------------------------------------------------
# n_alleles: scalar or length-n_loci vector giving alleles per locus.
# skew: Dirichlet-ish concentration. skew=1 -> even-ish; <1 -> a few common + rare
#       tail (realistic for microhaplotypes).
sim_allele_freqs <- function(n_loci = 200, n_alleles = 6, skew = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(n_alleles) == 1L) n_alleles <- rep(n_alleles, n_loci)
  stopifnot(length(n_alleles) == n_loci)
  lapply(n_alleles, function(k) {
    g <- rgamma(k, shape = skew, rate = 1)   # Dirichlet via normalized gammas
    g / sum(g)
  })
}


# -----------------------------------------------------------------------------
# 3. TRUE genotypes: HWE founders + Mendelian transmission
# -----------------------------------------------------------------------------
# Returns a simulator-private integer matrix: nrow = n_indiv, ncol = 2 * n_loci.
# Columns (2l-1, 2l) hold the two gene copies at locus l, coded 1..k_l.
drop_genotypes <- function(ped, freqs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_loci <- length(freqs)
  n_ind  <- nrow(ped)
  G <- matrix(0L, nrow = n_ind, ncol = 2L * n_loci)
  rownames(G) <- ped$id

  # process parents-before-offspring: birth year is a valid topological order
  ord    <- order(ped$birthyear)
  id2row <- setNames(seq_len(n_ind), ped$id)

  is_founder <- is.na(ped$dam) & is.na(ped$sire)

  for (i in ord) {
    cols_a <- 2L * seq_len(n_loci) - 1L   # first-copy columns
    cols_b <- 2L * seq_len(n_loci)        # second-copy columns
    if (is_founder[i]) {
      # sample both copies from HWE
      G[i, cols_a] <- vapply(freqs, function(p) sample.int(length(p), 1L, prob = p), 1L)
      G[i, cols_b] <- vapply(freqs, function(p) sample.int(length(p), 1L, prob = p), 1L)
    } else {
      dr <- id2row[[ped$dam[i]]]
      sr <- id2row[[ped$sire[i]]]
      # inherit one random copy from each parent, per locus
      pick_dam  <- ifelse(runif(n_loci) < 0.5, cols_a, cols_b)
      pick_sire <- ifelse(runif(n_loci) < 0.5, cols_a, cols_b)
      G[i, cols_a] <- G[dr, pick_dam]
      G[i, cols_b] <- G[sr, pick_sire]
    }
  }
  attr(G, "n_alleles") <- vapply(freqs, length, 1L)
  class(G) <- c("pap_simulator_internal", "matrix", "array")
  G
}


# -----------------------------------------------------------------------------
# 4. Genotyping error: allelic dropout, miscall, whole-locus missingness
# -----------------------------------------------------------------------------
# dropout : per-GENE-COPY dropout prob. One copy dropping in a het -> false
#           homozygote (the key microhap error mode). Both drop -> missing.
# epsilon : per-observed-copy substitution prob -> uniform random other allele.
# miss    : per-locus whole-genotype missingness.
add_errors <- function(G, dropout = 0.01, epsilon = 0.005, miss = 0.02, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k_all  <- attr(G, "n_alleles")
  n_loci <- length(k_all)
  n_ind  <- nrow(G)
  O <- G

  # dropout/epsilon/miss may be scalar (global) or length-n_loci (per-locus).
  dropout <- rep_len(dropout, n_loci)
  epsilon <- rep_len(epsilon, n_loci)
  miss    <- rep_len(miss,    n_loci)

  for (l in seq_len(n_loci)) {
    ca <- 2L * l - 1L; cb <- 2L * l
    a <- O[, ca]; b <- O[, cb]
    k <- k_all[l]

    # (i) whole-locus missing
    m <- runif(n_ind) < miss[l]
    # (ii) per-copy dropout
    da <- runif(n_ind) < dropout[l]
    db <- runif(n_ind) < dropout[l]
    both_drop <- da & db
    a_drop    <- da & !db          # copy a lost -> observed hom for b
    b_drop    <- db & !da          # copy b lost -> observed hom for a
    a[b_drop] <- a[b_drop]; b[b_drop] <- a[b_drop]
    b[a_drop] <- b[a_drop]; a[a_drop] <- b[a_drop]
    # (iii) substitution on surviving copies
    if (k > 1) {
      sa <- runif(n_ind) < epsilon[l]
      sb <- runif(n_ind) < epsilon[l]
      a[sa] <- misscall(a[sa], k)
      b[sb] <- misscall(b[sb], k)
    }
    # apply missing / both-dropped last
    gone <- m | both_drop
    a[gone] <- 0L; b[gone] <- 0L

    O[, ca] <- a; O[, cb] <- b
  }
  attr(O, "n_alleles") <- k_all
  class(O) <- c("pap_simulator_internal", "matrix", "array")
  O
}

# Convert the simulator-private two-copy matrix to PaPeR's public format.
as_public_genotypes <- function(G, locus_names = NULL) {
  if (!inherits(G, "pap_simulator_internal") || ncol(G) %% 2L != 0L)
    stop("G must be a simulator-private two-copy genotype matrix")
  n_loci <- ncol(G) / 2L
  if (is.null(locus_names)) locus_names <- sprintf("MH%03d", seq_len(n_loci))
  stopifnot(length(locus_names) == n_loci, !anyDuplicated(locus_names),
            all(nzchar(locus_names)))
  out <- matrix("0", nrow(G), n_loci,
                dimnames = list(rownames(G), locus_names))
  for (l in seq_len(n_loci)) {
    a <- G[, 2L * l - 1L]
    b <- G[, 2L * l]
    typed <- a > 0L & b > 0L
    out[typed, l] <- paste(pmin(a[typed], b[typed]),
                           pmax(a[typed], b[typed]), sep = "/")
  }
  attr(out, "n_alleles") <- attr(G, "n_alleles")
  out
}

# replace each allele with a uniformly random DIFFERENT allele in 1..k
misscall <- function(x, k) {
  if (length(x) == 0L) return(x)
  repl <- sample.int(k, length(x), replace = TRUE)
  bad  <- repl == x
  while (any(bad)) {                       # rejection to guarantee a change
    repl[bad] <- sample.int(k, sum(bad), replace = TRUE)
    bad <- repl == x
  }
  as.integer(repl)
}


# -----------------------------------------------------------------------------
# 5. Prune: mask individuals as "not genotyped" (for reconstruction tests)
# -----------------------------------------------------------------------------
# Keeps the FULL true pedigree (for scoring) but returns a GenoM containing only
# the retained individuals. `keep` is a logical or a vector of IDs to retain.
prune_samples <- function(G, keep) {
  if (is.logical(keep)) keep <- rownames(G)[keep]
  Gk <- G[rownames(G) %in% keep, , drop = FALSE]
  attr(Gk, "n_alleles") <- attr(G, "n_alleles")
  class(Gk) <- class(G)
  Gk
}


# -----------------------------------------------------------------------------
# 6. Collapse to Sequoia 0/1/2 (baseline: run stock Sequoia on best single allele)
# -----------------------------------------------------------------------------
# For each locus, reference = most frequent allele; code = count of reference
# copies (0/1/2), missing = -9. Lets you benchmark against real Sequoia.
microhap_to_snp <- function(G, freqs) {
  n_loci <- length(freqs)
  n_ind  <- nrow(G)
  S <- matrix(-9L, nrow = n_ind, ncol = n_loci,
              dimnames = list(rownames(G), sprintf("snp%03d", seq_len(n_loci))))
  for (l in seq_len(n_loci)) {
    ref <- which.max(freqs[[l]])
    a <- G[, 2L * l - 1L]; b <- G[, 2L * l]
    obs <- a != 0L & b != 0L
    S[obs, l] <- (a[obs] == ref) + (b[obs] == ref)
  }
  S
}


# -----------------------------------------------------------------------------
# 7. Convenience wrapper
# -----------------------------------------------------------------------------
simulate_dataset <- function(n_founders = 200, n_years = 8, n_per_year = 150,
                             n_loci = 200, n_alleles = 6, skew = 0.7,
                             dropout = 0.01, epsilon = 0.005, miss = 0.02,
                             seed = 1) {
  set.seed(seed)
  ped   <- sim_pedigree(n_founders, n_years, n_per_year)
  freqs <- sim_allele_freqs(n_loci, n_alleles, skew)
  Gtrue <- drop_genotypes(ped, freqs)
  Gobs  <- add_errors(Gtrue, dropout, epsilon, miss)
  list(pedigree     = ped,
       lifehist     = as_lifehist(ped),
       freqs        = freqs,
       geno_true    = as_public_genotypes(Gtrue),
       geno_obs     = as_public_genotypes(Gobs))
}
