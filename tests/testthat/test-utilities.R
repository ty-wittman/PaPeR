# Diagnostic utilities: duplicates, CalcOHLLR, and GetMaybeRel.

# Small on purpose: both utilities scan every pair of individuals.
make_data <- function(seed = 501, nloci = 60) {
  # Simulator functions come from helper-simulate.R.
  set.seed(seed)
  freqs <- lapply(sim_allele_freqs(nloci, sample(3:8, nloci, TRUE), 0.7, seed),
                  function(p) (p + 1e-3) / sum(p + 1e-3))
  ped <- sim_pedigree(30, 3, 40, seed = seed)
  Go <- as_public_genotypes(add_errors(drop_genotypes(ped, freqs, seed = seed + 1),
                                       0.008, 0.006, 0.01, seed = seed + 2))
  list(ped = ped, Go = Go, tables = build_tables(freqs, Err = 0.008))
}

test_that("find_duplicates detects an injected duplicate sample", {
  d <- make_data(511)
  dupid <- rownames(d$Go)[50]
  Gd <- rbind(d$Go, COPY = d$Go[dupid, ])
  dups <- .find_duplicates(Gd, d$tables, thresh = 10)
  expect_true(any((dups$ID1 == dupid & dups$ID2 == "COPY") |
                  (dups$ID2 == dupid & dups$ID1 == "COPY")))
})


test_that("get_maybe_rel flags a hidden parent-offspring link", {
  d <- make_data(513)
  ph <- d$ped[, c("id", "dam", "sire")]
  hid <- which(!is.na(ph$dam))[1]; truedam <- ph$dam[hid]; foc <- ph$id[hid]
  ph$dam[hid] <- NA
  mr <- .get_maybe_rel(d$Go, d$tables, Pedigree = ph, thresh = 5)
  expect_true(any((mr$ID1 == foc & mr$ID2 == truedam) |
                  (mr$ID2 == foc & mr$ID1 == truedam)))
})

