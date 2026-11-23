# Builds data/pap_example.rda, the simulated example dataset.
# Run from the repository root:  Rscript data-raw/make_example_data.R
#
# A population with overlapping generations: 150 founders born in 2015, then
# 110 offspring a year for five years (2016-2020). Parents breed from age 1
# to 4; sires can have several mates a year. 15% of individuals are not
# genotyped, so some true parents cannot be found. The same individuals are
# genotyped twice: on 100 microhaplotypes (3-8 alleles) and on 150 SNPs.

source("tests/testthat/helper-simulate.R")

set.seed(2026)
truth <- sim_pedigree(n_founders = 150, n_years = 5, n_per_year = 110,
                      mat_age = 1, max_age = 4)
truth$birthyear <- truth$birthyear + 2015L

# Readable IDs in birth order.
new_id <- setNames(sprintf("ID%04d", seq_len(nrow(truth))), truth$id)
truth$id   <- unname(new_id[truth$id])
truth$dam  <- unname(new_id[truth$dam])
truth$sire <- unname(new_id[truth$sire])

# Microhaplotypes: 100 loci with 3-8 alleles and skewed frequencies.
mh_freqs <- lapply(sim_allele_freqs(100, sample(3:8, 100, TRUE), skew = 0.7),
                   function(p) (p + 1e-3) / sum(p + 1e-3))
mh_true <- drop_genotypes(truth, mh_freqs)
# Per-copy dropout 0.003 (a heterozygote reads as a homozygote about 0.6% of
# the time), per-copy miscall 5e-4, 2% missing calls.
mh_obs  <- add_errors(mh_true, dropout = 0.003, epsilon = 5e-4, miss = 0.02)
GenoM <- as_public_genotypes(mh_obs, sprintf("MH%03d", seq_len(100)))

# SNPs: 150 biallelic loci, minor allele frequency 0.1-0.5, per-copy miscall
# 0.0025 (about 0.5% of genotypes wrong), 2% missing calls.
snp_freqs <- lapply(runif(150, 0.1, 0.5), function(q) c(1 - q, q))
snp_true <- drop_genotypes(truth, snp_freqs)
snp_obs  <- add_errors(snp_true, dropout = 0, epsilon = 0.0025, miss = 0.02)
GenoM_snp <- as_public_genotypes(snp_obs, sprintf("SNP%03d", seq_len(150)))

# 15% of individuals are never genotyped.
genotyped <- sort(sample(truth$id, round(0.85 * nrow(truth))))
GenoM     <- GenoM[genotyped, ]
GenoM_snp <- GenoM_snp[genotyped, ]
attr(GenoM, "n_alleles") <- NULL
attr(GenoM_snp, "n_alleles") <- NULL

LifeHistData <- data.frame(ID = genotyped,
                           Sex = truth$sex[match(genotyped, truth$id)],
                           BirthYear = truth$birthyear[match(genotyped, truth$id)])

pap_example <- list(
  GenoM = GenoM,
  GenoM_snp = GenoM_snp,
  LifeHistData = LifeHistData,
  Pedigree = data.frame(id = truth$id, dam = truth$dam, sire = truth$sire,
                        Sex = truth$sex, BirthYear = truth$birthyear,
                        genotyped = truth$id %in% genotyped))

save(pap_example, file = "data/pap_example.rda", compress = "xz")
