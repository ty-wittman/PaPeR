# Validation

This page summarizes the evidence behind PaPeR's defaults: what was tested,
on which data, and the known limitations. The studies were run by the author
during development, on simulated pedigrees.

## Scoring

- **Correct**: an assigned parent matches the true parent.
- **Wrong**: any other assigned parent, including an assignment where the
  true parent was never genotyped.
- **Missed**: an empty slot whose true parent was genotyped.

## The example data

`pap_example` ships with the package and is built by
[`data-raw/make_example_data.R`](../data-raw/make_example_data.R) (seed 2026).
It has 700 individuals over six birth years, 595 of them genotyped on 100
microhaplotypes with 3–8 alleles and on 150 SNPs, with dropout, miscalls and
2% missing calls. With default settings:

| Markers | Dams correct / wrong / missed | Sires correct / wrong / missed |
|---|---|---|
| 100 microhaplotypes (`Dropout = 0.005`, `Foreign = 5e-4`) | 387 / 0 / 0 | 400 / 0 / 0 |
| 150 SNPs (Sequoia error model, `Err = 0.005`) | 383 / 1 / 4 | 396 / 0 / 4 |

The microhaplotype run, with final LLRs, took 55 s on one core of a cloud
Linux machine (R 4.6.2, GCC 13, `-O2`); the SNP run took 16 s.

## PaPeR on microhaplotypes against Sequoia on SNPs

This was the main comparison: PaPeR on microhaplotypes against stock Sequoia
on SNPs, using the same pedigrees, the same number of loci and the same
missingness masks.

**Setup.**

- Sequoia 3.2.0 with `Module = "par"` and `Complex = "full"`; PaPeR with its
  public `assign_parents()` defaults.
- Pedigrees drawn from one simulated reservoir population with overlapping
  generations, at 600, 1,000 and 2,000 individuals.
- 100, 200 and 300 loci. Microhaplotype errors used Dropout = e and
  Foreign = 0.1e. SNPs were simulated with Sequoia's `SimGeno()` (minor allele
  frequency 0.3) and errors added with `MkGenoErrors()`.
- True error rates e of 1e-4, 0.005 and 0.01, each crossed with the same three
  assumed rates.
- Call missingness of 0, 5% and 20%, crossed with 0, 10% and 25% of
  individuals ungenotyped.

That gives 3 × 3 × 3 × 3 × 3 × 3 = 729 paired conditions. Each ran as one task
on one CPU of a SLURM cluster.

**Results.**

- In nearly every condition, PaPeR on microhaplotypes missed fewer
  recoverable parents and made fewer wrong assignments than Sequoia on SNPs.
  Wrong assignments were about zero at 200–300 loci. The advantage was
  largest at 100 loci and at 2,000 individuals.
- 26 of the 729 conditions favoured Sequoia on some count:
  - in about 9, PaPeR made one or two more wrong assignments, all to parents
    outside the recoverable set, while finding 16–138 more correct parents;
  - in about 17, at a true error of 1e-4 and mostly with 1,000 individuals and
    200 loci, PaPeR missed 2–6 more parents (at worst 6 of 814). This is the
    low-error limitation listed below.
- When the assumed error rate was too high, calibration lowered Dropout to
  2–6e-4, and PaPeR made no wrong assignments against Sequoia's 4–5.

**Spot checks** (correct / wrong / missed, and run time):

| Condition | Sequoia on SNPs | PaPeR on microhaplotypes | Time (Sequoia / PaPeR) |
|---|---|---|---|
| 600 individuals, assumed error = true error | 688 / 3 / 31 | 709 / 0 / 10 | 3 s / 20 s |
| 600 individuals, assumed error 50× too low | 657 / 0 / 62 | 707 / 0 / 12 | 2 s / 55 s |
| 2,000 individuals | 2,751 / 2 / 49 | 2,787 / 0 / 13 | 29 s / 102 s |

In the 50× too-low condition, calibration raised Dropout from 1e-4 to 0.0030
and then 0.0038; the first search alone found 682 correct parents. Peak memory
at 600 individuals was about 0.3 GiB.

## The same SNP data in both programs

Before error calibration was added, both programs ran on the SNP panels of
20 blind-test datasets, at the supplied error rates:

| | Correct | Wrong | Missed |
|---|---|---|---|
| Sequoia | 14,728 | 69 | 1,350 |
| PaPeR | 14,656 | 54 | 1,427 |
| PaPeR on the matching microhaplotypes | 15,637 | 5 | 464 |

On SNPs the two find about the same parents, and PaPeR makes fewer wrong
assignments. A later change to how evidence is gathered altered only one of
these results (one more correct dam).

## Blind test sets

Three rounds of blind tests were run with the truth sealed before any run and
released only after the results were frozen:

1. **Mating systems and error mismatches** (10 datasets).
2. **Ancestral against full evidence, and pedigree loops** (10 datasets).
3. **Error-rate policy** under matched, misstated and nonstandard error (55
   datasets). This round found the SNP error estimator's bias and its
   interval coverage, listed below.

## Error calibration

- **Dropout-only calibration** (the final policy): 25 multiallelic settings
  plus a foreign-only stress case, 27 calls and 54 searches. All converged,
  with no unavailable calculations. Against the first search, calibration
  added 170 correct and removed 39 wrong assignments.
- **Defaults with a 5× misstated rate**: on outbred datasets the calibrated
  results matched the results at the correct rate exactly.
- The calibration estimator for SNPs runs about 11% low, and its nominal 95%
  interval covered the true rate in 54% of cases. Trios share parents and the
  selection of confident trios is not modelled, so the intervals are
  approximate.

## Numerical checks

- Final LLRs on tiny pedigrees match brute-force enumeration over all latent
  genotypes (`tests/testthat/test-final-llr-enumeration.R`).
- Relationship likelihoods match CKMRsim's `make_matrix_X_l()` to 1e-12 for
  unrelated, parent–offspring, full-sibling, second-degree, half-avuncular and
  duplicate relationships (`tests/testthat/test-core.R`, run when CKMRsim is
  installed).
- A looped-pedigree test checks that incrementally updated, rebuilt,
  unpruned and move-by-move replayed evaluations agree
  (`tests/testthat/test-ancestral-cold-rebuild.R`). During development,
  rebuilding all evidence from scratch reproduced the cached values for 2,674
  questions to within 3.6e-12 (natural log), and pruning unobserved
  individuals changed 916 questions by at most 1.4e-12.


## Known limitations

- At very low true error (about 1e-4), PaPeR misses slightly more parents than
  Sequoia.
- The SNP error-rate estimator runs about 11% low.
- On weak SNP panels whose fitted rate is low, the error rate is not raised.
- `Foreign` is fixed at its supplied value, so a misstated `Foreign` is never
  corrected.
- Forest cuts approximate the correlations between ancestral lines outside
  each local question.
- Full evidence (`evidence = "full"`; includes descendant information) and the S37 and S45 hypothesis suites
  are experimental.
- Individuals with no genotype calls slow the search several-fold without
  adding assignments, and they are not filtered out automatically. Remove them
  before running PaPeR.


