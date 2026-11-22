# Test coverage

`Rscript tools/test.R` runs all 42 supported test files. Quick, standard and
extended currently select the same files; there are no slower oracle-only tiers.
Use `--jobs=1` for serial files and `--files=` for a focused selection.
Each scientific call uses one BLAS/OpenMP thread.

| Contract | Test files (under `tests/testthat/`) |
| --- | --- |
| Public exports and defaults | `test-public-api.R`, `test-engine-defaults.R`, `test-utilities.R` |
| Genotypes, missingness, emission models and numerical zeros | `test-genotypes.R`, `test-stats.R`, `test-missingness.R`, `test-errors.R`, `test-dropout-foreign.R`, `test-numerical-stability.R`, `test-likelihood-contrast-zero-contract.R` |
| Relationship likelihoods and catalogs | `test-core.R`, `test-relationship-catalog.R`, `test-relationship-universe.R`, `test-inference-context.R` |
| Candidate discovery, exclusion, typed ages, state and retention | `test-candidate-discovery.R`, `test-exclusion.R`, `test-pedigree-consistency.R`, `test-typed-age.R`, `test-search-initialization.R`, `test-retention-policy.R` |
| Ancestral boundaries, messages and cold rebuilds | `test-ancestral-evidence.R`, `test-ancestral-cold-rebuild.R`, `test-ancestral-survivor-roots.R`, `test-focal-role-cavity.R` |
| Pre-screen, support, pair rescue and hypothesis suites | `test-ancestral-prescreen-default.R`, `test-joint-support-gates.R`, `test-pair-rescue.R`, `test-hypothesis-suites.R`, `test-s7-reverse-po.R` |
| Calibration and fixed Foreign | `test-error-calibration-*.R`, `test-calibration-fixed-foreign.R`, `test-two-pass-error.R`, `test-review-*.R`, `test-confidence-dependency-records.R`, `test-confidence-no-raw.R` |
| Fixed-pedigree final LLRs and nonmutation | `test-fixed-pedigree-llr.R`, `test-final-llr-enumeration.R` |
| Duplicate detection and likely relatives | `test-duplicate-report.R`, `test-utilities.R` |

The enumeration test sums latent genotypes on tiny pedigrees and computes the
emissions independently inside the test. `helper-simulate.R` simulates the
pedigrees and genotypes used by the utility tests and by the example data
(`data-raw/make_example_data.R`).
