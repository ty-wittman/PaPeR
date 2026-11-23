# PaPeR 0.1.0

First public release.

* `assign_parents()` assigns genotyped dams and sires from microhaplotype or
  SNP genotypes and life-history data, with exact local pedigree likelihoods,
  the seven-relationship (S7) rival suite and Sequoia-style `Tassign` and
  `Tfilter` thresholds.
* Error models: dropout/foreign for microhaplotypes (`Dropout`, `Foreign`) and
  Sequoia's model for SNPs (`control$error_model = "sequoia"`, `Err`).
* Error-rate calibration from confident trios (`control$error_calibration`,
  default `"raise_guarded"`).
* Final LLRs for the assigned parents (`control$CalcLLR`) and for any fixed
  pedigree (`calc_llr()`).
* Diagnostics: `MHStats()`, `find_duplicates()`, `get_maybe_rel()`,
  `calc_pair_LL()`, `calc_trio_LL()` and `score_pedigree()`.
* Simulated example data `pap_example` and the getting-started vignette
  (`vignette("PaPeR")`).
