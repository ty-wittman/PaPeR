# Architecture

PaPeR assigns **genotyped parents to genotyped individuals**. The only
pedigree state is one dam slot and one sire slot per genotyped individual.
Sibship clustering, untyped dummy parents and inferred shared mates are
outside its scope.

The statistical model (genotype likelihoods, messages, scores, the error
models and calibration) is in [model.md](model.md), and results are in
[validation.md](validation.md). This page covers the interface, the order of
work in a run, the design rules the code keeps, and where each part lives.

## Interface

```r
assign_parents(GenoM, LifeHistData, Dropout, Foreign, control = list())
```

- `GenoM`: character calls `"a/b"`, `"0"` missing, one column per locus.
- `Dropout`, `Foreign`: one positive rate each for the whole panel
  (dropout/foreign model). SNPs can instead use Sequoia's model with
  `control = list(error_model = "sequoia", Err = e)`.
- `control`: every other setting. `pap_control()` returns the full list with
  defaults and rejects unknown names or invalid values.

The utilities `find_duplicates()`, `get_maybe_rel()`, `calc_pair_LL()` and
`calc_trio_LL()` take the same `GenoM, Dropout, Foreign, control` arguments;
`calc_llr()` scores a fixed pedigree with the same engine as the search.

Defaults: ancestral evidence, the seven-relationship (S7) hypothesis suite,
the candidate pre-screen (`prescreen_K = 10`, `prescreen_margin_ln =
log(100)`), `Tassign = 0.5`, `Tfilter = -2` and
`error_calibration = "raise_guarded"`. `evidence = "full"` and the S37 and S45
suites are experimental.

## What one call does

1. **Tables.** Allele frequencies come from `MHStats(GenoM)` unless
   `control$frequencies` is given. Each locus gets a compact emission table
   for the chosen error model.
2. **Bounds.** The life history gives birth-year intervals and parent age
   limits. Opposite-homozygote and Mendelian-error ceilings are calibrated
   from the tables, so that a true parent passes with probability `retain`.
3. **Candidates.** Each individual gets a Sequoia-style list of up to 50
   candidates of compatible sex and age that pass the ceilings and `Tfilter`.
4. **Search** (native code). Rounds visit every individual. A visit builds
   the local question for each candidate configuration: the focal, the
   proposed parents and their assigned ancestry, with deeper ancestry from
   cached forest messages. It screens the question, scores identity, support
   and the rival relationships, and accepts the best configuration that
   passes every gate. Accepted moves update the pedigree, the age intervals
   and the messages. Complete parent pairs are rechecked when their
   neighbourhood changes. The search stops at a fixed point or after
   `max_rounds`.
5. **Calibration.** The error rate is fitted on the confident trios of the
   final pedigree. Under `raise_guarded` it may be raised (or once lowered),
   and the search restarts from the original inputs, at most four searches
   in all. Foreign is never changed.
6. **Output.** `id`, `dam`, `sire`, plus final LLRs when `CalcLLR = TRUE`,
   and the attributes `calibration` (one row per search),
   `calibration_passes`, `converged`, `rounds` and `settings`.

## Design rules

- **Evidence scope.** A question contains the focal individual, the proposed
  parents and their assigned parents and grandparents (explicit to depth 2).
  Deeper ancestry enters only through forest messages, and in ancestral mode
  no descendant-to-ancestor message is used. Other candidates never enter a
  question, so adding an unrelated candidate cannot change a score.
- **Deterministic forest cuts.** Cuts depend only on the current pedigree:
  children in input order, families numbered by their first child, dam before
  sire before children. Union-find admits an incidence unless it closes a
  cycle. Small loops inside a question are contracted exactly; cuts
  approximate the correlations between ancestral lines outside it.
- **Barren pruning.** Unobserved leaves and their transmission factors are
  removed before planning; priors, normalizers and structural zeros are kept.
- **Numbers.** Natural logs internally, base-10 LLRs in the output. A
  structural zero is `-Inf`. A calculation above its resource limit (64 MiB
  and 400 million operations per program by default) is *unavailable*, never
  treated as zero evidence.
- **Error rates** are strictly positive and panel-wide.
- **Final LLRs** are computed on a fresh copy of the final pedigree and never
  change assignments.
- The search is single-threaded: one call uses one CPU core.

## Call paths

`assign_parents()` → `.pap_search()` → `.prepare_search()` → `.run_search()`
→ `run_parentage_cpp()` → `search_rounds_cpp()` → `run_visits()` →
configuration scoring → pedigree-conditioned genotype contractions.

`calc_llr()` → `.pap_final_llr()` → a fresh fixed-pedigree state and messages
→ `.prepare_search()` → `run_parentage_cpp()` → frozen `run_visits()` → the
same scorer and kernels, with only the assigned parents as candidates.
`assign_parents(CalcLLR = TRUE)` calls this path after calibration.

## R code (`R/`)

| File | Role |
| --- | --- |
| `pap.R` | Public settings, the `assign_parents()` wrapper, table preparation and search orchestration. |
| `genotypes.R` | Parse calls, missingness and allele encoding. |
| `stats.R` | `MHStats()`: allele frequencies and locus statistics. |
| `errors.R` | Dropout/foreign and Sequoia emission probabilities. |
| `tables.R` | Per-locus genotype, transmission and emission tables; `calc_pair_LL()`, `calc_trio_LL()`. |
| `lifehistory.R` | Validate sex and birth-year bounds. |
| `ageprior.R` | Parent age limits and admissible age intervals. |
| `exclusion.R` | Opposite-homozygote and Mendelian-error thresholds. |
| `candidate_pool.R`, `shortlist.R` | Discover compatible candidates and form the per-side lists. |
| `pedigree_state.R` | The authoritative parent slots, age state and graph updates. |
| `inference_context.R` | Genotype tables and observations shared by a call. |
| `relationship_catalog.R`, `relationships.R` | Relationship codes, orientations and question sets; pair and trio utilities. |
| `messages.R` | Initial forest messages, revisions and compiled snapshots. |
| `visit_inputs.R`, `search_engine.R`, `engine_controls.R`, `zz_scorer.R` | Prepare native inputs, call the engine and collect parents, moves and histories. |
| `error_calibration*.R` | Trio confidence, count likelihoods, rate fits, raise/lower decisions and restarts. |
| `final_llr.R` | `calc_llr()`: read-only final-pedigree LLRs, statuses and incompatibility counts. |
| `utilities.R` | `find_duplicates()` and `get_maybe_rel()`. |
| `score_pedigree.R` | Compare a pedigree with the truth. |
| `progress.R`, `log_units.R` | Progress display; natural-log and base-10 conversions. |

## Native code (`src/`)

The C++20 engine is organized in layers. Each `.cpp` compiles as its own
translation unit with a small declaration header; the remaining headers hold
shared value types (scores, topologies, evidence views and log arithmetic).

| Layer | Files | Role |
| --- | --- | --- |
| Driver | `driver.cpp`, `engine_state.*`, `search_inputs.*`, `execution_config.*` | Registered entry points; compile the per-call snapshot; validate inputs; orchestrate search rounds and collect results. |
| Visits | `visits.*`, `visit_graph.*`, `scoring.*`, `decisions.*`, `recheck.*`, `message_updates.*`, `graph_traversal.*`, `result_columns.*` | Visit focal individuals, score candidate configurations, choose and accept moves, recheck stale parent pairs, publish message updates and record results. |
| Scorer | `scorer.*`, `candidate_catalog.h`, `hypothesis_suite.h`, `relationship_topology.h`, `relationship_catalog.h`, `joint_identity.h`, `candidate_identity.h`, `predictive_models.h`, `partition_totals.h`, `panel_rejection.h`, `support_rejection.h`, `joint_support.h`, `pair_rescue.h` | Build each relationship question as an explicit local pedigree, reduce marker scores to identity, support and rival comparisons, and certify rejections. Supports are evaluated first, so a configuration that fails support is never scored against rivals. |
| Pedigree messages | `pedigree_forest.*`, `forest_messages.*` | Deterministic spanning forest and lazily updated sum-product genotype messages. |
| Numerical core | `plan_compiler.*`, `plan_evaluator.*`, `lane_kernels.*`, `likelihood_core.cpp`, `locus_tables.*`, `locus_workspace.*`, `resource_*.h`, `program_*.h` | Compile each question into a variable-elimination program and evaluate it per locus with vectorized kernels, within explicit memory and work limits. |
| Screens and life history | `exclusion.cpp`, `genotype_checks.*`, `typed_age.*`, `relative_age_cache.h` | Opposite-homozygote and trio checks, and birth-year intervals propagated through assigned links. |

## Tests

`tests/testthat/` has 42 test files. They cover the public interface,
genotype parsing and missingness, both error models, relationship likelihoods,
candidate discovery, age handling, messages and cold rebuilds, the decision
gates, calibration, and final LLRs. Exact final LLRs on tiny pedigrees are
checked against brute-force genotype enumeration
(`test-final-llr-enumeration.R`). Relationship likelihoods are optionally
cross-checked against CKMRsim.

```sh
R CMD INSTALL .
Rscript tools/test.R            # all test files, in parallel, about 30 seconds
```

`R CMD check` runs the same tests through `tests/testthat.R`.
