# The PaPeR model: genotypes, likelihoods, messages and errors

This is the scientific reference for what PaPeR computes. It covers the
genotype representation, the error models, Mendelian inheritance, the
pedigree messages, the local likelihood questions, the decision scores and
error calibration. Code locations are in [architecture.md](architecture.md).

PaPeR's search follows the heuristic parentage approach of Sequoia
(Huisman 2017, <https://doi.org/10.1111/1755-0998.12665>): rounds of greedy
assignment, parents tested against alternative relationships with the
`Tassign` and `Tfilter` thresholds, and Sequoia's SNP error model (section 2).
The likelihood calculations described here are PaPeR's own. The
[README](../README.md#relationship-to-sequoia) summarizes what is shared and
what differs.

Logs are natural logs internally. Public thresholds and LLRs are base 10,
so the internal threshold is `Tassign_ln = Tassign · ln 10`.

## 1. Genotypes and allele frequencies

Each locus *l* has *k* alleles (k = 2 for a SNP). An individual's true genotype
is an unordered pair *g = {a, b}*, with *a ≤ b*. A locus has *k(k+1)/2*
genotypes, indexed (1,1), (1,2), …, (1,k), (2,2), …. Calls are given as `"a/b"`,
and `"0"` means the locus is missing.

Allele frequencies *p* come from the called genotypes of all typed individuals
(`MHStats`). An allele that is never observed gets a pseudocount of 0.5:
*p_a = (n_a + 0.5·[n_a = 0]) / (2·n_typed + 0.5·n_zero)*, where *n_zero* is the
number of unobserved alleles.

Loci are treated as independent, with no linkage disequilibrium. Every
likelihood below is a product over loci, so its log is a sum over loci.

## 2. Observation (error) models

The emission *e_i(g) = P(observed call of i | true genotype g)* is applied to each
called individual. A missing call has emission 1 for every genotype, so it is
neutral. All error rates are **strictly positive and panel-wide**; zero rates
and per-locus rates are rejected.

**Sequoia model (SNPs, `Err = E`).** This matches the Sequoia model:

- *P(obs = g | g) = 1 − E*.
- A wrong call *h ≠ g* has probability *E·w(h)/Z(g)*, where *w(h)* depends on the
  class of *h*:
  - *w_Far = (E/2)²* if *h* shares no allele with *g*;
  - *w_Hom = E/2* if *h* is a homozygote sharing an allele;
  - *w_Het = E − (E/2)²* if *h* is a heterozygote sharing an allele.
- *Z(g)* normalizes the wrong calls of *g* to total *E*.

**Dropout/foreign model (microhaplotypes, `Dropout = D`, `Foreign = F`).** Two
mechanisms act in order on the true genotype:

1. *Dropout*: a true heterozygote *{a, b}* collapses to *{a, a}* or *{b, b}*,
   each with probability *D/2*. Homozygotes are unaffected.
2. *Foreign*: each of the two remaining copies independently is replaced, with
   probability *F*, by a draw from the foreign distribution *q*. By default
   *q = p*, so a draw can reproduce the original allele.

With *r(x→u) = (1 − F)[u = x] + F·q_u*, two precursor copies *{x, y}* give

- *R({x,y}→{u,v}) = r(x→u)r(y→v) + r(x→v)r(y→u)* for *u ≠ v*;
- *R({x,y}→{u,u}) = r(x→u)r(y→u)*.

The emission for a true heterozygote is then
*e({a,b}) = (1 − D)·R({a,b}) + (D/2)·R({a,a}) + (D/2)·R({b,b})*. For a true
homozygote it is *e({a,a}) = R({a,a})*.

## 3. Founders and Mendelian transmission

- **Founders.** An individual with no assigned parent has the Hardy–Weinberg
  prior *π({a,a}) = p_a²* and *π({a,b}) = 2·p_a·p_b*.
- **Transmission.** A parent *{a, b}* transmits allele *x* with probability
  *t(x) = ([x = a] + [x = b]) / 2*. A child *{u, v}* of dam genotype *d* and
  sire genotype *s* has probability *T({u,v} | d, s) = t_d(u)t_s(v) + t_d(v)t_s(u)*
  when *u ≠ v*, and *t_d(u)t_s(u)* when *u = v*.
- **Missing parent.** A missing (untyped, unassigned) parent transmits from the
  population: *t(x) = p_x*.

The likelihood of a pedigree at one locus sums, over all true genotypes, the
product of founder priors, transmission factors and emissions.

## 4. Genotype messages (the pedigree forest)

The pedigree is held as a graph of individuals and parental families. A family
is a (dam, sire) pair with its children. Each individual–family incidence
carries two genotype messages, computed by sum-product:

- **family → child:** *m(g_c) = Σ_{d,s} T(g_c | d, s) · μ_dam(d) · μ_sire(s)*;
- **parent → family:** *μ_x(g) = e_x(g) · (message from x's own parental
  family, or π(g) if x is a founder) · (messages from x's other families)*.

**Ancestral evidence, the default.** Every descendant-to-ancestor message is
set to 1: from a child up into its family, and from a family into its parents.
Information therefore flows only downward. The message arriving at an
individual summarizes the observations of its assigned ancestors only, never of
its descendants, mates or collateral relatives. The experimental full-evidence
mode keeps both directions.

**Loops and cuts.** A pedigree with loops (inbreeding, repeated matings) is not
a tree. A spanning forest is chosen deterministically from the current
pedigree:

- children are scanned in typed-ID (input row) order;
- each family is numbered by its first child;
- within a family, incidences come in the order dam, sire, children ascending.

Union-find admits an incidence unless it would close a cycle; the rest are cut.
Messages across a cut treat the two sides as independent, which is the one
approximation outside the local region. Cuts depend only on the current
pedigree, not on search history. Every accepted move rebuilds the forest and
invalidates the affected messages.

## 5. The local question

Each decision evaluates a configuration *c = (d, s)* for a focal individual *f*.
Either *d* or *s* may be absent.

**Explicit variables:**

- the genotypes of *f*, *d* and *s*;
- their assigned parents and grandparents (depth 2);
- any latent individuals a rival hypothesis needs, such as the unobserved
  shared parent of half-siblings.

**Factors:**

- emissions of the named, observed individuals;
- transmission factors among the explicit variables;
- founder priors;
- forest messages at the depth-2 boundary, which carry the deeper assigned
  ancestry.

**Excluded:** the focal's incoming assigned links (the question replaces them),
descendants and collateral relatives of the named individuals, and every other
candidate.

The question is contracted exactly, per locus, by variable elimination. Loops
inside the local region are therefore exact. Two simplifications are made, and
neither changes the likelihood:

- unobserved leaves with no evidence below them sum to one and are pruned
  recursively before planning;
- normalizers and structural zeros are kept.

Three special values have distinct meanings:

- **Structural zero.** An impossible configuration has likelihood 0, which is
  `-Inf` in logs.
- **Unavailable.** A calculation that exceeds the resource limits
  (64 MiB / 400M operations per program) is unavailable. It is never treated as
  zero evidence.
- **Range.** Factors whose range exceeds double precision are stored as logs.

## 6. Scores and decisions

**Identity.** The score of configuration *c* is

```
I(c) = ln P(E_c | c) − ln P(E_c | focal unrelated to d and s),     I(0,0) = 0
```

*E_c* is the evidence of section 5. The denominator keeps the chosen parents'
own ancestry and any relationship between them; only the focal's proposed
links are removed. Another candidate's data never enters *I(c)*. The chooser
takes the configuration with the highest identity. It must beat the best
configuration with a different parent identity by `margin` (base 10).

**Support.** A new parent must have positive support beyond `Tassign`:

```
dam support  = I(d,s) − I(0,s)
sire support = I(d,s) − I(d,0)
pair support = I(d,s) − max(I(0,0), I(d,0), I(0,s))
```

Each required support must exceed `Tassign_ln`. An unavailable comparison
never rejects on its own.

**Relationship gate (S7).** For each new parent, the parent-offspring
numerator is compared with every admitted nonparent rival on the same
evidence. The S7 rivals are:

- unrelated;
- full sibling;
- half sibling (2 orientations);
- grandparent (8);
- full avuncular (4);
- half avuncular (8);
- reverse parent-offspring (2), where the focal is the candidate's parent.

Each rival is an explicit topology with its latent individuals. The parent must
beat the best rival by more than `Tassign_ln`. Age, sex, role and cycle
constraints decide which rivals are admitted. S37 and S45 add compound
relationships and are opt-in experiments. In S37 the half-sibling plus
parent-offspring numerator requires both endpoints to have no assigned parents;
parent forms enter the parent numerator, never a nonparent rival.

**Screens before scoring:**

- **Opposite homozygotes.** Opposite-homozygote counts are bounded by `T`, and
  Mendelian-error counts in trios by `MaxMendelE`. Both are calibrated from
  the tables so that a true parent passes with retention 0.9999^(1/N).
- **Tfilter.** An unconditional pedigree-free parent-offspring versus unrelated
  LLR must exceed `Tfilter` (default −2).
- **Shortlist.** At most 50 candidates are kept per side.
- **Pre-screen.** Candidates are ranked by single-parent identity *I(d,0)* or
  *I(0,s)*. Up to `prescreen_K` (10) per side within `prescreen_margin_ln`
  (ln 100) below `Tassign_ln` are kept. Pair rescue then checks the winner
  against omitted pairs that share one parent.
- **Ages.** Birth years, or BY.min/BY.max intervals, must make every
  parent–child link feasible. Typed intervals propagate through assigned links.

## 7. Search

The search runs Sequoia-style rounds over the typed individuals in a
deterministic order:

- Each focal visit scores its candidate configurations and accepts the best
  one that passes every gate.
- Accepted moves update the pedigree immediately; later visits see them.
- Complete parent pairs are rechecked locally.
- The search stops at a fixed point (a round with no moves) or a detected
  cycle, with at most 41 rounds.

## 8. Error calibration (`raise_guarded`)

1. Run the search at the supplied rates.
2. Take the complete trios whose decision margins clear `Tassign` (and
   `margin`) by an extra one unit (base 10). These are the confident trios.
3. Fit the error rate to their Mendelian-inconsistency counts:
   - each called locus of trio *j* is inconsistent with a probability that
     depends on the rate;
   - the count distribution is computed exactly and truncated at the
     `MaxMendelE` ceiling that admitted the trio;
   - the composite log-likelihood is maximized on [1e-12, 0.25].
   For microhaplotypes there are two inconsistency categories, homozygous and
   heterozygous offspring, which identify Dropout and Foreign separately. Their
   probabilities are exact low-degree polynomials in (D, F).
4. Compute the 95% profile bounds: the rates within *χ²₁(0.95)/2* of the
   maximum.
5. Adjust and restart:
   - Raise the rate (SNP `Err`, or Dropout) to its lower bound when that bound
     is higher.
   - On pass 1 only, a supplied rate above twice its upper bound is lowered to
     the upper bound, if at least 70% of assigned full trios are confident.
   - Foreign is fitted and profiled as a nuisance parameter but never changed.
   - Each restart is a fresh search from the original inputs. There are at
     most four searches.

The intervals are approximate. Trios share parents, and confidence selection is
not modeled. In blind tests the SNP estimator ran about 11% low, with 54%
interval coverage.

## 9. Final LLRs

`CalcLLR = TRUE` scores the fixed final pedigree with the final calibrated
error rates and life history. `calc_llr()` scores a supplied pedigree at the
supplied rates, without running a search or calibration. Both use a fresh
pedigree-likelihood setup, `control$llr_evidence` (default `"ancestral"`, independent of the
search's `control$evidence`) and `control$hypotheses`.
Only assigned parents enter the candidate arrays; assignments never change.

For each assigned side, the score is the parent log likelihood minus the
largest nonparent alternative in the selected suite, keeping the other
assigned parent. Natural-log contrasts are converted once to base 10.
`LLRpair = min(LLRdam, LLRsire)` requires both parents and both side results;
it is a weakest-side summary, not a separate joint-pair hypothesis test.

`ancestral` uses the search's ancestral evidence boundary. `full` also uses
descendant evidence through forest messages. Forest cuts approximate external
correlations in loopy pedigrees. Program memory/work limits apply to each calculation. Unknown results stay `NA` with an
`unavailable` status; unassigned sides have `not_assigned` status. `OHdam`,
`OHsire`, and `MEpair` count observed incompatibilities independently of score
availability.

The file-by-file implementation path is in [architecture.md](architecture.md).
`pedigree_forest` defines deterministic cuts, `forest_messages` stores the
external evidence, `plan_compiler` builds each local contraction, and
`plan_evaluator` / `lane_kernels` evaluate it. This organization does not
change the question boundary, estimands or numerical limits.
