# Consolidated results — SSW surface-impact methodology

**Compiled 2026-08-03. This is the single source of truth for the project.**

It supersedes `paper/ng_manuscript.tex` (last modified 2026-07-28), which predates
every result below, contains none of them, and still leads with the snow-avalanche
application that the project abandoned. **Do not write from the manuscript. Write
from this file.**

Every number here is traced to a script and a JSON. Where a claim was withdrawn or
corrected, that is recorded in §7 rather than quietly dropped — nine headline
claims died during this project, and the tests that killed them are the reason the
survivors are trustworthy.

---

## 1. The claim, in one paragraph

The literature classifies sudden stratospheric warmings into "downward-propagating"
and not, and reports the between-group surface contrast as an effect size. That
contrast is mostly manufactured by the classification itself: the classifier window
overlaps the response window, so splitting events on it splits one distribution on
its own tail — 93–102% of the published contrast is reproduced by applying the same
criterion to dates with no SSW in them. Measured out of sample within a single
climate, **an SSW adds essentially nothing to the predictability of the surface
response** (+0.004 [−0.022, +0.027] pre-onset, +0.022 [−0.015, +0.060] post-onset,
against a size-matched null of 1,000 event-free draws), while the splits in use implicitly
claim up to **0.63**. SSWs do have a large, real surface effect — measured causally
in a nudged experiment at **+1.11 σ** across 9 models — and the all-event composite stands. What
carries little information is the sub-classification.

*(An earlier version of this paragraph put the predictable share at R² ≈ 0.11. That
figure was pooled across ensemble members and uncontrolled for what the same
pipeline achieves with no event present; it is withdrawn in §3 J-REVISED, and the
withdrawal now carries intervals.)*

---

## 2. The result chain

| # | result | system | status |
|---|---|---|---|
| A | design choices (estimator, catalogue) do **not** reshape estimates | obs, n=43 | established |
| B | the SSW→AO response is robust: −1.03 σ at +15..+29 d | obs, n=43 | established |
| C | R = 0.82 measures temporal **concentration**, not causal fraction | obs + 846 model events | **corrected — see §3C** |
| C2 | Granger vortex→surface asymmetry | obs n.s. (P=0.213); 4 models sig. | **not established in obs** |
| D | mediator measurement error under-attributes the stratosphere | obs, 6,958 days | established |
| E | selection on the outcome manufactures 30–44% of a reported impact | obs | established |
| F | the bias obeys a projection law, β × ΔS | obs R²=0.9945, CMIP6 R²=1.0000 | established |
| G | stratospheric classification is **not** exempt from it | obs + CMIP6 | established |
| H | the response distribution is one shifted population | CMIP6 n=1888 | established, power-bounded |
| I | forced between-event sd σ_f ≤ 0.26 σ (upper 95%) | CMIP6 n=1888 | established, 84% placebo correction |
| J | **out-of-sample predictable share R² = 0.115; 73% of diagnostic skill is window overlap** | CMIP6 n=1888 | **established, no correction** |
| K | causal stratosphere→surface effect S = +1.11 σ (sd 0.271; +1.36 incl. ECCC) | SNAPSI, 9 centres | established, 2 NH events; ECCC a Tukey outlier |
| **L** | **the DW/NDW contrast is identical with and without an SSW, in a designed experiment** | **SNAPSI, 9 centres, 61 ensembles** | **established — the decisive test, replicated at 3x scale 2026-09-17** |
| **M** | **an SSW shifts the surface distribution without inflating its variance** | **SNAPSI, 9 centres, 1798 vs 1805 members** | **established causally; ratio 0.955 [0.873, 1.051] after an estimator fix** |
| O | the field's two archetypes (Feb-2018 "propagating", Jan-2019 "not") have the same forced odds, and their observed difference is ordinary member-to-member noise | SNAPSI 9 centres + ERA5 | established for these 2 events; Jan-2019's NDW label flips with the pressure level |
| N | the week-2 100 hPa → surface coupling between members is the same with and without an SSW | SNAPSI, 9 centres, 36 matched ensembles | established; Δz +0.023 [−0.044, +0.090] |

**J is the headline.** It needs no correction of any kind, it answers the field's
own question in the field's own terms, and it independently confirms I.

---

## 3. Headline numbers

### A/B — robustness (`FINDING_design_robustness.md`)
- estimator choice (conventional vs within-winter): **0.002 σ**
- catalogue choice, common 1980–2019 window: **0.115 σ** (0.426 → 0.141 once the
  window is fixed; two-thirds of the apparent spread was record length, not event
  definition)
- post-onset AO, +15..+29 d: **−1.026 [−1.612, −0.422], p = 0.0013**, n=43
- robust to subsetting: isolated −1.226, satellite-era −0.934, both −1.020; every
  interval excludes zero
- pre-onset −30..−16 d: −0.820 (p=0.015), **fragile** — loses significance under
  every restriction, but every restricted interval still contains −0.820

### C — the R diagnostic: CORRECTED (`R_profile.py`, `R_gap_sensitivity.json`)
`R = within-winter effect / between-winter effect`. It was read as a causal
fraction: R≈1 ⇒ event-scale (causal downward influence), R≈0 ⇒ winter-scale
(common cause). **That reading is wrong, and this is now established at model
scale rather than suspected at n=43.**

R is **flat across lag** in both systems:

| lag bin | observations (n=43) | CanESM5 (10 members) |
|---|---|---|
| −30..−16 (**pre**-onset) | **+0.789** [0.34, 1.16] | **+0.754** [0.65, 0.84] (n=794) |
| −15..−1 | +0.714 | +0.744 |
| +0..+14 | +0.784 | +0.820 |
| +15..+29 (**post**-onset) | **+0.824** [0.49, 1.09] | **+0.825** [0.76, 0.89] (n=846) |

A precursor cannot be caused by the event that follows it, so R at a pre-onset
lag cannot be a causal fraction — yet it equals the post-onset value in both
systems, and in CanESM5 the intervals are ±0.10 wide. Flatness also survives
every baseline gap tested, 60/75/90 d (`R_gap_sensitivity.json`: "FLAT AT EVERY
GAP").

**Mechanism.** `within` subtracts a baseline that excludes |lag| ≤ 75 d, so *any*
anomaly living inside that window — before or after onset — is absent from the
baseline and returns R ≈ 1. **R measures temporal concentration of the anomaly
around the onset date.**

**Why the validation missed it — and the gap is now CLOSED.**
`validate_R_diagnostic.py` originally simulated only two structures: whole-winter
depression (f=0) and post-onset-only depression (f=1, `y[p:p+61] -= c`). It never
simulated a precursor, so it could not discover that a pre-onset anomaly also
scores R ≈ 1. A **precursor arm was added 2026-08-04** — the anomaly placed at
days −60..0, where its causal fraction is **zero by construction**:

| true f | R, anomaly AFTER onset | R, anomaly BEFORE onset |
|---|---|---|
| 0.00 | −0.032 | −0.062 |
| 0.25 | +0.203 | +0.214 |
| 0.50 | +0.481 | +0.479 |
| 0.75 | +0.703 | +0.749 |
| **1.00** | **+0.959** | **+1.012** |

**R returns the same value for an anomaly that cannot possibly have been caused by
the event.** Difference at f=1: +0.053. This settles it in simulation, under known
truth, and independently of the observational and CanESM5 evidence: *R is not an
estimator of the causal fraction.* What the grid parameter actually indexes is how
tightly the anomaly is concentrated around the onset date, not how much of it the
event caused.

> A first version of the precursor arm left the measurement window at the
> post-onset (15, 29) for **both** arms, so the pre-onset arm measured a window
> containing no anomaly: the between-winter denominator went to zero and R became
> a garbage ratio (sd = 26 at f=1), which looked like "R separates pre from post."
> The window must follow the anomaly. Fixed and recorded in the script.

**What survives:** the association is temporally **localised** around onset rather
than winter-scale. A short-timescale shared driver produces the same signature, so
this is not evidence of downward causation. For that, see K (SNAPSI), which is a
designed experiment.

> **Bug fixed 2026-08-03.** `R_profile.py` previously concatenated all 10 CanESM5
> members into one date-keyed frame and selected each event window by date alone —
> pulling rows from all ten members when only one had the SSW. That diluted the
> composite 3.5× (between −0.203 vs the correct −0.716 at +15..+29) and made the
> model look unlike observations. It now pools *pairs* per member, as
> `multimodel_R.py` always did; the two scripts now agree to 4 decimal places.
> The flatness conclusion is unaffected — it was always computed on the
> single-realisation observational path.

### C2 — Granger direction (`granger_direction.json`)
Vortex→surface minus surface→vortex, 10 lags, **2,000** winter-block bootstrap
replicates (re-run 2026-09-17; see the note below the table).

| dataset | spec | asymmetry | 95% CI | P(asym ≤ 0) |
|---|---|---|---|---|
| **observations** | conditional on v′T′ | **+0.00135** | [−0.0020, +0.0048] | **0.208** |
| **observations** | bivariate | **+0.00122** | [−0.0024, +0.0051] | **0.278** |
| CanESM5 | bivariate | +0.00335 | [+0.0021, +0.0047] | 0.000 |
| IPSL-CM5A2-INCA | bivariate | +0.00315 | [+0.0012, +0.0054] | 0.002 |
| CESM2-FV2 | bivariate | +0.00489 | [+0.0022, +0.0076] | 0.001 |
| MPI-ESM-1-2-HAM | bivariate | +0.00151 | [−0.0004, +0.0034] | 0.066 |

> **Re-run at 2,000 replicates on 2026-09-17, and nothing moved.** This was the
> only live result below the project's own N_BOOT ≥ 1000 floor — it shipped at
> 800 — and MPI-ESM-1-2-HAM at P = 0.063 sat close enough to 0.05 that the
> shortfall mattered. Point estimates are identical (the asymmetry is not
> resampled); every interval shifts in the fourth decimal; the p-values move
> 0.213→0.208, 0.296→0.278, 0.001→0.002, 0.003→0.001, and **0.063→0.066, still
> above 0.05**. No conclusion in this section changes. The earlier 800-replicate
> figures are superseded by the values in the table.

**Not significant in observations under either specification.** The observed CI
contains every model value, so this is consistent-but-underpowered at 47 winters,
not a null result. Effect sizes are tiny throughout — 0.1–0.5% of variance.

**Caveat the brief omits: the comparison is not matched.** Only observations were
run conditionally on v′T′; every model is bivariate. Conditioning is what removes
the shared wave-driving pathway, so the models' significance is measured against a
laxer specification than the observations'. Re-run the models conditionally before
this table is used to argue models and observations agree.

**This result was absent from earlier versions of this document.**

### D — mediation (`FINDING_mediation.md`)
Mediated share of the wave-driving→surface link, varying only the mediator:
single-day u10 **14.2%** → 10-day mean **45.2% [8, 88]** → disattenuated
**53.9% [−11, 140]**. Reliability of a 10-day u10 mean as a measure of the latent
vortex state: **0.592**, i.e. 41% error variance. Correlated errors bias toward
the naive estimate, so 53.9% is a **lower bound**.

### E/F/G — selection and the projection law (`FINDING_selection_on_outcome.md`, `FINDING_stratifier_law.md`)
- rate-matched selection bias under a true null: **AO 37%, NAO 30%, PNA 44%** of
  the reported dSSW impact
- decomposition closes numerically: residuals **0.010–0.024 σ** across three
  independent outcomes
- **the classification rate is uncontaminated**: 70% of real events vs 38% of
  pseudo-events satisfy the criterion, **P < 0.0001** — this is the quantity the
  field should report
- projection law `bias(Y) = β(Y-window on A-window) × bias(A)`, β the OLS **slope**:
  obs **R² = 0.9945** (7 stratifiers), CMIP6 **R² = 1.0000** (1888 events)
- **stratospheric classification is not exempt**: splitting random winter dates
  with no SSW by lower-stratospheric descent manufactures **+0.889 σ**; real SSWs
  split identically give +0.909 σ. The event adds 0.020 σ.
- **selection is not always inflationary** — for western-US snowpack the bias runs
  *opposite* to the effect. "Selection can only inflate, so this is a lower bound"
  is unsafe.
- body count on the published criterion (43 events, 30 DW / 13 NDW = 70%):

| outcome | contrast | share from selection | assumption-free null reproduces |
|---|---|---|---|
| AO | −1.782 | 99% | **96%** |
| NAO | −0.543 | 85% | **88%** |
| PNA | +0.222 | 76% | 67% |

- ERA5, in the literature's own fields: Karpechko 1000 hPa **96%** from selection
  (residual −0.026, p=0.589); ACP-2026 850 hPa **102%** (residual +0.013, p=0.906)

### H — one shifted population (`FINDING_no_class.md`)

| test | observations (n=43) | CMIP6 (n=1888) |
|---|---|---|
| dip, pseudo-calibrated | p = 0.964 | p = 0.888 |
| BIC k=1 vs k=2 | p = 0.140 | p = 0.835 |
| pure-shift KS | p = 0.326 | p = 0.279 |
| implied shift | −0.892 σ | −0.514 σ |

Power at n=1888, mixing fraction 2/3, total variance held fixed: **0.00** at 0.25σ
and 0.50σ, 0.04 at 0.75σ, **0.95 at 1.00σ**, 1.00 at 1.50σ. Nothing below ~1σ
is excluded and the claim must always carry that bound.

> **Re-run 2026-09-25 in the pinned environment:** the CMIP6 column moved (dip
> 0.846 → 0.888, BIC 0.865 → 0.835, KS 0.191 → 0.279, shift −0.511 → −0.514; power
> at 0.75σ 0.01 → 0.04, at 1.0σ 0.86 → 0.95); the observational column did not.
> Deterministic now (two runs identical); same pre-pin environment signature as
> I and J. The verdict is unchanged.

**"About two thirds propagate downward" is what one shifted population produces:**
observed pass rate on the criterion's surface conditions (1–2) 69.2% (53.8% with
all three conditions; 39 events in ERA5 coverage); unshifted pseudo-events 30.9%; pseudo-events displaced by
the measured −0.639 σ **74.4% [61.5, 87.2]**.

### I — forced variance (`FINDING_forced_variance_ceiling.md`)
`σ_f² = Var(Y|SSW) − Var(Y|pseudo)`.

| | CMIP6 (n=1888) | observations (n=43) |
|---|---|---|
| naive σ_f² | +0.0709 [+0.0354, +0.1074] | −0.0620 [−0.4332, +0.3136] |
| **placebo σ_f²** (pre-onset, true value 0) | **+0.0655** | +0.2832 |
| corrected σ_f² | **+0.0129 [−0.0401, +0.0656]** | −0.3453 |
| σ_f upper 95% | 0.256 (corrected) / 0.340 (raw) | 0.560 |

**84% of the raw excess variance is present before onset** — SSWs cluster in
disturbed winters whose variance is already elevated (pre-onset variance ratio
1.13). Without the placebo this was a clean false positive.

> **Re-run 2026-09-25 in the pinned environment (for the ACP constant) moved the
> CMIP6 column:** corrected σ_f² +0.0066 → +0.0129, upper 95% σ_f 0.239 → 0.256,
> placebo share 91% → 84%. Same 1,888 events; the pseudo-event set differs
> (37,747 → 37,734), the signature already seen in J (11,323 → 11,321) and
> attributed to package versions in the unpinned pre-2026-09-17 environment,
> which were never recorded. The new result is deterministic (two runs
> identical). The old one cannot be re-run, so it is superseded, not reconciled.
> The conclusion — σ_f² indistinguishable from zero, bounded near 0.25 σ — holds.

### J-REVISED (2026-08-03, after stress-testing) — READ THIS BEFORE THE TABLE BELOW

The figures in the original J are **pooled across ensemble members and uncontrolled
for what the same pipeline achieves with no event present**. Both were wrong to
report as "genuine predictability". `headline_stress_test.py` and
`within_model_check.py` correct them.

| tier | pooled (orig.) | within-model | pseudo-events, within-model | **event-specific** |
|---|---|---|---|---|
| P1 pre-onset | +0.0998 | +0.0664 | +0.0646 | **+0.0017** |
| P2 at-onset | +0.1121 | +0.0746 | +0.0767 | **−0.0021** |
| P3 + post-onset strat. | +0.4245 | +0.4013 | +0.3894 | **+0.0119** |

**Two corrections, both against the earlier claim.**
1. **A third of the pooled 0.115 was BETWEEN-MODEL signal.** GroupKFold by member
   prevents memorising a member but not exploiting cross-model structure — models
   with stronger vortex anomalies also have stronger surface responses. The
   within-climate value is 0.066, and the field's question is a within-climate one.
2. **Essentially all of what remains is not SSW-specific.** Applying the identical
   pipeline to pseudo-onsets on event-cleaned days gives 0.0646 against 0.0664 at
   real events. Against a size-matched null (below) **the event contributes
   +0.004 [−0.022, +0.027].**

**So "about 11% of between-event variance is genuinely predictable" is WRONG and is
withdrawn.** The correct statements are:
- there is real stratosphere→surface predictability in winter (~0.066 pre-onset,
  ~0.40 post-onset, within-climate);
- **an SSW adds essentially nothing to it** (+0.004 pre-onset, +0.022 post-onset,
  neither distinguishable from zero — intervals below);
- **94.5% of post-onset "diagnostic skill" is reproducible with no SSW present**
  (0.3792 of 0.4013), measured directly rather than inferred from a tier gap.

Against this, the DW/NDW split's implied R² of 0.63 stands against an event-specific
share whose largest 95% upper bound in any tier is 0.092 — **at least 6.8× smaller**
(12× for the pre-onset tier), not 5.5×.

**My earlier retraction of "there is nothing to classify" was itself premature.**
The defensible claim is narrower and sharper than either version: *the surface
response following an SSW is no more predictable from the stratosphere than that
of an ordinary winter day.*

Stress tests that this survived (`headline_stress_test.json`):
- **mediation objection**: P3 scores 0.386 on pseudo-events, 91% of its real value,
  so the post-onset gain is structural covariance and not downward mediation;
- **weak-coupling objection**: CV R² does not scale with a model's coupling
  strength (slope on \|coupling\|, r = +0.085), so the CMIP under-coupling bias
  does not suppress the estimate;
- **in-sample objection**: this model's optimism is only +0.015, so the published
  0.63 cannot be explained as ordinary overfitting.

#### J-REVISED with intervals (2026-09-24) — the numbers to quote

The table above compared the real arm against **6 pooled pseudo draws: 11,321
pseudo-events against 1,888 real ones**. A ridge fit on 6× the data scores higher
out of sample for that reason alone, which inflates the baseline and biases
"event-specific" toward zero. `within_model_check.py` now compares the real arm
against **1,000 independent pseudo draws, each the size of the real set** (same
members, same calendar days, event-cleaned days), and adds a member-cluster
bootstrap (1,000 resamples of the 20 members). Every draw and replicate has its own
seeded stream, so the numbers are identical at any worker count (checked: 4 vs 12).

| tier | real, within-model | null mean (sd) | **event-specific** [95% from null] | p(null ≥ real) | cluster bootstrap 95% |
|---|---|---|---|---|---|
| P1 pre-onset | 0.0667 | 0.0629 (0.0121) | **+0.004** [−0.022, +0.027] | 0.36 | [−0.048, +0.051] |
| P2 at-onset | 0.0738 | 0.0748 (0.0127) | **−0.001** [−0.028, +0.023] | 0.52 | [−0.060, +0.050] |
| P3 + post-onset strat. | 0.4013 | 0.3792 (0.0189) | **+0.022** [−0.015, +0.060] | 0.11 | [−0.054, +0.092] |

- **The withdrawal holds.** No tier's event-specific share differs from zero. The
  largest upper bound in any tier, P3's bootstrap +0.092, is 6.8× below the 0.63
  the DW/NDW split implies; the pre-onset tier's, +0.051, is 12× below it.
- **The 6-draw baseline was biased by size, by +0.005 to +0.010** — real, and
  smaller than feared. The old event-specific column was that much too low.
- **Two intervals, and they differ for a reason.** The null interval holds the 20
  simulation members fixed; the bootstrap also resamples WHICH members, so it
  carries between-member (largely between-model) heterogeneity and is about 2.5×
  wider. The 20 members come from 11 models, 10 of them CanESM5; "within-model"
  in the tables here means standardised within MEMBER. Quote the bootstrap as the
  conservative bound.
- **Imperfect matching errs against the claim.** Every draw has the real onset
  count, but in 1,017 of 20,000 member-draws some pseudo onsets' +8..+52 day windows
  lack enough in-season days, so their outcome is missing and the draw has fewer
  usable rows than the real set (all 1,888 real outcomes are present). Fewer rows
  lower the null, which raises "event-specific" — conservative for the conclusion.
- **Reproduction is partial, and the cause is known in part.** Re-running the
  original 6-draw arm gives within 0.0667 (was 0.0664) — reproduced — but pooled
  0.1026 (was 0.0998, the drift logged in `FAILURES.md`) and 6-draw pseudo
  0.0694 (was 0.0646) on a 2-event difference in the draws (11,321 vs 11,323).
  That a single 6-draw baseline moves by 0.005 is itself why the 1,000-draw null
  replaces it. Seeds, K, N_BOOT and package versions are now in the JSON.

### J — predictability ceiling, ORIGINAL POOLED NUMBERS (superseded by J-REVISED)
CMIP6, 1888 events, 20 members as CV groups, GroupKFold by member, ridge +
gradient boosting, permutation null, power curve, 86 features.

| tier | overlaps response window? | CV R² | p |
|---|---|---|---|
| **P1 pre-onset** | no | **+0.115** | 0.000 |
| **P2 at-onset** | no | +0.113 | 0.000 |
| **P3 + post-onset stratosphere** | **yes, 22 d** | **+0.422** | 0.000 |

**73% of apparent diagnostic skill is unavailable in advance.**
Power: detection 0.70 at injected R²=0.02, 1.00 at 0.05.

Implied R² of published splits, `q(1−q)C²/Var`, **matched variance**:

| split | C | q | implied R² | vs 0.115 |
|---|---|---|---|---|
| CMIP6 DW−NDW | 1.143 | 0.50 | **0.629** | 5.5× |
| Karpechko AO, obs | 1.782 | 0.70 | **0.597** | 5.2× |
| ACP 26,3723 NAO (published; different index, unmatched denominator — scale only) | 0.708 | 0.635 | 0.104 | 1.0× |
| ERA5 1000 hPa | 0.684 | 0.54 | 0.104 | 0.9× |
| ERA5 850 hPa | 0.639 | 0.59 | 0.088 | 0.8× |

**Not every published contrast is an overclaim on this metric** — the two ERA5
contrasts imply almost exactly the achievable variance. Only the large ones exceed
it. (They remain compromised for the *separate* reason that they are outcome-based;
β×ΔS puts them at 93–102% selection.)

> **The artifact now backs this table (re-run 2026-09-17).** Until today the table
> above existed only in this document: `predictability_ceiling.json` still held the
> pre-correction values, including the impossible **R² = 1.2852** for the
> observational Karpechko contrast, produced by dividing it by the CMIP6 variance.
> The script had been fixed and never re-run. Re-running it reproduces every figure
> in the table — 0.6294, 0.5972, 0.1041, 0.0885, 0.1565 — on Var(CMIP6) = 0.519 and
> Var(obs) = 1.117, and adds the two ERA5 rows, which the old JSON did not contain
> at all. The ACP NAO entry also moved, 0.3368 → 0.1565, for the same denominator
> reason.
>
> **The cross-validated R² moved at the third decimal and I could not establish
> why.** Same 1,888 events, same 38 features, same 20 groups, same fixed seed
> (20260803), same matched variances — but ridge P1 0.0998 → **0.1026**, P2 0.1121 →
> 0.1133, P3 0.4245 → **0.4226**, and gradient boosting up to 0.0081. The most
> likely cause is a library version difference: the original run's package versions
> were never recorded, so this cannot be demonstrated, only suspected. No conclusion
> moves — the tier contrast is (0.4226 − 0.1026)/0.4226 = **75.7%** of apparent
> diagnostic skill unavailable in advance, against the 73% quoted above for the
> 86-feature wave-driving run. From now on the environment is pinned
> (`environment/requirements.lock.txt`) so the next re-run is answerable.
>
> **2026-09-25:** the ACP row used C = 0.850, which subtracted the BOTH subtype
> alone (see §4, Lu & Rao). With the all-DW contrast 0.708 and their DW fraction
> 0.635 it implies 0.104 — and because their NAO index is not the AO whose
> variance is the denominator, the row is for scale only.

**Consistency.** The placebo-corrected σ_f upper bound (0.256, re-run 2026-09-25)
implies R² ≤ **0.126**; pooled cross-validation measures **0.103**. The two methods
are consistent — the measured value sits inside the bound — but the earlier
"agreeing to 0.01" (0.110 vs 0.100) does not survive the re-run and is withdrawn.

**Wave driving closed the last gap.** CMIP6 has no 3-D fields, so v'T' was proxied
via the TEM identity (div F enters through ∂ū/∂t): windowed tendency, 10−100 hPa
shear, 75N−45N vortex geometry. P1 0.0998 → **0.1148**; P2 +0.001; P3 −0.003. Wave
driving helps **only before onset**, which is physically right.

**Observations cannot settle it.** No skill in any catalogue — primary n=42
(R²=−0.081), union n=46 (−0.011), consensus_half 41 (+0.012), era5 41 (+0.002),
jra_55 40 (−0.020), all p > 0.2. Power at n=42 is **0.35** against a true R² of
0.10. **The observational null is uninformative, not negative, and must never be
cited as evidence against predictability.**

### K — SNAPSI causal effect (`snapsi_causal_effect.py`)
`S = mean(nudged) − mean(control)`, the causal contrast by experimental design.
Both arms carry identical nudging machinery; only the target differs (observed
event vs 1979–2019 climatology). Differencing against `free` would confound the
anomaly with the act of nudging.

**Re-run 2026-09-17 on 9 usable centres** (10 downloaded, NRL excluded by the
duplicate guard at ratio 0.004), 2 NH events × 2 initialisations, 38–52 members
each; re-run 2026-09-23 with the 10 CNR-ISAC `nudged s20190108` members the
manifest had mis-versioned (40 → 50 in that case), and 2026-09-24 after the
time-origin fix (UKMO and Meteo-France files start at 06 UTC, and every lead had
been measured from the file's first step — their windows sat 6 h early in K, L,
M, N and O; figures below are current). The original 3-centre figure is reproducible: CCCma, KMA and SNU alone give
**+0.976 σ** against the published +0.972.

| centre set | mean S | sd | note |
|---|---|---|---|
| original 3 (CCCma, KMA, SNU) | **+0.976 σ** | — | reproduces the published +0.972 |
| all 9 | **+1.363 σ** | 0.808 | |
| **9 excluding ECCC** | **+1.107 σ** | **0.271** | the robust figure |

**ECCC is a statistical outlier and must not be averaged in silently.** Its
S/σ is **+3.407** when the next highest is ECMWF at +1.467 and the other eight
span +0.618 to +1.467. A Tukey test on the other eight puts the upper fence at
1.669, so ECCC sits well outside it, at 2.9× their median. Including it moves the
grand mean from +1.11 to +1.36 and triples the between-model sd, 0.271 → 0.808.
**Quote +1.11 σ (or the all-9 median +1.24) as the causal effect, with the all-9 mean of
+1.36 as the ECCC-inclusive sensitivity.** Whether ECCC's response is physical or
a submission problem is not established here; its duplicate-guard ratio is 2.501,
the highest of any centre, which is the opposite of NRL's failure mode and is not
by itself evidence of a defect.

**Between models and between events, on the same basis.** Excluding ECCC,
between-model sd 0.271 against a between-event difference of **0.142 σ**
(Feb-2018 1.178, Jan-2019 1.036): about 2×. Including ECCC, 0.808 against 0.054:
about 15×. *(An earlier version paired the ECCC-excluded sd with the
ECCC-included event difference and called it "far more"; found in review
2026-09-25.)* Models disagree about downward-coupling strength more than the two
events differ, but by ~2×, not by an order of magnitude, on the headline basis.

Implied NDW response given S = +1.363: 1 of 5 published contrasts forces NDW to
**+0.116** (effectively zero); the other four imply **+0.79 to +0.99 σ** — events
labelled "non-propagating" would still carry a large real surface response. That
is a continuum, not two populations.

### L — THE DECISIVE TEST: selection under experimentally fixed forcing (`snapsi_selection_test.py`)

Every other demonstration in this project that "the contrast is manufactured"
rests on a *modelled* null — pseudo-events, or CMIP6 with an assumed structure.
SNAPSI removes that dependence. In the `nudged` runs **every member's
stratosphere is nudged to the same observed evolution**, so members differ only
in tropospheric noise: the causal driver is identical *by construction*, not by
assumption. Apply Karpechko conditions 1–2 (the surface conditions) and split.

| arm | what is true of it | DW−NDW contrast | DW rate (unconditional) |
|---|---|---|---|
| **nudged** | stratospheric forcing **identical** across members | **−1.622 σ** (sd 0.303, n=25) | **0.842** |
| **control** | **no SSW forcing at all** | **−1.590 σ** (sd 0.089, n=36) | **0.450** |

> **SOUTHERN HEMISPHERE, added 2026-09-17 — L replicates there too.**
> Using `s20190829`, the one SH initialisation whose forecasts span the +8..+25
> window (`s20191001` starts 13 days after the central date and is excluded by
> the coverage guard rather than averaged over a shorter stretch):
>
> | arm | contrast | DW rate |
> |---|---|---|
> | nudged | **−1.958 σ** (n=7 of 8) | 0.786 |
> | control | **−1.584 σ** (n=8 of 8) | 0.413 |
>
> Same signature as the NH: a large contrast in both arms, and the rate doing
> the separating. **Two caveats travel with every SH number.** The central date
> is 18 September 2019 from the SNAPSI protocol, where the 10 hPa 60°S wind
> *"did not reverse"* but reached its minimum — so this is an **austral MINOR
> warming, not an SSW under the WMO definition**. It tests generality to a
> different hemisphere *and* to a weaker class of event, which is more than the
> NH result covers but is not the same event mirrored. And it rests on one
> initialisation, so it is a replication, not an independent second sample.
>
> NH and SH are never pooled: they are summarised in separate blocks.

> **Re-run 2026-09-17 on the COMPLETE archive: 9 usable centres, 4 NH
> initialisations, 36 ensembles per arm, against the 3 centres and 9-12
> ensembles the original rested on. The result holds and tightens.** The
> contrast is the same with the SSW held identical and with no SSW at all
> (−1.622 vs −1.590), while the rate separates as before.

> **The control contrast is the value a threshold gives on pure noise.** Control
> members are standardised by their own ensemble, so they are close to a
> unit-variance Gaussian with mean 0. Splitting such a distribution at its mean
> gives a difference of group means of exactly 2√(2/π) = **1.596 σ**; the
> measured control contrast is **−1.590 σ** (sd 0.089 across 36 ensembles). The
> DW/NDW contrast is therefore not an effect size that the atmosphere supplies —
> it is the half-normal mean difference, a property of splitting a distribution
> at a threshold, and it appears at full size with no SSW. (Fig. 2a. The nudged
> arm cannot be read off the same line exactly: its members are shifted, and the
> 11 ensembles with < 3 NDW members, DW rate 0.99, cannot form a contrast.)
>
> **One correction the larger sample exposed.** The contrast needs ≥3 members in
> BOTH groups, which is asymmetric by construction: only an arm with a strong
> forced shift can push nearly every member into one class. **11 of 36 nudged
> ensembles are dropped that way, with a mean DW rate of 0.993; 0 of 36 control
> ensembles are.** Taking the rate over surviving ensembles understated it at
> 0.776. The unconditional rate is **0.842 against 0.450**, so the rate
> separation is *larger* than previously reported. The contrast remains
> necessarily conditioned on the ensembles where both groups exist, which are
> the weaker-responding ones — stated rather than hidden.

**The contrast is the same whether or not an SSW happened** (−1.52 vs −1.56).
It cannot be reporting a physical difference between groups: in the nudged arm
there is none to report, and in the control arm there is no event. The contrast
is the cut.

**The rate is where the signal lives** — 0.78 with forcing against 0.44 without.
This independently reproduces, in a designed experiment, what
`FINDING_selection_on_outcome.md` §6 found in observations (70% real vs 38%
pseudo, P < 0.0001). *The classification rate is real and uncontaminated; the
selected group's magnitude is not.*

**NRL excluded, and the guard was fixed to catch it.** Its nudged-minus-control
ensemble-mean effect is identically zero at all four initialisations, and for
`s20181213` its nudged submission is an exact duplicate of control. A first
version of the guard compared members *by label* and passed NRL, whose members
differ while the ensemble mean does not. The guard now tests the effect itself.

**Units caveat, stated because this is exactly the trap that produced an
impossible R²=1.285 earlier:** these contrasts are standardised by the control
ensemble's *member-to-member* spread for one event; the published −0.708 is
standardised by the *between-event* spread in observations. Different
yardsticks — the comparison indicates scale, not a like-for-like fraction.

### M — the SSW shifts the distribution, it does not split it (`snapsi_distribution_test.py`)

Result **I** bounded the forced between-event variance at σ_f² ≈ 0 in CMIP6, but
only after an **84% placebo correction** — 84% of the raw excess variance was
already present before onset. That correction is the obvious line of attack.
SNAPSI measures the same quantity **causally, with no correction at all**:
`control` has no SSW, `nudged` has the same SSW in every member.

Pooled over case-ensembles:

| quantity | 3 centres (original) | **9 centres, 36 ensembles (2026-09-17)** |
|---|---|---|
| causal shift | −0.889 σ | **−1.572 σ** (sd across cases 1.081) |
| **variance ratio** nudged/control | 0.939 [0.781, 1.130] | **0.955 [0.873, 1.051]** |
| pure translation (KS) | p = 0.805 (451 v 450) | **p = 0.371** (1798 v 1805) |
| excluding the short-lead init | — | **0.995 [0.895, 1.106]**, KS p = 0.650 |

> **The conclusion survives a fourfold expansion — but only after an estimator
> defect was fixed, and the naive scale-up would have overturned it.**
>
> Each ensemble is standardised by its OWN control mean and sd, so the control
> arm lands at ~0 in every case while the nudged arm lands on that case's causal
> shift. Concatenating without re-centring therefore adds the between-ensemble
> spread of those shifts to the nudged variance and to nothing else. At 3
> centres that spread was 0.295 σ, contributing 0.087 to a ratio of 0.941 — the
> published 0.939, reproduced. At 9 centres the spread is **1.071 σ**,
> contributing 1.15, and the uncorrected ratio reads **2.085 [1.922, 2.264]**
> with KS p = 0.000: a flat contradiction of the result, produced entirely by
> the pooling.
>
> Re-centring each ensemble before pooling — which is what the within-ensemble
> quantity requires — gave **0.952 [0.869, 1.047]** across all 9 centres and
> **0.995** excluding the short-lead initialisation. The arithmetic closed on that
> run: 0.952 within + 1.147 between ≈ 2.085 pooled. (With the 10 recovered
> CNR-ISAC members and the 06 UTC time-origin fix: **0.955 [0.873, 1.051]**.)
>
> **`shift_sd_across_cases` = 1.081 σ is NOT σ_f.** It mixes between-model with
> between-event variation, since the 36 ensembles span 9 centres and only 2
> events. It must never be quoted as an event-to-event forced spread.

**Why this matters.** If coupling strength genuinely varied between events — if
some SSWs were "downward-propagating" and others not — the forced ensemble would
be *more dispersed* than the unforced one. It is not. With the stratospheric
forcing identical across members by construction, the added between-member
variance is bounded at **+0.13 σ² at 95%**. Two independent routes now agree —
one needing an 84% placebo correction (I), one needing none (M).

> **POWER LIMITATION, stated not buried.** The Gaussian-mixture BIC test in this
> script has **essentially no power**: at a proper 5% critical value it detects a
> 2/3 mixture with probability 0.04 at 0.25 σ separation and **0.00 at 1.0–1.5 σ**,
> because the within-component spread stays near 1 σ while the separation grows.
> Its p = 0.415 is **not** evidence for one population and must never be quoted as
> such. **The variance bound and the KS test carry this result; the mixture test
> does not.** A first version of the power code compared the statistic against a
> *single* null draw — a coin flip under the null — and produced power
> *decreasing* with separation, which is impossible and is what exposed it.

> **The guard that decided the answer.** Including the `s20190108` initialisation
> flips every conclusion: variance ratio 1.644 [1.41, 1.93], KS p = 0.010. That
> initialisation begins **6 days after onset**, so the +8..+25 window sits at lead
> 2–19 days, before the ensemble has diverged; its control spread is ~138 Pa
> against 205–466 Pa elsewhere. Standardising by forecast-spread growth instead of
> climatological spread inflates both statistics. The anomaly was flagged from the
> spread *before* the run and confirmed by the lead arithmetic afterwards — and it
> is a property of the **initialisation**, present at both centres, not of a centre.

---

### N — Loeffel head-to-head: is the coupling a property of the event? (`snapsi_loeffel_test.py`)

Loeffel et al. (2026, WCD 7, 895–913) correlate week-2 100 hPa polar-cap GPH with
the weeks 3–7 surface response **across 18 events, on ensemble means** (r = 0.85)
and conclude that SSWs differ in their capacity to couple downward. In a SNAPSI
`nudged` ensemble the event is identical in every member, so the same relation can
be measured **between members**, where no event-to-event difference exists; the
`control` arm measures it with no SSW at all.

**Final run 2026-09-25 on the complete archive:** 9 centres, 4 NH
initialisations, both arms, 36 matched ensembles (5,210 zg members; leads on
the corrected 00 UTC origin). Supersedes the 8-centre run of 2026-09-23
(Δz −0.001 on 28 ensembles) and the interim of 2026-09-24.

**Design gate (pre-specified, sd ratio > 0.5): passes**, median 0.964 over 36
ensembles — nudging above 50 hPa does not pin 100 hPa.

Primary window (post-onset days 15–25; days 15–49 is covered by only 10 of 36
ensembles), Fisher-z pooled:

| | nudged (SSW imposed) | control (no SSW) |
|---|---|---|
| all 36 matched ensembles | **+0.111** [+0.063, +0.158], 1,797 members | **+0.088** [+0.041, +0.135], 1,805 members |
| **difference Δz** | **+0.023 [−0.044, +0.090], p = 0.50** | |
| excl. `s20190108`, 27 matched | +0.109 vs +0.108, **Δz +0.002 [−0.076, +0.079], p = 0.97** | |

Loeffel's own days 15–49 window, on the 10 matched ensembles that reach it (9 of
them the short-lead `s20190108`): +0.020 vs +0.070, Δz −0.049 [−0.178, +0.079].
The max-available window now resolves to days 15–25, identical to the primary.

**Reading.** Between members, the 100 hPa → surface relation is weak (r ≈ 0.1)
and **the same with no SSW present**: a general lower-stratosphere–surface
persistence, not something an SSW adds. The pre-specified falsifier — the control
arm reproducing the nudged correlation — is met. The per-ensemble r vary more
than sampling alone allows **in both arms** (Cochran Q p = 0.009 nudged, 0.025
control; 0.036 and 0.052 excluding `s20190108`): the coupling depends somewhat on
model or initialisation, but that dependence is not specific to the SSW.

**What it does NOT test.** Loeffel's between-event coefficient: SNAPSI has 2 NH
events. Across the 36 ensembles (centre × initialisation) the correlation of
ensemble means is **+0.641 nudged and +0.834 in control** — correlations of
Loeffel's size arise across **models**, and are larger with **no SSW at all**.
Descriptive only (raw means, not anomalies; the spread is between models); never
to be quoted as a test of their claim. It shows that a large across-ensemble-mean
r does not by itself require event-to-event differences.

**Power.** At n = 50, a within-ensemble r = 0.3 is detected with probability 0.56
per ensemble; pooled over 36 matched ensembles the Δz interval half-width is 0.067.

### O — the two archetypes: different events, or two draws? (`snapsi_archetype_test.py`)

Rao, Garfinkel & White (2020, JGR-Atmos 125, e2019JD031919) present Feb-2018 as a
strong downward-propagating SSW and Jan-2019 as one with a weak or opposite
surface response, attributing the difference chiefly to SSW strength; Nebel et
al. (2024, GRL 51, e2024GL110529) frame such outcomes as potential forecast
busts. SNAPSI nudges each model's stratosphere to each observed event, so the
forced odds of each outcome can be read off directly.

**1. The forced odds are the same; the forced shifts are comparable.** On the
headline basis (8 models, ECCC excluded), averaged over each event's two
initialisations: causal shift (K) 1.178 σ (Feb-2018) vs 1.036 σ (Jan-2019); DW
probability (L) 0.817 vs 0.828. At the primary pair plotted in Fig. 4
(s20180125 / s20181213, onset at lead 18 / 20 d): 1.014 vs 0.938 σ, DW odds
0.778 vs 0.715, and each event has the larger shift in 4 of 8 models. Jan-2019's
forced response is ~7-12% smaller — nowhere near the opposite outcome the
archetype reading attributes to it. *(The first version of this point said
Jan-2019 had the LARGER shift, +1.35 vs +1.19 σ: that was the all-9 mean, and
ECCC alone (4.65 vs 2.63) produced the ordering. Found in review 2026-09-25.)*

**2. The observed outcomes sit inside their own ensembles.** ERA5 polar-cap psl
(60–90N including the 60N row, 6-hourly, at exactly the members' forecast
times), placed in each model's nudged distribution for days +8..+25. A < 0 is
the downward sign, so a "weak" outcome sits in the TOP tail:

| init | event | obs percentile in nudged, median (range) | models with obs < 2.5% / > 97.5% |
|---|---|---|---|
| s20180125 | Feb-2018 | 0.43 (0.00–0.76) | 1 / 0 (0 / 0 after offset adjustment) |
| s20181213 | Jan-2019 | **0.68** (0.12–0.88) | **0 / 0** |
| s20190108 (short lead) | Jan-2019 | 0.96 (0.00–1.00) | 2 / 3 |

The falsifier — observed Jan-2019 above the 97.5th percentile (weaker than its
forced distribution) in ≥5 of 9 models — is met in 0 of 9 at the primary init
and 3 of 9 at the short-lead init. *(The first write-up named the 2.5th
percentile, the wrong tail for this sign convention; corrected after review.)*

**3. Their observed difference is ordinary.** A_2018 − A_2019, observed, against
every pairing of a 2018 member with a 2019 member of the same model (primary
pair s20180125 / s20181213, the two inits with comparable lead to onset): the
observed difference (−1.05 to +1.03 σ across the 9 model frames) lies inside
the central 95% in **9 of 9 models** (median percentile 0.33), against
member-pair sd 1.2–1.6 σ. The short-lead pair (s20180208 / s20190108) is outside
in 1 of 9 (2 of 9 after offset adjustment).

**4. The "non-propagating" label is a threshold call.** The project's ERA5
implementation of Karpechko et al. (2017), conditions 1–3:

| Jan-2019 | NAM mean | days negative | class |
|---|---|---|---|
| 1000 hPa, days +8..+52 (Karpechko) | **+0.019 σ** | 62% | **NDW** — condition 1 fails by 0.019 σ |
| 850 hPa, days +8..+52 (ACP 26, 3723) | −0.041 σ | 62% | DW |
| 1000 hPa, days +8..+25 | −0.214 σ | 83% | DW |
| 850 hPa, days +8..+25 | −0.199 σ | 83% | DW |

Feb-2018 is DW under all four. **The field's canonical non-propagating event is
non-propagating in one of four reasonable variants of the criterion, by 0.02 σ.**

**Reading.** The pair the field uses to illustrate that SSWs differ in their
downward impact received the same forced response and the same odds; their
observed outcomes are both ordinary draws; and the one label that separates them
sits on the threshold. Nothing here requires an event property.

**Limits.** Two events. Model–ERA5 representation offsets are large for some
centres (ECCC −1,213 Pa at s20180125 init; the adjusted column removes the
initial-time offset but not lead-dependent drift), which is why the event-pair
difference — where a model's bias largely cancels — is the robust test. One
observation per event, so the 9 per-model percentiles are not independent.

## 4. Prior art — verified from source, all must be cited

1. **Karpechko et al. (2017)**, QJRMS 143, 1459, doi 10.1002/qj.3017 — the
   criterion. Days +8..+52: mean 1000 hPa NAM negative; fraction of days negative
   > 0.5; fraction with negative 150 hPa NAM > 0.7. Conditions 1–2 **are** the
   surface response. 229 citations and live.
2. **Coughlin & Gray (2009)**, JAS 66, 531, "A Continuum of Sudden Stratospheric
   Warmings" — continuum on the SSW-*definition* axis.
3. **Maury et al. (2016)**, JGR-Atmos, doi 10.1002/2015JD024226 — *"the idea of a
   'warming continuum'"*; *"there is no statistical difference between SWEs with
   regard to their feedbacks on planetary waves and hence their potential influence
   into the troposphere."*
4. **White et al. (2019)**, J. Climate 32, 85, §3b — DW−NDW positive-lag
   differences *"entirely there by construction"*. **The closest prior art.**
5. **Baldwin et al. (2021)**, Rev. Geophys. 59, e2020RG000708, §7.2 — *"about two
   thirds … of SSW events are characterized as having a visible downward impact"*,
   with no uncertainty range and no null comparison.
6. **Loeffel et al. (2026)**, WCD 7, 895–913 (preprint EGUsphere 2025-4164), ICON ensemble
   re-forecasts, case-to-case variability. **Their r = 0.85 survives our control**
   (P(null ≥ 0.85) = 0.003), though the missing baseline is real and large: with no
   event at all, window overlap alone gives r = **+0.530** [+0.162, +0.789].

7. **Rao, Garfinkel & White (2020)**, JGR-Atmos 125, e2019JD031919 (read in
   full via PMC7507786) — Feb-2018 vs Jan-2019 in S2S models; *"the strength of
   the SSW is more important than the vortex morphology in determining the
   magnitude of its downward impact."* Result O tests this pair directly. Also
   verbatim from the full text: for Feb-2018 *"North Eurasia is ~4 °C colder than
   normal"*; for Jan-2019 *"The 2-m temperature in North Eurasia is anomalously
   warm"*.
8. **Nebel, Garfinkel, Cohen, Domeisen, Rao & Schwartz (2024)**, GRL 51,
   e2024GL110529 — 7 S2S models, 16 SSWs 1998–2022: models predict which SSWs
   have a stronger downward response to 100 hPa but *"struggle to predict which
   have a stronger tropospheric response"*; frames non-propagating outcomes as
   forecast busts. **Only the abstract has been read** (the publisher and
   preprint server refuse automated fetches); do not cite its internals until
   the full text is read.

9. **Sigmond, Scinocca, Kharin & Shepherd (2013)**, Nature Geoscience 6, 98–102,
   doi 10.1038/ngeo1698 (abstract read from the publisher) — forecasts initialised
   at SSW onset *"faithfully reproduce the observed mean tropospheric conditions in
   the months following"*, with enhanced skill against forecasts not initialised
   during SSWs, *"for atmospheric circulation patterns, surface temperatures over
   northern Russia and eastern Canada and North Atlantic precipitation"*. The skill
   claim this project decomposes.
10. **Hitchcock et al. (2022)**, GMD 15, 5073–5092, doi 10.5194/gmd-15-5073-2022
   (read in full) — SNAPSI protocol: zonal-mean nudging, *"tapering gradually from
   infinite (i.e., no nudging) below a lower limit of p_b = 90 hPa, to full
   strength at p_t = 50 hPa"*, 6-h timescale at full strength, all latitudes; the
   control ensemble amounts to *"a 'climatological' stratospheric forecast"*.

11. **Lu & Rao (2026)**, ACP 26, 3723–3742, doi 10.5194/acp-26-3723-2026 (read
   from the journal page) — ERA5 1940–2022: 52 SSWs, 33 DW (13 BOTH, 14 EA, 6 NA)
   and 19 NDW, criterion of White et al. (2019) with 850 hPa. 60-day NAO means:
   BOTH −0.762, EA −0.567, NA −0.435, NDW +0.088. **All-DW mean −0.620, so the
   DW−NDW contrast is −0.708, DW fraction 33/52 = 0.635.** The −0.850 this project
   used until 2026-09-25 subtracted the BOTH subtype alone.
12. **Loeffel et al. (2026)** details used in N, from the published text (recorded
   in `snapsi_loeffel_test.py`): predictor = week-2 (days 8–14) 100 hPa polar-cap
   standardised GPH anomaly; response = weeks 3–7 1000 hPa polar-cap GPH; r = 0.85,
   18 events, across events on ensemble means; conclusion *"pronounced and robust
   event-to-event differences in the tropospheric response to SSWs"*.

**What remains novel:** the quantification (β×ΔS, R²=0.99–1.00 in two systems); the
extension to stratospheric stratifiers; the power-bounded distributional test; the
shift predicting "two thirds"; and above all **the out-of-sample predictability
ceiling and the pre/post tier contrast**, which nobody has measured.

**Never lead with "unrecognised."** The field has circled this for 15+ years.

---

## 5. Open gaps

| gap | status |
|---|---|
| **forecast/operational circularity** | **UNANSWERED** — died twice on agent session limits. The one lead that could change the venue category. |
| **R re-validation** | `validate_R_diagnostic.py` must add a precursor arm (anomaly concentrated *before* onset) before any R-based claim is made again. |
| true v'T' in CMIP6 | proxy only; needs hundreds of GB of daily 3-D va/ta |
| observational power | 35% at R²=0.10; n=46 is the whole record |
| SNAPSI NH events | 2. Not fixable. σ_f not estimable there. |
| manuscript n inconsistency | manuscript uses 47 (union), analyses use 43 primary / 42 usable — **must be one choice** |
| `ncep_ncar` catalogue variant | crashes the CV loop (empty permutation array) |

### 5b. Event counts are NOT inconsistent — they are named catalogue sets

Every count that appears in this project traces to a set frozen in
`02_event_catalogues/FREEZE_RECORD.md` and materialised in `catalogue_sets.json`.
There is no unreconciled discrepancy:

| n | set | note |
|---|---|---|
| **43 / 36 winters** | **`primary`** | the frozen preregistered catalogue — **use this** |
| 47 | `union` | ≥1 reanalysis; permissive. What the stale manuscript uses |
| 42 | `consensus_half`, and `era5` | also the usable n in the observational CV arm (43 − 1958-01-30, whose predictors predate the record) |
| 41 | `primary_compendium_only`, `jra_55` | |
| 39 | `ncep_ncar` | **also the pre-re-freeze `primary`** — the collision that looks like an error |
| 36 | `consensus_strict` | equals primary's *winter* count by coincidence |
| 30 / 29 / 27 | AAO arm / satellite-era & isolated subsets / `merra2` | record-length limited |

**The one real problem is stale documentation, not stale data.** The catalogue was
re-frozen 2026-07-30 (39 → 43 events) after the consensus rule was found to be
era-dependent. The JSONs were all regenerated; several markdown files were not, and
still describe "frozen primary catalogue (39 events / 33 winters)":

**All six were resolved on 2026-09-17.** The list is kept because two of them were
not merely stale — their conclusions had been reversed by the re-freeze:

| document | what was wrong | now |
|---|---|---|
| `05_corrected_estimators/FINDING_canonical.md` | headline **reversed** | rewritten from the JSON |
| `05_corrected_estimators/FINDING_multi_index.md` | two conclusions **reversed** | rewritten from the JSON |
| `07_physical_decomposition/FINDING_gate8.md` | 39-event numbers, "30 of 39" | refreshed; no conclusion changed |
| `06_simulation_validation/FINDING_gate3_contamination.md` | quoted the withdrawn 0.947 / 0.053 | corrected to 0.930 / 0.070 |
| `06_simulation_validation/GATE3_CLOSURE.md`, `GATE3_FINAL.md` | already carried their own corrections | verified, left as they stand |

Three further documents quoted the withdrawn 0.947 / 0.053 pair as current and were
corrected at the same time: `GO_NO_GO.md` condition 3, `FINDING_design_robustness.md`,
and `06_simulation_validation/FINDING_spatial.md`.

**`FINDING_canonical.md` was reversed by its own JSON, and has been rewritten from
it (2026-09-17).** It had stated "under the preregistered primary catalogue the AO
pre-onset bins are not significant (−30..−16: −0.57, n.s.)" and built a finding on
pre-onset significance being catalogue-dependent. The current
`canonical_event_study.json` (43 events, 3 harmonics, 1,200 replicates) gives
**−0.8195, p = 0.005 — significant**.

**The corrected finding inverts the old one:** pre-onset significance is *not*
catalogue-dependent. AO at −30..−16 is significant in **8 of 9** catalogue sets
(p = 0.005 under `primary`), and NAO at the same bin in 8 of 9 (p = 0.003–0.018).
The sole exception in both outcomes is `merra2`, which is **record length, not
event definition** — MERRA-2 starts in 1980, giving 27 events / 23 winters against
43 / 36, and its post-onset bins remain significant and of ordinary size.
`consensus_strict` (36 events) shows the same attrition more mildly. This agrees
with A/B: catalogue choice moves the AO estimate 0.115 σ once the window is fixed.

A second stable feature was not previously recorded: **AO peaks at +15..+29 in 8 of
9 catalogues while NAO peaks at +30..+44 in all 9** — the NAO response lags the AO
response by about two weeks, in every set including the two shortest.

**`FINDING_multi_index.md` was reversed the same way, on two counts, and has been
rewritten (2026-09-17).** On the 43-event catalogue, under the *same* preregistered
32-test BH family:

| outcome | old (39 events) | current (43 events) |
|---|---|---|
| AO | 3/8 survive, no pre-onset | **5/8**, incl. −30..−16 (q=0.027) |
| NAO | 3/8 survive, no pre-onset | **4/8**, incl. −30..−16 (q=0.049) |
| PNA | **0/8** | **3/8**, all **positive**: −30..−16 +0.354 (q=0.030), −15..−1 +0.354, +15..+29 +0.283 |
| AAO (control) | 0/8 | **0/8**, own family, min q = 0.24 — still passes |

Two statements are therefore withdrawn: *"no pre-onset bin survives FDR in any
outcome"* and *"the circulation response is annular/Atlantic, not hemisphere-wide —
PNA shows nothing"*. The cause is the catalogue re-freeze, not the test: PNA
survives under the unchanged preregistered family. Consequently **GO_NO_GO
condition 6 flips from NOT MET to MET** — the response is not driven by AO/NAO
alone — and it closes on circulation indices, without the non-circulation outcomes
that were blocked behind the struck census.

This does **not** revive the precursor as causal. A precursor cannot be caused by
the event that follows it, and §7.3 stands: CanESM5 gives −0.387 and MIROC6 +0.380,
both p < 0.0001, so models disagree in sign. What is established is narrower — the
pre-onset depression is not an artifact of catalogue choice and not an artifact of
multiple testing.

### 5c. Artifact integrity

- **Gate 3 — RESOLVED 2026-08-04, and the documentation was the stale part.**
  The gap was misdiagnosed twice. First, `finalise_gate3.py` is the script
  `FINDING_gate3_contamination.md` line 11 labels **flawed** (its `excluded` draw
  leaves real events in the baseline, giving null +0.29); re-running it would have
  written a `gate3_final.json` that contradicts the closure. Second, the corrected
  script `gate3_clean_null.py` shipped configured at **N_BOOT = 150** — the very
  setting the project's own standing rule forbids — so every artifact on disk said
  `passes = False` while three documents recorded Gate 3 as CLOSED on 0.947/0.053.
  **No artifact supported the closure.**

  Recomputed at N_BOOT = 1000, N_COVER = 300 on the 43-event catalogue
  (8,862 clean days / 73 winters):

  | spec | null mean | coverage | FPR | pass |
  |---|---|---|---|---|
  | **3 harmonics (primary)** | **+0.0052** [−0.0062, +0.0166] | **0.930** (0.895, 0.956) | 0.070 | **yes** |
  | 6 harmonics | +0.0227 | 0.930 (0.895, 0.956) | 0.070 | yes |
  | 12-knot spline | +0.0397 | 0.937 (0.903, 0.961) | 0.063 | yes |

  **Gate 3 closes — but three quoted facts were wrong.** (a) The disclosed
  "residual positive bias of +0.033 σ, interval excluding zero" is **+0.0052 with
  an interval containing zero**; it was an artefact of the 39-event build and is
  withdrawn. (b) The 12-knot spline was rejected for bias > 0.05; at +0.0397 it
  passes, so that rejection is withdrawn. (c) Coverage/FPR are **0.930/0.070**,
  not 0.947/0.053 — the latter are not reproducible and must not be requoted.
  The ordering that justified 3 harmonics as primary survives (+0.005 < +0.023 <
  +0.040).
- **`predictability_ceiling.json` still holds the superseded implied-R² values**,
  including `Karpechko AO, observations: 1.2852` — an impossible R² > 1 caused by
  dividing an observational contrast by the CMIP6 variance. The script was
  corrected afterwards but never re-run, so the matched-variance numbers in §3J
  exist only in this document. **Re-run `predictability_ceiling.py`.**
- Genuinely stale outputs (script newer than its JSON by more than a few minutes):
  `stratifier_bias_law.py` (2 days), `ensemble_precursor.py` (1 day),
  `marginal_events.py` (6 h). `ao_null_perbin.py` and `predictability_ceiling.py`
  differ by 2–6 minutes, i.e. the script was touched just after it wrote its output.

---

## 6. Provenance

| result | script | output |
|---|---|---|
| A/B | `07_physical_decomposition/design_sensitivity.py`, `05_corrected_estimators/clean_subset.py` | — |
| C | `07_physical_decomposition/` causal-timescale scripts | `FINDING_causal_timescale.md` |
| D | `07_physical_decomposition/mediation_disattenuated.py` | `FINDING_mediation.md` |
| E | `06_simulation_validation/selection_on_outcome.py` | `FINDING_selection_on_outcome.md` |
| F/G | `06_simulation_validation/stratifier_bias_law.py`, `stratifier_law_cmip6.py` | `stratifier_bias_law{,_1958}.json` |
| H | `07_physical_decomposition/is_downward_propagation_a_class.py` | `is_downward_propagation_a_class.json` |
| H | `08_literature_audit/era5_recompute_and_two_thirds.py` | ERA5 body count + two-thirds |
| I | `07_physical_decomposition/forced_variance_ceiling.py` | `forced_variance_ceiling.json` |
| **J** | `07_physical_decomposition/predictability_ceiling.py`, `predictability_wave_driving.py` | `predictability_ceiling.json` |
| K | `03_data_ingestion/acquire_snapsi_surface.py`, `07_physical_decomposition/snapsi_causal_effect.py` | `snapsi_causal_effect.json` |
| Loeffel control | `08_literature_audit/loeffel2025_missing_control.py` | `loeffel2025_missing_control.json` |
| N | `03_data_ingestion/acquire_snapsi_zg.py`, `07_physical_decomposition/snapsi_loeffel_test.py` | `snapsi_loeffel_test.json` |
| O | `03_data_ingestion/acquire_era5_psl_cap.py`, `07_physical_decomposition/snapsi_archetype_test.py` | `snapsi_archetype_test.json` |
| era gate | `06_simulation_validation/era_dependence_gate.py` | `era_dependence_gate.json` |

**Layout changed 2026-08-04.** There is now ONE results tree for the whole
repository, at the repo root, replacing results scattered across six stage
directories plus a separate `archive/superseded_results/`:

```
results/
  current/<theme>/       53 live results, grouped by the question they answer
                         1_catalogue 2_event_study 3_calibration 4_selection_bias
                         5_mechanism 6_predictability 7_ensemble 8_experiment
                         9_literature
  superseded/<family>/   193 retired results (review_rounds, numbered, synthesis,
                         snowpack, stratosphere, sensitivity, limitations, …)
  _ALL_RESULTS.json      all 53 live payloads in a single file — the read path
  _PROVENANCE.json       producer + sha256 baseline, drives the staleness check
```

Scripts stayed in their stage folders (147 import sites resolve modules by stage
path); only the 50 write sites were repointed. `python tools/catalog_repo.py`
rebuilds `CATALOG.md` and flags stale results by producer sha256, not mtime.

Data: NCEP 1958–2024 via THREDDS NCSS (`/thredds/ncss/grid/`, bit-identical to
whole-file, 60× cheaper; 1975–76 return HTTP 500 and need the whole-file fallback).
ERA5 via WeatherBench2 public GCS zarr, anonymous. CMIP6 zonal means, 20 usable
members. SNAPSI via CEDA bearer JWT.

---

## 7. What died, and what it protects

Nine claims were withdrawn or corrected. This section is why the survivors should
be believed.

1. **"Study design reshapes SSW impact estimates"** — the original thesis.
   Estimator 0.002 σ, catalogue 0.115 σ. The project set out to indict the field's
   design choices and found the field's headline result robust to them.
2. **The snow-avalanche application** — abandoned. Still in the stale manuscript.
3. **"The pre-onset anomaly is real"** — withdrawn. CanESM5 gives −0.387, MIROC6
   gives **+0.380**, both p<0.0001. Models disagree in sign; the question is
   model-structural and adding members cannot resolve it.
4. **"Stratospheric (PJO) classification is exempt"** — withdrawn. It manufactures
   +0.889 σ under a true null.
5. **`u10_pre` as a forecast-usable positive** — killed. Window sign flip (−0.706 at
   (−45,−16) vs **+0.178** at (−60,−31)), 29.6-day seasonal-timing confound,
   BH-FDR q = 0.558.
6. **"σ_f ≈ 0"** — corrected. The truth sits at the *upper* end of the bound
   (~0.26 σ after the 2026-09-25 re-run; 0.24 before), consistent with
   cross-validation. The point estimate of 0.081 (0.113 after re-run) was
   misleadingly low.
7. **"Bootstrap under-covers by 20–27%"** — phantom. A sample-size confound:
   √(8638/5121) = 1.299 explained it exactly. **The same confound was re-created in
   CMIP6 and caught again.**
8. **"The response distribution is BIMODAL"** — a broken dip test. The statistic
   measures distance to a *monotone* density, so uniform calibration guarantees
   p=0 for unimodal data. Recalibrating against pseudo-events flipped p from 0.0000
   to 0.846. **Three tests disagreeing is what exposed it.**
9. **"The AAO negative control fails"** — a mis-specified FDR family. Pooling 32
   tests across four indices let the strong positives raise the BH threshold and
   make the *control* easier to call significant. Within its own 8-bin family,
   q = 0.133.

Plus a reproducibility defect: three scripts seeded bootstraps with `hash(str)`,
which Python randomises per process — no p-value from them was reproducible. Fixed
to `zlib.crc32`.

And in this session: **NRL's SNAPSI submission is corrupt** — `nudged` and
`control` are byte-identical for three of four initialisations. An automatic
duplicate guard now blocks it. This is worth reporting to the SNAPSI team.

---

## 8. Venue

Calibrated empirically, not from memory: Nature Geoscience published essentially no
SSW surface-impact work in 2024–26. The field's Nature-family outlets are **Nature
Communications**, **npj Climate and Atmospheric Science**, and **Communications
Earth & Environment**. NGeo's recent climate output is broad observational impact
work (global hail climatology, agricultural drought).

Against that: the headline is model-based, the observational arm has 35% power, and
the framing is corrective over prior art that includes White et al. (2019) saying
the differences are "entirely there by construction."

**Realistic best targets:** npj Climate and Atmospheric Science or Communications
Earth & Environment; JGR-Atmospheres / J. Climate / WCD as strong specialist homes.
This is not a downgrade — it is where this field's readers are.

**What would change the answer:** resolving whether the correction bites on any
operational or forecast-skill claim (§5). If an S2S skill claim is contaminated,
the significance changes category.
