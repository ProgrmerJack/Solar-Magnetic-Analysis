# Consolidated results — SSW surface-impact methodology

**Compiled 2026-08-03. This is the single source of truth for the project.**

It supersedes `archive/superseded_paper/ng_manuscript.tex` (last modified 2026-07-28), which predates
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
its own tail — 93–103% of the published contrast is reproduced by applying the same
criterion to dates with no SSW in them. Measured out of sample within a single
climate, **an SSW adds essentially nothing to the predictability of the surface
response** (−0.006 [−0.033, +0.020] pre-onset, +0.015 [−0.028, +0.059] post-onset,
against a size-matched null of 1,000 event-free draws), while the published
criterion's splits implicitly claim **0.34–0.61**. SSWs do have a large, real surface effect — measured causally
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
| C | R = 0.82 measures temporal **concentration**, not causal fraction | obs + 675 model events | **corrected — see §3C** |
| C2 | Granger vortex→surface asymmetry | obs n.s. (P=0.213); 4 models sig. | **not established in obs** |
| D | mediator measurement error under-attributes the stratosphere | obs, 6,958 days | established |
| E | selection on the outcome manufactures 31–45% of a reported impact | obs | established |
| F | the bias obeys a projection law, β × ΔS | obs R²=0.9963, CMIP6 R²=0.9999 | established |
| G | stratospheric classification is **not** exempt from it | obs + CMIP6 | established |
| H | the response distribution is one shifted population | CMIP6 n=1,517 | established, power-bounded |
| I | forced between-event sd σ_f ≤ 0.27 σ (upper 95%) | CMIP6 n=1,517 | established, 76% placebo correction |
| J | **an SSW adds no out-of-sample predictability: event-specific R² −0.006 pre-onset, +0.015 post-onset; 96% of diagnostic skill occurs with no SSW** | CMIP6 n=1,517 | **established, no correction** |
| K | causal stratosphere→surface effect S = +1.11 σ (sd 0.271; +1.36 incl. ECCC) | SNAPSI, 9 centres | established, 2 NH events; ECCC a Tukey outlier |
| **L** | **the DW/NDW contrast is identical with and without an SSW, in a designed experiment** | **SNAPSI, 9 centres, 61 ensembles** | **established — the decisive test, replicated at 3x scale 2026-09-17** |
| **M** | **an SSW shifts the surface distribution without inflating its variance** | **SNAPSI, 9 centres, 1798 vs 1805 members** | **established causally; ratio 0.955 [0.873, 1.051] after an estimator fix** |
| O | the field's two archetypes (Feb-2018 "propagating", Jan-2019 "not") have the same forced odds, and their observed difference is ordinary member-to-member noise | SNAPSI 9 centres + ERA5 | established for these 2 events; Jan-2019's NDW label flips with the pressure level |
| N | the week-2 100 hPa → surface coupling between members is the same with and without an SSW | SNAPSI, 9 centres, 36 matched ensembles | established; Δz +0.023 [−0.044, +0.090] |
| P | operational forecasts show no significant event-specific skill at SSWs; they underestimate the common shift | ECMWF S2S reforecasts, 12 SSWs 2003–2021 | established for one system; conditional p 0.115 / 0.366, power 0.88–0.90 |
| P′ | confirmed in 9 other systems: no event-specific discrimination (r +0.23 vs +0.52 event-free, p 0.91); outcomes on the downward side of the ensembles (mean PIT 0.40 vs 0.54, p 0.005); shift-size correction not significant (p 0.16, low power) | S2S reforecasts, 9 systems, 17 SSWs 1998–2021 | established (H1/H2 pre-registered; D, H1, H2 calibrated, D power 1.00 against strong skill) |
| Q | knowing the event adds no probabilistic forecast skill beyond the shift (CRPSS +0.008 [−0.003, +0.014]), less than on ordinary days (+0.037) | CMIP6, 1,517 events, leave-one-model-out | established |
| R | an imposed SSW raises the chance of a cold fortnight over northern Eurasia from 0.10 to 0.32 (0.26 without ECCC); the polar-cap shift through the ordinary relation predicts it (0.31; 0.27), residual −0.10σ [−0.36, +0.06]; the label adds nothing there (it does in mid-latitude N America, +0.29σ). Europe has extra cooling, N America an unexplained warming | SNAPSI 36 pairs (9 centres); ERA5 39 events; 10 S2S systems, 17 events | established (preregistered); operational: N Eurasia colder than forecast, rank 0.39 vs 0.53, p 0.014 |
| L′ | paired label test: contrast identical with and without an SSW on the same pairs (−0.03σ [−0.12, +0.05]); threshold model exact in controls | SNAPSI, 25 pairs, 8 centres | established |

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

| lag bin | observations (n=43) | CanESM5 (10 members, CP07 events) |
|---|---|---|
| −30..−16 (**pre**-onset) | **+0.789** [0.34, 1.16] | **+0.809** [0.70, 0.91] (n=633) |
| −15..−1 | +0.714 | +0.757 |
| +0..+14 | +0.784 | +0.847 |
| +15..+29 (**post**-onset) | **+0.824** [0.49, 1.09] | **+0.859** [0.79, 0.93] (n=675) |

A precursor cannot be caused by the event that follows it, so R at a pre-onset
lag cannot be a causal fraction — yet it equals the post-onset value in both
systems, and in CanESM5 the intervals are about ±0.10 wide. Flatness also survives
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
- rate-matched selection bias under a true null: **AO 38%, NAO 31%, PNA 45%** of
  the reported dSSW impact
- decomposition closes numerically: residuals **0.010–0.024 σ** across three
  independent outcomes
- **the classification rate is uncontaminated**: 70% of real events vs 39% of
  pseudo-events satisfy the criterion, **P < 0.0001** — this is the quantity the
  field should report
- projection law `bias(Y) = β(Y-window on A-window) × bias(A)`, β the OLS **slope**:
  obs **R² = 0.9963** (7 stratifiers; slope 1.023, intercept −0.034), CMIP6
  **R² = 0.9999** (1,517 events; slope 0.992, intercept +0.011)
- **stratospheric classification is not exempt**: splitting random winter dates
  with no SSW by lower-stratospheric descent manufactures **+0.795 σ**; real SSWs
  split identically give +0.979 σ (n=29). The residual after β × ΔS is +0.292,
  95% [−0.585, +0.870] — indistinguishable from zero.
- **selection is not always inflationary** — for western-US snowpack the bias runs
  *opposite* to the effect. "Selection can only inflate, so this is a lower bound"
  is unsafe.
- body count on the published criterion (43 events, 30 DW / 13 NDW = 70%):

| outcome | contrast | share from selection | assumption-free null reproduces |
|---|---|---|---|
| AO | −1.782 | 98% | **94%** |
| NAO | −0.543 | 88% | **87%** |
| PNA | +0.222 | 78% | 66% |

- ERA5, in the literature's own fields: Karpechko 1000 hPa **97%** from selection
  (residual −0.023, p=0.664), 97% reproduced by event-free dates; ACP-2026 850 hPa
  **103%** (residual +0.017, p=0.860), 93% reproduced

### H — one shifted population (`FINDING_no_class.md`)

| test | observations (n=43) | CMIP6 (n=1,517) |
|---|---|---|
| dip, pseudo-calibrated | p = 0.959 | p = 0.824 |
| BIC k=1 vs k=2 | p = 0.160 | p = 0.965 |
| pure-shift KS | p = 0.374 | p = 0.346 |
| implied shift | −0.865 σ | −0.476 σ |

Power at n=1,517, mixing fraction 2/3, total variance held fixed: **0.00** at 0.25σ
and 0.50σ, 0.005 at 0.75σ, **0.85 at 1.00σ**, 1.00 at 1.50σ. Nothing below ~1σ
is excluded and the claim must always carry that bound.

**"About two thirds propagate downward" is what one shifted population produces:**
observed pass rate on the criterion's surface conditions (1–2) 69.2% (53.8% with
all three conditions; 39 events in ERA5 coverage); unshifted pseudo-events 31.3%;
pseudo-events displaced by the measured −0.636 σ **75.4% [61.5, 87.2]**.

### I — forced variance (`FINDING_forced_variance_ceiling.md`)
`σ_f² = Var(Y|SSW) − Var(Y|pseudo)`.

| | CMIP6 (n=1,517) | observations (n=43) |
|---|---|---|
| naive σ_f² | +0.0783 [+0.0401, +0.1167] | −0.0236 [−0.3927, +0.3546] |
| **placebo σ_f²** (pre-onset, true value 0) | **+0.0595** [+0.0171, +0.1011] | +0.2933 |
| corrected σ_f² | **+0.0189 [−0.0332, +0.0722]** | −0.3169 [−1.158, +0.464] |
| σ_f upper 95% | **0.269** (corrected) / 0.342 (raw) | 0.596 (raw) / 0.681 (corrected) |

**76% of the raw excess variance is present before onset** — SSWs cluster in
disturbed winters whose variance is already elevated (pre-onset variance ratio
1.12). Without the placebo this was a clean false positive. σ_f² is
indistinguishable from zero and bounded near 0.27 σ.

**Ceiling against the published criterion, matched.** The DW−NDW contrast the
Karpechko criterion produces on the same outcome is 1.8× the largest contrast
the uncorrected σ_f bound allows in observations (|C| 1.782 at q 0.70, ceiling
0.985) and 2.0× in CMIP6 (conditions 1–2, |C| 1.154 at q 0.70, ceiling 0.567). The uncorrected bound is NOT always the wider one: in observations the placebo-corrected upper bound is larger (σ_f 0.681 against 0.596), and against it the observed contrast is **1.6×** the ceiling (1.127); in CMIP6 the uncorrected bound is the wider, so 2.0× stands (`over_wider_bound_ratio`, added 2026-09-30).

### J — predictability: an SSW adds nothing (`within_model_check.py`, `predictability_ceiling.py`, `headline_stress_test.py`, `predictability_wave_driving.py`)

CMIP6, 1,517 events (Charlton & Polvani detection, validated: §5c), 20 members
from 10 models (10 of them CanESM5). Ridge, GroupKFold(5) by member, three
predictor tiers (P1 before onset, 20 features; P2 to onset, 26; P3 adds the
stratosphere over days 0..+30, 38). Response: annular mode over the in-season
days of +8..+52 (18% of events, mostly March onsets, have < 44 of 45 days).

**The numbers to quote — within-member standardisation, against 1,000
size-matched event-free draws** (same members, onset counts and calendar days;
pseudo-onsets > 135 d from every real onset) and a member-cluster bootstrap
(1,000 resamples). Every draw has its own seeded stream.

| tier | real | null mean (sd) | **event-specific** [95% from null] | p(null ≥ real) | cluster bootstrap 95% |
|---|---|---|---|---|---|
| P1 pre-onset | 0.0544 | 0.0601 (0.0138) | **−0.006** [−0.033, +0.020] | 0.64 | [−0.061, +0.047] |
| P2 at-onset | 0.0583 | 0.0720 (0.0144) | **−0.014** [−0.043, +0.013] | 0.82 | [−0.076, +0.045] |
| P3 + post-onset strat. | 0.3815 | 0.3663 (0.0222) | **+0.015** [−0.028, +0.059] | 0.25 | [−0.061, +0.095] |

- **No tier's event-specific share differs from zero.** 96% of the post-onset
  "diagnostic skill" (0.3663 of 0.3815) is reproduced with no SSW present.
- **Against the split's implied share.** q(1−q)C²/Var(Y), computed on the same
  events as the contrast: CMIP6 Karpechko conditions 1–2 **0.528** (q 0.70;
  `predictability_ceiling.json`), observed AO **0.614** (q 0.70;
  `recompute_published_criterion.json`), ERA5 NAM 1000 hPa **0.402** (q 0.54) and
  850 hPa **0.341** (q 0.59) (`era5_recompute_and_two_thirds.json`). The largest
  event-specific upper bound, P3's bootstrap +0.095, is 3.6× below the smallest
  of these and 6.5× below the largest; P1's +0.047 is 7× to 13× below.
- **Imperfect matching errs against the claim**: in 1,022 of 20,000 member-draws
  some pseudo-onsets lack enough in-season outcome days, which lowers the null.
- **Pooled across members** the same tiers score 0.084 / 0.092 / 0.407 (ridge;
  gradient boosting 0.026 / 0.038 / 0.351); about a third of the pooled P1 skill is
  between-model structure (within-member 0.054). The within-member value is the
  one that answers the field's question.
- **Stress tests** (`headline_stress_test.json`): on pseudo-events the pooled
  pipeline retains 75% (P1) and 89% (P3) of its real skill; per-member CV R² does
  not scale with a member's coupling strength (r = +0.016); in-sample optimism
  is +0.020 (P1), so a split's 0.5–0.6 is not ordinary overfitting.
- **Wave driving** (TEM-identity proxies, 86 features): P1 0.084 → 0.097, P2
  +0.001, P3 −0.016, against the same pipeline without them on the same events.
- **Power** (pooled P2 features): injected R² 0.02 detected in 75% of trials,
  0.05 in 100%.
- **Consistency with I**: the placebo-corrected σ_f upper bound 0.269 implies a
  classifier R² ≤ 0.269²/Var(Y) ≈ 0.14 (Var(Y|SSW) 0.526); the pooled
  cross-validated 0.084 sits inside it.

**Observations cannot settle it.** n = 42 (the heat-flux record bounds the
window): vortex only R² = −0.126, vortex + wave −0.112; over 200 equivalent
assignments of winters to folds the 5–95% ranges are −0.21 to −0.04 and −0.36 to
−0.03, positive in 0% and 2.5% of assignments. Skill requires p < 0.05 **and**
R² > 0. With the training-fold mean as baseline instead of the pooled mean, the
two values are −0.022 and −0.009 (methods audit 2026-09-25); neither is
positive. **The observational null is uninformative, not negative, and must
never be cited as evidence against predictability.**

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

**Paired test (added 2026-09-25, `paired_NH` in the same JSON).** On the 25
centre × initialisation pairs where both arms form a contrast (8 centres):
nudged −1.622 vs control −1.594, **paired difference −0.028σ [−0.122, +0.054]**
(10,000 centre-cluster bootstrap). Against the one-population expectation (a cut
at zero on a Gaussian with the ensemble's own mean and sd): control residual
**+0.002 [−0.024, +0.023]** — the threshold model is exact with no SSW; nudged
residual +0.195 [+0.107, +0.302] — the forced contrast is SMALLER than one
population predicts, the opposite of what two populations would give (probably
Karpechko condition 2 or skewness; not tested). Falsifier (|paired difference|
> 0.2σ) not met.

**The contrast is the same whether or not an SSW happened** (−1.622 vs −1.590).
It cannot be reporting a physical difference between groups: in the nudged arm
there is none to report, and in the control arm there is no event. The contrast
is the cut.

**The rate is where the signal lives** — 0.842 with forcing against 0.450 without.
This independently reproduces, in a designed experiment, what
`FINDING_selection_on_outcome.md` §6 found in observations (70% real vs 38%
pseudo, P < 0.0001). *The classification rate is real and uncontaminated; the
selected group's magnitude is not.*

**NRL excluded, and the guard was fixed to catch it.** At three of its four NH
initialisations (s20180208, s20181213, s20190108) its nudged and control members
differ by at most 0.023 Pa (7, 10 and 2 of 80 members differ at all); at
s20180125 by up to 10.4 Pa (mean 3.1 Pa), against a member spread of ~110 Pa. A first
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
| excluding the short-lead init | — | **0.996 [0.895, 1.106]**, KS p = 0.650 |

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
variance is bounded at **+0.05 σ² at 95%** (variance-ratio upper bound 1.051). Two independent routes now agree —
one needing an 84% placebo correction (I), one needing none (M).

> **POWER LIMITATION, stated not buried.** The Gaussian-mixture BIC test in this
> script has **essentially no power**: at a proper 5% critical value it detects a
> 2/3 mixture with probability 0.07 at 0.25 σ separation, **0.00 at 1.0 σ** and
> 0.06 at 1.5 σ, because the within-component spread stays near 1 σ while the
> separation grows. Its p = 0.985 is **not** evidence for one population and must never be quoted as
> such. **The variance bound and the KS test carry this result; the mixture test
> does not.** A first version of the power code compared the statistic against a
> *single* null draw — a coin flip under the null — and produced power
> *decreasing* with separation, which is impossible and is what exposed it.

> **The short-lead initialisation.** `s20190108` begins **6 days after onset**, so
> the +8..+25 window sits at lead 2–19 days, before the ensemble has fully
> diverged. The headline (0.955, KS p = 0.37, 36 ensembles) includes it; without
> it the ratio is 0.996 [0.895, 1.106], KS p = 0.650. Both support a pure shift.

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

### P — operational forecasts: no event-specific skill; the shift is underestimated (`s2s_forecast_test.py`)

Plan approved 2026-09-25 before any data were retrieved. ECMWF S2S reforecasts from
the ECMWF Data Store (model year 2022, one version; hindcasts 2002–2021; 11 members;
all 35 Monday/Thursday starts Dec–Mar; `acquire_s2s_reforecasts.py`), polar-cap
MSLP at 00 UTC (instantaneous; GRIB `stepType=instant`), reduced with the ERA5
`cap()`. 12 primary-catalogue SSWs 2003–2021; ~2,200 zone-free control starts.
Outcome: days +8..+25 NAM proxy; DW = Karpechko conditions 1–2.

| | short (starts 2–9 d before) | mid (10–17 d before) |
|---|---|---|
| r(ensemble mean, observed) across events | +0.765 (null mean +0.563) | +0.466 (+0.385) |
| **conditional p** (primary) | **0.115** | **0.366** |
| Brier skill of P(DW); conditional p | +0.294; 0.069 | +0.005; 0.629 |
| forecast P(DW), SSW vs control | 0.60 vs 0.35 | 0.48 vs 0.37 |
| observed DW rate, SSW vs control | 0.75 vs 0.29 | 0.75 vs 0.30 |
| ensemble-mean vs observed A (Pa), SSW | −166 vs −269 | −25 vs −269 |
| mean PIT at SSW starts (0.5 = calibrated) | 0.39 | 0.30 |
| outer-bin share (0.167 expected); binomial p | 0.04; 0.11 | 0.14; 1.00 |

- **Discrimination:** not significantly better at SSWs than at event-free dates of
  the same calendar period and lead. Short-lead point estimate is above the null;
  "no strong event-specific skill", not "zero". Consistent with Nebel et al. 2024
  (their SSW skill at the 91st percentile of random dates).
- **No second population:** the rank histogram at SSW starts is not U-shaped (the
  pre-stated falsifier). The deviation is a one-sided offset: outcomes more
  downward than the members.
- **The shift is underestimated:** the −166 / −269 Pa above are anomalies against
  a climatology that includes the SSW years, not shifts. Relative to ordinary
  (control) starts the forecast shift is **−270 Pa against −423 Pa observed (64%)**
  at short lead and **−113 against −428 Pa (26%)** at mid lead
  (`shift_forecast`, `shift_observed`, added 2026-09-26). The forecast raises
  P(DW) 0.35 → 0.60 while reality goes 0.29 → 0.75.

**Validation** (`validate_s2s_null.py`, 40 synthetic datasets per case on the real
layout; binomial s.e. 0.034): conditional-r rejection rate 0.03/0.05 under no
skill and 0.05/0.05 under equal skill everywhere; power against skill present only
before SSWs 0.88/0.90. Brier test: power 0.28/0.15, rejection 0.10 under no skill at
short lead: secondary only. Two design changes made during validation, before the
real forecasts were analysed: the first null (moving each event to another year)
had only 8 event-free winters and rejected 20% under noise in the first check; the
unconditional correlation test was conservative (mean p 0.67 under equal skill),
biased toward our thesis, so the null is conditioned on the across-event spread of
the outcome.

**Limits:** 12 events; one model version; starts of one event share an outcome, so
start-level reliability tests are liberal (the event is the unit); spread/error
1.36 (short) indicates an over-dispersed ensemble at SSW starts.

### P′ — nine more systems: no event-specific skill; outcomes lie low in the ensembles (`s2s_multimodel_test.py`)

Plan approved 2026-09-25 (item A); H1 and H2 written into the script before any
non-ECMWF forecast was retrieved (D was run as in P but is not among the registered
hypotheses). ECMWF is the discovery set; the confirmatory set is nine of the other
eleven S2S systems in the ECMWF Data Store (BoM: only a 2014 version; UKMO, IAP-CAS:
no msl), one model version each (`acquire_s2s_reforecasts.py`; members per start as
archived; KMA has Jan–Mar starts only):

| system | version (model date) | hindcasts | members | starts used |
|---|---|---|---|---|
| ECCC | on the fly, 2025 | 2001–2020 | 4 | 700 |
| CMA | on the fly, 2022 | 2007–2021 | 4 | 525 |
| HMCR | on the fly, 2025 | 1991–2020 | 11 | 510 |
| KMA | on the fly, 2026 | 1993–2016 | 7 | 288 |
| CNRM | fixed, 2025-06-01 | 1999–2024 | 11 | 416 |
| JMA | fixed, 2022-09-30 | 1991–2020 | 5 | 240 |
| CNR-ISAC | fixed, 2023-10-16 | 2001–2020 | 8 | 480 |
| NCEP | fixed, 2011-03-01 | 1999–2010 | 4 | 516 |
| CPTEC | fixed, 2023-01-04 | 1999–2018 | 11 | 123 |

Starts 2–9 d before onset; each system's event value is the mean over its starts;
the multi-model value is the mean over the systems covering the event (2–8 per
event). Null: 10,000 draws moving every event to a zone-free date within ±21
calendar days (any hindcast year) at which at least half of the event's systems
have starts; D conditioned on the across-event variance of the outcome (factor
1.25).

| | D: r (event-free mean) | D p (cond.) | H1: mean PIT (event-free mean) | H1 p | H2: CRPS gain (event-free mean) | H2 p |
|---|---|---|---|---|---|---|
| **multi-model, 9 systems, 17 SSWs** | **+0.23 (+0.52)** | **0.91** | **0.398 (0.544)** | **0.005** | +121 (−39) | 0.16 |
| incl. ECMWF | +0.24 (+0.53) | 0.91 | 0.399 (0.545) | 0.004 | +106 (−40) | 0.17 |

Per system (secondary): mean PIT below its event-free mean in 9 of 10 systems
(not NCEP, whose null itself sits at 0.36, 1999–2010 hindcasts, 4 members);
individual H1 p < 0.05 for ECMWF 0.035, CMA 0.001, HMCR 0.004, KMA 0.005,
CNR-ISAC 0.020 (CNRM 0.050). No system shows discrimination significantly above
its event-free null (smallest conditional p: ECMWF 0.11). Ensemble-mean vs
observed A at SSWs, multi-model: −116 vs −241 Pa (anomalies against a
climatology that includes SSW years; not a shift, see P).

- **Discrimination** of event-to-event differences after SSWs is not better than
  on event-free dates; the multi-model point estimate is below it.
- **H1 confirmed:** outcomes after SSWs fall on the downward side of the
  ensembles, more than on event-free dates — the systems under-forecast the
  common shift. (The rank-histogram / second-population test exists only for
  ECMWF, in P; it was not run for the nine systems.)
- **H2 not confirmed**; calibration gives it little power (0.11 under damping), so
  this is absence of evidence.
- ECMWF here (event-weighted PIT 0.401, p_cond 0.110) differs slightly from P
  (start-weighted 0.391, 0.115): different weighting, min-other-years rule and
  random stream; same data.

**Validation** (`validate_s2s_multimodel.py`, 100 synthetic datasets per case on
every system's real layout; binomial s.e. 0.022): equal skill everywhere — D 0.06,
H1 0.03, H2 0.03; no information — D 0.07, H1 1.00 (correctly: a forecast that
knows nothing misses the shift), H2 0.00 (mean p 0.19: low power); damped 0.6×
truth — H1 power 0.66, H2 0.11; skill only at starts 2–9 d before SSWs (another
year's outcome elsewhere; added after the real result, because D's multi-model power
had not been measured) — D power 1.00. H1c (sensitivity below): 0.02 under equal
skill, 1.00 no information, 0.55 damped. One
design change, found on synthetic data before the real result: the first null
required ALL covering systems at a pseudo-onset and had no candidate for 4 events.

**Sensitivities (added 2026-09-26, after the primary result, from review):**
| | D p (cond.) | H1 p | H1c p (system-centred PIT) | H2 p |
|---|---|---|---|---|
| confirmatory (primary) | 0.91 | 0.005 | 0.012 | 0.16 |
| excluding CPTEC | 0.90 | 0.005 | 0.009 | 0.15 |
| incl. ECMWF | 0.91 | 0.004 | 0.009 | 0.17 |

H1c removes each system's mean PIT over its candidate pseudo-onsets before
averaging: a null draw averages only the systems present (quorum), so a system with
an unusual PIT level (NCEP 0.36) otherwise enters events and null unequally, which
favours H1. CPTEC's forecast climatology rests on 1–3 other years per start date.

**Limits:** 17 events, many shared across systems (not independent tests);
per-system versions differ in age (NCEP 2011); starts of one event share an outcome;
ERA5 starts 1998-11, so the observed climatology omits HMCR/JMA 1991–98 and KMA
1993–98 hindcast years (a constant offset per start date and lead, applied alike
to events and null; not quantified here).

### Q — forecast value: knowing the event adds no probabilistic skill (`forecast_value_test.py`)

Plan approved 2026-09-25. CMIP6, 1,517 events, 20 members, 10 models. Gaussian
forecasts of the days +8..+52 annular-mode response (member sigma units, not
demeaned): SHIFT (calendar regression on other models' SSWs) vs EVENT-AWARE
(ridge on stratospheric predictors + calendar; out-of-fold residual spread).
Leave-one-MODEL-out, closed-form CRPS, 2,000 model-bootstrap replicates; the same
pipeline on 200 size-matched sets of event-free dates.

| predictors | CRPSS after SSWs [95%] | CRPSS ordinary days [2.5–97.5%] | SSW-specific | p(pseudo ≥ SSW) |
|---|---|---|---|---|
| P2 up to onset | **+0.0076 [−0.0028, +0.0142]** | +0.0374 [+0.0234, +0.0530] | −0.030 | 1.000 |
| P1 before onset | +0.0065 [−0.0028, +0.0139] | +0.0314 [+0.0185, +0.0437] | −0.025 | 1.000 |

Event-aware information improves on the shifted distribution by <1% of CRPS
after SSWs (interval includes 0; upper bound 1.4%), and is worth ~4x more on
ordinary winter days; in none of 200 matched samples was it worth as little as
after SSWs. Falsifier (SSW-specific value > 0 at p < 0.05) not met — the data
point the other way. Interpretation (not tested): after an SSW the stratosphere's
information is mostly the shift.

### R — regional cold risk: the shift, acting through the ordinary circulation–temperature relation (`snapsi_regional_test.py`, `era5_regional_test.py`)

Plan approved 2026-09-26, written into both scripts before any regional
temperature was analysed. Regions from the literature: NEURASIA 50–65N 10–130E
(Kretschmer et al. 2018; PRIMARY), HI_EUROPE 55–70N 0–60E, MID_EASIA 35–55N
90–150E, MID_NAMER 35–55N 120–60W (Huang et al. 2021). Window: daily means of
post-onset days +8..+24; T in sigma of the control (SNAPSI) or event-free (ERA5)
distribution. Mediation test: fit T on the member's/date's polar-cap NAM WITHOUT an
SSW, predict with the SSW, residual R. Label test: does the downward indicator add
to N? SNAPSI tas: `acquire_snapsi_tas.py` (3,603 members, 9 centres, 2.5° daily);
ERA5: `acquire_era5_t2m_regions.py` (WB2, 1959–Jan 2023; 39 of 43 events).

**SNAPSI (causal; 36 centre × initialisation pairs, ~450 members per arm and init;
centre-cluster bootstrap):**

| region | SSW effect S_T (σ) | residual beyond NAM, R (σ) | share explained | P(cold fortnight): control / nudged / predicted from shift | label coef (nudged) |
|---|---|---|---|---|---|
| **NEURASIA** | **−0.85 [−1.40, −0.49]** | **−0.10 [−0.36, +0.06]** | 88% | 0.10 / 0.32 / 0.31 | −0.04 [−0.19, +0.09] |
| HI_EUROPE | −0.87 [−1.22, −0.63] | −0.27 [−0.47, −0.14] | 69% | 0.10 / 0.32 / 0.26 | +0.14 [−0.03, +0.36] |
| MID_EASIA | −0.51 [−0.80, −0.27] | −0.21 [−0.50, +0.03] | 60% | 0.10 / 0.27 / 0.20 | −0.10 [−0.23, +0.03] |
| MID_NAMER | +0.49 [+0.37, +0.63] | +0.80 [+0.63, +1.04] | — | 0.10 / 0.08 / 0.18 | +0.29 [+0.06, +0.53] |

Excluding ECCC (post-hoc, as in K): NEURASIA S_T −0.60, R +0.02 [−0.02, +0.07];
HI_EUROPE R −0.19 [−0.26, −0.12]; MID_NAMER R +0.69 [+0.61, +0.76].
Primary falsifier (R interval excludes 0 AND |R| > 0.2σ in NEURASIA): **not met**.
Self-test: the 8-centre bootstrap interval covered 0 in 180/200 null draws (90%,
nominal 95%: slightly liberal) and detected a −0.5σ direct effect in 200/200.

**ERA5 (39 events; calendar-window null of event-free dates):** NEURASIA −1.00 K,
R −0.33 K (−0.12 sd), p 0.38, 67% explained; HI_EUROPE −1.04 K, R −0.42 K, p 0.27;
MID_EASIA −0.19 K, R +0.15 K, p 0.53; MID_NAMER +0.04 K, R +0.55 K (+0.27 sd),
p 0.09. The null's 95% range is about ±0.7 K: residuals below ~0.27 sd are not
detectable here. Label coefficients uninformative (intervals ±1.7 K).

- Over northern Eurasia an imposed SSW raises the chance of a cold fortnight
  from 0.10 to 0.32 (0.26 without ECCC), and the shift of the polar-cap
  circulation, acting through the ordinary SSW-free relation, predicts it (0.31;
  0.27 without ECCC). The downward label adds nothing there; in mid-latitude
  North America it does (+0.29σ [+0.06, +0.53]), where the polar-cap index does
  not describe the response.
- Regionally the polar-cap index is not the whole story: extra cooling over
  high-latitude Europe (−0.19 to −0.27σ) and a warming of mid-latitude North
  America that the NAM relation does not predict (+0.7 to +0.8σ; same sign in
  ERA5). The forecast product after an SSW should therefore be the shifted
  distribution of the regional variable itself, not a label or the NAM alone.
**Operational (`s2s_regional_test.py`; S2S daily-mean 2 m temperature of the ten
systems of P′, `acquire_s2s_reforecasts.py --var t2m`; 17 SSWs; P′ machinery; H1-T
registered as primary before the data were analysed):**

| region | mean rank (event-free mean) | H1-T p | system-centred p | leave-one-out max p | ens-mean / obs anomaly |
|---|---|---|---|---|---|
| **NEURASIA** | **0.390 (0.527)** | **0.014** | 0.021 | 0.036 | −0.3 / −1.3 K |
| HI_EUROPE | 0.426 (0.529) | 0.069 | 0.083 | 0.134 | −0.8 / −1.3 K |
| MID_EASIA | 0.539 (0.497) | 0.74 | 0.79 | 0.88 | +0.1 / +0.2 K |
| MID_NAMER | 0.484 (0.532) | 0.22 | 0.27 | 0.34 | +0.3 / +0.2 K |

Holm over the four regions: NEURASIA 0.054. Over NEURASIA all ten systems lie below
their event-free value. Discrimination D is not interpretable here: its calibration
rejected 10% under equal skill (below). Anomalies are against climatologies that
include the SSW years (not shifts).
**Calibration** (`validate_s2s_multimodel.py --var t2m`, 100 datasets per case, real
t2m layouts and ERA5 NEURASIA outcome, noise scaled to the pressure case): equal
skill — D 0.10 (liberal), H1 0.03, H2 0.02, H1c 0.01; no information — D 0.05, H1
1.00; damped 0.6× — H1 0.32 (low power), H2 0.03 (noise scaled by the anomaly-s.d.
ratio after the 2026-09-29 review; the raw-s.d. scaling had given 0.23); skill only before SSWs — D 1.00.
- The operational systems under-forecast the northern-Eurasian cold after SSWs as
  they under-forecast the polar-cap shift (P′), with a test of low power.

### S — revision robustness, 2026-09-29 (referee concerns; plan approved 2026-09-29)

All post hoc unless marked; every script's design is in its docstring, and the
two marked "registered" were committed (57019f6) before their data existed.

**Operational rank deficit (`s2s_baseline_test.py`, `s2s_multimodel_robustness.py`;
nine confirmatory systems, 17 events):**

| check | polar-cap NAM | N-Eurasian T |
|---|---|---|
| mean rank after SSWs, 95% event bootstrap | 0.40 [0.29, 0.52] | 0.39 [0.29, 0.50] |
| all-winter baseline (SSW winters in, SSW windows out): null mean, p | 0.53, p 0.012 | 0.51, p 0.028 |
| vortex-state regression (NCEP 60-90N u10): slope per m/s; SSW residual, p | 0.0006; −0.10, p 0.028 | 0.0006; −0.09, p 0.057 |
| mixed model, (1 given date) + system: SSW coefficient [95%] | −0.14 [−0.25, −0.04], p 0.009 | −0.13 [−0.25, −0.01], p 0.031 |
| events covered by ≥5 systems (15): rank vs null, p | 0.41 vs 0.54, p 0.012 | 0.41 vs 0.53, p 0.039 |
| P(outcome<0), forecast vs observed after SSWs | 0.57 vs 0.82 (reliability 0.079) | 0.53 vs 0.59 (0.019) |
| same on event-free dates | 0.41 vs 0.29 (0.014) | 0.41 vs 0.42 (0.003) |
| discrimination Δr [95%]; MDE (α 0.05, power 0.8); power at Δr 0.2/0.3 | −0.29 [−0.76, +0.22]; 0.43; 0.42/0.63 | −0.17 [−0.82, +0.30]; 0.33; 0.53/0.76 |

The deficit is SSW-specific: rank barely depends on vortex state across winter
dates and SSWs fall below that relation. The absolute test against 0.5 is marginal
(interval reaches 0.52); the systems sit slightly above 0.5 on ordinary dates.

**SNAPSI contrast shortfall explained (`snapsi_contrast_diagnosis.py`, 25 NH pairs):**
Gaussian-expectation residual +0.195 [0.108, 0.304] (as L); condition 1 only +0.090
[0.047, 0.145]; EMPIRICAL shifted-control null (control members shifted by the
imposed effect, full criterion) +0.040 [−0.026, +0.097] (23 pairs); skewness
difference +0.02 [−0.12, +0.16]; s.d. ratio 0.98 [0.94, 1.03]. One shifted
population reproduces the nudged contrast. **SH variance ratio 1.74 [1.41, 2.15]**
(8 ensembles, 402+402 members): the austral minor warming widens the distribution
as well as shifting it — an exception to the NH translation.

**Event differences bounded (`event_difference_bounds.py`, from I):** CMIP6 forced
between-event sd 0.137 (upper 0.269) in member units, shift −0.48, noise sd 0.67:
forced P(negative outcome) across events 0.62–0.87 (0.47–0.93 at the upper bound);
share of a single event's label variance due to forced differences **2.2% (≤7.8%)**.

**Observed cold anchor (`era5_regional_test.py`, added):** cold northern-Eurasian
fortnight after 39 SSWs **0.26 [0.13, 0.41]** vs 0.10 event-free (HI_EUROPE same;
MID_EASIA 0.21 [0.08, 0.33]; MID_NAMER 0.08 [0.00, 0.18]).

**Vortex geometry (`vortex_geometry_test.py`, registered 57019f6; area weights and
validation period corrected to Seviour et al. 2013 before the first complete run):**
NCEP Z10 moments; the paper's event rule gives 17 displaced / 13 split events in
their 52 winters (paper, ERA: 17 / 18 — NCEP's 2.5° grid smooths splits).
Catalogued SSWs, registered any-day rule: 27 split, 10 displaced (skewed to split);
7-day rule: 9 / 10. Split minus displaced, ERA5 NAM1000 days 8–25: −0.25σ [−0.61,
+0.13], p 0.28 (7-day: −0.40, p 0.12); days 8–52: −0.02, p 0.92; N-Eurasian T
−0.42 K, p 0.70; DW rate 0.78 vs 0.60, Fisher p 0.41 (7-day 0.89 vs 0.50, p 0.14).
Same direction as the literature (splits stronger early), not significant at 37
events.

**CMIP6 forecast value, robustness (`forecast_value_robustness.py`; 1,514 events
scored in every variant after the precursor-coverage fix; 200 draws per null):**

| set / tier | CRPSS after SSWs [95%] | balanced | no CanESM5 | ordinary dates (p) | weak-vortex dates (p) |
|---|---|---|---|---|---|
| base P1 | +0.65% [−0.33, 1.34] | 0.44% | 0.49% | 3.26% (1.0) | 1.13% (0.82) |
| base P2 | +0.76% [−0.24, 1.38] | 0.36% | 0.62% | 3.91% (1.0) | 1.62% (0.94) |
| base P3 (post-onset) | 17.2% [10.8, 21.4] | 14.7% | 16.8% | 20.7% (1.0) | 18.9% (0.92) |
| plus P1 | +0.62% | 0.51% | 0.95% | 4.38% (1.0) | 2.41% (0.99) |
| plus P2 | +0.76% [−0.26, 1.92] | 0.52% | 1.14% | 4.83% (1.0) | 2.93% (1.0) |
| plus P3 | 17.6% | 15.2% | 16.0% | 21.2% (0.99) | 19.1% (0.87) |

"plus" = base + wind tendencies (10/50/100 hPa, 60N and cap), 10-100 hPa shear,
75N-45N geometry and the pre-onset surface annular mode. Predictor s.d. at SSWs /
nulls: 1.02 (ordinary), 1.03 (weak). Event information is worth less after SSWs
than on ordinary or weak-vortex days in every variant; richer predictors help the
comparison days, not SSWs.

**Missed warming or weak response? (`s2s_postonset_test.py`; registered 57019f6;
the matched null (a201458) was adopted AFTER a code review had emulated the
polar-cap hit/miss test with the all-starts null on leads 10+ (misses p 0.0095);
the change raised p; nine systems):**

| starts | polar-cap NAM rank vs null, p | N-Eurasian T rank vs null, p |
|---|---|---|
| 2-9 d before onset, all (17 events) | 0.40 vs 0.55, p 0.004 | 0.39 vs 0.52, p 0.022 |
| — HITS: >=50% of members reverse u10(60N) within +-3 d of onset (125 starts, 17 events) | 0.41 vs 0.55 (matched), p 0.014 | 0.40 vs 0.52, p 0.059 |
| — MISSES (43 starts, 15 events) | 0.41 vs 0.55, p 0.029 | 0.41 vs 0.53, p 0.11 |
| — hit minus miss (15 events) | +0.003 [−0.047, +0.053] | −0.002 [−0.121, +0.109] |
| 0-7 d AFTER onset (observed SSW in initial state; 17 events) | 0.50 vs 0.56, p 0.14 | 0.46 vs 0.52, p 0.18 |

Not a missed warming: forecasts that predict the reversal are as far off as those
that miss it. With the observed SSW in the initial state the deficit shrinks and is
not significant (cf. Dai et al. 2025). The reversal criterion does not measure the
strength or persistence of the predicted warming, so "forecast warming too weak"
and "coupling too weak" cannot be separated here. Neither pre-stated reading
matched exactly (T2-hit low but T1 not significant); reported as found.

### T — second revision round, 2026-09-30 (plans approved 2026-09-30; designs committed in `4c8b3f2` before runs/data)

**T-a Criterion sweep (`criterion_sweep.py`, 9_literature).** 108 versions of the
criterion in ERA5 (39 events; window end 25/35/52 d, 1000/850 hPa, mean below
0/−0.25/−0.5σ, >50/60/70% of days negative, 150 hPa condition off/on) against
event-free dates displaced by each version's own measured shift (400 sets).
Conditions 1–2 (primary, 54 versions): max|z| 1.15, **p 0.965**; mean z −0.04,
p(excess) 0.50; **no version outside its shifted-null 95% interval**; event-free
dates reproduce **76–117%** of the class contrast (median 94%). With the 150 hPa
condition (54): p 0.87 / 0.48; 89–154%. Published version: 53.8% observed against
54.1% [38.5, 69.3] shifted (share of contrast reproduced 95.7% in this run's 400
sets). CMIP6 (1,517 events, 27 versions, conditions 1–2, 20 sets): p 0.50 / 0.85.
**Registered reference too lenient (code review):** the events' shift is measured
on them, so under the null the test never rejects (0/200 simulated). **Calibrated
(post hoc, 200 replays of the whole procedure): surface versions p 0.845 (max|z|),
0.96 (mean z); with condition 3, 0.745 / 0.995. Power against a planted two-class
structure (2/3 at −0.6σ, 1/3 at +1.2σ beyond the shift): 0.815 (max|z|), 0.51 (mean
z); with condition 3, 0.37 / 0.10.** Per-version 95% intervals are uninformative
(0.2% of versions outside under the null); the Fig. 1b interval has the same
construction (a consistency check, not a sized test). CMIP6 arm not recalibrated.
**No version of the criterion separates SSWs from one shifted population.**

**T-b Observations vs CMIP6 (`obs_model_compatibility.py`, 6_predictability).**
Observed statistics (43 events, CPC AO days 8–52): KS 0.133, variance ratio 0.96,
shift −0.81 s.d., conditions-1–2 rate 0.72, contrast share 0.95. Against 10,000
43-event draws from all CMIP6 events: two-sided p 0.22 (KS), 0.40 (variance
ratio), 0.59 (shift), 1.00 (rate), 0.74 (share); **combined p 0.62: the observed
record is a typical CMIP6 draw on every statistic** (compatibility, which is weaker
than exchangeability). **Per-model draws are INVALID as registered** for models
with 53–152 events: drawing 43 of 59 events without replacement leaves almost no
spread, so p collapses toward 0 as an artefact (FAILURES 2026-09-30); only CanESM5
(683 events: combined p 0.64) is interpretable. Model shifts genuinely differ
(full-population −1.02 to +0.10 s.d.; NorESM2-LM, 53 events, is the one positive); the observed −0.81 is inside that range.
**Perfect model** (each model as truth, 1,000 draws of 43): the pure-shift KS test
rejects in **0.0–0.1%** of draws (nominal 5%) even where the model's full
population has variance ratio 1.27 — **at n = 43 the KS test has almost no size or
power; the observed KS p = 0.37 must not be cited as support.** The two-thirds
test's error rate is 0.0–1.4% (conservative); the variance-ratio interval covers
the model's value in 91–100% of draws. Registered pass rule failed in every model
only through PM1/PM2 being below 1% (conservative, not liberal).

**T-c Forced-variance identification audit (`forced_variance_ceiling.py` C4).**
Earlier numbers reproduced exactly (separate seeded stream). C4a (1,000 synthetic
data sets per cell, real pseudo-event pairs): with the assumption true, the
placebo-corrected upper bound covers the true s_f in 0.945–0.975; the uncorrected
bound in 0.956–0.973 (CMIP6) but **0.896–0.898 in observations, just under the
pre-set 0.90** (slightly liberal; the falsifier uses the wider, corrected bound).
Forcing tied to the pre-onset state (corr(pre, post) 0.07 obs, 0.12 CMIP6) barely
matters; damping (λ² 0.873) and forcing anti-correlated with the concurrent
internal anomaly destroy coverage, as expected. C4b, damping at the SNAPSI lower
bound λ² = 0.873: s_f ≤ 0.78 (obs) / 0.42 (CMIP6); **contrast/ceiling 1.38 (obs),
1.67 (CMIP6)** — still above 1; forced share of one label's variance **22% (obs),
18% (CMIP6)** (below the pre-set 25%); forced P ranges 0.26–0.99 / 0.30–0.98.
C4c tipping: through the pre-onset state unreachable (ρ* −4.1 obs, −3.2 CMIP6);
**through the concurrent internal anomaly reachable at ρ* −0.30 (obs) / −0.40
(CMIP6), −0.26 / −0.36 together with the SNAPSI-bounded damping** (the
concurrent case and the combined case are not in the registered design; the
combined case was added after review) — a dependence that neither observations nor these models can test
(within a SNAPSI ensemble the forcing is identical, so it cannot vary). Reading:
the falsifier survives the audit but is weaker than stated: forced differences can
produce at most ~60–72% of the contrast, and up to ~20% of a label, under the worst
damping SNAPSI allows; it fails only if forcing anti-correlates with concurrent
internal variability at |ρ| ≥ 0.3.

**T-d Shift vs continuous vs two regimes (`class_model_comparison.py`, 6_predictability).**
CMIP6, 1,517 events, leave-one-model-out, predictors up to onset. CRPS: M0 shift
0.4093, M1 continuous (mean and log-spread linear in predictors) 0.4059, **M2 two
regimes with state-dependent membership 0.4047**, M2b state-blind two populations
0.4092. **Registered test: M2 beats M1 after SSWs by more than on ordinary days —
G = +0.0012 against 200 event-free sets' −0.0023 [−0.0053, +0.0003], p 0.01
(ignorance p 0.01). By the pre-set reading this SUPPORTS more than a pure shift.**
M2b vs M0: p 0.18 (no evidence of two populations irrespective of state). The
registered test's own controls: planted regimes (2:1, 1 s.d. apart) detected in
29–33%; **no-regime synthetic data exceed the event-free 95th percentile in 13–15%
(miscalibrated)**. POST HOC (added after seeing the result, labelled): a skewed
one-population model (M1 + constant skew-normal shape) barely improves on M1
(0.00015) and M2 still beats it beyond ordinary days (p 0.02 CRPS, 0.005 ign);
calibrated against the negative controls themselves, p 0.02 (CRPS) / 0.00 (ign).
**What M2 finds after SSWs: components 0.49 s.d. apart in mean with spreads 0.58
and 0.72 — CLOSER in mean than the same fit finds on ordinary days (0.71 [0.58,
0.84]), different mainly in spread**; membership s.d. 0.32 (ordinary 0.26 [0.22,
0.31]); response skew −0.25 (ordinary −0.16 [−0.26, −0.06]). Gain 0.3% of CRPS.
Reading: the response after an SSW has a state-dependent SHAPE (spread) that a
linear log-spread model misses; it is not a separation into two kinds of mean
response, and its forecast value is negligible. SNAPSI arm: mixture gain lower
with the SSW imposed (−0.029 [−0.043, −0.014]) but the planted-mixture positive
control was detected in 0% — the arm has no power and is not evidence either way.

**T-f Held-out 2023–24 SSWs (`s2s_heldout_test.py`, 6_predictability; registered in
`4c8b3f2` before any of these forecasts or post-Jan-2023 ERA5 were retrieved).**
Events 2023-02-16, 2024-01-16, 2024-03-04 (primary catalogue); systems CNRM (on disk,
never scored at these dates), ECMWF model year 2025 (2005–2024), CMA model year 2025
(2010–2024), all three covering every event. ERA5 extended with ARCO-ERA5 (overlap
correlation 1.000, offset −0.01 Pa / ≤0.04 K). **HO1 polar cap: mean rank 0.72
against null 0.53 [0.29, 0.77], p 0.93 — ABOVE the null mean: by the registered
reading a FAILURE TO REPLICATE** (per event 0.68, 0.83, 0.64; every system above 0.5
in every event; observed outcomes on the upward, positive-NAM side: A_obs +166,
+446, −35 Pa). **HO2 N-Eurasian T: 0.57 vs 0.51, p 0.64** (observed warmer than
forecast in 2023 and Jan 2024). HO3 pooled (secondary, 20 events): polar cap 0.445
vs 0.541, p 0.032; N-Eurasia 0.414 vs 0.525, p 0.031; the 17 P' events alone with
the spliced observations 0.397 (p 0.004) and 0.387 (p 0.013). Reading: the
registered under-forecast is not seen in the three later events (different model
versions for ECMWF/CMA; Lee et al. 2025 describe the 2024 events as non-canonical
and weakly coupled; March 2024 was mildly downward); with 3
events power is low, but the direction is the opposite, so the operational claim
must be stated as holding for 1998–2021 and weakened overall.

**T-f diagnosis (`s2s_heldout_diagnosis.py`; EXPLORATORY, POST HOC, requested
2026-09-30).** D1 no defect: the held-out code path reproduces CNRM's P' rank
exactly (0.4271). **D2 not the model versions: the 2025 versions also under-forecast
the 1998–2021 events** — ECMWF 2025 0.416 vs null 0.526 (10 events, p 0.08; ECMWF
2022 on the same events 0.397); CMA 2025 0.24 vs 0.54 (5 events, p 0.002); the
four 2025–26 versions already in P' (ECCC, HMCR, KMA, CNRM) were all below their
null. D2b paired: ECMWF CY49R1 minus CY47R3 on the same 10 events +0.019 [−0.055,
+0.076]. "Model year" is the archive's production label: ECMWF changed cycle
(CY47R3 → CY49R1), CMA's configuration matches the same BCC-CPS-S2Sv2 in both sets
(research notes, model_versions.md; version date not retained in our files). **D4 the rank follows the realised outcome:** events with negative-NAM
outcomes average rank 0.31 (n 14), positive 0.75 (n 6), in both periods; 76% of the
1998–2021 outcomes were negative, 1 of 3 held-out; three events drawn from the
earlier 17 average ≥ 0.72 with probability 0.007 (the earlier sample was
outcome-lucky relative to later events). **D5 all three held-out events lie inside
the 95% prediction interval of the 1998–2021 relation** between observed and
ensemble-mean response (slope 0.38 [−0.57, 1.49]; residual s.d. 395 Pa). **D3 winter
effect:** ranks on ordinary dates move together across systems from winter to winter
(~0.2 in 2013 to ~0.7; winter-mean correlations between systems 0.55–0.87);
2022/23 was low (0.25–0.54) and 2023/24 high on only 14 ordinary dates, so a
2023/24 winter-level error is not established. D4b: given the outcome sign the
periods agree — upward outcomes rank 0.76 (4 events, 1998–2021) and 0.75 (2 held
out); downward 0.29 (13) and 0.64 (1). **D6 (added
after D1–D5):** the forecasts gave the held-out negative outcomes 0.64, 0.55, 0.85;
only the third (March 2024) was negative; the three signs are 3.8 times
likelier under the forecasts as issued than under the 1998–2021 bias (logit shift
+1.02) — moderate evidence that the earlier bias was overestimated. Reading: the
failure is chiefly outcome sampling (a different mix of upward and downward
outcomes), not a model change; the under-forecast claim must be reported as uncertain in size and
persistence. ED Fig. 5.

**T-g — the 4 March 2026 SSW (`s2s_heldout2026_test.py`; registered `a511fbc`, AO disclosure `a04a4a8`).** Real-time ECDS forecasts, starts 23 Feb–2 Mar with a same-month-day reforecast: ECMWF 4 starts × 101 members, ECCC 3 × 21, HMCR 1 × 41, KMA 2 × 8, NCEP 3 × 16 (primary); JMA 1 × 5 (secondary); CMA, CNRM, CNR-ISAC have no matching reforecast month-day and are not used. Observed outcome upward: polar-cap A_obs +156 Pa (positive NAM), N-Eurasia +3.0 K against the ECMWF hindcast mean. **Primary polar cap: mean rank 0.67 vs null 0.61 [0.28, 0.93], p 0.60 — above the null mean: by the registered reading not consistent with the 1998–2021 under-forecast; one event, reported as one more event.** Per system 0.62/0.66/0.32/0.78/0.95. With JMA 0.65, p 0.60. **Four held-out events pooled: 0.70 vs 0.55 [0.35, 0.75], p 0.93** (2023 0.68, Jan 2024 0.83, Mar 2024 0.64, 2026 0.67). **Secondary N-Eurasian temperature: 0.93 vs 0.56 [0.15, 0.91], p 0.98** (every system 0.89–0.98: much warmer than forecast); with JMA 0.93, p 0.98; four events pooled 0.66 vs 0.53, p 0.86. Reading: like the 2023–24 events, an upward outcome with a high rank (D4: upward outcomes ranked 0.76 in 1998–2021 too); the operational under-forecast is a 1998–2021 property not seen in any of four later events.

### U — the SSW itself is a threshold on a continuum (`vortex_threshold_continuity.py`, 5_mechanism; registered `2a0c40b` before any outcome was related to the wind)

Units: 114 NCEP vortex-weakening episodes 1958–2025 (42 reversals; rule reproduces all
39 NCEP compendium onsets), running variable = minimum u(10 hPa, 60N), cutoff 0.
**Dose (primary for size): the surface response scales smoothly with how far the
wind falls** — per 10 m/s weaker minimum: AO −0.36 [−0.51, −0.21] (sign: slope +0.36
in X), NAM1000 0.22 [0.15, 0.30], N-Eurasian T 0.36 K [0.07, 0.66] colder, and the
probability of a "downward" outcome +0.076 [0.009, 0.142] higher. **Global jump at
reversal (continuity): none** — AO −0.28 [−0.99, +0.47] (90%: [−0.87, +0.35]), NAM1000
−0.18 [−0.55, +0.18], T +0.09 [−1.87, +2.01], downward label +0.12 [−0.21, +0.43]; Holm
1.0 for all. Local-linear RD (rdrobust) and local randomisation at 0 also null
(AO +0.29 [−0.77, +1.76]; permutation p 0.65/0.62), BUT the local fits are unstable:
4 of 16 placebo-cutoff tests have p < 0.05 (largest jumps 3–6 units at ±5–10 m/s),
so by the pre-set falsification rule the LOCAL RD is uninterpretable; conclusions rest
on the global estimates. Placebo outcome (pre-deceleration AO) and balance (day of
year, era) pass; deceleration size borderline (p 0.07). **CMIP6 (5,726 episodes, 1,483
reversals): slope 0.173 [0.159, 0.187] σ per 10 m/s; jump +0.05 [−0.02, +0.13] σ
(90% upper 0.12)** — the reversal adds at most a small fraction of the SSW effect.
Calibration: the observational global-jump test rejects 2.5% of 68-winter model blocks
(no planted jump) and detects a planted −0.8σ jump (the whole-SSW shift) in 78%.
**Reading: "SSW" is itself a threshold on a continuum of vortex weakenings, as the
downward label is a threshold on a continuum of outcomes; shown in observations
without a model, and in CMIP6 with power.**

### V — the shift rule as a forecast (`shift_rule_forecast.py`, 6_predictability; registered `bf453bd`)

Model probabilities only (SNAPSI, 36 pairs): P(N-Eurasian days 8–24 below the no-SSW
q-quantile) = 0.24 [0.15, 0.37] (q 0.05), 0.32 [0.22, 0.47] (0.10), 0.45 [0.37, 0.57]
(0.20); CMIP6 P(negative annular mode, days 8–52) 0.74. **Verified untuned on 39
observed SSWs: 0.21 (8), 0.26 (10), 0.38 (15) — consistent with the rule at every
threshold (binomial p 0.71, 0.40, 0.42) and rejecting climatology (p 0.0006, 0.004,
0.008)**; Brier skill vs climatology +0.12, +0.09, +0.11; relative economic value
at cost/loss = q: 0.75 [0.23, 0.87], 0.68 [0.24, 0.84], 0.60 [0.17, 0.79]. **The value is positive only for cost/loss ratios between q and the observed frequency; between the observed frequency and the rule probability it is negative (the rule slightly over-states the risk, e.g. q 0.10: −0.36 at 0.32); zero elsewhere** (ED Fig. 7b). Do not quote 0.60–0.75 without that range. Risk
ratios rise with severity as a shift implies (rule ×4.8/3.2/2.3; observed
×4.1/2.6/1.9). With the 2023–24 events (42): same. Negative AO: rule 0.74, observed
0.70 (30/43), climatology 0.42; p 0.49 vs rule, 0.0003 vs climatology; BSS +0.26.
Only the predictor is out of sample (observed frequencies at q 0.10 had been seen).
**Prospective 2026 (run after T-g): N-Eurasian days 8–24 +4.7 K, not cold at any q** (rule gave 0.24/0.32/0.45, so 'not cold' had probability 0.76/0.68/0.55); one event, uninformative; CPC AO window not complete in the file (ends 31 March 2026). Retrospective V1/V2/AO numbers identical on re-run (parsed comparison).

### W — where the forecasts lose the shift (`s2s_leadlag_test.py`, 6_predictability; registered `4c8b3f2` before any 100 hPa forecast was retrieved)

17 P' events (1998–2021), starts 2–9 d before onset, eight confirmatory systems
(CNRM excluded: its ECDS "100 hPa" fields are 10.8–11.5 km, not 16 km — deviation
recorded before the first run), days +8..+25, P' calendar-window null. **L1 the
lower-stratospheric anomaly is under-forecast: observed 100 hPa polar-cap height
ranks 0.446 in the ensembles against null 0.595 [0.485, 0.701], p 0.005 (Holm
0.01).** **L2 given the forecast 100 hPa state the surface is not detectably
under-forecast: 0.447 vs 0.477 [0.371, 0.586], p 0.30.** L4 (secondary) the 10 hPa
wind is not under-forecast: 0.505 vs 0.435, p 0.89. ECMWF (discovery, 12 events)
the same: L1 0.420 vs 0.582 p 0.008; L2 p 0.30; L4 p 0.89. L3 (descriptive) mean
surface error −131 Pa [−312, +59], of which −67 [−167, +34] through the 100 hPa
error and −64 [−202, +82] remainder (within-ensemble slope 3.5 Pa/gpm) — roughly
half, intervals wide. **Reading by the pre-set rule (L1 low, L2 ~ null): the shift
is lost in the lower stratosphere — a persistence deficit below a correctly forecast
10 hPa warming (cf. Garfinkel et al. 2025), not a surface-coupling deficit.** Scope:
the 1998–2021 events, whose surface under-forecast did not replicate in 2023–24
(T-f); L2 not rejecting is absence of evidence for a coupling deficit, not proof of
calibration.

**V revision 3 (registered `add7ba9`; `shift_rule_forecast.py --revision3`, JSON key `revision3`):**
R1 strict independence (37 events, without 2018-02-12 and 2019-01-02): 7, 9, 14 cold
fortnights; p under rule 0.57/0.38/0.41, under climatology 0.002/0.009/0.012 — HOLDS.
R2 severity N-Eurasia: q 0.025/0.05/0.10/0.20/0.33 → observed 0.10/0.21/0.26/0.38/0.51
vs SNAPSI 0.18/0.24/0.32/0.45/0.60 (p rule 0.30–0.71; p clim 0.016/0.0006/0.004/0.008/
0.025); RR observed 4.1/4.1/2.6/1.9/1.55. The rule is ABOVE the observed frequency at
every q (slight over-prediction). R3 regions: HI_EUROPE moderate cold consistent and
beats climatology (q 0.10–0.33, p clim ≤ 0.004) but EXTREMES over-predicted (q 0.025:
0/39 vs 0.20, p < 0.001; q 0.05: 3/39 vs 0.24, p 0.014); MID_EASIA consistent with
rule, not with a detectable change vs climatology; MID_NAMER rule ≈ climatology, observed
consistent with both. R4 weeks: days 15–28 (SNAPSI 22 pairs, 9 centres, only the later
initialisations of the two events): N-Eurasia observed 0.08/0.15/0.28 vs rule
0.20/0.27/0.42 (p 0.07–0.11), NOT distinguishable from climatology; days 29–42 (12
pairs): consistent with the rule; HI_EUROPE beats climatology (p 0.004–0.045).
R5 second predictor (CMIP6 days 8–24 NAM, mean −0.49σ, through the ERA5 event-free
T–NAM regression, 2.0 K per σ): p 0.11/0.19/0.33; observed 0.21/0.26/0.38 lies
BETWEEN the two independent model routes (p under R5 0.07/0.31/0.50). **R6 OUT OF SAMPLE (1940–58, outcomes never examined): FAILED AT q 0.10.**
ERA5 (CDS) detector finds 12 onsets 1941-02-07 … 1958-02-01 (the last = the catalogue's
1958-01-30, outcome not examined before); on 1958–2024 it recovers all 43 catalogued onsets
within 3 d (+1 extra, 2002-02-17). N-Eurasia days 8–24 below the period's event-free
q-quantile: 0/12 (q 0.05; rule 0.24, p 0.08), 0/12 (q 0.10; rule 0.32, p 0.012 —
REJECTED), 2/12 (q 0.20; rule 0.45, p 0.08); vs climatology p ≥ 0.62 (no detectable
rise). From 1946: 0/9, 0/9 (p 0.036), 2/9. Negative NAM proxy days 8–52: 7/12 (0.58) vs
CMIP6 0.74 (p 0.20) and period climatology 0.42 (p 0.38). Pooled with V1 (51): 8, 10, 17 —
p under rule 0.19 / 0.052 / 0.092; vs climatology 0.004 / 0.03 / 0.02.
Exploratory (post hoc): the early reversals are weaker (median minimum u −3.8 vs −8.5 m/s
for 1959–2022; 58% vs 26% above −4 m/s); mean N-Eurasian anomaly −0.2 K vs −1.0 K, NAM
proxy shift −0.36σ; but the dose slope (~0.37 K per 10 m/s) explains only ~0.2 K of the
0.8 K gap — weak events do NOT explain the failure. Early ERA5 stratosphere weakly
constrained (few radiosondes reached 10 hPa); 12 events. Report as an out-of-sample
failure of a single per-SSW probability.

**U revision 3 — ERA5 running variable 1940–2025** (`--era5`, JSON key `era5_1940_2025`;
registered `add7ba9`): 144 episodes (54 reversals), 111 with outcomes; NCEP vs ERA5 u_min
r 0.98 on 96 common units, 5 switch side of 0 (1958-12-04, 1968-01-07, 1977-01-11,
1981-03-04, 2002-02-17). Dose per 10 m/s: NAM proxy 0.16σ [0.10, 0.22], AO 0.29 [0.15,
0.43] (from 1950), T 0.38 K [0.14, 0.62], downward label −0.12 [−0.16, −0.06] (all p <
0.003). Global step: NAM proxy −0.28 [−0.67, +0.11], AO −0.73 [−1.51, −0.06] (raw p 0.028),
T −0.66 K [−2.15, +0.77], label +0.26 [−0.06, +0.59]; **Holm 0.34 / 0.11 / 0.39 / 0.34 →
no step by the registered rule**; local RD null (AO −0.02 [−1.40, +1.58]) but 5 placebo
cutoffs p < 0.05 → local estimates uninterpretable again; balance on deceleration size
FAILS (p 0.002: reversals have bigger drops). From 1946: the same. The AO global step is
the largest hint of a step anywhere in the project; CMIP6 (5,726 episodes) bounds it at
+0.05 [−0.02, +0.13]σ.

### Y — within-start ensemble comparison (`s2s_member_experiment.py`, 5_mechanism; registered `be804a7`)

10 S2S systems, ~4,500 starts; 461 mixed starts (≥2 reversing, ≥2 not); primary window
(anchor+8..+25, k* ≤ 9) only 40 starts → N1 shift −0.05σ [−0.29, +0.17] (NOT significant);
N2 variance ratio 0.97 [0.73, 1.36]; N3 DW rate 0.54 vs 0.45, contrast −0.875 vs −0.872
(diff −0.004 [−0.33, +0.30]); secondary window (+8..+20, 162 starts): shift −0.15σ [−0.24,
−0.03], variance ratio 0.96 [0.82, 1.18]. N4 (11,481 members, start fixed effects): dose
−0.05σ/10 m/s [−0.10, +0.01], step −0.05σ [−0.19, +0.10] — NO within-ensemble dose or step
(cf. SNAPSI between-member coupling weak). N5 balance −0.01σ (clean). N6 temperature shift
−0.08σ [−0.34, +0.13]. Quantile shift lower tail −0.58σ vs upper +0.12σ. Reading: translation
and unchanged contrast replicate in a third model set; the shift is small and imprecise.

### Z1 — SNAPSI residual quantiles (`snapsi_residual_quantiles.py`, 8_experiment; registered `113c7d3`, tolerance |R| ≤ 0.25σ, |E| ≤ 0.05 declared before run)
Pooled NH NAM (36 pairs): ALL intervals inside tolerance → APPROXIMATE TRANSLATION (R from
+0.10 at q05 to −0.11 at q95; E within [−0.03, +0.02]). Feb-2018 unresolved (R(0.05) +0.24
[+0.02, +0.38], compressed lower tail); Jan-2019 unresolved. N-Eurasian T (days 8–24):
unresolved (R(0.05) −0.11 [−0.45, +0.13]; E(0.05) +0.018 [+0.002, +0.046]).

### Z2 — transport cross-validation of the shifted null (`criterion_transport_cv.py`, 9_literature; registered `bbfeb9c`; RETROSPECTIVE)
Whole-winter folds, 39 events: published surface DW 27 observed vs 29.2 expected (p 0.52);
with 150 hPa 21 vs 21.3; 0/108 versions p < 0.05. Brier constant-rate minus shifted −0.003
[−0.024, +0.015] (label adds nothing beyond the shift). CRPS days 8–52 NAM: climatology minus
shifted +0.19 [+0.07, +0.30]; mean PIT 0.50; outer deciles 26%.

### Z3 — ICON 18 spin-off ensembles, summary level (`icon_event_heterogeneity.py`; registered `6ad6e8a` before values read)
Between-event variance of ensemble-mean polar-cap response (NAM sign, ~1000 hPa, daily-sd
units), days 8–25: 0.41–0.43 (≈80% of within-ensemble daily variance); days 26–42: 0.10–0.12
(13–15%). Corr with week-2 100 hPa anomaly r 0.80 [0.63, 0.91]. Includes tropospheric
initial-state memory and selected events → NOT comparable to the CMIP6 forced variance
(0.019, net of pre-onset); but RETRACT "events differ little": realised responses differ
substantially. Shape/labels not assessable (no member data).

### Z4 — issue-time state forecast (`state_forecast.py`, 6_predictability; registered `4f5c53a`; RETROSPECTIVE whole-winter CV)
E1, 57 ERA5 SSWs 1940–2026 (incl. a one-day ERA5-only reversal 2025-11-28 not in NCEP): mean CRPS
clim 1.65, SNAPSI rule 1.66, CMIP6 rule 1.55, STATE 1.49, STATE+SSW 1.51. State vs clim
0.0942 [-0.0328, 0.2092]; vs SNAPSI rule 0.103 [-0.0157, 0.2117]; positive
both eras & LOWO, but CI includes 0 → DECISION: forecasting = LIMITATION. SSW indicator adds
nothing (-0.0081 [-0.0184, 0.0028]); oracle min-wind adds nothing. Cold periods 13 obs vs
expected 5.7/18.1/12.2/11.9. E2 (17 events, starts 0–7 d after onset): raw ens 1.18, calibrated
1.09, state 1.58, clim 1.98 (event-free) — dynamical ensembles far better; calibration +8% [5, 11].
State model and climatology are REISSUED at each start date (ERA5 state to that date), same valid
window; every one of the 198 starts is scored by all four methods (checked 2026-10-02, `n_scored_per_method`).

### CODE-REVIEW CORRECTIONS (2026-10-02; fresh-context review of the revision-3 scripts; all re-run)
- **Y/N4** (`s2s_member_experiment.py`): fixed effects were keyed on start DATE across systems
  (1,450 "starts"); keyed on (system, start) = 2,526 groups → dose **+0.12σ/10 m/s [−0.02, +0.27]**
  (was −0.05), step **+0.04 [−0.16, +0.30]** (was −0.05). N1–N3, N5, N6 unchanged. 4,498 starts total.
- **Z3 ICON**: the review's "6/18 have no onset in the record" was WRONG and is reverted (same day):
  Loeffel et al. 2026 (Sect. 2) start every spin-off at onset or the day before, so t = 0 is the onset
  for the 6 already-easterly ensembles. Registered 18-event result stands (Z3 above); the 12-crossing
  subset (V_b 0.59–0.60, r 0.92) is a sensitivity only (JSON key `sensitivity_crossing_only`).
- **V-R5** (`shift_rule_forecast.py`): CMIP6 index in all-year s.d. units vs ERA5 Nov–Mar event-free;
  threshold from a subset of windows. Fixed: **p_rule2 0.099/0.184/0.323** (was 0.112/0.194/0.327),
  CMIP6 mean −0.42σ; q 0.05 binomial p **0.052** (borderline). Member Nov–Mar s.d. 1.03–1.27.
- **Z4** (`state_forecast.py`): E2 climatology not event-free → fixed: clim CRPS **1.98** (was 1.83),
  state vs clim **0.20 [−0.02, 0.38]**, calibrated vs clim 0.45; F2 unit fix: F2 CRPS 1.546, cold
  expected 12.2. Decision unchanged.
- **U-ERA5** (`vortex_threshold_continuity.py`): covariate lacked Sep–Oct → 32/144 units dropped.
  Autumn pressure added: n = 143 (AO 124). Dose replicates (AO 0.31 [0.18, 0.43]; T 0.37 K [0.16,
  0.59]; label +0.12 [0.07, 0.16]); **AO step −0.54 [−1.23, +0.09], raw p 0.10, Holm 0.40 — the
  earlier "hint of a step" (−0.73, p 0.03) was partly the dropped sample.** Added fixed-bandwidth E3
  and doy/era balance (pass); drop balance still fails (p 0.002).
- **Z2** transport: registered Brier vs unshifted climatology added: +0.145 [+0.001, +0.273]
  (+0.151 [+0.014, +0.296] with 150 hPa).
- **Z1** residual quantiles: per-centre summaries added: Meteo-France NAM upper tail compressed
  (q95 −0.42σ [−0.75, −0.04]); ECCC temperature widened; others unresolved.

### X — the state-dependent shape is a continuous dose (`response_shape_dose.py`, 6_predictability; registered `7d8ed6f`)

S3 fitted M2 (P2, all SSWs): component means −0.66 / −0.16, sd 0.72 / 0.58 — the
DOWNWARD component is the WIDER one (weight 0.62); one mode at every weight
(Stephenson et al. 2004: "unimodal but not multinormal" = single regime). S4
quantile shift SSW minus event-free: −0.40 to −0.63σ across q 0.05–0.95; demeaned,
lower tail stretched (−0.15σ [−0.18, −0.015] at q 0.05), rest within intervals.
**S1 (primary) with the realised 100 hPa wind days 0–30 as predictor: G −0.0146
vs event-free −0.0109 [−0.0160, −0.0062], p 0.92 (ign 0.78); S2 all post-onset
levels p 0.84. The pre-set rule's criterion "accounted for by the dose" is met,
BUT the test has NO POWER (below) → reported as INCONCLUSIVE (underpowered).** Continuous-model skill 13.8% after SSWs / 16.1% event-free; M2 10.3% /
13.2%. S5 the dose itself is unimodal (ΔBIC +25 for two components; fitted
mixture one mode). Registered controls INVALID here (both rates 1.0: noise
outcomes lack the dose relation); post-hoc power check (synthetic y = real in-sample ridge mean incl. dose, R² 0.32, + noise s.d. 0.60; planted 2:1 regimes 1σ apart): detection **0.08** vs 0.05 nominal → underpowered by the registered <0.2 rule; real G at the 7th percentile of the no-regime synthetic sets. Do NOT claim the dose explains the shape; claim only unimodality, wider downward component, unimodal dose.

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
   (P(null ≥ 0.85) = 0.0005), though the missing baseline is real and large: with no
   event at all, window overlap alone gives r = **+0.489** [+0.107, +0.766].

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
   forecast busts. **Full text read 2026-09-25** (open access, CC BY; local copy
   in archive/, git-ignored). Findings we use: (i) Zcap100 (days 0–14) from
   initialisations 2–9 d before onset explains ~2/3 of the inter-event spread,
   Zcap500 only ~25%; (ii) the SSW-date skill for Zcap500 over days 5–35 (r = 0.52,
   13 hindcast SSWs, predictor model Zcap100) sits at the **91st percentile** of
   10,000 sets of 13 random Dec–Mar dates (mean r 0.15) — *not* significant at 95%,
   in their words indistinguishable from sampling variability: an operational
   analogue of J; (iii) Zcap500 response significantly stronger after **split**
   SSWs (and after Kodera-absorbing events) in models and ERA5 — a non-zonal event
   property that neither SNAPSI (zonal-mean nudging) nor our zonal-mean CMIP6
   predictors can test; stated as a limit in the manuscript.

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

13. **Spaeth, Rupp, Garny & Birner (2024)**, Commun. Earth Environ. 5, 126,
   doi 10.1038/s43247-024-01292-z (abstract read via OpenAlex) — S2S ensemble spread
   in near-surface geopotential height is *reduced* after weak polar vortex states
   over the Norwegian Sea and Scandinavia, by up to 20%. **Scope of M:** M's
   "shape unchanged" is for the POLAR-CAP MEAN; regional spread can change. Say so.
14. **Roberts & Vitart**, QJRMS (arXiv 2411.17694, abstract read) — ECMWF
   subseasonal wintertime NAO, days 16–45, 2001–2020: RPC ≈ 1.5 (signal-to-noise
   paradox, possible unreliability); well calibrated over 1959–2023. **Any
   forecast-reliability test must check calibration on non-SSW start dates.**

**What remains novel:** the quantification (β×ΔS, R²=0.99–1.00 in two systems); the
extension to stratospheric stratifiers; the power-bounded distributional test; the
shift predicting "two thirds"; and above all **the out-of-sample predictability
ceiling and the pre/post tier contrast**, which nobody has measured.

**Never lead with "unrecognised."** The field has circled this for 15+ years.

---

## 5. Open gaps

| gap | status |
|---|---|
| forecast/operational test | answered for ECMWF (result P); other S2S systems not tested |
| true v'T' in CMIP6 | proxy only; needs hundreds of GB of daily 3-D va/ta |
| observational power | not quantified: at n=42 the CV R² moves by up to 0.3 with the fold assignment (§3 J); n=46 is the whole record |
| SNAPSI NH events | 2. Not fixable. σ_f not estimable there. |
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
- **CMIP6 event detection is validated** (`validate_ssw_detector.json`, 2026-09-25).
  `detect_ssw` implements Charlton & Polvani (2007) as stated by Butler et al.
  (2017): first easterly day at 10 hPa 60N in Nov–Mar, 20 consecutive westerly
  days between events, final warmings (no 10 consecutive westerly days before
  30 April) excluded. On NCEP R1 1958–2024 it reproduces all 39 NCEP-NCAR
  compendium dates to the day, with no extra events in the compendium's
  coverage; after it, 2024-01-17 and 2024-03-04. The detector used until then
  (any easterly day, then a 20-day skip) gave 1,888 CMIP6 events; every CMIP6
  result was re-run on the 1,517 CP07 events, all callers detecting on the full
  daily series.
- **Pseudo-events are zone-free** (2026-09-25): onsets > 135 days from every real
  onset (`zone_free_index`), so no pseudo-event window overlaps a real event's
  −60..+75-day zone. Cleaning only the onset day had left 3.6% of observational
  draws (and 5–9% of CMIP6 windows) contaminated. All affected producers re-run.
- **Implied R² is computed where the contrast is**, from the same outcome and
  events (`implied_r2_eta2` in the two literature JSONs, `implied` in
  `predictability_ceiling.json`); no contrast is typed into another script.
- **January 2019's label survives a true daily mean**
  (`era5_nam_daily_mean_check.json`): +0.025 σ (±0.003) against +0.019 σ from
  00 UTC; still NDW at 1000 hPa over days 8–52.

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
| P | `03_data_ingestion/acquire_s2s_reforecasts.py`, `07_physical_decomposition/s2s_forecast_test.py` | `s2s_forecast_test.json` |
| P′ | `07_physical_decomposition/s2s_multimodel_test.py`, `06_simulation_validation/validate_s2s_multimodel.py` | `s2s_multimodel_test.json`, `validate_s2s_multimodel.json` |
| Q | `07_physical_decomposition/forecast_value_test.py` | `forecast_value_test.json` |

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

Against that: the headline is model-based, the observational arm is uninformative at n=42, and
the framing is corrective over prior art that includes White et al. (2019) saying
the differences are "entirely there by construction."

**Realistic best targets:** npj Climate and Atmospheric Science or Communications
Earth & Environment; JGR-Atmospheres / J. Climate / WCD as strong specialist homes.
This is not a downgrade — it is where this field's readers are.

**What would change the answer:** resolving whether the correction bites on any
operational or forecast-skill claim (§5). If an S2S skill claim is contaminated,
the significance changes category.
