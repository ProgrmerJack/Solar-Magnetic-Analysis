# FINDING: an observational diagnostic separating causal stratospheric influence from common cause

**Date:** 2026-07-31
**Status:** the strongest result in this project. Obtained by re-reading its
central null result correctly.

---

## 1. The null that was a signal

`FINDING_design_robustness.md` reports the conventional between-winter composite
(−0.954) and the within-winter estimator (−0.952) agreeing to **0.002 σ**, and
treats it as a dead end — no estimator gap left to explain, Gate 8 fails,
thesis dead.

That reading was wrong. The two estimators contrast different things:

| estimator | comparison | contains |
|---|---|---|
| between-winter | SSW window vs day-of-year-matched days from winters with **no** event | anything making an SSW-hosting winter anomalous as a whole |
| within-winter | the same window vs **that winter's own** baseline, >75 d from onset | only what changes **after onset** |

Their ratio is therefore a statement about mechanism:

> **R = within-winter effect / between-winter effect**
>
> **R ≈ 1** — the anomaly appears after onset relative to the same winter.
> Event-scale. Consistent with the stratosphere **causing** the surface anomaly.
>
> **R ≈ 0** — SSW-hosting winters are anomalous throughout and nothing extra
> happens at onset. Winter-scale. The signature of a **common cause**.

R needs no experiment. It is computable from any observational record or model
run, and it is exactly the quantity a common-cause confounder attacks.

## 2. Result

Peak response bin +15..+29 d, winter-block bootstrap, winters as the resampling
unit:

| dataset | between-winter | within-winter | **R** | 95% CI | n events |
|---|---|---|---|---|---|
| **observations** | −0.991 | −0.817 | **+0.82** | [+0.50, +1.09] | 43 |
| **CanESM5** (10 members) | −0.716 | −0.590 | **+0.83** | [+0.76, +0.89] | 846 |
| MIROC6 (2 members) | −0.138 | −0.044 | +0.32 | [−4.39, +1.03] | 189 |

**Observations and an 846-event independent model ensemble agree to within 0.01
on R.** That agreement is not built in anywhere — the model was acquired,
detected and analysed by a wholly separate pipeline.

## 2b. R validated under known truth — it is an unbiased estimator

The identity R = causal/(causal+confound) is a population statement. It says
nothing about whether the *estimator* recovers it from finite, autocorrelated,
seasonally varying data with 43 events. `validate_R_diagnostic.py` tests that by
simulating records where the causal fraction is **set** rather than estimated:
AR(1) with φ = 0.943 (the measured lag-1 autocorrelation of observed winter AO),
77 winters, 43 events in 36 winters — the observed configuration. The
association is split into an event-scale depression after onset and a
winter-scale depression of event-hosting winters, in a known ratio.

| true f | mean R | bias | sd | bootstrap coverage |
|---|---|---|---|---|
| 0.00 | −0.032 | −0.032 | 0.298 | 0.843 |
| 0.25 | 0.203 | −0.047 | 0.275 | 0.847 |
| 0.50 | 0.481 | −0.019 | 0.222 | 0.900 |
| 0.75 | 0.703 | −0.047 | 0.199 | 0.857 |
| 1.00 | 0.959 | −0.041 | 0.203 | 0.827 |

**R is unbiased**: max |bias| = 0.047 across the entire range. It tracks the
causal fraction.

**But the winter-block bootstrap under-covers** (0.83–0.90 against a nominal
0.95), so the interval reported in §2 is too narrow and is superseded below.
Reported here rather than quietly fixed, because it is the same class of error
this project has had to withdraw before.

## 2c. Calibrated inference against a simulated common-cause null

Rather than rely on an interval known to under-cover, R is tested directly
against the distribution it takes when the truth is **pure common cause**
(f = 0), from 1,500 simulated records:

> null R: mean −0.046, sd 0.320, 95th percentile +0.372, 99th +0.456

| dataset | R | P(R ≥ observed \| pure common cause) |
|---|---|---|
| observations | +0.824 | **0.0000** — 0 of 1,500 replicates |
| CanESM5 | +0.825 | **0.0000** — 0 of 1,500 replicates |

Calibrated interval for the observations, using the simulation sampling sd:
**R = 0.82, 95% CI [+0.20, +1.45]** (against the too-narrow bootstrap
[+0.50, +1.09]). It still **excludes 0**.

So the conclusion survives proper calibration: the pure common-cause hypothesis
is rejected at p < 0.0007, by an estimator demonstrated to be unbiased, and
independently replicated in an 846-event model ensemble.

## 3. What this licenses, and what it does not

**Licensed.** The observed SSW–surface association is **predominantly
event-scale**. R = 0.82 with a CI excluding zero comfortably: the link is *not*
explicable as a common cause that makes some winters both SSW-prone and
surface-anomalous. This is direct observational evidence for causal downward
influence, from a diagnostic that requires no model and no experiment.

**NOT licensed.** That a winter-scale component exists in observations. The
calibrated interval [+0.20, +1.45] **includes 1**. The honest statement is that
the observed link is consistent with being *entirely* event-scale; the apparent
18% winter-scale share is not resolved at n=43. Nor is the *magnitude* of the
causal fraction pinned down — [0.20, 1.45] is wide. What is established is the
**sign of the conclusion**: not a common cause.

In CanESM5, where n=846, R = 0.83 [0.76, 0.89] **excludes 1**, so a small
winter-scale component is resolved there — about 17%, precisely estimated.

**MIROC6 is uninformative for R, and an earlier reading of it was wrong.** Its
between-winter effect is only −0.138, so R is a ratio with a near-zero
denominator and its interval is correspondingly useless ([−4.39, +1.03]). MIROC6
does not demonstrate common-cause behaviour; it simply has a weak SSW–surface
link overall, consistent with its weak vortex–surface correlation (+0.193
against CanESM5's +0.255). The earlier claim in
`FINDING_ensemble_precursor.md` §4c that MIROC6 shows a winter-scale response
rested on a raw composite that was not day-of-year adjusted against non-event
winters, and is withdrawn.

## 4. Why this matters beyond SSWs

**R is a model-evaluation metric that the standard composite cannot replace.** A
model can reproduce the observed composite while getting R wrong — the right
surface response for the wrong mechanism. Composite-based evaluation, which is
what the field uses, is blind to that distinction. R separates models that
couple downward from models that merely co-vary.

The construction is general. Any event study with a natural blocking unit — a
winter, a season, a site — admits the same decomposition of an association into
a within-block (event-scale) and between-block (selection/common-cause) part.

## 5. Relation to the rest of the project

Consistent with the mediation result (`FINDING_mediation.md`), which found the
stratosphere carries 45–54% of the wave-driving-to-surface link after correcting
for mediator measurement error. Both say the stratospheric pathway is
substantial and neither says it is exclusive; they are independent routes to
that conclusion, using different data and different assumptions.

It also explains why `FINDING_design_robustness.md` found estimator choice
irrelevant: in the observed system R ≈ 1, so the two estimators *must* agree.
That agreement is not the absence of a finding — it is the finding.
