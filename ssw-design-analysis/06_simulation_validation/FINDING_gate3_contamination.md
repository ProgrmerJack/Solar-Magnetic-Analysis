> **⚠ PRE-RE-FREEZE DOCUMENT.** Written before 2026-07-30, when the event
> catalogue was re-frozen from **39 events / 33 winters** to **43 / 36** after the
> consensus rule was found to be era-dependent
> (see `../02_event_catalogues/FREEZE_RECORD.md`). Every count below is the old
> build. The JSONs in `results/` were all regenerated; this document was not.
> Qualitative conclusions are mostly unaffected, but **no number here should be
> quoted without checking it against `results/`.**

# Gate 3, check 1: the requested exclusion introduced its own bias

Scripts: `finalise_gate3.py` (flawed), `gate3_clean_null.py` (corrected).

## What the three draws gave

Frozen primary catalogue (39 events / 33 winters), 5,000 draws, tested bin +15..+29:

| pseudo-onset draw | null mean | coverage | FPR |
|---|---|---|---|
| `full_pool` — any winter | **−0.006** [−0.015, +0.002] | — | — |
| `excluded` — ≥60 d before / ≥75 d after any real onset | **+0.234** [+0.225, +0.243] | 0.878 | 0.122 |
| `nonssw_winters` — only SSW-free winters | **+0.149** [+0.140, +0.159] | — | — |

The three were required to agree. They do not, and the strictest draw is biased
**positive** by about a quarter of a standard deviation. Identical for 6 harmonics and the
12-knot spline (+0.2337 vs +0.2339), so it is not a seasonal-model artefact.

## Cause: the baseline, not the onsets

In an event study the baseline is "days far from any onset" — and under the pseudo-onset
draw, *onset* means the **fake** onset. Pushing fake onsets away from real events therefore
leaves the real, genuinely depressed SSW periods sitting **inside the baseline**. The fake
bins are then measured against a reference contaminated with the real negative response, so
they read positive. **+0.234 is the mirror image of the real effect leaking into the
reference category.**

This is the same class of error as the original single-window flaw that started this
rebuild: the contaminated group was the comparison group, not the treatment group. It is
worth recording that the error recurred in a test specifically designed to detect it.

## Correction

Remove real-event influence from the **data** before drawing anything, so neither the fake
bins nor the fake baseline can contain it:

1. drop every day within [−60, +75] d of any **real** onset — 4,516 of 13,864 winter days
   (33%), leaving **9,348 clean days across 74 winters**;
2. draw pseudo-onsets among the surviving days, preserving the observed day-of-year set;
3. fit on the cleaned series.

Every day in the calibration, treated and control alike, is then free of real SSW
influence, so a non-zero result can only come from the estimator or the seasonal model.

Sanity check on the corrected design: **+0.012** (30 draws), against +0.234 under the
flawed one. Full run with 4,000 draws and Clopper–Pearson intervals in progress.

## Consequence for the gate

Gate 3 stays **provisionally cleared, not closed**. The `full_pool` result (−0.006) was
correct but was not by itself sufficient evidence, because the check meant to confirm it was
itself biased. Closure now rests on the clean-null run.

---

## Corrected clean-null result: Gate 3 does NOT pass

4,000 draws on the cleaned series (9,348 days, 74 winters), Clopper–Pearson intervals:

| spec | null mean | coverage (95% CI) | FPR (95% CI) | pass |
|---|---|---|---|---|
| 6 harmonics | +0.046 [0.034, 0.057] | 0.922 (0.902, 0.940) | 0.077 (0.060, 0.098) | **no** |
| cyclic spline, 12 knots | +0.063 [0.052, 0.075] | 0.919 (0.898, 0.937) | 0.081 (0.063, 0.102) | **no** |
| 3 harmonics | +0.033 [0.021, 0.044] | 0.929 (0.909, 0.946) | 0.071 (0.054, 0.091) | **no** |

The diagnosis of the previous failure is confirmed directly: **43% of baseline days under
the old `excluded` draw lay inside real-event influence**, which is what produced the
+0.234 bias.

Two problems remain on the corrected design, both smaller than before but real:

1. **A small positive bias survives.** +0.03 to +0.06 σ, with intervals excluding zero.
   Under the pass criterion (|mean| ≤ 0.05) harm3 and harm6 are inside, the spline is not —
   but none of them is centred.
2. **Under-coverage.** Intervals cover 0.92 rather than 0.95, and the false-positive rate is
   0.07–0.08 rather than 0.05. The winter-block bootstrap intervals are too narrow, so the
   estimator is mildly anti-conservative.

Note the ordering: the *least* flexible seasonal model (3 harmonics) has the smallest bias
and the best coverage, and the most flexible (12-knot spline) the worst. That is the
opposite of what residual-seasonality-as-cause would predict, and points away from the
seasonal model as the explanation.

### Open question being tested

Whether the under-coverage is a property of the estimator or an artefact of using only 150
bootstrap replicates for a percentile interval. At 150 replicates, coverage is 0.920
(CI 0.889–0.945). A 1,000-replicate run is in progress; if coverage moves to ~0.95 the
cause is bootstrap resolution and the fix is simply more replicates. If it does not, the
interval construction itself needs changing (BCa, or a studentised interval).

**Gate 3 remains open.** It must not be recorded as passed on the strength of the
`full_pool` number alone.

---

## Resolved: the under-coverage was bootstrap resolution, not the estimator

Same data, same draws, same specification (6 harmonics) — only the number of bootstrap
replicates changed:

| replicates | coverage (95% CI) | FPR |
|---|---|---|
| 150 | 0.920 (0.889, 0.945) | 0.080 |
| **1000** | **0.947 (0.915, 0.969)** | **0.053** |

At 1,000 replicates both criteria are met: the coverage interval contains 0.95 and the
false-positive rate contains 0.05. A percentile interval built from 150 replicates has
2.5th/97.5th percentiles estimated from roughly 4 order statistics in each tail, which is
far too coarse; the resulting intervals are systematically too narrow.

**Standing requirement added: no reported interval may use fewer than 1,000 winter-block
bootstrap replicates.** Every earlier number in this project used 150–800, and any of them
that was near a significance boundary should be treated as unreliable until recomputed.

## Residual bias

The small positive null mean is unaffected by replicate count (it comes from point
estimates, not intervals): +0.033 (3 harmonics), +0.046 (6 harmonics), +0.063 (12-knot
spline), against a criterion of |mean| ≤ 0.05. Three and six harmonics satisfy it; the
spline does not.

The ordering is informative and runs opposite to the residual-seasonality hypothesis: the
**least** flexible seasonal model has the smallest bias and the best coverage. Adding
seasonal flexibility makes calibration worse, not better. On this evidence the earlier
recommendation of 6 harmonics over 3 was not justified — **3 harmonics is the better
primary candidate**, pending its own 1,000-replicate confirmation.

The residual +0.03 σ is about 3% of the response being measured (≈ −0.95 σ at +15..+29 d).
It is small but it is not zero, and it must be disclosed rather than rounded away.

## Gate 3 status

Coverage and false-positive criteria: **met** at ≥1,000 replicates (confirmed for 6
harmonics; 3-harmonic confirmation running).
Null-centring criterion: **met** for 3 and 6 harmonics, **not met** for the 12-knot spline.
Spatial outcomes: **not addressed** — a shared seasonal curve is inadmissible there and
requires separate calibration.
