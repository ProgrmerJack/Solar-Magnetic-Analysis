# Gate 8 — and a result that overturns the paper's framing

> **Refreshed 2026-09-17 from `design_sensitivity.json` on the 43-event catalogue.**
> This document previously carried the 39-event build behind a warning banner. Every
> number below is now read from the artifact. **No conclusion changed** — R² moved
> 0.210 → 0.2167 and the permutation P 0.30 → 0.2808, both still nothing. The event
> subset is still 30, but it is now 30 of 43 rather than 30 of 39, which makes the
> coverage caveat at the foot of this file stronger, not weaker.

Script: `design_sensitivity.py` · Output: `design_sensitivity.json`
Data acquired for this: 100 hPa eddy heat flux 45–75°N, 1958–2024, from NCEP/NCAR via
NOAA PSL OPeNDAP (`03_data_ingestion/acquire_heatflux.py`, 24,472 daily values,
winter mean +13.8 ± 10.3 K m/s).

## The headline: the two estimators agree

Per event, same frozen catalogue, same window (+15..+29 d), both DOY-adjusted against the
same non-event seasonal cycle (n = **30 of 43** events with sufficient baseline):

| estimator | mean AO response |
|---|---|
| conventional (day-of-year-matched, non-event winters) | **-0.9536** |
| corrected (own winter, >75 d from onset) | **-0.9518** |
| **difference** | **-0.0018** |

**They agree to 0.002 σ.** Once the catalogue is fixed and the seasonal adjustment is
applied consistently to both, the between-winter / within-winter distinction makes
essentially no difference to the estimated response.

## What does move the answer: the catalogue

| contrast | shift in AO estimate |
|---|---|
| estimator (conventional vs corrected, same catalogue) | 0.008 |
| catalogue (union vs primary, same estimator), +15..+29 | 0.143 — **18×** larger |
| catalogue (union vs primary, same estimator), −30..−16 | 0.379 — **47×** larger |

And the catalogue contrast is not merely quantitative: pre-onset significance appears under
the union catalogue and disappears under the preregistered consensus one.

## Consequence for the paper

The working title is *"Study design reshapes estimates of sudden stratospheric warming
impacts."* On this evidence, **the estimator half of that claim does not hold.** What
reshapes the estimate is **which events you call events**, not how you compare them.

That is still a real and useful result — arguably a sharper one, because event catalogues
are chosen casually in this literature (the superseded work here carried four mutually
incompatible ones) while estimator choice is where methodological attention usually goes.
But the paper has to be rebuilt around it, and the current framing must be dropped.

Revised claim, supported by what is now in hand:

> Estimated SSW surface responses are far more sensitive to the event catalogue than to the
> compositing design. A union-across-reanalyses catalogue yields a significant pre-onset
> anomaly that a consensus catalogue does not, while conventional and within-winter
> estimators applied to the same catalogue agree to within 0.01 σ.

## Gate 8 itself: not met

Design sensitivity is not predictable from dynamics — largely because there is almost none
to predict.

| predictor (standardised) | β | robust SE | permutation p |
|---|---|---|---|
| precursor v′T′ (100 hPa, −45..−1 d) | −0.084 | 0.261 | 0.67 |
| vortex u10 (−5..+5 d) | -0.041 | 0.211 | 0.82 |
| polar-cap Z100 (0..+30 d) | -0.265 | 0.201 | 0.18 |
| onset day-of-year | -0.246 | 0.179 | 0.22 |
| **winter AO anomaly** | **+0.261** | 0.199 | 0.13 |

R² = 0.2167, **permutation P = 0.2808** (5,000 shuffles). Nothing survives.

The one suggestive signal is the winter AO anomaly — the direct selection variable —
with pairwise r = +0.353, p = 0.0557. Its sign is what winter selection predicts, and it is
the largest coefficient, but at n = 30 it is not significant and must not be reported as a
mechanism.

## Honest note on scope

Only 30 of 43 events had enough within-winter baseline for the corrected estimator. The
comparison is therefore on a subset, and events in winters densely packed with warmings are
under-represented — exactly the winters where the two designs might diverge most. That is a
real limitation of this test, not a reason to discount it, and it should be probed with a
relaxed baseline gap before the finding is finalised.
