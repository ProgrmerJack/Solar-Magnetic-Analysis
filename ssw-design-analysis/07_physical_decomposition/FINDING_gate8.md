> **⚠ PRE-RE-FREEZE DOCUMENT.** Written before 2026-07-30, when the event
> catalogue was re-frozen from **39 events / 33 winters** to **43 / 36** after the
> consensus rule was found to be era-dependent
> (see `../02_event_catalogues/FREEZE_RECORD.md`). Every count below is the old
> build. The JSONs in `results/` were all regenerated; this document was not.
> Qualitative conclusions are mostly unaffected, but **no number here should be
> quoted without checking it against `results/`.**

# Gate 8 — and a result that overturns the paper's framing

Script: `design_sensitivity.py` · Output: `design_sensitivity.json`
Data acquired for this: 100 hPa eddy heat flux 45–75°N, 1958–2024, from NCEP/NCAR via
NOAA PSL OPeNDAP (`03_data_ingestion/acquire_heatflux.py`, 24,472 daily values,
winter mean +13.8 ± 10.3 K m/s).

## The headline: the two estimators agree

Per event, same frozen catalogue, same window (+15..+29 d), both DOY-adjusted against the
same non-event seasonal cycle (n = 30 of 39 events with sufficient baseline):

| estimator | mean AO response |
|---|---|
| conventional (day-of-year-matched, non-event winters) | **−0.902** |
| corrected (own winter, >75 d from onset) | **−0.910** |
| **difference** | **+0.008** |

**They agree to 0.008 σ.** Once the catalogue is fixed and the seasonal adjustment is
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
| vortex u10 (−5..+5 d) | −0.037 | 0.212 | 0.84 |
| polar-cap Z100 (0..+30 d) | −0.247 | 0.200 | 0.21 |
| onset day-of-year | −0.246 | 0.181 | 0.22 |
| **winter AO anomaly** | **+0.265** | 0.201 | 0.12 |

R² = 0.210, **permutation P = 0.30** (5,000 shuffles). Nothing survives.

The one suggestive signal is the winter AO anomaly — the direct selection variable —
with pairwise r = +0.352, p = 0.057. Its sign is what winter selection predicts, and it is
the largest coefficient, but at n = 30 it is not significant and must not be reported as a
mechanism.

## Honest note on scope

Only 30 of 39 events had enough within-winter baseline for the corrected estimator. The
comparison is therefore on a subset, and events in winters densely packed with warmings are
under-represented — exactly the winters where the two designs might diverge most. That is a
real limitation of this test, not a reason to discount it, and it should be probed with a
relaxed baseline gap before the finding is finalised.
