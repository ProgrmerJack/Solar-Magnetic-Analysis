# Gate 3 — seasonality calibration: PASSED, and the premise was wrong

Script: `calibrate_seasonality.py` · Output: `calibration_results.json`
Data: CPC daily AO 1950–2026, 13,864 winter days, 77 winters, 47 canonical events.

## The premise this gate was built on does not hold

The plan required fixing a pseudo-onset null centred at **−0.258** instead of zero, on the
assumption that three annual harmonics left residual seasonality. That diagnosis was wrong.

The old test drew fake onsets **only from winters that actually contained an SSW**. Those
winters really do carry a depressed post-onset period, so a fake onset placed inside one
lands, on average, partly on the genuine response. The "null" was therefore a diluted copy
of the real signal, not a null.

Reproduced directly, same estimator, same seasonal model, only the winter pool changed:

| fake onsets drawn from | null mean | sd |
|---|---|---|
| all 77 winters | **+0.006** | 0.280 |
| the 38 SSW-hosting winters (what the old test did) | **−0.267** | 0.275 |

That recovers the −0.258 exactly. The estimator was never miscalibrated in the way assumed.

**Consequence for the superseded work:** r113's "observed −0.999 vs null −0.258, P = 0.001"
compared the real effect against a null that already contained part of that effect. The
comparison was therefore **conservative**, and the reported significance if anything
understated.

## Calibration under a correct null

Pseudo-onsets permuting the observed onset day-of-year distribution across the full winter
pool; 400 draws; winter-block bootstrap intervals; tested bin +15..+29 d.

| seasonal specification | null mean | coverage | false-positive rate | verdict |
|---|---|---|---|---|
| 3 annual harmonics | −0.024 | 0.933 | 0.067 | pass (both at the edge) |
| **6 annual harmonics** | **−0.018** | **0.950** | **0.050** | **clean pass** |
| cyclic spline, 8 knots | −0.045 | 0.933 | 0.067 | pass (edge) |
| **cyclic spline, 12 knots** | **−0.012** | **0.950** | **0.050** | **clean pass** |
| leave-one-winter-out climatology | −0.024 | 0.917 | 0.083 | **fails coverage and FPR** |

Criteria: |mean| ≤ 0.05, coverage 0.93–0.97, FPR 0.03–0.07.

## Decisions

- **Primary seasonal specification: 6 annual harmonics.** Clean pass, fewest parameters of
  the two clean passes, and no material difference from the 12-knot spline (−0.018 vs
  −0.012) — satisfying the "no dependence on spline complexity" criterion.
- **Sensitivity: cyclic spline, 12 knots.**
- **Leave-one-winter-out climatology is rejected**, despite being the plan's suggested fix:
  it over-fits, degrading coverage to 0.917 and inflating the false-positive rate to 0.083.
- Any future pseudo-onset null must draw from the **full winter pool**, never only from
  event-hosting winters, and must reproduce the observed calendar distribution.

## Still open

This calibration is for a single station-free index (AO). The plan correctly requires that
**spatial outcomes use region- or station-specific seasonality**; a shared curve is
inadmissible there. That must be calibrated separately before any gridded or
multi-station outcome is analysed, and Gate 3 is only closed for index-type outcomes.
