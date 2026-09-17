# Gate 3 for spatial outcomes — phase 1 (null means)

Script: `calibrate_spatial.py` · Output: `spatial_calibration.json`
Test case: SNOTEL daily mean air temperature aggregated to (state, day) —
**81,836 state-days, 10 units, 47 winters, 1980–2026**, 29 events in range.
Real-event influence removed first (32,342 rows, 40%), leaving 49,494 clean rows.
Temperature was chosen deliberately as the hardest case: strongest and most
unit-dependent seasonal cycle in the archive.

## Result: unit-specific seasonality makes no difference here

| seasonal specification | null mean | in sd units | sd |
|---|---|---|---|
| shared 3 harmonics (one curve for all units) | **+0.157 °C** | +0.026 | 0.598 |
| unit-specific 2 harmonics | +0.186 °C | +0.031 | 0.595 |
| unit-specific 3 harmonics | **+0.157 °C** | +0.026 | 0.598 |

Shared and unit-specific 3-harmonic specifications differ by **0.0002 °C** — they are tied,
not ranked. The script's "least-biased" selection picked `unit_harm3` by that margin, which
is meaningless; the honest statement is that the three are indistinguishable.

**This contradicts the working assumption.** The plan stated that for spatial data "a shared
global seasonal curve is inadequate" because seasonality differs by elevation, latitude and
continentality. On this test it is not inadequate: interacting the seasonal curve with unit
changes the null bias by less than 0.03 °C.

The likely reason is that the design already carries **unit × winter fixed effects**, which
absorb each unit's level in each winter. What remains for the seasonal term to explain is
only the *shape* of the within-winter cycle, and at this aggregation the shapes are similar
enough that sharing one curve costs nothing.

## Limitation that materially qualifies this

The units are **US states** — aggregates of many stations. Averaging within a state already
removes most of the elevation and aspect variation that motivated the concern. **A
station-level test, where a 3,000 m and a 1,200 m site sit in the same unit, is the harder
case and has not been run.** This result should therefore be read as: *shared seasonality is
adequate for spatially aggregated outcomes*, and as saying nothing yet about station-level
or fine-gridded outcomes.

## Residual bias

+0.026 sd, close to the +0.033 sd found for index outcomes, so the spatial case is not
worse-behaved than the index case on this measure. As there, the bias is small, real, and
its mechanism is unexplained.

## Status

Phase 2 (coverage and false-positive rate at 1,000 replicates) is running for the tied
best specification. Gate 3 for spatial outcomes stays **open** until it reports.

---

## Phase 2 — coverage and false-positive rate

Unit-specific 3 harmonics, 60 draws × 1,000 winter-block bootstrap replicates:

| metric | value | 95% CI | criterion | verdict |
|---|---|---|---|---|
| null mean | +0.026 sd | [+0.023, +0.029] | \|mean\| ≤ 0.05 | pass |
| coverage | **0.983** | (0.911, 1.000) | CI contains 0.95 | pass, but see below |
| false-positive rate | **0.017** | — | ≈ 0.05 | conservative |

The criteria are met, but the point estimates say something the pass/fail line hides:
**the intervals are too wide, not too narrow.** Coverage 0.983 against a nominal 0.95, and a
false-positive rate of 0.017 against 0.05, mean the estimator is *conservative* on spatial
data — the opposite of the index case, which sits at 0.930 / 0.070 on the 43-event
catalogue (the 0.947 / 0.053 once quoted here was the 39-event build).

Conservative is far safer than anti-conservative: it costs power, not validity. A real
effect is harder to detect, but a spurious one is not easier to claim.

Two reasons for caution about how firm this is:

1. **n = 60 draws.** The Clopper–Pearson interval spans 0.911–1.000. It contains 0.95, but
   it would contain almost anything above 0.91. This is weak evidence for "nominal
   coverage" and is quoted as such, not as a demonstration.
2. **29 events, not 39.** SNOTEL begins in 1980, so ten catalogue events fall outside the
   record. Fewer events with strong cross-unit correlation is exactly the regime where a
   winter-block bootstrap tends to over-cover.

## Gate 3 (spatial): conditionally closed

**Closed for spatially aggregated outcomes**, with the seasonal specification a free choice
(shared and unit-specific are indistinguishable here) and the estimator conservative rather
than anti-conservative.

**Not closed for station-level or fine-gridded outcomes.** The concern that motivated this
gate — seasonality differing with elevation, latitude and aspect — is largely averaged away
by aggregating to states. A station-level calibration is still required before any
station-resolved or gridded outcome is analysed, and nothing in this result licenses that.
