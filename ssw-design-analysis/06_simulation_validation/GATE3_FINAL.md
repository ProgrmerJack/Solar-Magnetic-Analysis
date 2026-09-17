> **⚠ PRE-RE-FREEZE DOCUMENT.** Written before 2026-07-30, when the event
> catalogue was re-frozen from **39 events / 33 winters** to **43 / 36** after the
> consensus rule was found to be era-dependent
> (see `../02_event_catalogues/FREEZE_RECORD.md`). Every count below is the old
> build. The JSONs in `results/` were all regenerated; this document was not.
> Qualitative conclusions are mostly unaffected, but **no number here should be
> quoted without checking it against `results/`.**

# Gate 3 — final status across all three outcome regimes

| regime | required seasonal specification | null mean | coverage | FPR | draws |
|---|---|---|---|---|---|
| **index** (AO, NAO, PNA, AAO) | 3 annual harmonics | **+0.005 sd** | **0.930 (0.895–0.956)** | **0.070** | 300 |
| **spatially aggregated** (state-level) | shared *or* unit-specific — indistinguishable | +0.026 sd | 0.983 (0.911–1.000) | 0.017 | 60 |
| **station-level** | **per-station LOWO day-of-year climatology** | **+0.009 sd** | 0.950 (0.831–0.994) | 0.050 | 40 |

All three regimes meet the criteria. Gate 3 is **closed**. Index-regime figures recomputed 2026-08-04 on the 43-event
catalogue; the spatial and station regimes have **not** been re-run since the
re-freeze and their numbers are still from the 39-event build.

## The rule that emerged, which was not obvious in advance

Seasonal treatment must scale with the heterogeneity of the units, and the requirement is
**not** monotone in flexibility:

* **index outcomes** — 3 harmonics beat 6 harmonics, a 12-knot spline, and a
  leave-one-winter-out climatology. More flexibility made calibration *worse*.
* **aggregated outcomes** — shared and unit-specific are identical; averaging within a
  state has already removed the between-station variation.
* **station-level outcomes** — a shared curve leaves 4× the bias; elevation banding helps
  not at all; only per-station climatology works.

So the same estimator needs the *least* flexible seasonal model for a single index and the
*most* flexible for a station network, and the intermediate case needs neither. A single
blanket rule ("always use splines", "always use climatologies") would have been wrong in
two of the three regimes.

## Standing requirements from this gate

1. **N_BOOT ≥ 1,000** for any reported interval. At 150 replicates the same estimator gives
   coverage 0.920 / FPR 0.080. (The 0.947/0.053 pair once quoted here was the
   39-event build; on the 43-event catalogue 1,000 replicates give 0.930/0.070.)
2. **Seasonal specification by regime**, as tabled above.
3. **Pseudo-onset nulls must remove real-event influence from the DATA first**, not merely
   push fake onsets away from real ones — the latter contaminates the baseline and biased
   the null by +0.234 (43% of baseline days affected).
4. Pseudo-onsets must be drawn from the **full winter pool**, never only from event-hosting
   winters, and must preserve the observed day-of-year distribution.

## Honest limitations carried forward

* Residual positive bias remains in the spatial and station regimes (+0.009 to
  +0.026 sd) but **not in the index regime**, where it is now +0.005 with an
  interval containing zero (recomputed 2026-08-04 on the 43-event catalogue).
  Where it remains it is small, real,
  mechanism unexplained. It is ~1–3% of the response being measured.
* Coverage estimates rest on 40–300 draws; the station and aggregated intervals are wide.
* The seasonal specifications were chosen by exploratory calibration conducted before the
  confirmatory analysis was frozen. Disclosed, not presented as a priori.
* Station-level was calibrated on 103 of 887 stations, stratified by elevation.
