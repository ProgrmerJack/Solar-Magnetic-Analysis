> **STATUS.** The index-regime calibration was **recomputed 2026-08-04 on the
> 43-event catalogue** and the current table is below. Text outside that table
> may still carry pre-re-freeze wording; where the two disagree, the recomputed
> table and `results/current/3_calibration/gate3_clean_null.json` are correct.

# Gate 3 — closure record (index outcomes)

Final calibration on the cleaned series: **8,862** days free of real-event influence,
**73 winters**, frozen primary catalogue (**43 events / 36 winters**), pseudo-onsets
preserving the observed day-of-year distribution, **1,000 winter-block bootstrap
replicates**.

> **⚠ RECOMPUTED 2026-08-04 ON THE 43-EVENT CATALOGUE. The table below is the
> current, reproducible one; the figures previously recorded here were from the
> superseded 39-event build and are struck through underneath.**
>
> The calibration was never re-run after the 2026-07-30 catalogue re-freeze
> (39 → 43 events). More days now fall inside the real-event influence zone
> (5,002 of 13,864 = 36%, against 4,516 = 33%), leaving **8,862 clean days across
> 73 winters** rather than 9,348 across 74. Every number moved.

### Current (43-event catalogue, `gate3_clean_null.py`, N_BOOT = 1000)

| specification | null mean (4,000 draws) | coverage (300 draws, 1,000 reps) | FPR | verdict |
|---|---|---|---|---|
| **3 harmonics** | **+0.0052** [−0.0062, +0.0166] | **0.930** (0.895, 0.956) | **0.070** (0.044, 0.105) | **PASS — primary** |
| 6 harmonics | +0.0227 [0.0113, 0.0342] | 0.930 (0.895, 0.956) | 0.070 (0.044, 0.105) | PASS — sensitivity |
| cyclic spline, 12 knots | +0.0397 [0.0285, 0.0510] | 0.937 (0.903, 0.961) | 0.063 (0.039, 0.097) | **PASS** |

Criteria: |mean| ≤ 0.05; coverage interval contains 0.95; FPR interval contains
0.05. **All three specifications now pass.**

**Three substantive changes from the superseded table.**

1. **The disclosed residual bias is gone.** 3 harmonics gave +0.033 with an
   interval excluding zero, carried here and in `GATE3_FINAL.md` as a real,
   unexplained limitation. On the current catalogue it is **+0.0052, interval
   containing zero**. The bias was a property of the 39-event build, not of the
   estimator.
2. **The 12-knot spline no longer fails.** Its bias falls from +0.063 to +0.0397,
   inside the criterion. It was rejected on that basis; that rejection no longer
   stands on the evidence.
3. **Coverage is 0.930, not 0.947**, and FPR 0.070, not 0.053. Both still pass on
   their intervals, but the 0.947/0.053 figures quoted throughout this project
   are **not reproducible on current data** and must not be requoted.

The ordering that motivated the primary choice **survives**: the least flexible
seasonal model still has the smallest bias (+0.005 < +0.023 < +0.040), so
3 harmonics remains the primary specification.

Coverage counts for 3 and 6 harmonics are again identical (279/300). Both share
the coverage seed and therefore the same pseudo-onset draws; the null means differ
substantially (+0.005 vs +0.023), so the estimates are genuinely different and the
match is shared randomness, not a configuration error. The spline, on the same
draws, gives 281/300 — showing the specifications can and do differ.

### Superseded (39-event build) — retained so the change is visible

| specification | null mean | coverage | FPR | verdict |
|---|---|---|---|---|
| ~~3 harmonics~~ | ~~+0.033 [0.021, 0.044]~~ | ~~0.947 (0.915, 0.969)~~ | ~~0.053~~ | ~~PASS — primary~~ |
| ~~6 harmonics~~ | ~~+0.046 [0.034, 0.057]~~ | ~~0.947 (0.915, 0.969)~~ | ~~0.053~~ | ~~PASS — sensitivity~~ |
| ~~cyclic spline, 12 knots~~ | ~~+0.063 [0.052, 0.075]~~ | ~~—~~ | ~~—~~ | ~~**FAIL** (bias > 0.05)~~ |
| leave-one-winter-out climatology | +0.024 | 0.917 | 0.083 | **FAIL** (coverage, FPR) — not re-run |

## Decisions

- **Primary seasonal specification: 3 annual harmonics.** Smallest residual bias, best
  calibration, fewest parameters.
- **Sensitivity: 6 annual harmonics.**
- **Rejected: leave-one-winter-out climatology** (coverage 0.917, FPR 0.083).
  ~~Rejected: 12-knot cyclic spline (bias +0.063).~~ **The spline rejection is
  WITHDRAWN 2026-08-04** — on the 43-event catalogue its bias is +0.0397, inside
  the |mean| ≤ 0.05 criterion, and it passes coverage and FPR. It is not adopted
  (3 harmonics still has the smallest bias and fewest parameters) but it is no
  longer excluded by the evidence, and any statement that it fails is wrong.
  The *ordering* still holds: added seasonal flexibility monotonically worsens the
  bias (+0.005 → +0.023 → +0.040), which is the opposite of what a
  residual-seasonality explanation predicts.
- **Standing requirement: no reported interval may use fewer than 1,000 winter-block
  bootstrap replicates.** At 150 replicates the same estimator gives coverage 0.920 and
  FPR 0.080. The rule stands; the specific 0.947/0.053 pair it was justified
  with came from the 39-event build and is superseded by 0.930/0.070.

## Disclosed limitations

1. ~~**A residual positive bias of ≈ +0.033 σ remains** and its interval excludes zero.~~
   **WITHDRAWN 2026-08-04.** On the 43-event catalogue the 3-harmonic null mean is
   **+0.0052 [−0.0062, +0.0166]** — an interval containing zero. There is no
   residual bias left to disclose for the primary specification. The old text: that
   is roughly 3% of the response being measured (≈ −0.95 σ at +15..+29 d). It is small but
   real and is reported, not rounded away. Its mechanism is not established.
2. **The seasonal specification was chosen by an exploratory calibration exercise**
   conducted before the confirmatory analysis was frozen. This is disclosed openly rather
   than presented as an a priori choice.
3. **Gate 3 is closed for index-type outcomes only.** Spatial outcomes (gridded fields,
   station networks) require region- or station-specific seasonality; a shared curve is
   inadmissible and has not been calibrated. Gate 3 must be re-run for those before any
   spatial outcome is analysed.

## Recomputation required

`05_corrected_estimators/canonical_event_study.py` was run with 800 replicates and
6 harmonics. Under the rules just fixed it must be re-run with **1,000+ replicates and
3 harmonics** before its numbers are used. Its qualitative conclusions (shape stable across
catalogues; pre-onset significance catalogue-dependent) are unlikely to change, but the
coefficients and p-values are provisional until then.
