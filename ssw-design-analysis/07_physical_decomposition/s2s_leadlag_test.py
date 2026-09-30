#!/usr/bin/env python3
"""
s2s_leadlag_test.py
===================
WHERE IN THE FORECAST CHAIN IS THE SHIFT LOST: THE LOWER STRATOSPHERE, OR BELOW IT?

Plan approved 2026-09-30. This design was committed BEFORE any forecast of
100 hPa geopotential height or of 10 hPa wind beyond lead 15 was retrieved.

QUESTION
  Operational ensembles started 2-9 days before SSWs place the observed days
  +8..+25 polar-cap surface outcome too far on the downward side (P', mean rank
  0.40 against 0.54). Is that because the forecast lower-stratospheric anomaly is
  too weak or too short-lived (a stratospheric persistence deficit, as reported
  for most S2S systems by Garfinkel et al. 2025), or because, given the forecast
  lower stratosphere, the surface responds too little (a coupling deficit)?

DATA
  Forecasts: the ten P' systems, model versions, starts and members; polar-cap
  (60-90N, cos-lat, the 60N row included) geopotential height at 100 hPa, 00 UTC,
  leads 1-34 (acquire_s2s_reforecasts.py --var gh100); zonal-mean u at 10 hPa, 60N,
  leads 16-34 (--var u10_long), joined to the cached leads 1-15. Surface outcome as
  in P'. Observations: ERA5 (WeatherBench 2) 100 hPa geopotential height, 00 UTC,
  the same cap (acquire_era5_psl_cap.py --z100); ERA5 polar-cap msl as in P'.
  Anomalies: leave-one-year-out, as in P'. Z100 anomalies are sign-reversed
  (B = -anomaly), so that NEGATIVE is the weak-vortex direction, like A.

TESTS (starts 2-9 d before onset; the 17 P' events; multi-model mean of the nine
confirmatory systems; the P' calendar-window null with the quorum rule;
one-sided p = P(null <= observed))
  L1 (primary)  rank of the observed days +8..+25 mean B within the ensemble: is
                the lower-stratospheric anomaly after SSWs under-forecast?
  L2 (primary)  conditional surface rank: within each ensemble, regress members' A
                (surface, days +8..+25) on their B (same window); the conditional
                forecast at the OBSERVED B is a + b*B_obs plus the members'
                residuals; rank of the observed A in it. Is the surface
                under-forecast even given the lower stratosphere?
  L3 (descriptive) decomposition of the ensemble-mean surface error per event:
                A_obs - A_ens = b (B_obs - B_ens) + remainder; multi-model means of
                both parts with 10,000-resample event-bootstrap intervals.
  L4 (secondary) the same as L1 for u(10 hPa, 60N) over days +8..+25.
  Holm correction over L1 and L2 (this diagnostic family only; they are not added
  to the paper's primary family, whose register is unchanged).
  Reading, fixed now:
    L1 low, L2 ~ null   -> the shift is lost in the stratosphere (persistence);
    L1 ~ null, L2 low   -> it is lost below 100 hPa (coupling);
    both low            -> both;  neither -> the diagnostic does not locate it.
  ECMWF (discovery system) is reported separately and is not in the primary mean.

Output: results/current/6_predictability/s2s_leadlag_test.json
"""
raise SystemExit("design registered 2026-09-30; implementation follows the data")
