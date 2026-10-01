#!/usr/bin/env python3
"""
s2s_heldout2026_test.py
=======================
A FOURTH OUT-OF-SAMPLE SSW: 4 MARCH 2026, SCORED WITH REAL-TIME FORECASTS

Plan approved 2026-09-30. This design was committed BEFORE any forecast for 2026,
any 2026 reforecast of ECCC/HMCR or model year 2026 of ECMWF, and any surface
observation (sea-level pressure, 2 m temperature) after 30 April 2024 was
retrieved or examined.

WHAT HAD BEEN SEEN
  Only the stratospheric wind: the project's Charlton-Polvani detector
  (ensemble_precursor.detect_ssw, validated on the NCEP compendium) run on NCEP
  u(10 hPa, 60N) daily means to 17 March 2026 (where the NCEP series ends),
  continued with ERA5 (ARCO) 00 UTC values (correlation 0.976 with NCEP over the
  overlap), finds one onset in winter 2025-26: 4 March 2026 (easterly 4-10 March,
  westerly 11-23 March, final warming about 10 April). None in winter 2024-25. The
  onset is marginal (a late, short reversal), which is stated. The earlier held-out
  test (s2s_heldout_test.py) had failed to replicate; this design was written after
  that result and its exploratory diagnosis were known.

DISCLOSURE ADDED 2026-10-01 (after registration, before any test was run): while
  checking the date range of data/processed/atmospheric/ao_daily_cpc.txt, its last
  two lines were displayed: the CPC daily AO for 30 and 31 March 2026 (+2.53 and
  +1.96), i.e. days 26-27 after the 4 March onset, inside the outcome window. No
  other surface observation for 2026 has been examined. The design is unchanged.

FORECASTS (ECDS dataset s2s-forecasts, real time; starts 2-9 days before onset,
i.e. 23 February - 2 March 2026, each start used only if a reforecast of the same
system and model version exists for the same start month-day)
  PRIMARY systems, whose 2026 real-time model is the one of available reforecasts:
    ECMWF (reforecasts model year 2026, on the fly, 2006-2025), ECCC (model year
    2026), HMCR (model year 2026), KMA (model year 2026), NCEP (CFSv2, the only
    version, fixed hindcasts 1999-2010).
  SECONDARY (version match not verifiable): CMA (latest reforecasts model year
    2025), CNRM (2025), JMA (2022), CNR-ISAC (2023).
  Real-time ensembles are larger than reforecast ensembles; the rank of the observed
  outcome is (below + 0.5 equal + 0.5)/(n + 1) with the real-time members, whose
  expectation under calibration is 0.5 whatever n; this is stated as a difference
  from the null draws, which use reforecast members.

ANOMALIES AND OUTCOMES (as P' and the regional test)
  Forecast anomaly: real-time value minus the mean of the reforecasts of the same
  system, model version, start month-day and lead over all their hindcast years.
  Observed anomaly: ERA5 minus its mean over those hindcast years at the same
  calendar dates (WeatherBench 2 series spliced with ARCO-ERA5 exactly as in
  s2s_heldout_test.py, extended to winters 2024-25 and 2025-26 with the same
  reduction; overlap rule unchanged). Outcomes: polar-cap NAM proxy days +8..+25
  (primary) and northern-Eurasian 2 m temperature days +8..+24 (secondary); sign
  so that a low rank means downward / colder than forecast.

TESTS
  HO2026 (primary): multi-system mean rank of the primary systems at 4 March 2026
    against the calendar-window null of their reforecasts (pseudo-onsets within
    +-21 days in any hindcast year, > 135 days from every catalogued onset, quorum
    of half), 10,000 draws; one-sided p = P(null <= observed).
  Secondary: the same for temperature; with the secondary systems added; and the
    four held-out events (2023, Jan 2024, Mar 2024, 2026) together, each against its
    own systems' null.
  Reading, fixed now: a rank below the null mean is consistent with the 1998-2021
  under-forecast, above it inconsistent; with one event neither is decisive, and
  the result is reported as one more event, not as a replication or its failure.

Output: results/current/6_predictability/s2s_heldout2026_test.json
"""
raise SystemExit("design registered 2026-09-30; implementation follows the data")
