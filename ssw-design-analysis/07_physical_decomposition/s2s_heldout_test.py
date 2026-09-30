#!/usr/bin/env python3
"""
s2s_heldout_test.py
===================
HELD-OUT EVALUATION: DO FORECASTS UNDER-PREDICT THE SHIFT AFTER THE 2023-2024 SSWs?

Plan approved 2026-09-30. This design was committed BEFORE any forecast of model
year 2025 (ECMWF, CMA) or any ERA5 value after 2023-01-10 was retrieved, and
before any statistic at these dates was computed from the CNRM files already on
disk. Whatever it shows is reported.

WHAT HAD BEEN SEEN
  - Published descriptions of the events: Lee, Butler & Manney (2025, Weather 80,
    45-53) describe January 2024 as not exerting a canonical surface influence and
    March 2024 as weakly coupled.
  - The CNRM reforecast files (s2s_cnrm_psl_cap.parquet, s2s_cnrm_t2m_regions.parquet;
    hindcast years 1999-2024) were retrieved for the confirmatory test; their
    2023/2024 starts entered only the leave-one-year-out forecast climatology of
    other years. No outcome, rank or ensemble statistic at the dates below was
    computed or looked at.

EVENTS (the frozen primary catalogue, load_catalogue("primary"); no new detection)
  2023-02-16, 2024-01-16, 2024-03-04: the catalogued onsets after the end of the
  ERA5 series used by the operational tests (polar cap to 2022-04-30, regional
  temperature to 2023-01-10). They were never part of any operational test.

FORECASTS (one model version per system, Dec-Mar starts, as in the P' design)
  cnrm       model version 1 June 2025, hindcast years 1999-2024 (on disk)
  ecmwf2025  ECMWF model year 2025, hindcast years 2005-2024, odd-day starts, 11 members
  cma2025    CMA model year 2025, hindcast years 2010-2024, Jan-Mar starts (none in
             December), 4 members
  Polar-cap msl leads 10-34 and daily-mean 2 m temperature windows 10-33, reduced
  exactly as in P' and the regional test (acquire_s2s_reforecasts.py).

OBSERVATIONS
  The existing ERA5 series (WeatherBench 2) wherever they exist; after their end,
  ARCO-ERA5 reduced to the same grid and regions (acquire_era5_arco_extension.py).
  Before splicing, the mean ARCO-minus-WB2 difference over the overlap (00 UTC
  polar cap: Nov-Apr 2020/21 and 2021/22; regional temperature: those winters plus
  2022-11-01..2023-01-10) is removed, one constant per series; the constant, the
  daily RMS difference and the correlation are reported. If the correlation of
  daily values over the overlap is below 0.99, the splice is reported as failed and
  the held-out test is not run.

STATISTICS (identical to P' and the regional test: s2s_forecast_test.start_stats,
leave-one-year-out anomalies, starts 2-9 days before onset, outcome days +8..+25
(polar cap) and +8..+24 (temperature), sign so that a LOW rank means the observed
outcome lies on the downward / cold side of the ensemble)
  HO1 (primary)   mean over the held-out events of the mean rank over the systems
                  that have starts 2-9 d before the event; null: 10,000 draws of
                  pseudo-onsets within +-21 calendar days of each event in any
                  hindcast year of those systems, > 135 d from every catalogued
                  onset, a pseudo-onset qualifying if at least half of the event's
                  systems have starts before it (the P' quorum rule);
                  p = P(null mean rank <= observed), one-sided.
  HO2 (secondary) the same for northern-Eurasian (NEURASIA) temperature.
  HO3 (secondary) pooled: the 17 P' events (nine confirmatory systems) together
                  with the held-out events (their systems), each event against its
                  own systems' null; p as above; for the polar cap and NEURASIA.
  Per-system and per-event ranks are reported, with the number of starts.
  Reading, fixed now: HO1 below the null mean is consistency with the registered
  result; p < 0.05 is an independent replication. With three events and about
  three systems the power is low, and a non-significant result is reported as
  uninformative, not as a failure to replicate, unless the mean rank lies above
  the null mean.

Output: results/current/6_predictability/s2s_heldout_test.json
"""
raise SystemExit("design registered 2026-09-30; implementation follows the data")
