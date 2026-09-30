#!/usr/bin/env python3
"""
obs_model_compatibility.py
==========================
IS THE OBSERVED RECORD A TYPICAL 43-EVENT DRAW FROM THE CMIP6 SSW POPULATION, AND
DOES THE OBSERVATIONAL PIPELINE RECOVER THE TRUTH WHEN A MODEL PLAYS REALITY?

Plan approved 2026-09-30. Written before it was run.

DATA (as is_downward_propagation_a_class.py)
  Observations: CPC AO daily, primary catalogue (43 events), per-event mean over
  days +8..+52, anomalies from the event-free day-of-year climatology; 300
  day-of-year-matched event-free sets (gate3_clean_null).
  CMIP6: each member's annular-mode index, Charlton-Polvani events (1,517), the
  same outcome; per member 20 event-free sets (stratifier_law_cmip6.draw_pseudo).
  The observed index is an EOF of 1000 hPa height and the CMIP6 one of zonal-mean
  sea-level pressure; both are standardised; this is stated as a limit.

STATISTICS (the observational pipeline, computed identically on any event set)
  T1 pure-shift KS statistic (events against the event-free pool plus the shift)
  T2 variance ratio, events / event-free pool
  T3 implied shift, mean(events) - mean(pool)
  T4 downward rate: conditions 1-2 of the criterion on the same series over days
     8-52 (mean negative, more than half of days negative)
  T5 share of the DW-NDW contrast (same outcome) reproduced by event-free sets

A. COMPATIBILITY. 10,000 draws of 43 events without replacement from the pooled
   CMIP6 events (with the event-free sets of the members drawn) and, for each
   model with at least 43 events, 2,000 draws from that model. For each statistic,
   the two-sided percentile of the observed value; combined by the minimum of the
   five two-sided p values, calibrated on the same draws. Reading: p < 0.05 means
   the observed record is not a typical CMIP6 draw on that statistic. Passing
   shows compatibility, which is weaker than exchangeability, and is described so.

B. PERFECT MODEL. Each model in turn is "reality": 1,000 draws of 43 of its events,
   the observational pipeline applied, and the answers compared with that model's
   full-population values:
   PM1 size of the pure-shift KS test at 5% (rejection rate);
   PM2 how often the observed downward rate of the draw falls outside the 95%
       interval of shifted event-free sets (the two-thirds test's error rate);
   PM3 coverage of the model's full-population variance ratio by the 43-event
       bootstrap 95% interval (2,000 resamples).
   Pass, fixed now: PM1 and PM2 rates within [0.01, 0.10]; PM3 coverage >= 0.90.

Output: results/current/6_predictability/obs_model_compatibility.json
"""
raise SystemExit("design registered 2026-09-30; implementation follows")
