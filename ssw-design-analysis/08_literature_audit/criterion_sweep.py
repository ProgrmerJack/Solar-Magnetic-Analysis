#!/usr/bin/env python3
"""
criterion_sweep.py
==================
DOES THE "DOWNWARD" RATE AND CONTRAST TRACK ONE SHIFTED POPULATION FOR EVERY
REASONABLE VERSION OF THE CRITERION, OR ONLY FOR THE ONE WE CHOSE?

Plan approved 2026-09-30. Written before the sweep was run.

QUESTION
  era5_recompute_and_two_thirds.py showed, for the Karpechko et al. (2017)
  criterion as published, that event-free dates displaced by the measured surface
  shift reproduce the observed downward rate, and that event-free dates reproduce
  97% of the class contrast. A referee can answer that the criterion was badly
  tuned. If the label reflected two kinds of event, some versions of the criterion
  would separate SSWs from the shifted population (a plateau: a rate insensitive
  to the threshold, or a contrast the null cannot reach). If it is a threshold on
  one shifted population, every version should track the shifted null smoothly.

GRID (fixed now; 108 specifications in ERA5)
  classification window days 8..W, W in {25, 35, 52}
  level of conditions 1-2: NAM 1000 hPa or 850 hPa (ERA5, era5_nam_daily.parquet)
  condition 1 threshold: window mean below tau, tau in {0, -0.25, -0.5} (sigma)
  condition 2: share of days negative above f, f in {0.5, 0.6, 0.7}
  condition 3 (150 hPa NAM negative on more than 70% of days): off / on
  The published criterion is (52, 1000, 0, 0.5, on).

EVENTS AND NULL (as era5_recompute_and_two_thirds.py)
  The 39 primary-catalogue events inside the ERA5 NAM record; event-free dates
  from gate3_clean_null (more than 135 days from every catalogued onset), drawn
  day-of-year matched, 400 sets. Shift per specification: mean of the window-mean
  NAM at events minus that at event-free dates (each level measured on its own;
  with condition 3 on, the 150 hPa field is shifted by its own measured shift).
  Outcome for the contrast: the per-event 1000 hPa NAM of selection_on_outcome
  (its outcome window), as in the published-criterion audit.

STATISTICS PER SPECIFICATION
  rate at SSWs; rate of shifted event-free sets (mean, 2.5-97.5%);
  z = (rate_SSW - mean) / s.d. over shifted sets;
  contrast DW - NDW at SSWs and the share reproduced by unshifted event-free sets.

TESTS
  S1 (primary, conditions 1-2 only, 54 specifications): max |z| over
     specifications, referred to the same statistic computed for each shifted
     event-free set against the others (the joint null, correlation between
     specifications kept); p = P(null max|z| >= observed).
  S2 (primary): mean z over the 54 (a systematic excess), same reference.
  S3 (secondary): S1 and S2 with condition 3 on (54 specifications); the
     stratospheric shift is not a surface shift, so a failure here speaks to the
     stratospheric condition, not the surface classification.
  Descriptive: the share of the contrast reproduced, per specification.
  CMIP6 (secondary): conditions 1-2 on each member's annular mode, W x tau x f
     (27 specifications), 1,517 events, 20 pseudo sets per member; S1 and S2.
  Reading, fixed now: S1 or S2 at p < 0.05 in ERA5 means some version of the
  criterion separates SSWs from one shifted population, and the paper says so.

Output: results/current/9_literature/criterion_sweep.json
"""
raise SystemExit("design registered 2026-09-30; implementation follows")
