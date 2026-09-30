# Preregistration — Study design and estimates of sudden stratospheric warming impacts

**Status: DRAFT, not yet frozen.** Freeze and timestamp (OSF) before any published effect
is re-estimated. Nothing in Sections 3–7 may be changed after freezing; changes go in a
dated amendment log with reasons.

**Prepared:** 2026-07-28

---

## 1. Background and motivation

Surface impacts of major sudden stratospheric warmings (SSWs) are usually estimated by
compositing days near an event against a climatology or against days from other winters.
Preliminary work by the author (superseded, retained only as motivation) found that
estimates of the same quantity can change substantially with the comparison design, and
that a specification fitting one lag window against all remaining winter days is unsound
because the reference category contains the post-onset response.

That preliminary work is **not** evidence about the published literature. It examined one
hazard application and five circulation variables, and its own central interpretation
(that pre-onset anomalies reflect winter selection) was contradicted when the design was
corrected: pre-onset anomalies persisted under winter fixed effects, implying genuine
tropospheric precursors. This study is designed to answer the question properly.

## 2. Research question

> How strongly, and in which direction, do event definition, seasonal baseline, precursor
> treatment and comparison design alter estimated surface impacts following major SSWs?

**Explicitly not assumed:** that conventional estimates are too large. The bias may
amplify, attenuate or reverse an association, and the direction is an outcome of this
study, not a premise. Terms "inflate", "artefact" and "upper bound" are barred from the
title, abstract and conclusions unless the sign is established empirically.

## 3. Literature census

### 3.1 Databases and queries
Web of Science Core Collection, Scopus, ADS, and Google Scholar (first 200 hits).
Query (adapted per database syntax):

```
("sudden stratospheric warming" OR "stratospheric sudden warming" OR
 "major warming" OR "polar vortex disruption" OR "vortex weakening")
AND (surface OR tropospheric OR temperature OR precipitation OR "cold air"
     OR wind OR "sea ice" OR snow OR "air quality" OR hazard OR impact)
AND PY = 2000-2026
```

### 3.2 Inclusion criteria
1. Peer-reviewed, English.
2. Reports at least one **observational or reanalysis-based** surface or tropospheric
   quantity estimated relative to SSW events.
3. Uses a defined event set (any catalogue).
4. Reports an effect estimate, composite anomaly, or significance statement.

### 3.3 Exclusion criteria
1. Purely model-experiment studies with no observational estimate (recorded but not coded).
2. Case studies of a single event.
3. Reviews without new estimates.
4. Stratosphere-only outcomes with no surface or tropospheric quantity.

Target corpus **75–120 studies**. If the search returns more, all are screened; if fewer
than 75 pass, the year window is widened to 1990 and this is logged as an amendment.

### 3.4 Coding sheet
For each study: surface outcome; spatial domain; event catalogue and definition; analysis
period; lag window(s); baseline/control construction; climatology calculation; treatment of
pre-onset periods; seasonal adjustment; handling of multiple SSWs per winter; inference and
clustering level; **number of events and number of unique winters**; causal vs associative
language; data and code availability.

### 3.5 Coder reliability
Two independent coders on a random **20% subset**; Cohen's κ reported per field. Fields
with κ < 0.6 are re-defined and the subset re-coded before full coding proceeds.
*The author is a solo independent researcher; a second coder must be recruited before this
section can be executed, and this is a hard prerequisite, not an optional refinement.*

## 4. Replication sample

10–15 studies selected by a **frozen stratified procedure** across domains: AO/NAO and
circulation; regional temperature; precipitation; cold extremes; wind/renewable energy; air
quality; sea ice or snow cover; one hazard application. Within each stratum, studies are
ordered by a seeded pseudo-random permutation of their DOIs and taken in order subject to
data availability.

**Selection is fixed before any re-estimation.** Studies whose published result cannot be
reproduced remain in the sample and are reported as reproduction failures; they are not
replaced.

## 5. Estimators

Primary: mutually exclusive lead–lag bins fitted simultaneously,

    y_it = α_i + s_i(doy_t) + Σ_{k≠k0} β_k · 1[t − T_i ∈ B_k] + ε_it

with bins [−60,−46], [−45,−31], [−30,−16], [−15,−1], [0,14], [15,29], [30,44], [45,60]
and an explicit omitted baseline of days more than 75 d from any onset. Days are assigned
to the nearest central date. Inference by winter-block bootstrap resampling whole winters
(and, for multi-region outcomes, all regions of a winter together).

Comparators, all applied to the same outcome and window: conventional climatological
composite; non-SSW-winter control; single-window winter fixed effects; calendar-matched
design; case-crossover.

Seasonal specification is chosen by the calibration in Section 7 — **not** by fit to the
outcome.

## 6. Event catalogues

- **Primary:** ERA5-detected major mid-winter SSWs from the NOAA CSL compendium
  (Butler et al. 2017), extended by peer-reviewed events through 2024-03-04. Final warmings
  excluded by the compendium's own rule. Analysis period **1958–2024** (the event period;
  outcome records may extend later).
- **Sensitivity:** JRA-55, MERRA-2 and NCEP-NCAR catalogues separately; a consensus
  catalogue requiring detection in ≥4 of 6 reanalyses; and each reanalysis's own central
  date for date uncertainty.

The catalogue is fixed once and never varied by outcome. Both the number of events and the
number of unique winters are reported everywhere.

## 7. Calibration requirement (blocking)

Before any published effect is re-estimated, the estimator must pass, on pseudo-onset dates
that reproduce the **observed calendar distribution** of real onsets:

1. mean estimate within ±0.05 σ of zero;
2. 95% interval coverage 0.93–0.97;
3. false-positive rate 0.03–0.07;
4. no material dependence on spline complexity.

Seasonal specifications compared: 3 and 6 annual harmonics; cyclic cubic splines (two knot
counts); leave-one-winter-out day-of-year climatology. For spatial outcomes a shared
seasonal curve is inadmissible; seasonality must be region- or station-specific.

## 8. Multiplicity

Primary outcome: the +15..+29 d bin of the primary estimator for each replicated study.
All other bins and comparators are secondary. Benjamini–Hochberg FDR at 0.05 across the
family of primary outcomes; the full family size is reported with the results.

## 9. Physical decomposition

Design sensitivity (change between conventional and primary estimates) is regressed on
100-hPa eddy heat flux, polar-cap geopotential height, 10-hPa zonal wind, blocking
frequency, AO/NAO evolution, vortex displacement vs split, downward-propagation strength,
event timing and tropical background state. This is **exploratory** and labelled as such.

## 10. What would falsify the study's thesis

If corrected and conventional estimates agree within their intervals across the replication
sample, the conclusion is that design choice does **not** materially affect SSW impact
estimates. That result will be reported with the same prominence as the alternative.

## 11. Deviations

Any departure from this document is logged in `AMENDMENTS.md` with date, reason, and
whether it was made before or after seeing the affected results.
