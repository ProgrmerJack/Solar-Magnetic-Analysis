# Finding: "downward propagation" is a threshold on a continuum, not a class

**Status: established in a 1888-event CMIP6 ensemble with an explicit power bound.
Consistent with, but not established by, 43 observed events.**

Script: `07_physical_decomposition/is_downward_propagation_a_class.py`
Data: `is_downward_propagation_a_class.json`

## The implicit physical model, never tested

The literature treats "downward-propagating" (DW) SSWs as a CLASS. Karpechko et
al. (2017) split events in two; roughly two thirds are said to propagate; papers
compare DW against NDW composites and report the contrast as an effect size. That
framing carries a physical claim: **there are two populations of SSWs, one that
couples to the surface and one that does not.**

Two populations predict a mixture-shaped distribution of per-event surface
response. One population predicts a unimodal one, in which case DW/NDW is a cut
through a continuum, "two thirds propagate" is the fraction below a threshold
rather than a property of nature, and comparing the groups compares one
distribution's tails with itself.

Nobody appears to have plotted that distribution or tested it.

## Three tests, two systems, one answer

Per-event day-of-year-adjusted surface anomaly over days +8..+52.

| test | observations (n=43) | CMIP6 (n=1888, 10 models) |
|---|---|---|
| dip statistic, pseudo-event calibrated | p = 0.964 | p = 0.846 |
| Gaussian mixture, BIC k=1 vs k=2 | p = 0.140 | p = 0.865 |
| pure-shift sufficiency (KS) | p = 0.326 | p = 0.191 |
| implied shift | **-0.892 sigma** | **-0.511 sigma** |
| skew / excess kurtosis | +0.274 / -0.691 | -0.222 / -0.039 |

The post-SSW surface response distribution is a **pure translation of the ordinary
winter distribution**. One component suffices; no second mode; a shift alone
reproduces the whole shape.

## The power bound, which is what makes this a claim rather than a null

Failing to reject a mixture is not evidence against one. This project has already
mistaken one absence of evidence for a finding (the AAO negative control) and
nearly a second (the phantom bootstrap under-coverage). So the exclusion is
quantified: two-component mixtures were simulated at n=1888 with mixing fraction
2/3 (the literature's own figure) and total variance held at the observed value,
so the mixture cannot be detected by spread alone.

| component separation | power to detect |
|---|---|
| 0.25 sigma | 0.00 |
| 0.50 sigma | 0.00 |
| 0.75 sigma | 0.01 |
| **1.00 sigma** | **0.86** |
| 1.50 sigma | 1.00 |

**The argument.** The DW-minus-NDW contrast actually reported is +1.143 sigma in
this ensemble and -1.782 sigma (AO) in observations. A genuine two-population
structure producing a contrast of that size implies a component separation of at
least ~1 sigma. That is exactly where the test has 86-100% power. It is not there.

**What is NOT excluded:** any mixture with separation below ~0.75 sigma, where
power is essentially zero. The claim is bounded accordingly and must be stated
that way — not as "there are no classes of SSW".

## The shift predicts "about two thirds", the field's headline number

Baldwin et al. (2021), *Reviews of Geophysics* 59, e2020RG000708, section 7.2 --
the discipline's consensus statement -- says verbatim:

> "Most studies agree that about two thirds (Charlton-Perez et al., 2018;
> Domeisen, 2019; White et al., 2019) of SSW events are characterized as having a
> visible downward impact (e.g., persistent negative phase of the NAM or NAO in
> the lower troposphere and/or the lower stratosphere, (e.g., Domeisen, 2019;
> Karpechko et al., 2017)."

quoted with no uncertainty range, no null-model comparison, and no statement of
what fraction of random winter dates would pass the same test. The review adds:

> "it is still impossible to predict which individual SSW will have a visible
> downward impact—meaning that the tropospheric anomalies (e.g., NAM index or
> pressure) are of the same sign as as those in the stratosphere."

(the "as as" is in the published text). The label is defined by the SIGN of the
post-onset tropospheric anomaly, and the review treats it as a detection problem,
never as a selection problem.

The one-population model makes a hard prediction here: displacing event-free dates
by the measured shift must reproduce the observed classification rate. Tested on
ERA5 (`08_literature_audit/era5_recompute_and_two_thirds.py`, 39 events):

| | surface conditions (1 AND 2) pass rate |
|---|---|
| observed at real SSWs | **69.2%** |
| pseudo-events, no shift | 30.9% |
| pseudo-events shifted by the measured -0.639 sigma | **74.4%  [61.5, 87.2]** |

The interval contains the observed rate. **"About two thirds of SSWs propagate
downward" is what one shifted population produces.** No classes are required to
explain the field's headline statistic.

## Body count in ERA5, in the literature's own fields

`era5_nam_daily.parquet`: polar-cap NAM at 1000/850/150 hPa from WeatherBench2's
public ERA5, 1959-2023, validated at corr = +0.839 against the CPC AO on 11,609
shared days.

| variant | DW rate (null) | DW | NDW | contrast | from selection | reproduced under null |
|---|---|---|---|---|---|---|
| Karpechko, 1000 hPa | 54% (11%) | -0.932 | -0.248 | **-0.684** | **96%** | **97%** |
| ACP 2026, 850 hPa | 59% (13%) | -0.878 | -0.239 | **-0.639** | **102%** | **93%** |

Residuals -0.026 (p=0.589) and +0.013 (p=0.906). ACP 26, 3723 (2026) publishes
-0.762 / +0.088, contrast -0.850, from the same criterion; the values here are the
same sign and comparable magnitude on a different catalogue and NAO construction,
so they corroborate scale without claiming to reproduce their digits.

## What is NOT being claimed

**SSWs do have a surface effect.** The shift is real and large: -0.639 sigma in
ERA5, and the full criterion separates real events from random dates decisively
(54-59% pass versus 11-13% under the null). The all-event composite stands.

What carries no information is the SUB-CLASSIFICATION. Splitting SSWs into
propagating and non-propagating tells you where an event sits in one shifted
distribution, not which of two mechanisms produced it, and the between-group
contrast that split generates is 93-97% reproducible with no event present.

## Why this and the bias result are the same result

`FINDING_stratifier_law.md` showed the DW/NDW surface contrast is 85-99%
reproducible by applying the identical criterion to random winter dates. This
explains why: if there is only one population, then splitting it on its own tail
must manufacture a contrast, and the size of that contrast is fixed by the
regression slope beta and the separation the cut creates. The methods result and
the physical result are one result seen from two sides.

It also disposes of the natural objection to the methods result — "but the groups
really are different" — on its own terms. In a sample thirty times larger than the
observational record, they are not different in the way "class" requires.

## Two errors caught in producing this

1. **A broken dip test nearly produced the opposite conclusion.** The first run
   reported `dip = 0.2468, p = 0.0000, BIMODAL` at n=1888. The implemented
   statistic measures the distance from the ECDF to a globally convex or concave
   CDF — i.e. to a MONOTONE density — not to the nearest unimodal one. Any
   Gaussian sample scores ~0.25 on it; a uniform scores ~0. Calibrating against
   the uniform (correct for Hartigan's exact dip, wrong for this statistic)
   therefore guarantees p=0 for perfectly unimodal data. Recalibrating against
   resampled pseudo-events — a single population by construction, matched n and
   shape — makes the test valid whatever the statistic measures, and flips the
   verdict to p=0.846. **The contradiction between the three tests is what exposed
   it; a single-test design would have published the artefact.**
2. Observations cannot settle this. n=43 gives essentially no power against any
   mixture, and the observational row above is reported as consistency, not
   evidence.

## Limits

- Established in MODELS. Ten CMIP6 models reproduce the SSW surface composite with
  varying fidelity (three barely at all — see `REVIEW_BRIEF.md` result E), so a
  model-derived distributional claim needs the observational record to be
  extended before it can be asserted for the real atmosphere.
- Mixing fraction fixed at 2/3 for the power analysis; other fractions shift the
  detectable separation and were not swept.
- The annular mode in models is EOF-based from sea level pressure, not the CPC
  index; the two constructions are not identical.
