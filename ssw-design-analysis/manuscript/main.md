# Sudden stratospheric warmings shift rather than split the surface response

*Draft for Nature Geoscience (Analysis). Every number is taken from
`ssw-design-analysis/CONSOLIDATED_RESULTS.md` or the result JSONs it indexes; display items are built by
`09_figures/fig*.py` from `results/current/`.*

---

## Abstract

Sudden stratospheric warmings (SSWs) shift the odds of cold, blocked winter
weather, and forecasts initialised at an SSW are more skilful. Research and
forecasting nonetheless sort SSWs into the roughly two thirds that propagate
downward to the surface and the rest. Here we show that this label does not
identify a kind of event and that knowing the event adds no detectable forecast
skill beyond a shifted distribution. In a multi-model experiment imposing one
observed SSW on every member, the contrast between downward and non-downward
members (−1.62σ) is no larger than in paired ensembles with no SSW (−1.59σ;
difference −0.03σ [−0.12, 0.05]), where it is what a threshold on one Gaussian
gives. The SSW triples the chance of a cold fortnight over northern Eurasia, and
the circulation shift predicts it. In 1,517 simulated SSWs, stratospheric
information up to onset adds 0.8% [−0.3, 1.4] to probabilistic skill beyond the
shift, less than on ordinary days (3.7%). Ten operational systems predict
differences between events not significantly better than on ordinary days and
under-forecast the shift and its northern-Eurasian cold. After an SSW, forecasts
should issue the shifted distribution; the label records where a threshold cuts
it.

*(194 words by `wc`; limit 200.)*

---

## Main

Sudden stratospheric warmings, the rapid breakdown of the winter polar vortex,
are often followed by weeks of anomalous surface weather that projects onto a
negative Northern Annular Mode (NAM)^{baldwin21}^. Forecasts initialised at SSW
onset are more skilful than forecasts initialised at other times^{sigmond13}^,
and SSWs are one of the main windows of opportunity for subseasonal prediction.

The prevailing account of this skill is event-specific. SSWs are classified by
whether their signal "propagates downward" to the surface, most often with the
criterion of Karpechko et al.^{karpechko17}^, and about two thirds are said to
have a visible downward impact^{baldwin21}^. Differences between the classes are
attributed to the events' strength, morphology or wave forcing^{rao20,lu26}^,
forecast systems are evaluated on whether they predict the class before
onset^{nebel24}^, and non-downward outcomes are framed as potential forecast
busts^{nebel24}^.

Much of the underlying physics is not in dispute. Nudged-ensemble experiments
show that an imposed SSW drives a negative NAM^{hitchcock14,hong26}^ and shifts
weather regimes towards negative North Atlantic Oscillation states^{leeRW25}^,
and relaxation experiments attribute the contrasting outcomes of the 2018 and
2019 events to tropical rather than stratospheric influences^{knight20}^. The
tropospheric response scales with the strength of the lower-stratospheric
anomaly^{white20}^, and ensembles of distinct events differ in their mean
responses in proportion to their post-onset lower stratosphere^{loeffel26}^.
Composites of downward and non-downward events are known to be partly built by
the classification itself^{white19}^, and large seasonal ensembles find little
in the pre-onset state that distinguishes one event's surface response from
another's^{bett23}^.

What has not been tested directly is whether the label itself carries information: whether
downward and non-downward outcomes are distinct kinds of response, and whether
knowing the event helps a probabilistic forecast beyond the shift that every SSW
imposes. We test both, in the Stratospheric Nudging And Predictable Surface
Impacts (SNAPSI) experiment^{hitchcock22}^, in 1,517 SSWs in 20 CMIP6 members of
10 models, in 43 observed events, and in operational reforecasts, for the
polar-cap circulation and for regional temperature.

### An imposed SSW shifts the surface distribution without splitting it

In SNAPSI's `nudged` ensembles the zonal-mean stratosphere of every member is
relaxed to the observed evolution of an SSW above 50 hPa, with no nudging below
90 hPa; in the `control` ensembles it is relaxed to climatology. We analyse nine
models (4,407 members), the February 2018 and January 2019 events from two
initialisations each, and the September 2019 Southern Hemisphere minor warming.
Averaged over days 8–25 after onset, the causal effect of the SSW on the
polar-cap NAM is −1.11σ of each model's control spread (s.d. 0.27 across eight
models; one model, ECCC, is a Tukey outlier at −3.41σ, and including it gives
−1.36σ; Fig. 1a).

The effect is a translation. Re-centring each ensemble on its own mean and
pooling 36 ensembles, the nudged members (1,798) have the same variance as the
controls (1,805): ratio 0.955 [0.873, 1.051], which excludes an added forced
variance larger than about 0.05σ² (Fig. 1b). A shape test does not reject a pure
translation (Kolmogorov–Smirnov p = 0.37), and neither does the per-event
response in observations (p = 0.37, n = 43) or CMIP6 (p = 0.35, n = 1,517),
although two populations separated by less than about 1σ cannot be excluded at
these sample sizes. Once the variance already present before onset is removed,
the forced between-event standard deviation in CMIP6 is at most 0.27σ (upper
95% bound). This concerns the polar-cap mean; regional forecast spread can still
change after weak-vortex states^{spaeth24}^.

### The downward label is a threshold on one population

We apply the surface conditions of the Karpechko criterion^{karpechko17}^ to
every member over days 8–25 and split each ensemble into downward (DW) and
non-downward (NDW) members. With the same SSW imposed on every member, the DW
minus NDW contrast is −1.62σ; in the same models and initialisations with no SSW
it is −1.59σ (Fig. 2a). On the 25 centre–initialisation pairs where both arms
can form a contrast, the paired difference is −0.03σ [−0.12, +0.05]
(centre-cluster bootstrap). With no SSW the contrast is exactly what a cut at
zero produces on one Gaussian population (within 0.002σ [−0.02, +0.02]). With
the SSW imposed the distribution is shifted, so a single population would give a
slightly larger contrast than a cut through its middle; the observed contrast is instead
smaller than that expectation (by 0.19σ [0.11, 0.30]; in 22 of 25 ensembles),
the opposite of what two distinct populations would give. In the Southern
Hemisphere minor warming the nudged ensembles are wider than the controls, so
the raw contrast is larger (−1.96σ against −1.58σ); against the threshold
expectation for each ensemble's own mean and spread it is again smaller (by
0.20σ, in six of seven ensembles).

What the SSW changes is the rate: 84% of nudged members are classified DW
against 45% of controls (Fig. 2b; 79% against 41% in the Southern Hemisphere).
The rate carries the causal shift; the contrast carries the threshold. The same
holds in observations. Applying the published criterion to ERA5, 97% of the
resulting DW–NDW contrast is reproduced by event-free dates, and the fraction of
events meeting the surface conditions, 69.2%, is reproduced by event-free dates
displaced by the measured shift (75.4% [61.5, 87.2]): "about two thirds" is what
one shifted population produces.

### Knowing the event adds no detectable forecast skill

If downward coupling were a property of events, knowing the event should improve
a forecast of how strongly it couples. We compare two probabilistic forecasts of
each CMIP6 event's surface response (days 8–52): the shifted distribution alone
(the mean and spread of SSW responses in other models, as a function of calendar
date), and an
event-aware forecast that adds the stratospheric state up to onset. Scored with
the continuous ranked probability score on events of held-out models, the
event-aware forecast improves on the shift by 0.8% [−0.3, 1.4] (Fig. 3a). The
same information improves forecasts on size-matched ordinary winter days by 3.7%
[2.3, 5.3], and in none of 200 such samples was it worth as little as after
SSWs. The stratospheric state is informative in general; after an SSW, what it
adds beyond the shift is not detectable. The label itself cannot be a
predictor — it is defined from the outcome — so the test is of the event
information the label is meant to summarise.

Point prediction agrees. Before onset, the out-of-sample R² of the surface
response within a model is 0.054 at SSWs and 0.060 on event-free dates: the
event adds −0.006 [−0.033, +0.020] (95% from the event-free null; [−0.061,
+0.047] when simulation members are also resampled; Fig. 3b). With the post-onset
stratosphere, R² rises to 0.38 at SSWs and to 0.37 on event-free dates: 96% of
that diagnostic skill is present without an SSW, and the event adds +0.015
[−0.028, +0.059] ([−0.061, +0.095]; Fig. 3c).

### Regional cold follows the shift

The same holds for the weather an SSW is feared for. Over northern Eurasia
(50–65° N, 10–130° E), where cold extremes become about twice as frequent after
weak-vortex states^{kretschmer18}^, the imposed SSW lowers the days 8–24
temperature by 0.85σ [0.49, 1.40] and raises the chance of a fortnight colder
than the control's 10th percentile from 0.10 to 0.32 (Fig. 5a,b). The shift of
the polar-cap circulation, acting through the relation between that circulation
and regional temperature that holds without any SSW, predicts 0.31: it accounts
for 88% of the cooling and leaves −0.10σ [−0.36, +0.06] (+0.02σ [−0.02, +0.07]
without ECCC). Within the SSW ensembles the downward label adds nothing to a
member's regional temperature beyond its circulation (−0.04σ [−0.19, +0.09]). In
39 observed SSWs the same relation accounts for two thirds of the 1.0 K
northern-Eurasian cooling, and the remainder (−0.33 K) lies within the range of
event-free dates (p = 0.38; Fig. 5c). The polar-cap index is not the whole of the
regional response: high-latitude Europe cools by a further 0.19–0.27σ, and
mid-latitude North America warms (+0.49σ) where the circulation relation implies
cooling, with the same sign of residual in observations (+0.55 K, p = 0.09) and
in line with the mixed North American signal after weak-vortex
states^{kretschmer18,huang21}^. What an SSW changes regionally is therefore also a
distribution, of the regional variable itself, and the label adds nothing to it.

### Operational forecasts miss the size of the shift, not the class of event

Operational forecasts give the same answer. In seven subseasonal systems, the
skill with which the 500 hPa polar-cap response after 13 SSWs was predicted
before onset lay at the 91st percentile of the skill at random winter dates, not
significantly higher^{nebel24}^. We tested ten systems of the S2S reforecast
archive^{vitart17}^, ECMWF as a discovery set and nine others as a confirmatory
set (Methods). ECMWF forecasts started 2–9 days before 12 SSWs track the
event-to-event differences in the surface response with r = 0.77, against 0.52
on average for event-free dates of the same calendar period and lead (p = 0.11;
p = 0.37 for starts 10–17 days before onset; Fig. 4a). Across 17 SSWs the mean
of the nine further systems correlates with the observed response at r = 0.23,
against 0.52 on event-free dates (p = 0.91), and no single system does
significantly better than on event-free dates (Fig. 4b). Both tests detect
strong event-specific skill in synthetic forecasts (in about 90% of 40 and in
all of 100 trials).

What the forecasts get wrong is the size of the shift. If they missed a second
population of outcomes, the observed values would scatter into both tails of the
ensembles; in ECMWF they do not (4% in the two outer rank bins against 17%
expected). They fall on the downward side, and so they do in the nine further
systems, as registered before their data were retrieved: mean rank 0.40 against
0.54 on event-free dates (p = 0.005; p = 0.012 with each system's own rank level
removed; at most p = 0.022 with any one event left out), lower than the event-free
value in nine of ten systems (Fig. 4c). The
ECMWF forecasts raise the probability of a downward outcome from 0.35 on
ordinary dates to 0.60 after SSWs, whereas the observed rate rises from 0.29 to
0.75: relative to ordinary dates they capture 64% of the observed shift
(−270 Pa against −423 Pa) from starts 2–9 days before onset, and 26% from starts
10–17 days before. Correcting the size of the shift by its leave-one-event-out
mean error did not improve the systems' probabilistic scores significantly
(p = 0.16), a test with little power at 17 events.

Whether models under-represent the surface response to SSWs is disputed:
initialised at onset, one seasonal model overestimated it^{sigmond13}^; ECMWF
over-persists the negative North Atlantic Oscillation after weak-vortex
states^{kolstad20}^; subseasonal systems couple too strongly from the lower
stratosphere to the surface at short lags^{garfinkel25}^; and nudged forecasts of
the 2018 event reproduce the predictable part of its observed
response^{dai25}^. Our test does not separate a missed SSW from weak coupling:
starts 2–9 days before onset include forecasts that do not predict the warming
itself, which would also explain the fall to 26% at longer lead. What it
establishes is the forecast consequence: before onset, operational ensembles place
too little probability on the shifted state, whatever the cause.

The same holds for the regional cold. Registered before the data were analysed,
the primary regional test asks whether observed northern-Eurasian temperature
after SSWs lies on the cold side of the forecast ensembles more than on
event-free dates. It does: mean rank 0.39 against 0.53 across the nine further
systems (p = 0.014; 0.021 with each system's own rank level removed; at most
0.036 with any one event left out), below the event-free value in all ten
systems (Fig. 5d). The ensemble-mean anomaly after these SSWs is −0.3 K where
the observed is −1.3 K. High-latitude Europe points the same way (0.43 against
0.53, p = 0.07); mid-latitude East Asia and North America, where the imposed SSW
cools less or warms, show no such error. On synthetic forecasts this test is
conservative (2% false positives) and weak (it detected a 40% under-prediction of
every anomaly in 23% of trials), so the absence of a signal elsewhere is not
evidence of calibration.

### Event differences follow a general stratosphere–surface relation

Events do differ, and the lower stratosphere after onset tracks how much: across
18 ensembles of distinct events the week-2 100 hPa anomaly correlates with the
later surface response at r = 0.85^{loeffel26}^. Two findings place that relation.
Between members of one SNAPSI ensemble, where the zonal-mean forcing is shared,
the same week-2 relation is weak and is equally present with no SSW (r = 0.11
against 0.09; difference in Fisher z +0.02 [−0.04, +0.09], 36 matched ensembles;
Extended Data Fig. 1). And in CMIP6 the post-onset stratosphere predicts the
surface as well on ordinary winter days as after SSWs (above). The relation is a
general property of stratosphere–troposphere coupling, measurable continuously
after onset, not a feature that sorts SSWs into kinds. Consistently, the field's
archetype pair, February 2018 and January 2019, received comparable forced
shifts and similar forced probabilities of a downward outcome (0.78 and 0.71 in
eight models excluding ECCC, at the primary initialisations), and at those
initialisations the observed outcomes lie within the central 95% of the model
distributions in 17 of 18 cases of nine models and two events (Extended Data Fig. 2).

### Discussion

The downward label measures where a threshold cuts a shifted distribution. Under
an imposed SSW, members fall into both classes in proportions set by the shift,
the contrast between the classes is no larger than a threshold on one population
produces, and knowledge of the event adds nothing detectable to a probabilistic
forecast beyond the shift. None of this denies that SSWs matter or that events
differ: the shift is large and causal, and the post-onset lower stratosphere
carries information, as it does on any winter day.

Two consequences follow. For forecasting, the product after an SSW is the shifted
distribution of the circulation — for the polar-cap NAM, the probability of a
negative state — rather than a categorical prediction of whether this event will
propagate — and, for regional weather, the shifted distribution of the regional
variable itself, which the polar-cap index does not fully carry; the quantity to
verify and correct is the size of that shift,
which operational systems under-predict, not the class of the event, which none
predicts significantly better than it predicts ordinary winter days. For
research, wherever "downward-propagating SSW" is used to stratify impacts, the
rate at which events meet the criterion, set against a matched null, is an
honest quantity; the contrast between the classes is not, because the threshold
manufactures most of it (97% in ERA5; see also ref.^{white19}^).

Five limits bound these conclusions. First, the decisive evidence is from models;
with 42 usable observed events the observational out-of-sample tests are
uninformative, and the operational tests rest on 12 to 17 events. Second, SNAPSI
contains two Northern Hemisphere events, so statements about how events differ
rest on n = 2; the result that the class contrast is no larger than the
threshold's does not depend on n, because it is measured within each ensemble.
Third, SNAPSI nudges only the zonal-mean stratosphere and our CMIP6 fields are
zonal means, so non-zonal vortex geometry, which is associated with stronger
responses^{nebel24}^ and differs between nudged members^{feng25}^, is neither
held fixed nor among our predictors. Fourth, "not detectable" is bounded, not
zero: after SSWs the event-aware forecast's gain is at most 1.4% of the
probabilistic score (upper 95% bound). Fifth, the regional analysis covers mean
temperature in four literature-defined boxes, not precipitation, wind or local
extremes, and observed regional samples of 39 events can exclude only residuals
larger than about a quarter of a standard deviation.

---

## Methods

**Event catalogue.** The 43 observed SSWs over 36 winters 1958–2024 (primary
catalogue, frozen; file checksum in the repository) are the 41 major events of
the NOAA CSL SSW compendium^{butler17}^ detected in at least two thirds of the
reanalyses that cover their date, and the two major warmings of 16 January and
4 March 2024 that post-date it^{leeSH25}^.

**Reanalysis.** ERA5^{hersbach20}^ from the WeatherBench 2 public copy^{rasp24}^
(1.5°, 6-hourly). Polar-cap sea-level pressure is the cos-latitude-weighted mean
over 60–90° N (60–90° S for the SH case), including the 60° row. NAM indices at
1000, 850 and 150 hPa: cos-latitude-weighted geopotential height over the rows
poleward of 65° N (66–90° N on this grid) from the 00 UTC field of each day, as
an anomaly from its day-of-year mean, negated and standardised by day of year;
with true daily means the one label that depends on it (Extended Data Fig. 2c)
is unchanged.

**SNAPSI.** Nudged and control ensembles for CCCma, CNR-ISAC, ECCC, ECMWF, KMA,
Meteo-France, NCAR, SNU and UKMO (38–52 members each) at initialisations
s20180125, s20180208 (February 2018; central date 12 February), s20181213,
s20190108 (January 2019; central date 2 January) and s20190829 (SH, eight
centres without CCCma; central date 18 September 2019). NRL is excluded: at three of its four Northern Hemisphere
initialisations its nudged and control members differ by less than 0.03 Pa. Lead
is measured from 00 UTC on the initialisation date; the time origin of every
ensemble was measured from its files (UKMO and Meteo-France start at 06 UTC).
"Days +8..+25" denotes the 69 six-hourly steps from 00 UTC on day 8 to 00 UTC on
day 25 after the central date, covered in full by every Northern Hemisphere
ensemble and s20190829. s20190108 initialises six days after onset; s20191001
starts after the SH central date and is not used.

**Causal shift and distribution (Fig. 1).** S = mean(nudged) − mean(control) of
the days +8..+25 polar-cap sea-level pressure, divided by each model's control
spread pooled over its initialisations, reported with the NAM sign; standard
errors from the two ensemble variances. For the distribution, members are
standardised by their own ensemble's control mean and s.d. and each ensemble is
re-centred on its own mean before pooling (pooling without re-centring adds the
between-ensemble spread of shifts to the nudged variance only). Variance-ratio
interval from 4,000 bootstrap resamples of members within ensembles.

**Classification and paired test (Fig. 2).** Karpechko et al.^{karpechko17}^
conditions 1–2 over days +8..+25 on the member's NAM proxy −(polar-cap sea-level
pressure − control mean)/control s.d.: window mean negative and more than half of
6-hourly values negative. A contrast needs at least three members in each class,
which removes 11 of 36 nudged ensembles (DW rate 0.99) and no control ensemble;
the DW rate is reported for all ensembles. The paired test uses the 25 centre ×
initialisation pairs in which both arms form a contrast, with 10,000 bootstrap
resamples of the eight centres. The one-population expectation for each ensemble
is the contrast a cut at zero gives on a Gaussian with that ensemble's mean and
s.d., −s φ(a)[1/Φ(a) + 1/(1 − Φ(a))] with a = −m/s (−1.596s at m = 0).

**CMIP6 events and predictors.** CMIP6 zonal-mean fields for 20 members of 10
models (10 members are CanESM5). SSWs follow Charlton and Polvani^{charlton07}^:
the first day from November to March on which the daily zonal-mean wind at 10 hPa,
60° N turns easterly; a further event requires 20 consecutive westerly days in
between; a reversal not followed by 10 consecutive westerly days before 30 April
is a final warming and is excluded. On NCEP–NCAR reanalysis 1958–2023 this
detector reproduces all 39 NCEP–NCAR central dates of the compendium^{butler17}^
to the day, with no additional events. Members with fewer than 15 events are not
used (one member of an eleventh model). Predictors: 10, 50 and 100 hPa zonal wind
at 60° N and over the cap in windows before, at and after onset, plus
seasonality. Response: each member's annular-mode index (leading EOF of
zonal-mean sea-level pressure, 20–90° N) over the November–April days of days
+8..+52 (18% of events, mostly March onsets, have fewer than 44 of the 45 days).
Pseudo-onsets are taken only from days more than 135 days from every real onset,
so no pseudo-event's −60..+75-day neighbourhood overlaps a real event's.

**Forecast value (Fig. 3a).** Two Gaussian forecasts of each event's response,
in the member's own σ units and not demeaned: the shifted distribution (mean from
a calendar regression, day-of-year sine and cosine, fitted on the SSWs of the
other models; spread its residual s.d.) and the event-aware forecast (ridge
regression on the predictors up to onset plus calendar; spread the out-of-fold
residual s.d. within the training models). Leave-one-model-out cross-validation;
continuous ranked probability score in closed form; skill 1 − CRPS(event-aware)/
CRPS(shift) with a 2,000-replicate bootstrap over models. The identical pipeline
is applied to 200 size-matched sets of pseudo-onsets (same members, onset counts
and calendar days).

**Predictability (Fig. 3b).** Ridge regression, five-fold cross-validation grouped
by member, all variables standardised within member, R² on pooled out-of-fold
predictions; null from 1,000 size-matched pseudo-onset sets; a 1,000-replicate
member-cluster bootstrap gives the second interval. In 1,022 of 20,000
member-draws some pseudo-onsets lack enough in-season outcome days, which lowers
the null and so errs against our conclusion. Each draw and replicate has its own
seeded stream.

**Operational reforecasts (Fig. 4).** S2S reforecasts^{vitart17}^ from the ECMWF
Data Store, one model version per system (members per start; hindcast years).
Discovery set: ECMWF (model year 2022; 11; 2002–2021; Monday and Thursday
starts). Confirmatory set: ECCC (2025; 4; 2001–2020), CMA (2022; 4; 2007–2021),
HMCR (2025; 11; 1991–2020), KMA (2026; 7; 1993–2016), CNRM (1 June 2025; 11;
1999–2024), JMA (30 September 2022; 5; 1991–2020), CNR-ISAC (16 October 2023; 8;
2001–2020), NCEP (1 March 2011; 4; 1999–2010) and CPTEC (4 January 2023; 11;
1999–2018). Starts from December to March (KMA January to March), leads to 34
days (ECMWF 42). BoM (whose only version is from 2014), UKMO and IAP-CAS (no
sea-level pressure) were not used. Polar-cap sea-level pressure at 00 UTC
(instantaneous), reduced as for ERA5. Forecast anomalies are departures from the
mean of the other hindcast years at the same start date and lead; observed
anomalies are departures from the mean of the other years at the same calendar
dates (our ERA5 extraction spans November 1998 to April 2022, so it omits the
earliest hindcast years of HMCR, JMA and KMA and CNRM's after 2021). The outcome is the days +8..+25 NAM proxy; an SSW start lies
2–9 days (for ECMWF also 10–17 days) before a catalogued onset, and an event's
value for a system is the mean over its starts.

**Reforecast tests.** Discrimination (D): correlation across events between the
ensemble-mean and observed outcomes. Rank (H1): mean over events of the
probability integral transform of the observed outcome within the ensemble.
Shift correction (H2): the gain in ensemble CRPS from adding to every member the
leave-one-event-out mean of (observed − ensemble mean). Each is compared with
10,000 sets of pseudo-onsets, one per event, within ±21 days of the event's
calendar date in any hindcast year and more than 135 days from every catalogued
onset; for D the sets are restricted to those whose across-event variance of the
observed outcome is within a factor of 1.25 of the events' (s.d. within about
1.12). p values are one-sided: P(null r ≥ observed) for D, P(null rank ≤
observed) for H1 and P(null gain ≥ observed) for H2. H1 and H2 were registered
before the confirmatory data were retrieved; their primary statistic is the mean
over the confirmatory systems covering each event (2–8 systems; 17 events,
1998–2021), and a pseudo-onset enters the null for an event if at least half of
that event's systems have starts before it, those systems being averaged. The
design was validated on synthetic forecasts with every system's real layout and
the real observations (100 datasets per case): rejection rates at p < 0.05 were
0.06 (D), 0.03 (H1) and 0.03 (H2) for forecasts equally skilful at all dates;
0.07 (D) for forecasts without information; 0.66 (H1) and 0.11 (H2) for forecasts
that under-predict every anomaly by 40%; and 1.00 (D) for forecasts skilful only
before SSWs (elsewhere the outcome of another year; a case added after the real
result). For ECMWF alone (40 datasets per case) the corresponding power of D was
0.88–0.90. Two null designs were replaced on
synthetic data before the real forecasts were analysed: moving each event to
another year (too few event-free winters to be calibrated), and requiring all of
an event's systems at each pseudo-onset (no candidates for four events). Two
sensitivities were added after the primary result. Because a null draw averages
only the systems with starts at that pseudo-onset, a system's overall rank level
can enter events and null unequally; removing each system's mean rank over its
pseudo-onsets before averaging gives H1 p = 0.012 (rejection rate 0.02 under
equal skill). CPTEC's forecast climatology rests on one to three other years per
start date; without CPTEC, H1 p = 0.005 and D p = 0.90.

**Regional temperature (Fig. 5).** Regions from the literature: northern
Eurasia 50–65° N, 10–130° E (primary)^{kretschmer18}^, high-latitude Europe
55–70° N, 0–60° E, mid-latitude East Asia 35–55° N, 90–150° E and mid-latitude
North America 35–55° N, 120–60° W^{huang21}^; cos-latitude means over all grid
points. SNAPSI near-surface temperature (6-hourly) of all nine centres, both arms,
four Northern Hemisphere initialisations (3,603 members), interpolated to a
2.5° grid and averaged to daily means from 00 UTC on the initialisation date;
days +8..+24 after the central date, each day complete (four steps). T is
standardised by the control ensemble of the same centre and initialisation. For
each pair, T = a + bN is fitted on control members (N, the member's polar-cap
NAM proxy as above) and applied to nudged members; the residual R is their mean
departure. The cold-fortnight probability predicted from the shift averages, over
nudged members, a Gaussian with mean a + bN and the control residual s.d. The
label test regresses T on N and the downward indicator within each nudged
ensemble (both demeaned). Intervals: 10,000 bootstrap resamples of centres; on
synthetic data the interval covered a true zero in 180 of 200 draws (slightly
liberal) and detected a −0.5σ direct effect in all 200. ERA5 2 m temperature
(WeatherBench 2, 1959 to January 2023; 39 of the 43 events), anomalies from a
31-day-smoothed day-of-year climatology; N is the 1000 hPa NAM over days +8..+25;
a and b are fitted on all November–March dates more than 135 days from every
onset, and R is compared with 10,000 sets of such dates within ±21 calendar days
of each event. The design, regions and falsifier (an interval excluding zero and
|R| > 0.2σ in northern Eurasia) were fixed before any regional temperature was
analysed; the analysis without ECCC was added after.
Operational: daily-mean 2 m temperature ("2t", averaged over each 24 h lead
window; the definition of the daily mean differs between centres) for the same
ten systems, model versions and starts, over the same regions from the S2S
reforecast archive (1.5° grid), and the ERA5 regional series above as the
observation; anomalies, event and pseudo-onset definitions, null and statistics
as for the polar-cap tests, with the outcome the mean over days +8..+24 and the
sign chosen so that a low rank means colder than forecast. The primary test (rank
in northern Eurasia, multi-model mean of the nine confirmatory systems) was
registered before the regional data were analysed; the other regions are
secondary (Holm-adjusted p for northern Eurasia over the four regions, 0.054).
On 100 synthetic datasets per case with the real temperature layouts and
outcome, rejection rates for the rank test were 0.02 for equally skilful
forecasts, 1.00 for forecasts without information and 0.23 for forecasts that
under-predict every anomaly by 40%; the discrimination test rejected 12% of
equally skilful datasets, so it is not interpreted for temperature.

**Precursor coupling (Extended Data Fig. 1).** Week-2 (days 8–14) 100 hPa
polar-cap geopotential height against days 15–25 polar-cap sea-level pressure,
correlated across members within each ensemble; Fisher-z pooling weighted by
n − 3; nudged-minus-control difference on matched ensembles with a normal-theory
interval. A pre-specified gate required the nudged/control spread ratio at
100 hPa to exceed 0.5 (median 0.96).

**Archetypes (Extended Data Fig. 2).** ERA5 polar-cap sea-level pressure at exactly
the members' forecast times, placed in each model's nudged distribution with the
same control base and spread; the pair test compares the observed 2018-minus-2019
difference with all member pairings of the same model.

**Reproducibility.** Every result is produced by one script writing one JSON;
random streams are seeded in the code, producer hashes are recorded and checked,
and the environment is pinned. Defects found during the work are logged
in the repository.

---

## References

[[REFERENCES]]

---

## Figure legends

**Fig. 1 | An imposed SSW shifts the surface distribution without splitting it.**
**a**, Causal polar-cap NAM shift (nudged minus control, days +8 to +25) for nine
SNAPSI models and four Northern Hemisphere initialisations, in units of each
model's control spread; bars, ±1.96 s.e.; black ticks, model means; grey band,
mean ± s.d. over models excluding ECCC. **b**, Members of all 36 ensembles, each
re-centred on its own mean and pooled: control (grey), nudged (vermilion, placed at the
mean shift) and the control translated by the same shift (dashed), in units of
each ensemble's control spread; the shift here is the mean over all 36 ensembles
including ECCC in these per-ensemble units, and so is larger than in **a**.

**Fig. 2 | The downward label is a threshold on one population.** **a**, DW minus
NDW contrast with the SSW imposed (nudged) against the same model and
initialisation with no SSW (control), for the 25 pairs in which both arms form a
contrast; the dashed line is equality and the dotted lines mark the value a cut at
zero gives on a unit Gaussian, −2√(2/π). **b**, Fraction of members classified DW,
all 36 Northern Hemisphere ensembles per arm.

**Fig. 3 | Knowing the event adds no detectable forecast skill after an SSW.** **a**,
Continuous ranked probability skill of the event-aware forecast over the shifted
distribution, on held-out models, after SSWs (point, 95% model-bootstrap
interval) and on 200 size-matched sets of ordinary winter days (grey; mean and
2.5–97.5% range), for
predictors before onset (P1) and up to onset (P2). **b**,**c**, Within-model
out-of-sample R² of the surface response at SSWs (vermilion lines) and on 1,000
event-free sets with the same onset counts (grey histograms), for predictors
before onset (**b**) and including the post-onset stratosphere (**c**).

**Fig. 4 | Operational reforecasts after SSWs.** **a**, ECMWF reforecasts
started 2–9 days before 12 SSWs: ensemble-mean against observed days +8..+25
polar-cap sea-level-pressure response (sign reversed, so negative is the
downward, negative-NAM outcome); dotted, equality. **b**, Correlation across
events between the ensemble-mean and observed responses for each system and for
the mean of the nine confirmatory systems (diamond), against the central 95% of
event-free dates of the same calendar period and lead with matched spread of the
outcome (grey bars); numbers, events per system. **c**, Mean rank (probability
integral transform) of the observed outcome within the ensemble, against the
central 95% of the same event-free dates without the spread restriction; 0.5, a
calibrated forecast.

**Fig. 5 | Regional cold after an imposed SSW is the circulation shift.** **a**,
SNAPSI: effect of the imposed SSW on days 8–24 temperature in four regions
(blue; 95% centre-bootstrap interval), the part implied by each member's
polar-cap circulation through the relation that holds without an SSW (grey), and
the residual (vermilion, with interval), in units of the control spread.
**b**, Probability of a fortnight colder than the control's 10th percentile: no SSW
(0.10 by construction), SSW imposed, and predicted from the circulation shift
alone (hatched). **c**, ERA5, 39 observed SSWs: the same decomposition in kelvin;
whiskers, the central 95% of the residual on event-free dates of the same
calendar period. **d**, Ten S2S reforecast systems, 17 SSWs: mean rank
(probability integral transform) of the observed days 8–24 regional temperature
within the forecast ensembles, for the mean of the nine confirmatory systems
(diamonds) and each system (dots), against the central 95% of event-free dates
of the same calendar period and lead (grey); p, the registered one-sided test
that outcomes are colder than forecast more often than on event-free dates.

**Extended Data Fig. 1 | The lower-stratospheric precursor couples to the surface
as strongly with no SSW.** **a**, For each of 36 SNAPSI ensembles, the correlation
across members between week-2 (days 8–14) 100 hPa polar-cap geopotential height
and the days 15–25 polar-cap surface response, with the SSW imposed (nudged)
against the same model and initialisation with no SSW (control); open squares,
the initialisation that starts after onset. **b**, Fisher-pooled correlations with
95% intervals and the nudged-minus-control difference, for all matched ensembles
and excluding the short-lead initialisation.

**Extended Data Fig. 2 | The field's two archetype events are two draws.** **a**,
Nudged members of nine models for February 2018 (s20180125) and January 2019
(s20181213) as the surface NAM proxy over days +8 to +25 (negative,
downward-propagating); bars, interquartile range; filled diamonds, ERA5; open
diamonds, ERA5 with the model-minus-ERA5 initial offset removed, drawn where
removing it moves the value by more than 0.25σ. **b**, Percentile of the observed
2018-minus-2019 difference among all member pairings of the same model; grey,
central 95%. **c**, Observed NAM means for both events under four variants of the
Karpechko criterion^{karpechko17}^; the only NDW label is January 2019 at
1000 hPa over days 8–52.

---

## Acknowledgements

This work is based on S2S data. S2S is a joint initiative of the World Weather
Research Programme (WWRP) and the World Climate Research Programme (WCRP). The
original S2S database is hosted at ECMWF as an extension of the TIGGE database.

## Author contributions

[To be completed by the authors.]

## Competing interests

[To be completed by the authors.]

## Data availability

- SNAPSI: CEDA archive, `https://dap.ceda.ac.uk/badc/snap/data/post-cmip6/SNAPSI/`
  (free registration for download).
- ERA5: WeatherBench 2 public copy,
  `gs://weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr`
  (anonymous access).
- NCEP–NCAR reanalysis 1: NOAA PSL THREDDS server (`ncep.reanalysis.dailyavgs`).
- CMIP6 zonal means: ESGF; the member list is in the repository
  (`03_data_ingestion/cmip6_ensemble_manifest.csv`).
- SSW compendium: NOAA CSL; the frozen copy's checksum is in the repository.
- S2S reforecasts: ECMWF Data Store, dataset `s2s-reforecasts`
  (`https://ecds.ecmwf.int/`); Copyright © 2026 ECMWF; licence CC BY-NC 4.0
  (`https://creativecommons.org/licenses/by-nc/4.0/legalcode`); ECMWF does not
  accept any liability whatsoever for any error or omission in the data, their
  availability, or for any loss or damage arising from their use; the data were
  reduced here to polar-cap means.

## Code availability

All code, reduced data products, result files and run logs:
`https://github.com/ProgrmerJack/Solar-Magnetic-Analysis` (directory
`ssw-design-analysis/`).
