# Stratospheric warmings shift winter cold risk rather than creating two kinds of event

*Draft for Nature Geoscience (Article). Every number is taken from
`ssw-design-analysis/CONSOLIDATED_RESULTS.md` or the result JSONs it indexes;
display items are built by `09_figures/fig*.py` and `10_tables/table*.py` from
`results/current/`.*

---

## Abstract

Sudden stratospheric warmings (SSWs) raise the risk of cold, blocked winter
weather, and research and forecasting sort them into the roughly two thirds
whose influence propagates downward to the surface and the rest. Here we show
that this label is a threshold on one shifted distribution, and that what
forecasts miss after SSWs is the size of the shift. Imposing one observed SSW on
every member of multi-model ensembles raises the share of downward outcomes from
45% to 84% without changing the spread of Northern Hemisphere outcomes; forced
differences between events account for a few per cent of any event's label. Over
northern Eurasia the chance of a cold fortnight rises from 10% to about 30% in
the experiment, as the circulation shift alone predicts, and to 26% after 39
observed SSWs. In 1,517 simulated SSWs, knowing the event adds no detectable
probabilistic skill beyond the shift. In ten operational forecast systems,
outcomes after SSWs fall on the negative-annular-mode side of the ensembles more
often than on ordinary winter dates, as much when forecasts predicted the
warming as when they missed it. After an SSW, forecasts should issue, verify and
correct the shifted distribution, not a class.

---

## Main

Cold-air outbreaks over Europe and Asia damage health, energy systems and
transport, and the stratosphere is one of their few sources of predictability
weeks ahead^{domeisenbutler20}^. The largest stratospheric disturbances, sudden
stratospheric warmings (SSWs), are followed by weeks of anomalous surface weather
that projects onto a negative Northern Annular Mode (NAM)^{baldwin21}^; after
weak-vortex states the coldest days over northern Eurasia become about twice as
frequent^{kretschmer18}^, and forecasts initialised at an SSW are more skilful
than forecasts initialised at other times^{sigmond13}^.

Not every SSW is followed by the canonical response, and the field has organised
this variability into two classes. Events are labelled by whether their signal
"propagates downward", most often with the criterion of Karpechko et
al.^{karpechko17}^; about two thirds are said to have a visible downward
impact^{baldwin21}^. Differences between the classes are attributed to the
events' strength, morphology or wave forcing and are used as effect sizes and
stratifiers^{rao20,lu26}^, forecast systems are assessed on whether they
anticipate which events will propagate downward^{nebel24}^, and non-downward
outcomes are framed as potential forecast busts^{nebel24}^.

Several pieces of this picture are already known to be fragile. The warmings
themselves form a continuum, with no clear threshold between major and minor or
between split and displaced events^{coughlin09,maury16}^; composite differences
between downward and non-downward events are partly built by the
classification^{white19}^; large ensembles find little in the pre-onset state
that separates one event's surface response from another's^{bett23}^; and
operational systems predict before onset which events will propagate to
100 hPa but not which will reach the troposphere^{nebel24}^. Nudged experiments
show that an imposed SSW drives a negative NAM^{hitchcock14,hong26}^ and that the
contrasting 2018 and 2019 outcomes owed much to the tropics^{knight20}^, while
ensembles of distinct events differ in their mean response in proportion to
their post-onset lower stratosphere^{loeffel26}^.

What has not been tested directly is whether the label identifies anything in
the Earth system: whether downward and non-downward outcomes are distinct kinds
of response, how much of an individual event's label reflects the event, what
the label means for regional cold, and what, if anything, forecasts get wrong
after SSWs. We answer these with an experiment that imposes the same observed
SSW on every ensemble member in nine models (SNAPSI^{hitchcock22}^), 1,517 SSWs
in 20 CMIP6 simulations of 10 models, 43 observed events, and reforecasts of ten
operational systems, for the polar-cap circulation and for regional temperature.
The reforecast tests and the regional mediation test were checked on synthetic
data; tests fixed before their data were
analysed, and those added afterwards, are listed in Extended Data Table 1.

### One shifted population, cut by a threshold

In SNAPSI's `nudged` ensembles the zonal-mean stratosphere of every member is
relaxed to the observed evolution of an SSW above 50 hPa, with no nudging below
90 hPa; in the `control` ensembles it is relaxed to climatology. We analyse nine
models (4,407 members), the February 2018 and January 2019 events from two
initialisations each, and the September 2019 Southern Hemisphere minor warming.
Averaged over days 8–25 after onset, the imposed SSW moves the polar-cap NAM by
−1.36σ of each model's control spread across the nine models (−1.11σ, s.d. 0.27,
without ECCC, whose −3.41σ is a statistical outlier; Extended Data Fig. 1).

The effect is a translation (Fig. 1a). Re-centring each ensemble on its own mean
and pooling 36 ensembles, the members with the SSW imposed (1,798) have the same
variance as those without it (1,805): ratio 0.955 [0.873, 1.051], which excludes
an added forced variance larger than about 0.05σ². A shape test does not reject
a pure translation (Kolmogorov–Smirnov p = 0.37), and neither does the per-event
response in observations (p = 0.37, n = 43) or CMIP6 (p = 0.35, n = 1,517),
although two populations separated by less than about 1σ cannot be excluded at
these sample sizes. The Southern Hemisphere minor warming is the exception: there
the imposed warming also widens the distribution (variance ratio 1.74 [1.41,
2.15]).

Applying the surface conditions of the Karpechko criterion to every member
splits each ensemble into downward (DW) and non-downward (NDW) members. What the SSW changes is the rate: 84% of members are DW with the SSW imposed
against 45% without it (Fig. 1a; 82% against 45% without ECCC, and 79% against
41% in the Southern Hemisphere minor warming). The same holds in observations:
the share of the 39 observed events with ERA5 outcomes that meet the surface
conditions, 69%, is reproduced by event-free dates displaced by the measured
shift (75% [61, 87]; Fig. 1b), so "about two thirds" is what one shifted
population produces. In the Northern Hemisphere the contrast between the classes
does not change. With the SSW imposed, the DW minus
NDW contrast is −1.62σ; in the same models and initialisations without it,
−1.59σ, and on the 25 centre–initialisation pairs where both arms can form a
contrast the paired difference is −0.03σ [−0.12, +0.05] (centre-cluster
bootstrap; Fig. 2a), while the rate separates the arms in every model
(Fig. 2b). Without an SSW the contrast is what a cut at zero gives on
one Gaussian population (within 0.002σ [−0.02, +0.02]). With it, the contrast is
0.19σ [0.11, 0.30] smaller than a Gaussian cut predicts, but that shortfall
belongs to the Gaussian shortcut, not to the ensembles: half of it comes from the
criterion's second condition (a fraction of days negative), and when the
members without an SSW are shifted by the imposed effect and classified with the
full criterion, they reproduce the contrast of the members with it (difference
+0.04σ [−0.03, +0.10]). In ERA5, likewise, 97% of the DW–NDW contrast that the
published criterion produces is reproduced by event-free dates. In the Southern
Hemisphere minor warming, where the imposed warming also widens the
distribution, the raw contrast is larger (−1.96σ against −1.58σ without it), and
against the threshold expectation for each ensemble's own mean and spread it is
again smaller (by 0.20σ, in six of seven ensembles).

Events do differ, but by little, and not in kind. In CMIP6 the forced variance
between events, after removing the variance already present before onset, is
0.019σ² [−0.033, 0.072]: across events the forced probability of a negative
outcome ranges from about 0.62 to 0.87 (0.47 to 0.93 at the upper bound), and
those differences account for 2% of the variance of a single event's label (at
most 8%). Any one label is therefore almost entirely the shift plus chance. The
field's archetype pair shows it (Fig. 1c): February 2018 and January 2019
received comparable forced shifts and similar forced probabilities of a downward
outcome (0.80 and 0.75 across nine models; 0.78 and 0.71 without ECCC), and at
the primary initialisations the observed outcomes lie within the central 95% of
the model distributions in 17 of 18 model–event cases (Extended Data Fig. 2), although they
are conventionally taken as opposite kinds. The lower stratosphere after onset does
track how much events differ: across 18 ensembles of distinct events the week-2
100 hPa anomaly correlates with the later surface response at
r = 0.85^{loeffel26}^, a relation that survives a control for window overlap
(p = 0.0005), although overlap alone yields r = 0.49 [0.11, 0.77]. Between
members of one SNAPSI ensemble, where the forcing is shared, the same relation
is weak and is equally present with no SSW (r = 0.11 against 0.09; difference in
Fisher z +0.02 [−0.04, +0.09]; Extended Data Fig. 3), and in CMIP6 the
post-onset stratosphere predicts the surface as well on ordinary winter days as
after SSWs (next section): the relation is a general property of
stratosphere–troposphere coupling, measurable after onset, not a feature that
sorts SSWs into kinds. The difference most often proposed, vortex geometry,
cannot be tested in SNAPSI or our zonal-mean CMIP6 fields, so we test it in
reanalysis^{seviour13}^. Events classified as vortex splits tend towards a
stronger early surface response than displacements, as earlier studies
report^{nebel24}^, but not significantly: −0.25σ [−0.61, +0.13] over days 8–25
(p = 0.28; −0.40σ, p = 0.12, with a stricter, persistent classification), with
no difference over days 8–52 (−0.02σ, p = 0.92) and none in the downward rate
(0.78 against 0.60, p = 0.41; 37 events).

### Regional cold follows the shift

The same holds for the weather an SSW is feared for. Over northern Eurasia
(50–65° N, 10–130° E), the imposed SSW lowers the days 8–24 temperature by 0.85σ
[0.49, 1.40] and raises the chance of a fortnight colder than the control's 10th
percentile from 0.10 to 0.32 (Fig. 1d, Fig. 3a,b). The shift of the polar-cap
circulation, acting through the relation between that circulation and regional
temperature that holds without any SSW, predicts 0.31: it accounts for 88% of the
cooling and leaves −0.10σ [−0.36, +0.06]. Without ECCC, the outlier of the
circulation shift, the cooling is 0.60σ [0.44, 0.71] and the cold-fortnight
probability 0.26 against 0.27 predicted, with a residual of +0.02σ [−0.02,
+0.07]. Within the SSW ensembles the downward label adds nothing to a member's
regional temperature beyond its circulation (−0.04σ [−0.19, +0.09]). The
observations agree. After 39 observed SSWs a cold northern-Eurasian fortnight
occurred in 26% [13, 41] of cases, against 10% on event-free dates; the same
circulation–temperature relation accounts for two thirds of the 1.0 K
northern-Eurasian cooling, and the remainder (−0.33 K) lies within the range of
event-free dates (p = 0.38; Fig. 3c).

The polar-cap index is not the whole of the regional response. High-latitude
Europe cools by a further 0.19–0.27σ (its cold-fortnight probability rises to 0.32 against 0.26 predicted from the shift;
0.29 against 0.23 without ECCC), mid-latitude East Asia cools by 0.51σ, of which the circulation accounts for 60%
(0.41σ and 79% without ECCC) (0.27 against 0.20 predicted;
residual −0.21σ [−0.50, +0.03]), and mid-latitude North America warms (+0.49σ)
where the circulation relation implies cooling, with the same sign of residual in
observations (+0.55 K, p = 0.09), in line with the mixed North American signal
during weak-vortex states^{kretschmer18,huang21}^. Over northern Eurasia the
label adds nothing; in mid-latitude North America, where the polar-cap index does
not describe the response, downward members are warmer than their circulation
implies (+0.29σ [+0.06, +0.53]). What an SSW changes regionally is therefore also
a distribution, of the regional variable itself.

### Knowing the event adds no detectable skill; forecasts miss the size of the shift

If downward coupling were a property of events, knowing the event should improve
a forecast of how strongly it couples. We compare two probabilistic forecasts of
each CMIP6 event's surface response (days 8–52): the shifted distribution alone
(the mean and spread of SSW responses in other models, as a function of calendar
date), and an event-aware forecast that adds the stratospheric state up to
onset. Scored with the continuous ranked probability score on events of held-out
models, the event-aware forecast improves on the shift by 0.8% [−0.3, 1.4]
(Fig. 4a). The same information improves forecasts on size-matched ordinary
winter days by 3.7% [2.3, 5.3], and in none of 200 such samples was it worth as
little as after SSWs. The result holds with the post-onset stratosphere among the predictors (17.2% after SSWs against 20.7% on ordinary days), against weak-vortex days that are not SSWs (0.8% against 1.6%, p = 0.94), with models weighted equally (0.4% after SSWs) or CanESM5 excluded (0.6%; the
null was not recomputed without it), and when wave-driving proxies and the pre-onset surface annular mode are added: the richer predictors raise the skill to 4.8% on ordinary days and 2.9% on weak-vortex days but leave it at 0.8% after SSWs. The predictors span the same range after SSWs as on the comparison days (s.d. ratios 1.02 and 1.03), so the gap is not an artefact of a restricted range. The label itself cannot be a
predictor, since it is defined from the outcome, so the test is of the event
information the label is meant to summarise. Point prediction agrees. Before
onset, the out-of-sample R² of the surface response within a model is 0.054 at
SSWs and 0.060 on event-free dates: the event adds −0.006 [−0.033, +0.020]
(Fig. 4b). With the post-onset stratosphere, R² rises to 0.38 at SSWs and to
0.37 on event-free dates, and the event adds +0.015 [−0.028, +0.059] (Fig. 4c):
96% of that diagnostic skill is present without an SSW.

Operational forecasts give the same answer on event differences. In seven
subseasonal systems the skill with which the 500 hPa polar-cap response after 13
SSWs was predicted before onset lay at the 91st percentile of the skill at random
winter dates, not significantly higher^{nebel24}^. We tested ten systems of the
S2S reforecast archive^{vitart17}^, ECMWF as a discovery set and nine others as a
confirmatory set (Methods). ECMWF forecasts started 2–9 days before 12 SSWs track
the event-to-event differences in the surface response with r = 0.77, against
0.52 on average for event-free dates of the same calendar period and lead
(p = 0.11; p = 0.37 for starts 10–17 days before onset; Fig. 5a). Across 17 SSWs
the mean of the nine further systems correlates with the observed response at
r = 0.23, against 0.52 on event-free dates (p = 0.91), and no single system does
significantly better than on event-free dates (Fig. 5b). These tests rule out
large event-specific skill, not small: with 17 events the smallest gain in
correlation detectable with 80% power is 0.43, and the interval for the nine-system gain is −0.29 [−0.76, +0.22].

What the forecasts get wrong is the size of the shift. If they missed a second
population of outcomes, the observed values would scatter into both tails of the
ensembles; in ECMWF they do not (4% in the two outer rank bins against 17%
expected). They fall on the negative-NAM side, and so they do in the nine further
systems, as registered before the confirmatory test was run: mean rank 0.40
against 0.54 on event-free dates (p = 0.005; Holm-adjusted over all primary tests, 0.03), lower than the event-free value in nine of ten systems
(Fig. 5c). The deficit is specific to SSWs. It holds against event-free dates
drawn from all winters rather than from winters without SSWs (0.53, p = 0.012),
it is not a general error that grows as the vortex weakens (the rank barely
depends on vortex strength across winter dates, and SSWs fall below that
relation; p = 0.028), and it survives treating events as the unit (mixed model,
−0.14 [−0.25, −0.04]), restricting to the 15 events covered by at least five
systems (p = 0.012), removing each system's own rank level (p = 0.012) and
leaving out any one event (largest p = 0.022). In probability terms, the systems
give a negative-NAM fortnight after SSWs a mean probability of 0.57, whereas it
occurred 82% of the time (on event-free dates, 0.41 against 0.29). The ECMWF
forecasts raise the probability of a downward outcome from 0.35 on ordinary dates
to 0.60 after SSWs, whereas the observed rate rises from 0.29 to 0.75: they
capture 64% of the observed shift from starts 2–9 days before onset and 26% from
starts 10–17 days before. Correcting the size of the shift by its
leave-one-event-out mean error did not improve the probabilistic scores
significantly (p = 0.16), a test with little power at 17 events.

The same error reaches the regional cold. Registered before its data were
assembled, the primary regional test asks whether observed northern-Eurasian
temperature after SSWs lies on the cold side of the ensembles more than on
event-free dates. It does: mean rank 0.39 against 0.53 across the nine further
systems (p = 0.014; Holm-adjusted over all primary tests 0.07, so not significant after
correction; mixed model −0.13 [−0.25, −0.01]; at most
0.036 with any one event left out), below the event-free value in all ten systems (Fig. 5d). The ensemble-mean anomaly after these SSWs is −0.3 K where the
observed is −1.3 K. High-latitude Europe points the same way (0.43 against 0.53,
p = 0.07); mid-latitude East Asia and North America, where the imposed SSW cools
less or warms, show no such error. On synthetic forecasts this test rejected 3% of datasets with no error and
detected a 40% under-prediction of every anomaly in 32% of trials, so the absence of a signal elsewhere is not evidence of
calibration.

Whether models under-represent the surface response to SSWs is disputed:
initialised at onset, one seasonal model overestimated it^{sigmond13}^; ECMWF
over-persists the negative North Atlantic Oscillation after weak-vortex
states^{kolstad20}^; subseasonal systems couple too strongly from the lower
stratosphere to the surface at short lags^{garfinkel25}^; and nudged forecasts of
the 2018 event reproduce the predictable part of its observed
response^{dai25}^. The fall from 64% to 26% with lead suggests that part of the
error might simply be a missed warming. It is not. Registered before the data
were retrieved, the test splits the starts 2–9 days before onset into those
whose ensemble reversed the 10 hPa, 60° N wind within three days of the observed
onset (125 starts, all 17 events) and those that did not (43 starts, 15 events).
Outcomes lie as low in the ensembles that caught the warming as in those that
missed it: mean rank 0.41 against 0.55 on matched event-free draws (p = 0.014; Holm-adjusted 0.07) and 0.41 against
0.55 (p = 0.029), a difference of +0.003 [−0.05, +0.05]; for
northern-Eurasian temperature 0.40 (p = 0.06) and 0.41 (p = 0.11). Forecasts
started 0–7 days after onset, with the observed warming in their initial state,
are closer to calibrated (0.50 against 0.56, p = 0.14; temperature 0.46 against
0.52, p = 0.18), consistent with nudged forecasts that reproduce the predictable
response once the observed stratosphere is imposed^{dai25}^. Before onset, then,
forecasts that predict an SSW still under-predict its surface consequence; our
reversal criterion does not measure how strong or persistent the predicted
warming is, so a warming forecast too weak and too weak a coupling to it cannot
be told apart here.

### Implications for forecasting and research

The downward label measures where a threshold cuts a shifted distribution. Under
an imposed SSW, members fall into both classes in proportions set by the shift,
the contrast between the classes is what one shifted population produces, the
event's own forced odds explain a few per cent of its label, and knowledge of the
event adds nothing detectable to a probabilistic forecast beyond the shift. None
of this denies that SSWs matter or that events differ: the shift is large and
causal, it roughly triples the chance of a cold northern-Eurasian fortnight, and
the post-onset lower stratosphere carries information, as it does on any winter
day.

Two consequences follow. For forecasting, the product after an SSW is the shifted
distribution — for the polar-cap NAM, the probability of a negative state, and
for regional weather, the shifted distribution of the regional variable itself,
which the polar-cap index does not fully carry — rather than a categorical
prediction of whether this event will propagate. The quantity to verify and
correct is the size of that shift, which operational systems under-predict, not
the event-to-event differences the class is meant to summarise, which, for the polar-cap circulation, none predicts significantly better than it
predicts differences between ordinary winter days. For research, wherever "downward-propagating SSW" is used to
stratify impacts, the rate at which events meet the criterion, set against a
matched null, is an honest quantity; the contrast between the classes is not,
because the threshold manufactures most of it (97% in ERA5; see also
ref.^{white19}^).

The results also say what should be verified. A system's value after an SSW is
not whether it named the class of the event, a class fixed in large part by the
outcome it is scored against, but whether its ensemble placed the right
probability on the shifted state. That can be checked routinely with the test
used here: rank the outcomes after SSWs within the ensembles and compare the
ranks with those on event-free dates of the same season and lead. Across ten current systems the ensembles started before onset are, over the
weeks after SSWs, not negative enough in the annular mode (significantly in their multi-model
mean, and as point estimates in nine of the ten systems) and too warm over
northern Eurasia (as point estimates in all ten, not significantly after
correction for multiple tests).

Five limits bound these conclusions. First, the decisive evidence is from models;
with 42 usable observed events the observational out-of-sample tests are
uninformative, and the operational tests rest on 12 to 17 events. Second, SNAPSI
contains two Northern Hemisphere events, so statements about how events differ
rest on n = 2 there and on the CMIP6 bound; the result that the class contrast is
what one shifted population produces does not depend on n, because it is measured
within each ensemble. Third, SNAPSI nudges only the zonal-mean stratosphere and
our CMIP6 fields are zonal means, so non-zonal vortex geometry, which is
associated with stronger responses^{nebel24}^ and differs between nudged
members^{feng25}^, is neither held fixed nor among our model predictors; we test
it only in reanalysis. Fourth, "not detectable" is bounded, not zero: after SSWs
the event-aware forecast's gain is at most 1.4% of the probabilistic score (1.9%
with the extended predictors; upper 95% bounds), and operational gains in correlation below about 0.4 could not be
detected. Fifth, the regional analysis covers mean temperature in four
literature-defined boxes, not precipitation, wind or local extremes, and observed
regional samples of 39 events can exclude only residuals larger than about a
quarter of a standard deviation.

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

**Causal shift and distribution (Fig. 1a; Extended Data Fig. 1).** S = mean(nudged) − mean(control) of
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

**Forecast value (Fig. 4a).** Two Gaussian forecasts of each event's response,
in the member's own σ units and not demeaned: the shifted distribution (mean from
a calendar regression, day-of-year sine and cosine, fitted on the SSWs of the
other models; spread its residual s.d.) and the event-aware forecast (ridge
regression on the predictors up to onset plus calendar; spread the out-of-fold
residual s.d. within the training models). Leave-one-model-out cross-validation;
continuous ranked probability score in closed form; skill 1 − CRPS(event-aware)/
CRPS(shift) with a 2,000-replicate bootstrap over models. The identical pipeline
is applied to 200 size-matched sets of pseudo-onsets (same members, onset counts
and calendar days).

**Predictability (Fig. 4b,c).** Ridge regression, five-fold cross-validation grouped
by member, all variables standardised within member, R² on pooled out-of-fold
predictions; null from 1,000 size-matched pseudo-onset sets; a 1,000-replicate
member-cluster bootstrap gives the second interval. In 1,022 of 20,000
member-draws some pseudo-onsets lack enough in-season outcome days, which lowers
the null and so errs against our conclusion. Each draw and replicate has its own
seeded stream.

**Operational reforecasts (Fig. 5).** S2S reforecasts^{vitart17}^ from the ECMWF
Data Store, one model version per system (members per start; hindcast years).
Discovery set: ECMWF (model year 2022; 11; 2002–2021; Monday and Thursday
starts). Confirmatory set: ECCC (2025; 4; 2001–2020), CMA (2022; 4; 2007–2021),
HMCR (2025; 11; 1991–2020), KMA (2026; 7; 1993–2016), CNRM (1 June 2025; 11;
1999–2024), JMA (30 September 2022; 5; 1991–2020), CNR-ISAC (16 October 2023; 8;
2001–2020), NCEP (1 March 2011; 4; 1999–2010) and CPTEC (4 January 2023; 11;
1999–2018). Starts from December to March (KMA January to March), leads to 34
days (ECMWF 42). BoM was excluded when the design was fixed, before the confirmatory data were
retrieved, and
UKMO and IAP-CAS provide no sea-level pressure. Polar-cap sea-level pressure at 00 UTC
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
an event's systems at each pseudo-onset (no candidates for four events). Three sensitivities were added after the primary result. Because a null draw averages
only the systems with starts at that pseudo-onset, a system's overall rank level
can enter events and null unequally; removing each system's mean rank over its
pseudo-onsets before averaging gives H1 p = 0.012 (rejection rate 0.02 under
equal skill). CPTEC's forecast climatology rests on one to three other years per
start date; without CPTEC, H1 p = 0.005 and D p = 0.90. Leaving out each
event in turn, H1 p is at most 0.022.

**Regional temperature (Fig. 3).** Regions from the literature: northern
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
secondary (Holm-adjusted p for northern Eurasia over the four regions, 0.054);
the system-centred and leave-one-event-out variants were added after the
primary result.
On 100 synthetic datasets per case with the real temperature layouts and
outcome, rejection rates for the rank test were 0.03 for equally skilful
forecasts, 1.00 for forecasts without information and 0.32 for forecasts that
under-predict every anomaly by 40% (noise scaled by the ratio of the anomaly s.d.
of temperature and pressure); the discrimination test rejected 10% of
equally skilful datasets, so it is not interpreted for temperature.

**Precursor coupling (Extended Data Fig. 3).** Week-2 (days 8–14) 100 hPa
polar-cap geopotential height against days 15–25 polar-cap sea-level pressure,
correlated across members within each ensemble; Fisher-z pooling weighted by
n − 3; nudged-minus-control difference on matched ensembles with a normal-theory
interval. A pre-specified gate required the nudged/control spread ratio at
100 hPa to exceed 0.5 (median 0.96).

**Archetypes (Extended Data Fig. 2).** ERA5 polar-cap sea-level pressure at exactly
the members' forecast times, placed in each model's nudged distribution with the
same control base and spread; the pair test compares the observed 2018-minus-2019
difference with all member pairings of the same model.

**Bounds on event differences.** From the CMIP6 forced between-event variance
(the variance of the response across events minus that across size-matched
event-free dates, corrected by the same difference before onset), an event's
forced shift is modelled as Gaussian around the mean shift and an outcome as that
shift plus the event-free noise; the forced probability of a negative outcome is
evaluated for 1,000,000 simulated events at the point estimate and at the upper
95% bound, and the share of a single event's label variance due to forced
differences is Var(P)/[P̄(1 − P̄)].

**Contrast diagnosis and Southern Hemisphere spread.** On the paired Northern
Hemisphere ensembles (24 for the condition-1 check and 23 for the empirical null,
where both classes keep at least three members), the nudged contrast is compared with (i) the Gaussian
expectation for condition 1 alone and (ii) an empirical null: the same pair's
control members with every 6-hourly value shifted by the imposed effect,
classified with conditions 1–2. Skewness and s.d. of member means are compared
between arms. The Southern Hemisphere variance ratio pools the eight s20190829
ensembles re-centred on their own means, with 10,000 resamples of members within
ensembles.

**Reforecast checks.** Against all-winter dates: pseudo-onsets from every
December–March date of the hindcast years more than 30 days from every
catalogued onset. Vortex state: the NCEP 60–90° N mean zonal wind at 10 hPa on
the (pseudo-)onset date; ranks regressed on it over all such dates, and the SSW
residual compared with the same residual on calendar-matched draws. Mixed model:
rank of each system on each date, with a random intercept per date and systems
as fixed effects, on the events and all event-free candidate dates (971 for the
polar cap, 1,495 for temperature). Reliability:
the probability of a negative outcome as the share of members with a negative
window mean, averaged over starts and systems; Brier decomposition with bins of
width 0.2. Discrimination power: the change in correlation relative to the mean
of the conditional null, with an event bootstrap, and the minimum detectable
change (95th percentile minus mean plus 0.84 s.d. of the conditional null).
Post-onset and SSW-hit tests: forecast leads 1–9 days were retrieved for
polar-cap sea-level pressure and 2 m temperature and joined to the leads used
above (same starts and members), and the zonal-mean zonal wind at 10 hPa, 60° N
for leads 1–15. Starts 0–7 days after a catalogued onset are compared with the
same calendar-window null. A start 2–9 days before onset is a hit if at least
half of its members have a negative 10 hPa, 60° N wind on some lead within ±3
days of the observed onset. Hits and misses are compared with a matched null
(the same events and systems, and in each system a random subset of the
pseudo-onset's starts of the same size), and with each other by an event
bootstrap. The matched null replaced a comparison with the all-starts null, which was too
liberal; the change was made after a code review had emulated the polar-cap
version of the test with the all-starts null (misses p = 0.0095 then, 0.029
now), so it is not independent of that result, and it raised the p values.

**CMIP6 robustness.** The forecast-value pipeline repeated with (i) the post-onset
stratosphere (days 0–30); (ii) a weak-vortex null, zone-free dates whose member
10 hPa, 60° N zonal wind lies in that member's lowest 15% of November–March
zone-free days, one per real onset within ±30 calendar days; (iii) models
weighted equally, and CanESM5 (ten of the 20 simulations) excluded; and (iv) an
extended predictor set adding zonal-wind tendencies at 10, 50 and 100 hPa (60° N
and cap), 10–100 hPa shear and 75° N–45° N wind difference at 10 hPa (wave-driving
proxies; the archive holds no eddy fluxes) and the pre-onset surface annular mode
as a tropospheric precursor. 200 draws per null.

**Vortex geometry.** Following ref.^{seviour13}^: NCEP–NCAR daily 10 hPa
geopotential height north of 20° N (1958–2024); vortex edge, the
December–March zonal-mean height at 60° N; two-dimensional moments of (edge − Z)
over the region Z < edge, on the polar-stereographic plane with planar area
weights; centroid latitude and aspect ratio. The paper's event rule (centroid
below 66° N or aspect ratio above 2.4 for at least 7 days, events at least 30 days
apart) gives 17 displaced and 13 split events over their 52 winters, against 17
and 18 with ERA-40/ERA-Interim; the coarser NCEP grid likely smooths splits. Each
catalogued SSW is classified in a window of ±10 days around its central date:
split if the aspect ratio exceeds 2.4 on any day, otherwise displaced if the
centroid falls below 66° N on any day (registered rule: of the 43 events, 30 split, 11 displaced and 2 unclassified;
27 and 10 of those with ERA5 outcomes are compared),
and, as a sensitivity, with the 7-day persistence (9 split, 13 displaced, 21 unclassified; 9 and 10
compared). Outcomes from ERA5:
the 1000 hPa NAM over days +8..+25 and +8..+52, the northern-Eurasian
temperature anomaly over days +8..+24 and the downward label; split minus
displaced with within-class bootstrap intervals and two-sided permutation p
(10,000); Fisher's exact test for the label. Area weights and the validation
period were corrected to the paper before the test was first run on the complete
record.

**Test register.** Extended Data Table 1 lists every test, whether it is primary,
secondary or a sensitivity, whether its design was committed to the repository
before its data existed, before its run, or with its result, the commit, and the
raw and Holm-adjusted p values (over the registered primary tests with p values).

**Use of AI tools.** An AI assistant (Anthropic's Claude, in Claude Code) was used
to write and review analysis code, retrieve and check literature metadata, and
draft text under the authors' direction. Every number in the paper is produced by
a script in the repository and was checked against that script's output.
[[AUTHORS: confirm or revise this statement.]]

**Reproducibility.** Every result is produced by one script writing one JSON;
random streams are seeded in the code, producer hashes are recorded and checked,
and the environment is pinned. Defects found during the work are logged
in the repository.

---

## References

[[REFERENCES]]

---

## Figure legends

**Fig. 1 | One shifted population, cut by a threshold.** **a**, SNAPSI, 36
Northern Hemisphere ensembles of nine models: members without an SSW (control,
grey) and with the observed SSW imposed (nudged, vermilion), each ensemble
re-centred on its own mean and the nudged members placed at the mean imposed
shift, in units of each ensemble's control spread; vertical line, the threshold
of the downward criterion's first condition. The shares quoted are of members
meeting both surface conditions. **b**, ERA5: share of the 39 observed SSWs
meeting the surface conditions of the criterion, against event-free dates as they
are and displaced by the measured shift (95% interval). **c**, The archetype pair:
forced probability of a downward outcome in each model (nudged members; lines
join models) for February 2018 and January 2019 at the primary initialisations;
vermilion, means over nine models. **d**, Probability of a fortnight (days 8–24)
colder than the 10th percentile of the no-SSW distribution over northern Eurasia
(50–65° N, 10–130° E): SNAPSI without and with the imposed SSW and as predicted
from the circulation shift alone (black ticks, without ECCC), and ERA5 on
event-free dates and after 39 observed SSWs (95% event-bootstrap interval).

**Fig. 2 | The downward label is a threshold on one population.** **a**, DW minus
NDW contrast with the SSW imposed (nudged) against the same model and
initialisation with no SSW (control), for the 25 pairs in which both arms form a
contrast; the dashed line is equality and the dotted lines mark the value a cut at
zero gives on a unit Gaussian, −2√(2/π). **b**, Fraction of members classified DW,
all 36 Northern Hemisphere ensembles per arm.

**Fig. 3 | Northern-Eurasian cold follows the circulation shift.** **a**,
SNAPSI: effect of the imposed SSW on days 8–24 temperature in four regions
(blue; 95% centre-bootstrap interval), the part implied by each member's
polar-cap circulation through the relation that holds without an SSW (grey), and
the residual (vermilion, with interval), in units of the control spread.
**b**, Probability of a fortnight colder than the control's 10th percentile: no SSW
(0.10 by construction), SSW imposed, and predicted from the circulation shift
alone (hatched). **c**, ERA5, 39 observed SSWs: the same decomposition in kelvin;
whiskers, the central 95% of the mean residual over sets of event-free dates
(one per event, same calendar period).

**Fig. 4 | Knowing the event adds no detectable forecast skill after an SSW.** **a**,
Continuous ranked probability skill of the event-aware forecast over the shifted
distribution, on held-out models, after SSWs (point, 95% model-bootstrap
interval) and on 200 size-matched sets of ordinary winter days (grey; mean and
2.5–97.5% range), for
predictors before onset (P1) and up to onset (P2). **b**,**c**, Within-model
out-of-sample R² of the surface response at SSWs (vermilion lines) and on 1,000
event-free sets with the same onset counts (grey histograms), for predictors
before onset (**b**) and including the post-onset stratosphere (**c**).

**Fig. 5 | Operational reforecasts after SSWs.** **a**, ECMWF reforecasts
started 2–9 days before 12 SSWs: ensemble-mean against observed days +8..+25
polar-cap sea-level-pressure response (sign reversed, so negative is the
downward, negative-NAM outcome); dotted, equality. **b**, Correlation across
events between the ensemble-mean and observed responses for each system and for
the mean of the nine confirmatory systems (diamond), against the central 95% of
event-free dates of the same calendar period and lead with matched spread of the
outcome (grey bars); numbers, events per system; open symbols, systems with fewer
than six events (CPTEC, 4), whose values are not informative. **c**, Mean rank (probability
integral transform) of the observed outcome within the ensemble, against the
central 95% of the same event-free dates without the spread restriction; 0.5, a
calibrated forecast. **d**, Regional 2 m temperature: mean rank
(probability integral transform) of the observed days 8–24 regional temperature
within the forecast ensembles, for the mean of the nine confirmatory systems
(diamonds) and each system (dots), against the central 95% of event-free dates
of the same calendar period and lead (grey); p, the registered one-sided test
that outcomes are colder than forecast more often than on event-free dates.

**Extended Data Fig. 1 | An imposed SSW shifts the surface distribution without splitting it.**
**a**, Causal polar-cap NAM shift (nudged minus control, days +8 to +25) for nine
SNAPSI models and four Northern Hemisphere initialisations, in units of each
model's control spread; bars, ±1.96 s.e.; black ticks, model means; grey band,
mean ± s.d. over models excluding ECCC. **b**, Members of all 36 ensembles, each
re-centred on its own mean and pooled: control (grey), nudged (vermilion, placed at the
mean shift) and the control translated by the same shift (dashed), in units of
each ensemble's control spread; the shift here is the mean over all 36 ensembles
including ECCC in these per-ensemble units, and so is larger than in **a**.

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

**Extended Data Fig. 3 | The lower-stratospheric precursor couples to the surface
as strongly with no SSW.** **a**, For each of 36 SNAPSI ensembles, the correlation
across members between week-2 (days 8–14) 100 hPa polar-cap geopotential height
and the days 15–25 polar-cap surface response, with the SSW imposed (nudged)
against the same model and initialisation with no SSW (control); open squares,
the initialisation that starts after onset. **b**, Fisher-pooled correlations with
95% intervals and the nudged-minus-control difference, for all matched ensembles
and excluding the short-lead initialisation.

**Extended Data Table 1 | Test register.** Every test in the paper: its role
(primary, secondary, discovery, sensitivity or diagnostic), how its design was
registered (committed to the repository before its data existed, before its run,
or with its result), the commit, the statistic, the raw p value and the p value
Holm-adjusted over the registered primary tests that report one. Built by
`10_tables/tableED1_tests.py` from the result files.

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
