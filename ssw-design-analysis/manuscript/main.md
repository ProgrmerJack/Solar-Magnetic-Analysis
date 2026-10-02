# Stratospheric warmings shift winter cold risk rather than creating two kinds of event

*Draft for Nature Geoscience (Article). Every number is taken from
`ssw-design-analysis/CONSOLIDATED_RESULTS.md` or the result JSONs it indexes;
display items are built by `09_figures/fig*.py` and `10_tables/table*.py` from
`results/current/`.*

---

## Abstract

Earth-system extremes are often understood by sorting them into kinds. Sudden
stratospheric warmings (SSWs), which raise the risk of cold winter weather, are
sorted into the two thirds whose influence propagates to the surface and the
rest, by a label defined from the outcome it is then used to explain. Here we
show that the label behaves as a threshold on one shifted distribution. Imposing
observed SSWs on nine models raises the share of downward outcomes from 45% to
84% without widening the northern annular mode; event-free dates, once shifted,
reproduce the observed two thirds; and in 1,517 simulated SSWs event information
adds about 1% to probabilistic skill. The warming is itself a threshold: across
114 observed vortex weakenings, the surface response grows smoothly with how far
the wind falls, with no resolved step where it reverses. Cold-risk probabilities
taken from the models alone match 39 SSWs of 1959–2022, after which cold 17-day
northern-Eurasian periods were 2.6 times as frequent as otherwise, but
over-predict after twelve earlier SSWs held out. Operational forecasts
under-predicted the shift after SSWs of 1998–2021, losing it in the lower
stratosphere, but not after four later ones. Forecasts should issue the shifted
distribution, not a class.

---

## Main

Extreme events in the Earth system are often understood by sorting them into
kinds, such as El Niño flavours, weather regimes or climatic epochs. Whether the
kinds exist, or are thresholds placed on one continuous population, decides what
can be predicted and how forecasts should be judged; tested against a
one-population null, several typologies have proved consistent with a single
regime or with stochastic variability^{stephenson04,neukom19}^, or reproducible
as a continuum^{thual23}^. Sudden stratospheric warmings (SSWs) are a case with
direct stakes. Cold-air outbreaks over Europe and Asia damage health, energy
systems and transport, and the stratosphere is one of their few sources of
predictability weeks ahead^{domeisenbutler20}^; in the United Kingdom alone, the
cold weather after an SSW has been linked to about 620 additional deaths per
event^{charltonperez21}^, and in February 2021 more than 4.5 million people in
Texas lost power in a cold outbreak^{ferc21}^ whose attribution to the
stratosphere is still argued both ways^{cohen21,davis22}^. SSWs are followed by
weeks of anomalous surface weather that projects onto a negative Northern
Annular Mode (NAM)^{baldwin21}^; after weak-vortex states the coldest days over
northern Eurasia become about twice as frequent^{kretschmer18}^, and forecasts
initialised at an SSW are more skilful than forecasts initialised at other
times^{sigmond13}^.

Not every SSW is followed by the canonical response (over the British Isles and
central Europe only about 45% are followed by a cold, negative-NAO
period^{hall23}^), and the field has organised this variability into two
classes. Events are labelled by whether their signal "propagates downward", most
often with the criterion of Karpechko et al.^{karpechko17}^; about two thirds
are said to have a visible downward impact^{baldwin21}^. Differences between the
classes are attributed to the events' strength, morphology, wave forcing or
stratospheric persistence and are used as effect sizes and
stratifiers^{runde16,rao20,lu26}^; related typologies sort events by wave
absorption or reflection^{kodera16}^ and by the weather regime at
onset^{domeisen20}^. The label was introduced to describe outcomes, often in
probabilistic terms^{karpechko17}^; the question is whether it can carry the
causal and predictive weight it is given. Forecast systems are assessed on
whether they anticipate which events will propagate downward, and non-downward
outcomes are framed as potential forecast busts^{nebel24}^. The label is defined
from the surface outcome it is then used to explain, and a cut placed on a
continuous outcome makes cases on either side of it look different whether or
not two kinds exist^{altman06}^.

Parts of this picture are known to be fragile, and the field has dissolved one
dichotomy before. The warmings form a continuum, with no clear threshold between
major and minor or between split and displaced events^{coughlin09,maury16}^;
splits and displacements, distinct in the stratosphere, differ little at the
surface, and their differences reduce to the mean lower-stratospheric
wind^{charlton07,maycock15}^. Composite differences between downward and
non-downward events are partly built by the classification^{white19}^; in an
idealised model the later tropospheric response is set linearly by the strength
of the lower-stratospheric warming^{white20}^; large ensembles find little in
the pre-onset state that separates one event's surface response from
another's^{bett23}^; no distinct tropospheric configuration is systematically
linked to vortex decelerations^{gallego26}^; and operational systems predict
before onset which events will propagate to 100 hPa but not which will reach the
troposphere^{nebel24}^. Nudged experiments show that an imposed SSW drives a
negative NAM^{hitchcock14,hong26}^ and that the contrasting 2018 and 2019
outcomes owed much to the tropics^{knight20}^.

What has not been tested directly is whether the label, or the warming itself,
identifies a kind of event. Three questions are easily conflated: whether
downward and non-downward outcomes are distinct kinds of response or one
distribution cut by a threshold; how much events differ in the response they
force; and whether knowing the event improves a forecast. Events can differ
without forming classes. The downward label differs from the split–displacement
typology in a way that matters for testing it: it is defined on the surface
outcome itself, so its classes differ by construction. In one model, composites
of the two classes already differ before onset, back to 40 days
earlier^{white19}^, as a label defined on the later outcome selects on anomalies
that persist; a test therefore needs an intervention that fixes the forcing, and
a null that does not depend on how the criterion is set. Attribution studies in
the same experiment quantified the stratosphere's role in individual
extremes^{seviour26}^; our subject is the classification. We use an experiment
that imposes the same observed SSW on every ensemble member in nine models
(SNAPSI^{hitchcock22}^), 1,517 SSWs in 20 CMIP6 simulations of 10 models,
observed SSWs and vortex weakenings since 1940 in two reanalyses, and
reforecasts and real-time forecasts of ten operational systems; Extended Data
Table 1 lists every test and whether it was fixed before its data.

### One shifted population, cut by thresholds

In SNAPSI's `nudged` ensembles every member's zonal-mean stratosphere is relaxed
to an observed SSW, in `control` ensembles to climatology (Methods). We analyse
nine models (4,407 members), the February 2018 and January 2019 events from two
initialisations each, and the September 2019 Southern Hemisphere minor warming.
Averaged over days 8–25 after onset, the imposed SSW moves the polar-cap NAM by
−1.36σ of each model's control spread across the nine models (−1.11σ without
ECCC, an outlier; Extended Data Fig. 1).

The effect is a translation (Fig. 1a). Re-centring each ensemble on its own mean
and pooling 36 ensembles, the members with the SSW imposed (1,798) have the same
variance as those without it (1,805) within each ensemble: ratio 0.955 [0.873,
1.051], which excludes an added forced variance larger than about 0.05σ². Tested
directly against a tolerance declared before the test (±0.25σ in any quantile,
±0.05 in a tail probability), the nudged distribution equals the translated
control at every quantile from the 5th to the 95th, pooled over the 36
ensembles; for single events, and for northern-Eurasian temperature, the
intervals are too wide to decide (Supplementary Note 14). This concerns the
polar-cap annular mode; regional spread can change^{spaeth24}^. The Southern
Hemisphere minor warming is the exception: there the imposed warming also widens
the distribution (variance ratio 1.74 [1.41, 2.15]), in an event whose
tropospheric response arrived more than a month after onset, underestimated even
with nudging^{feng25}^.

Applying the surface conditions of the Karpechko criterion to every member
splits each ensemble into downward (DW) and non-downward (NDW) members
(Methods). What the SSW changes is the rate: 84% of members are DW with the SSW
imposed against 45% without it (Fig. 1a; 79% against 41% in the Southern
Hemisphere). The same holds in observations: the share of the 39 observed events
with ERA5 outcomes that meet the surface conditions, 69% (54% with all three
published conditions), is reproduced by event-free dates displaced by the
measured shift (75% [61, 87]; Fig. 1b), so "about two thirds" is what one
shifted population produces. The shifted null also transports: with the shift
estimated from other winters only, it predicts 29.2 downward events where 27
occurred (whole-winter cross-validation; calibration p = 0.52, and p > 0.05 for
all 108 criterion versions), and a constant class rate predicts the labels no
better (Supplementary Note 15). In the Northern Hemisphere the contrast between
the classes does not change: −1.62σ with the SSW imposed and −1.59σ without it
(paired difference −0.03σ [−0.12, +0.05]; Fig. 2a), while the rate separates the
arms in every model (Fig. 2b); members without an SSW, shifted by the imposed
effect and classified identically, reproduce the contrast (Methods). A
within-start comparison in ten operational systems, which is associational
rather than an intervention (a member's own tropospheric waves can drive both
its reversal and its later state), shows the same signature: ensemble members
that reverse the vortex are more often downward than members started from the
same initial state that do not (0.54 against 0.45), with an unchanged class
contrast (difference −0.004σ [−0.33, +0.30]) and no added spread (variance ratio
0.97 [0.73, 1.36]; 0.96 [0.82, 1.18] over a shorter window with 162 starts),
although their mean shift is small (Supplementary Note 8). In ERA5, event-free
dates reproduce 97% of the surface DW–NDW contrast of the published criterion.
Nor does this depend on how the criterion is set: across 108 versions varying
its window, level, thresholds and stratospheric condition, the observed downward
rates are those of shifted event-free dates (calibrated joint test p = 0.85, and
0.75 with the stratospheric condition; two planted classes 1.8σ apart are
detected with probability 0.82), and those dates reproduce 76–154% of the class
contrast (Extended Data Fig. 4).

The same reasoning applies one step earlier, without a model. An SSW is itself
defined by a threshold, the reversal of the zonal wind at 10 hPa, 60° N; if
reversal marked a distinct kind of event, surface outcomes would jump at zero
among otherwise similar vortex disruptions. In a test registered before any
outcome was related to the wind, across 114 vortex weakenings since 1958 (42
SSWs), the outcome instead grows steadily with how far the wind falls (Extended
Data Fig. 6): per 10 m s⁻¹ weaker minimum wind the Arctic Oscillation over days
8–52 falls by 0.36 [0.21, 0.51] and northern Eurasia cools by 0.36 K [0.07,
0.66], with no step at reversal for any outcome. ERA5 winds for 1940–2025 (144
weakenings), which place five NCEP episodes on the other side of zero, give the
same dose and no step after correction; the largest, on the Arctic Oscillation
(−0.73 [−1.51, −0.06], Holm-adjusted p = 0.11), coincides with larger
decelerations before reversals (Supplementary Note 6). The observational record
excludes only steps about as large as the whole SSW effect; in CMIP6, 5,726
episodes bound the step at reversal to at most 15% of it (Methods). The tests
thus resolve no surface-response discontinuity at the reversal that defines an
SSW, within these limits; like "downward", "SSW" behaves here as a threshold on
a continuum, which does not deny that reversal may matter dynamically.

### Events differ continuously; no kinds are resolved

How much the forced response differs between events can be estimated only under
assumptions. In CMIP6 the forced variance between events, after removing the
variance already present before onset, is 0.019σ² [−0.033, 0.072]: the forced
probability of a negative outcome ranges across events from about 0.62 to 0.87,
accounting in this estimate for 2% of the variance of a single event's label (at
most 8%): small as a share of one realisation, but a range wide enough to change
the best action for users whose cost/loss ratio lies within it. Under the
latent-effect model tested, classes that produced the criterion's full class
contrast would have to exceed this assumption-dependent bound (subtypes with
smaller effects are not excluded), and an audit finds that statement robust
unless forcing is anti-correlated with concurrent internal variability (Methods;
Supplementary Note 1). Event-conditioned responses differ much more. In 18 model
SSWs re-run as 40-member ensembles^{loeffel26}^, the ensemble-mean polar-cap
response over days 8–25 varies between events with a variance of 0.41–0.43 (in
units of the daily climatological s.d., net of finite-ensemble noise), about 80%
of the variance within an ensemble, falling to 0.10–0.12 over days 26–42, and it
tracks the week-2 100 hPa anomaly (r = 0.80 [0.63, 0.91]; Supplementary Note
16). Those means also carry each event's tropospheric state at onset, which the
CMIP6 estimate removes, and the 18 events were selected, so the two numbers
answer different questions: realised responses to different SSWs differ
substantially, and how much of that difference the stratosphere forces is not
settled by either.

If the event information available at onset determined how strongly an event
couples, adding it should improve a forecast of that coupling. Scored on events
of held-out CMIP6 models, a forecast that adds the stratospheric state up to
onset improves on the shifted distribution alone by 0.8% [−0.3, 1.4] of the
continuous ranked probability score (Fig. 4a), against 3.7% [2.3, 5.3] on
size-matched ordinary winter days, and in none of 200 such samples was the
information worth as little as after SSWs; the result holds in five robustness
checks, and point prediction agrees (Fig. 4b,c; Methods). The realised label
cannot be an issue-time predictor, because it is defined from the later outcome
(a model can still predict its probability); the test is of the event
information available at onset, with these predictors and methods, and limited
gains do not bound physical predictability. Operational systems agree: forecasts
started 2–9 days before onset track event-to-event differences no better than on
event-free dates^{nebel24}^ (Fig. 5a,b; Supplementary Note 4).

Event differences show in shape as well as size, but not as kinds. A registered
comparison of statistical models finds that a two-component model with
state-dependent weights predicts CMIP6 responses after SSWs slightly better than
a continuous one, more so than on ordinary days (p = 0.01; Methods). Its
components lie 0.49σ apart in mean, closer than on ordinary days (0.71σ), and
the fitted density has one mode at every weight: it does not resolve
well-separated modes (a single regime in the sense of ref.^{stephenson04}^,
unimodal but not Gaussian), although statistical components need not correspond
to physical mechanisms and overlapping mechanisms are not excluded. Its lower,
"downward" component is also the wider one (0.72σ against 0.58σ), so the
predicted spread grows with the predicted shift, as when a forced shift varies
continuously in size between events^{white20}^, and the realised
lower-stratospheric dose is itself unimodal; whether that dose accounts for the
whole shape could not be tested with useful power (Methods). Among operational
ensemble members started from the same initial state, those that reverse the
vortex show no added spread either, and their lower tail moves more than their
upper (Supplementary Note 8).

The field's archetype pair shows the same (Fig. 1c): February 2018 and January
2019 received comparable forced shifts and similar forced probabilities of a
downward outcome (0.80 and 0.75 across nine models), and at the primary
initialisations the observed outcomes lie within the central 95% of the model
distributions in 17 of 18 model–event cases (Extended Data Fig. 2), although
they are conventionally taken as opposite kinds. Across events, the post-onset
lower stratosphere tracks the ensemble-mean response^{loeffel26}^; within one
SNAPSI ensemble, where the forcing is shared, the member-to-member relation is
weak and equally present without an SSW (Extended Data Fig. 3). These are
different relations (averaging removes noise from the first; the second asks
whether members differ by their own stratosphere), and in CMIP6 the post-onset
stratosphere predicts the surface as well on ordinary winter days: a general
property of stratosphere–troposphere coupling, not a feature that sorts SSWs
into kinds. Vortex splits, tested in reanalysis^{seviour13}^, tend towards a
stronger early surface response than displacements^{lehtonen16,nebel24}^, but
not significantly (−0.25σ [−0.61, +0.13]; 37 events), as expected when more than
50 events are needed to separate the two^{maycock15}^ (Supplementary Note 2).

### Consequences for cold risk and forecasts, and their limits

Over northern Eurasia (50–65° N, 10–130° E), the imposed SSW lowers the days
8–24 temperature by 0.85σ [0.49, 1.40] and raises the chance of a 17-day period
(days 8–24) colder than the control's 10th percentile from 0.10 to 0.32 (Fig.
1d, Fig. 3a,b), a risk ratio of 3.3. The shift of the polar-cap circulation,
acting through the relation between that circulation and regional temperature
that holds without any SSW, predicts 0.31: it reconstructs 88% of the cooling
and leaves −0.10σ [−0.36, +0.06]. This does not prove mediation: the regional
mean response is statistically consistent with that reconstruction, and other
pathways are not excluded. Without ECCC the cold-period probability is 0.26
against 0.27 predicted (Methods). Because nudged and control members share their
initial tropospheric states, pre-onset cold, often stronger than the cold after
onset^{lehtonen16}^, cannot enter these differences. At the earlier of two
initialisations before the February 2018 warming, imposing the observed
stratosphere roughly doubled the predicted risk of the ensuing Eurasian cold
spell relative to free-running forecasts^{seviour26}^. Within the SSW ensembles
the downward label adds nothing to a member's regional temperature beyond its
circulation (−0.04σ [−0.19, +0.09]). After 39 observed SSWs a cold
northern-Eurasian period occurred in 26% [13, 41] of cases, against 10% on
event-free dates (risk ratio 2.6); the same circulation–temperature relation
reconstructs two thirds of the 1.0 K cooling, and the remainder (−0.33 K) lies
within the range of event-free dates (p = 0.38; Fig. 3c).

The shift provides a fixed, model-derived probabilistic benchmark, and its
transfer is the test. In a test registered before it was run, the probability of
a cold northern-Eurasian period after an SSW, taken from the nudged
experiment alone (0.24, 0.32 and 0.45 for periods below the event-free 5th,
10th and 20th percentiles), matches the 39 observed SSWs of 1959–2022 (0.21,
0.26 and 0.38; p = 0.40–0.71) and not climatology (p = 0.0006–0.008), at every
threshold from the 33rd to the 2.5th percentile and without the two events
SNAPSI imposes; an independent rule built from the CMIP6 circulation shift and
the ERA5 event-free relation between circulation and temperature brackets the
observations from below (0.11, 0.19 and 0.33; Extended Data Fig. 7;
Supplementary Tables 1 and 3). Only the predictor is out of sample there, and
the product has limits that it must carry: it overstates extreme cold over
high-latitude Europe, is not confirmed over East Asia, carries no cold-risk
signal over North America (nor do observations), and in weeks 3–4
northern-Eurasian cold was rarer than predicted (Supplementary Table 2). Out of
sample it failed: after twelve SSWs of 1941–1958 in ERA5, whose surface outcomes
no one had examined, no northern-Eurasian period fell below the 10th
percentile, against 32% predicted (p = 0.012; Supplementary Table 4). Those
reversals were weaker and their data older, but neither explains the failure
(Supplementary Note 7); a single probability per SSW held on the events it was
checked on and did not transfer.

The polar-cap index is not the whole of the regional response: high-latitude
Europe and East Asia cool more than the circulation implies, and mid-latitude
North America warms where it implies cooling, in observations
too^{kretschmer18,huang21}^ (Supplementary Note 3). What an SSW changes
regionally is a distribution of the regional variable itself.

Operational forecasts show the same structure, and the same limit. After the 17
SSWs of 1998–2021, outcomes fell on the negative-NAM side of the ensembles of
ten systems (registered: mean rank 0.40 against 0.54 on event-free dates, p =
0.005, Holm-adjusted 0.03; Fig. 5c,d): the systems under-predicted the shift. A
follow-up registered before its data places the loss in the lower stratosphere:
the forecast 100 hPa anomaly was too weak (p = 0.005) while the 10 hPa wind was
not, and given it the surface was not detectably under-predicted (p = 0.30), the
persistence deficit reported for nearly all systems^{garfinkel25}^. But the bias
did not recur after four later SSWs scored as registered held-out tests
(2023–2026; mean rank 0.70 against 0.55, p = 0.93; Extended Data Fig. 5);
exploratory checks point to the mix of outcomes rather than to model changes
(Supplementary Note 5). Whether models under-represent the surface response to
SSWs is disputed^{sigmond13,kolstad20,afargan24,dai25}^; in this historical
sample and these diagnostics, forecast errors after SSWs concerned the size of
the shift rather than event-to-event differences, but neither a general
mechanism nor its persistence is established.

### Implications

The downward label measures where a threshold cuts a shifted distribution: the
shift sets the class proportions and the contrast, an event's own forced odds
explain a few per cent of its label, and knowing the event adds little to a
probabilistic forecast beyond the shift. The warming that defines an SSW behaves
the same way. None of this denies that SSWs matter, that events differ, or that
the label usefully describes what happened: the shift is large and causal, it
more than doubles the chance of a cold northern-Eurasian period, and the
post-onset lower stratosphere carries information, as on any winter day. The
lesson is general: a class defined on the outcome it is then said to predict
always shows a contrast in that outcome, and like regimes and epochs
elsewhere^{stephenson04,neukom19}^ it has to be tested against one shifted
population before it is read as a kind. Whether a particular SSW "caused" a
particular cold outbreak, as argued for Texas in 2021^{cohen21,davis22}^, is a
question of event-specific counterfactual attribution, distinct from
classification; for forecasting and risk the useful quantity is how much the
warming changed the odds, the sense in which forecasters already say that a
major warming "loads the dice"^{butler24}^.

For forecasting, the product after an SSW is the shifted distribution (for the
polar-cap NAM, the probability of a negative state; for regional weather, the
regional variable itself, which the polar-cap index does not fully carry), not a
categorical prediction of whether this event will propagate; such a product,
taken from models alone, matched 39 events and had value for decisions but
over-predicted on twelve held out, so it must be verified as events accrue. What
should be verified is whether an ensemble placed the right probability on the
shifted state, not whether it named a class fixed largely by the outcome it is
scored against. Two of the three 2023–24 events fail the criterion's surface
conditions and would count as non-downward outcomes, the kind framed as
potential forecast busts^{nebel24}^, yet the forecasts had given those outcomes
36% and 45% probability (an exploratory check), as a shifted distribution
should. The quantity to verify, and with more events to correct, is the size of
the shift and its persistence in the lower stratosphere. For research, wherever
"downward-propagating SSW" stratifies impacts, the rate at which events meet the
criterion, set against a matched null, is an honest quantity; the contrast
between the classes is not, because matched event-free dates reproduce 97% of
the selected surface composite contrast under the same procedure in ERA5 (see
also ref.^{white19}^); that is a property of the composite, not a measure of how
much causal influence is artefactual.

Five limits bound these conclusions. First, the causal evidence and the bound on
event differences come from models, which may share biases; observations are
decisive where no model is needed (event-free dates reproduce the observed class
contrast and, once shifted, the downward rate; surface outcomes vary smoothly
across the reversal), and the observed record is a typical 43-event draw from
the CMIP6 events on each of five statistics (combined p = 0.62), which shows
compatibility, not exchangeability. CMIP6 models reproduce the surface evolution
after SSWs particularly well over Siberia^{hall22}^, and the direction of a
shared bias is not obvious: the observed Eurasian cold may be inflated by
sampling^{kolstad22}^, while operational systems under-predicted the shift in
1998–2021. Second, SNAPSI contains two Northern Hemisphere events, so statements
about how events differ rest on n = 2 there and on the CMIP6 bound; the result
that the class contrast is what one shifted population produces does not depend
on n, because it is measured within each ensemble. Third, SNAPSI nudges only the
zonal-mean stratosphere and our CMIP6 fields are zonal means, so non-zonal
vortex geometry, associated with stronger responses^{nebel24}^ and varying
between nudged members^{feng25}^, is tested only in reanalysis. Fourth, event
information after SSWs is worth little, not nothing (at most 1.4–2.0% of the
probabilistic score, upper 95% bounds), and operational gains in correlation
below about 0.4 could not be detected; the operational tests rest on 4 to 20
events per system or set. Fifth, the regional analysis covers mean temperature
in four literature-defined boxes, not precipitation, wind or local extremes, and
regional samples of 39 events exclude only residuals larger than about a quarter
of a standard deviation.

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

**SNAPSI.** Nudging relaxes the zonal-mean stratosphere to the observed evolution
above 50 hPa, with none below 90 hPa^{hitchcock22}^. Nudged and control ensembles for CCCma, CNR-ISAC, ECCC, ECMWF, KMA,
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
spread pooled over its initialisations, reported with the NAM sign (ECCC −3.41σ, a statistical outlier; without it the
mean is −1.11σ, s.d. 0.27 across models); standard
errors from the two ensemble variances. For the distribution, members are
standardised by their own ensemble's control mean and s.d. and each ensemble is
re-centred on its own mean before pooling (pooling without re-centring adds the
between-ensemble spread of shifts to the nudged variance only). Variance-ratio
interval from 4,000 bootstrap resamples of members within ensembles.
Shape: the pooled re-centred members do not reject a pure translation
(Kolmogorov–Smirnov p = 0.37), nor does the per-event CMIP6 response (p = 0.35,
n = 1,517); at 43 observed events the same test almost never rejects even where a
model's response is wider (Observation–model compatibility), so it is not evidence
there. A mixture test on the CMIP6 events detects a two-to-one mixture of
populations 1σ apart with probability 0.85 but one 0.75σ apart with probability
0.005, so closer populations are not excluded this way.

**Classification and paired test (Fig. 2).** The published criterion^{karpechko17}^
classifies an event with the NAM at 1000 hPa (conditions 1–2) and 150 hPa
(condition 3) over days 8–52 after onset; the SNAPSI forecasts (45–60 days long)
do not reach day 52 after onset at most initialisations, so we use its conditions 1–2 over days +8..+25 on the member's NAM proxy −(polar-cap sea-level
pressure − control mean)/control s.d.: window mean negative and more than half of
6-hourly values negative. A contrast needs at least three members in each class,
which removes 11 of 36 nudged ensembles (DW rate 0.99) and no control ensemble;
the DW rate is reported for all ensembles (without ECCC, 82% against 45%). The paired test uses the 25 centre ×
initialisation pairs in which both arms form a contrast, with 10,000 bootstrap
resamples of the eight centres. The one-population expectation for each ensemble
is the contrast a cut at zero gives on a Gaussian with that ensemble's mean and
s.d., −s φ(a)[1/Φ(a) + 1/(1 − Φ(a))] with a = −m/s (−1.596s at m = 0). Without an
SSW the contrast equals that expectation within 0.002σ [−0.02, +0.02]; with it,
control members shifted by the imposed effect and classified identically reproduce
the nudged contrast (difference +0.04σ [−0.03, +0.10]; centre-cluster bootstrap).
The ERA5 contrast is that of the published criterion's surface index (1000 hPa
NAM, days 8–52).

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
predictions (before onset 0.054 at SSWs and 0.060 on event-free dates; with the
post-onset stratosphere 0.38 and 0.37); with 42 usable observed events the
same out-of-sample tests are uninformative in observations; null from 1,000 size-matched pseudo-onset sets; a 1,000-replicate
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
2–9 days (for ECMWF also 10–17 days; discrimination there p = 0.37) before a catalogued onset, and an event's
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
departure. The cold-period probability predicted from the shift averages, over
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
analysed; the analysis without ECCC was added after. Cold-period
probabilities with the SSW imposed against those predicted from the circulation
shift: high-latitude Europe 0.32 against 0.26, mid-latitude East Asia 0.27
against 0.20 (East Asia cools by 0.51σ, 60% of it reconstructed). Without ECCC, the
northern-Eurasian cooling is 0.60σ [0.44, 0.71] with a residual of +0.02σ
[−0.02, +0.07]; the
high-latitude European cold-period probability is 0.29 against 0.23
predicted, and mid-latitude East Asia cools by 0.41σ, 79% of it reconstructed
by the circulation. With all models, high-latitude Europe and mid-latitude East Asia cool more
than the circulation implies (by 0.19–0.27σ, and by 0.21σ [−0.50, +0.03]); over
mid-latitude North America downward members are warmer than their circulation
implies (+0.29σ [+0.06, +0.53]).
Operational: daily-mean 2 m temperature ("2t", averaged over each 24 h lead
window; the definition of the daily mean differs between centres) for the same
ten systems, model versions and starts, over the same regions from the S2S
reforecast archive (1.5° grid), and the ERA5 regional series above as the
observation; anomalies, event and pseudo-onset definitions, null and statistics
as for the polar-cap tests, with the outcome the mean over days +8..+24 and the
sign chosen so that a low rank means colder than forecast. The primary test (rank
in northern Eurasia, multi-model mean of the nine confirmatory systems) was
registered before the regional data were analysed; the other regions are
secondary (Holm-adjusted p for northern Eurasia over the four regions, 0.054;
mixed model −0.13 [−0.25, −0.01]; high-latitude Europe 0.43 against 0.53,
p = 0.07);
the system-centred and leave-one-event-out variants were added after the
primary result (with any one event left out, p is at most 0.036).
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
interval (difference in Fisher z, +0.02 [−0.04, +0.09]). A pre-specified gate
required the nudged/control spread ratio at 100 hPa to exceed 0.5 (median
0.96). The published correlation across 18
ensembles of distinct events^{loeffel26}^ (r = 0.85) survives a control for the
overlap of its predictor and outcome windows (p = 0.0005), although overlap
alone yields r = 0.49 [0.11, 0.77].

**Archetypes (Extended Data Fig. 2).** ERA5 polar-cap sea-level pressure at exactly
the members' forecast times, placed in each model's nudged distribution with the
same control base and spread; the pair test compares the observed 2018-minus-2019
difference with all member pairings of the same model. Without ECCC the
forced probabilities of a downward outcome are 0.78 (2018) and 0.71 (2019).

**Bounds on event differences.** From the CMIP6 forced between-event variance
(the variance of the response across events minus that across size-matched
event-free dates, corrected by the same difference before onset), an event's
forced shift is modelled as Gaussian around the mean shift and an outcome as that
shift plus the event-free noise; the forced probability of a negative outcome is
evaluated for 1,000,000 simulated events at the point estimate and at the upper
95% bound, and the share of a single event's label variance due to forced
differences is Var(P)/[P̄(1 − P̄)]. The excess variance equals Var(f) +
2 Cov(f, ε) for forced part f and internal part ε, so the bound assumes
Cov(f, ε) = 0. The ceiling on the class contrast is what a classifier that knew
every event's forced response exactly would produce by splitting Gaussian forced
responses at the criterion's downward fraction q, s_f φ(z_q)[1/q + 1/(1 − q)],
evaluated at the wider of the uncorrected and placebo-corrected upper 95%
bounds of s_f (0.342 and 0.269 in CMIP6; 0.596 and 0.681 in observations), and
compared with the contrast the criterion produces on the same outcome
(conditions 1–3 with the daily AO as the 1000 hPa index in observations, where
the uncorrected bound gives a ratio of 1.8; conditions 1–2 in CMIP6, whose fields
include no 150 hPa height). At the upper limit the criterion's class contrasts are 2.0 (CMIP6) and
1.6 (observations) times that ceiling. Combining the wider, uncorrected bound with
the strongest damping of internal variability the experiment allows (Identification
audit) raises the share of a single label's variance due to forced differences to
about a fifth; at that damping the falsifier is still met (ceilings 60% and 72% of
the contrasts).

**Contrast diagnosis and Southern Hemisphere spread.** On the paired Northern
Hemisphere ensembles (24 for the condition-1 check and 23 for the empirical null,
where both classes keep at least three members), the nudged contrast is compared with (i) the Gaussian
expectation for condition 1 alone and (ii) an empirical null: the same pair's
control members with every 6-hourly value shifted by the imposed effect,
classified with conditions 1–2. Skewness and s.d. of member means are compared
between arms. With the SSW imposed, the contrast is 0.19σ [0.11, 0.30] smaller
than a Gaussian cut on each ensemble's mean and spread predicts; half of that
shortfall comes from the criterion's second condition (a fraction of days
negative), and the empirical null reproduces the nudged contrast, so the
shortfall belongs to the Gaussian shortcut, not to the ensembles. The Southern
Hemisphere variance ratio pools the eight s20190829 ensembles re-centred on their
own means, with 10,000 resamples of members within ensembles. In that minor
warming, where the imposed warming also widens the distribution, the raw
contrast is larger (−1.96σ against −1.58σ without it), and against the threshold
expectation for each ensemble's own mean and spread it is again smaller (by
0.20σ, in six of seven ensembles).

**Reforecast checks.** Against all-winter dates: pseudo-onsets from every
December–March date of the hindcast years more than 30 days from every
catalogued onset (null mean rank 0.53, p = 0.012). Vortex state: the NCEP
60–90° N mean zonal wind at 10 hPa on the (pseudo-)onset date; ranks regressed on
it over all such dates, and the SSW residual compared with the same residual on
calendar-matched draws (the rank barely depends on vortex strength, and SSWs fall
below that relation, p = 0.028). Mixed model: rank of each system on each date,
with a random intercept per date and systems as fixed effects, on the events and
all event-free candidate dates (971 for the polar cap, 1,495 for temperature);
SSW coefficient −0.14 [−0.25, −0.04]. Restricted to the 15 events covered by at
least five confirmatory systems, H1 p = 0.012. Reliability:
the probability of a negative outcome as the share of members with a negative
window mean, averaged over starts and systems; Brier decomposition with bins of
width 0.2. ECMWF shift capture: forecast probability of a downward outcome 0.35 on ordinary
dates and 0.60 after SSWs, against observed rates 0.29 and 0.75 (starts 2–9 days
before onset). Discrimination power: the change in correlation relative to the mean
of the conditional null, with an event bootstrap (nine-system gain −0.29 [−0.76,
+0.22]), and the minimum detectable
change (95th percentile minus mean plus 0.84 s.d. of the conditional null).
Post-onset and SSW-hit tests: forecast leads 1–9 days were retrieved for
polar-cap sea-level pressure and 2 m temperature and joined to the leads used
above (same starts and members), and the zonal-mean zonal wind at 10 hPa, 60° N
for leads 1–15. Starts 0–7 days after a catalogued onset (northern-Eurasian temperature 0.46
against 0.52, p = 0.18) are compared with the
same calendar-window null. A start 2–9 days before onset is a hit if at least
half of its members have a negative 10 hPa, 60° N wind on some lead within ±3
days of the observed onset (125 starts covering all 17 events; misses, 43 starts
covering 15 events; for northern-Eurasian temperature, hits 0.40, p = 0.06, and
misses 0.41, p = 0.11). Hits and misses are compared with a matched null
(the same events and systems, and in each system a random subset of the
pseudo-onset's starts of the same size), and with each other by an event
bootstrap. The matched null replaced a comparison with the all-starts null, which was too
liberal; the change was made after a code review had emulated the polar-cap
version of the test with the all-starts null (misses p = 0.0095 then, 0.029
now), so it is not independent of that result, and it raised the p values. With 17 events the
smallest detectable gain in correlation is 0.43. Correcting the size of the shift by
its leave-one-event-out mean error (H2) did not improve the probabilistic scores
significantly (p = 0.16), a test with little power at 17 events. Relative to
ordinary dates, ECMWF's forecast polar-cap shift is 64% of the observed from starts
2–9 days before onset (−270 against −423 Pa) and 26% from starts 10–17 days before
(−113 against −428 Pa); in the downward rate, 54% and 24%. Forecasts started 0–7
days after onset, with the warming in their initial state, are closer to calibrated
for the polar cap (0.50 against 0.56, p = 0.14), consistent with nudged forecasts
that reproduce the predictable response^{dai25}^. A low mean rank could also come
from a skewed ensemble or one too narrow on its negative side; the forecast
probabilities and the ensemble-mean responses show that the centre of the forecast
distribution is displaced. Whether models under-represent the surface response to
SSWs is disputed: one seasonal model initialised at onset overestimated
it^{sigmond13}^; ECMWF over-persists the negative North Atlantic Oscillation after
weak-vortex states^{kolstad20}^; many subseasonal systems couple 100 to 850 hPa
polar-cap height too strongly at lags under a week (realistically in the
multi-model mean) while sustaining lower-stratospheric anomalies too
briefly^{garfinkel25}^; reforecasts may be overconfident in the canonical
storm-track response^{afargan24}^. Of the other regions, high-latitude Europe
pointed the same way as northern Eurasia; East Asia and North America did not.

**Held-out SSWs.** The three catalogued onsets after the end of the ERA5 series
used above (16 February 2023, 16 January 2024, 4 March 2024) were never part of
an operational test. Registered before their forecasts were retrieved: CNRM (the
confirmatory version, whose files hold these winters but had not been scored at
these dates), ECMWF model year 2025 (hindcast years 2005–2024, 11 members) and
CMA model year 2025 (2010–2024, 4 members, January–March starts). "Model year"
is the archive's reforecast-production label: ECMWF's changed cycle (CY47R3 in the
earlier set, CY49R1 here), whereas CMA's configuration is consistent with the same
BCC-CPS-S2Sv2 model in both sets (unconfirmed: the version date was not retained
in our files). Each system covers
every event with one to four starts 2–9 days before onset. Observations after the
end of our WeatherBench 2 extractions (April 2022 for the polar cap, January 2023
for temperature) are ARCO-ERA5 averaged onto the same 1.5° grid and reduced
identically; over the overlapping winters the daily polar-cap and regional series
agree with correlation ≥ 0.9998 and mean differences ≤ 0.04 K (0.01 Pa), which
are removed. Statistics, null and quorum rule as for the polar-cap and regional
tests; the registered reading treats a non-significant result as uninformative
unless the mean rank lies above the null mean, which it does for the primary
polar-cap outcome (and for the secondary temperature outcome). Secondary,
northern-Eurasian temperature: 0.57 against 0.51 (p = 0.64). January 2024 followed
a rapid vortex recovery^{leeSH25}^. CNRM, whose model
version is unchanged, also placed every held-out outcome above the ensemble median
(ranks 0.79, 0.71, 0.63).
Per event (polar cap): 0.68, 0.83 and 0.64. Pooling the 17 earlier events (nine
confirmatory systems, with the spliced observations, whose climatology now
includes 2022–2024) and the three held-out events, each against its own systems'
null: polar cap 0.45 against 0.54 (p = 0.032), northern Eurasia 0.41 against 0.53
(p = 0.031); the earlier events alone reproduce their result (0.40, p = 0.004).
Lee et al.^{leeSH25}^ describe the January 2024 event as without canonical
surface influence and the March 2024 event as weakly coupled to the troposphere;
those descriptions had been read before registration.

**Held-out diagnosis (Extended Data Fig. 5).** Exploratory and added after the
held-out result was seen. (i) The held-out code reproduces CNRM's earlier rank
exactly (0.4271). (ii) With the earlier design, ECMWF CY49R1 on the ten 1998–2021
SSWs in its hindcast years gives 0.42 against a null of 0.53 (p = 0.08), CY47R3
on the same events 0.40 (paired difference +0.019 [−0.055, +0.076], 10,000 event
resamples), and CMA model year 2025 on its five 0.24 (p = 0.002). (iii) Across all
20 events (multi-system means over the nine confirmatory systems for 1998–2021
and the three held-out systems for 2023–24) the rank averages 0.31 when the
observed outcome was negative (14 events) and 0.75 when positive (6); by period,
upward outcomes rank 0.76 (four events, 1998–2021) and 0.75 (two, held out), and
downward ones 0.29 (13) and 0.64 (one); 76% of the 1998–2021 outcomes
were negative against one of three held-out. (iv) Forecast probabilities of a
negative outcome (share of members) were 0.64, 0.55 and 0.85 for the held-out
events; a constant logit shift fitted to the earlier events (+1.02) makes the
observed held-out signs 3.8 times less likely than the forecasts as issued do.
(v) The held-out events lie within the 95% prediction interval of the 1998–2021
regression of observed on ensemble-mean response. (vi) Ranks on dates more than
30 days from any SSW vary coherently across systems from winter to winter
(correlations of winter means 0.55–0.87 between pairs of systems), so forecast
errors are shared within a winter and events of one winter are not independent.

**A fourth held-out SSW.** The project's detector, run on NCEP–NCAR u(10 hPa,
60° N) to 17 March 2026 (where that series ends) and continued with ERA5, finds no
SSW in winter 2024–25 and one onset in 2025–26, on 4 March 2026 (a late, short
reversal: easterly 4–10 March). A test was registered before any 2026 forecast
or 2026 surface observation was retrieved; the CPC Arctic Oscillation for 30–31
March 2026 was seen afterwards and is disclosed in the register. Forecasts are the
real-time ensembles of the S2S archive started 2–9 days before onset, each start
used only where a reforecast of the same system and model version exists for the
same start month-day: primary ECMWF (four starts, 101 members), ECCC (three, 21),
HMCR (one, 41), KMA (two, 8; all model year 2026) and NCEP (three, 16); JMA (one,
5) secondary; CMA, CNRM and CNR-ISAC had no matching start. Forecast anomalies are
the real-time values minus the mean of those reforecasts over their hindcast
years; observed anomalies use the same hindcast years, from ERA5 continued with
ARCO-ERA5 to April 2026 with the offsets above. The rank uses the real-time
members (its expectation under calibration is 0.5 for any ensemble size); the null
draws pseudo-onsets from the reforecasts, as before. The observed outcome was
upward (polar-cap NAM positive): mean rank 0.67 against a null of 0.61
[0.28, 0.93] (p = 0.60; with JMA 0.65, p = 0.60); secondary, northern Eurasia was
warmer than forecast in every system (0.93 against 0.56, p = 0.98). By the registered reading this lies above the null mean and is not
consistent with the 1998–2021 under-forecast; it is one event. All four held-out
events together: 0.70 against 0.55 [0.35, 0.75] (p = 0.93); temperature 0.66
against 0.53 (p = 0.86).

**Where the shift is lost.** Registered before any 100 hPa forecast was retrieved.
For the 17 events and starts of the confirmatory test, the polar-cap (60–90° N)
geopotential height at 100 hPa, days 8–25, from eight confirmatory systems
(CNRM's archived "100 hPa" fields have polar-cap values of 10.8–11.5 km, not about
16 km, and were excluded before the first run), against ERA5; sign reversed so that
low is the weak-vortex direction. L1: the rank of the observed 100 hPa anomaly in
the ensembles, 0.45 against 0.59 [0.49, 0.70] on matched event-free draws
(p = 0.005; Holm over L1 and L2, 0.01). L2: within each ensemble, the members'
surface anomaly regressed on their 100 hPa anomaly gives a forecast conditional on
the observed 100 hPa state; the rank of the observed surface outcome in it is 0.45
against 0.48 [0.37, 0.59] (p = 0.30). L4 (secondary): u(10 hPa, 60° N), 0.50
against 0.44 (p = 0.89). By the pre-set reading (L1 low, L2 near the null) the shift
is lost in the lower stratosphere, not in its coupling to the surface; L2 not
rejecting is absence of evidence for a coupling deficit, not proof of calibration.
ECMWF, the discovery system, agrees (L1 p = 0.008, L2 p = 0.30). Of the mean
ensemble-mean surface error (−131 Pa [−312, +59]), about half (−67 Pa [−167, +34])
passes through the 100 hPa error (within-ensemble slope 3.5 Pa per gpm); the
intervals are wide.

**CMIP6 robustness.** The forecast-value pipeline repeated with (i) the post-onset
stratosphere (days 0–30); (ii) a weak-vortex null, zone-free dates whose member
10 hPa, 60° N zonal wind lies in that member's lowest 15% of November–March
zone-free days, one per real onset within ±30 calendar days; (iii) models
weighted equally, and CanESM5 (ten of the 20 simulations) excluded; and (iv) an
extended predictor set adding zonal-wind tendencies at 10, 50 and 100 hPa (60° N
and cap), 10–100 hPa shear and 75° N–45° N wind difference at 10 hPa (wave-driving
proxies; the archive holds no eddy fluxes) and the pre-onset surface annular mode
as a tropospheric precursor. 200 draws per null. Skill after SSWs against
ordinary days: (i) 17.2% against 20.7%; (iii) 0.4% after SSWs with models
weighted equally and 0.6% without CanESM5 (the null was not recomputed without
it); (ii) 0.8% after SSWs against 1.6% on weak-vortex days (p = 0.94); (iv) 0.8% after SSWs against 4.8% on ordinary days and 2.9% on weak-vortex
days. The predictors span the same range after SSWs as on the comparison days
(s.d. ratios 1.02 and 1.03), so the gap is not an artefact of a restricted
range.

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
compared). Split minus displaced: −0.25σ [−0.61, +0.13] over days +8..+25
(p = 0.28), −0.40σ (p = 0.12) with the persistent classification, −0.02σ
(p = 0.92) over days +8..+52, and downward rates 0.78 against 0.60 (p = 0.41).
Outcomes from ERA5:
the 1000 hPa NAM over days +8..+25 and +8..+52, the northern-Eurasian
temperature anomaly over days +8..+24 and the downward label; split minus
displaced with within-class bootstrap intervals and two-sided permutation p
(10,000); Fisher's exact test for the label. Area weights and the validation
period were corrected to the paper before the test was first run on the complete
record. In a 1,000-year simulation the sign of the split–displacement difference
depends on the detection algorithm^{maycock15}^.

**Criterion sweep (Extended Data Fig. 4).** 108 versions of the criterion in ERA5
(39 events): classification window days 8 to 25, 35 or 52; conditions 1–2 on the
1000 or 850 hPa NAM; condition 1 threshold 0, −0.25 or −0.5σ; condition 2 share of
negative days above 0.5, 0.6 or 0.7; condition 3 (150 hPa NAM negative on more
than 70% of days) off or on. For each version, 400 day-of-year-matched sets of
event-free dates (more than 135 days from every catalogued onset) are displaced
by that version's measured shift of each level and classified identically. The
joint tests take, across the 54 versions with conditions 1–2, the maximum |z| and
the mean z of the observed rate against the shifted sets, and refer them to the
same statistics of each set against the others (leave-one-out moments); p = 0.965
and 0.50 (with condition 3: 0.87 and 0.48). A code review found this registered
reference too lenient: the events' own shift is measured on them, so their window
mean matches the shifted null by construction, and in simulation the test never
rejected under its null. The calibrated reference, added afterwards, replays the
whole procedure with each of 200 event-free sets playing the events (displaced by
the measured shift, its own shift then re-estimated): p = 0.85 for the maximum
|z| and 0.96 for the mean z (with condition 3, 0.75 and 0.995). With a planted
two-class structure (two thirds of events displaced by a further −0.6σ, one third
by +1.2σ, keeping the mean shift) the calibrated maximum-|z| test detects it in
82% of replays (37% with condition 3). Per-version 95% intervals of the shifted
sets are not informative on their own: under the null only 0.2% of versions fall
outside them. The same construction underlies the interval of Fig. 1b, which is
therefore a consistency check rather than a test of nominal size. Event-free dates reproduce 76–117% of the class contrast with
conditions 1–2 and 89–154% with condition 3; the published version gives 53.8%
observed against 54.1% [38.5, 69.3] shifted. In CMIP6 (27 versions, conditions
1–2, 1,517 events of which 1,513–1,516 are classified per version, 20 sets per
member), p = 0.50 and 0.85 (registered reference, not recalibrated). An event is
classified when at least 90% of its window days exist, a rule fixed after
registration and before the run; the design was registered before the sweep was
run.

**Observation–model compatibility.** Five statistics of the observational
pipeline (43 events, CPC AO, days 8–52): the pure-shift Kolmogorov–Smirnov
statistic, the variance ratio against event-free dates, the shift in units of
the event-free s.d., the rate of conditions 1–2 and the share of the class
contrast reproduced within the event-free pool. Each is ranked among 10,000 draws
of 43 of the 1,517 CMIP6 events (two-sided p 0.22, 0.40, 0.59, 1.00 and 0.74;
minimum-p combination calibrated on the draws, 0.62). Per-model draws, registered
alongside, are not interpretable for models with 53–152 events, because 43 of so
few events overlap too much to represent sampling spread. In a perfect-model
check (each model as reality, 1,000 draws of 43 events) the registered pass rule
(error rates within 1–10%, coverage at least 90%) failed in every model, through
tests that erred too rarely: the Kolmogorov–Smirnov test with an estimated shift
rejected in 0–0.1% of draws and the two-thirds test erred in 0–1.4%. Only CanESM5
(683 events) is large enough for 43-event draws to represent sampling spread;
there the KS test rejected in 0.1% of draws although the model's full population
has a variance ratio of 1.22, so at n = 43 it has almost no size or power, and
the variance-ratio interval covered the model's value in 91%. The two-thirds
test's low error rate is partly built in, because its shift is measured on the
events it tests (see Criterion sweep). Registered before it was run; the
implementation choices above (the pooled event-free reference, T3 in s.d. units)
were fixed after registration and before the run.

**Identification audit.** The forced-variance bound assumes Cov(f, ε) = 0. On
synthetic data built from the real pairs of pre- and post-onset event-free
anomalies (1,000 data sets per case), the placebo-corrected upper bound covers a
known s_f in 94.5–97.5% of data sets when the assumption holds; the uncorrected
bound covers it in 89.6–89.8% at 43 events, just below the registered 90%.
Forcing tied to the pre-onset state lowers coverage at 1,517 events (to 73% at
ρ = −0.5 and s_f = 0.4) but moves the estimate little (s_f² 0.127 against a true
0.160), because pre- and post-onset anomalies correlate weakly (0.07–0.12). A forced response that damps internal
variability, or that anti-correlates with the concurrent internal anomaly, hides
forced variance. Damping is bounded by experiment: within SNAPSI ensembles, where
the forcing is shared, the variance ratio is at least 0.873 (lower 95% bound).
At that damping, the wider upper bound on s_f becomes 0.78 (observations) and
0.42 (CMIP6), the class contrast is 1.38 and 1.67 times the ceiling, and forced
differences account for at most 22% and 18% of one label's variance. The ceiling
reaches the published contrast only if forcing anti-correlates with the
concurrent internal anomaly at ρ ≤ −0.30 (observations) or −0.40 (CMIP6), or
−0.26 and −0.36 together with that damping (combined case added after review);
through the pre-onset state it cannot (required ρ beyond −2.8). Transferring a
SNAPSI days 8–25 polar-cap ratio to annular-mode outcomes over days 8–52 is an
assumption. Registered before it was run; the concurrent-anomaly case was added
before the run but is not in the registered design.

**Model comparison.** CMIP6 events and predictors as for Fig. 4a (information up
to onset), leave-one-model-out. M0: Gaussian with calendar mean (the shift). M1:
Gaussian with ridge mean and log-spread linear in the predictors (penalty by
inner grouped cross-validation). M2: two Gaussian components with their own
means (shared calendar terms) and spreads, membership logistic in the
predictors, fitted by generalised EM from 20 starts. M2b: M2 with constant
membership. Scores: closed-form CRPS and ignorance. The registered test compares
G = score(M1) − score(M2) after SSWs with 200 size-matched sets of event-free
dates: G = 0.0012 against −0.0023 [−0.0053, +0.0003], p = 0.01 (ignorance
p = 0.01); M2b against M0, p = 0.18. Synthetic controls with the real predictors
detected planted regimes (two to one, 1σ apart) in 29–33% of data sets and
exceeded the event-free 95th percentile without regimes in 13–15%. Added after
the result, and labelled so: against those no-regime controls p = 0.02 (0.00 for
ignorance); a skewed one-population model (M1 with a constant skew-normal shape)
improves on M1 by 0.00015 and M2 still beats it beyond ordinary days (p = 0.02);
fitted to all SSW events, M2's components are 0.49σ apart in mean with spreads
0.58 and 0.72, against 0.71σ [0.58, 0.84] apart on event-free sets. In SNAPSI,
five-fold cross-validated two-component mixtures fit ensembles with the SSW
imposed worse than those without (difference −0.029 [−0.043, −0.014], centre
bootstrap), but a planted 1σ mixture was never detected, so that arm has no
power. The implementation choices above (EM details, penalty grid, control
construction) were fixed after registration and before the run. Also added afterwards: skill against M0 after SSWs, with a model-cluster
bootstrap, is 0.8% [−0.2, 1.4] for M1, 1.1% [0.4, 2.0] for M2 and 0.9% [0.1, 1.5]
for the skewed model, against 3.7%, 3.1% and 3.7% on the event-free sets (M2:
p = 0.015).

**Response shape and dose.** Registered before it was run. The fitted
two-component model (all 1,517 events, predictors up to onset) has component
means −0.66 and −0.16 with spreads 0.72 and 0.58: the lower ("downward")
component is the wider one, with mean weight 0.62, and the mixture has one mode
at every membership weight from 0.01 to 0.99. Quantiles 0.05–0.95 of the
responses after SSWs lie 0.40–0.63σ below those on the 200 size-matched
event-free sets; with each sample's mean removed, the lower tail is stretched
(−0.15σ [−0.18, −0.015] at the 5th percentile, model-cluster bootstrap) and the
rest is unchanged within its intervals. The registered comparison (G and its
event-free reference as above) was repeated with the realised 100 hPa zonal wind
over days 0–15 and 15–30 after onset (60° N and polar cap) added to the
predictors, on the events and on every event-free set alike. This is a mechanism
check, not a forecast: the dose is observed after onset. G after SSWs is
−0.0146 against −0.0109 [−0.0160, −0.0062] on the event-free sets (p = 0.92;
ignorance p = 0.78; with all post-onset levels p = 0.84); the continuous model's
skill over the shift is 13.8% after SSWs and 16.1% on event-free sets, the
two-component model's 10.3% and 13.2%. The pre-set rule's
criterion for "accounted for by the dose" is met, but the test has no power: the
registered synthetic controls, drawn as noise around the calendar mean without the
dose relation, are not comparable to the real reference with these predictors
(both rates 1.0), and a check added afterwards, with synthetic outcomes built on
the real dose relation (in-sample R² 0.32, residual s.d. 0.60), detects planted
two-to-one regimes 1σ apart in 8% of data sets against 5% expected without them.
By the registered rule for detection below 0.2, the test is underpowered and does
not show whether the dose accounts for the shape. The dose itself (polar-cap 100 hPa wind, days 0–30, after
SSWs) is unimodal: one Gaussian component is preferred to two by ΔBIC = 25, and
the fitted two-component density has one mode.

**Continuity at the wind reversal (Extended Data Fig. 6).** Registered before
any outcome was related to the wind. Units are vortex-weakening episodes in
NCEP–NCAR u(10 hPa, 60° N), 1958 to April 2025: a day of November–March that is
the minimum over ±20 days, preceded within 30 days by a fall of at least 15 m s⁻¹,
and followed before 30 April by 10 consecutive days above max(u_min, 0) (for
reversals, the final-warming rule of the CP07 definition). This yields 114 units,
42 with u_min < 0; the rule recovers every NCEP onset of the compendium. Running
variable u_min, cutoff 0, key date the minimum. Outcomes, days 8–52: CPC Arctic
Oscillation (primary), ERA5 1000 hPa NAM, northern-Eurasian 2 m temperature, and
the downward label (conditions 1–2 on the AO). Dose: ordinary least squares slope
on u_min with the pre-deceleration AO, day of year and era as covariates;
global jump: the coefficient of 1(u_min < 0) in the same regression with separate
slopes on each side; winter-cluster bootstrap intervals (10,000), Holm over the
four outcomes. Per 10 m s⁻¹ weaker minimum the AO falls by 0.36 [0.21, 0.51], the
NAM by 0.22 [0.15, 0.30], northern Eurasia cools by 0.36 K [0.07, 0.66] and the
probability of a downward label rises by 0.08 [0.01, 0.14]. Jumps at reversal: AO
−0.28 [−0.99, +0.47] (90%: [−0.87, +0.35]), NAM −0.18 [−0.55, +0.18],
temperature +0.09 K [−1.87, +2.01], label +0.12 [−0.21, +0.43] (Holm-adjusted p 1.0
for all). Local polynomial regression discontinuity (rdrobust^{calonico14}^;
triangular kernel, MSE-optimal bandwidth, robust bias-corrected intervals) and
local randomisation within ±2.5 and ±5 m s⁻¹ are also null, but four of 16
placebo-cutoff tests (±5, ±10 m s⁻¹) give p < 0.05, so by the registered
falsification rule the local estimates are not interpretable with these numbers
of units and the conclusion rests on the global estimates; placebo outcome and
covariate balance at 0 pass (deceleration size p = 0.07). In the 20 CMIP6
members, the same rule gives 5,726 units (1,483 reversals): dose 0.17σ [0.16,
0.19] per 10 m s⁻¹ and jump +0.05σ [−0.02, +0.13] (90% upper limit 0.12σ). On
1,000 draws of 68 consecutive model winters the observational jump test rejects
in 2.5% (no jump in the model) and detects a planted jump of −0.8σ, the size of
the whole SSW effect, in 78%.

**The shift rule as a forecast (Extended Data Fig. 7).** Registered before the
verification was run. Probabilities are taken from models only: for a
northern-Eurasian period (days 8–24) colder than the q-quantile of the
no-SSW distribution, the share of nudged SNAPSI members below their control
ensemble's q-quantile, averaged over the 36 centre–initialisation pairs (0.24
[0.15, 0.37], 0.32 [0.22, 0.47] and 0.45 [0.37, 0.57] for q = 0.05, 0.10, 0.20;
centre bootstrap; without ECCC 0.18, 0.26, 0.40); for a negative annular mode
(days 8–52), the share of the 1,517 CMIP6 events, 0.74. Verification, untuned, on
the 39 observed SSWs of the regional test (ERA5 anomalies from a smoothed
1959–2022 climatology, quantiles from event-free dates): 8, 10 and 15 cold
periods (0.21, 0.26, 0.38); two-sided binomial p under the rule 0.71, 0.40,
0.42 and under climatology 0.0006, 0.004, 0.008; Brier skill against climatology
+0.12, +0.09, +0.11. With the three 2023–24 events (42): p under the rule 0.59,
0.25, 0.22. Negative AO: 30 of 43 (0.70) against 0.74 (p = 0.49) and a
climatological 0.42 (p = 0.0003); Brier skill +0.26. Relative economic value
(cost–loss model, one decision per event) is positive for users whose cost/loss
ratio lies between the climatological rate and the observed rate after SSWs
(0.75 [0.23, 0.87], 0.68 [0.24, 0.84], 0.60 [0.17, 0.79] at a ratio equal to q),
negative between the observed rate and the rule's probability, where the rule
over-states the risk, and zero elsewhere. Only the predictor is out of sample: the
observed frequencies at q = 0.10 had been reported before. Prospectively, after
the 4 March 2026 SSW (run after its forecast test) northern Eurasia was 4.7 K
warmer than normal over days 8–24, not cold at any threshold, an outcome the rule
gave probability 0.55–0.76; one event does not test it.

**Shift rule: independence, severity, regions, lead and a second predictor
(Supplementary Tables 1–3).** Registered together, before these outcomes were
computed. Probabilities for every region, threshold and window are the SNAPSI
shares as above, from the pairs whose members fully cover the window (days 8–24:
36 pairs, nine centres; days 15–28: 22 pairs, nine centres, both later
initialisations at every centre plus both earlier ones of CCCma and UKMO; days
29–42: 12 pairs, nine centres, mainly the 8 January 2019 initialisation).
Observed anomalies and thresholds are as for the 39 events, with a November–May
series for the later windows. Without the two SSWs imposed in SNAPSI (12
February 2018, 2 January 2019) the counts are 7, 9 and 14 of 37. In northern
Eurasia the observed frequencies are 0.10, 0.21, 0.26, 0.38 and 0.51 at the
2.5th, 5th, 10th, 20th and 33rd percentiles against 0.18, 0.24, 0.32, 0.45 and
0.60 from SNAPSI (binomial p 0.30–0.71) and against climatology p =
0.0006–0.025. The second predictor averages, over the 1,517 CMIP6 events, the
ERA5 probability of a cold window given each event's days 8–24 annular-mode
anomaly (mean −0.49σ), from the regression of northern-Eurasian temperature on
the 1000 hPa NAM over 3,912 event-free windows (2.0 K per σ, residual s.d. 2.4
K); no observation after an SSW enters it. Its probabilities, 0.11, 0.19 and
0.33, are consistent with the observations (p = 0.07, 0.31, 0.50). Over
high-latitude Europe the counts at the 2.5th and 5th percentiles (0 and 3 of 39)
reject the rule (p < 0.001 and 0.014); at the 10th–33rd percentiles they are
consistent with it and reject climatology (p ≤ 0.004). Over days 15–28 the
northern-Eurasian counts (3, 6, 11) are below the rule (p = 0.07–0.11) and
consistent with climatology; over days 29–42 they are consistent with the rule,
and over high-latitude Europe they reject climatology (p = 0.004–0.045).

**Out-of-sample SSWs, 1941–1958 (Supplementary Table 4).** Registered with the
tests above, before any ERA5 surface value before 1959 was retrieved. SSWs are
detected with the project's CP07 detector on the ERA5 zonal-mean wind at 10 hPa,
60° N (daily mean of four synoptic hours, Copernicus Climate Data Store); on
1958–2024 it recovers all 43 catalogued onsets within three days, with one
addition (17 February 2002). Twelve onsets fall between November 1940 and March
1958; the last coincides with the catalogue's first event, whose outcome had not
been examined. Outcomes are ARCO-ERA5 northern-Eurasian 2 m temperature reduced
as above; because of the warming trend, anomalies are taken from the 1941–1958
smoothed day-of-year mean and thresholds from that period's event-free dates.
Counts below the 5th, 10th and 20th percentiles are 0, 0 and 2 of 12 against
0.24, 0.32 and 0.45 (binomial p = 0.08, 0.012, 0.08) and are consistent with the
period's climatology (p ≥ 0.62); from 1946, when ERA5's upper-air record
improves, 0, 0 and 2 of 9 (p = 0.13, 0.036, 0.20). A negative polar-cap annular
mode over days 8–52 followed 7 of 12 (CMIP6 0.74, p = 0.20; climatology 0.42, p
= 0.38). Pooled with the 39 later events, 8, 10 and 17 of 51 (p under the rule
0.19, 0.052, 0.092; under climatology 0.004, 0.03, 0.02). Added afterwards and
exploratory: the early reversals are weaker (median minimum wind −3.8 against
−8.5 m s⁻¹; 58% against 26% above −4 m s⁻¹), their mean northern-Eurasian
anomaly is −0.2 K against −1.0 K, and their polar-cap shift −0.36σ; the dose
slope of the continuity test accounts for about 0.2 K of that difference.

**Residual quantiles of the imposed SSW (Supplementary Note 14).** Registered
with a tolerance declared before it was run. In each Northern Hemisphere
nudged/control pair (36), the residual R(q) = Q_nudged(q) − Q_control(q) − Δ,
with Δ the difference in means, at q = 0.05–0.95, and the tail-probability error
E(q) = P(nudged < Q_control(q) + Δ) − q at q = 0.05, 0.10, 0.20; intervals from
2,000 two-stage resamples (centres, then members within arms, Δ re-estimated in
each). Tolerance: |R| ≤ 0.25σ (a 0.25σ shift moves the probability below a
Gaussian 10th percentile by 0.044) and |E| ≤ 0.05. Pooled polar-cap NAM: every
interval within the tolerance (largest residual −0.11σ [−0.22, +0.05] at the
95th percentile; tail errors within [−0.03, +0.02]), approximate translation by
the registered reading. February 2018 alone: unresolved, with a compressed lower
tail (R(0.05) = +0.24σ [+0.02, +0.38]); January 2019 alone: unresolved.
Northern-Eurasian temperature, days 8–24: unresolved (R(0.05) = −0.11σ [−0.45,
+0.13]; E(0.05) = +0.018 [+0.002, +0.046]).

**Transport of the shifted null (Supplementary Note 15).** Registered before it
was run; retrospective, because the outcomes had been examined. Whole-winter
cross-validation over the 39 ERA5 events: in each fold the event-free
candidates, the shift of every level and the constant class rate are estimated
from the other winters; each held-out event is predicted by training-winter
event-free days within ±10 days of its date, displaced by the training shift.
Published surface conditions: 27 downward events observed, 29.2 expected
(Poisson-binomial p = 0.52); with the 150 hPa condition 21 and 21.3; no version
of the 108 has p < 0.05. Brier score, constant class rate minus shifted null,
−0.003 [−0.024, +0.015] (winter resamples). CRPS of the days 8–52 NAM,
climatology minus shifted null, +0.19 [+0.07, +0.30]; mean PIT 0.50, 26% in the
outer deciles.

**Within-start ensemble comparison (Supplementary Note 8).** Registered before
any member split was computed. In every December–March reforecast start of the
ten S2S systems (about 4,500 starts), a member reverses if its zonal wind at 10
hPa, 60° N is westerly at leads 1–2 and first changes to easterly at lead 3–20,
and does not reverse if it stays westerly through lead 20. Starts with at least
two of each are compared (461); the anchor is the median onset lead of the
reversing members, and outcomes are polar-cap NAM-proxy means over 8–25 days
after it (40 starts whose window fits within lead 34) and, secondarily, 8–20
days (162 starts), in units of the system's s.d.; intervals resample hindcast
winters (2,000). Shift (reversing minus non-reversing, mean over starts): −0.05σ
[−0.29, +0.17] and −0.15σ [−0.24, −0.03]; pooled within-start variance ratio
0.97 [0.73, 1.36] and 0.96 [0.82, 1.18]; downward rate 0.54 against 0.45, class
contrast −0.875σ against −0.872σ (difference −0.004σ [−0.33, +0.30]); before
onset the groups do not differ (−0.01σ [−0.07, +0.04]). Across 11,481 members of
1,450 starts with start fixed effects, the outcome over 8–25 days after each
member's minimum wind shows neither a dose relation (−0.05σ per 10 m s⁻¹ [−0.10,
+0.01]) nor a step at reversal (−0.05σ [−0.19, +0.10]). Northern-Eurasian
temperature: shift −0.08σ [−0.34, +0.13], variance ratio 1.29 [0.79, 2.66]. With
few starts in the primary window the shift is imprecise; the comparison shares
SNAPSI's design (one initial state, SSW or not) but not its forcing, which here
arises spontaneously and is weaker.

**Continuity with ERA5 winds (Supplementary Note 6).** Registered with the tests
above. The episode rule, outcome windows, estimators, Holm family and
falsification rule are those of the NCEP test, applied to the ERA5 wind for
minima from November 1940 to March 2025, with a linear year term among the
covariates. Outcomes: a polar-cap annular-mode proxy (ERA5 polar-cap mean
sea-level pressure at 00 UTC, WeatherBench 2 from 1959 and ARCO-ERA5 otherwise,
anomaly from the smoothed day-of-year mean, sign reversed, in units of its
November–March daily s.d.), the CPC Arctic Oscillation (from 1950),
northern-Eurasian temperature with a linear trend (0.38 K per decade) removed,
and the downward label on the proxy. Because the CPC index begins in 1950, the
pre-deceleration covariate and placebo outcome use the proxy (an implementation
note written before the first run). 144 episodes, 54 reversals, 111 with
outcomes; NCEP and ERA5 minima agree with r = 0.98 on 96 common episodes, and
five lie on opposite sides of zero (December 1958, January 1968, January 1977,
March 1981, February 2002). Dose per 10 m s⁻¹: proxy 0.16σ [0.10, 0.22], Arctic
Oscillation 0.29 [0.15, 0.43], temperature 0.38 K [0.14, 0.62], downward label
−0.12 [−0.16, −0.06]. Global steps: −0.28 [−0.67, +0.11], −0.73 [−1.51, −0.06],
−0.66 K [−2.15, +0.77] and +0.26 [−0.06, +0.59]; Holm-adjusted p 0.34, 0.11,
0.39 and 0.34. The local estimates are null (Arctic Oscillation −0.02 [−1.40,
+1.58]) but five placebo cutoffs give p < 0.05, so they are not interpretable,
and deceleration size is unbalanced at zero (p = 0.002): reversals follow larger
falls in the wind. Results from 1946 are the same.

**Test register.** Extended Data Table 1 lists every test, whether it is primary,
secondary or a sensitivity, whether its design was committed to the repository
before its data existed, before its run, or with its result, the commit, and the
raw and Holm-adjusted p values (over the registered primary tests with p values).
Tests added in the second revision (robustness checks and the held-out
replication) each carry their own registered reading; they are not
added to that family, which would change the adjusted p values of tests
registered earlier.

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
vermilion, means over nine models. **d**, Probability of a period (days 8–24)
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
**b**, Probability of a period colder than the control's 10th percentile: no SSW
(0.10 by construction), SSW imposed, and predicted from the circulation shift
alone (hatched). **c**, ERA5, 39 observed SSWs: the same decomposition in kelvin;
whiskers, the central 95% of the mean residual over sets of event-free dates
(one per event, same calendar period).

**Fig. 4 | Knowing the event adds little forecast skill after an SSW.** **a**,
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

**Extended Data Fig. 4 | Every version of the criterion tracks one shifted
population.** **a**, ERA5: downward rate after 39 SSWs against that of event-free
dates displaced by each version's measured shift (mean and central 95% of 400
sets), for 54 versions with the surface conditions (filled) and 54 with the
150 hPa condition added (open); star, the published criterion; dashed, equality.
**b**, Share of the class contrast reproduced by event-free dates, per version,
sorted. **c**, CMIP6, 27 versions of the surface conditions, 1,517 events. Panel
titles give the calibrated joint-test p (surface versions; with the 150 hPa
condition in parentheses), a calibration added after review (Methods); the
intervals in **a** are not tests of nominal size.

**Extended Data Fig. 5 | Why the held-out replication failed (exploratory).**
**a**, Every event (17 of 1998–2021 with the nine confirmatory systems, grey; three
held out with their three systems, vermilion): mean rank of
the observed polar-cap outcome in the forecast ensembles against the observed
response (negative, downward); the rank follows the realised outcome. **b**,
ECMWF reforecasts of cycles CY47R3 (model year 2022) and CY49R1 (model year 2025)
on the same ten 1998–2021 SSWs; vermilion, means. **c**, Mean rank on dates more
than 30 days from every SSW, by winter, for the three held-out systems.

**Extended Data Fig. 6 | The warming itself is a threshold on a continuum.**
**a**, The 114 vortex-weakening episodes in NCEP–NCAR u(10 hPa, 60° N), 1958–2025
(filled, wind reversed: SSWs): CPC Arctic Oscillation over days 8–52 after the
wind minimum against the minimum wind; titles, the registered dose (per
10 m s⁻¹) and step at reversal with 95% winter-cluster bootstrap intervals.
**b**, The same for the northern-Eurasian 2 m temperature anomaly. **c**, Dose
(blue, per 10 m s⁻¹ stronger minimum wind; positive means a weaker vortex gives a
lower index or colder Eurasia) and step at reversal (vermilion), in s.d. of each
outcome across episodes, for three observed outcomes and for the annular mode in
5,726 CMIP6 episodes.

**Extended Data Fig. 7 | The shift alone as a forecast.** **a**, Observed
frequency after SSWs (Wilson 95% interval) against the probability taken from the
models alone (filled) and climatology (open): northern-Eurasian periods (days
8–24) below the event-free 5th, 10th and 20th percentiles (SNAPSI probabilities;
39 ERA5 events) and a negative annular mode over days 8–52 (CMIP6 probability;
CPC AO, 43 events); dashed, equality. **b**, Relative economic value of the rule
against climatology for users with a given cost/loss ratio: positive between the
climatological and the observed rate, negative between the observed rate and the
rule probability, zero elsewhere, where rule and climatology lead to the same
decision. **c**, Northern Eurasia, days 8–24: observed frequency of a colder
period (shading, Wilson 95% interval) against the event-free percentile that
defines it, with the SNAPSI rule, the independent rule from the CMIP6
circulation shift and the ERA5 event-free circulation–temperature relation, and
climatology (dashed).

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
- ERA5 after 2022 (held-out events, to April 2026): ARCO-ERA5,
  `gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3`
  (anonymous access).
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
- S2S real-time forecasts (the 2026 event): ECMWF Data Store, dataset
  `s2s-forecasts`, under the same terms.
- Daily Arctic Oscillation index: NOAA Climate Prediction Center.
- ERA5 zonal wind at 10 hPa, 60° N, 1940–2026: Copernicus Climate Data Store,
  `reanalysis-era5-pressure-levels` (licence CC BY 4.0).

## Code availability

All code, reduced data products, result files and run logs:
`https://github.com/ProgrmerJack/Solar-Magnetic-Analysis` (directory
`ssw-design-analysis/`).
