# Supplementary Information

*Stratospheric warmings shift winter cold risk rather than creating two kinds of event.* Supplementary Notes hold text moved from the main text to keep it within the article length; every number is from `ssw-design-analysis/CONSOLIDATED_RESULTS.md` or the result JSONs it indexes. Supplementary Tables are built by `10_tables/tableS*.py` from `results/current/`.

## Supplementary Note 1 | Event differences: the bound, its assumption and its falsifier

How much the forced response differs between events can be estimated only under assumptions. In CMIP6 the forced variance between events,
after removing the variance already present before onset, is 0.019σ²
[−0.033, 0.072]: the forced probability of a negative outcome ranges across
events from about 0.62 to 0.87 (0.47 to 0.93 at the upper bound), accounting for
2% of the variance of a single event's label (at most 8%). Any one label is
therefore almost entirely the shift plus chance. The bound assumes that an
event's forced response is uncorrelated with the internal variability it adds
to. It also gives a falsifier: at its upper limit, forced differences could
produce at most about half the class contrast of the criterion in CMIP6 and two
thirds in observations, so, under the latent-effect model tested, classes producing those full contrasts would have to exceed the measured, assumption-dependent bound; subtypes with smaller effects are not excluded. An audit finds the falsifier robust to the
strongest damping of internal variability the experiment allows and to forcing
that depends on the pre-onset state; it fails only if forcing anti-correlated
with concurrent internal variability (ρ ≤ −0.26 to −0.40), which no available
data can test (Methods).

## Supplementary Note 2 | Archetypes, post-onset coupling and vortex geometry

The field's archetype pair shows the same (Fig. 1c): February 2018 and January
2019 received comparable forced shifts and similar forced probabilities of a
downward outcome (0.80 and 0.75 across nine models), and at the primary
initialisations the observed outcomes lie within the central 95% of the model
distributions in 17 of 18 model–event cases (Extended Data Fig. 2), although
they are conventionally taken as opposite kinds. The lower stratosphere after
onset does track how much events differ (week-2 100 hPa anomaly against the
later surface response, r = 0.85 across 18 ensembles of distinct
events^{loeffel26}^; Methods), but between members of one SNAPSI ensemble, where
the forcing is shared, the relation is weak and equally present with no SSW
(r = 0.11 against 0.09; Extended Data Fig. 3), and in CMIP6 the post-onset
stratosphere predicts the surface as well on ordinary winter days: it is a
general property of stratosphere–troposphere coupling, not a feature that sorts
SSWs into kinds. The difference most often proposed, vortex geometry, cannot be
tested in SNAPSI or our zonal-mean CMIP6 fields, so we test it in
reanalysis^{seviour13}^. Events classified as vortex splits tend towards a
stronger early surface response than displacements^{lehtonen16,nebel24}^, but
not significantly: −0.25σ [−0.61, +0.13] over days 8–25 (p = 0.28; 37 events;
Methods). With 27 splits and 10 displacements the interval still admits a split
response stronger by 0.6σ, as expected when a 1,000-year simulation needs more
than 50 events to separate the two^{maycock15}^.

## Supplementary Note 3 | Regional residuals beyond the polar-cap circulation

The polar-cap index is not the whole of the regional response: high-latitude
Europe and mid-latitude East Asia cool more than the circulation implies, and
mid-latitude North America warms (+0.49σ) where the circulation relation implies
cooling, with the same sign of residual in observations (+0.55 K, p = 0.09), in
line with the mixed North American signal during weak-vortex
states^{kretschmer18,huang21}^ (Methods). What an SSW changes regionally is
therefore a distribution of the regional variable itself.

## Supplementary Note 4 | Operational discrimination

Operational forecasts agree on event differences: in seven subseasonal systems,
pre-onset skill for the 500 hPa polar-cap response after 13 SSWs lay at the 91st
percentile of that at random winter dates, not significantly
higher^{nebel24}^. In ten systems of the S2S reforecast archive^{vitart17}^
(ECMWF as a discovery set, nine others confirmatory; Methods), forecasts started
2–9 days before onset track event-to-event differences in the surface response
no better than on event-free dates of the same season and lead (ECMWF r = 0.77
against 0.52, p = 0.11, 12 SSWs; the nine others 0.23 against 0.52, p = 0.91,
17 SSWs; Fig. 5a,b), which rules out large event-specific skill, not small
(Methods).

## Supplementary Note 5 | Why the operational bias did not recur

Exploratory checks added afterwards point to the events,
not the systems: the newer ECMWF cycle ranks the earlier SSWs much as the older
one did (paired difference +0.019 [−0.055, +0.076]); given the kind of outcome,
the periods agree (upward outcomes ranked 0.76 in 1998–2021 and 0.75 in 2023–24);
what changed is the mix, 13 of 17 earlier outcomes downward against one of four
since; and the signs of the 2023–24 outcomes are 3.8 times likelier under the
forecasts as issued than under the 1998–2021 bias (Methods). Whether models under-represent
the surface response to SSWs is disputed^{sigmond13,kolstad20,afargan24,dai25}^;
our results bound such a bias as one of the size of the shift, of uncertain
persistence, rather than of event-to-event differences.

## Supplementary Note 6 | Continuity at the wind reversal with ERA5 winds, 1940–2025

The NCEP test of the main text uses 114 weakening episodes of 1958–2025. Its local
regression-discontinuity estimates were unstable (four of 16 placebo cutoffs with p < 0.05).
Several NCEP episodes just above zero are reversals in other reanalyses, so measurement error sits
exactly at the cutoff. The ERA5 version (registered before the ERA5 data were retrieved) uses the
ERA5 10 hPa, 60° N wind from 1940. Its minima agree with NCEP's (r = 0.98 on 96 shared episodes),
but five lie on the other side of zero: December 1958 (NCEP −5.8, ERA5 +6.4 m s⁻¹), January 1968
(+1.1, −4.7), January 1977 (0.0, −2.6), March 1981 (+0.9, −0.7) and February 2002 (+1.8, −0.4).
The dose response replicates for all four outcomes. No outcome shows a step after correction. The
largest step, −0.73 on the Arctic Oscillation (raw p = 0.03; Holm 0.11), is the strongest hint of a
step at reversal anywhere in this work. In CMIP6, 5,726 episodes bound the step at +0.05σ
[−0.02, +0.13]. In ERA5, reversals follow larger decelerations (balance p = 0.002), so a step
estimated at zero partly reflects how much faster the wind fell. The local estimates again fail their
placebo cutoffs, so the observational conclusion rests on the global estimates, as registered.
All numbers are in Methods and in `results/current/5_mechanism/vortex_threshold_continuity.json`
(key `era5_1940_2025`).

## Supplementary Note 7 | The shift rule out of sample

The 39 events of 1959–2022 verify the SNAPSI probabilities, but only the predictor was out of sample.
Twelve SSWs of 1941–1958 in ERA5 provided a test whose outcomes had not been seen. The rule failed
it at the registered 10th-percentile threshold (0 of 12 against 0.32; p = 0.012).

The early events differ in three ways:
- **They are weaker reversals.** The median minimum wind is −3.8 m s⁻¹, against −8.5 m s⁻¹ in
  1959–2022.
- **The data are older.** ERA5's stratosphere in the 1940s rests on few soundings reaching 10 hPa.
- **There are only twelve of them.**

None of these explains the failure on its own. The dose response accounts for about a quarter of
the smaller cooling, and the failure persists from 1946, when the upper-air record improves.

Two readings follow, and with twelve events they cannot be separated:
- the single SNAPSI probability, which rests on two strong imposed events, overstates the risk
  after weak reversals;
- the 1941–1958 sample was unusually mild.

Either way, a per-SSW probability verified on one set of events did not transfer to another. The
forecast product should therefore scale with the size of the disturbance and be verified
prospectively. Full counts are in Supplementary Table 4.

## Supplementary Note 8 | Members that do and do not reverse the vortex from one initial state

The ten S2S reforecast systems give a comparison with SNAPSI's design, in which every member shares one initial state
and members either have an SSW or do not. Here the SSW is not imposed: some members reverse the 10 hPa wind on
their own. The test was registered before any split was computed.

Of 461 starts with at least two members of each kind, the onset comes early enough for the 8–25-day window in 40
(162 for the 8–20-day window). Results are in Methods.

The signature of a threshold on one population appears:
- **Rate:** members that reverse are more often downward.
- **Contrast:** their class contrast is unchanged.
- **Spread:** they are no more spread out.
- **Shape:** their lower tail moves more than their upper. The quantile shift is −0.58σ at the 5th percentile
  against +0.12σ at the 95th, the asymmetry also seen in CMIP6.

Two results are weaker than in SNAPSI:
- **The mean shift is small** (−0.15σ over the shorter window).
- **Within starts there is no dose relation and no step.** This matches SNAPSI's own finding that, between members
  sharing one forcing, the post-onset stratosphere predicts the surface only weakly.

Spontaneous reversals in these ensembles mostly come late in the forecast (only 40 of the 461 mixed starts have the median onset by lead 9). The comparison therefore adds
evidence on the form of the change (translation, unchanged contrast), not on its size.

## Supplementary Note 9 | Extended results: the label as a threshold (rates, contrasts, criterion sweep)

Applying the surface conditions of the Karpechko criterion to every member
splits each ensemble into downward (DW) and non-downward (NDW) members (adapted
to the length of the runs; in observations we also apply it as published;
Methods). What the SSW changes is the rate: 84% of members are DW with the SSW
imposed against 45% without it (Fig. 1a; 79% against 41% in the Southern
Hemisphere). The same holds in observations: the share of the 39 observed events
with ERA5 outcomes that meet the surface conditions, 69% (54% with all three
published conditions), is reproduced by event-free dates displaced by the
measured shift (75% [61, 87]; Fig. 1b), so "about two thirds" is what one
shifted population produces. In the Northern Hemisphere the contrast between the
classes does not change: −1.62σ with the SSW imposed and −1.59σ without it
(paired difference −0.03σ [−0.12, +0.05] on 25 centre–initialisation pairs; Fig.
2a), while the rate separates the arms in every model (Fig. 2b). Without an SSW
the contrast is what a cut at zero gives on one Gaussian population; with it,
members without an SSW that are shifted by the imposed effect and classified
identically reproduce it (Methods). In ERA5, event-free dates reproduce 97% of
the surface DW–NDW contrast of the published criterion (93% of that of an 850
hPa variant^{lu26}^). Nor does this depend on how the criterion is set: across
108 versions varying its window, level, thresholds and stratospheric condition,
the observed downward rates are those of shifted event-free dates (calibrated
joint test p = 0.85, and 0.75 with the stratospheric condition; two planted
classes 1.8σ apart are detected with probability 0.82), and those dates
reproduce 76–154% of the class contrast (Extended Data Fig. 4).

What has not been tested directly is whether the label, or the warming itself,
identifies a kind of event. Three questions are easily conflated: whether
downward and non-downward outcomes are distinct kinds of response or one
distribution cut by a threshold; how much events differ in the response they
force; and whether knowing the event improves a forecast. Events can differ
without forming classes. The downward label differs from the split–displacement
typology in a way that matters for testing it: it is defined on the surface
outcome itself, so its classes differ by construction, and a test needs an
intervention that fixes the forcing and a null that does not depend on how the
criterion is set. We also ask what the shift means for regional cold risk and
what forecasts get wrong. Attribution studies in the same experiment quantified
the stratosphere's role in individual extremes^{seviour26}^; our subject is the
classification. We use an experiment that imposes the same observed SSW on every
ensemble member in nine models (SNAPSI^{hitchcock22}^), 1,517 SSWs in 20 CMIP6
simulations of 10 models, 43 observed SSWs and 114 observed vortex weakenings,
and reforecasts and real-time forecasts of ten operational systems. Extended
Data Table 1 lists every test and whether it was fixed before its data.

## Supplementary Note 10 | Extended results: continuity at the wind reversal

The same reasoning applies one step earlier, without a model. An SSW is defined
by a threshold, the reversal of the zonal wind at 10 hPa, 60° N. If reversal
marked a distinct kind of event, surface outcomes would jump at zero among
otherwise similar vortex disruptions. In a test registered before any outcome
was related to the wind, across 114 vortex weakenings since 1958, 42 of them
SSWs, the outcome instead grows smoothly with how far the wind falls (Extended
Data Fig. 6): per 10 m s⁻¹ weaker minimum wind the Arctic Oscillation over days
8–52 falls by 0.36 [0.21, 0.51], northern Eurasia cools by 0.36 K [0.07, 0.66]
and a downward label becomes 0.08 [0.01, 0.14] more likely, with no step at
reversal for any outcome (Arctic Oscillation −0.28 [−0.99, +0.47]). With 114
episodes this excludes only steps about as large as the whole SSW effect, and
local discontinuity estimates fail their placebo checks (Methods). In CMIP6 the
same rule yields 5,726 episodes and a step of +0.05σ [−0.02, +0.13], whose upper
90% limit is 15% of the 0.8σ shift observed after SSWs, against a dose of 0.17σ
per 10 m s⁻¹. With ERA5 winds for 1940–2025 instead (144 weakenings, 54 of them
reversals), which place five NCEP episodes on the other side of zero, the dose
replicates (Arctic Oscillation 0.29 [0.15, 0.43] per 10 m s⁻¹; northern Eurasia
0.38 K) and no outcome shows a step after correction for the four tested; the
largest, on the Arctic Oscillation, is −0.73 [−1.51, −0.06] (Holm-adjusted p =
0.11), and reversals there also follow larger decelerations (Supplementary Note
6). The tests resolve no surface-response discontinuity at reversal within these limits; "SSW" behaves here as a threshold on a continuum, which does not deny that reversal may matter dynamically.

## Supplementary Note 11 | Extended results: the shift rule as a forecast

The shift provides a fixed, model-derived probabilistic benchmark; on these 39 events it verifies. In a test registered
before it was run, the probability of a cold northern-Eurasian period after
an SSW was taken from the nudged experiment alone (0.24, 0.32 and 0.45 for
periods below the event-free 5th, 10th and 20th percentiles) and verified,
untuned, on the 39 observed SSWs: the observed frequencies, 0.21, 0.26 and 0.38,
are consistent with those probabilities (p = 0.40–0.71) and not with climatology
(p = 0.0006–0.008), and the risk ratio rises with severity, as a shift implies
(observed 1.9, 2.6 and 4.1 from the 20th to the 5th percentile; Extended Data
Fig. 7). A negative annular mode, given 0.74 by the CMIP6 events, followed 30 of
43 observed SSWs (climatology 0.42). A user who protects whenever the risk
exceeds their cost/loss ratio recovers 60–75% of the value of a perfect forecast
at ratios near the climatological rate, and loses value between the observed
rate and the slightly higher model probability (Methods). Only the predictor is
out of sample: the observed frequencies had been reported before. The
verification holds without the two SSWs that SNAPSI imposes (37 events: p =
0.38–0.57 under the rule, 0.002–0.012 under climatology) and at every threshold
from the 33rd to the 2.5th percentile, and it is bracketed by a second rule
built independently from the CMIP6 circulation shift and the ERA5 relation
between circulation and temperature on event-free dates (0.11, 0.19 and 0.33),
which under-predicts where SNAPSI slightly over-predicts (Supplementary Tables 1
and 3). The product has limits that it must carry: over high-latitude Europe the
models overstate the most extreme cold (none of 39 events below the event-free
2.5th percentile, against 20% predicted); over mid-latitude East Asia the
observed rise is not distinguishable from climatology with 39 events, and over
North America neither models nor observations show raised cold risk; and in
weeks 3–4, where SNAPSI rests mainly on the later initialisations of its two
events, northern-Eurasian cold was rarer than predicted (0.15 against 0.27) and
not distinguishable from climatology (Supplementary Table 2). Out of sample the
rule did not verify: after twelve SSWs of 1941–1958 in ERA5, whose surface
outcomes no one had examined, no northern-Eurasian period fell below the
event-free 10th percentile, against 32% predicted (p = 0.012), and cold risk did
not rise detectably above that period's climatology (Supplementary Table 4).
These reversals were weaker (an exploratory check; median minimum wind −3.8
against −8.5 m s⁻¹ in 1959–2022), but by the dose response that explains only
about a quarter of their smaller mean cooling (0.2 against 1.0 K), and 1940s
stratospheric winds are weakly constrained; a single probability per SSW
therefore holds for the events it was verified on and over-predicted for this
earlier, weaker set.

## Supplementary Note 12 | Extended results: event information and operational skill

If the event information available at onset determined how strongly an event couples, adding it should improve a forecast of that coupling; limited gains with these predictors and methods do not bound physical predictability. We compare two probabilistic forecasts of
each CMIP6 event's surface response (days 8–52): the shifted distribution alone
(SSW responses in other models, by calendar date) and an event-aware forecast
that adds the stratospheric state up to onset. Scored with the continuous ranked
probability score on events of held-out models, the event-aware forecast
improves on the shift by 0.8% [−0.3, 1.4] (Fig. 4a). The same information
improves forecasts on size-matched ordinary winter days by 3.7% [2.3, 5.3], and
in none of 200 such samples was it worth as little as after SSWs; the
two-component model gains 1.1% [0.4, 2.0], against 3.1% on ordinary days (a
comparison added afterwards). The result holds in five robustness checks
(Methods). The realised label cannot be an issue-time predictor, because it is defined from the later outcome (a model can still predict its probability); the test is of the event information available at onset. Point prediction agrees: the event adds −0.006 [−0.033, +0.020] to
the out-of-sample R² of the surface response before onset (Fig. 4b) and +0.015
[−0.028, +0.059] with the post-onset stratosphere, when 96% of that diagnostic
skill is present without an SSW (Fig. 4c; Methods).

Operational forecasts agree on event differences^{nebel24}^: in ten systems of
the S2S reforecast archive^{vitart17}^, forecasts started 2–9 days before onset
track event-to-event differences in the surface response no better than on
event-free dates of the same season and lead (Fig. 5a,b; Supplementary Note 4),
which rules out large event-specific skill, not small.

## Supplementary Note 13 | Extended results: operational forecasts in 1998–2021 and after

Over 1998–2021, what the forecasts got wrong was the size of the shift. In ECMWF, 4% of observed outcomes fell in the two outer rank bins against 17% expected (a non-symmetric rank histogram does not by itself exclude omitted structure). They fall on the negative-NAM side, as registered before
the confirmatory test was run: mean rank 0.40 against 0.54 on event-free dates
(p = 0.005; Holm-adjusted over all primary tests, 0.03), lower than the
event-free value in nine of ten systems (Fig. 5c), a deficit specific to SSWs
that survives the further checks in Methods. The systems gave a negative-NAM
period after these SSWs a mean probability of 0.57, whereas it occurred 82%
of the time, and ensemble-mean polar-cap responses were about half those
observed (116 against 241 Pa). Observed northern-Eurasian temperature lay on the
cold side in all ten systems, not significantly after correction (mean rank 0.39
against 0.53, p = 0.014; Holm-adjusted 0.07; Fig. 5d). The error was not merely
a missed warming: starts whose ensembles reversed the 10 hPa wind near the
observed onset under-predicted as much as those that did not (difference +0.003
[−0.05, +0.05]). A follow-up registered before the 100 hPa forecasts were
retrieved places the loss in the lower stratosphere: the forecast 100 hPa
anomaly was too weak (p = 0.005) while the 10 hPa wind was not, and given the
forecast 100 hPa state the surface was not detectably under-predicted (p = 0.30;
Methods) — the lower-stratospheric persistence deficit reported for nearly all
systems^{garfinkel25}^, not a coupling deficit.

The bias did not recur. After the three SSWs of 2023 and 2024, held out and
registered before they were scored, and after the SSW of 4 March 2026, scored
with real-time forecasts of five systems under a design registered before any
2026 forecast was retrieved, each observed polar-cap outcome lay, on average
over systems, on the upward side of the ensembles (four events: mean rank 0.70
against 0.55, p = 0.93; Methods, Extended Data Fig. 5). Exploratory checks added
afterwards point to the events, not the systems: given the kind of outcome, the
periods agree, and what changed is the mix, 13 of 17 earlier outcomes downward
against one of four since (Supplementary Note 5). Whether models under-represent
the surface response to SSWs is disputed^{sigmond13,kolstad20,afargan24,dai25}^;
our results bound such a bias as one of the size of the shift, of uncertain
persistence, rather than of event-to-event differences.

## Supplementary Note 14 | Is the imposed SSW a translation? Residual quantiles against a declared tolerance

A pooled variance ratio near 1 and a non-rejected shape test do not show that a distribution is translated, so the
translation was tested directly. In every nudged/control pair the nudged quantiles were compared with the control
quantiles moved by the difference in means, against a tolerance fixed before the test: ±0.25σ in any quantile and
±0.05 in a tail probability.

Results:
- **Polar-cap annular mode, all 36 ensembles pooled:** every interval falls inside the tolerance. This is a
  positive equivalence result, not a failure to reject.
- **Single events:** unresolved. In February 2018 the lower tail is compressed: the most negative states are less
  extreme than a pure translation implies.
- **Northern-Eurasian temperature:** unresolved. There is a small excess of extreme cold below the 5th
  percentile, inside the tolerance, and lower-quantile intervals that cross it.

So translation is established for the circulation index at the resolution of these ensembles. For regional
temperature it is neither established nor refuted. Figures are in Methods; values per pair are in
`results/current/8_experiment/snapsi_residual_quantiles.json`.

## Supplementary Note 15 | Does the shifted null transport? Whole-winter cross-validation

The shift used to build the null had been estimated on the same events it was compared with. Here it is
estimated without each held-out winter; climatology and the constant class rate are also estimated without that
winter. The test is retrospective, because the 39 outcomes had been seen, but nothing about a held-out winter
enters its own prediction.

Results:
- **Calibrated label rates out of fold:** 27 downward events observed against 29.2 expected; no version of the
  criterion is miscalibrated.
- **Better than climatology for the outcome itself** (CRPS).
- **No better from a constant class rate:** the label carries no information beyond the shift, as one shifted
  population predicts.

Values are in Methods and in `results/current/9_literature/criterion_transport_cv.json`.

## Supplementary Note 16 | How much do event-conditioned responses differ? The 18 ICON ensembles

Loeffel et al. (2026) re-ran 18 SSWs from a long ICON simulation as 40-member spin-off ensembles. Their public
archive gives the ensemble mean and spread of the standardised polar-cap geopotential anomaly. It does not give
member trajectories, so only the differences between the ensemble-mean responses can be assessed. Distribution
shape and label rates cannot.

The test was registered before any value was read. Onsets are where the ensemble-mean 10 hPa wind first turns
negative. The variance of the ensemble means between events, net of finite-ensemble noise, is bounded by
assuming the window-mean spread lies between zero and the daily spread.

Results:
- **Days 8–25:** between-event variance 0.41–0.43, about 80% of the within-ensemble variance.
- **Days 26–42:** 0.10–0.12, about 13–15%.
- **Link to the lower stratosphere:** the event means track the week-2 100 hPa anomaly (r = 0.80).

These conditional means differ greatly from the CMIP6 forced-variance estimate (0.019) because they estimate
different things:
- **They keep each event's initial tropospheric state,** which persists into the first weeks. The fall from days
  8–25 to days 26–42 is consistent with that. The CMIP6 estimate subtracts the differences present before onset.
- **The 18 events were selected,** so their spread is not a population variance.

What survives is the distinction the main text now draws: realised responses to different SSWs differ
substantially and in proportion to the lower-stratospheric anomaly. How much of that the stratosphere itself
forces is not settled by these summaries or by the assumption-dependent CMIP6 bound.

Values are in `results/current/5_mechanism/icon_event_heterogeneity.json`.
