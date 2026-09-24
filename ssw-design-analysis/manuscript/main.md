# Surface effects of sudden stratospheric warmings are a common shift, not a property of the event

*Draft for Nature Geoscience (Article). Every number is taken from
`ssw-design-analysis/CONSOLIDATED_RESULTS.md`; display items are built by
`09_figures/fig*.py` from `results/current/`.*

---

## Abstract

Sudden stratospheric warmings (SSWs) are a major source of extended-range
predictability in the winter extratropics, and forecasts initialised at an SSW
are measurably more skilful. The field's account of that skill rests on a
classification: about two thirds of events are said to propagate downward to the
surface and the rest not, and the difference is read as a property of the event.
Here we test that reading in a designed experiment in which the stratosphere of
every ensemble member is nudged to the same observed SSW, across nine forecast
models. The SSW shifts the surface distribution by −1.11 σ (s.d. 0.27 across
eight models; one outlier) without changing its shape (variance ratio 0.955 [0.873, 1.051]).
Splitting members into "downward" and "non-downward" produces the same −1.6 σ
contrast whether the SSW is imposed or absent, the value a threshold yields on
pure noise. In 1,888 simulated events, an SSW adds no out-of-sample
predictability over an ordinary winter day (R² +0.004; 95% upper bound 0.05). The field's two archetype events receive the same
forced odds. The skill after an SSW is the common shift; the
classification adds nothing to it.

*(183 words)*

---

## Main

Sudden stratospheric warmings — the rapid breakdown of the winter polar vortex —
are often followed by weeks of anomalous surface weather, projecting onto a
negative Northern Annular Mode (NAM) (Baldwin et al. 2021). Forecasts initialised
at SSW onset reproduce the observed mean tropospheric conditions of the following
months and are more skilful than forecasts initialised at other times, for
circulation patterns, surface temperatures over northern Russia and eastern
Canada, and North Atlantic precipitation (Sigmond et al. 2013).

The dominant account of that predictability is event-specific. SSWs are divided
into those whose signal "propagates downward" to the surface and those whose
signal does not, commonly with the criterion of Karpechko et al. (2017), and
"about two thirds" of events are said to have a visible downward impact (Baldwin
et al. 2021). The difference between the two groups is then read as a property of
the events — of their strength, morphology or wave forcing — and recent work asks
which SSWs will couple downward, and how well forecast systems can tell in advance
(Rao et al. 2020; Nebel et al. 2024; Loeffel et al. 2026); non-downward outcomes
are framed as potential forecast busts (Nebel et al. 2024).

That reading has a known weakness. The classification is made on the surface
response itself, so a contrast between classes is partly guaranteed by the
selection: White et al. (2019) noted that differences between downward and
non-downward composites at positive lags are "entirely there by construction".
Others have argued that SSWs form a continuum rather than classes (Coughlin &
Gray 2009; Maury et al. 2016). What has been missing is a test that separates the
two readings — an experiment
in which the stratospheric forcing is identical by design, so that any difference
between classes cannot come from the event.

The Stratospheric Nudging And Predictable Surface Impacts (SNAPSI) experiment
provides one (Hitchcock et al. 2022). In its `nudged` ensembles the zonal-mean
stratosphere of every member is relaxed to the observed evolution of a real SSW
above 50 hPa, with no nudging below 90 hPa; in its `control` ensembles it is
relaxed to climatology. Members of a nudged ensemble therefore share the same
zonal-mean stratospheric evolution and differ in tropospheric weather. We analyse the complete archive
for nine models (4,407 members: 3,603 in the Northern Hemisphere, 804 in the
Southern; a tenth submission is excluded as corrupt), the two Northern
Hemisphere events it contains — February 2018 and January 2019, each
from two initialisations — and the September 2019 Southern Hemisphere minor
warming (eight models), together with 43 observed events and 1,888 events in 20
CMIP6 members of 11 models.

### An SSW shifts the surface distribution without splitting it

The causal surface effect of an SSW is the difference between nudged and control
ensembles initialised on the same date. Averaged over days 8–25 after onset, it
is a negative polar-cap NAM anomaly of −1.11 σ of each model's control spread
(s.d. 0.27 across eight models; Fig. 1a). One model, ECCC, responds three times
as strongly (−3.41 σ) and lies well outside the Tukey fence of the other eight;
including it gives −1.36 σ. Across the eight models the spread (s.d. 0.27) is
about twice the difference between the two events (0.14 σ): models disagree
about the strength of downward coupling more than these two SSWs differ.

The effect is a translation. Re-centring each ensemble on its own mean and pooling
across 36 ensembles, the nudged distribution (1,798 members) has the same variance
as the control (1,805 members) — ratio 0.955 [0.873, 1.051] — and the same shape
(Kolmogorov–Smirnov p = 0.37; Fig. 1b). The SSW moves the whole distribution of surface
outcomes; it does not divide it into responders and non-responders. The
observational record and CMIP6 say the same: the per-event surface response is
unimodal and consistent with a pure shift (observations p = 0.33, n = 43; CMIP6
p = 0.19, n = 1,888), although a two-population mixture with separation below
about 0.75 σ cannot be excluded at this sample size; and once the variance
already present before onset is removed, the forced between-event standard
deviation is at most 0.26 σ (CMIP6, upper 95% bound).

### The downward/non-downward contrast needs no SSW

We apply the surface conditions of the Karpechko et al. (2017) criterion to every
member and split each ensemble into downward-propagating (DW) and non-downward
(NDW) members over days 8–25. In the nudged ensembles, where the same SSW is imposed on every
member, the DW-minus-NDW surface contrast is −1.62 σ (25 ensembles; Fig. 2a). In
the control ensembles, where there is no SSW at all, it is −1.59 σ (36
ensembles). The contrast is the same with the event imposed and with the event
absent. It is, moreover, the value a threshold produces on noise: control members
are standardised by their own ensemble, and splitting a unit-variance Gaussian at
its mean gives a difference of group means of 2√(2/π) = 1.60 σ. The Southern
Hemisphere minor warming shows the same signature (−1.96 σ nudged, −1.58 σ
control).

The classification therefore does not measure an attribute of events. What does
separate the arms is the *rate*: on average 84% of nudged members are classified
DW against 45% of control members (Fig. 2b; 79% against 41% in the Southern Hemisphere). The
rate carries the causal shift; the contrast carries the threshold. In
observations the same decomposition holds: applying the published criterion to
ERA5 ourselves, 97% of the resulting DW–NDW contrast is reproduced by event-free
dates, and the fraction of the 39 events in ERA5 coverage that "propagate
downward" — 69.2% — is reproduced by event-free dates displaced by the measured
shift (74.4% [61.5, 87.2]). "About two thirds"
is what one shifted population produces.

### An SSW adds no predictability beyond the shift

If downward coupling were an event property, knowing the event should predict
how strongly it couples. We measure that directly, out of sample, in CMIP6 (1,888
events in 20 members). Ridge regressions predict each event's surface response
(days 8–52) from stratospheric winds at three sets of predictors: strictly before
onset, up to onset, and including the stratosphere after onset. Every variable is
standardised within each simulation member, so that the skill measured is within
a single climate, which is the question forecasters face. The same pipeline is
then applied to 1,000 sets of event-free winter dates with the same onset counts
and calendar days (Fig. 3a–c).

Real SSWs are not more predictable than ordinary winter days. Before onset the
out-of-sample R² is 0.067 at SSWs and 0.063 on event-free dates: the event adds
+0.004 [−0.022, +0.027] (p = 0.36). Including the post-onset stratosphere raises R²
to 0.40 — but to 0.38 on event-free dates as well, so 94.5% of that apparent
diagnostic skill is reproduced with no SSW present, and the event adds +0.022
[−0.015, +0.060] (p = 0.11). Resampling simulation members as well as dates
widens the intervals (pre-onset [−0.048, +0.051]; post-onset [−0.054, +0.092])
without changing the conclusion.

DW/NDW contrasts obtained with the published criterion imply much larger
event-specific shares. Converted to the fraction of between-event variance a
two-class split explains, the CMIP6 DW–NDW contrast and the observational
Karpechko AO contrast imply R² of 0.63 and 0.60 — six to seven times the largest upper bound we measure (Fig. 3d). Not
every published split overstates: the ERA5 NAM contrasts imply 0.09–0.10, near
the upper bound of the post-onset event-specific share.

### The field's two archetypes are two draws

February 2018 and January 2019 are the standard pair of opposites: the first a
strong downward-propagating event with a cold northern Eurasia, the second an
event whose surface response was weak or reversed, a difference attributed
chiefly to the strength of the SSW (Rao et al. 2020). SNAPSI nudges each model's
stratosphere to each event, so their forced consequences can be compared
directly. Averaged over both initialisations and the eight models, the two
events receive comparable forced shifts (−1.18 σ for February 2018, −1.04 σ for
January 2019) and the same probability of a DW outcome (0.82 and 0.83). At the
two initialisations with comparable lead to onset (18 and 20 days; Fig. 4a) the
shifts are −1.01 σ and −0.94 σ, and each event has the larger one in half of
the models.

The observed outcomes are ordinary draws from these forced distributions. Placed
among each model's nudged members at the primary initialisation, observed
January 2019 falls at a median percentile of 0.68 and outside the central 95% in
none of nine models (at the short-lead initialisation, which starts after onset,
it falls outside in five). The
observed difference between the two events lies within the central 95% of all
pairings of a 2018 member with a 2019 member in nine of nine models (Fig. 4b).
Finally, the label that separates them sits on the threshold (Fig. 4c): under
the Karpechko criterion at 1000 hPa, January 2019 is non-downward because its
mean NAM over days 8–52 is +0.019 σ rather than negative; at 850 hPa, or over
days 8–25, it is downward. The pair the field uses to show that SSWs differ in
their surface impact received the same forced odds, produced ordinary
outcomes, and is separated by a label that flips with the choice of pressure
level.

The same holds for the lower-stratospheric precursor. Loeffel et al. (2026) find
that week-2 geopotential height at 100 hPa predicts the weeks 3–7 surface
response across 18 events (r = 0.85) and conclude that SSWs differ in their
capacity to couple downward. Between members of a single nudged ensemble, where
no event-to-event difference exists, the same relation is weak (r = 0.11) and is
equally present with no SSW (r = 0.09; difference Δz = +0.02 [−0.04, +0.09],
36 matched ensembles in nine models; Extended Data Fig. 1). Across ensembles, the
correlation of ensemble means is 0.64 with the SSW imposed and 0.83 with no SSW,
driven by differences between models. A large correlation of
ensemble means does not require the events to differ.

### Discussion

After an SSW, the surface distribution is displaced by about one standard
deviation of the model's own spread, and it keeps its shape. That single number
is where the skill documented by Sigmond et al. (2013) lives: a forecast that
knows the stratosphere knows the shift. What the downward-propagation
classification adds on top of it — a claim about which events will couple and
which will not — is not supported by the experiment. With the event held
identical, members fall into both classes in proportions set by the shift, the
contrast between classes is the threshold's own, and the event-specific share of
predictability is indistinguishable from zero. An outcome in the "non-downward"
tail is an expected draw from a shifted distribution: under identical forcing,
about one member in six lands there.

The practical consequence is a change of forecast product. After an SSW the
appropriate statement is a shifted probability distribution — the probability of
a negative NAM, of a cold-air outbreak, of a given precipitation anomaly — rather
than a categorical prediction of whether this event will propagate. The same
applies wherever "downward-propagating SSW" is used to stratify impacts: the
classification rate against a matched null is an honest quantity; the contrast
between classes is not.

Three limits bound these conclusions. First, the evidence is model-based: with
42 usable observed events, the observational arm has 35% power at a true R² of
0.10 and
cannot settle the question on its own, although every observational test we can
make agrees with the models. Second, SNAPSI contains two Northern Hemisphere
events, so statements about how events differ rest on n = 2; the conclusion that
the classification contrast is manufactured does not depend on n, because it
holds in every ensemble in which a contrast can be formed. Third, the Southern Hemisphere case is a minor
warming, not an SSW. None of this implies that SSWs do not matter: the shift they
cause is large, causal, and the foundation of the skill that follows them.

*(Main text: ~1,860 words; limit 3,000.)*

---

## Methods

**Event catalogue.** Observed SSWs are the 43 major events of the NOAA CSL
compendium over 36 winters 1958–2024 (primary catalogue, frozen; file checksum
recorded in the repository). No inline date lists are used.

**Reanalysis.** ERA5 is read from the WeatherBench2 public copy (1.5°, 6-hourly).
Polar-cap sea-level pressure is the cos-latitude-weighted mean over 60–90° N
(60–90° S for the SH case), including the 60° row. NAM indices at 1000, 850 and
150 hPa follow the polar-cap definition (65–90° N geopotential height anomaly,
negated and standardised by day of year).

**SNAPSI.** Nudged and control ensembles for CCCma, CNR-ISAC, ECCC, ECMWF, KMA,
Meteo-France, NCAR, SNU and UKMO (38–52 members each) at initialisations
s20180125, s20180208 (February 2018; central date 12 February), s20181213,
s20190108 (January 2019; central date 2 January) and s20190829 (SH, central date
18 September 2019). NRL is excluded because its nudged-minus-control effect is
identically zero (duplicate submission). Forecast lead is measured from 00 UTC on
the initialisation date; the time origin of every ensemble was measured from its
files (UKMO and Meteo-France start at 06 UTC). Post-onset windows are applied only
to ensembles whose forecasts span them. s20190108 initialises six days after
onset and is reported separately where its short lead matters.

**Causal shift (Fig. 1a).** S = mean(nudged) − mean(control) of the day +8..+25
polar-cap sea-level pressure, divided by each model's control spread pooled over
its initialisations; reported in the NAM sign convention (−S, negative =
downward). Standard errors from the two ensemble variances.

**Distribution shape (Fig. 1b).** Members standardised by their own ensemble's
control mean and s.d.; each ensemble re-centred on its own mean before pooling,
because pooling without re-centring adds the between-ensemble spread of shifts to
the nudged variance only. Variance-ratio interval from 1,000 member bootstraps
within ensembles (4,000 resamples). A two-component mixture test was also run
but has no power at this sample size and is not used as evidence.

**Classification (Fig. 2).** Karpechko et al. (2017) conditions 1–2, over post-onset
days +8..+25 (the window every SNAPSI initialisation covers), on the member's NAM
proxy −(polar-cap sea-level pressure − control mean)/control s.d.:
window mean negative and more than half of 6-hourly values negative. The contrast
requires at least three members in each class; the DW rate is reported for all
ensembles.

**Predictability (Fig. 3).** CMIP6 zonal-mean fields for 20 members; SSWs detected
by the Charlton–Polvani reversal of the 10 hPa, 60° N zonal-mean wind
(November–March, final warmings excluded). Predictors: 10, 50 and 100 hPa zonal wind at
60° N and over the cap in windows before, at and after onset, plus seasonality;
38 features. Response: each model's annular-mode index (leading EOF of zonal-mean
sea-level pressure, 20–90° N) averaged over days +8..+52. Ridge regression, five-fold
cross-validation grouped by member, all variables standardised within member.
The null repeats the pipeline on 1,000 draws of event-free days (days influenced
by any real event removed), each with the same members and calendar days as the
real set; a cluster bootstrap (1,000 resamples of members) gives the second
interval. The 20 members come from 11 models (10 are CanESM5), so
standardisation and resampling are by member. Each draw and replicate has its own seeded stream, so the results are
identical at any degree of parallelism.

**Archetypes (Fig. 4).** Observed polar-cap sea-level pressure from ERA5 at exactly
the members' forecast times, placed in each model's nudged distribution using the
same control base and spread. The pair test compares the observed 2018-minus-2019
difference with all member pairings of the same model; a model's bias largely
cancels in the difference. Initial-time model-minus-ERA5 offsets are reported and
removed in a sensitivity.

**Loeffel comparison (Extended Data Fig. 1).** Week-2 (days 8–14) 100 hPa polar-cap
geopotential height against days 15–25 polar-cap sea-level pressure, correlated
across members within each ensemble; Fisher-z pooling; nudged-minus-control
difference on matched ensembles; permutation p-values. A pre-specified design
gate required the nudged/control spread ratio at 100 hPa to exceed 0.5 (median
0.96).

**Reproducibility.** Every result is produced by one script writing one JSON with
its seed, resample count and input counts; producer hashes are recorded and
checked. The environment is pinned. Defects found during the work, including
several that would have changed a result, are logged in the repository.

---

## References (verified from source; metadata from Crossref)

1. Sigmond, M., Scinocca, J. F., Kharin, V. V. & Shepherd, T. G. Enhanced seasonal
   forecast skill following stratospheric sudden warmings. *Nat. Geosci.* **6**,
   98–102 (2013). doi:10.1038/ngeo1698
2. Karpechko, A. Y., Hitchcock, P., Peters, D. H. W. & Schneidereit, A.
   Predictability of downward propagation of major sudden stratospheric warmings.
   *Q. J. R. Meteorol. Soc.* **143**, 1459–1470 (2017). doi:10.1002/qj.3017
3. Baldwin, M. P. et al. Sudden stratospheric warmings. *Rev. Geophys.* **59**,
   e2020RG000708 (2021). doi:10.1029/2020RG000708
4. Rao, J., Garfinkel, C. I. & White, I. P. Predicting the downward and surface
   influence of the February 2018 and January 2019 sudden stratospheric warming
   events in subseasonal to seasonal (S2S) models. *J. Geophys. Res. Atmos.*
   **125**, e2019JD031919 (2020). doi:10.1029/2019JD031919
5. Nebel, D. M., Garfinkel, C. I., Cohen, J., Domeisen, D. I. V., Rao, J. &
   Schwartz, C. The predictability of the downward versus non-downward
   propagation of sudden stratospheric warmings in S2S hindcasts. *Geophys. Res.
   Lett.* **51**, e2024GL110529 (2024). doi:10.1029/2024GL110529
6. Loeffel, S., Rupp, P., Kiefer, S., Pinto, J. G., Birner, T. & Garny, H.
   Quantifying the tropospheric response to individual sudden stratospheric
   warmings revealed by an ensemble simulation strategy. *Weather Clim. Dynam.*
   **7**, 895–913 (2026). doi:10.5194/wcd-7-895-2026
7. White, I., Garfinkel, C. I., Gerber, E. P., Jucker, M., Aquila, V. & Oman, L. D.
   The downward influence of sudden stratospheric warmings: association with
   tropospheric precursors. *J. Clim.* **32**, 85–108 (2019).
   doi:10.1175/JCLI-D-18-0053.1
8. Hitchcock, P. et al. Stratospheric Nudging And Predictable Surface Impacts
   (SNAPSI): a protocol for investigating the role of stratospheric polar vortex
   disturbances in subseasonal to seasonal forecasts. *Geosci. Model Dev.* **15**,
   5073–5092 (2022). doi:10.5194/gmd-15-5073-2022
9. Coughlin, K. & Gray, L. J. A continuum of sudden stratospheric warmings.
   *J. Atmos. Sci.* **66**, 531–540 (2009). doi:10.1175/2008JAS2792.1
10. Maury, P., Claud, C., Manzini, E., Hauchecorne, A. & Keckhut, P. Characteristics
    of stratospheric warming events during Northern winter. *J. Geophys. Res.
    Atmos.* **121**, 5368–5380 (2016). doi:10.1002/2015JD024226
11. Lu, R. & Rao, J. Sorting sudden stratospheric warmings with the downward
    tropospheric influence using ERA5 and CESM2-WACCM. *Atmos. Chem. Phys.* **26**,
    3723–3742 (2026). doi:10.5194/acp-26-3723-2026

*Before submission: convert author–year citations to Nature numbering; read
Nebel et al. (2024) in full — only its abstract has been read, and the text
cites only what the abstract states.*

---

## Figure legends

**Fig. 1 | An SSW shifts the surface distribution without splitting it.**
**a**, Causal surface NAM shift (nudged minus control, days +8 to +25) for nine
SNAPSI models and four Northern Hemisphere initialisations, in units of each
model's control spread; bars, ±1.96 s.e.; black ticks, model means; grey band,
mean ± s.d. over models excluding ECCC. **b**, Members of all 36 ensembles, each
re-centred on its own mean and pooled: control (grey), nudged (orange) and the
control translated by the mean shift (dashed), in units of each ensemble's control
spread.

**Fig. 2 | The downward/non-downward contrast needs no SSW.** **a**, DW-minus-NDW
surface contrast per ensemble with the SSW imposed (nudged) and absent (control);
circles, Northern Hemisphere; triangles, Southern Hemisphere minor warming;
black bars, means. 11 of 36 nudged ensembles have fewer than three NDW members
(their DW rate 0.99) and cannot form a contrast. Green dashed line, the contrast a
threshold at the mean produces on a unit Gaussian, −2√(2/π); grey dotted line,
the published ERA5 NAO contrast of Lu & Rao (2026), −0.708, for scale only (a
different index and yardstick). **b**, Fraction of members classified
DW, all ensembles.

**Fig. 3 | An SSW adds no out-of-sample predictability.** **a–c**, Within-model
cross-validated R² (variables standardised within simulation member) of the
surface response in CMIP6 at real SSWs (orange line)
and on 1,000 event-free sets of the same size (histogram), for three predictor
sets. **d**, Event-specific R² (real minus null mean; thick bars, 95% from the
null; thin bars, 95% from a member-cluster bootstrap) against the R² implied by
DW/NDW contrasts: four obtained by applying the published criterion to this
study's data (filled) and one published contrast in a different index (open, for
scale). The implied values measure variance explained by a diagnostic split and
are shown for scale, not as the same estimand.

**Fig. 4 | The field's two archetype events are two draws.** **a**, Nudged members
of nine models for February 2018 (s20180125) and January 2019 (s20181213) as the
surface NAM proxy over days +8 to +25 (negative, downward-propagating); bars,
interquartile range; filled diamonds, ERA5; open diamonds, ERA5 with the
model-minus-ERA5 initial offset removed; legend values are means over the eight
models excluding ECCC. **b**, Percentile of the observed
2018-minus-2019 difference among all member pairings of the same model; grey,
central 95%. **c**, Observed NAM means for both events under four variants of the
Karpechko et al. (2017) criterion; the only NDW label is January 2019 at 1000 hPa
over days 8–52 (+0.019 σ).

---

**Extended Data Fig. 1 | The lower-stratospheric precursor couples to the surface
as strongly with no SSW.** **a**, For each of 36 SNAPSI ensembles, the
within-member correlation of week-2 (days 8–14) 100 hPa polar-cap geopotential
height with the days 15–25 polar-cap surface response, with the SSW imposed
(nudged) against the same model and initialisation with no SSW (control); open
squares, the initialisation that starts after onset. **b**, Fisher-pooled
correlations with 95% intervals and the nudged-minus-control difference, for all
matched ensembles and excluding the short-lead initialisation.

## Data and code availability

SNAPSI data: CEDA archive (free registration). ERA5: WeatherBench2 public
cloud copy. CMIP6: ESGF. NOAA CSL SSW compendium: public. All code, result files
and run logs: the project repository.
