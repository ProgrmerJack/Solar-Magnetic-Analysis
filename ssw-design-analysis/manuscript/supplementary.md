# Supplementary Information

*Stratospheric warmings shift winter cold risk rather than creating two kinds of event.* Supplementary Notes hold text moved from the main text to keep it within the article length; every number is from `ssw-design-analysis/CONSOLIDATED_RESULTS.md` or the result JSONs it indexes. Supplementary Tables are built by `10_tables/tableS*.py` from `results/current/`.

## Supplementary Note 1 | Event differences: the bound, its assumption and its falsifier

Events do differ, but by little. In CMIP6 the forced variance between events,
after removing the variance already present before onset, is 0.019σ²
[−0.033, 0.072]: the forced probability of a negative outcome ranges across
events from about 0.62 to 0.87 (0.47 to 0.93 at the upper bound), accounting for
2% of the variance of a single event's label (at most 8%). Any one label is
therefore almost entirely the shift plus chance. The bound assumes that an
event's forced response is uncorrelated with the internal variability it adds
to. It also gives a falsifier: at its upper limit, forced differences could
produce at most about half the class contrast of the criterion in CMIP6 and two
thirds in observations, so latent classes producing those contrasts would have
to exceed the measured bound. An audit finds the falsifier robust to the
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
