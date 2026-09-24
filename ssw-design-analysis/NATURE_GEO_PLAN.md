# What would make this Nature Geoscience — and what would not

**Rewritten 2026-09-17, after the SNAPSI archive was completed and all four
experimental results re-run on it.** The 2026-08-04 version declared a go/no-go
before running the test. That was the point of declaring it, and this file now
records the outcome rather than moving the line.

---

## 1. The go/no-go, and the answer

The condition declared on 2026-08-04, before any of the data existed:

> | outcome | verdict |
> |---|---|
> | L holds across **≥8 centres AND the SH case** | **submit to Nature Geoscience** |
> | L holds in NH only, or in <8 centres | npj Clim Atmos Sci or Comms Earth Environ |
> | L fails to replicate | back to WCD as a methods note |

**It holds across 9 centres and it replicates in the Southern Hemisphere.**

| | contrast, nudged (SSW identical in every member) | contrast, control (no SSW at all) | DW rate |
|---|---|---|---|
| **NH**, 9 centres, 36 ensembles/arm | **−1.622 σ** (sd 0.303, n=25) | **−1.590 σ** (sd 0.089, n=36) | 0.842 vs 0.450 |
| **SH**, 8 centres, `s20190829` | **−1.958 σ** (n=7 of 8) | **−1.584 σ** (n=8 of 8) | 0.786 vs 0.413 |

The condition is met. **That does not mean the paper will clear the desk**, and
§5 states the odds honestly. It means the scientific precondition the project set
itself is satisfied and the remaining risk is framing and fit, not evidence.

---

## 2. What the archive changed, including what it nearly broke

The published SNAPSI results rested on **3 usable centres**. The archive is now
complete: 10 centres, 6 initialisations, 6,161 members, 1.17 M rows, of which
**9 centres are usable** (NRL excluded — its nudged-minus-control effect is
identically zero, duplicate-guard ratio 0.004).

| result | 3 centres (published) | 9 centres (now) | verdict |
|---|---|---|---|
| **K** causal effect S | +0.972 σ | **+1.107 σ** (sd 0.271) excl. ECCC; +1.363 incl. | holds, ECCC is an outlier |
| **L** contrast nudged vs control | −1.518 vs −1.558 | **−1.622 vs −1.590** | **replicates, 3× the ensembles** |
| **M** variance ratio | 0.939 [0.781, 1.130] | **0.955 [0.873, 1.051]** | holds, tighter — *after an estimator fix* |

**Three things the larger sample exposed. All three are in the record.**

1. **M's pooled estimator was contaminated.** Each ensemble is standardised by
   its own control mean and sd, so control lands at ~0 while nudged lands on that
   case's causal shift. Pooling without re-centring added the between-ensemble
   spread of shifts to the nudged variance and to nothing else. At 3 centres that
   spread was 0.295 σ and the bias was invisible; at 9 it is 1.071 σ and the
   uncorrected ratio reads **2.085** with KS p = 0.000 — a flat contradiction of
   M, produced entirely by the pooling. Re-centring restores 0.952 (0.955 with the 10 recovered members). **The naive
   scale-up would have overturned the result and been wrong.**
2. **L's DW rate was understated.** The contrast needs ≥3 members in both groups,
   which only an arm with a strong forced shift can fail. 11 of 36 nudged
   ensembles are dropped that way (their rate: 0.993); 0 of 36 control are. The
   unconditional rate is **0.842 vs 0.450**, not 0.776 vs 0.450.
3. **ECCC is a K outlier** at +3.407 σ against a Tukey fence of 1.669 on the other
   eight. Quote **+1.11 σ**, with +1.36 as the ECCC-inclusive sensitivity.

---

## 3. The framing that could clear the desk

The corrective framing — *"a classification in the SSW literature manufactures
its own contrast"* — is a methods paper about a subfield with **no Nature
Geoscience footprint**. Measured, not assumed: of **987 NGeo research articles
published 2023-01-01 to 2026-09**, titles matching sudden stratospheric warming,
polar vortex or stratosphere–troposphere coupling number **zero**. NGeo's
stratospheric output is Brewer–Dobson circulation, radiative cooling and aerosol
transport.

A full-text search of the journal's history returns **one** SSW surface-impact
paper, and it is the right one:

> **Sigmond, Scinocca, Kharin & Shepherd (2013), "Enhanced seasonal forecast
> skill following stratospheric sudden warmings", Nature Geoscience 6, 98–102**
> (doi 10.1038/ngeo1698, ~370 citations). Forecasts initialised at SSW onset show
> enhanced skill for circulation, surface temperature over northern Russia and
> eastern Canada, and North Atlantic precipitation.

**Nobody has decomposed where that skill lives.** This project can, and the answer
is a positive physical claim rather than a complaint:

> ### Forecast skill after a sudden stratospheric warming comes from a common shift, not from the individual event

1. **The shift is real, large and causal** — S = **+1.11 σ**, measured by nudging
   across 9 models (K). This is emphatically not "SSWs don't matter".
2. **The shift is all there is** — with the stratospheric driver held identical by
   experiment, the surface distribution translates without dispersing: variance
   ratio **0.955 [0.873, 1.051]**, KS p = 0.371, on 1,798 vs 1,805 members (M).
3. **So nothing event-specific is recoverable** — an SSW adds +0.004
   [−0.022, +0.027] to out-of-sample pre-onset R² over an ordinary winter day
   (size-matched null, 1,000 draws), and 95% of apparent
   post-onset diagnostic skill reproduces with no SSW present (J).
4. **Which explains the field's organising number** — displacing event-free dates
   by the measured shift gives a 74.4% [61.5, 87.2] "downward propagation" pass
   rate against 69.2% observed. *"About two thirds propagate downward"* is what
   one shifted population produces (H).
5. **The classification carries no extra information, proven by experiment** —
   the DW/NDW contrast is the same with the SSW held identical and with no SSW at
   all, in both hemispheres (L).
6. **Neither is the lower-stratospheric precursor** — between members of one
   event, week-2 100 hPa GPH correlates with the surface at r = +0.100; with no
   SSW at all, +0.101 (Δz −0.001 [−0.077, +0.076], 8 models). Löffel's r = 0.85
   across events does not require events to differ in coupling (N).
7. **The field's own showcase pair is two draws from one distribution** — Feb-2018
   ("propagating") and Jan-2019 ("not") receive the same DW odds (0.82 vs 0.83)
   and comparable forced shifts (1.18 vs 1.04 σ) in 8 models; their observed
   difference is inside the central 95% of member-pair differences in 9 of 9;
   and Jan-2019's NDW label holds in one of four criterion variants, by 0.02 σ (O).
8. **The constructive replacement** — report the classification **rate** against a
   matched null. It is uncontaminated and it separates cleanly: 0.842 vs 0.450.

**Why this is Earth-system significance rather than methodology:** it turns a
13-year-old result in this journal into a quantified forecast statement. After an
SSW the correct product is a **fixed shifted PDF**, not a conditional forecast —
which bears on subseasonal prediction of cold-air outbreaks, winter energy demand,
and the air-quality literature that already uses "downward-propagating SSW" as an
exposure variable (ACP 24, 1389, 2024; the practice continues in ACP 26, 3723, 2026).

It also **reconciles rather than attacks**: Sigmond is right, Karpechko's criterion
does select real events at a real rate, and Löffel's r = 0.85 is a real
correlation. What none of them establishes is that the coupling is a property of
the event.

---

## 4. Work remaining, ranked

| # | task | state |
|---|---|---|
| 1 | **Löffel head-to-head** (WCD 7, 895–913, 2026) | **done on 8 of 9 centres** (result N): within-ensemble r = +0.100 nudged vs +0.101 control, Δz −0.001 [−0.077, +0.076] on 28 matched ensembles. Outstanding: Meteo-France (`snap34` axis) and ECCC control, both need a valid CEDA token |
| 2 | **Intervals on the J-REVISED withdrawal** | **done**: event-specific +0.004 [−0.022, +0.027] pre-onset, +0.022 [−0.015, +0.060] post-onset, against 1,000 size-matched null draws; the old baseline was size-biased by +0.005–0.010 |
| 3 | **The operational/S2S bite** | **partly answered** (result O): the Feb-2018/Jan-2019 contrast that S2S studies (Rao 2020; Nebel 2024 frames NDW outcomes as forecast busts) treat as an event property is member-to-member noise under identical forcing. Still open: whether a published S2S *skill* number is conditioned on the observed outcome; read Nebel 2024 in full |
| 4 | Re-run H and I on the 9-centre archive | not started; both are CMIP6-based so unaffected, but the SNAPSI cross-checks should match |
| 5 | Figures and manuscript | **four main figures built** (`09_figures/fig1-4`); **draft** `manuscript/main.md` (183-word abstract, ~1,860-word main text); Extended Data Fig. 1 (N) waits on the Meteo-France zg; reference list to complete |

Items 6–8 of the old list are **done**: the impossible R² is cleared, Gate 3's
conflicting artifact is resolved, and all six stale documents are rewritten —
two of which had been asserting the opposite of their own JSONs.

---

## 5. The odds, stated honestly

**Nature Geoscience sends 15–20% of submissions to review and rejects about half
of those: roughly 7–10% overall.** No paper has a "high probability" of
acceptance there. With zero subfield precedent since 2013, the dominant risk is
the desk, not the referees.

| route | realistic odds |
|---|---|
| NGeo, corrective framing | ~0 — desk reject |
| NGeo, §3 framing, items 1–3 done | perhaps 15–25% of reaching review; ~10% overall |
| Nature Communications / npj Clim Atmos Sci / Comms Earth Environ | good — this is where the field's Nature-family work appears |
| WCD / JGR-Atmos / J. Climate | high — where SNAPSI work lands, and where Löffel published |

**The strategy that maximises expected outcome is not to pick NGeo. It is to build
the paper in NGeo's shape and submit there first with a pre-planned cascade.** A
desk rejection costs 1–3 weeks and the Nature Portfolio transfer service carries
the file onward. The §3 framing is also the strongest framing for the fallback
venues, so nothing is wasted.

**What would most move the needle is item 3.** If a published S2S or operational
skill claim is contaminated by the same cut, the paper stops being about a
classification and becomes about forecasts people actually issue.

---

## 6. What must be disclosed, not buried

- **The evidence is model-based.** The observational arm has 35% power at true
  R² = 0.10 and cannot carry the claim. Say so before a referee finds it.
- **Prior art is close.** White et al. (2019, J. Climate 32, 85 §3b) already wrote
  that DW−NDW positive-lag differences are *"entirely there by construction."*
  Cite it prominently. Never lead with "unrecognised."
- **SNAPSI has two NH events.** Anything about how events *differ* rests on n=2.
  Excluding ECCC, the between-model spread (sd 0.271) is about twice the
  between-event difference (0.142 σ) — models disagree about downward coupling
  more than these two events differ, but not by an order of magnitude.
- **The SH case is an austral MINOR warming**, not an SSW: the 10 hPa 60°S wind
  never reversed. It tests generality to a different hemisphere *and* a weaker
  class of event, on one initialisation.
- **NRL's SNAPSI submission is corrupt** and worth reporting to the SNAPSI team
  independently of this paper.
- **Bug count.** `FAILURES.md` logs **16 distinct recurrable defects** found and
  fixed in this project's own code, **6 of them during the 2026-09 expansion** —
  including one that would have overturned result M had the pooling not been
  checked, and one where a cache-filename change made two analyses silently use a
  third of the data. `CONSOLIDATED_RESULTS.md` §7 records the withdrawn claims
  separately. Disclose the count: it is evidence the audit works, not against it.

---

## 7. What to drop

- the literature census, and any prevalence claim
- "study design reshapes SSW impact estimates" — dead
- result C as a causal claim — R measures temporal concentration
- the snow-avalanche application, and `paper/ng_manuscript.tex` entirely
