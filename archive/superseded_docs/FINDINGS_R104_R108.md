# Verification round R104–R108: corrected inputs and re-tested claims

Date: 2026-07-26. All numbers below are reproducible from the scripts named.

This round did three things: it **corrected the SSW catalog** against the authoritative
published source, it **acquired the full multi-region Alpine incident record**, and it
**re-tested each of the manuscript's load-bearing claims** under designs that control the
confounders the original analyses did not.

---

## 0. Input correction: the SSW catalog was wrong

Script: `scripts/download_extra/download_ssw_compendium_canonical.py`
Source: NOAA CSL SSW Compendium (Butler et al. 2017, ESSD 9, 63–76) — per-reanalysis
central dates for all six reanalyses.
Output: `data/processed/atmospheric/ssw_canonical.csv`

The manuscript's event list (`data/results/ssw_event_catalog.csv`, 16 events):

| problem | detail |
|---|---|
| **spurious event** | `2012-01-11` is **not** a major SSW in **any** of the six reanalyses |
| **missing events** | `2000-03-20` (MAR 2000) and `2010-03-24` (MAR 2010) are major SSWs inside the Davos record era |

Canonical Davos-era count is **17**, not 16. A locally written Charlton–Polvani detector
used in R102/R104/R105 over-detected badly (40 events 1979–2024 vs 28–30 authoritative),
admitting final warmings and double counts; results depending on it are flagged below.

---

## 1. Headline suppression — survives, and improves, but is design-dependent

Script: `scripts/analysis_extra/r106_canonical_catalog_headline.py`
Observable: Davos natural dry-slab avalanche counts.

| catalog | between-winter (manuscript design) | within-winter (winter FE) |
|---|---|---|
| manuscript 16 (with spurious event) | gmRR 0.305, 14/16, P = 0.0021 | IRR 0.684 [0.30, 1.32], P = 0.119 |
| **canonical 17** | **gmRR 0.241, 15/17, P = 0.0012** | **IRR 0.622 [0.24, 1.14], P = 0.060** |
| canonical, March events excluded (15) | gmRR 0.303, 12/15, P = 0.018 | IRR 0.641 [0.23, 1.28], P = 0.091 |

Correcting the catalog makes the published result **stronger**. Good news, honestly obtained.

**But the design matters more than the catalog.** The manuscript compares SSW windows with
day-of-year-matched days drawn from **other winters**, and the control pool is only
**7 non-SSW winters**. When each winter serves as its own control, the effect roughly
halves (70% → 38% reduction) and significance becomes marginal (P = 0.060).

---

## 2. PWL deepening / "loaded gun" — a seasonal artifact

Script: `scripts/analysis_extra/r108_pwl_seasonal_control.py`
Data: 285,157 SNOWPACK station-days, 130 stations, 1997–2020.

The original analysis (`scripts/analysis/r82_pwl_depth_validation.py`) compares SSW-window
station-days with **all** other winter station-days — with no day-of-year, winter, or
station control. Snowpack deepens monotonically through the season and SSWs cluster in
mid-winter, so season alone predicts the reported sign.

| metric | uncontrolled (original design) | season + winter + station controlled | event-level (own-winter, DOY-matched) |
|---|---|---|---|
| `Pen_depth` (PWL depth) | **+1.68 cm**, P ≈ 0 | −0.11 [−1.49, 1.19], **P = 0.80** | **8/17**, P = 1.0 |
| `min_ccl_pen` | −0.211, P ≈ 0 | −0.008 [−0.046, 0.033], P = 0.75 | 9/17, P = 1.0 |
| `sk38_pwl` | −0.781, P ≈ 0 | −0.121 [−0.355, 0.077], P = 0.29 | 8/17, P = 1.0 |
| **`HS_mod` (total snow depth — CONTROL)** | **+4.39 cm, P = 1.6e-53** | −4.12 [−10.8, 3.7], P = 0.31 | 7/17 |

The last row is the decisive diagnostic: the uncontrolled design "detects" a hugely
significant increase in **total snow depth** during SSW windows. That is simply the
snowpack deepening through winter. **The design manufactures significance from seasonality.**

Every event-level test is a coin flip (7–9 out of 17). The structural-instability pillar
of the manuscript does not survive.

---

## 3. Opposing-direction / human-trigger increase — contradicted *in realised accidents*

> **Superseded in part by R109 (below).** This section shows no rise in realised
> *accidents*. R109 finds that forecaster-*assessed* human triggerability **does** rise in
> Norway with a clean placebo. Both can be true; read the two together.

Scripts: `r104_accident_case_crossover.py`, `r107_pooled_within_winter.py`
Data acquired this round: full LAWIS record (`scripts/download_extra/download_lawis_full.py`)
— **4,859 incidents, 11 countries, 1992–2026**, of which Austria 4,423 (Tirol 3,510 plus
913 from Steiermark/Salzburg/Vorarlberg/Kärnten/Nieder-/Oberösterreich); 4,504 with full
detail records. Previously the repo held Tirol only (3,060).

No increase in human-involved incidents is detectable anywhere:

* Austria (3,771 winter incidents, 1992–2024): IRR ≈ 1.0–1.1, all windows
* US CAIC (787 winter accidents, 1979–2024): no increase; post-onset 1.21 [0.92, 1.57]
* Tirol lag profile is *suppressive* after onset, not elevated

**Artifact caught:** a single ±15 d window in Tirol gives IRR = 1.335, P = 0.004 — an
apparent increase. It is a **reference-category artifact**: ±15 d is the least-suppressed
slice of an entirely suppressed post-onset period, so against a reference contaminated with
strongly suppressed days it looks elevated. The lag decomposition dissolves it. This is a
warning about every single-window contrast in the manuscript.

---

## 4. The venue-deciding test: pooled multi-region, within-winter — NULL

Script: `scripts/analysis_extra/r107_pooled_within_winter.py`
Method: conditional Poisson (`scripts/analysis_extra/condpois.py`, verified against
dummy-variable GLM to 1e-5) with region × winter strata and region-specific DOY harmonics;
bootstrap resamples **whole winters across all regions together**, preserving the fact that
Swiss and Austrian records share the same Alpine weather.

**ALPINE (primary — the manuscript's mechanism is explicitly *Alpine blocking*), 30 region-winters:**

| window | IRR | 95% CI | P(suppression) |
|---|---|---|---|
| placebo (−45,−16) | 0.97 | 0.57–1.45 | 0.41 ✔ placebo correctly null |
| **onset (−15,+15) — manuscript's window** | **1.01** | **0.64–1.50** | **0.52** |
| early post (0,+14) | 1.02 | 0.63–1.34 | 0.54 |
| mid post (+15,+29) | 1.11 | 0.61–2.06 | 0.59 |
| late post (+30,+44) | 1.00 | 0.62–1.51 | 0.52 |

US (secondary, out-of-domain), 49 region-winters: onset 0.81 [0.55, 1.16]; post 1.21 [0.92, 1.57].
ALL regions, 79 region-winters: onset 0.99 [0.64, 1.38].

Adding an independent Alpine region did not strengthen the effect — **it removed it.**

### The structural reason, and it is important
**Only one region in the world's accessible data has natural-avalanche-count records: Davos.**
Every other record acquired here (Austria, US, Italy, Slovakia) is *accident/incident* based
and therefore driven by human exposure. So:

* the suppression hypothesis concerns **natural releases**, and is testable in **exactly one region**;
* the independent regions test a **different observable** and show nothing;
* the pooled null is therefore not a clean refutation of the Davos result — it is evidence
  that the claim cannot currently be replicated, because the data to replicate it do not exist
  outside Davos.

---

## 5. Where this leaves the paper

| pillar | status after this round |
|---|---|
| Novel question (first stratosphere→avalanche test) | ✔ intact; literature check found no prior work |
| Natural-avalanche suppression (Davos) | ~ real in direction; gmRR 0.24 between-winter, but IRR 0.62 / P = 0.060 within-winter; single region **— but see R109: replicated in Norway** |
| PWL deepening / structural "loaded gun" | ✘ **falsified** — seasonal artifact (and again in R109) |
| Opposing-direction compound hazard | ✘ contradicted in realised accidents **— but see R109: supported in assessed problem type** |
| Multi-region replication | ✘ null; and not currently possible for the right observable |
| Forecast skill | ✘ BSS = −0.95 (unchanged) |

**Nature Geoscience: not reachable.** Two of the three pillars are gone, and the surviving
one is single-region and marginal under the correct design.
**Nature Communications: not reachable** on this evidence, for the same reason.

**What is honestly publishable, and worth publishing:**
a rigorous single-region result plus a genuinely useful methodological warning — that the
standard "SSW-window vs all other winter days" comparison used throughout this literature
manufactures significance from seasonality (the `HS_mod` control proves it), and that
single-window contrasts invert sign against a contaminated reference. Natural target:
*The Cryosphere*, *NHESS*, or *Weather and Climate Dynamics*.

**What would change the verdict:** natural-avalanche activity counts (not accidents) from
additional regions — the full SLF Swiss record 1970–present, French EPA/CLPA, Italian
AINEVA. That is a data-access problem, not an analysis problem.

---

# R109 — a NEW PILLAR found in data already in the repository

Script: `scripts/analysis_extra/r109_norway_problem_mechanism.py`
Data: `data/processed/cryosphere/norway_avalanche.parquet` — **53,907 Norwegian
avalanche bulletins, 65 regions, 2013–2026** (45,135 winter region-days), each carrying
*structured, forecaster-assessed avalanche problems*: weak-layer cause, trigger
sensitivity, avalanche type. This dataset was already in the repository and unused.

It provides an operational observable for the manuscript's mechanism that is independent
of the Swiss data, the SNOWPACK model, ERA5 and the accident records.

Pre-specified mapping (from the manuscript's own mechanism):
* `natural_trigger` = "Spontaneous release" problem cited → predicted **DOWN**
* `human_trigger` = "Low additional load" / "Easy to trigger" → predicted **UP**
* `PWL_persistent` = faceted-snow or surface-hoar cause → predicted **UP**

## Result: the trigger dissociation REPLICATES, with a clean placebo

| outcome | placebo (−45,−16) | onset (±15) | **mid post (+15,+29)** | **post (+15,+44)** |
|---|---|---|---|---|
| `natural_trigger` (pred. ↓) | 1.17 [0.39, 3.10] ✔ null | 0.90 | **0.54 [0.20, 0.87], P = 0.002** | **0.60 [0.29, 1.19], P = 0.031** |
| `human_trigger` (pred. ↑) | 0.92 [0.83, 1.75] ✔ null | 0.82 | **1.14 [1.04, 1.28], P = 0.004** | **1.13 [1.03, 1.24], P < 0.001** |

Both limbs move in the **predicted opposite directions**, in the **same lag window**, and
that window (+15 to +44 d) matches the downward-coupling timescale found independently in
the Tirol lag profile (R104). The pre-onset placebo is null for both — the requirement
that failed nowhere else in this project.

**This is the first independent support for the manuscript's central dissociation claim.**

## But two limbs FAIL, and one diagnostic invalidates a third

| outcome | placebo | verdict |
|---|---|---|
| `PWL_persistent` | **0.58 [0.49, 0.85], P = 0.023 — placebo NOT null** | **uninterpretable**; a significant "effect" *before* onset means residual confounding. Consistent with R108, where the Swiss structural claim was a seasonal artifact. |
| `danger_ge3` | 0.64 [0.40, 1.06], P = 0.045 — marginal fail | uninterpretable |
| `dry_slab` | 1.01 ✔ | no effect (1.03 [0.92, 1.14]) |

So the **structural limb of the mechanism fails in Norway exactly as it failed in
Switzerland**, while the **trigger-dissociation limb succeeds**.

## Honest limitations
1. **n = 4 SSW events** (2018-02-12, 2019-01-01, 2021-01-05, 2023-02-16). The archive
   begins in 2013. Region-day power is large; event-level power is not. Four events cannot
   distinguish an SSW effect from four unusual winters, and the winter-block bootstrap
   (13 winters, 4 treated) only partly reflects that.
2. These are **forecaster assessments**, not observed avalanches. They record what experts
   judged the problem to be — informative and operationally meaningful, but not an
   independent physical count.
3. This **does not** contradict the accident results: Austrian and US *accidents* show no
   increase, while Norwegian *assessed triggerability* rises. Both can hold — forecasters
   may correctly flag rising triggerability while realised accidents do not rise (effective
   warning, or reduced exposure in cold-blocked weather). Earlier phrasing that the
   opposing-direction claim was simply "contradicted" was too strong: it is **contradicted
   in realised accidents and supported in assessed problem type.**

## Revised status of the paper's pillars

| pillar | status |
|---|---|
| Natural-release suppression | Davos gmRR 0.24 (between-winter) / IRR 0.62, P = 0.060 (within-winter) **+ Norway RR 0.54–0.60, clean placebo** → now **two independent systems** |
| Human triggerability rises | **supported in Norway** (RR 1.13–1.14, clean placebo); **not** in realised accidents |
| Structural / PWL "loaded gun" | **falsified twice** — Swiss seasonal artifact (R108); Norwegian placebo failure (R109) |
| Multi-region Alpine replication | null (R107) — but the observable does not exist outside Davos |

---

# R109 (updated, 6 events) and R110 — Alpine replication FAILS

## Catalog extended with properly-sourced post-compendium events
The NOAA CSL compendium table ends at FEB 2023. Two further major mid-winter SSWs are
documented in the peer-reviewed literature and are now in `ssw_canonical.csv`, flagged
`date_source = Lee2025_Weather`:
**2024-01-16** and **2024-03-04** (Lee et al. 2025, *Weather*, doi:10.1002/wea.7656).
The March 2025 event is **excluded** — documented as that winter's early *final* warming,
outside the Charlton–Polvani major mid-winter definition.
Norway coverage 4 → **6 events**; ALBINA **5 events**.

## R109 Norway at 6 events — the pillar strengthens and the placebo clears

| outcome | placebo (−45,−16) | mid post (+15,+29) | post (+15,+44) |
|---|---|---|---|
| `natural_trigger` (pred. ↓) | 1.40 [0.48, 2.46] ✔ null | **0.68 [0.30, 0.91], P = 0.005** | **0.65 [0.44, 0.98], P = 0.016** |
| `human_trigger` (pred. ↑) | 0.96 [0.84, 1.07] ✔ null | **1.12 [1.06, 1.23], P = 0.001** | — |
| `PWL_persistent` (pred. ↑) | 0.81 [0.49, 1.18] ✔ **now null** (failed at n=4) | 0.93 | late post **1.23 [1.07, 1.36], P = 0.001** |
| `danger_ge3` (pred. ↑) | 1.00 ✔ null | 0.83 | 0.93 [0.61, 1.36] — **null** |
| `dry_slab` | 1.02 ✔ | 1.04 | 1.05 [1.00, 1.13] |

Adding two events strengthened every limb and cleared the PWL placebo failure. In Norway
the dissociation is solid.

## R110 ALBINA (Alps) — does NOT replicate

Data: 117,047 EAWS problems, **97 regions** (Tyrol / South Tyrol / Trentino),
2018-12-04 .. 2024-04-30, 62,714 region-days, **5 SSWs**, 6 winters.
Script: `scripts/analysis_extra/r110_albina_alpine_problems.py`
Downloader: `scripts/download_extra/download_albina_caaml_problems.py`

| outcome | placebo | post (+15,+44) | verdict |
|---|---|---|---|
| `skier_triggerable` (human channel, pred. ↑) | 0.98 [0.93, 1.06] ✔ null | **0.97 [0.80, 1.34]** | **NO effect — contradicts Norway's 1.12** |
| `spontaneous_types` (natural channel, pred. ↓) | 1.30 [0.24, 4.21] ✔ null (wide) | 1.47 [0.87, 2.42] — **wrong sign** (onset 0.72 ↓, mid/post ↑) | inconsistent |
| `persistent_weak_layers` (pred. ↑) | 0.79 [0.32, 1.66] ✔ null | 1.36 [0.81, 1.73]; late post 1.52, P = 0.050 | weak, marginal |
| `poor_stability` | **0.98 [0.96, 0.98], P = 0.000 — PLACEBO FAILS** | — | uninterpretable |
| `high_frequency` | **0.90 [0.66, 0.90], P = 0.000 — PLACEBO FAILS** | — | uninterpretable |
| `danger_ge3` | 0.72 [0.43, 1.77], P = 0.058 marginal | **1.67 [1.13, 2.24], P = 0.002** | danger rises, but placebo marginal |

**The Norway trigger-dissociation does not replicate in the Alps** — and the Alps are the
manuscript's own mechanism domain ("planetary-wave *Alpine* blocking"). The human channel
is flat (0.97 vs Norway's 1.12) and the natural channel changes sign with window.

Two of six ALBINA placebos fail outright, which is itself the diagnosis: **six winters is
not enough winters to identify these effects**, whatever the region-day count. Region-days
(62,714) create an illusion of power that the winter-block bootstrap correctly removes.

## Net position after R109 + R110

| claim | status |
|---|---|
| Natural-release suppression | Davos (gmRR 0.24 between-winter / IRR 0.62, P=0.060 within-winter) **+ Norway (0.65–0.68, clean placebo)** — two systems, two climates |
| Human triggerability rises | **Norway yes** (1.12, clean placebo); **Alps no** (0.97); realised accidents no |
| Structural / PWL | Swiss SNOWPACK falsified (R108); Norway late-post 1.23 (P=0.001) but not at the main window; Alps marginal (1.52, P=0.050) — **still not established** |
| Danger ratings rise | Alps yes (1.67, P=0.002, marginal placebo); Norway no (0.93) — **contradictory** |

The genuinely new, defensible result is the **natural-trigger suppression replicating in an
independent country and observing system**. The opposing-direction ("loaded gun") claim
remains unsupported outside Norway.

---

# R111 — the instrumental, human-free test. NULL at 31 events.

Script: `scripts/analysis_extra/r111_kinetic_growth_instrumental.py`
Output: `data/results/r111_kinetic_growth_instrumental.json`

Every observable tested before this is mediated by people: accidents (exposure),
bulletins and danger ratings (judgement), observed counts (observer effort). This test
uses **instruments only**.

**Physics.** The manuscript's structural limb is kinetic-growth metamorphism, which is
driven by the snowpack temperature gradient, `TG = (T_base − T_surface)/HS`. Under snow the
base sits near 0 °C and the surface tracks air temperature, so `TG ≈ (0 − T_air)/HS` is
computable from any automatic station reporting air temperature and ultrasonic snow depth.
The classical kinetic-growth threshold is 10 K/m, strong growth / depth hoar 20 K/m
(Colbeck 1982; Akitaya 1974).

**Data.** SNOTEL, 945 automatic stations, 1980–2026, aggregated to 84,110 state-days across
10 states and **47 winters**, covering **31 canonical SSWs** — the largest event sample
anywhere in this project, roughly five times the bulletin-era tests. Rate model with offset
`log(reporting stations)`, state × winter strata, DOY harmonics, winter-block bootstrap.

## Continental US (the manuscript's own "continental-specificity" domain)

| outcome | predicted | placebo (−45,−16) | post (+15,+44) |
|---|---|---|---|
| `kinetic` TG > 10 K/m | ↑ | 0.95 [0.78, 1.15] ✔ | **0.95 [0.76, 1.23]** |
| `strong_kinetic` TG > 20 K/m | ↑ | 0.87 [0.64, 1.26] ✔ | **0.92 [0.62, 1.40]** |
| `melt` Tmax > 0 | ↓ | 1.05 [0.91, 1.18] ✔ | **0.98 [0.91, 1.06]** |
| `new_snow` | ↓ | 0.99 [0.95, 1.03] ✔ | **0.99 [0.93, 1.06]** |
| `ros` rain-on-snow | ↓ | 1.06 [0.90, 1.21] ✔ | **0.98 [0.90, 1.06]** |

Maritime: equally null (`kinetic` post 0.96 [0.81, 1.12]).

**Every placebo is clean. Every effect is null.** At 31 events, with instrumental data and
tight intervals, there is no detectable shift in either the structural regime or the natural
trigger channels.

## This overturns R102

R102 reported SNOTEL as *independent confirmation* of trigger suppression
(melt d = −0.47, P = 0.006; rain-on-snow d = −0.50, P = 0.004, n = 39). R111 uses the same
network and finds melt 0.98 [0.91, 1.06] and ROS 0.98 [0.90, 1.06].

The difference is entirely **design and catalog**: R102 used the over-detecting local
detector (40 events vs 31 canonical) and compared SSW windows with day-of-year-matched days
from *other winters*. R111 uses the canonical catalog and makes each winter its own control.
The earlier SNOTEL "win" was the same between-winter artifact as the Davos headline.

---

# THE UNIFYING RESULT

Across four independent observing systems, two continents and four observable classes, one
pattern holds without exception:

| design | result |
|---|---|
| SSW window vs **other winters** (uncontrolled / between-winter) | **significant** — Davos gmRR 0.24 (P = 0.0012); SNOWPACK PWL +1.68 cm (P ≈ 0); SNOTEL d ≈ −0.5 (P ≈ 0.005) |
| **each winter its own control**, with pre-onset placebo | **null** — Davos IRR 0.62 (P = 0.060); PWL −0.11 (P = 0.80); SNOTEL 0.95–0.99; pooled Alpine 1.01 |

The single exception is the Norwegian *forecaster-assessed* problem archive (R109), which is
human judgement rather than instrument.

**Interpretation.** The apparent SSW–avalanche association is consistent with a
between-winter confound: winters containing SSWs differ from winters that do not (ENSO/QBO
phase, decadal circulation state), and any design that compares across winters absorbs that
difference into the "SSW effect". Once the comparison is made *within* winters — where the
event is the only thing that changes — the effect disappears in every instrumental dataset.

This is a real scientific result, and a genuinely useful one for the field, because the
between-winter design is standard in the SSW-impacts literature. But it is a **negative**
result, and it is not what a Nature Geoscience paper claims.

---

# R112 — THE POSITIVE CONTROL. Design validated, and the standard method shown to be confounded.

Script: `scripts/analysis_extra/r112_design_positive_control.py`
Output: `data/results/r112_design_positive_control.json`
Data: CPC daily AO/NAO 1950–2026 (**all 47 canonical SSWs**), NCEP daily Z500/SLP/U850 1979–2024 (32 SSWs).

Every null in R106–R111 rests on the within-winter design. That design had never been
validated. Two explanations were possible and they have opposite consequences:
(A) the avalanche effect really is a between-winter confound, or (B) the within-winter design
is over-conservative — SSW impacts last ~60 days, a third of a snow season, so winter fixed
effects might absorb the signal being tested. If (B), every null is uninformative and the
manuscript's original results stand.

The test: apply the identical design to an SSW impact that is beyond dispute — the shift to
negative Arctic Oscillation persisting ~60 days (Baldwin & Dunkerton, *Science* 294, 581, 2001).

## Result 1 — the design is VALIDATED

| AO window | **within-winter** | between-winter |
|---|---|---|
| **placebo (−45,−16)** | **−0.057 [−0.446, +0.343], P = 0.75 ✔ null** | **−0.656, 33/47 neg, P = 0.0025 ✗** |
| post (0,+60) — Baldwin–Dunkerton | **−0.614 [−1.038, −0.179], P = 0.0065 ✔** | −0.829, 36/47, P < 0.0001 |
| post (+15,+44) | −0.476 [−0.929, −0.033], P = 0.034 ✔ | −0.907, 34/47, P < 0.0001 |

The within-winter design **recovers the canonical negative-AO response at the textbook
magnitude (≈ −0.6 σ) with a clean pre-onset placebo.** Explanation (B) is dead. The design
detects real SSW surface impacts when they exist — so the avalanche nulls in R106–R111 are
informative evidence of absence, not blindness.

Corroborated on U850 (within-winter mid-post −0.173 [−0.317, −0.046], P = 0.007).

## Result 2 — the STANDARD method fails its own placebo

This is the finding with reach beyond avalanches. The between-winter composite — the design
the manuscript uses, and the field's standard SSW compositing approach — reports a
**significant AO anomaly 45 to 16 days BEFORE the SSW central date** (−0.656, 33/47 events,
P = 0.0025). NAO likewise (−0.180, P = 0.041). SLP post-onset P = 0.0025 while within-winter
is null.

An event cannot cause an anomaly six weeks before it happens. The pre-onset signal is the
**whole-winter offset**: winters that contain SSWs are dynamically active, weak-vortex
winters with a lower mean AO. Between-winter compositing attributes that standing difference
to the event. Within-winter compositing does not — its placebo is clean (P = 0.75).

*Honest caveat:* part of a pre-onset tropospheric anomaly can be a genuine precursor
(blocking often precedes SSWs). The diagnostic is not the pre-onset value alone but the
**contrast between the two designs at the same window**, together with the clean
within-winter placebo.

## What this makes possible

The manuscript's headline is produced by exactly the design that fails its own placebo:

| observable | between-winter (standard) | within-winter (validated) |
|---|---|---|
| AO (undisputed physics) | −0.83 | −0.61 ✔ real, and recovered |
| Davos natural avalanches | gmRR 0.24, P = 0.0012 | IRR 0.62, P = 0.060 |
| SNOWPACK weak-layer depth | +1.68 cm, P ≈ 0 | −0.11 cm, P = 0.80 |
| SNOTEL trigger/kinetic physics (31 SSWs) | d ≈ −0.5, P ≈ 0.005 | 0.95–0.99, all null |
| Pooled Alpine avalanche activity | — | 1.01 |

**The defensible high-impact paper is no longer about avalanches.** It is:

> *Between-winter compositing inflates the surface impacts of sudden stratospheric warmings.*

with (i) a positive control proving the corrected design recovers the canonical AO response,
(ii) a placebo test showing the standard design reports significant anomalies before onset,
(iii) five circulation variables and 47 events, and (iv) a worked application in which a
large, publishable-looking hazard effect (68% avalanche suppression) dissolves entirely.

That has genuine broad significance: SSW compositing underpins a large literature and
operational subseasonal forecasting. The avalanche work becomes the demonstrating case
study rather than the claim.

---

# R113 — peer review was right, and the fix OVERTURNS the R112 claim

Script: `scripts/analysis_extra/r113_event_study.py` · Output: `data/results/r113_event_study.json`

Review identified that R112 fitted **one lag window at a time**, so the control group for
any window was "all other days in that winter" — meaning the placebo estimate was compared
against a control group containing the central date and the entire post-onset response. A
null placebo therefore did not establish what R112 claimed.

R113 implements the requested fix: mutually exclusive lead/lag bins fitted **simultaneously**
against an explicit omitted baseline of winter days more than 75 days from *any* onset, with
days assigned to the nearest central date so bins stay exclusive when a winter holds two
SSWs (2023/24).

## Result: the pre-onset anomaly is present WITHIN winters too

| bin | AO | P | NAO | P |
|---|---|---|---|---|
| −60..−46 | −0.269 | 0.44 | −0.189 | 0.19 |
| −45..−31 | −0.455 | 0.13 | −0.113 | 0.38 |
| **−30..−16** | **−0.955** | **0.0010** | **−0.381** | **0.003** |
| **−15..−1** | **−0.603** | **0.018** | **−0.293** | **0.009** |
| +0..+14 | −0.864 | 0.005 | −0.283 | 0.002 |
| +15..+29 | −1.089 | <0.001 | −0.304 | 0.018 |
| +30..+44 | −1.046 | 0.002 | −0.490 | <0.001 |
| +45..+60 | −0.919 | <0.001 | −0.356 | 0.002 |

**The R112 headline is withdrawn.** The within-winter estimator does *not* report nothing
before onset; it reports a significant negative AO from −30 days onward. The pre-onset
signal therefore is **not** purely a between-winter selection artefact — it survives winter
fixed effects, which means it is substantially a genuine tropospheric precursor, exactly as
the review argued and as the literature on SSW precursors would predict.

## Validation of the estimator itself (which does pass)

* **Effect injection:** adding a known −0.5 synthetic post-onset response shifts the estimate
  from −0.999 to −1.487, i.e. **−0.488 recovered against −0.500 injected (97.6%)**. The
  estimator is unbiased. (The ratio reported by the script is mis-specified — it must be
  computed as a difference against the un-injected fit, not as a raw ratio.)
* **Pseudo-onset randomisation (1000 draws):** observed post-onset mean −0.999 versus a null
  mean of −0.258 (sd 0.221), **P = 0.001**. The post-onset response is genuinely stronger
  than random mid-winter dates.
* **Caution:** that null is centred at −0.258, not zero, so three annual harmonics do **not**
  fully remove the seasonal structure — mid-winter days differ from the baseline (which is
  weighted to early/late winter) even after adjustment. Any bin estimate carries this
  residual.
* **Baseline-gap sensitivity (60/75/90 d):** pre-onset −0.35/−0.45/−0.54, post-onset
  −1.00/−1.09/−1.14 — magnitudes drift with the baseline definition, so absolute values are
  not robust; the profile shape is.

## Consequence for the manuscript

`paper/ng_manuscript.tex` must **not** be submitted. Its central claims are now contradicted
by its own corrected analysis:

| manuscript claim | status after R113 |
|---|---|
| "within-winter placebo is null" | **false** — significant at −30..−16 and −15..−1 |
| "pre-onset signal is the standing winter offset" | **not supported** — survives winter FE |
| "conventional composites inflate; read as upper bounds" | **unsupported**, and SLP already flipped sign |
| "compositing manufactures significance" | **overstated** |

What survives is weaker and narrower: between-winter and within-winter estimators give
materially different magnitudes; the difference is *not* cleanly attributable to selection;
and the SSW surface response has a real precursor component that event-window composites
cannot separate from the post-onset response without an explicit event-study design.

The avalanche nulls (R106–R111) were computed with the same single-window specification and
must be recomputed in the event-study form before being relied upon.

---

# R114 — avalanche results recomputed as an event study. Nulls survive.

Script: `scripts/analysis_extra/r114_avalanche_event_study.py` · Output: `data/results/r114_avalanche_event_study.json`

R106–R111 all used the single-window specification that R113 showed to be unsound. All
avalanche arms were therefore recomputed in event-study form: mutually exclusive lead/lag
bins fitted simultaneously against an explicit baseline of winter days >75 d from any onset.

| arm | winters | pre-onset bins | post-onset bins |
|---|---|---|---|
| CH-Davos natural counts | 21 | 1.13, 0.73, 0.35 | 0.46, 0.63, 0.54 — **all CIs enormous** (e.g. 0.16–13.84) |
| Pooled Alpine | 33 | 1.22, 0.86, 0.92 | 0.86, 0.86, 0.80 — all null |
| SNOTEL kinetic growth | 33 | 1.00, 0.97, 0.97 | 1.00, 0.93, 1.15 — all null |
| SNOTEL melt | 47 | 1.04, 0.95, 0.97 ✔ | 0.94, 0.95, **0.88 (P=0.016)** |
| SNOTEL new snow | 47 | 1.07, 1.08, 1.10 | **1.13 (P=0.038)**, 1.08, **1.18 (P=0.016)** |
| SNOTEL rain-on-snow | 47 | 1.04, 0.96, 0.96 ✔ | 0.92, 0.95, **0.86 (P=0.012)** |

## Multiplicity: nothing survives

This family contains **36 tests** (6 arms × 6 bins). Correcting across it:

| | survivors |
|---|---|
| Benjamini–Hochberg FDR 0.05 | **0 of 36** (smallest q = 0.192) |
| Bonferroni | **0 of 36** |

The three nominally significant late-post bins (melt 0.88, rain-on-snow 0.86, new snow 1.18,
all at +30..+44 d) are within chance expectation for 36 tests (~1.8 expected at P<0.05; three
observed). They are also not mutually coherent: new-snow loading rises while melt and
rain-on-snow fall, and the manuscript predicted loading to fall.

**Conclusion: the avalanche nulls survive the corrected design.** Reassuringly, the melt and
rain-on-snow arms now have clean pre-onset bins, so their late-post point estimates are not
obviously contaminated — but they do not survive multiplicity and cannot be claimed.

The Davos arm deserves a specific caveat: under the event-study specification its intervals
span an order of magnitude (0.16–13.84). With 21 winters, one region and a baseline of only
1,391 day-observations, that arm is **uninformative in either direction** — it does not
support the manuscript's original claim, but neither does it refute it. Earlier statements
that the Davos effect was "null" under the within-winter design overstated what the data can
say; the honest statement is that it is unresolved.
