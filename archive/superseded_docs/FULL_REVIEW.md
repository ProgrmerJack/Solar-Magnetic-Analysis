# Full review: SSW–avalanche project, verification rounds R104–R114

**Date:** 2026-07-28
**Scope:** everything done in this audit, what it found, what it broke, and where the
project actually stands.
**Status of the science:** no manuscript in this repository is currently submittable.

---

## 1. Executive summary

This project set out to establish that sudden stratospheric warmings (SSWs) modulate snow
avalanche hazard. Over this audit it was tested about as thoroughly as the available data
permit, from four independent observing systems, two continents, five observable classes,
and up to 47 events.

**The central claim did not survive.** Neither did the methodological claim that replaced
it. Three successive headlines were each killed by a stricter specification of the same
data:

| version | headline claim | killed by |
|---|---|---|
| 1. Solar–magnetic | solar activity modulates polar vortex → surface hazard | abandoned before this audit |
| 2. Opposing-direction avalanche hazard | 68% natural-avalanche suppression + weak-layer deepening ("loaded gun") | wrong event catalogue (R106); seasonal confound (R108); no replication (R107, R110, R111) |
| 3. Compositing artefact | between-winter composites inflate SSW surface impacts | its own corrected analysis (R113) |

What the audit did produce is real, but smaller than any of those: a set of demonstrated
design defects, a verified estimation toolkit, a corrected event catalogue, ~150 MB of
newly acquired data, and a fully catalogued repository.

**The most reliable finding of this work is the pattern itself**: every positive result in
this project depended on a design choice that a stricter specification dissolved.

---

## 2. Starting position

At the beginning of the audit the manuscript (`paper/main.tex`, 10,771 words, 17 figures)
claimed:

- 68% reduction in Swiss natural dry-slab avalanche counts during SSW windows
  (gmRR = 0.32; 14/16 events; P = 0.004)
- buried weak layers deepening in 15/16 events (P = 0.0005)
- danger bulletins elevated in 13/16 events
- a "loaded-gun" compound hazard state: natural triggers suppressed while structural
  instability intensifies
- an operational blind spot for count-based monitoring

Target: *Nature Geoscience*.

---

## 3. Defects found

### 3.1 Input defects

| defect | detail | consequence |
|---|---|---|
| **Wrong event catalogue** | The manuscript's 16-event list contains `2012-01-11`, which is **not a major SSW in any of the six reanalyses** in the NOAA CSL compendium, and omits `2000-03-20` (MAR 2000) and `2010-03-24` (MAR 2010), both major SSWs inside the Davos record era. Canonical count for that era is **17**, not 16. | every result computed on the wrong exposure variable |
| **Over-detecting detector** | A locally implemented Charlton–Polvani scheme found 40 events 1979–2024 against ~31 canonical, admitting final warmings and double counts. | R102, R104, R105 event sets invalid |
| **Event period mislabelled** | Stated as 1958–2026 in places; events actually end 2024-03-04. | reproducibility |
| **Catalogue definition ambiguous** | 47 events is the *union* across reanalyses, while Methods say "ERA5 where available" (ERA5-only = 44). | not a reproducible definition |

Correcting the catalogue *strengthened* the original headline (gmRR 0.305 → 0.241;
14/16 → 15/17; P = 0.0021 → 0.0012). That correction was honestly obtained and is the one
piece of good news in this section.

### 3.2 Design defects

| defect | demonstration |
|---|---|
| **Seasonal confounding** in "SSW window vs all other winter days" | The same design applied to **total snow depth** — a quantity with no candidate mechanism for responding within days to a stratospheric event — returns +4.4 cm at **P = 1.6×10⁻⁵³**. That is the seasonal deepening of the snowpack, recovered because SSWs cluster mid-winter while the control pool is weighted to early/late season. |
| **Winter selection** in between-winter compositing | SSW-hosting winters are dynamically distinct. The composite therefore mixes event response with a standing winter offset. |
| **Single-window specification** | Fitting one lag window at a time makes the control group "all other days in that winter" — which contains the central date and the entire post-onset response. A "placebo" estimated this way is not a placebo. |
| **Reference-category inversion** | Tirol accidents, ±15 d window: IRR = 1.335, P = 0.004 — an apparent *increase*. Lag decomposition shows the whole post-onset period is suppressed and ±15 d is merely its least-suppressed slice. Sign inverts against a contaminated reference. |
| **Pseudo-replication** | The weak-layer result treats 2.85×10⁵ station-days as independent. The reported +1.68 cm is **2.8%** of a ~60 cm mean depth. |
| **Anti-conservative standard errors** | Austrian incidents: model-based P = 0.0004 versus winter-block bootstrap **P = 0.11**. Roughly an order of magnitude. |
| **Residual seasonality** (still unresolved) | Pseudo-onset randomisation returns a null centred at **−0.258, not 0**. Three annual harmonics do not fully remove the seasonal cycle. |
| **Baseline dependence** (still unresolved) | Varying the baseline gap 60/75/90 d moves post-onset AO from −1.00 to −1.14, pre-onset from −0.35 to −0.54. Only profile shape is robust. |

### 3.3 Code defects found and fixed

| bug | impact | fix |
|---|---|---|
| `%y` parses `58` as **2058** | corrupted compendium dates | century resolved from event-name year; assertion added |
| LAWIS API timestamps carry time-of-day | Austrian arm silently reduced from 4,423 events to **3** | `.normalize()`; assert ≥50 events per region |
| `Categorical.codes` is int8; `code + j*n_strata` overflows | crash / stratum collision in bootstrap | explicit int64 cast |
| Duplicate background processes writing one log | interleaved, unreliable output | distinct logs; verified single-process reruns |
| Injection-recovery reported as a raw ratio | falsely suggested the estimator was 3× biased | computed as a difference: **97.6% recovery**, estimator unbiased |

The last one was my own error in interpreting a validation test, and it would have caused
a correct estimator to be discarded.

---

## 4. Chronological account

### R104 — accident case-crossover
CAIC (US, 1951–2025) and LAWIS Tirol (1992–2024) accident records; winter fixed effects
absorb the ~10× growth in backcountry recreation. Tirol lag profile: flat before onset,
then IRR 0.82 → 0.74 → 0.70 across days 0→44. **Caught the reference-category inversion.**
Model P = 0.0004 collapsed to bootstrap **P = 0.11**.

### R105 — multi-region random-effects meta-analysis
Pooled Davos, Tirol, US-continental, US-maritime: **IRR 1.03 [0.79, 1.35], P = 0.80**,
I² = 39.5%, with Alpine and US arms disagreeing in sign.

### R106 — canonical catalogue vs manuscript catalogue
Between-winter (manuscript design), canonical 17 events: **gmRR 0.241, 15/17, P = 0.0012**.
Within-winter: **IRR 0.622 [0.24, 1.14], P = 0.060**. First clear evidence that the
headline was design-dependent. Control pool: only **7 non-SSW winters**.

### R107 — pooled within-winter, conditional Poisson
Alpine primary domain, 30 region-winters: onset **IRR 1.01 [0.64, 1.50]**, placebo null.
Adding an independent Alpine region *removed* the effect rather than strengthening it.

### R108 — weak layers with seasonal control
285,157 SNOWPACK station-days.

| metric | uncontrolled | controlled | event-level |
|---|---|---|---|
| Pen_depth | +1.68 cm, P≈0 | −0.11, **P = 0.80** | 8/17 |
| min_ccl_pen | −0.211, P≈0 | −0.008, P = 0.75 | 9/17 |
| sk38_pwl | −0.781, P≈0 | −0.121, P = 0.29 | 8/17 |
| **HS_mod (control)** | **+4.39 cm, P = 1.6e-53** | −4.12, P = 0.31 | 7/17 |

The structural pillar was a seasonal artifact.

### R109 — Norwegian bulletin archive (found unused in the repo)
53,907 bulletins, 65 regions, 2013–2026, structured EAWS problems. With six events:
natural (spontaneous-release) problems **0.65, P = 0.016**; human-triggerable **1.12,
P = 0.001**; placebos null. The only positive replication found — and it is
forecaster judgement, not instrument.

### R110 — ALBINA Alpine replication (117,047 problems downloaded)
97 Alpine regions, 2018–2024, five events. Human channel **0.97 [0.80, 1.34] — does not
replicate Norway's 1.12**. Two of six placebos fail outright (`poor_stability`,
`high_frequency`, both P = 0.000), diagnosing six winters as insufficient regardless of
62,714 region-days.

### R111 — instrumental test, no humans in the chain
Snowpack temperature gradient `TG ≈ (0−T_air)/HS`, kinetic-growth threshold 10 K m⁻¹;
945 automatic SNOTEL stations, 47 winters, **31 events**. Kinetic growth 0.95 [0.76, 1.23];
melt 0.98; new snow 0.99; rain-on-snow 0.98. All placebos clean. **Overturned R102's
earlier SNOTEL "confirmation"**, which had used the over-detecting catalogue and a
between-winter comparison.

### R112 — positive control *(later withdrawn)*
Applied the within-winter estimator to the Baldwin–Dunkerton AO response. Reported:
within-winter recovers −0.614 (P = 0.0065) with a null placebo, while the conventional
composite fires at −0.656 (P = 0.0025) **before** onset. This became the basis of
`paper/ng_manuscript.tex`.

### R113 — event study; **withdrew R112**
External review identified that R112 fitted one window at a time, so the placebo's control
group contained the post-onset response. Refitted with mutually exclusive bins against an
explicit baseline (>75 d from any onset):

| bin | AO | P | NAO | P |
|---|---|---|---|---|
| −45..−31 | −0.455 | 0.13 | −0.113 | 0.38 |
| **−30..−16** | **−0.955** | **0.0010** | **−0.381** | **0.003** |
| **−15..−1** | **−0.603** | **0.018** | **−0.293** | **0.009** |
| +0..+14 | −0.864 | 0.005 | −0.283 | 0.002 |
| +15..+29 | −1.089 | <0.001 | −0.304 | 0.018 |
| +45..+60 | −0.919 | <0.001 | −0.356 | 0.002 |

**The within-winter estimator also shows a significant pre-onset anomaly.** The pre-onset
signal is therefore not a selection artifact — it survives winter fixed effects and is
substantially a genuine tropospheric precursor. R112's headline is false.

Validation that *did* pass: injection recovery **97.6%**; pseudo-onset randomisation
observed −0.999 vs null −0.258 ± 0.221, **P = 0.001**.

### R114 — avalanche results recomputed as event studies
All arms refitted in the corrected form. Three bins nominally significant (melt 0.88,
rain-on-snow 0.86, new snow 1.18, all +30..+44 d). Across the **36-test family**:

| correction | survivors |
|---|---|
| Benjamini–Hochberg FDR 0.05 | **0 of 36** (smallest q = 0.192) |
| Bonferroni | **0 of 36** |

~1.8 were expected by chance; three appeared, and they are mutually incoherent (loading
rises while melt and rain fall — opposite to the hypothesis).

**Correction to earlier reporting:** the Davos arm was previously called "null." Under the
event study its intervals span **0.16–13.84**. That is *uninformative*, not null — it
neither supports nor refutes the original claim.

---

## 5. What is established

| finding | confidence |
|---|---|
| The manuscript's 16-event catalogue is wrong | **high** — checked against all six reanalyses |
| Uncontrolled window-vs-rest designs manufacture significance from seasonality | **high** — the snow-depth control has no candidate mechanism |
| Model-based SEs are anti-conservative here | **high** — 0.0004 vs 0.11 |
| Single-window specifications are unsound | **high** — demonstrated on AO |
| The estimation machinery is correct | **high** — reproduces dummy GLM to 1e-5, OLS to 1e-7, injection recovery 97.6% |
| SSW → negative AO days 0–60 | **high** — textbook, and recovered |
| AO is *also* significantly negative before onset | **high** — implies real precursors |
| No avalanche/snowpack response survives multiplicity | **moderate–high** |

## 6. What is NOT established

- That between-winter composites are **upper bounds**. SLP flips sign; bias can amplify,
  attenuate or reverse.
- That the pre-onset anomaly is a **selection artefact**. It survives winter fixed effects.
- That published SSW impacts are **generally** inflated. No systematic survey was done.
- That the Davos avalanche effect is **absent**. Uninformative, not null.
- Anything from the Norwegian bulletin arm — forecaster judgement, six events, unresolved.

---

## 7. Data

**54 GB, 8,739 files.** Catalogued in `docs/DATA_INVENTORY.md`.

| group | size | relevance |
|---|---|---|
| POES particle flux + PSP/GOES/OMNI solar | **43.9 GB** | **none** — legacy from the abandoned solar hypothesis |
| Aura MLS O₃/Temperature | 5.1 GB | stratospheric verification |
| ERA5 polar stratosphere | 2.3 GB | **1979–2014 only** — does not cover 2018/19/21/23 SSWs |
| SNOTEL daily (945 stations) | 80 MB | the instrumental test |
| **all avalanche data** | **446 MB = 0.8%** | the entire evidence base |

**Structural fact:** of ~600,000 avalanche-related records, **exactly one region (Davos)
has observed natural avalanche activity counts.** Everything else is a danger rating
(forecaster prediction), an accident (human exposure), or model output.

### Acquired during this audit (~150 MB)

| dataset | content |
|---|---|
| LAWIS full | 4,859 incidents, 11 countries, 1992–2026 (was Tirol-only, 3,060) |
| ALBINA CAAML | 117,047 EAWS problems, 97 Alpine regions, 2018–2024 |
| NOAA CSL compendium | per-reanalysis central dates, six reanalyses |
| Catalogue extension | 2024-01-16, 2024-03-04 (Lee et al. 2025, *Weather*); March 2025 excluded as a final warming |

Checked and rejected: SLF CAAML API (serves only current bulletins); EAWS pan-European
archive (starts winter 2021/22 — too few events); Swiss `danger_descriptions` (dry/wet
only, not EAWS taxonomy).

---

## 8. Repository state

| status | n |
|---|---|
| CURRENT (corrected pipeline) | **9** |
| SUPERSEDED (old catalogue / over-detector / single-window) | **25** |
| UNREVIEWED (predate the audit) | **300** |
| result files | 188 |
| orphan results (no producing script) | 75 |

Structure: `docs/` (5 live reference documents), `archive/superseded_docs/` (4 withdrawn
narratives, retained), `archive/2026-07_root_cleanup/` (12 build artefacts), clean root.
Both manuscripts carry DO-NOT-SUBMIT banners.

Re-runnable audit tooling: `scripts/utilities/catalog_results.py`,
`scripts/utilities/inventory_datasets.py`.

---

## 9. Methodological lessons

1. **Run the placebo before believing the effect.** Every headline here would have been
   caught earlier by a pre-onset window.
2. **Validate the estimator on a known signal.** R112 looked decisive precisely because its
   positive control was mis-specified.
3. **A control variable with no candidate mechanism is the cheapest confound test there is.**
   Total snow depth at P = 1.6e-53 settled the weak-layer question instantly.
4. **Model SEs are not usable for clustered geophysical data.** Order-of-magnitude error.
5. **Region-days are not events.** 62,714 region-days across six winters gave failing
   placebos; the winter-block bootstrap correctly removed the illusory power.
6. **Correct multiplicity before claiming a bin.** Three of 36 tests at P<0.05 is chance.
7. **Check the exposure variable first.** The catalogue was wrong for the entire project.

---

## 10. Honest assessment and options

The avalanche hypothesis has been tested with the best available data — instrumental,
multi-continental, up to 47 events — and has not survived. The compositing hypothesis that
replaced it was withdrawn by its own corrected analysis.

**Options, in order of cost:**

1. **Short methods note.** "Single-window and between-winter SSW composites conflate
   post-onset response with precursors and winter selection." Defensible today on
   R113/R114. Target *Weather and Climate Dynamics* or *GRL*. This is the honest
   deliverable from the work as it stands.
2. **The systematic study a referee asked for.** Classify 50–100 published SSW-impact
   papers by design, reproduce 8–12, re-estimate under the corrected specification, report
   the full distribution of changes *including those that grow*, pre-register selection
   rules. Only this would support a field-level claim. It is a substantial project.
3. **Stop.** A well-tested negative is a legitimate outcome, and this one is unusually well
   documented.

**What I would not do:** submit either manuscript, or run further specification searches on
the avalanche data. Three headlines have now died the same way, and a fourth found by
continued searching would carry no more credibility than the first three.

**Before any of it:** the six residual problems in `docs/STATE_OF_EVIDENCE.md` — incomplete
seasonal adjustment, baseline dependence, ambiguous catalogue definition, mislabelled event
period, no multi-reanalysis sensitivity, unvalidated TG proxy — must be resolved. They
affect the corrected pipeline too, not only the superseded work.
