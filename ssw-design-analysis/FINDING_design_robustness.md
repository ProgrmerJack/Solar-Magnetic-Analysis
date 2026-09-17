# FINDING: the design choices do not reshape the estimate

**Date:** 2026-07-30
**Supersedes:** the working thesis "Study design reshapes estimates of sudden
stratospheric warming impacts"
**Status:** the planned headline is not supported. What replaces it is a
robustness result plus one honest null.

---

## 1. What was being claimed

The preregistered plan asserted that study-design choices — which event
catalogue, which estimator — materially change published SSW impact estimates,
and that this heterogeneity could be explained dynamically (Gate 8). Every
component of that claim has now been tested on the corrected 43-event
catalogue. None of it survives at the strength required.

## 2. Estimator choice: no effect

`07_physical_decomposition/design_sensitivity.py`, 43-event primary:

| quantity | value |
|---|---|
| mean conventional composite | −0.954 |
| mean corrected (within-winter) | −0.952 |
| **mean design sensitivity** | **−0.002** |
| R² of 5 dynamical predictors | 0.217, permutation **P = 0.28** |
| significant predictors | **0 of 5** |

Gate 8 fails, and it fails for an uninteresting reason: there is no estimator
gap left to explain. The conventional composite and the within-winter estimator
agree to 0.002 σ. The earlier figure of 0.008 σ was already negligible; the
catalogue correction made it smaller.

## 3. Catalogue choice: small, and two-thirds of it was record length

`05_corrected_estimators/common_period.py`. The spread across eight catalogues,
before and after fixing the analysis window to the span all six reanalyses
actually cover:

| AO bin | full record 1950–2026 | common window 1980–2019 |
|---|---|---|
| pre-onset −30..−16 | 0.426 | **0.141** |
| post-onset +15..+29 | 0.183 | **0.115** |

MERRA2 contributed 27 events against NCEP-NCAR's 39 not because it defines
events differently but because its record starts in 1980. Comparing their
profiles compared two periods as much as two event lists. Once the window is
common, catalogue choice moves the post-onset estimate by 0.115 σ against a
response of ≈1.0 σ. That is a robustness result, not a design problem.

## 4. The post-onset response is robust to everything tested

`05_corrected_estimators/clean_subset.py`, 8,000 replicates:

| subset | n | +15..+29 | 95% CI | p |
|---|---|---|---|---|
| all | 43 | −1.026 | [−1.612, −0.422] | 0.0013 |
| isolated only | 29 | −1.226 | [−1.845, −0.574] | 0.0005 |
| satellite era only | 29 | −0.934 | [−1.584, −0.279] | 0.0035 |
| isolated **and** satellite | 21 | −1.020 | [−1.765, −0.262] | 0.0088 |

Every interval excludes zero. Across the eight catalogues it runs −0.96 to
−1.14 on the full record and −0.94 to −1.05 on the common window. This is the
one solid quantity in the project.

## 5. The pre-onset anomaly is unresolved — NOT shown to be an artefact

This is where the temptation to overclaim was strongest, and it must be
resisted.

The pre-onset anomaly is −0.82 (p=0.015) on the full 43-event record and loses
significance under every restriction:

| subset | n | −30..−16 | 95% CI | p | CI contains −0.82? |
|---|---|---|---|---|---|
| all | 43 | −0.820 | [−1.413, −0.202] | 0.015 | — |
| isolated only | 29 | −0.548 | [−1.219, +0.143] | 0.117 | **yes** |
| satellite era only | 29 | −0.388 | [−1.041, +0.302] | 0.277 | **yes** |
| isolated **and** satellite | 21 | −0.275 | [−1.019, +0.501] | 0.496 | **yes** |

Two candidate mechanisms were tested formally, each as an interaction with its
own bootstrap interval rather than as a comparison of two separate fits:

- **era** (`era_split.py`): pre-1980 −1.735 vs 1980+ −0.468, difference −1.268,
  **p = 0.115**
- **event clustering** (`isolation_interaction.py`): isolated −0.564 vs
  clustered −2.940, difference +2.375, **p = 0.234**

Neither reaches significance. They are also not the same confound in disguise —
clustered events split 6 pre-1980 / 8 post-1980, Fisher **P = 0.49** — so these
are two independent tests that both fail to resolve the question.

Critically, **every restricted subset's interval contains the full-record
−0.82**. Nothing here contradicts the full-record estimate. The intervals simply
widen as n falls from 43 to 21. The correct statement is that the pre-onset
anomaly is *fragile and underpowered*, not that it is an artefact. 43 events
cannot settle it.

Supporting context, not evidence: the CPC daily AO index shows no
discontinuity at 1979 (Levene P = 0.26 on variance; lag-1 autocorrelation 0.943
vs 0.945), so an early-record data-quality explanation has no independent
support. And the large-ensemble literature already reports that the pre-SSW NAO
"has little bearing on its post-SSW state" (r = 0.19), so a large surface
precursor was never an established claim to overturn.

## 6. A reproducibility defect found in the process

`canonical_event_study.py`, `multi_index_event_study.py` and
`marginal_events.py` seeded the winter-block bootstrap with `hash(<str>)`.
Python randomises string hashing per process, so the seed changed on every run:
`hash("aao") % 10000` returned 3533, 7014, 6074 on three consecutive
invocations. No bootstrap p-value from those three scripts was reproducible.

This was caught because a re-run "changed" the AAO negative control from pass
to fail while its event set was provably identical — the four events the
catalogue fix added are all pre-1979 and the daily AAO record begins in 1979.
Fixed to `zlib.crc32(name.encode())` in all three files. At p≈0.05 with 1,200
replicates the Monte Carlo SE is ≈0.009, which is on its own enough to move a
result across the threshold; threshold-crossing claims need several thousand
replicates, and the numbers in §4 and §5 use 4,000–8,000.

`seed_stability.py` measured how much damage this did — same data, same
estimator, 40 fixed seeds:

| case | p range over 40 seeds | significant in | verdict |
|---|---|---|---|
| AO / primary, −30..−16 | 0.0050–0.0167 | 40/40 | stable |
| AO / **consensus_strict**, −30..−16 | 0.0317–**0.0617** | **32/40** | **seed-dependent** |
| AO / primary, +15..+29 | 0.0000–0.0067 | 40/40 | stable |
| NAO / primary, both bins | 0.0067–0.0333 | 40/40 | stable |

The headline quantities are stable. But the `consensus_strict` pre-onset result
straddles the threshold and changed significance on 8 of 40 seeds — it was
reported at p=0.048 in an earlier run, which was a coin flip presented as a
finding. Any conclusion that rested on that catalogue's pre-onset bin is void.

## 7. The negative control: passes, but exposes an inference problem

The 43-event multi-index run reported the AAO negative control failing at
+45..+60 (−0.659, surviving FDR at q=0.036). That was not seed noise:
`negcontrol_stability.py` found it significant on **60 of 60** fixed seeds at
6,000 replicates, p range 0.011–0.023. PNA was stably null (0/60).

`aao_null.py` then ran the clean pseudo-onset null on the AAO — real-event
influence deleted from the data first, then fake onsets drawn from surviving
days at the observed day-of-year set:

| quantity | value |
|---|---|
| null mean at +45..+60 | **−0.018** (estimator does not manufacture signal) |
| null SD | 0.368 |
| null 95% range | [−0.733, +0.691] |
| observed | −0.659 — **inside** the null range |
| randomisation two-sided p | **0.081** |
| winter-block bootstrap p | 0.0147 |
| bootstrap SE vs null SD | 0.268 vs 0.368 — **27% too small** |

**That 27% comparison is invalid, and so is the conclusion I first drew from
it.** The null SD came from 5,121 clean days while the bootstrap SE used all
8,638: √(8638/5121) = 1.299 rescales the null SD of 0.368 to 0.283, against a
bootstrap SE of 0.268. They agree. The apparent disagreement was data volume,
nothing else. The same confound voided `ao_null_perbin.py`, where a median SE
ratio of 0.80 across AO bins was exactly the 1/√(13864/8862) = 0.80 predicted by
sample size alone.

`boot_calibration.py` settles it properly — bootstrap SE and null SD both from
identical clean data, under a known truth:

| series | SE ratio (median) | coverage (median) | at the disputed bin |
|---|---|---|---|
| AO | 0.96 | 0.923 | — |
| AAO | 1.00 | 0.940 | +45..+60: ratio **1.06**, coverage **0.945** |

The winter-block bootstrap is calibrated for both series, and specifically at
the bin in question. Therefore:

So the raw −0.659 at p=0.0147 is a real number, not an estimator or seed
artefact. Two further tests found what it actually is.

**ENSO is not the mechanism.** ONI was acquired for this (NOAA PSL, 917 monthly
values 1950–2026) and entered as a daily covariate. It has a genuine direct
effect on the AAO (β = −0.383), but controlling for it moves the SSW–AAO
association not at all:

| series | without ONI | with ONI | attenuation |
|---|---|---|---|
| AAO +45..+60 | −0.659 (p=0.0147) | −0.661 (p=0.0137) | **0%** |
| AO +15..+29 | −1.027 (p=0.0013) | −1.023 (p=0.0013) | 0% |

**The failure was a mis-specified FDR family.** The q=0.036 that raised the
alarm came from pooling all 32 tests across four indices. Doing that lets the
AO and NAO results' very small p-values raise the BH threshold, which makes the
*negative control* easier to call significant — precisely backwards for a
control. Within the AAO's own 8-bin family:

| bin | p | q |
|---|---|---|
| +45..+60 | 0.0167 | **0.133** |
| +0..+14 | 0.0817 | 0.327 |
| all others | 0.21–0.97 | 0.55–0.97 |

**Gate 6 passes.** One bin of eight at raw p<0.05 is what chance predicts —
P(≥1 of 8 at α=0.05) = 0.34. The lesson is about family definition, not about
the estimator: a negative control must be judged against its own tests, never
pooled with the positive results it exists to check.

Note for the AO results: coverage of 0.923 against a nominal 0.95 is mildly
anti-conservative (true type-I error ≈0.077), though with 200 draws the Monte
Carlo SE on coverage is ±0.019 and 0.95 sits just inside the interval. AO
p-values near a threshold should be read with that slack in mind; the §4
post-onset results are far enough from 0.05 to be unaffected.

## 8. What this leaves

Supported:
- the SSW→AO surface response, −1.0 σ at +15..+29 d, robust to catalogue
  (0.115 σ), estimator (0.002 σ), era, and event isolation
- a calibrated, preregistered estimator with validated coverage (Gate 3, 0.947)
- a demonstration that the record period, not the event definition, drives most
  apparent disagreement between reanalysis-derived catalogues

Not supported:
- that design choices materially reshape SSW impact estimates — they do not
- that the estimator gap is dynamically predictable (Gate 8: P = 0.28)
- that the pre-onset anomaly is a design artefact — underpowered, unresolved

The honest summary is that this project set out to indict the field's design
choices and instead found the field's headline result to be robust to them.

## 9. The route to settling the open question

The pre-onset anomaly (§5) is the one live scientific question, and observations
cannot settle it at n=43. SNAPSI can: `nudged` minus `control` is the causal
stratospheric contribution by experimental design, with 50 members per centre
per case.

An earlier round of this project concluded SNAPSI was unreachable behind CEDA
authentication. **That was wrong** — it probed `/badc/snapsi`, which does not
exist. The real root, `/badc/snap/data/post-cmip6/SNAPSI`, is fully browsable
anonymously; only the download needs a free account.

`acquire_snapsi.py` was rewritten against the verified layout and has enumerated
the archive without credentials: `snapsi_manifest.csv`, **21,486 files**, 11
centres, 5 experiments, 6 cases, 50 members each, with pattern-constructed URLs
sample-verified (21,386/21,486; the 2 failing nodes flagged, not trusted).

Per-centre file size varies 34× for the same variable (psl: 2.0 MB at NRL,
67.2 MB at UKMO), so the download is tiered rather than all-or-nothing:

| scope (nudged+control, psl, 6hrPt) | files | size |
|---|---|---|
| CCCma pilot | 385 | 2.1 GB |
| 4 smallest-grid centres | 2,547 | 16.8 GB |
| all 11 centres | 6,121 | 147 GB |

Blocked only on a free CEDA account.
