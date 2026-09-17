> **⚠ PRE-RE-FREEZE DOCUMENT.** Written before 2026-07-30, when the event
> catalogue was re-frozen from **39 events / 33 winters** to **43 / 36** after the
> consensus rule was found to be era-dependent
> (see `../02_event_catalogues/FREEZE_RECORD.md`). Every count below is the old
> build. The JSONs in `results/` were all regenerated; this document was not.
> Qualitative conclusions are mostly unaffected, but **no number here should be
> quoted without checking it against `results/`.**

# Gate 6 — multi-index event study and the negative control

Script: `multi_index_event_study.py` · Output: `multi_index_event_study.json`
Frozen primary catalogue (39 events / 33 winters), 3 harmonics, 1,200 bootstrap
replicates. Family = **32 tests** (4 outcomes × 8 bins), Benjamini–Hochberg at 0.05 as
preregistered.

| outcome | role | raw p<0.05 | **survives BH-FDR** |
|---|---|---|---|
| AO | canonical response | 3/8 | **3/8** — +15..+29, +30..+44, +45..+60 |
| NAO | Atlantic sector | 4/8 | **3/8** — +15..+29, +30..+44, +45..+60 |
| PNA | Pacific sector | 3/8 | **0/8** |
| AAO | **negative control** | 1/8 | **0/8** |

## The negative control passes

A Northern Hemisphere mid-winter SSW has no mechanism to shift the Southern annular mode at
0–60 day lag, so the AAO is exposed to identical dates, identical seasonality and an
identical estimator with no signal to find.

The script's built-in criterion — "zero post-onset bins significant at raw p<0.05" — flagged
a **failure** on one bin (+45..+60, −0.66, raw p = 0.017). **That criterion was wrong.** It
ignores multiplicity: with 4 post-onset bins, P(≥1 significant by chance at α = 0.05) = 0.19,
so a single hit is unremarkable. Under the correction the preregistration actually
specifies, **no AAO bin survives (0/8)**.

The negative control therefore **passes**, and this is the strongest estimator validation in
the project so far — stronger than any pseudo-onset null, because it uses the *real* event
dates and asks whether a signal appears where none can exist.

## Two substantive findings

**1. The circulation response is annular/Atlantic, not hemisphere-wide.** PNA shows nothing
surviving FDR. Its three raw-significant bins are all *pre*-onset and positive, which is the
pattern chance produces, and they vanish under correction.

**2. No pre-onset bin survives FDR in any outcome, including AO.** The AO −30..−16 bin
(−0.57, raw p = 0.070) does not even reach nominal significance under the primary catalogue.
Taken with the catalogue-dependence already documented, the "significant precursor" result
should be regarded as an artefact of the union catalogue and of uncorrected testing. Earlier
statements in this project that the pre-onset anomaly was a firm, genuine precursor are not
supported.

## Gate 6 status: NOT met

The condition is "the result is not driven by AO/NAO alone". The honest reading of this
table is that, among circulation indices, **it is** driven by AO/NAO — PNA contributes
nothing. That is a legitimate scientific result (the response is annular), but it means
Gate 6 cannot be closed on circulation indices.

Closing it requires **non-circulation outcomes** — regional temperature, precipitation,
snow, wind — which come from the replication sample and are therefore blocked behind the
literature census. Gate 6 remains open, and its dependency is now explicit.
