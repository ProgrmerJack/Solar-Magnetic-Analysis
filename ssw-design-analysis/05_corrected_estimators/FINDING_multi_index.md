# Gate 6 — multi-index event study and the negative control

Script: `multi_index_event_study.py` · Output:
`results/current/2_event_study/multi_index_event_study.json`
Frozen primary catalogue (**43 events / 36 winters**), 3 harmonics, 1,200 winter-block
bootstrap replicates. Preregistered family = **32 tests** (4 outcomes × 8 bins),
Benjamini–Hochberg at 0.05.

> **Rewritten 2026-09-17 from the current artifact. Two of its conclusions were
> reversed by the 2026-07-30 catalogue re-freeze (39 → 43 events)** and had survived
> behind a warning banner. The document previously concluded that *no pre-onset bin
> survives FDR in any outcome* and that *PNA shows nothing surviving FDR*. Both are
> false on the current data, and the cause is the catalogue, not the test: PNA
> survives under the preregistered 32-test family, unchanged. The old tables are in
> git history.

## Effect by bin (σ), primary catalogue, with BH over the preregistered 32-test family

| outcome | role | raw p<0.05 | **survives BH-FDR** |
|---|---|---|---|
| AO | canonical response | 6/8 | **5/8** — −30..−16, +0..+14, +15..+29, +30..+44, +45..+60 |
| NAO | Atlantic sector | 6/8 | **4/8** — −30..−16, +0..+14, +30..+44, +45..+60 |
| PNA | Pacific sector | 3/8 | **3/8** — −30..−16, −15..−1, +15..+29 (all **positive**) |
| AAO | **negative control** | 1/8 | **0/8**, judged on its own 8-bin family, min q = **0.24** |

Key values: AO −30..−16 = **−0.818, p = 0.0050, q = 0.027**; NAO −30..−16 = −0.325,
p = 0.018, q = 0.049; PNA −30..−16 = **+0.354**, p = 0.008, q = 0.030.

**Family-definition sensitivity, stated because it changes the count.** The
preregistration specified 32 tests. The project later established that a negative
control must be judged on its own family and never pooled with the positives
(`../CONSOLIDATED_RESULTS.md` §7.9). That correction was about how to judge the
*control*; it does not by itself license shrinking the *positive* family. If the
positive family is taken as 24 tests (AO, NAO, PNA only), BH is less strict and the
counts rise to AO 6/8, NAO 6/8, PNA 3/8, adding −15..−1 for AO and NAO.
**The preregistered 32 is reported as primary here**, because moving to 24 post hoc
would be choosing the family that yields more significance — the project's own
documented failure mode. Nothing below depends on the choice: every claim holds
under both.

## The negative control passes

A Northern Hemisphere mid-winter SSW has no mechanism to shift the Southern annular
mode at 0–60 day lag, so the AAO is exposed to identical dates, identical
seasonality and an identical estimator with no signal to find. On its own 8 bins,
**0 survive and min q = 0.24**. One bin is raw-significant (+45..+60, −0.659,
p = 0.030) against 0.4 expected by chance from 8 tests.

The script's original built-in criterion — "zero post-onset bins significant at raw
p<0.05" — would flag that single bin as a failure. That criterion is wrong: with 4
post-onset bins, P(≥1 significant by chance at α = 0.05) = 0.19. The control passes,
and it remains the strongest estimator validation in the project, because it uses
the *real* event dates and asks whether a signal appears where none can exist.

## Two substantive findings, both the opposite of what this file used to say

**1. The circulation response is hemisphere-wide, and the Pacific signal is
opposite in sign.** PNA is not null: 3 of 8 bins survive FDR, all **positive**
(−30..−16 +0.354, −15..−1 +0.354, +15..+29 +0.283). The earlier reading — that
PNA's raw-significant bins were "the pattern chance produces" and "vanish under
correction" — was true of the 39-event build and is not true now. A positive PNA
alongside a negative AO/NAO is a sign contrast, not an absence, and it is the kind
of directional result the project has committed to reporting rather than smoothing
over (`../GO_NO_GO.md` condition 7).

**2. Pre-onset bins DO survive FDR, in all three positive outcomes.** AO −30..−16
survives at q = 0.027, NAO at q = 0.049, PNA at q = 0.030. The previous statement
that the AO −30..−16 bin "does not even reach nominal significance under the
primary catalogue" described the 39-event build; the current value is −0.818 at
p = 0.005. The withdrawal of the precursor as "an artefact of the union catalogue
and of uncorrected testing" is therefore **itself withdrawn**.

**This does not make the precursor causal, and the reason it was doubted still
stands on other grounds.** A precursor cannot be caused by the event that follows
it, and the model evidence remains contradictory — CanESM5 gives −0.387 while
MIROC6 gives **+0.380**, both p < 0.0001, so the question is model-structural
(`../CONSOLIDATED_RESULTS.md` §7.3). What is now established is narrower: the
observed pre-onset depression is not an artifact of catalogue choice and not an
artifact of multiple testing. It is a real feature of the observational record
whose interpretation is open.

## Gate 6 status: **MET on circulation indices**

The condition is "the result is not driven by AO/NAO alone". On the 39-event build
the answer was no — PNA contributed nothing and the response looked purely
annular/Atlantic. On the current catalogue PNA carries three FDR-surviving bins, so
the response is **not** confined to AO/NAO.

This closes Gate 6 on circulation indices without needing the non-circulation
outcomes that were blocked behind the struck literature census. `../GO_NO_GO.md`
condition 6 records "NOT MET ... PNA 0/8" and must be updated.
