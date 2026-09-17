# What would make this Nature Geoscience — and what would not

**Written 2026-08-04, after reading the full corpus and running the SNAPSI
selection test (result L).**

This supersedes the venue section of `GO_NO_GO.md` **only for the new framing**.
That document's verdict — *"NO GO for Nature Geoscience"* — was correct **for the
paper it was assessing**: a field-wide *Analysis* of the literature, whose
conditions 1, 2 and 10 required a systematic census with a second coder and were
struck as unachievable by a solo author. That paper is still dead. This is a
different paper.

---

## 1. Why the old framing could never reach NGeo

| blocker | status |
|---|---|
| conditions 1, 2, 10 — literature census, second coder | **struck**, unachievable solo (`REMOVING_THE_HUMAN_DEPENDENCY.md`) |
| the thesis "study design reshapes estimates" | **dead** — estimator moves the AO estimate 0.002 σ, catalogue 0.115 σ |
| Gate 8, design sensitivity predictable from dynamics | **failed**, permutation P = 0.28, 0/5 predictors |
| framing | corrective/methodological — a WCD paper, correctly judged |

A paper whose claim is *"the field's method is biased"* is a methods paper. NGeo
publishes claims about the Earth system. That is the whole gap.

---

## 2. What changed: result L is a designed experiment, not a modelled null

Every earlier demonstration that the DW/NDW contrast is manufactured rested on a
null the author built — pseudo-events, or CMIP6 under an assumed structure. A
referee can always argue the null is wrong.

`snapsi_selection_test.py` removes that. In SNAPSI's `nudged` runs **every member's
stratosphere is nudged to the same observed evolution**, so members differ only in
tropospheric noise. The causal driver is identical *by construction*.

| arm | what is true of it | DW−NDW contrast | DW rate |
|---|---|---|---|
| **nudged** | stratospheric forcing **identical** across members | **−1.518 σ** (sd 0.304, n=10) | **0.78** |
| **control** | **no SSW forcing at all** | **−1.558 σ** (sd 0.103, n=12) | **0.44** |

**The contrast is the same whether or not an SSW happened.** It cannot be
reporting a physical difference: in the nudged arm there is none, and in the
control arm there is no event. Meanwhile the **rate** separates cleanly, 0.78 vs
0.44 — independently reproducing the observational 70% vs 38% (P < 0.0001).

This is the experiment the literature never ran, and it needs no census, no
second coder, and no new data acquisition.

---

## 3. The paper that could be NGeo

> ### Downward propagation is not a property of individual sudden stratospheric warmings

A **positive physical claim**, not a methods complaint. The chain, all in hand:

1. **L** — with the stratospheric driver held identical by experiment, the standard
   classification still manufactures its full contrast; identical contrast with no
   SSW at all.
2. **H** — the post-SSW surface response distribution is **one shifted population**,
   power-bounded: a two-component mixture with the separation the literature's own
   contrast implies (~1 σ) would be detected at 86–100% power in n = 1888. It is
   not there.
3. **the shift predicts the field's headline number** — displacing event-free dates
   by the measured −0.639 σ gives a 74.4% [61.5, 87.2] pass rate against 69.2%
   observed. *"About two thirds propagate downward"* is what one shifted population
   produces.
4. **J** — an SSW adds ~nothing to out-of-sample predictability within a single
   climate (+0.002 pre-onset, +0.012 post-onset). You cannot predict which event
   couples because there is nothing to predict.
5. **K** — and yet the causal effect is **large and real**: S = +0.97 σ from nudging.
   This is emphatically *not* "SSWs don't matter."
6. **the constructive replacement** — report the classification **rate** against a
   matched null. It is uncontaminated, strong, and costs nothing in rigour.

**Why this is Earth-system significance rather than methodology:** it dissolves a
15-year organizing concept. "Downward-propagating SSW" is used as an exposure
variable in air quality, energy demand, cold extremes and S2S forecasting. If it
is a threshold on one continuum rather than a class of event, every study
stratifying on it is measuring its own selection rule — and there is a better
statistic available.

---

## 4. The go/no-go, stated in advance

**The SNAPSI expansion is the decision point.** Currently 3 of 11 centres are
usable (NRL corrupt — see §6). The archive is enumerated: 21,486 files, 11 centres,
6 cases, ~50 members each; `nudged`+`control` psl for all 11 centres is ~147 GB and
blocked on nothing but disk and time.

| outcome | verdict |
|---|---|
| L holds across **≥8 centres AND the SH case** | **submit to Nature Geoscience** with the framing above |
| L holds in NH only, or in <8 centres | **npj Climate and Atmospheric Science** or **Communications Earth & Environment** |
| L fails to replicate | back to WCD as a methods note; the correction still stands |

Declaring this before running it is the point. The project's documented failure
mode is choosing a target and then finding a result to fit it
(`AUDIT_PROTOCOL.md`).

---

## 5. Work required, ranked

| # | task | why | cost |
|---|---|---|---|
| 1 | **download remaining 8 SNAPSI centres** | converts "3 models" into the community protocol; the go/no-go | ~147 GB, no permission needed |
| 2 | **run L on the SH 2020 case** | hemispheric generality; kills "NH-specific artefact" | trivial once (1) lands |
| 3 | **re-run `predictability_ceiling.py`** | its JSON still holds `Karpechko AO: 1.2852`, an impossible R² > 1 | minutes |
| 4 | **re-run `finalise_gate3.py`** | `gate3_final.json` has never existed; Gate 3 closure is unverifiable from artifacts | minutes |
| 5 | **intervals on the J-REVISED withdrawal** | three point estimates, one seed, no error bars — currently not publishable either way | hours |
| 6 | **precursor arm in `validate_R_diagnostic.py`** | the validation never simulated a precursor, which is why R's causal reading survived so long | hours |
| 7 | **fix 6 stale documents** | they cite the 39-event catalogue; `FINDING_canonical.md`'s headline is *reversed* by its own JSON | hours |
| 8 | **independent re-implementation of the D pipeline** | `AUDIT_PROTOCOL.md` already commits to this as the substitute for a co-author | days |

Items 3–7 are credibility, not science. A referee who finds an impossible R² > 1
in a supplementary file stops reading.

---

## 6. What must be disclosed, not buried

- **NRL's SNAPSI submission is corrupt.** Nudged-minus-control ensemble-mean effect
  is identically zero at all four initialisations; for `s20181213` the nudged
  submission is an exact duplicate of control. Worth reporting to the SNAPSI team
  independently of this paper.
- **Two bugs found in this project's own code on 2026-08-03/04**: `R_profile.py`
  pooled ensemble members by date (diluting the CanESM5 composite 3.5×), and the
  first version of the corruption guard compared members by label and passed NRL.
  Both are fixed and recorded. Seven substantive bugs total — disclose the count;
  it is evidence of the audit working, not against it.
- **The observational arm cannot carry the claim.** n = 43, power 0.35 at true
  R² = 0.10. Say so; do not let a referee find it.
- **Prior art is close.** White et al. (2019, *J. Climate* 32, 85 §3b) already
  wrote that DW−NDW positive-lag differences are *"entirely there by
  construction."* **Cite it prominently and state precisely what is new**: the
  designed-experiment demonstration (L), the power-bounded distributional test (H),
  the shift that predicts "two thirds", the β×ΔS correction law, and the
  out-of-sample predictability ceiling. Never lead with "unrecognised."

---

## 7. What to drop

- the literature census, and with it any prevalence claim
- "study design reshapes SSW impact estimates" — dead
- result C as a causal claim — R measures temporal concentration, corrected
- the snow-avalanche application, and `paper/ng_manuscript.tex` entirely
