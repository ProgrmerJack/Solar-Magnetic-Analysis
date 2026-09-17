# Response to Referee & Revision Plan

**Manuscript:** "Stratospheric sudden warming episodes are associated with opposing-direction avalanche hazard signals through planetary-wave Alpine blocking"
**Prepared:** 2026-06 · Addresses the senior-referee report point-by-point.

> **UPDATE (this revision):** M1, M2, and M3 have now been **implemented and run** on the
> local data (scripts `r99`/`r100`/`r101` in `scripts/analysis_extra/`, outputs in `data/results/`,
> full numbers in `M1_M2_M3_RESULTS.md`). Headlines: **M2** — no statistically resolved SSW
> suppression beyond blocking (IRR_S = 0.88 [0.28, 1.43]); **M3** — 3 post-2019 SSWs verified, OOS
> direction matches in-sample (accidents 2/3 ↑, gmRR 1.49) but null at n = 3, and the count-channel
> OOS is blocked on post-2019 SLF data; **M1** — confirmation needs ≈57–70 events, and a manuscript
> power error ("<25% for |ρ|=0.5", actually ≈51%) was corrected. All results are integrated into
> `main.tex`/SI, which recompile cleanly.

Each item is tagged:
- ✅ **DONE (text)** — edit already applied to `paper/main.tex` or `paper/supplementary_information.tex` in this revision.
- ✍️ **REFRAME** — wording/positioning change applied or specified.
- 🔬 **NEEDS ANALYSIS** — requires running code on the restricted data; cannot be done by editing. A concrete, runnable protocol is given so the author can execute it.
- 🧭 **DECISION** — a strategic choice only the author can make; options laid out.

---

## A. Major comments

### M1 — Headline novelty (compound sub-type) is not statistically established
**Status: 🧭 DECISION + 🔬 NEEDS ANALYSIS + ✍️ REFRAME**

This is the load-bearing issue. The triple-positive joint test (P = 0.166) and cross-arm correlation (ρ = −0.33, P = 0.21, <25% power) cannot, at n = 16, establish an *irreducible* opposing-direction sub-type. Three honest routes:

- **Route 1 (recommended for NG attempt): demote the typology from "finding" to "framework/hypothesis."** Make the paper's *result* the robust descriptive dissociation (counts down, instability/danger up) and present the compound-event typology purely as an interpretive proposal. The manuscript already says "candidate … requiring larger-sample confirmation" — make that consistent in **title, abstract first sentence, and the one-line contribution claim**. Concretely: change the abstract's opening rhetorical framing ("Whether a single meteorological forcing can…is unknown. Here we show…") so that "Here we show" attaches only to the dissociation, not to the typology.
- **Route 2: power the irreducibility test by pooling events across independent regions/forcings** (see M3 protocol) to raise effective n. This is the only route that converts the candidate into a result.
- **Route 3 (venue): move to a specialist journal** where an exploratory mechanism paper at n = 16 is in-scope (see §D).

**Protocol to strengthen the joint claim without overclaiming** (`scripts/review_rounds/`):
1. Pre-register the three arms (suppression, bulletin elevation, PWL deepening) and the cross-arm statistic *before* re-running.
2. Report the joint test as a **power analysis**, not a p-value: "to detect the observed enrichment at 80% power requires N ≈ ___ events" — turning the negative result into a quantified sample-size requirement (constructive, honest).
3. If Route 2 data are assembled, re-estimate cross-arm ρ on the pooled set with region as a blocking factor.

### M2 — SSW-specificity beyond ordinary blocking is weak
**Status: 🔬 NEEDS ANALYSIS + ✍️ REFRAME (partly already in text)**

The text already reports IRR = 0.89 (P = 0.06) after regime conditioning and that non-SSW blocking suppresses *more* (RR 0.41 vs 0.95). The referee's core question — *what is the SSW-specific increment, with uncertainty* — is not yet answered as a single headline number.

**Protocol (run on existing Swiss daily series + ERA5 Z500 blocking index):**
1. Fit a single negative-binomial GLM on daily natural dry-slab counts with terms: `blocking_index`, `SSW_window`, `blocking_index × SSW_window`, winter FE, day-of-winter spline. The **interaction coefficient** is the SSW-specific increment beyond blocking. Report IRR and winter-block-bootstrap CI.
2. If the interaction CI spans 1 (likely), state explicitly in the Discussion: *"We find no statistically resolved avalanche suppression attributable to SSWs beyond their association with Alpine blocking (interaction IRR = ___, 95% CI [__, __]); the SSW's role is as a predictable 2–4-week precursor marker of blocking, not an independent suppressor."* This is a defensible, novel-enough framing and removes the overclaim.
3. Reuse `scripts/review_rounds/r98_*` (continuous-vortex GLM) as the code scaffold.

### M3 — Exploratory origin / multiple comparisons / mechanism switched between drafts
**Status: 🔬 NEEDS ANALYSIS + ✍️ REFRAME**

The proposed surface mechanism changed from **NAO** (earlier draft, `notebooks/nature_geoscience_paper.md`) to **Z500 blocking** (current), with NAO now a "rejected pathway." Referees who find the draft history will read this as mechanism-fitting. Two actions:

1. ✍️ **Disclose the evolution honestly** in Methods "Study origin": add one sentence — *"An earlier version of this analysis examined annular-mode (NAO) mediation; on fuller analysis the event-level annular-mode correlation was null (R² < 0.01) and the regional Z500 blocking pathway carried the signal, which is the mechanism reported here."* Turning a hidden pivot into a documented model-selection step defuses the HARKing concern.
2. 🔬 **Out-of-sample confirmation (the single highest-value new analysis).** Freeze the entire pipeline and apply it, unchanged, to **SSW events after the 2018/19 training cutoff** (2019/20–2024/25; the catalogue and Swiss data exist). Pre-commit the prediction (RR < 1) and the analysis script hash. Even 4–6 genuinely out-of-sample events providing directional confirmation would do more for credibility than any in-sample robustness check. Script scaffold: `scripts/analysis/35_prospective_2021_test.py` already exists — generalise it to the full post-2019 window.
3. ✍️ Keep the worst-case Bonferroni bound but **add the garden-of-forking-paths caveat explicitly**: the bound covers forcing choice only, not window/lag/outcome/regime-definition choices.

### M4 — Pseudo-replication inflates the most striking p-values
**Status: ✅ DONE (text)**

- SI "Aggregate surface weather changes" now carries an explicit caveat that P < 10⁻⁴⁷ etc. are pooled station-day descriptive statistics over >10⁴ autocorrelated days and that the valid inferential unit is the event (n = 16).
- **Remaining author action (light):** sweep the SI once more and ensure every sub-10⁻¹⁰ p-value derived from station-days carries the word "pooled/descriptive." A grep for `10^{-` in the SI will list them.

### M5 — SNOWPACK/ERA5 circularity limits the mechanistic chain
**Status: ✅ Already well-handled; ✍️ minor tightening**

The Methods circularity caveat is exemplary. One tightening: in the Results "Snowpack mechanism" subsection, add a half-sentence reminding the reader that the independent legs (Rutschblock, bulletin, French S2M) support **suppression/loading**, but do **not** independently establish the *opposing-direction structure* — that structure rests on combining the (shared-forcing) SNOWPACK instability arm with the independent suppression arm. This pre-empts the referee's exact objection.

### M6 — Operational/forecasting framing overstated
**Status: ✅ DONE (text)**

- `main.tex` Discussion: "most actionable forecast pathway" → now "contingent on future skill validation … no calibrated forecast skill is demonstrated here (LOO BSS = −0.95) … at most a 2–4-week categorical prior shift."
- SI plain-language summary: forecasting sentence now states no calibrated skill is demonstrated.
- Abstract already contained the BSS = −0.95 disclosure; no change needed there.

### M7 — Accident/fatality channel adds little
**Status: ✍️ REFRAME + 🧭 DECISION**

Already demoted to "supportive consequence-context" in text. Recommended further action (author decision): **remove the accident analysis from the abstract-adjacent narrative and Fig. count_rating_dissociation panel (c)**, relocating it to a single SI subsection titled "Consequence context (exploratory)." Rationale: event-level direction is null (15/29, P = 0.50) and the secular-trend confound invites suspicion of the lone Jan-2021 RR = 2.52 highlight. Cutting it *strengthens* the paper by removing the weakest limb. Left as a decision because it touches a figure you may want to keep.

### M8 — Authorship, data access, submission completeness
**Status: ✅ DONE (cover letter + title) + 🧭 DECISION**

- **Cover letter written** (`paper/COVER_LETTER.md`) — was entirely missing; NG requires it. Covers novelty, significance, exploratory-origin transparency, openness, suggested reviewers, and a direct question to the editors on venue fit.
- **Title consistency:** `main.tex` header comment now matches `\title{}`. (SI retains its own descriptive subtitle, which is acceptable.)
- **Byline/account:** reconcile the submitting email/account with the manuscript byline (`Jack00040008@outlook.com` / ORCID 0009-0003-5482-5526) before upload — 🧭 author action.
- **Single-author + restricted data:** consider whether an SLF or atmospheric-science collaborator could co-author or formally vouch for data QC; this materially helps editor confidence (decision).

---

## B. Minor comments

| # | Item | Status |
|---|------|--------|
| 1 | Human-trigger numbers (1.45 vs 1.81 vs 1.70) | ✅ SI extended abstract now states "increase in absolute terms (RR = 1.45; human-to-natural *differential* 1.81×)". Verify the same disambiguation in any other SI mention via grep `1.81`. |
| 2 | BF "approaching decisive" promotional | ✅ `main.tex` limitations now: "'strong' evidence; the wide range reflects prior sensitivity at n = 16, and we rely on the conservative lower bound." |
| 3 | MLS validates catalogue, not mechanism | ✅ Already explicitly stated in Results and Discussion; no overreach found. No change needed. |
| 4 | Norway shown as clean replication in SI summary | ✅ SI extended abstract + plain-language summary now flag n_eff ≈ 4 and "inconclusive / weaker." |
| 5 | Wet-avalanche increase fragile (3 events = 78%) | ✅ Already labelled hypothesis-generating in SI "Dry slab specificity"; ensure any main-text mention carries "speculative/spring-rebound." (Author: grep `1.62` in main.tex.) |
| 6 | Fig. 1 DAG: promote common-cause arrow | 🧭 Figure edit — the timing evidence (52% pre-onset) supports the dashed common-cause arrow; consider rendering it solid/secondary-supported rather than "hypothesised, not tested." Author decision (TikZ in `main.tex` ~line 367). |
| 7 | Spec-curve = sign stability only | ✅ `main.tex` now: "a demonstration of *sign* stability across analyst choices, not of effect-size or inferential independence, because the variants share the same 16 events." |
| 8 | Climate-change Discussion is speculative scope-padding | ✍️ Recommend trimming the two "competing processes" paragraphs to ~4 sentences; flagged, left to author to preserve preferred content. |
| 9 | Compound-typology literature thin | 🔬 Add 2–3 citations engaging whether "opposing-direction" is genuinely outside Zscheischler's "preconditioned" class (it may be a special case). Strengthens the conceptual claim or honestly bounds it. |
| 10 | Repo working-file clutter in submission bundle | ✅ See "Submission packaging" below — list of files to exclude. |

---

## C. Questions for the author (drafted responses to include in rebuttal)

- **Q1 (out-of-sample test):** See M3 protocol #2. Draft answer: *"We have applied the frozen pipeline to N post-2019 SSW events; M/N show RR < 1, consistent with the training-period direction (binomial P = __)."* — run before resubmission.
- **Q2 (SSW-specific increment):** See M2 protocol. Provide the single interaction-IRR with CI.
- **Q3 (NAO vs Z500):** See M3 #1 disclosure. Draft answer: both are reanalysis-diagnosed; we do not claim to mechanistically separate them and now present Z500 as the empirically dominant *carrier*, with the top-down vs common-cause attribution explicitly left open and supported by published GCM experiments, not by our correlations.
- **Q4 (can any n=16 test distinguish the sub-type?):** Honest answer: **no** — state that the typology claim is therefore presented as a framework/hypothesis (ties to M1 Route 1), and the distinguishing test is specified for future larger-N work.
- **Q5 (data provenance/reproducibility):** Document the SLF data-request reference and deposit the **event-level (n = 16) processed series** underlying the primary sign test as an open aggregated CSV on Zenodo, so the headline inference is reproducible without a data request.

---

## D. 🧭 Strategic decision: target venue

My referee assessment was that the structural limits (n = 16; candidate-only typology; exploratory origin; no skill) make a *Nature Geoscience* acceptance unlikely **even after these fixes**, because the novel element cannot be powered by revision alone. Two honest paths:

1. **Attempt NG with maximal honesty** (this revision + the M3 out-of-sample test). Best case: editor finds the conceptual framing and multi-system convergence compelling enough to send out. The cover letter is written for this and explicitly invites the editor's venue judgment.
2. **Pivot to a specialist venue** — *Weather and Climate Dynamics* (EGU; ideal for the stratosphere–troposphere coupling + exploratory mechanism), *The Cryosphere*, *JGR-Atmospheres*, or *Natural Hazards and Earth System Sciences*. There the robustness work and candid uncertainty treatment are strong assets and n = 16 is in-scope. I can produce a reframed abstract + cover letter for any of these in ~minutes if you choose this path.

**Recommendation:** run the M3 out-of-sample test first. If post-2019 events confirm directionally, attempt NG (path 1) with a genuine confirmatory result in hand; if not, pivot (path 2) — which is the higher-probability publication route regardless.

---

## E. Submission packaging (exclude from upload)

These working artefacts should not be in the submission bundle: `paper/main_backup.tex`, `paper/main_temp.tex`, `paper/main_full.txt`, `paper/main_for_review.txt`, `paper/main_text_dump.txt`, `paper/r59_*.txt`, `paper/si_key_for_review.txt`, `paper/supplementary.tex`, `paper/supplementary_information_backup.tex`, root `temp_tables.tex`, `temp_abstract.tex`, `move_text*.py`, `r98*_out.txt`, `ANALYSIS_SUMMARY_R60.txt`. Submit only: `main.tex` (+ compiled PDF), `supplementary_information.tex` (+ PDF), figure files referenced in `main.tex`, `COVER_LETTER.md` (as PDF), and `REVIEWER_INDEX.md` if the journal allows a reproducibility appendix.

---

## F. Summary of edits applied in this pass

**`paper/main.tex`** — (1) header-comment title aligned to `\title{}`; (2) spec-curve sign-stability caveat; (3) Bayes-factor prior-sensitivity hedge; (4) operational-claim qualifier tying to BSS = −0.95.

**`paper/supplementary_information.tex`** — (5) pooled-vs-event-level inferential caveat on the extreme surface-weather p-values; (6) human-trigger 1.45-vs-1.81 disambiguation; (7) Norway n_eff ≈ 4 / inconclusive-replication caveat in extended abstract and plain-language summary; (8) plain-language forecasting sentence now states no calibrated skill.

**New files** — `paper/COVER_LETTER.md`; this plan.

**Not done by editing (require author):** M1 typology demotion in title/abstract (decision); M2 interaction-IRR analysis; M3 out-of-sample test + origin disclosure; M7 accident relocation (decision); minors 6, 8, 9 (figure/trim/citations); venue decision (§D); packaging (§E).
