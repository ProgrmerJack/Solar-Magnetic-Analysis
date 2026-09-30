# M1 / M2 / M3 — Implemented analyses and results

New scripts (all runnable from repo root), outputs in `data/results/`:

| Fix | Script | Output JSON |
|-----|--------|-------------|
| M2 | `scripts/analysis_extra/r99_blocking_ssw_interaction.py` | `r99_blocking_ssw_interaction.json` |
| M3 | `scripts/analysis_extra/r100_oos_post2019.py` | `r100_oos_post2019.json` |
| M1 | `scripts/analysis_extra/r101_irreducibility_power.py` | `r101_irreducibility_power.json` |

---

## M2 — Is there suppression *beyond* ordinary Alpine blocking?  ✅ ANSWERED

Negative-binomial GLM of daily natural dry-slab counts (`dry_natural_size_1234`, the
primary endpoint) on **blocking × SSW interaction**, with winter fixed effects, a
day-of-season quartic, and a winter-block bootstrap (B = 1000). Blocking defined as the
paper's own operational index (standardized NH Z500 anomaly > 1σ; same as script 56).
Sample: 2,702 winter-days, 21 winters, 477 blocking-days, 465 SSW-window days,
**only 25 days that are both**.

**Headline (full binary-blocking model):**
- **IRR_S = 0.88, 95% CI [0.28, 1.43]** — SSW effect controlling for blocking; point
  estimate ~12% suppression but **CI includes 1 → not statistically resolved beyond blocking.**
- IRR_B (blocking | non-SSW) = 0.71 [0.30, 1.37].
- IRR_BS (interaction) = 1.97 [0.14, 4.65] — **unidentifiable** (only 25 co-occurrence days).
- NB dispersion α = 14.2 (strong overdispersion → NB warranted).

**Continuous-blocking sensitivity** (uses full Z500 strength, not a 1σ cut):
- IRR per +1σ blocking = 0.91 [0.75, 1.11]; **IRR_S (SSW | mean blocking) = 1.09 [0.66, 1.79]** — essentially null.

**Nested progression of IRR_S** (shows where the apparent effect lives):
raw 0.66 → +seasonal quartic 0.54 → +winter FE **0.88**. Between-winter structure absorbs most of it.

**Cross-check:** IRR_S = 0.88 reproduces the manuscript's independently-derived regime-conditioned
IRR = 0.89 (P = 0.06) — confirms the model is correctly specified.

**Conclusion → manuscript:** there is **no statistically resolved SSW-specific avalanche
suppression beyond Alpine blocking**; SSW onset is best read as a *predictable precursor
marker of the blocking regime*, not an independent suppressor. (Integrated into Discussion/limitations.)
*Caveat:* blocking proxy is NH-mean Z500 (the only Z500 in the panel); the within-winter daily
design is conservative. The event-level sign test (14/16) remains the primary evidence.

---

## M3 — Frozen-pipeline out-of-sample test on post-2019 SSWs  ⚠️ RUN, BUT DATA-LIMITED

Generalises `35_prospective_2021_test.py` to all post-study SSWs, same matched-control design,
reference = post-2019 **non-SSW** winters (controls the secular trend; out-of-sample vs pre-2019 training).

**SSW verification (NCEP 10 hPa wind reversal) — all three confirmed:**
- 5 Jan 2021 (u₁₀ min −21.2 m/s), 16 Feb 2023 (−25.0), 16 Jan 2024 (−9.4).
- Independent DJFM reversal scan found exactly these three mid-winter events (2022-03-19 is a March final-warming, correctly excluded).

**Results (n_OOS = 3):**
- **Accidents** (predict ↑, loaded-gun): 2/3 increase — Jan 2021 RR = 2.84, Feb 2023 RR = 0.79,
  Jan 2024 RR = 1.48; **gmRR = 1.49** (matches in-sample 1.40); sign P = 1.0.
- **Bulletin danger** (predict ↑): 2/3 increase — RR = 1.28 / 0.65 / 1.26; gmRR = 1.02; sign P = 1.0.
- **Direction consistent with the in-sample loaded-gun prediction, but statistically null at n = 3.**

**The decisive test is blocked on data, not code:**
- The **primary endpoint (natural dry-slab counts) ends May 2019** in all local data, so the
  count-channel OOS **cannot be run** — confirms the manuscript's own note. `oos_counts()` is
  implemented and will execute the instant post-2019 SLF activity data is provided at
  `data/processed/cryosphere/slf_activity_post2019.parquet`.
- **Hard physical limit:** only ~3 major mid-winter SSWs have occurred since 2019, so *no*
  post-2019 OOS test can be powered. This is supportive context, honestly framed as such.

**Conclusion → manuscript:** OOS direction matches in-sample on the available (weak) channels;
powered prospective confirmation awaits more events or post-2019 SLF count data. (Integrated into limitations.)

---

## M1 — Sample-size requirement instead of an underpowered p-value  ✅ COMPUTED + REFRAMED

Replaces the "joint test P = 0.166 / non-significant ρ" framing with the events-needed-to-confirm.

**Events for 80% power (α = 0.05):**
- Cross-arm correlation: **N = 70** for the observed |ρ| = 0.33; N = 47 for 0.40; N = 30 for 0.50.
- Joint triple-positive enrichment (0.81 vs product-of-marginals null 0.666): **N ≈ 57**.
- **Headline: ≈57–70 SSW events (3–4× the 16 available).**

**Power at n = 16:** ρ = 0.33 → 23.5%; enrichment → 16.9%. (Consistent with the observed non-significance.)

**⚠️ Manuscript numerical error found and corrected:** the text claimed "<25% power to detect
|ρ| = 0.5." The correct figure for ρ = 0.5 at n = 16 is **≈51%** (Fisher-z). The <25% applies to the
*observed* |ρ| ≈ 0.33. All instances corrected in `main.tex` (×2) and SI (×2).

**Reframing applied:** title kept (already descriptive); abstract, Results §"loaded-gun", and Discussion
now call the sub-type a **hypothesis** and state the ≈57–70-event requirement explicitly.

---

## Manuscript edits made from these results
- **Abstract:** typology → "hypothesised … requires ≈57–70 events"; "reveal" → "suggest."
- **Results (joint/irreducibility):** p-value framing → sample-size requirement; power figures corrected.
- **Discussion (compound-event extension):** "candidate" → "hypothesised"; corrected power numbers.
- **Limitations:** added the M2 interaction-GLM result (IRR_S = 0.88 [0.28, 1.43]) and the M3 OOS result (2/3 accidents, gmRR 1.49, null at n = 3; count-channel data gap).
- **SI:** joint-irreducibility caption + synthesis corrected and reframed.

Both `main.tex` and `supplementary_information.tex` recompile cleanly.

---

## R102 — Breakthrough attempt: extended SSW catalog × SNOTEL (n=39, independent)

**Goal:** break the n=16 ceiling for the *mechanism* by exploiting that SSWs are hemispheric.
Computed the full Charlton–Polvani major-SSW catalog **1979–2024 from NCEP 10 hPa winds (40 events)**
and tested the continental-US **SNOTEL network (945 stations, 1980–2026)** across the 39 events in that era.
`scripts/analysis_extra/r102_snotel_extended_dissociation.py` → `data/results/r102_snotel_extended_dissociation.json`.

**What replicated at n=39 (genuine win — independent of ERA5/SNOWPACK):**
- **Melt-day frequency ↓: d = −0.47, t-P = 0.006**
- **Rain-on-snow frequency ↓: d = −0.50, t-P = 0.004**
- Solid-precip-day fraction ↑ (sign P = 0.012); SWE accumulation rate ↑ (d ≈ +0.4).
→ The **trigger-suppression arm** of the loaded gun is confirmed at >2× the Alpine event sample,
with a physically independent observing system. This directly answers the SNOWPACK/ERA5 circularity critique (M5).

**What did NOT replicate (claims now corrected in the manuscript):**
- **"SWE increases 15/15 events" does not hold at n=39** — bulk SWE is null for continental stations;
  **total precipitation DECREASES (d = −0.49)**. SSW windows are cold *and dry* → preservation, not extra loading.
- **Continental-specificity of the trigger meteorology is not supported** — triggers drop in maritime too.
  The continental focus must come from snowpack *structure* (faceting), not trigger meteorology.

**Honest verdict:** this strengthens one arm (trigger suppression, independently, n=39) and forces honest
softening of the loading and continental-specificity claims. It does **not** make the avalanche-*outcome*
claim decisive — that remains Alpine and n=16, because no more SSW events exist in the avalanche-count era.

**Bottom line across M1/M2/M3/R102:** every higher-power test points the same way — the natural-trigger
suppression is real and now independently confirmed, but the SSW-specific / novel-compound / strong-loading
framing is not supported at power. The decisive, honest paper is a **stratosphere→circulation→avalanche-trigger
teleconnection** with an independently-confirmed trigger-suppression mechanism — not a new compound-hazard type.

---

## R103 — Monitoring blind spot + Davos trigger fractions (NG-killer attempt)

Tried to relocate the decisive/novel claim to the "monitoring blind spot" (danger exceeding
count-predicted danger) and the natural/human trigger-fraction shift, both testable as sharper questions.
`scripts/analysis_extra/r103_monitoring_blindspot.py` → `data/results/r103_monitoring_blindspot.json`.

- **Blind spot** (danger residual above count-predicted, non-SSW isotonic mapping): **11/16 events, sign P=0.11, mean +0.076 danger levels** — directionally right, NOT decisive. Once you control for the count level, the "hidden hazard" is small.
- **Davos triggers** (13,918 avalanches): natural fraction down **9/16 (P=0.40)**, human fraction up **10/16 (P=0.23)** — directionally right, NOT significant.

**Verdict:** another honest non-decisive result. Controlling properly for activity level, the decoupling is small; the trigger-fraction shift is weak. No NG-decisive claim emerges.

---

## OVERALL CONCLUSION (after M1, M2, M3, R102, R103 + extended catalog + SNOTEL + Davos)

Every rigorous, higher-power, or sharper test converges on the same answer:
- **Decisive + reducible-to-blocking:** natural-trigger suppression (14/16) and PWL deepening (15/16, P=0.0005) are real; trigger-suppression mechanism is independently confirmed at **n=39** (SNOTEL, d=-0.5, P=0.004).
- **NOT decisive at any achievable power:** SSW effect beyond blocking (M2 null), compound-type irreducibility (needs ~57 events), out-of-sample outcome (n≤3), the blind spot (11/16 P=0.11), trigger fractions (9–10/16).

The "novel + broad + decisive" combination NG requires is **not present in the data and is not reachable by honest means**, because (a) the strong effect is mediated by ordinary blocking and (b) the genuinely novel parts are capped by SSW rarity (16 events in the avalanche-count era). Further hunting would be garden-of-forking-paths p-mining — the exact pathology to avoid.
