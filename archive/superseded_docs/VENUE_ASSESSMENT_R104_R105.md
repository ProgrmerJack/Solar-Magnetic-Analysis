# Venue assessment: what would make this Nature Geoscience / Nature Communications

Date: 2026-07-26. Scripts: `r104_accident_case_crossover.py`, `r105_multiregion_meta.py`.

## Headline: the n=16 ceiling was NOT physical — and breaking it changed the answer

Earlier conclusion was that the compound-hazard claim was capped at 16 SSW events and
therefore unconfirmable. That was wrong about the *cause*: the cap came from the Swiss
daily-count record (1998/99–2018/19), not from SSW rarity. Longer avalanche records exist
**already in this repo**:

| dataset | span | SSWs covered | n records |
|---|---|---|---|
| Swiss Davos counts (used in paper) | 1998–2019 | 22 | 3,325 winter days |
| **CAIC US accidents** | **1951–2025** | **40 (all NCEP-era)** | 802 winter (1979+) |
| **LAWIS Tirol incidents** | **1992–2024** | **29** | 2,189 |
| US danger ratings (multi-center) | 2011–2025 | 10 | 69,871 |
| Norway NVE danger | 2018–2025 | 5 | 22,848 |

So the ceiling was breakable. It was broken. The results below are what came out.

---

## Finding 1 — the "opposing-direction / loaded gun" claim is NOT supported at higher power

The paper's most novel claim is that human-trigger hazard *rises* while natural activity falls.
Tested on the two long accident records (r104), lag-resolved, with winter fixed effects
(absorbing the ~10× growth in backcountry recreation) and DOY harmonics (absorbing seasonality):

**LAWIS Tirol — lag profile (IRR vs non-SSW winter days)**

| lag (days from onset) | IRR | 95% CI | P |
|---|---|---|---|
| −45 to −31 | 1.14 | 0.92–1.40 | 0.23 |
| −30 to −16 | 0.90 | 0.73–1.11 | 0.33 |
| −15 to −1 | 1.05 | 0.86–1.27 | 0.62 |
| **0 to +14** | **0.82** | 0.68–1.00 | 0.048 |
| **+15 to +29** | **0.74** | 0.61–0.90 | 0.002 |
| **+30 to +44** | **0.70** | 0.58–0.84 | 0.0004 |

Accidents go **down** after onset, not up — flat before, monotonically suppressed after.
That is a textbook downward-coupling signature, but it points the *opposite* way to the
paper's human-trigger claim.

**Critical artifact caught:** a single ±15 d window in Tirol gives IRR = 1.335, P = 0.004
(apparent *increase*). That is a **reference-category artifact** — the ±15 d window is the
*least* suppressed part of an entirely suppressed post-onset period, so against a reference
contaminated with strongly suppressed days (+15 to +44) it looks elevated. The lag
decomposition dissolves it. Any single-window analysis of this kind is unsafe.

## Finding 2 — model p-values here are badly anti-conservative

Tirol post-onset [+15,+44]:

| inference | IRR | 95% CI | P |
|---|---|---|---|
| model SE | 0.78 | 0.68–0.89 | **0.0004** |
| 3 harmonics / month FE | 0.76–0.77 | 0.66–0.88 | 0.0001–0.0002 |
| **winter-block bootstrap** | 0.78 | **0.52–1.14** | **0.11** |
| 2000+ (complete reporting) | 0.90 | 0.77–1.04 | 0.15 |
| fatal only (n=280) | 0.54 | 0.34–0.85 | 0.008 |
| leave-one-winter-out | 0.65–0.88 | — | all < 1 |

Avalanches cluster massively within winters (one cycle = many events), so model SEs are
wrong. **P = 0.0004 → P = 0.11.** Direction is consistent (LOWO all < 1); significance is not.

## Finding 3 — pooled across regions, the effect is null and sign-inconsistent

Random-effects meta-analysis, post-onset [+15,+44] d, winter-block-bootstrap SEs (r105):

| region | IRR | 95% CI |
|---|---|---|
| CH-Davos natural counts | 0.73 | 0.48–2.00 |
| AT-Tirol incidents | 0.78 | 0.53–1.14 |
| US-continental accidents | 1.26 | 0.85–1.73 |
| US-maritime accidents | 1.23 | 0.84–1.70 |
| **POOLED** | **1.03** | **0.79–1.35**, P = 0.80, I² = 39.5% |

Placebo (pre-onset −30 to −1): pooled 1.07, P = 0.59 — placebo behaves correctly.

**The Alps go down (~0.75), the US goes up (~1.25).** Either that geographic split is real and
mechanistically explained (SSW-driven blocking is Atlantic/European-centred), or the pooled
claim collapses. It is currently unexplained.

## Finding 4 — the paper's headline is a BETWEEN-winter effect, not a within-winter one

This is the most important result for reviewability. Using the **paper's own 16-event catalog**
and the **paper's own ±15 d window** on Davos natural dry-slab counts:

| design | estimate | inference |
|---|---|---|
| **between-winter** (paper's: SSW windows vs DOY-matched days in non-SSW winters) | gmRR = **0.305**, 14/16 down | **P = 0.0021** ✔ reproduces the paper exactly |
| **within-winter** (winter fixed effects, same window) | IRR = **0.68** | bootstrap CI 0.30–1.31, **P = 0.124** |

The control pool is **7 non-SSW winters** (1,100 days). Winter-to-winter variability in
avalanche activity is enormous, so the headline rests on a thin between-winter contrast.
When each winter serves as its own control, the effect **halves** (70% → 32% reduction) and
loses significance. SSW winters are not quieter overall (ratio 1.35, Mann-Whitney P = 0.80).

The direction survives everywhere. The magnitude and the p-value do not.

---

## Venue verdict

**Novelty is real and confirmed.** A literature check found no existing SSW→avalanche work;
the closest analogue (Arctic vortex collapse → South China rainfall, 2025) went to
*npj Climate and Atmospheric Science*, not to NG.

| criterion | NG | Nat Comms |
|---|---|---|
| Novelty | ✔ first stratosphere→avalanche link | ✔ |
| Decisive support | ✘ headline is between-winter only; pooled null; sign-inconsistent | ✘ same |
| Broad significance | ✘ mechanism reduces to known SSW→blocking (M2 null) | ~ borderline |
| Format | ✘ 10,771 words / 17 figures vs ~3,000 / 4 | ✘ vs ~5,000 |

**Neither is reachable with the current evidence.** Not because the science is bad — the
direction is consistent across four independent observing systems — but because the decisive
claim rests on a 7-control-winter contrast that halves under within-winter controls, and the
most novel claim (opposing-direction) is contradicted by the better-powered data.

## The concrete path to Nature Communications

The binding constraint is **control winters**, not SSW events. Fix that and the paper changes class.

1. **Get the full SLF Swiss record (1970–present).** 21 → ~55 winters; 7 → ~25 non-SSW control
   winters. This directly attacks the identification weakness in Finding 4. Request template
   already drafted in `scripts/download/06_download_caic_avalanche.py`.
2. **Add the remaining Alpine regions** — LAWIS full Austria (public API; currently only Tirol),
   France EPA/CLPA (1900–present, longest avalanche record in existence), Italy AINEVA,
   South Tyrol/Trentino. Pooling ~5 Alpine regions is what buys decisiveness: single-region CI
   half-width ~0.73 in log space → ~0.33 pooled → CI ≈ 0.52–0.94, i.e. significant.
3. **Pre-register the within-winter design** (winter FE + DOY harmonics + winter-block bootstrap
   + lag profile + pre-onset placebo). Make the within-winter estimate the headline, not the
   between-winter one.
4. **Drop the opposing-direction / loaded-gun framing.** Finding 1 contradicts it. Keep PWL
   deepening as a mechanism observation, not a hazard-type claim.
5. **Explain or drop the Alps–US sign split** using the ERA5 blocking composite region by region.
6. **Cut to ~5,000 words / 6 figures.**

Honest probability: with steps 1–3 delivering a pooled Alpine IRR whose CI excludes 1,
this is a credible Nature Communications submission. Without new data, it is not — and no
reanalysis of the current data will change that, because the limit is control winters.

**Publishable right now, without any new data:** *Journal of Glaciology*, *The Cryosphere*,
*Weather and Climate Dynamics*, or *npj Climate and Atmospheric Science* — reframed around
the consistent-direction, multi-system suppression signal with the design caveats stated openly.
