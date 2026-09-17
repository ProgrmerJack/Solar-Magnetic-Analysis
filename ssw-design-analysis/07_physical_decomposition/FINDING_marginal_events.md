> **SUPERSEDED — see FINDING_marginal_RETRACTION.md.** Monotonicity fails,
> the mechanism rests on three events, and the strong/marginal split is
> confounded with pre-satellite era (Fisher P = 0.011). Do not cite.

# The strongest result in the project: apparent SSW precursors come from marginally-detected events

Scripts: `marginal_events.py`, `design_sensitivity.py` · Output: `marginal_events.json`

## The double dissociation

Frozen catalogue split by how many of six reanalyses detected each event.
AO event study, 3 harmonics, 1,200 bootstrap replicates, **BH-FDR across all 24 tests**
(3 subsets × 8 bins):

| subset | n | pre-onset bins q<0.05 | post-onset bins q<0.05 |
|---|---|---|---|
| **strong** (≥4/6 reanalyses) | 37 | **0 / 4** | **3 / 4** |
| **marginal** (1–3/6) | 8 | **1 / 4** (−30..−16: **−1.38**, q = 0.022) | **0 / 4** |
| union (≥1/6) | 45 | 2 / 4 | 4 / 4 |

Strong events show the canonical post-onset response and **no** pre-onset anomaly.
Marginal events show a pre-onset anomaly and **no** post-onset response at all.
Pooling them — which the union catalogue does — produces both.

This is a double dissociation, not a difference in noise level. Marginal events are not
simply weaker versions of strong ones; their signal sits on the *other side* of the nominal
date.

## The physical cause

| quantity at onset | strong (n=37) | marginal (n=8) |
|---|---|---|
| 10 hPa 60°N zonal wind | **−0.62 m/s** | **+8.09 m/s** |
| 100 hPa eddy heat flux, −45..−1 d | 19.6 K m/s | 18.2 K m/s |

A major mid-winter warming is *defined* by reversal of the 10 hPa 60°N wind to easterly.
Strong events average −0.62 m/s — reversed, as required. **Marginal events average +8.09 m/s
— still westerly.** Averaged across reanalyses they do not satisfy the defining criterion.
Their precursor wave driving is indistinguishable from strong events (18.2 vs 19.6 K m/s),
so this is not a wave-forcing difference: it is a classification difference.

The natural reading is that marginal detections are vortex disturbances whose central date
is assigned by whichever one or two reanalyses happened to cross the threshold, so their
surface anomaly is not aligned with the nominal date. Including them shifts weight into the
pre-onset window.

## Why this matters beyond bookkeeping

Tropospheric precursors to SSWs are an active subfield, and "is there a significant signal
before onset?" is a question people answer by compositing. On these data the answer depends
on whether the catalogue admits marginally-detected events:

* consensus catalogue → **no** significant pre-onset AO anomaly
* union catalogue → **yes**, significant at −30..−16 d

and the difference is carried by 8 events whose vortex, on average, never reversed.

This also explains the earlier result that catalogue choice moves the estimate 18–47× more
than estimator choice: the estimator contrast is 0.008 σ, while admitting or excluding these
8 events moves the pre-onset estimate by 0.38 σ.

## What this does NOT establish — limits to state plainly

1. **n = 8 marginal events.** The −1.38 estimate has a 95% interval of [−2.25, −0.37]. It
   survives FDR, but it rests on eight events and must never be quoted without that number.
2. **The physical contrast is not statistically significant.** u10 +8.09 vs −0.62 m/s gives
   Mann–Whitney P = 0.166. The magnitude is decisive-looking; the test, at n = 8 vs 37, is
   not. It is offered as the likely mechanism, not a demonstrated one.
3. **One index.** Demonstrated on AO. NAO has not been decomposed this way.
4. **One catalogue pair.** Whether the same holds for other consensus thresholds, or in
   other reanalysis families, is untested.

## What would make this publishable at the highest level

In rough order of value:

1. **Replicate the decomposition on NAO** and on the two circulation fields — cheap, and it
   turns one index into a pattern.
2. **Test monotonicity in detection count** (1–2 vs 3 vs 4–5 vs 6). If pre-onset weight
   falls monotonically as consensus rises, that is far stronger than a two-way split.
3. **Identify the 8 events against the literature.** If they are independently known as
   ambiguous or as minor warmings, the classification argument becomes documented rather
   than inferred.
4. **Add reanalysis-specific catalogues.** If the effect appears whichever single reanalysis
   is trusted, it generalises beyond this compendium.
5. **SNAPSI**, when CEDA returns, to check the same decomposition against nudged truth.

Items 1–3 need no new data and no collaborator.
