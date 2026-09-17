> # ⚠ SUPERSEDED — AND ITS HEADLINE IS NOW REVERSED. DO NOT CITE AS WRITTEN.
>
> Written 2026-07-28, **before** the catalogue was re-frozen on 2026-07-30
> (39 events → 43; see `../02_event_catalogues/FREEZE_RECORD.md`). Every table
> below is the 39-event / 6-harmonic / 800-replicate build.
>
> The conclusion below — *"under the preregistered primary catalogue the AO
> pre-onset bins are not significant (−30..−16: −0.57, n.s.)"* — **does not hold
> on the current data.** `results/canonical_event_study.json` (43 events, 3
> harmonics, 1,200 replicates) gives **−0.8195, p = 0.005 — significant**, and
> that bin is significant in every catalogue except `merra2`.
>
> The "pre-onset significance is catalogue-dependent" finding built on it is
> therefore not supported. See `../CONSOLIDATED_RESULTS.md` §5b.

# Canonical AO/NAO event study — frozen catalogue, frozen estimator

Script: `canonical_event_study.py` · Output: `canonical_event_study.json`
**Re-run under the closed Gate 3 rules: 3 annual harmonics, 1,200 winter-block bootstrap
replicates.** Supersedes the earlier run (6 harmonics, 800 replicates), which was
provisional because it predated the replicate floor and the seasonal-model decision.

## AO — effect by bin (σ), frozen primary catalogue

| catalogue | n / winters | −60..−46 | −45..−31 | −30..−16 | −15..−1 | +0..+14 | +15..+29 | +30..+44 | +45..+60 |
|---|---|---|---|---|---|---|---|---|---|
| **primary** | **39 / 33** | −0.10 | −0.27 | −0.57 | −0.39 | −0.51 | **−0.97** | **−0.82** | **−0.78** |
| compendium-only | 37 / 32 | −0.07 | −0.24 | −0.58 | −0.40 | −0.52 | **−0.99** | **−0.88** | **−0.76** |
| union (old 47-event set) | 47 / 38 | −0.31 | −0.43 | **−0.95** | **−0.62** | **−0.85** | **−1.11** | **−1.04** | **−0.93** |
| ERA5 | 42 / 35 | −0.13 | −0.31 | **−0.86** | **−0.60** | **−0.84** | **−1.05** | **−1.02** | **−0.82** |

Bold = P < 0.05. NAO primary: −0.21, −0.14, −0.24, **−0.26**, **−0.25**, **−0.33**, **−0.53**, **−0.35**.

## Stability against the estimator change

Changing from 6 to 3 harmonics and 800 to 1,200 replicates moved the AO primary profile by
at most 0.03 σ in any bin (e.g. +15..+29: −0.95 → −0.97) and changed no significance
verdict. The result is not sensitive to the seasonal-model choice within the calibrated
family — which is the reassuring outcome, given the choice was made on calibration rather
than on fit to this outcome.

## The finding, unchanged and now on frozen settings

**Pre-onset significance is catalogue-dependent.** Under the preregistered primary
catalogue the AO pre-onset bins are not significant (−30..−16: −0.57, n.s.). Under the
union catalogue — the 47-event set the superseded work used — the same bin is significant
(−0.95, P < 0.05).

The union admits events detected by a single reanalysis; if such events are marginal or
mis-dated, their true response can fall inside the window labelled "pre-onset", which would
produce exactly this pattern. Testable, not established.

NAO differs in detail: its −15..−1 bin **is** significant under primary while −30..−16 is
not, so some pre-onset structure survives the strict catalogue in NAO but not in AO.

## Gate 5 condition 7

Post-onset bins negative in **every** catalogue, both outcomes; pre-onset bins negative in
every catalogue (sign stable); profile deepest at +15..+29 everywhere. The substantive
lead–lag structure is catalogue-independent, so condition 7 is **met for the shape**. The
significance of the pre-onset bins is not stable, and that is reported as a result rather
than smoothed over.
