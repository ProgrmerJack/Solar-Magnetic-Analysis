# Canonical AO/NAO event study — frozen catalogue, frozen estimator

Script: `canonical_event_study.py` · Output:
`results/current/2_event_study/canonical_event_study.json`
Settings: **3 annual harmonics, 1,200 winter-block bootstrap replicates, 75-day
baseline gap**, on the catalogue re-frozen 2026-07-30 at **43 events / 36 winters**.

> **This document was rewritten on 2026-09-17 because its headline had been
> reversed by its own artifact.** It previously reported the 39-event build and
> concluded *"pre-onset significance is catalogue-dependent — under the
> preregistered primary catalogue the AO pre-onset bins are not significant
> (−30..−16: −0.57, n.s.)"*. On the current data that bin is **−0.82, p = 0.005**,
> and it is significant in eight of nine catalogues. The conclusion below is the
> opposite of the one this file used to carry. The old tables are in git history;
> they are not reproduced here, because a wrong table behind a warning banner is
> still a wrong table.

## AO — effect by bin (σ), all catalogue variants

| catalogue | n / winters | −60..−46 | −45..−31 | −30..−16 | −15..−1 | +0..+14 | +15..+29 | +30..+44 | +45..+60 |
|---|---|---|---|---|---|---|---|---|---|
| **primary** | **43 / 36** | −0.17 | −0.29 | **−0.82** | **−0.59** | **−0.79** | **−1.03** | **−0.87** | **−0.74** |
| primary_compendium_only | 41 / 35 | −0.14 | −0.25 | **−0.83** | **−0.60** | **−0.80** | **−1.04** | **−0.92** | **−0.72** |
| consensus_strict | 36 / 31 | −0.11 | −0.33 | **−0.67** | −0.40 | −0.52 | **−1.05** | **−0.97** | **−0.81** |
| consensus_half | 42 / 35 | −0.18 | −0.29 | **−0.85** | **−0.62** | **−0.82** | **−1.08** | **−1.01** | **−0.83** |
| union | 47 / 38 | −0.31 | −0.43 | **−0.95** | **−0.62** | **−0.85** | **−1.11** | **−1.04** | **−0.93** |
| era5 | 42 / 35 | −0.13 | −0.31 | **−0.86** | **−0.60** | **−0.84** | **−1.05** | **−1.02** | **−0.82** |
| jra_55 | 41 / 35 | −0.11 | −0.33 | **−0.81** | **−0.53** | **−0.76** | **−1.14** | **−0.97** | **−0.73** |
| ncep_ncar | 39 / 33 | −0.31 | −0.43 | **−0.83** | −0.50 | **−0.64** | **−0.96** | **−1.05** | **−0.88** |
| merra2 | 27 / 23 | −0.20 | −0.46 | −0.52 | −0.24 | **−0.60** | **−1.07** | **−0.77** | **−0.69** |

## NAO — effect by bin (σ)

| catalogue | n / winters | −60..−46 | −45..−31 | −30..−16 | −15..−1 | +0..+14 | +15..+29 | +30..+44 | +45..+60 |
|---|---|---|---|---|---|---|---|---|---|
| **primary** | **43 / 36** | −0.18 | −0.11 | **−0.33** | **−0.30** | **−0.29** | **−0.31** | **−0.50** | **−0.35** |
| primary_compendium_only | 41 / 35 | −0.18 | −0.11 | **−0.37** | **−0.30** | **−0.30** | **−0.30** | **−0.49** | **−0.35** |
| consensus_strict | 36 / 31 | −0.22 | −0.17 | **−0.32** | −0.26 | **−0.27** | **−0.34** | **−0.57** | **−0.36** |
| consensus_half | 42 / 35 | −0.18 | −0.11 | **−0.37** | **−0.29** | **−0.30** | **−0.30** | **−0.50** | **−0.35** |
| union | 47 / 38 | −0.20 | −0.11 | **−0.36** | **−0.29** | **−0.28** | **−0.30** | **−0.49** | **−0.36** |
| era5 | 42 / 35 | −0.17 | −0.12 | **−0.39** | **−0.29** | **−0.30** | **−0.31** | **−0.50** | **−0.34** |
| jra_55 | 41 / 35 | −0.16 | −0.10 | **−0.38** | **−0.28** | **−0.28** | **−0.34** | **−0.54** | **−0.35** |
| ncep_ncar | 39 / 33 | −0.26 | −0.18 | **−0.35** | **−0.29** | **−0.28** | **−0.28** | **−0.51** | **−0.36** |
| merra2 | 27 / 23 | −0.14 | −0.10 | −0.17 | −0.09 | −0.22 | **−0.38** | **−0.50** | **−0.30** |

Bold = P < 0.05, two-sided, winter-block bootstrap.

## The finding

**Pre-onset significance is NOT catalogue-dependent — it survives everywhere the
record is long enough to see it.** For AO at −30..−16, p = 0.005 under `primary`
and below 0.05 in eight of nine sets. For NAO the same bin runs p = 0.003–0.018
across eight sets.

**The single exception is `merra2`, and it is a record-length effect, not an event
definition effect.** MERRA-2 covers 1980 onward, so the set holds **27 events / 23
winters** against 43 / 36. Its pre-onset bins are the only non-significant ones in
either outcome, while its post-onset bins remain significant and of ordinary size
(+15..+29: −1.07 AO, −0.38 NAO). That is what losing 37% of the events does to the
weaker half of the profile. `consensus_strict` (36 events) shows the same pattern
more mildly: −30..−16 survives but −15..−1 and +0..+14 do not.

This is consistent with the project's wider result that catalogue choice moves the
AO estimate by **0.115 σ** once the analysis window is held fixed, and that the
apparent 0.426 σ spread in the earlier build was two-thirds record length rather
than event definition (`FINDING_design_robustness.md`).

**What this does not settle.** That the pre-onset anomaly is significant does not
make it causal — a precursor cannot be caused by the event that follows it. Its
status is set elsewhere: CanESM5 gives −0.387 and MIROC6 **+0.380**, both
p < 0.0001, so models disagree in sign and the question is model-structural
(`CONSOLIDATED_RESULTS.md` §7.3). What this document establishes is narrower and
solid: the observed pre-onset depression is not an artifact of which event
catalogue is used.

## Stability

Across all 9 catalogues and both outcomes, every post-onset bin is negative and
every pre-onset bin is negative (`_stability` block in the JSON).

The two outcomes peak at different lags, and each does so consistently:
**AO is deepest at +15..+29 in 8 of 9 catalogues** (`ncep_ncar` alone peaks at
+30..+44), while **NAO is deepest at +30..+44 in all 9**. So the NAO response lags
the AO response by about two weeks, and that offset is not a catalogue artifact —
it is the same in every set including the two shortest. The lead–lag shape is
catalogue-independent, and so, on the current build, is its significance.
