# Finding: how much of SSW downward coupling is actually predictable

**Status: established in a 1888-event CMIP6 ensemble, out of sample, with no
correction of any kind. Observations (n=42) have 35% power and cannot settle it.**

Scripts: `07_physical_decomposition/predictability_ceiling.py`,
`07_physical_decomposition/predictability_wave_driving.py`

## Why this supersedes the variance route

`FINDING_forced_variance_ceiling.md` estimated the between-event forced variance
`sigma_f` by differencing variances against pseudo-events and subtracting a
placebo. It works, but the placebo removes **91%** of the raw signal in CMIP6, and
that correction was the weakest joint in the project.

Cross-validated prediction needs no correction. Fit on training events, score on
events never seen, and the out-of-sample R^2 **is** the genuinely predictable
share of between-event variance. Nothing to subtract, nothing to calibrate. It
also asks the field's own question in the field's own terms: *which SSWs propagate
downward?*

## The result

CMIP6, 1888 events, 20 ensemble members as CV groups, GroupKFold by member, ridge
and gradient boosting, permutation null, power curve. 86 features after the
wave-driving closure below.

| predictor tier | overlaps response window? | CV R^2 | p |
|---|---|---|---|
| **P1 pre-onset** (strictly before day 0) | no | **+0.115** | 0.000 |
| **P2 at-onset** (adds days −5..0) | no | +0.113 | 0.000 |
| **P3 + post-onset stratosphere** (days 0–30) | **yes, 22 d** | **+0.422** | 0.000 |

**There IS real skill.** About **11–12% of between-event variance in the surface
response is genuinely predictable in advance**. An earlier claim in this project
that "there is nothing to classify" was too strong and is retracted. Power:
detection rate 0.70 at injected R^2 = 0.02 and 1.00 at 0.05, so 0.115 sits far
above the floor.

**The tier contrast is the finding.** Adding post-onset stratospheric information
— never the surface — nearly quadruples apparent skill, 0.115 → 0.422.
**73% of the apparent diagnostic skill is not available in advance**, and comes
from the 22-day overlap between the classifier window and the response window.
That is this project's thesis tested by prediction rather than by projection.

## What the published splits implicitly claim

Implied R^2 of a binary split is `q(1-q)C^2/Var`, using the variance of the **same
system** the contrast was measured in (CMIP6 0.519, observations 1.117). An
earlier version divided observational contrasts by the CMIP6 variance and produced
R^2 = 1.285 — an impossible value, and the giveaway that the denominator was
wrong.

| reported contrast | C | q | implied R^2 | vs 0.115 achievable |
|---|---|---|---|---|
| CMIP6 DW−NDW | 1.143 | 0.50 | **0.629** | **5.5x** |
| Karpechko AO, observations | 1.782 | 0.70 | **0.597** | **5.2x** |
| ACP 26,3723 published NAO | 0.850 | 0.59 | 0.157 | 1.4x |
| ERA5 1000 hPa NAM | 0.684 | 0.54 | 0.104 | 0.9x |
| ERA5 850 hPa NAM | 0.639 | 0.59 | 0.088 | 0.8x |

**Not every published contrast is an overclaim on this metric, and the paper must
say so.** The two ERA5 contrasts imply almost exactly the genuine predictable
variance. Only the large ones (CMIP6 DW−NDW, observational AO) imply ~0.60, which
exceeds even the post-onset diagnostic value.

The ERA5 splits are still not vindicated — they are *outcome-based*, so they
explain the response tautologically rather than predictively, and the separate
`beta x dS` decomposition puts them at 93–102% selection. Implied-R^2 and the
projection law answer different questions; both belong in the paper.

## Independent convergence

| method | R^2 |
|---|---|
| placebo-corrected `sigma_f` upper 95% (0.2388) | 0.110 |
| cross-validated, pre-onset | 0.100 → 0.115 with wave proxies |
| cross-validated, at-onset | 0.113 |

Two routes sharing no assumptions — one needing a 91% placebo correction, one
needing none — agree to within 0.01. The correction that was the obvious referee
target is independently corroborated, and `sigma_f` sits at the **upper end** of
its bound (~0.24 sigma), not at the low point estimate of 0.081.

## The wave-driving gap, and how it was closed

CMIP6 here holds zonal-mean `u` and `psl` only (21 `*_zm.nc`, no 3-D fields), so a
true eddy heat flux would need hundreds of GB of daily `va`/`ta`. Closed from two
directions instead.

**1. Proxies from zonal-mean u.** The TEM momentum budget puts the wave forcing
`div F` directly into `du/dt` — the identity the Charlton–Polvani wind-tendency
SSW definition already rests on. Added: windowed `du/dt`, 10−100 hPa vertical
shear, and 75N−45N vortex geometry. Features 38 → 86.

| tier | without wave | with wave | gain |
|---|---|---|---|
| P1 pre-onset | +0.0998 | **+0.1148** | +0.015 |
| P2 at-onset | +0.1121 | +0.1133 | +0.001 |
| P3 post-onset | +0.4245 | +0.4220 | −0.003 |

Wave driving helps **only before onset**, which is physically right: once the
at-onset vortex state is known, the wave history is already encoded in it. The
headline moves to ~0.115 and the conclusion is unchanged.

**2. Real v'T' in observations** (100 hPa, 45–75N, 1958–2024), nested A/B at low
dimension: vortex-only CV R^2 = −0.030 (p = 0.395), vortex + wave = −0.081
(p = 0.185). No skill either way — **but the arm is uninformative, not negative.**
Power at n=42: detection rate **0.35 at a true R^2 of 0.10**, 0.58 at 0.20.

**The observational null must never be cited as evidence against predictability.**
It has 35% power against the effect CMIP6 establishes.

## Limits

- The model arm carries the claim. n=42 observations cannot settle it, and that is
  a property of the record rather than of the method.
- Wave driving in CMIP6 is a proxy, not a measured eddy flux, and is labelled as
  one throughout.
- Established in ten CMIP6 models with known spread in SSW surface-composite
  fidelity.
