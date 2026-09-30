# FINDING: mediator measurement error under-attributes the stratosphere

**Date:** 2026-07-31
**Status:** solid methodological result; the underlying physics question is NOT
settled by it, and this document says so.

---

## 1. The question

The live debate in stratosphere–troposphere coupling is whether the stratosphere
is **causal** or whether tropospheric wave activity drives both the vortex
anomaly and the surface response — the "common cause" reading. It is genuinely
open: recent work argues surface anomalies are "just as consistent with ocean,
sea ice, and land forcing," with tropospheric circulation playing "an equal if
not more important role."

Every design problem this project documented — catalogue choice, consensus
rules, era dependence, compositing, selection on the outcome — comes from
discretising a continuous process into "events." So this analysis drops events
entirely and works on the continuous record: **6,958 winter days across 47
winters**, against 43 events in 36 winters.

## 2. The chain

| | variable | timing |
|---|---|---|
| X | 100 hPa eddy heat flux v′T′, 45–75°N | days t−15..t−1 |
| M | stratospheric vortex state | day t onward |
| Y | AO | days t+15..t+45 |

The ordering is physical, not assumed from correlation: v′T′ at 100 hPa is a
tropospheric flux of wave activity *into* the stratosphere and necessarily
precedes the stratospheric response it produces.

Observed structure: corr(X,M) = −0.437 (waves decelerate the vortex),
corr(M,Y) = +0.185 (weak vortex → negative AO), corr(X,Y) = −0.149.

## 3. The result that looked like a breakthrough, and was not

With the vortex measured as **10 hPa 60°N wind on a single day t**:

| | estimate | CI | p |
|---|---|---|---|
| total X→Y | −0.152 | [−0.224, −0.076] | <0.0001 |
| direct X→Y \| M | −0.130 | [−0.218, −0.044] | 0.003 |
| indirect (via vortex) | −0.022 | [−0.087, +0.041] | 0.50 |
| **mediated share** | **14.2%** | [−35%, +60%] | — |

Read naively: 86% direct, stratosphere a bystander, common-cause confirmed.
That would be a major claim.

**It is an artefact.** Classical measurement error in a *mediator* attenuates the
indirect path toward zero and inflates the apparent direct path — exactly the
observed pattern. A one-day wind snapshot is a noisy proxy for an anomaly that
persists for weeks and descends over 10–20 days.

Varying only the mediator definition:

| mediator | mediated share | overlaps outcome window? |
|---|---|---|
| u10, single day | 14.2% [−35, +60] | no |
| u10, mean t..t+10 | **45.2% [+8, +88]** | no |
| u10, mean t..t+20 | 74.4% | yes — not interpretable |
| u10, mean t..t+30 | 89.7% | yes — not interpretable |
| 10 hPa polar-cap height, t..t+20 | 103.8% | yes — not interpretable |

Only the first two are interpretable; windows extending past t+15 overlap the
outcome and induce mechanical correlation. Between them the share triples.

## 4. Correcting it properly

Two imperfect measures of one latent vortex state, both ending at t+10 so
neither overlaps the outcome window:

- M1 = 10 hPa 60°N zonal wind, mean t..t+10
- M2 = 100 hPa polar-cap height (sign-flipped), mean t..t+10
- corr(M1,M2) = +0.722

Instrumenting M1 with M2 under classical errors-in-variables:

| quantity | estimate | CI95 |
|---|---|---|
| **reliability of M1** | **0.592** | [0.521, 0.663] |
| latent vortex → surface (IV) | +0.169 | [−0.028, +0.369] |
| mediated share, naive | 45.2% | [8%, 87%] |
| mediated share, disattenuated | 53.9% | [−11%, 140%] |

Even a 10-day mean of the standard vortex index carries **41% error variance**
as a measure of the latent state that couples to the surface.

**Assumption, not waved away:** the correction needs the measurement errors of
M1 and M2 to be uncorrelated. u10 and polar-cap height are linked by thermal
wind, so that is not guaranteed. Correlated errors bias the IV estimate *toward*
the naive one — so 53.9% is a **lower bound** on stratospheric mediation, which
is the conservative direction for the argument made here.

## 5. What is and is not established

**Established.** Mediator measurement error systematically under-attributes the
stratosphere. Moving from a single-day index to a 10-day mean moves the answer
by **31 percentage points**; disattenuating moves it a further ~9. Any study
estimating stratospheric influence with an instantaneous or noisy vortex index
will understate it, and the size of that understatement is large enough to flip
a qualitative conclusion. Published findings of weak stratospheric influence
should be checked for this before being read as evidence for common cause.

**Not established.** The causality question itself. The best-supported estimate
is 45–54% mediated, but the disattenuated interval spans −11% to 140%. With 47
winters of reanalysis this cannot be narrowed usefully — the binding constraint
is that one realisation of history contains a limited number of independent
winters, not the estimator.

Settling it needs experimental data where the stratospheric state is *set*
rather than observed: SNAPSI's nudged-to-observed minus nudged-to-climatology
contrast, 50 members × 11 centres. See `GO_NO_GO.md`; blocked on CEDA
authentication, not on method.

## 6. Honest placement

This is a real contribution to an active debate and a concrete warning with a
quantified correction. It is a strong *WCD* / *JGR-Atmospheres* / *J. Climate*
result.

It is **not** a Nature Geoscience breakthrough, and should not be dressed as
one. It does not settle the causal question — it shows that one common way of
answering it is biased, and in which direction.
