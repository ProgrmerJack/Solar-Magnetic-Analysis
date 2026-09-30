# FINDING: the pre-onset anomaly is real, and the observational estimate is too large

**Date:** 2026-07-31
**Status:** major result, one model so far. Multi-model confirmation in progress
— the headline claim is **not** safe until a second model agrees.

---

## 1. The question observations could not answer

The pre-onset AO anomaly at −30..−16 d before an SSW is −0.82 (p=0.015) on the
frozen 43-event catalogue. It loses significance under every restriction —
isolated events −0.55 (p=0.117), satellite era −0.39 (p=0.277), both −0.28
(p=0.496) — but **every one of those intervals still contains −0.82**, and
neither mechanism test separates (era p=0.115, clustering p=0.234).

Nothing was resolved because the binding constraint is 47 winters of one
realisation of history. No estimator fixes that.

## 2. Breaking the constraint without SNAPSI

SNAPSI is behind CEDA authentication, which has been unavailable. It is not the
only source of many independent realisations. **CMIP6 daily output is on Google
Cloud as consolidated zarr with no authentication at all**, and daily `ua` is
archived on plev8 = 1000/850/700/500/250/100/50/**10 hPa** — 10 hPa is present,
which is exactly what Charlton–Polvani detection needs.

`acquire_cmip6_ensemble.py` pulls CanESM5 historical members, keeping only zonal
means (~140 MB retained per member against ~18 GB read), threaded at 35–42 MB/s.
Each member is 165 years.

Two corrections were needed to make the data usable:

- **Access pattern.** CESM2-LENS on AWS chunks `U` as (1 member, 10 days, all 32
  levels, all lat/lon), so a 10 hPa/60°N time series would require reading
  ~850 GB per member. CMIP6 on Google Cloud is chunked far more favourably.
- **Calendar.** CanESM5 uses `calendar="noleap"`. The raw time values are integer
  day counts, and casting them to `datetime64` silently produced
  nanoseconds-since-epoch — collapsing 165 years into a single "winter" and
  yielding 1 detected SSW per member. Decoding through `cftime` with the declared
  units and calendar fixed it.

Detection sanity check after the fix: **0.51–0.55 SSWs per winter**, against
~0.6 observed. The detector is behaving.

## 3. Result — 176 SSWs in 332 winters (2 members)

| AO bin | effect | 95% CI | p |
|---|---|---|---|
| −60..−46 | −0.120 | [−0.348, +0.119] | 0.320 |
| −45..−31 | −0.295 | [−0.486, −0.070] | 0.010 |
| **−30..−16** | **−0.314** | **[−0.528, −0.084]** | **0.008** |
| −15..−1 | −0.457 | [−0.658, −0.251] | <0.001 |
| +0..+14 | −0.818 | [−1.006, −0.617] | <0.001 |
| +15..+29 | −0.881 | [−1.067, −0.697] | <0.001 |
| +30..+44 | −0.502 | [−0.705, −0.293] | <0.001 |
| +45..+60 | −0.317 | [−0.525, −0.113] | 0.002 |

**The pre-onset anomaly is real.** It is significant at p=0.008 on 176 events,
with an interval 2.7× tighter than the observational one. The profile shape
matches observations: anomaly builds from −45 d, peaks at +15..+29, decays.

## 3b. Final numbers — 856 SSWs in 1,660 winters (10 members)

| AO bin | effect | 95% CI | p |
|---|---|---|---|
| −60..−46 | −0.130 | [−0.222, −0.036] | 0.006 |
| −45..−31 | −0.258 | [−0.355, −0.167] | <0.0001 |
| **−30..−16** | **−0.358** | **[−0.455, −0.261]** | **<0.0001** |
| −15..−1 | −0.376 | [−0.463, −0.292] | <0.0001 |
| +0..+14 | −0.656 | [−0.750, −0.564] | <0.0001 |
| +15..+29 | −0.609 | [−0.696, −0.523] | <0.0001 |
| +30..+44 | −0.469 | [−0.554, −0.387] | <0.0001 |
| +45..+60 | −0.332 | [−0.422, −0.240] | <0.0001 |

20× the observed events, 46× the winters, interval **6.3× tighter**. The estimate
is stable as members are added: −0.314 (2 members) → −0.358 (10 members).

## 4. The magnitude comparison — and a correction

**An earlier version of this document claimed the observed −0.820 is "too large"
and not explicable as sampling noise (subsample test, P=0.007). That claim was
wrong, and it was wrong for a reason worth recording: it compared an observation
against a model without asking whether the model is unbiased.**

CanESM5 understates the surface response at **every** lag, not just before onset:

| bin | observed | CanESM5 | obs/model |
|---|---|---|---|
| −30..−16 | −0.820 | −0.358 | 2.29 |
| +0..+14 | −0.789 | −0.656 | 1.20 |
| +15..+29 | −1.026 | −0.609 | 1.69 |
| +30..+44 | −0.867 | −0.469 | 1.85 |
| +45..+60 | −0.743 | −0.332 | 2.24 |
| | | **mean over post-onset** | **1.74** |

This is the documented CMIP weak stratosphere–troposphere coupling bias. Scaling
the model precursor by the ratio estimated from the *post-onset* bins — a
quantity that does not involve the precursor at all:

| | estimate | 95% CI |
|---|---|---|
| bias-adjusted model pre-onset | −0.624 | [−0.792, −0.455] |
| observed pre-onset | −0.820 | [−1.413, −0.202] |
| | | **intervals overlap** |

So the observed precursor is **consistent** with the ensemble once the model's
general weak-coupling bias is accounted for. The precursor is somewhat more
understated than the response (2.29 against 1.74), but not enough to claim
anything specific to the precursor.

The subsample test in §4a below remains valid *conditional on the model being
unbiased*, which it is not. It is retained because the conditional statement is
still informative, not because it supports the withdrawn claim.

## 4a. Subsample test (conditional on an unbiased model)

The observed −0.820 lies **outside** the ensemble interval [−0.528, −0.084].
That comparison has two readings, and they must be separated rather than
assumed:

- (a) the observation is noisy — 43 events is few
- (b) the model understates the coupling

**Subsampling settles it.** Drawing 36-winter blocks from the ensemble, where
the truth is known to be −0.314:

| | value |
|---|---|
| subsample mean | −0.277 |
| subsample sd | 0.193 |
| 5th–95th percentile | [−0.599, +0.037] |
| **P(estimate ≤ −0.820 \| truth = −0.314)** | **0.007** |

Reading (a) is rejected: −0.820 is not ordinary 43-event sampling noise. So
either CanESM5 understates the precursor — CMIP models are documented to have
weak stratosphere–troposphere coupling — or the observational estimate is
inflated by something beyond sampling variance.

## 4b. WITHDRAWN: the second model disagrees in sign

**The claim in §5 below that "the precursor is real" is withdrawn.** It rested on
one model. MIROC6 — independent, 128×256 against CanESM5's 64×128, acquired and
analysed identically — gives the opposite sign:

| | pre-onset −30..−16 | post-onset +15..+29 | n events |
|---|---|---|---|
| CanESM5 (10 members) | **−0.387** [−0.483, −0.292] | −0.574 [−0.663, −0.489] | 856 |
| MIROC6 (2 members) | **+0.380** [+0.194, +0.573] | +0.055 [−0.133, +0.259] | 192 |

Both are significant at p<0.0001 in their own ensembles, in **opposite
directions**. Between-model spread vastly exceeds within-model sampling
uncertainty. CMIP models do not settle the precursor question; they disagree
about it.

Checks performed before accepting this as a real disagreement rather than a bug:

- **Detection validated in both.** Composite u10 at 10 hPa/60°N across onsets:
  CanESM5 +19.4 → +1.0 → +9.0 m/s; MIROC6 +25.3 → +1.2 → +4.3 m/s. Both are
  genuine vortex reversals, at 0.51 and 0.58 events/winter against ~0.6 observed.
- **Fill contamination fixed.** MIROC6's first download had 26% of cells
  corrupted by `_FillValue` averaged into the zonal mean (a 10 hPa wind of
  2.6×10¹⁸ m/s). Re-acquired with masking; u10 now +10.5 m/s mean. The
  disagreement **survives** the fix, so it is not the fill bug.
- **Index definition ruled out.** The fixed-band annular index gave MIROC6 a
  vortex–surface correlation of only +0.073. Replacing it with the leading EOF of
  zonal-mean SLP, computed per model, raised that to +0.193 and left CanESM5's
  headline essentially unchanged (−0.358 → −0.387). MIROC6 still gives +0.380.
- **SLP climatology sane in both**, with the expected subpolar minimum.

## 4c. A separate finding: the estimator contrast is system-dependent

MIROC6's **raw** post-onset composite is −0.535, but the within-winter estimator
returns +0.055. The winter fixed effect removes the entire response, which means
that in MIROC6 the SSW–AO association lives at the **winter** timescale, not the
event timescale: SSW-hosting winters are anomalous throughout rather than after
onset.

This qualifies `FINDING_design_robustness.md`. Estimator choice moves the
observed estimate by 0.002 σ, and CanESM5 behaves the same way — but in MIROC6
the conventional and within-winter estimators disagree completely. So "estimator
choice does not matter" is true **of the observed system**, not as a general
statement. Where the signal is winter-scale, the choice is decisive.

## 5. What is and is not established

**SUPERSEDED BY §4b — read that first.** The paragraph below was written from
CanESM5 alone and its conclusion does not survive the second model.

~~**Established, and this is the result.** The pre-onset surface anomaly before
SSWs is a **real physical feature**, not a design artefact. 856 model events give
−0.358 [−0.455, −0.261] at p<0.0001, with a coherent profile that builds from
−60 d, peaks just after onset, and decays. This project could not settle it from
43 observed events across many attempts; the ensemble settles it.~~

**What actually stands:** within CanESM5 the precursor is negative, large, and
precisely estimated; within MIROC6 it is positive and precisely estimated. The
observational question remains open, and the ensemble work has established *why*
it is hard rather than answering it — the quantity is model-structural, so
adding realisations of one model cannot resolve it.

**Established.** The observational estimate of −0.820 is consistent with the
ensemble once CanESM5's general weak-coupling bias (1.74×, estimated from
post-onset bins only) is applied. There is no evidence the observational
precursor estimate is anomalous.

**Not established.** Whether the 1.74× bias factor is a CanESM5 property or
shared across models. MIROC6 (independent, 128×256 vs 64×128) is downloading and
provides that check. If MIROC6's obs/model ratio is also ≈1.7, the weak-coupling
bias is shared and the adjustment stands; if MIROC6 matches observations
directly, the bias is CanESM5-specific and the adjustment above is unnecessary
(the conclusion — precursor real, observation consistent — is unchanged either
way, only the route to it).

MIROC6 (5 members, independent model, different resolution 128×256) is
downloading. Ten CanESM5 members will also tighten the ensemble estimate from
2 members' worth.

**Caveat that does not go away.** These are models. Model coupling strength is a
model property. What the ensemble establishes cleanly is the *existence* and
*shape* of the precursor under a design with abundant realisations; the
*magnitude* comparison to observations inherits every model bias.
