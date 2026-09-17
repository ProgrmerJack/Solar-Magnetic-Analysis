# Finding: the stratifier bias law, and the collapse of the S1/S2 distinction

**Status: established on observations AND on a 1888-event CMIP6 ensemble.
Two errors of my own are corrected here; one candidate positive result is killed.**

Scripts: `06_simulation_validation/stratifier_bias_law.py`,
`stratifier_inference.py`, `stratifier_law_cmip6.py`.
Data: `stratifier_bias_law.json`, `stratifier_inference.json`,
`stratifier_law_cmip6.json`.

## What was previously claimed

`selection_projection_law.py` derived, **in simulation only**, that selecting SSWs
on the surface outcome biases the composite by

    bias(Y) = beta(Y_window on A_window) x bias(A)                          (1)

with beta the ordinary least-squares SLOPE. R^2 = 0.999 across five outcomes,
generalising to ENSO at R^2 = 0.985.

`AUDIT_PROTOCOL.md` then split the literature into **S1** (outcome-classified,
predicted biased) and **S2** (stratospheric or pre-onset classification, predicted
**unbiased**), and declared S2 the falsification control arm.

## Correction 1 — S2 is not a control arm. The dichotomy was wrong.

Equation (1) never required A to be the outcome. beta is defined for **any**
stratifier, and a stratospheric stratifier is not independent of the surface --
that correlation is the entire subject of the field. So stratospheric
classification carries the same bias with a smaller beta.

Measured on observations, under a true null (pseudo-onsets on event-cleaned days,
so the differential effect is exactly zero by construction):

| stratifier | class | beta | contrast manufactured under a TRUE NULL |
|---|---|---|---|
| `z100_post` lower-stratospheric descent (PJO-like) | stratospheric | +0.775 | **+0.889 sigma** |
| `u10_post45` reversal persistence | stratospheric | +0.777 | +0.732 |
| `u10_post30` reversal depth | stratospheric | +0.723 | +0.693 |
| `AO_post` the published surface criterion | outcome | +1.824 | +1.814 |
| `u10_pre` vortex preconditioning | pre-onset | +0.089 | +0.011 |

Splitting **random winter dates with no SSW at all** by lower-stratospheric vortex
strength produces a +0.889 sigma surface contrast. Real SSWs split the same way
give +0.909 sigma. The event adds 0.020 sigma.

This means the ACP 24, 1389, 2024 result — that the applied literature uses the
Hitchcock PJO *stratospheric* classification rather than a surface one — does
**not** exempt that literature. My previous conclusion that "D's scope collapses"
was wrong and is withdrawn.

## The law, validated outside simulation

Across stratifiers, under the null, predicted vs observed contrast:

| system | n stratifiers | slope | intercept | R^2 |
|---|---|---|---|---|
| observations (29-43 events) | 7 | 0.961 | +0.002 | 0.9945 |
| CMIP6 (1888 events, 10 models) | 4 | 1.001 | +0.002 | 1.0000 |

Two independent systems, structurally different stratifiers, beta spanning
-0.97 to +1.82. The law is not an artefact of the simulation that produced it.

## Correction 2 — a sample-size confound I re-created and caught

The first CMIP6 run compared a residual pooled over 1888 events against a null
built from ~90-event draws. That null is sqrt(1888/90) ~ 4.6x too wide, making
every p-value spuriously conservative: `AM_post` appeared bounded to +/-0.099 when
the matched bound is +/-0.008. This is the **same** confound that produced the
phantom "20-27% bootstrap under-coverage" earlier in this project. Null
realisations are now pooled across members at matched draw index. All CMIP6
numbers below are the matched ones.

## What the ensemble buys: the residual is small but REAL

| stratifier | contrast | selection term | residual | p | share explained |
|---|---|---|---|---|---|
| `u10_pre` | +0.162 | +0.192 | -0.030 | 0.267 | zero to +/-0.051 |
| `u10_post30` | +0.352 | +0.251 | **+0.100** | <1e-4 | 72% |
| `u10_post45` | +0.320 | +0.221 | **+0.099** | <1e-4 | 69% |
| `AM_post` | +1.143 | +1.097 | **+0.046** | <1e-4 | 96% |

So the honest statement is **not** "the differential effect is zero". It is:
roughly 70% of the between-group surface contrast for a stratospheric
classification, and 96% for an outcome-based one, is the stratifier-surface
regression rather than a differential effect of the event. The remainder is small
and, in models, statistically real.

### Two limits stated in the open, not buried

1. **`AM_post` is near-tautological.** Its classification window (8,52) IS the
   outcome window, so "96% explained" is close to an identity. It is faithful to
   the published practice and it is the correct diagnosis of that practice, but it
   is not a discovery. The load-bearing results are the stratospheric stratifiers,
   where the stratifier and outcome are different variables in different windows.
2. **Post-onset stratospheric residuals are an upper bound on the artefact, not a
   causal effect.** The event partly *causes* those stratifiers, so beta x dS may
   contain a mediated causal path; subtracting it can remove signal as well as
   bias. The artefact share for `u10_post*` is therefore **at most** 69-72%. Only
   `u10_pre` (which the event cannot have caused) and `AM_post` (where the
   stratifier is the outcome) are cleanly interpretable.

## Killed: `u10_pre` as a forecast-usable positive result

On observations `u10_pre` looked like the constructive half — contrast -0.775
against a null of +0.011, the only pre-onset and therefore forecast-usable
stratifier behaving that way. It does not survive:

- **window sign flip**: -0.706 at (-45,-16) but **+0.178** at (-60,-31);
- **seasonal confound**: high-vortex events onset 29.6 days later in the season
  than low-vortex events;
- **q_residual = 0.558, q_slope = 0.598** under BH-FDR over the pre-specified
  7-stratifier family;
- CMIP6 independently bounds its residual to **+/-0.051** at n=1745.

Observations and the ensemble agree it is nothing. Recorded as exploratory-only
and not carried forward. This is the fifth headline of this project to die on a
robustness test, and the test is why.

## The premise, VERIFIED from source 2026-08-01

Karpechko et al. (2017), QJRMS 143, 1459, doi 10.1002/qj.3017. An SSW is
"downward propagating" if, over days +8..+52 after onset: (1) mean NAM at
**1000 hPa** is negative; (2) the fraction of those 45 days with negative
**1000 hPa** NAM exceeds 0.5; (3) the fraction with negative 150 hPa NAM exceeds
0.7. Conditions 1 and 2 are the surface response. The abstract is explicit:
events are distinguished "based on whether they were followed by an anomalous
Northern Annular Mode response."

**229 citations and live**: 30 (2020), 37 (2021), 54 (2022), 22 (2023), 28 (2024),
14 (2025), 9 (2026 partial).

My earlier worry that the criterion might be stratospheric — raised because ACP 24,
1389 (2024) uses the Hitchcock PJO classification — is resolved. Both exist; the
Karpechko surface criterion is the one in active use, and the PJO one is not exempt
either (see Correction 1).

## Stratum S1 is not empty

| study | classification | presents surface composite as an effect size? |
|---|---|---|
| Karpechko et al. 2017, QJRMS | 1000 hPa NAM, +8..+52 | origin of the criterion |
| J. Climate 32, 85, 2019 | follows Karpechko, 850/100 hPa | yes — DW anomalies "around twice that of the total" |
| **ACP 26, 3723, 2026** | *"follows the method by Karpechko et al. (2017)"*, 850 hPa | yes — **NAO −0.762 (DW) vs +0.088 (NDW)**, contrast **−0.850** |

ACP 2026 additionally sub-classifies its DW events using **2 m temperature over 40
days post-onset** — conditioned on the surface twice.

## The body count

`08_literature_audit/recompute_published_criterion.py` applies conditions 1-3 as
published to the frozen 43-event catalogue (30 DW / 13 NDW, 70%):

| outcome | DW | NDW | contrast | beta x dS | residual | p | share from selection |
|---|---|---|---|---|---|---|---|
| NAO | -0.464 | +0.078 | **-0.543** | -0.464 | -0.079 | 0.443 | **85%** |
| AO | -1.433 | +0.349 | **-1.782** | -1.756 | -0.026 | 0.595 | **99%** |
| PNA | +0.186 | -0.036 | +0.222 | +0.169 | +0.053 | 0.698 | 76% |

**The assumption-free statement.** Applying the identical criterion to pseudo-onsets
on event-cleaned days, where the true DW-minus-NDW difference is exactly zero by
construction, reproduces:

- **88%** of the observed NAO contrast (null -0.476 vs observed -0.543)
- **96%** of the observed AO contrast (null -1.712 vs observed -1.782)
- 67% of the PNA contrast

This involves no counterfactual model, no beta, and no correction — two composites
computed identically, one on real onsets and one on random winter dates.

Against the protocol's pre-specified "conclusion change" rule (significance lost,
or effect size changing by >33%, or a sign flip), this is a change of 85-99%.

**Scope limit, stated plainly:** this is NOT a reproduction of ACP 2026's -0.850.
They use ERA5 850 hPa NAM and a 1000 hPa-height NAO on their own catalogue; this
uses the CPC AO and NAO on the frozen 43-event catalogue. The published value is
quoted for scale. The claim is about the decomposition of a contrast produced by
this criterion, not about any paper's digits.

## What this does and does not license

**Licensed:** that between-group surface contrasts in SSW stratification are
quantitatively predicted by beta x dS in two independent systems; that
stratospheric classification is not exempt; that the criterion in active use
produces a contrast 85-99% reproducible under a true null.

**Not licensed:** that the differential effect is zero. At n=43 the residual bound
is loose (p=0.44-0.70), and the 1888-event ensemble finds the residual small but
genuinely NONZERO. The correct claim is "mostly artefact", never "no effect".

## Next, in order

1. **Extend NCEP stratospheric fields back to 1958** to match the v'T' record.
   Observational n goes 29 -> 42 (+45%) for every stratospheric stratifier, which
   is the binding constraint on the observational arm.
2. **Run the literature audit** (`08_literature_audit/AUDIT_PROTOCOL.md`), with
   the S1/S2 strata now redefined: the question is not "did they classify on the
   surface?" but "what is beta for the stratifier they used?"
3. **Verify the Karpechko et al. (2017) criterion from the primary source.** Still
   unverified. It now matters less than it did — the bias applies either way — but
   the paper cannot describe the criterion without having read it.
