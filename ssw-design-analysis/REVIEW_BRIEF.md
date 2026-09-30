# Review brief — SSW design/causality project

Single-author observational study. I need a strategic verdict, not encouragement.

## History: four headlines tested and withdrawn
1. Solar-magnetic forcing of avalanche hazard — abandoned.
2. Opposing-direction avalanche hazard — killed by catalogue error + seasonal confound.
3. "Conventional composites inflate SSW impacts" — killed by my own corrected analysis.
4. "Study design reshapes SSW impact estimates" — killed by direct test: estimator
   choice moves the AO estimate 0.002 sigma; catalogue choice 0.115 sigma on a
   common 1980-2019 window (the larger 0.426 was two-thirds record-length, not
   event definition). Gate 8 permutation P=0.28 with 0/5 predictors.

## Data
- Observations: CPC daily AO/NAO/PNA/AAO 1950-2026; NCEP stratosphere 1979-2024;
  100 hPa eddy heat flux v'T' 1958-2024 (own OPeNDAP pull); ONI.
- Frozen event catalogue: 43 events / 36 winters, era-independent consensus rule
  (>=2/3 of *covering* reanalyses; the original ">=4 of 6" demanded unanimity
  pre-1979, Fisher P=0.011 for era confounding). Rebuilds byte-for-byte.
- CMIP6 daily ua(10 hPa)+psl from Google Cloud zarr, no auth: 11 models,
  20 members, ~180 GB read, zonal means retained.

## Current results

**A. R diagnostic (main).** R = within-winter effect / between-winter effect at
+15..+29 d. Within-winter compares to the same winter's own baseline (>75 d from
onset); between-winter to day-of-year-matched non-event winters.
- Observations R = 0.824. Simulated pure-common-cause null (1500 reps, AR(1)
  phi=0.943, 43 events/36 winters): mean -0.046, sd 0.320. P(R>=0.824 | f=0) =
  0.0000. Winter-scale confounding rejected, p<0.0007.
- Validated under known truth: bias <= 0.047 across causal fraction f in [0,1];
  f=0 vs f=1 separable at n=43. BUT winter-block bootstrap under-covers
  (0.83-0.90 vs 0.95), so I use the simulated null, and the calibrated CI is
  [0.20, 1.45].
- Replicated: 6 CMIP6 models give R = 0.80-1.00 (CanESM5 0.83 n=846,
  IPSL-CM5A2-INCA 0.89, MPI-ESM-1-2-HAM 1.00, CESM2-FV2 0.83, INM-CM5-0 0.80,
  MPI-ESM1-2-LR 0.88).
- **PROBLEM I FOUND**: R is FLAT across lag in observations — 0.79 at -30..-16
  (pre-onset) vs 0.82 at +15..+29. A precursor cannot be caused by its own event,
  so R is NOT a causal fraction. It measures temporal CONCENTRATION of the
  anomaly around the event. So the valid claim narrows to "the association is
  temporally localised, not winter-scale". A short-timescale shared driver
  produces the same signature.

**B. Granger direction.** Vortex->surface minus surface->vortex, 10 lags,
winter-block bootstrap, conditional on v'T'.
- Observations: asymmetry +0.0014, P(asym<=0)=0.213 — NOT significant.
- Models: CanESM5 +0.0034 (P=0.0000), IPSL +0.0032 (P=0.0013), CESM2-FV2 +0.0049
  (P=0.0025), MPI-ESM-1-2-HAM +0.0015 (P=0.063).
- Observed CI contains every model value => consistent but underpowered at 47
  winters. Effect sizes tiny (0.3-0.5% of variance).

**C. Mediation, wave driving -> vortex -> surface.** Continuous daily, 6958
winter days, 47 winters. Total X->Y = -0.152 (p<0.0001).
- Single-day mediator gives 14.2% mediated => looks like "common cause".
  That is measurement-error attenuation: 10-day mean gives 45.2% [8,88];
  reliability of the 10-day u10 mean as a measure of the latent vortex state is
  0.592 (two-indicator IV using 100 hPa polar-cap height); disattenuated 53.9%
  but CI [-11%, 140%].
- Claim: published weak-stratospheric-influence findings may be mediator
  measurement-error artefacts. Direction of the bias is known and one-sided.

**D. Selection-on-outcome law.** Karpechko et al. (2017) classify "downward
propagating" SSWs using criteria that ARE the surface response (mean surface NAM
over +8..+52 negative; fraction negative > 0.5). Applied to pseudo-onsets with
true effect zero, selection reproduces 30-44% of the reported anomaly for
AO/NAO/PNA. Bias is additive and subtracts cleanly (residuals 0.010-0.024 sigma).
- Derived law: bias(Y) = beta(Y_window on A_window) x bias(selection variable),
  where beta is the ordinary regression SLOPE (not correlation — dropping
  sigma_Y/sigma_A over-predicts 3x). R^2 = 0.999 across 5 outcomes.
- Generalises to ENSO (structurally unrelated): slope +0.976, intercept -0.002,
  R^2 = 0.985 on circulation outcomes. FAILS for persistent non-Gaussian impact
  variables (snow depth, R^2 = 0.004).
- Uncontaminated alternative: classification RATE, 70% real vs 38% pseudo,
  P<0.0001.

**E. Model evaluation side-finding.** Of 11 CMIP6 models, three barely reproduce
the observed SSW surface composite at all (observed between-winter -0.991;
NorESM2-LM +0.072 i.e. wrong sign, MIROC6 -0.138, GFDL-CM4 -0.244).
corr(composite, R) across 8 usable models = -0.634, P=0.091 — not significant.

## Errors I caught and corrected (context for how much to trust the above)
- Bootstrap seeded from hash(str): Python randomises per process, so no p-value
  in 3 scripts was reproducible. Fixed (crc32), now byte-identical across runs.
- Negative control "failure": AAO flagged at q=0.036 pooled across 32 tests;
  pooling a control with positive outcomes raises the BH threshold and makes the
  control easier to flag. On its own 8 bins q=0.133 — passes. Chased seed noise,
  bootstrap calibration and ENSO before finding it was family definition.
- A "20-27% bootstrap under-coverage" that was purely a sample-size confound
  (sqrt(13864/8862)=1.251, observed ratio 0.80=1/1.251). Nearly reported.
- CMIP6 _FillValue (1e20) averaged into zonal means: MIROC6 gave a 10 hPa wind of
  2.6e18 m/s, 26% of cells corrupt. CanESM5 had no fill, which hid it for 10
  members.
- CanESM5 noleap calendar cast to datetime64 collapsed 165 years into 1 "winter"
  (1 SSW/member instead of ~90).

## What I need decided
1. Given R is concentration not causation, what is the strongest *honest*
   headline? Is "stratospheric influence is directional and temporally localised
   but weak, and the observational record alone cannot establish direction" a
   paper, or is it a null dressed up?
2. Is Nature Geoscience realistic for any of A-E, or is that self-deception?
   I currently think no. Am I wrong?
3. Which of A-E is the actual lead result? I keep switching.
4. Is there a decisive test I am missing that works with 47 winters of
   observations + 11 CMIP6 models + no intervention data? SNAPSI (nudged minus
   control) is the obvious answer and is blocked behind a CEDA account.
5. Anything above that you think is wrong or overclaimed.
