| id | test | role | registration | commit | statistic | p | p_Holm_primary |
|---|---|---|---|---|---|---|---|
| M-var | variance ratio nudged/control, re-centred (SNAPSI NH) | secondary | in script, run same day | 0b6886d | 0.9554 [0.8732, 1.0507] |  |  |
| M-KS | pure translation, Kolmogorov-Smirnov (SNAPSI NH) | secondary | in script, run same day | 0b6886d | KS | 0.3714 |  |
| N | week-2 100 hPa coupling between members, nudged minus control (Fisher z) | secondary | in script with a design gate | 07324a4 | dz 0.023 | 0.5037 |  |
| J-P1 | within-model out-of-sample R2 added by the SSW, before onset (CMIP6) | secondary | in script, run same day | 0b6886d | -0.0057 | 0.6444 |  |
| J-P3 | within-model R2 added by the SSW, with post-onset stratosphere | secondary | in script, run same day | 0b6886d | 0.0151 | 0.2517 |  |
| O | archetype pair: observed 2018-2019 difference within model pairings | secondary | in script, run same day | 79007c2 | median percentile 0.331; outside 95%: 0 of 9 |  |  |
| L | paired DW-NDW contrast, nudged minus control (SNAPSI) | primary | in script, run same day | 0eda6d9 | -0.0283 [-0.1219, 0.0544] |  |  |
| L-d | nudged contrast vs empirical shifted-control null | diagnostic | post hoc | this revision | 0.0401 [-0.0263, 0.0974] |  |  |
| R | regional residual beyond the NAM, N Eurasia (SNAPSI) | primary | in script, run same day | 103a2c1 | -0.099 [-0.36, 0.0592] |  |  |
| R-obs | regional residual beyond the NAM, N Eurasia (ERA5) | secondary | in script, run same day | 103a2c1 | -0.3255 K | 0.3814 |  |
| Q | SSW-specific CRPS skill of event-aware forecasts (CMIP6, P2) | primary | in script, run same day | 0eda6d9 | -0.0298 | 1.0 | 1.0 |
| P-D | event discrimination, ECMWF (conditional null) | discovery | in script, run same day | d9d8243 | r 0.7655 | 0.1147 |  |
| P'-D | event discrimination, nine-system mean | secondary | committed before run | 0eda6d9 | r 0.2262 | 0.9129 |  |
| P'-H1 | outcomes low in ensembles, nine-system mean (polar cap) | primary | committed before run | 0eda6d9 | rank 0.3979 vs 0.5437 | 0.005 | 0.03 |
| P'-H2 | shift-size correction improves CRPS | primary | committed before run | 0eda6d9 | gain 120.894 | 0.1604 | 0.4344 |
| R-op | N-Eurasian temperature low in ensembles (nine-system mean) | primary | committed before data | 103a2c1 | rank 0.3896 vs 0.5273 | 0.0135 | 0.0675 |
| B1 | rank deficit vs all-winter baseline (polar cap) | sensitivity | post hoc | this revision | 0.3979 vs 0.5288 | 0.0124 |  |
| B2 | rank deficit beyond vortex-state regression (polar cap) | sensitivity | post hoc | this revision | residual -0.1016 | 0.0284 |  |
| R3 | rank deficit, mixed model with date random effects (polar cap) | sensitivity | post hoc | this revision | -0.1446 [-0.2538, -0.0355] | 0.00938 |  |
| T1 | rank deficit, starts after onset | primary (revision) | committed before data | 57019f6 | 0.4966 | 0.1448 | 0.4344 |
| T2 | rank deficit, pre-onset starts that caught the SSW | primary (revision) | committed before data; matched null adopted after a review had emulated the polar-cap test (it raised p) | 57019f6, a201458 | 0.4113 | 0.0136 | 0.0675 |
| T2-miss | rank deficit, pre-onset starts that missed the SSW | secondary (revision) | as T2 | 57019f6, a201458 | 0.4137 | 0.029 |  |
| T2-d | hit minus miss rank (polar cap) | secondary (revision) | as T2 | 57019f6, a201458 | 0.0028 [-0.0473, 0.0529] |  |  |
| G1 | split minus displaced surface response (ERA5) | secondary (revision) | committed before data | 57019f6 | -0.245 | 0.2768 |  |
| Q-P3 | CRPS skill with post-onset stratosphere: after SSWs vs ordinary days | sensitivity | post hoc | this revision | 0.1717 vs 0.2068 | 1.0 |  |
| T-a1 | criterion sweep, 54 surface versions: max abs z of SSW rate vs shifted null | robustness (revision 2) | committed before run | 4c8b3f2 | max abs z 1.154 | 0.965 |  |
| T-a2 | criterion sweep, 54 surface versions: systematic excess (mean z) | robustness (revision 2) | committed before run | 4c8b3f2 | mean z -0.036 | 0.4975 |  |
| T-a1-cal | criterion sweep, surface versions: max abs z, calibrated reference | robustness (revision 2) | post hoc calibration of a registered test | this revision | max abs z 1.154 | 0.845 |  |
| T-a2-cal | criterion sweep, surface versions: mean z, calibrated reference (two-sided) | robustness (revision 2) | post hoc calibration of a registered test | this revision | mean z -0.036 | 0.96 |  |
| T-b | observed statistics within 43-event CMIP6 draws (min-p combination) | robustness (revision 2) | committed before run | 4c8b3f2 | min p 0.222 | 0.623 |  |
| T-c | identification audit: contrast / ceiling at SNAPSI-bounded damping (obs; CMIP6) | robustness (revision 2) | committed before run | 4c8b3f2 | 1.381; 1.671 |  |  |
| T-d | two regimes beat continuous model after SSWs vs ordinary days (CRPS) | robustness (revision 2) | committed before run | 4c8b3f2 | G 0.00124 vs -0.0023 | 0.01 |  |
| T-d-PH | two regimes vs one skewed population, after SSWs vs ordinary days | sensitivity | post hoc | this revision | G 0.00109 | 0.02 |  |
| T-d-PH2 | two regimes vs continuous, calibrated on no-regime synthetic data | sensitivity | post hoc | this revision | G 0.00124 | 0.02 |  |
| T-d-PH4 | two-regime model skill over the shift: after SSWs vs ordinary days | sensitivity | post hoc | this revision | 0.0112 vs 0.0309 | 0.015 |  |
| T-d-PH3 | two-regime components after SSWs: mean separation (ordinary days) | descriptive | post hoc | this revision | 0.4937 (0.7142) |  |  |
| HO1 | held-out 2023-24 SSWs: polar-cap rank deficit | replication (revision 2) | committed before data | 4c8b3f2 | pending |  |  |
| HO2 | held-out 2023-24 SSWs: N-Eurasian temperature rank deficit | replication (revision 2) | committed before data | 4c8b3f2 | pending |  |  |
| L1 | 100 hPa polar-cap height rank after SSWs | diagnostic (revision 2) | committed before data | 4c8b3f2 | pending |  |  |
| L2 | surface rank conditional on forecast 100 hPa anomaly | diagnostic (revision 2) | committed before data | 4c8b3f2 | pending |  |  |
