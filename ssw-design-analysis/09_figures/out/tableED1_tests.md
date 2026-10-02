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
| HO1 | held-out 2023-24 SSWs: polar-cap rank deficit | replication (revision 2) | committed before data | 4c8b3f2 | 0.7167 vs 0.532 | 0.9308 |  |
| HO2 | held-out 2023-24 SSWs: N-Eurasian temperature rank deficit | replication (revision 2) | committed before data | 4c8b3f2 | 0.5704 vs 0.5144 | 0.6383 |  |
| HD-2a | ECMWF CY49R1 (model year 2025) on the 1998-2021 SSWs: rank deficit | diagnostic | post hoc (exploratory) | this revision | 0.4162 vs 0.5259 | 0.0813 |  |
| HD-2b | ECMWF CY49R1 minus CY47R3 rank, same 10 events | diagnostic | post hoc (exploratory) | this revision | 0.019 [-0.055, 0.0762] |  |  |
| HD-2c | CMA model year 2025 on the 1998-2021 SSWs: rank deficit | diagnostic | post hoc (exploratory) | this revision | 0.24 vs 0.5376 | 0.002 |  |
| HD-6 | held-out signs: likelihood ratio, forecasts as issued vs 1998-2021 bias | diagnostic | post hoc (exploratory) | this revision | 3.83 |  |  |
| L1 | 100 hPa polar-cap height rank after SSWs | diagnostic (revision 2) | committed before data | 4c8b3f2 | 0.4458 vs 0.5945 | 0.005 |  |
| L2 | surface rank conditional on forecast 100 hPa anomaly | diagnostic (revision 2) | committed before data | 4c8b3f2 | 0.4473 vs 0.4765 | 0.2999 |  |
| HO2026 | held-out 4 March 2026 SSW, real-time forecasts: polar-cap rank deficit | replication (revision 3) | committed before data | a511fbc | 0.6664 vs 0.6057 | 0.598 |  |
| HO2026-T | held-out 4 March 2026 SSW: N-Eurasian temperature rank deficit | replication (revision 3) | committed before data | a511fbc | 0.9323 vs 0.5644 | 0.9815 |  |
| HO-4 | four held-out SSWs (2023-2026) pooled: polar-cap rank deficit | replication (revision 3) | committed before data | a511fbc | 0.7041 vs 0.5503 | 0.9268 |  |
| U-dose | AO days 8-52 per 10 m/s of minimum 10 hPa wind, 114 observed episodes | revision-3 primary | committed before outcomes | 2a0c40b | 0.3582 [0.213, 0.5104] |  |  |
| U-jump | AO step at wind reversal, 114 observed episodes (Holm over 4 outcomes: 1.0) | revision-3 primary | committed before outcomes | 2a0c40b | -0.2755 [-0.9867, 0.465] | 0.4598 |  |
| U-cmip6 | CMIP6 step at wind reversal (sigma), 5,726 episodes | revision-3 secondary | committed before outcomes | 2a0c40b | 0.0512 [-0.02, 0.1278] | 0.157 |  |
| V-0.10 | shift rule (SNAPSI) vs observed cold fortnights below 10th pct., 39 SSWs: binomial under rule | revision-3 primary | committed before run | bf453bd | 10/39 vs p_rule 0.3243 | 0.3986 |  |
| V-clim | observed cold fortnights below 10th pct. vs climatology | revision-3 primary | committed before run | bf453bd | 10/39 vs 0.10 | 0.0042 |  |
| V-R1 | shift rule, strict independence (37 events, no SNAPSI events), q 0.10: binomial under rule | revision-3 secondary | committed before outcomes | add7ba9 | 9/37 | 0.3801 |  |
| V-R2 | shift rule, 2.5th percentile, N Eurasia: binomial under rule | revision-3 secondary | committed before outcomes | add7ba9 | 4/39 vs 0.1786 | 0.295 |  |
| V-R3 | shift rule, 2.5th percentile, high-latitude Europe: binomial under rule (rule rejected) | revision-3 secondary | committed before outcomes | add7ba9 | 0/39 vs 0.201 | 0.0002 |  |
| V-R4 | shift rule, days 15-28, N Eurasia, q 0.10: binomial under rule | revision-3 secondary | committed before outcomes | add7ba9 | 6/39 vs 0.2724 | 0.107 |  |
| V-R5 | second rule (CMIP6 x ERA5 event-free), N Eurasia q 0.10: binomial under rule | revision-3 secondary | committed before outcomes | add7ba9 | 10/39 vs 0.1843 | 0.2985 |  |
| V-R6 | shift rule OUT OF SAMPLE: 12 ERA5 SSWs 1941-58, N Eurasia q 0.10 (rule rejected) | revision-3 primary | committed before outcomes | add7ba9 | 0/12 vs 0.3243 | 0.0122 |  |
| U-ERA5 | AO step at reversal, ERA5 winds 1940-2025 (Holm over 4 outcomes) | revision-3 secondary | committed before data | add7ba9 | -0.5364 [-1.2269, 0.093] | 0.4008 |  |
| N-shift | within-start ensemble comparison: shift, reversing minus non-reversing members (8-25 d, 40 starts) | revision-3 primary | committed before computation | be804a7 | -0.0535 [-0.2864, 0.1671] |  |  |
| N-var | within-start ensemble comparison: variance ratio (no widening) | revision-3 primary | committed before computation | be804a7 | 0.974 [0.7317, 1.3602] |  |  |
| N-contrast | within-start ensemble comparison: class contrast difference | revision-3 primary | committed before computation | be804a7 | -0.0036 [-0.3335, 0.2963] |  |  |
| N-step | within-ensemble step at reversal, initial state fixed (11,481 members) | revision-3 primary | committed before computation | be804a7 | 0.0439 [-0.1631, 0.2994] | 0.652 |  |
| RQ-A | SNAPSI residual quantiles, polar-cap NAM, pooled: inside declared tolerance (+-0.25 sigma, +-0.05)? | revision-3 primary | committed before run | 113c7d3 | approximate translation |  |  |
| RQ-T | SNAPSI residual quantiles, N-Eurasian temperature, pooled | revision-3 secondary | committed before run | 113c7d3 | unresolved |  |  |
| TCV | shifted null, whole-winter cross-validation: observed vs expected downward count (retrospective) | revision-3 primary | committed before run | bbfeb9c | 27 vs 29.21 | 0.5156 |  |
| SF-F0 | state forecast vs climatology, 57 ERA5 SSWs, CRPS skill (retrospective CV) | revision-3 primary | committed before run | 4f5c53a | 0.0942 [-0.0328, 0.2092] |  |  |
| SF-F1 | state forecast vs fixed SNAPSI rule, CRPS skill | revision-3 primary | committed before run | 4f5c53a | 0.103 [-0.0157, 0.2117] |  |  |
| SF-SSW | SSW indicator beyond the continuous state, CRPS skill | revision-3 primary | committed before run | 4f5c53a | -0.0081 [-0.0184, 0.0028] |  |  |
| SF-ops | calibrated operational ensembles vs state forecast, 17 events, CRPS skill | revision-3 secondary | committed before run | 4f5c53a | state vs calibrated -0.4601 [-0.836, -0.0756] |  |  |
| ICON | ICON event ensembles (18, each started at onset or the day before): between-event variance of mean responses, days 8-25 (bounds) | revision-3 secondary | committed before values read | 6ad6e8a | 0.4145-0.4277 |  |  |
| RC | N-Eurasian DW-NDW temperature contrast (days 8-24) vs matched shifted null, winter CV (power 0.30) | revision-4 primary | committed before run | 8c324ed | obs -2.2049 K vs null -2.2524 K | 0.9389 |  |
| GEN-1 | back-shifted SSW-imposed contrast vs control, the 11 pairs with < 3 NDW members (tolerance 0.25 sigma); inconclusive | revision-5 primary | committed before run | 892377c | 0.1734 [0.0404, 0.3304]; net of the spread change 0.0299 [-0.0268, 0.0635] (post hoc); forced variance of half the control variance added: -0.1485 [-0.2992, 0.0354] (not detected) |  |  |
| GEN-2 | dose-response: back-shift difference on imposed shift, 36 pairs (slope per sigma) | revision-5 secondary | committed before run | 892377c | -0.0311 [-0.2048, 0.0399] |  |  |
| GEN-3 | SH minor warming: nudged contrast vs location-scale null (and shift-only null) | revision-5 primary | committed before run | 892377c | -0.0365 [-0.1376, 0.0794] (shift only -0.136 [-0.2216, -0.033]) |  |  |
| GEN-4 | SNAPSI N-Eurasian DW-NDW temperature contrast vs matched shifted null (sigma; tolerance 0.25) | revision-5 primary | committed before run | 892377c | -0.0861 [-0.2398, 0.057]; planted 0.5 sigma -0.5861 |  |  |
| PRO | prospective: shifted null on every NH SSW with onset Nov 2026 - Mar 2036 (CRPS, days 8-52 NAM) | prospective primary | predictions committed before any onset | d2b54a8 (amended before any onset: 4cfba00, b2416e0) | 0 events scored |  |  |
| S-dose | two regimes vs continuous with the realised 100 hPa dose (G, p vs ordinary days); inconclusive | revision-3 primary | committed before run | 7d8ed6f | G -0.0146 vs -0.01094; power vs planted regimes 0.08 (post hoc): underpowered | 0.92 |  |
