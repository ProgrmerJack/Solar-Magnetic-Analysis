# results

One tree for every computed result in this repository. Regenerated indexes come from `python tools/catalog_repo.py`; **the JSON files themselves are written by the analysis scripts and must not be hand-edited.**

| file | what it is |
|---|---|
| `_ALL_RESULTS.json` | all live result payloads in ONE file — read this when you want the numbers |
| `_PROVENANCE.json` | producer + sha256 per result; drives the staleness check in `CATALOG.md` |

## `current/` — the live analysis

53 results from `ssw-design-analysis/`, grouped by the question they answer.

### `current/1_catalogue/` — 1

the frozen event catalogue and its nine sensitivity sets

| result | produced by |
|---|---|
| `catalogue_sets.json` | `build_catalogue.py` |

### `current/2_event_study/` — 10

the core AO/NAO lead-lag estimates and their robustness subsets

| result | produced by |
|---|---|
| `canonical_event_study.json` | `canonical_event_study.py` |
| `clean_subset.json` | `clean_subset.py` |
| `common_period.json` | `common_period.py` |
| `enso_confound.json` | `enso_confound.py` |
| `era_split.json` | `era_split.py` |
| `isolated_events.json` | `isolated_events.py` |
| `isolation_interaction.json` | `isolation_interaction.py` |
| `multi_index_event_study.json` | `multi_index_event_study.py` |
| `negcontrol_stability.json` | `negcontrol_stability.py` |
| `seed_stability.json` | `seed_stability.py` |

### `current/3_calibration/` — 8

estimator behaviour under a known-zero truth (Gate 3)

| result | produced by |
|---|---|
| `aao_null.json` | `aao_null.py` |
| `ao_null_perbin.json` | `ao_null_perbin.py` |
| `boot_calibration.json` | `boot_calibration.py` |
| `calibration_results.json` | `calibrate_seasonality.py` |
| `gate3_clean_null.json` | `gate3_clean_null.py` |
| `spatial_calibration.json` | `calibrate_spatial.py` |
| `station_calibration.json` | `calibrate_station_level.py` |
| `validate_R_diagnostic.json` | `validate_R_diagnostic.py` |

### `current/4_selection_bias/` — 8

selection-on-outcome bias and the beta x dS projection law

| result | produced by |
|---|---|
| `selection_law_enso.json` | `selection_law_enso.py` |
| `selection_on_outcome.json` | `selection_on_outcome.py` |
| `selection_projection_law.json` | `selection_projection_law.py` |
| `selection_snotel.json` | `selection_snotel.py` |
| `stratifier_bias_law.json` | `stratifier_bias_law.py` |
| `stratifier_bias_law_1958.json` | `stratifier_bias_law.py` |
| `stratifier_inference.json` | `stratifier_inference.py` |
| `stratifier_law_cmip6.json` | `stratifier_law_cmip6.py` |

### `current/5_mechanism/` — 10

the R diagnostic, Granger direction, and mediation

| result | produced by |
|---|---|
| `R_gap_sensitivity.json` | `R_gap_sensitivity.py` |
| `R_null_calibration.json` | **orphan** |
| `R_profile.json` | `R_profile.py` |
| `causal_timescale_ratio.json` | `causal_timescale_ratio.py` |
| `continuous_mediation.json` | `continuous_mediation.py` |
| `design_sensitivity.json` | `design_sensitivity.py` |
| `granger_direction.json` | `granger_direction.py` |
| `mediation_disattenuated.json` | `mediation_disattenuated.py` |
| `mediator_specification.json` | `mediator_specification.py` |
| `multimodel_R.json` | `multimodel_R.py` |

### `current/6_predictability/` — 6

out-of-sample skill, forced variance, and the class test

| result | produced by |
|---|---|
| `forced_variance_ceiling.json` | `forced_variance_ceiling.py` |
| `headline_stress_test.json` | `headline_stress_test.py` |
| `is_downward_propagation_a_class.json` | `is_downward_propagation_a_class.py` |
| `predictability_ceiling.json` | `predictability_ceiling.py` |
| `predictability_wave_driving.json` | `predictability_wave_driving.py` |
| `within_model_check.json` | `within_model_check.py` |

### `current/7_ensemble/` — 6

CMIP6 ensembles, marginal events, and the era gate

| result | produced by |
|---|---|
| `consensus_monotonicity.json` | `consensus_monotonicity.py` |
| `ensemble_precursor.json` | **orphan** |
| `ensemble_precursor_CanESM5.json` | `ensemble_precursor.py` |
| `ensemble_precursor_MIROC6.json` | `ensemble_precursor.py` |
| `era_dependence_gate.json` | `era_dependence_gate.py` |
| `marginal_events.json` | `marginal_events.py` |

### `current/8_experiment/` — 1

SNAPSI nudged-minus-control, the one causal estimate

| result | produced by |
|---|---|
| `snapsi_causal_effect.json` | `snapsi_causal_effect.py` |

### `current/9_literature/` — 3

recomputation of published criteria on our catalogue

| result | produced by |
|---|---|
| `era5_recompute_and_two_thirds.json` | `era5_recompute_and_two_thirds.py` |
| `loeffel2025_missing_control.json` | `loeffel2025_missing_control.py` |
| `recompute_published_criterion.json` | `recompute_published_criterion.py` |

## `superseded/` — the retired project

Kept for provenance. Nothing here is current; the avalanche and solar-magnetic theses were both abandoned. Grouped by family:

| family | files |
|---|---|
| `audit_and_meta/` | 6 |
| `event_catalogue/` | 3 |
| `fresh_analysis/` | 6 |
| `limitations/` | 7 |
| `mechanism/` | 7 |
| `numbered/` | 23 |
| `review_rounds/` | 102 |
| `revision/` | 4 |
| `sensitivity/` | 7 |
| `snowpack/` | 9 |
| `stratosphere/` | 6 |
| `synthesis/` | 13 |

## Orphans

Results no script writes. Both are real and were verified by grep:

- `R_null_calibration.json` — nothing writes **or** reads it, yet its calibrated CI [0.20, 1.45] is quoted in `REVIEW_BRIEF.md`. Result C is corrected anyway (R is concentration, not a causal fraction).
- `ensemble_precursor.json` — a leftover from before `ensemble_precursor.py` switched to per-model output (`ensemble_precursor_CanESM5.json` etc.).
