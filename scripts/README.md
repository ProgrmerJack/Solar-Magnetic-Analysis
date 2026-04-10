# Scripts Directory Structure

This directory contains all code used to produce the results in the manuscript. Scripts are organized into thematic subdirectories.

## Directory Layout

```
scripts/
├── download/           # Core data acquisition pipeline (numbered 01–53)
├── download_extra/     # Supplementary data downloads (ALBINA, EAWS, ERA5, Norway, etc.)
├── process/            # Raw data processing pipeline (numbered 00–11)
├── process_extra/      # Additional data processing (ERA5 extension, QC fixes, sintering model)
├── analysis/           # Core analysis pipeline (numbered 00–28)
├── analysis_extra/     # Standalone analyses (ERA5 mechanism, mediation, Norway, US, etc.)
├── review_rounds/      # Iterative analysis for peer review rounds (r11–r39)
├── utilities/          # QC checks, data inspection, verification, tests
└── README.md           # This file
```

## Subdirectory Descriptions

### `download/` — Core Data Acquisition (53 scripts)
Numbered pipeline for downloading all primary datasets: solar indices, OMNI, GOES X-ray, flare catalogs, SSW catalog, CAIC avalanche, SNOTEL, DSCOVR/ACE, POES/MEPED, ERA5, NCEP, Norwegian NVE, Kyoto geomagnetic indices, PSP, MLS, SuperMAG, and more. Run `download_all.py` or execute scripts in order.

### `download_extra/` — Supplementary Downloads (24 scripts)
Additional data downloads added during the research process: ALBINA bulletins, EAWS multi-country danger levels, EnviDat Swiss danger/stability data, ERA5 extended variables, French BRA data, Norwegian targeted data, US danger ratings, Austrian LAWIS, climate indices.

### `process/` — Core Data Processing (11 scripts)
Numbered pipeline converting raw data to analysis-ready formats: geomagnetic indices, solar catalogs, OMNI, GOES XRS/particle, POES, MLS, ERA5, MODIS, cryosphere, PSP/ACE. Run `00_run_pipeline.py` or execute in order.

### `process_extra/` — Supplementary Processing (7 scripts)
Additional processing: ERA5 extended field processing, NetCDF fixes, QC failure remediation, SNOWPACK sintering model implementation.

### `analysis/` — Core Analysis Pipeline (28 scripts)
Numbered from `00_preregistration.py` through `28_mechanism_upgrade.py`. Includes panel construction, event definition, primary endpoint analysis, figures, synthesis, and mechanism analysis. The key scripts for reproducing paper figures and tables are:
- `24_fresh_analysis.py` — Primary Swiss analysis + figures
- `25_revision_analysis.py` — Revised analysis + main paper figures
- `26_mechanism_analysis.py` — ERA5 mechanism composites + figures
- `27_planetary_wave_analysis.py` — Planetary wave timing
- `28_mechanism_upgrade.py` — Full mechanism chain

### `analysis_extra/` — Standalone Analyses (20 scripts)
Individual analyses addressing specific research questions: ERA5 mechanism pathways, mediation analysis, Norwegian phase analysis, US danger analysis, sintering timescale, SSW subtype stratification, threshold vs. continuous models, autocorrelation correction.

### `review_rounds/` — Peer Review Round Analyses (44 scripts)
Scripts created during iterative peer review, prefixed by round number (r11–r39). Each addresses specific reviewer concerns:
- **r20–r21**: Definitive multi-country analysis, paper-verified statistics
- **r22–r23**: SSW-type stratification, confounder independence
- **r25–r28**: Mechanism deepening, mediation, wave activity
- **r29–r32**: EAWS multi-country, EnviDat, SNOWPACK stability, trigger suppression
- **r34–r39**: Multi-country replication, Bayesian analysis, specification curve, ALBINA/Norway expansion

### `utilities/` — Quality Control & Verification (24 scripts)
Data QC checks, file inspections, result verification, API tests, and paper-number auditing.

## Reproduction

To reproduce the full analysis from raw data:
1. Run `download/download_all.py` (or scripts 01–53 in order)
2. Run `process/00_run_pipeline.py` (or scripts 01–11 in order)
3. Run `analysis/01_build_daily_panel.py` through `analysis/28_mechanism_upgrade.py`
4. Run review-round scripts as needed for specific claims

See `REVIEWER_INDEX.md` for a claim-by-claim mapping to generating code.
