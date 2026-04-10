# Planetary Wave Forcing Suppresses Natural Dry Slab Avalanche Activity via Stratospheric Sudden Warmings

[![DOI](https://img.shields.io/badge/DOI-pending-blue)](https://doi.org/)
[![Data: Zenodo](https://img.shields.io/badge/Data-Zenodo-blue)](https://zenodo.org/)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Target: Nature Geoscience](https://img.shields.io/badge/Target-Nature%20Geoscience-green.svg)](https://www.nature.com/ngeo/)

> **Multi-Country Evidence, Process-Model Validation, and Count–Rating Dissociation**

---

## Abstract

Stratospheric sudden warming (SSW) events produce persistent surface weather anomalies, yet their consequences for geophysical hazards have never been systematically quantified. Here we identify a previously unknown connection: SSW episodes drive a large, reproducible reduction in natural dry slab avalanche activity across four countries and two continents, mediated by a **threshold amplification** mechanism in which modest stratospheric-origin weather shifts produce disproportionately large hazard responses by switching binary trigger conditions across mountain networks. In Switzerland (21 winters, 16 SSW events), 14 of 16 events show reduced counts (geometric mean rate ratio RR = 0.32; 95% CI [0.20, 0.54]; *d* = −1.06). Norwegian avalanche danger levels decrease concordantly (*d* = −0.67; *P* < 10⁻⁶), and all four Utah SSW events show decreased counts (RR = 0.34). A specification curve of 180 analytical variants confirms robustness: 100% show RR < 1 (permutation *P* < 0.001). Five-country EAWS data reveal a predictable geographic gradient in hazard response (*r* = 0.69, *P* = 0.0004), confirmed out-of-sample by independent ALBINA bulletin data. A count–rating dissociation and a differential trigger response reveal a "loaded gun" mechanism: SSW-associated cold regimes suppress surface-energy triggers (−57% warming, −29% rain) while preserving snowpack instability.

---

## Scientific Hypothesis

### SOC Universality and Stratosphere–Surface Coupling

Snow avalanches are a canonical example of **Self-Organised Criticality (SOC)**: snowpack accumulates stress incrementally until a threshold is crossed, producing a power-law distributed release event. We hypothesise that stratospheric variability—specifically, the enhanced planetary wave forcing that triggers SSW events—modulates the *trigger conditions* for these threshold-governed hazards without altering the underlying criticality. Because natural dry slab release requires *any* of *k* independent trigger channels to exceed a firing threshold, a modest weather shift that simultaneously suppresses all channels produces a **super-linear** (threshold-amplified) reduction in release probability.

This framework predicts that SOC-governed geophysical hazards (avalanches, landslides, debris flows, flooding) should all exhibit amplified sensitivity to stratospheric-origin weather perturbations—a testable prediction that generalises beyond avalanche science.

### Mechanistic Chain

```
Planetary wave forcing (tropospheric)
        │
        ├──→ Stratospheric Sudden Warming (SSW)
        │         │
        │         ├──→ EPP → NOₓ / O₃ perturbation (chemical fingerprint)
        │         │
        │         └──→ Downward coupling (1–4 weeks)
        │                   │
        │                   └──→ Z500 depression / blocking
        │
        └──→ Direct tropospheric reorganisation (common cause)
                      │
                      └──→ Cold-dry weather regime shift
                                │
                                ├──→ Surface warming events −57%
                                ├──→ Rain-on-snow events −29%
                                ├──→ Shortwave radiation −17.5%
                                │
                                └──→ Trigger suppression
                                          │
                                          ├──→ Natural dry slab avalanches ↓ 68%
                                          ├──→ Human-triggered avalanches ↑ (loaded gun)
                                          └──→ Danger ratings ↑ (count–rating dissociation)
```

The temporal profile—suppression emerging **before** SSW onset—identifies planetary wave forcing as the **common cause** rather than purely top-down stratospheric forcing.

---

## Repository Structure

```
Solar-Magnetic-Analysis/
│
├── paper/                      # Manuscript and supplementary materials
│   ├── main.tex                # Main manuscript (Nature Geoscience Article format)
│   ├── supplementary_information.tex
│   └── main.pdf                # Compiled manuscript
│
├── scripts/                    # All analysis code (see below)
│   ├── download/               # Core data acquisition pipeline (53 scripts)
│   ├── download_extra/         # Supplementary data downloads (24 scripts)
│   ├── process/                # Raw → processed data pipeline (11 scripts)
│   ├── process_extra/          # Additional processing & sintering model (7 scripts)
│   ├── analysis/               # Core analysis pipeline (28 scripts)
│   ├── analysis_extra/         # Standalone mechanism & mediation analyses (20 scripts)
│   ├── review_rounds/          # Peer-review-driven analyses, r11–r39 (44 scripts)
│   └── utilities/              # QC checks, data inspection, verification (24 scripts)
│
├── data/                       # Data directory (not tracked in Git)
│   ├── raw/                    # Raw downloaded data
│   ├── processed/              # Analysis-ready parquet files
│   ├── results/                # JSON/CSV result files
│   ├── figures/                # Generated figures (PDF)
│   ├── solar/                  # Solar index data
│   ├── geomagnetic/            # Geomagnetic index data
│   ├── atmospheric/            # ERA5, NCEP, MLS reanalysis
│   └── cryosphere/             # Avalanche, snow, and danger-level data
│
├── notebooks/                  # Exploratory notebooks
│   └── nature_geoscience_paper.md
│
├── logs/                       # Download and processing logs
│
├── REVIEWER_INDEX.md           # Claim-by-claim traceability (every number → script)
├── ZENODO_DEPOSITS.json        # Zenodo deposit metadata (raw & processed data)
├── zenodo_upload.py            # Zenodo data upload script
├── requirements.txt            # Python dependencies
└── README.md                   # This file
```

### Key Directories

| Directory | Contents | Scripts |
|-----------|----------|---------|
| `scripts/download/` | Core data acquisition: solar indices, OMNI, GOES X-ray, flare catalogs, SSW catalog, SLF avalanche, SNOTEL, ERA5, NCEP, NVE, Kyoto indices, and more | 53 numbered scripts + `download_all.py` |
| `scripts/download_extra/` | Supplementary downloads: ALBINA, EAWS, EnviDat, ERA5 extended, French BRA, Norwegian targeted, US danger ratings, LAWIS | 24 scripts |
| `scripts/process/` | Raw → analysis-ready conversion: geomagnetic, solar, OMNI, GOES, POES, MLS, ERA5, MODIS, cryosphere, PSP/ACE | 11 numbered scripts + `00_run_pipeline.py` |
| `scripts/process_extra/` | Extended ERA5 processing, NetCDF fixes, SNOWPACK sintering model | 7 scripts |
| `scripts/analysis/` | Core analysis: panel construction, SSW event definition, primary endpoint, figures, mechanism analysis | 28 numbered scripts (00–28) |
| `scripts/analysis_extra/` | ERA5 mechanism pathways, mediation analysis, Norwegian phase analysis, US danger, sintering timescale, SSW subtypes | 20 standalone scripts |
| `scripts/review_rounds/` | Iterative peer-review analyses: multi-country replication, confounder independence, specification curve, Bayesian evidence, EAWS gradient, trigger suppression | 44 scripts (r11–r39) |
| `scripts/utilities/` | Data QC, file inspection, result verification, API tests | 24 scripts |

---

## Key Scripts

| Research Question | Primary Script |
|---|---|
| Swiss primary analysis (RR, CI, *P*-values) | `scripts/review_rounds/r21_paper_analysis.py` |
| Phase-resolved temporal structure | `scripts/analysis/24_fresh_analysis.py` |
| Norwegian danger-level analysis | `scripts/review_rounds/r38_norway_analysis.py` |
| Utah replication | `scripts/review_rounds/r20_definitive_analysis.py` |
| EAWS European geographic gradient | `scripts/review_rounds/r37_eaws_gradient.py` |
| ALBINA Austrian/Italian gradient | `scripts/review_rounds/r31c_albina_gradient.py` |
| Austrian LAWIS incident analysis | `scripts/review_rounds/r29_grand_analysis.py` |
| SNOWPACK stability indices | `scripts/review_rounds/r31_snowpack_stability.py` |
| Trigger suppression mechanism | `scripts/review_rounds/r31b_trigger_suppression.py` |
| Weather regime mediation (Blinder–Oaxaca) | `scripts/review_rounds/r37_mediation.py` |
| ERA5 surface meteorology composites | `scripts/analysis/26_mechanism_analysis.py` |
| Specification curve (180 variants) | `scripts/review_rounds/r37_spec_curve.py` |
| Bayesian evidence | `scripts/review_rounds/r34_wave_flux_bayesian.py` |
| Confounder independence (QBO, ENSO, PDO, NAO, F10.7) | `scripts/review_rounds/r23_reviewer_upgrades.py` |
| Figures (main paper) | `scripts/analysis/25_revision_analysis.py`, `26_mechanism_analysis.py` |

---

## Data Availability

### Zenodo Deposits

Raw and processed datasets are archived on Zenodo. See [`ZENODO_DEPOSITS.json`](ZENODO_DEPOSITS.json) for deposit metadata.

| Deposit | Zenodo ID | Status |
|---------|-----------|--------|
| Raw data | [19493865](https://zenodo.org/deposit/19493865) | Draft |
| Processed data | [19495580](https://zenodo.org/deposit/19495580) | Draft |

### Primary Data Sources

| Dataset | Source | Access |
|---------|--------|--------|
| Swiss WSL/SLF avalanche counts | WSL/SLF | Upon request |
| Swiss EnviDat danger levels, Rutschblock, SNOWPACK | [EnviDat](https://www.envidat.ch) | Public |
| Norwegian NVE danger forecasts | [NVE API](https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/) | Public |
| Utah Avalanche Center | [UAC](https://utahavalanchecenter.org) | Public |
| EAWS danger levels | [avalanches.org](https://avalanches.org) | Public |
| ALBINA danger levels | [avalanche.report](https://avalanche.report) | Public |
| Austrian LAWIS incidents | [LAWIS](https://lawis.at) | Public |
| ERA5 reanalysis | [Copernicus CDS](https://cds.climate.copernicus.eu) | Registration required |
| Aura/MLS Level 2 | [NASA Earthdata](https://earthdata.nasa.gov) | Registration required |
| SNOTEL | [NRCS](https://www.nrcs.usda.gov/wps/portal/wcc/home/) | Public |
| IMS Snow Cover | [NSIDC](https://nsidc.org) | Public |
| Butler et al. SSW catalog | Butler et al. (2015) | Published |

---

## Reproducibility

### Installation

```bash
# Clone the repository
git clone https://github.com/ProgrmerJack/Solar-Magnetic-Analysis.git
cd Solar-Magnetic-Analysis

# Create a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate   # Linux/macOS
.venv\Scripts\activate      # Windows

# Install dependencies
pip install -r requirements.txt
```

**Note:** Some data sources require credentials:
- **ERA5**: Register at [Copernicus CDS](https://cds.climate.copernicus.eu) and configure `~/.cdsapirc`
- **NASA Earthdata**: Register at [Earthdata](https://earthdata.nasa.gov) and configure `earthaccess`
- **HAPI**: Configure `hapiclient` credentials for CDAWEB access

### Full Reproduction Pipeline

```bash
# 1. Download all raw data
python scripts/download/download_all.py

# 2. Process raw data into analysis-ready formats
python scripts/process/00_run_pipeline.py

# 3. Build the master analysis panel
python scripts/analysis/01_build_daily_panel.py

# 4. Run the core analysis pipeline (scripts 02–28)
python scripts/analysis/02_define_events.py
# ... through ...
python scripts/analysis/28_mechanism_upgrade.py

# 5. Run review-round analyses for specific claims
python scripts/review_rounds/r21_paper_analysis.py   # Swiss primary
python scripts/review_rounds/r38_norway_analysis.py   # Norway
# ... etc.
```

### Tracing Claims to Code

**Every quantitative claim** in the manuscript is mapped to its generating script, input data, and output file in [`REVIEWER_INDEX.md`](REVIEWER_INDEX.md). To verify any number in the paper:

1. Find the claim in `REVIEWER_INDEX.md` (organised by manuscript section and line number)
2. Follow the link to the generating script
3. Run the script to reproduce the result
4. Inspect the output JSON/CSV file

Example: the primary Swiss rate ratio (RR = 0.32) is produced by:
- **Script:** `scripts/review_rounds/r21_paper_analysis.py`
- **Output:** `data/results/r21_paper_analysis.json` → `swiss_core.geo_rr`
- **Input:** `data/processed/analysis_panel_v2.parquet`

### Data Pipeline Overview

```
Raw downloads           →  Processed parquets       →  Analysis results    →  Paper claims
scripts/download/          data/processed/              data/results/          paper/main.tex
scripts/download_extra/    data/processed/              data/results/
```

| Key Intermediate File | Description | Generated By |
|---|---|---|
| `data/processed/analysis_panel_v2.parquet` | Master Swiss daily panel (7,518 rows × 110 cols) | `scripts/analysis/01_build_daily_panel.py` |
| `data/processed/era5_swiss_alps_extended.parquet` | ERA5 daily weather (1998–2019) | `scripts/process_extra/process_extended_era5.py` |
| `data/results/ssw_event_catalog.csv` | 16 SSW events with dates and types | `scripts/analysis/02_define_events.py` |

---

## Python Dependencies

Core analysis requires Python ≥ 3.10 with the following packages (see [`requirements.txt`](requirements.txt)):

| Category | Packages |
|----------|----------|
| **Data I/O** | `pandas`, `numpy`, `netCDF4`, `cdflib`, `h5py`, `xarray` |
| **Analysis** | `scipy`, `powerlaw`, `astropy` |
| **Visualisation** | `matplotlib` |
| **Web/API** | `requests`, `beautifulsoup4`, `lxml`, `tqdm` |
| **Registration-required** | `cdsapi` (ERA5), `earthaccess` (NASA), `hapiclient` (CDAWEB) |

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{Ashuraliyev2025_SSW_Avalanche,
  title   = {Planetary Wave Forcing Suppresses Natural Dry Slab Avalanche
             Activity via Stratospheric Sudden Warmings: Multi-Country
             Evidence, Process-Model Validation, and Count--Rating Dissociation},
  author  = {Ashuraliyev, Abduxoliq},
  year    = {2026},
  journal = {Nature Geoscience},
  note    = {Under review},
  url     = {https://github.com/ProgrmerJack/Solar-Magnetic-Analysis}
}
```

**Author:** Abduxoliq Ashuraliyev  
**Affiliation:** Independent Researcher, Tashkent, Uzbekistan  
**Contact:** Jack00040008@outlook.com  
**ORCID:** [0009-0003-5482-5526](https://orcid.org/0009-0003-5482-5526)

---

## License

This project is released under the [MIT License](LICENSE).

*© 2026 Abduxoliq Ashuraliyev. All rights reserved.*
