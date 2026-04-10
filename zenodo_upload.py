#!/usr/bin/env python3
"""
Zenodo Data Upload Script — Solar-Magnetic Avalanche Analysis
Uploads ALL data (~55 GB) across two deposits:
  1. Raw/Source Data   (~24.6 GB) — unprocessed files from original sources
  2. Processed Data    (~33.2 GB) — analysis-ready datasets, results, figures

Strategy:
  - Small-file directories → zipped by category
  - Large individual files (POES, ERA5, etc.) → uploaded directly
  - PSP magnetometer data → zipped by year (~2 GB each)
  - Aura/MLS data → zipped as one archive (~6 GB)

Deposits are created as DRAFTS. Run with --publish to finalise.
"""

import requests
import json
import os
import sys
import zipfile
import tempfile
import time
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
ZENODO_API_TOKEN = "v0vwEqX8u9dw6MUFZqAQJSGjwcqA3JImFA5zQbPJx4MIJrhlfQgVp77jJz7p"
ZENODO_API_URL = "https://zenodo.org/api/deposit/depositions"
BASE_DIR = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis")
DATA_DIR = BASE_DIR / "data"

MAX_RETRIES = 3
UPLOAD_TIMEOUT = 7200  # 2 hours for very large files

AUTHOR = {
    "name": "Ashuraliyev, Abduxoliq",
    "affiliation": "Independent Researcher, Tashkent, Uzbekistan",
    "orcid": "0009-0003-5482-5526",
}

GITHUB_REPO = "https://github.com/ProgrmerJack/Solar-Magnetic-Analysis"


# ── Zenodo API helpers ────────────────────────────────────────────────────────
def _auth():
    return {"Authorization": f"Bearer {ZENODO_API_TOKEN}"}


def create_deposit(metadata):
    """Create a new Zenodo draft deposit and return JSON response."""
    r = requests.post(
        ZENODO_API_URL,
        headers={**_auth(), "Content-Type": "application/json"},
        json={"metadata": metadata},
    )
    if r.status_code == 201:
        return r.json()
    print(f"  ✗ Error creating deposit: {r.status_code}\n  {r.text[:400]}")
    return None


def upload_file(bucket_url, filepath, filename=None):
    """Upload a single file to a Zenodo deposit bucket with retry logic."""
    filename = filename or os.path.basename(filepath)
    size_mb = os.path.getsize(filepath) / (1024 * 1024)
    size_str = f"{size_mb:.1f} MB" if size_mb < 1024 else f"{size_mb/1024:.2f} GB"
    print(f"  ↑ Uploading {filename} ({size_str}) …", end=" ", flush=True)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with open(filepath, "rb") as f:
                r = requests.put(
                    f"{bucket_url}/{filename}",
                    headers=_auth(),
                    data=f,
                    timeout=UPLOAD_TIMEOUT,
                )
            if r.status_code in (200, 201):
                print("✓")
                return True
            print(f"✗ HTTP {r.status_code}", end="")
            if attempt < MAX_RETRIES:
                wait = 10 * attempt
                print(f" (retry {attempt}/{MAX_RETRIES} in {wait}s)", end="")
                time.sleep(wait)
            else:
                print(f"\n    {r.text[:300]}")
        except requests.exceptions.RequestException as e:
            print(f"✗ {type(e).__name__}", end="")
            if attempt < MAX_RETRIES:
                wait = 15 * attempt
                print(f" (retry {attempt}/{MAX_RETRIES} in {wait}s)", end="")
                time.sleep(wait)
            else:
                print(f"\n    {e}")
    return False


def upload_direct(bucket_url, filepath, zenodo_name=None):
    """Upload a file directly from disk (no zipping)."""
    return upload_file(bucket_url, str(filepath), zenodo_name)


def publish_deposit(deposit_id):
    """Publish a draft deposit to mint a DOI."""
    r = requests.post(
        f"{ZENODO_API_URL}/{deposit_id}/actions/publish",
        headers=_auth(),
    )
    if r.status_code == 202:
        return r.json()
    print(f"  ✗ Error publishing: {r.status_code}\n  {r.text[:400]}")
    return None


# ── Zip helpers ───────────────────────────────────────────────────────────────
def zip_directory(src_dir, zip_path, base_name=None, exclude_dirs=None):
    """Zip entire *src_dir* into *zip_path* — NO size filtering."""
    exclude_dirs = set(exclude_dirs or [])
    count = 0
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_STORED) as zf:
        src = Path(src_dir)
        for fp in sorted(src.rglob("*")):
            if fp.is_dir():
                continue
            parts = fp.relative_to(src).parts
            if any(p.startswith(".") or p == "__pycache__" for p in parts):
                continue
            if parts and parts[0] in exclude_dirs:
                continue
            arcname = f"{base_name}/{fp.relative_to(src)}" if base_name else str(fp.relative_to(src))
            zf.write(fp, arcname)
            count += 1
    return zip_path, count


def zip_files(file_list, zip_path, base_name=None):
    """Zip a list of (src_path, archive_name) tuples."""
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_STORED) as zf:
        for src, arcname in file_list:
            if os.path.exists(src):
                full_arcname = f"{base_name}/{arcname}" if base_name else arcname
                zf.write(src, full_arcname)
    return zip_path


# ── Data README generators ────────────────────────────────────────────────────
def write_raw_readme(path):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"""# Raw/Source Data — Solar-Magnetic Avalanche Analysis

## Overview
This deposit contains **ALL** raw, unprocessed source data (~24.6 GB) for:

> **Planetary Wave Forcing Suppresses Natural Dry Slab Avalanche Activity
> via Stratospheric Sudden Warmings: Multi-Country Evidence, Process-Model
> Validation, and Count–Rating Dissociation**
>
> Abduxoliq Ashuraliyev (2026)

## Contents

| Archive/File | Size | Description |
|---|---|---|
| `raw_cryosphere.zip` | ~460 MB | Avalanche databases (WSL/SLF, UAC, NVE, EAWS, ALBINA, LAWIS, EnviDat, CAIC) |
| `raw_atmospheric_indices.zip` | ~90 MB | Climate indices (NAO/AO/PNA/QBO), ERA5, POES HPI, TIMED/SABER, volcanic aerosol |
| `raw_aura_mls.zip` | ~6.2 GB | NASA Aura/MLS Level 3 daily zonal mean O3 & Temperature (2004–2024) |
| `raw_solar_other.zip` | ~60 MB | ACE/WIND/DSCOVR, SWPC flare catalog, Parker Solar Probe (non-MAG) |
| `raw_psp_mag_YYYY.zip` | ~0.4–2.7 GB | Parker Solar Probe MAG L2 RTN data, one zip per year (2018–2025) |
| `raw_geomagnetic.zip` | ~11 MB | GFZ Kp/Ap, OMNI Dst |
| `raw_misc.zip` | ~4 MB | ERA5 extended Swiss regional, CAIC accident data |
| `DATA_README.md` | — | This file |

## Reproducing the Analysis
1. Clone: `git clone {GITHUB_REPO}`
2. Extract all archives into `data/` (preserving directory structure)
3. Install dependencies: `pip install -r requirements.txt`
4. Consult `REVIEWER_INDEX.md` for claim → script mapping

## Data Sources & Licences

| Source | URL | Licence |
|---|---|---|
| WSL/SLF (Swiss avalanche) | https://www.envidat.ch/ | CC-BY-4.0 |
| UAC (Utah avalanche) | https://utahavalanchecenter.org/ | Public domain |
| NVE (Norway) | https://api01.nve.no/ | Norwegian Open Data Licence |
| EAWS | https://www.avalanches.org/ | CC-BY-4.0 |
| ALBINA (Tyrol/South Tyrol/Trentino) | https://avalanche.report/ | CC-BY-4.0 |
| Austria LAWIS | https://lawis.at/ | CC-BY-4.0 |
| NOAA CPC (NAO/AO/PNA) | https://www.cpc.ncep.noaa.gov/ | Public domain |
| ECMWF ERA5 | https://cds.climate.copernicus.eu/ | Copernicus Licence |
| NASA Aura/MLS | https://mls.jpl.nasa.gov/ | Public domain |
| NOAA POES/MEPED | https://www.ngdc.noaa.gov/ | Public domain |
| NOAA SWPC | https://www.swpc.noaa.gov/ | Public domain |
| ACE/DSCOVR | https://cdaweb.gsfc.nasa.gov/ | Public domain |
| Parker Solar Probe | https://spdf.gsfc.nasa.gov/ | Public domain |
| GFZ Potsdam | https://www.gfz-potsdam.de/ | CC-BY-4.0 |
| CAIC | https://avalanche.state.co.us/ | Public domain |

## GitHub Repository
{GITHUB_REPO}

## Contact
Abduxoliq Ashuraliyev — Jack00040008@outlook.com
ORCID: 0009-0003-5482-5526
""")
    return path


def write_processed_readme(path):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"""# Processed/Analysis-Ready Data — Solar-Magnetic Avalanche Analysis

## Overview
This deposit contains **ALL** processed datasets (~33.2 GB), statistical
results, and figures for:

> **Planetary Wave Forcing Suppresses Natural Dry Slab Avalanche Activity
> via Stratospheric Sudden Warmings: Multi-Country Evidence, Process-Model
> Validation, and Count–Rating Dissociation**
>
> Abduxoliq Ashuraliyev (2026)

## Contents

### Individual large files (uploaded directly)
| File | Size | Description |
|---|---|---|
| `poes_noaa15_YYYY.parquet` | 0.7–1.6 GB × 13 | POES NOAA-15 MEPED particle precipitation (2013–2025) |
| `poes_noaa18_YYYY.parquet` | 0.5–1.5 GB × 8 | POES NOAA-18 MEPED particle precipitation (2013–2019) |
| `era5_polar_strat_gridded.nc` | 2.3 GB | ERA5 polar stratospheric fields, gridded |
| `goes_r_particle.parquet` | 1.1 GB | GOES-R combined particle data |
| `goes_r_particle_goes16.parquet` | 552 MB | GOES-16 particle data |
| `omni_1min.parquet` | 412 MB | OMNI 1-minute solar wind data |
| `goes_r_particle_goes17.parquet` | 318 MB | GOES-17 particle data |
| `goes_xrs.parquet` | 278 MB | GOES X-ray sensor flare data |
| `psp_mag.parquet` | 265 MB | Parker Solar Probe merged MAG data |
| `goes_r_particle_goes18.parquet` | 229 MB | GOES-18 particle data |

### Zipped collections
| Archive | Size | Description |
|---|---|---|
| `processed_atmospheric_small.zip` | ~200 MB | MLS gridded (O3, T, N2O, HNO3), MERRA-2, NCEP, climate indices, SSW catalog |
| `processed_cryosphere.zip` | ~100 MB | SLF, Utah, Norway avalanche data; SNOTEL SWE; SNOWPACK stability |
| `processed_solar_small.zip` | ~100 MB | Flare catalog, solar indices, SDO/HMI, OMNI hourly, ACE/DSCOVR, CME catalog |
| `processed_geomagnetic.zip` | ~8 MB | Dst and Kp indices |
| `analysis_panels.zip` | ~5 MB | Merged analysis panels, ERA5 Swiss Alps daily |
| `analysis_results.zip` | ~1 MB | 91 JSON statistical result files |
| `figures.zip` | ~7 MB | 76 publication-quality figures (PNG + PDF) |

## Reproducing All Claims
1. Clone: `git clone {GITHUB_REPO}`
2. Download this deposit and extract into `data/processed/`, `data/results/`, `data/figures/`
3. Install: `pip install -r requirements.txt`
4. Consult `REVIEWER_INDEX.md` for claim → script → data traceability

### Key analysis scripts
- `scripts/analysis/r21_paper_analysis.py` — Primary Swiss SSW analysis
- `scripts/analysis/r20_definitive_analysis.py` — Multi-country replication
- `scripts/analysis/r31_snowpack_stability.py` — SNOWPACK validation
- `scripts/analysis/r29_eaws_multicountry.py` — EAWS pan-European analysis

## File Formats
- `.parquet` — Apache Parquet columnar format (readable with pandas, pyarrow)
- `.json` — JSON (analysis results, statistical outputs)
- `.csv` — Comma-separated values
- `.nc` — NetCDF-4 (gridded atmospheric data)
- `.pdf` / `.png` — Figures

## GitHub Repository
{GITHUB_REPO}

## Contact
Abduxoliq Ashuraliyev — Jack00040008@outlook.com
ORCID: 0009-0003-5482-5526
""")
    return path


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    do_publish = "--publish" in sys.argv
    tmp = Path(tempfile.mkdtemp(prefix="zenodo_full_"))
    print(f"Working directory: {tmp}\n")
    results = {}

    # ══════════════════════════════════════════════════════════════════════════
    #  STEP 1 — Build zip archives (for small-file collections)
    # ══════════════════════════════════════════════════════════════════════════
    print("=" * 70)
    print("STEP 1: Building zip archives")
    print("=" * 70)

    def report_zip(label, path, count=None):
        sz = os.path.getsize(path)
        s = f"{sz/1e6:.1f} MB" if sz < 1e9 else f"{sz/1e9:.2f} GB"
        extra = f", {count} files" if count else ""
        print(f"  ✓ {label}: {s}{extra}")

    # ── RAW ZIPS ──────────────────────────────────────────────────────────────

    # 1. Cryosphere — all subdirs
    print("\n  [RAW] Zipping cryosphere …")
    p, n = zip_directory(DATA_DIR / "cryosphere", tmp / "raw_cryosphere.zip",
                         base_name="cryosphere")
    report_zip("raw_cryosphere.zip", p, n)

    # 2. Atmospheric (excluding aura_mls — that gets its own zip)
    print("  [RAW] Zipping atmospheric indices …")
    p, n = zip_directory(DATA_DIR / "atmospheric", tmp / "raw_atmospheric_indices.zip",
                         base_name="atmospheric", exclude_dirs={"aura_mls"})
    report_zip("raw_atmospheric_indices.zip", p, n)

    # 3. Aura/MLS — ~6.2 GB, 88 files
    print("  [RAW] Zipping Aura/MLS (6.2 GB, ~88 files) …")
    p, n = zip_directory(DATA_DIR / "atmospheric" / "aura_mls",
                         tmp / "raw_aura_mls.zip", base_name="atmospheric/aura_mls")
    report_zip("raw_aura_mls.zip", p, n)

    # 4. Solar (excluding psp_mag)
    print("  [RAW] Zipping solar (non-PSP) …")
    p, n = zip_directory(DATA_DIR / "solar", tmp / "raw_solar_other.zip",
                         base_name="solar", exclude_dirs={"psp_mag"})
    report_zip("raw_solar_other.zip", p, n)

    # 5. PSP MAG — zip by year (each ~0.4–2.7 GB)
    psp_dir = DATA_DIR / "solar" / "psp_mag"
    psp_zips = []
    for year_dir in sorted(psp_dir.iterdir()):
        if year_dir.is_dir() and year_dir.name.isdigit():
            yr = year_dir.name
            print(f"  [RAW] Zipping PSP MAG {yr} …")
            zname = f"raw_psp_mag_{yr}.zip"
            p, n = zip_directory(year_dir, tmp / zname,
                                 base_name=f"solar/psp_mag/{yr}")
            psp_zips.append(zname)
            report_zip(zname, p, n)

    # 6. Geomagnetic
    print("  [RAW] Zipping geomagnetic …")
    p, n = zip_directory(DATA_DIR / "geomagnetic", tmp / "raw_geomagnetic.zip",
                         base_name="geomagnetic")
    report_zip("raw_geomagnetic.zip", p, n)

    # 7. Misc (data/raw/)
    print("  [RAW] Zipping misc (data/raw/) …")
    p, n = zip_directory(DATA_DIR / "raw", tmp / "raw_misc.zip", base_name="raw")
    report_zip("raw_misc.zip", p, n)

    # ── PROCESSED ZIPS ────────────────────────────────────────────────────────

    # Small atmospheric files (everything except the large POES/ERA5/gridded files)
    print("\n  [PROC] Zipping atmospheric small files …")
    atm_dir = DATA_DIR / "processed" / "atmospheric"
    large_atm_prefixes = ("poes_", "era5_polar_strat_gridded", "merra2_polar_strat.nc",
                          "mls_ozone_gridded", "mls_temperature_gridded",
                          "mls_n2o_gridded", "mls_hno3_gridded")
    atm_small = []
    for fp in sorted(atm_dir.iterdir()):
        if fp.is_file() and not fp.name.startswith(large_atm_prefixes):
            atm_small.append((str(fp), f"processed/atmospheric/{fp.name}"))
    # Also include the medium gridded MLS/MERRA2 files (15–85 MB)
    for name in ["merra2_polar_strat.nc", "mls_ozone_gridded.nc",
                 "mls_temperature_gridded.nc", "mls_n2o_gridded.nc",
                 "mls_hno3_gridded.nc"]:
        fp = atm_dir / name
        if fp.exists():
            atm_small.append((str(fp), f"processed/atmospheric/{name}"))
    zip_files(atm_small, tmp / "processed_atmospheric_small.zip")
    report_zip("processed_atmospheric_small.zip", tmp / "processed_atmospheric_small.zip",
               len(atm_small))

    # Cryosphere
    print("  [PROC] Zipping cryosphere …")
    p, n = zip_directory(DATA_DIR / "processed" / "cryosphere",
                         tmp / "processed_cryosphere.zip",
                         base_name="processed/cryosphere")
    report_zip("processed_cryosphere.zip", p, n)

    # Solar — small files (exclude the large ones that upload individually)
    print("  [PROC] Zipping solar small files …")
    sol_dir = DATA_DIR / "processed" / "solar"
    large_solar = {"goes_r_particle.parquet", "goes_r_particle_goes16.parquet",
                   "omni_1min.parquet", "goes_r_particle_goes17.parquet",
                   "goes_xrs.parquet", "psp_mag.parquet",
                   "goes_r_particle_goes18.parquet", "goes_legacy_particle.parquet"}
    sol_small = []
    for fp in sorted(sol_dir.iterdir()):
        if fp.is_file() and fp.name not in large_solar:
            sol_small.append((str(fp), f"processed/solar/{fp.name}"))
    zip_files(sol_small, tmp / "processed_solar_small.zip")
    report_zip("processed_solar_small.zip", tmp / "processed_solar_small.zip",
               len(sol_small))

    # Geomagnetic
    print("  [PROC] Zipping geomagnetic …")
    p, n = zip_directory(DATA_DIR / "processed" / "geomagnetic",
                         tmp / "processed_geomagnetic.zip",
                         base_name="processed/geomagnetic")
    report_zip("processed_geomagnetic.zip", p, n)

    # Analysis panels
    print("  [PROC] Zipping analysis panels …")
    panel_files = []
    for fn in ["analysis_panel.parquet", "analysis_panel_v2.parquet",
               "era5_swiss_alps_daily.parquet", "era5_swiss_alps_extended.parquet",
               "manifest.json", "preregistration.json"]:
        fp = DATA_DIR / "processed" / fn
        if fp.exists():
            panel_files.append((str(fp), fn))
    zip_files(panel_files, tmp / "analysis_panels.zip", base_name="processed")
    report_zip("analysis_panels.zip", tmp / "analysis_panels.zip", len(panel_files))

    # Results
    print("  [PROC] Zipping results …")
    p, n = zip_directory(DATA_DIR / "results", tmp / "analysis_results.zip",
                         base_name="results")
    report_zip("analysis_results.zip", p, n)

    # Figures
    print("  [PROC] Zipping figures …")
    p, n = zip_directory(DATA_DIR / "figures", tmp / "figures.zip",
                         base_name="figures")
    report_zip("figures.zip", p, n)

    # ── READMEs ───────────────────────────────────────────────────────────────
    write_raw_readme(tmp / "DATA_README_RAW.md")
    write_processed_readme(tmp / "DATA_README_PROC.md")
    print("\n  ✓ README files written")

    # ══════════════════════════════════════════════════════════════════════════
    #  STEP 2 — Create Zenodo deposits and upload EVERYTHING
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("STEP 2: Creating Zenodo deposits and uploading")
    print("=" * 70)

    # ── DEPOSIT 1: Raw/Source Data (~24.6 GB) ─────────────────────────────────
    print("\n── DEPOSIT 1: Raw/Source Data (~24.6 GB) ──")
    meta_raw = {
        "title": (
            "Raw Source Data: Planetary Wave Forcing Suppresses Natural "
            "Dry Slab Avalanche Activity via Stratospheric Sudden Warmings"
        ),
        "upload_type": "dataset",
        "description": (
            "Complete raw source data (~24.6 GB) for the study linking "
            "Stratospheric Sudden Warming (SSW) events to reduced natural "
            "dry slab avalanche activity across four countries. Contains: "
            "avalanche databases (WSL/SLF Switzerland, UAC Utah, NVE Norway, "
            "EAWS pan-European, ALBINA Tyrol, LAWIS Austria), atmospheric "
            "reanalysis (ERA5, NCEP, Aura/MLS O3 & Temperature), climate "
            "indices (NAO, AO, PNA, QBO), solar wind and particle data "
            "(ACE/DSCOVR, GOES, Parker Solar Probe including full MAG L2 "
            "RTN ~17.8 GB), geomagnetic indices (Dst, Kp), and snow "
            "observations (SNOTEL, EnviDat SNOWPACK, MODIS). "
            "ALL source files are included — no external downloads required. "
            "Analysis code: " + GITHUB_REPO
        ),
        "creators": [AUTHOR],
        "keywords": [
            "stratospheric sudden warming", "avalanche", "SSW",
            "snow avalanche", "ERA5", "NCEP", "SNOTEL", "Aura MLS",
            "Parker Solar Probe", "WSL SLF", "EAWS", "solar wind",
            "geomagnetic", "cryosphere", "raw data", "reproducibility",
        ],
        "license": "cc-by-4.0",
        "access_right": "open",
        "related_identifiers": [
            {"identifier": GITHUB_REPO, "relation": "isSupplementTo",
             "scheme": "url"}
        ],
        "notes": (
            "COMPLETE raw data deposit — all files included. Data assembled "
            "from public archives. See DATA_README.md for individual source "
            "licences and extraction instructions."
        ),
    }

    dep_raw = create_deposit(meta_raw)
    if not dep_raw:
        print("FATAL: Could not create raw deposit.")
        sys.exit(1)

    print(f"  ✓ Created deposit ID: {dep_raw['id']}")
    bucket_raw = dep_raw["links"]["bucket"]

    # Upload all raw zips
    raw_zips = [
        "raw_cryosphere.zip", "raw_atmospheric_indices.zip",
        "raw_aura_mls.zip", "raw_solar_other.zip",
        "raw_geomagnetic.zip", "raw_misc.zip",
    ] + psp_zips  # raw_psp_mag_2018.zip .. raw_psp_mag_2025.zip

    for zname in raw_zips:
        zpath = tmp / zname
        if zpath.exists() and zpath.stat().st_size > 100:
            upload_file(bucket_raw, str(zpath), zname)
            time.sleep(1)

    upload_file(bucket_raw, str(tmp / "DATA_README_RAW.md"), "DATA_README.md")

    if do_publish:
        pub = publish_deposit(dep_raw["id"])
        results["raw"] = ({"doi": pub["doi"], "doi_url": pub["doi_url"],
                           "id": dep_raw["id"]} if pub
                          else {"status": "DRAFT", "id": dep_raw["id"]})
    else:
        results["raw"] = {"status": "DRAFT", "id": dep_raw["id"]}
        print(f"  → Draft: https://zenodo.org/deposit/{dep_raw['id']}")

    # ── DEPOSIT 2: Processed/Analysis-Ready Data (~33.2 GB) ───────────────────
    print("\n── DEPOSIT 2: Processed/Analysis-Ready Data (~33.2 GB) ──")
    meta_proc = {
        "title": (
            "Processed Data and Results: Planetary Wave Forcing Suppresses "
            "Natural Dry Slab Avalanche Activity via Stratospheric Sudden "
            "Warmings"
        ),
        "upload_type": "dataset",
        "description": (
            "Complete processed datasets (~33.2 GB), statistical results, "
            "and figures. Includes POES NOAA-15/18 MEPED particle "
            "precipitation data (21 yearly parquets, ~28 GB), ERA5 polar "
            "stratospheric fields (2.3 GB), GOES-R particle data, solar "
            "indices, avalanche daily counts, SNOWPACK stability metrics, "
            "91 JSON result files, and 76 publication-quality figures. "
            "ALL processed files included — no regeneration required. "
            "Rate ratio = 0.32, P = 0.004 (sign test), n = 16 SSW events "
            "over 21 winters across Switzerland, Utah, Norway, and "
            "pan-European EAWS regions. "
            "Analysis code: " + GITHUB_REPO
        ),
        "creators": [AUTHOR],
        "keywords": [
            "stratospheric sudden warming", "avalanche", "SSW",
            "processed data", "POES MEPED", "ERA5", "GOES",
            "analysis results", "reproducibility", "rate ratio",
            "sign test", "SNOWPACK", "sintering model",
            "specification curve", "meta-analysis",
        ],
        "license": "cc-by-4.0",
        "access_right": "open",
        "related_identifiers": [
            {"identifier": GITHUB_REPO, "relation": "isSupplementTo",
             "scheme": "url"}
        ],
        "notes": (
            "COMPLETE processed data deposit — all files included. "
            "Processed from raw source data (see companion deposit). "
            "Consult REVIEWER_INDEX.md in the GitHub repository for "
            "claim-to-script-to-data traceability."
        ),
    }

    dep_proc = create_deposit(meta_proc)
    if not dep_proc:
        print("FATAL: Could not create processed deposit.")
        sys.exit(1)

    print(f"  ✓ Created deposit ID: {dep_proc['id']}")
    bucket_proc = dep_proc["links"]["bucket"]

    # Upload large individual atmospheric files
    print("\n  ── Large atmospheric files (individual upload) ──")
    atm_dir = DATA_DIR / "processed" / "atmospheric"
    for fp in sorted(atm_dir.iterdir()):
        if fp.is_file() and (fp.name.startswith("poes_") or
                             fp.name == "era5_polar_strat_gridded.nc"):
            upload_direct(bucket_proc, fp,
                          f"processed_atmospheric_{fp.name}")
            time.sleep(1)

    # Upload large individual solar files
    print("\n  ── Large solar files (individual upload) ──")
    for name in sorted(large_solar):
        fp = sol_dir / name
        if fp.exists():
            upload_direct(bucket_proc, fp, f"processed_solar_{name}")
            time.sleep(1)

    # Upload all zipped collections
    print("\n  ── Zipped collections ──")
    proc_zips = [
        "processed_atmospheric_small.zip", "processed_cryosphere.zip",
        "processed_solar_small.zip", "processed_geomagnetic.zip",
        "analysis_panels.zip", "analysis_results.zip", "figures.zip",
    ]
    for zname in proc_zips:
        zpath = tmp / zname
        if zpath.exists() and zpath.stat().st_size > 100:
            upload_file(bucket_proc, str(zpath), zname)
            time.sleep(1)

    upload_file(bucket_proc, str(tmp / "DATA_README_PROC.md"), "DATA_README.md")

    if do_publish:
        pub = publish_deposit(dep_proc["id"])
        results["processed"] = (
            {"doi": pub["doi"], "doi_url": pub["doi_url"],
             "id": dep_proc["id"]} if pub
            else {"status": "DRAFT", "id": dep_proc["id"]})
    else:
        results["processed"] = {"status": "DRAFT", "id": dep_proc["id"]}
        print(f"  → Draft: https://zenodo.org/deposit/{dep_proc['id']}")

    # ══════════════════════════════════════════════════════════════════════════
    #  SUMMARY
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for name, info in results.items():
        status = info.get("doi", info.get("status", "UNKNOWN"))
        print(f"  {name.upper():12s}: {status}  (id={info.get('id', '?')})")
        if "doi_url" in info:
            print(f"               {info['doi_url']}")

    out_path = BASE_DIR / "ZENODO_DEPOSITS.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")
    if not do_publish:
        print("\n⚠  Deposits are DRAFTS. Review at https://zenodo.org/deposit/")
        print("   Re-run with --publish to finalise and mint DOIs.")

    return results


if __name__ == "__main__":
    main()
