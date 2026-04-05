"""
download_all.py — Master Download Orchestrator
Solar-Magnetic-Analysis Research Project

Runs all download scripts in sequence, logging progress.
Scripts that can run immediately (no credentials):
  01 — Solar indices (SILSO, GFZ Kp, SWPC)
  02 — OMNIWeb (solar wind + Kp/Dst/AE, 1963–present)
  03 — GOES XRS (GOES-16/18 X-ray flux, 2017–present)
  04 — Flare catalog (SWPC events, HEK, NGDC)
  05 — SSW catalog (Butler 2015, hard-coded + fetch attempts)
  06 — Avalanche / CAIC / EAWS / EnviDat
  07 — SNOTEL snowpack network (NRCS API)
  08 — DSCOVR/ACE solar wind (NOAA archive + CDAWeb)
  09 — POES/MEPED particle precipitation (NOAA)
  12 — Geomagnetic indices (Dst, AE supplemental)
  14 — Norway avalanche (Varsom, Regobs, NVE)
  13 — PSP FIELDS data (NASA SPDF, no auth)

Setup-only scripts (credentials required):
  10 — ERA5 (Copernicus CDS API)
  11 — NASA Earthdata (MERRA-2, Aura/MLS, MODIS, TIMED/SABER)
"""
import sys
import subprocess
import time
import json
from pathlib import Path
from datetime import datetime

BASE_DIR = Path(__file__).resolve().parents[2]
SCRIPTS  = Path(__file__).parent
LOG_DIR  = BASE_DIR / "logs"
LOG_DIR.mkdir(exist_ok=True)

DOWNLOAD_SCRIPTS = [
    ("01_download_solar_indices.py",   "Solar Indices (SILSO, GFZ Kp, SWPC)",          True),
    ("02_download_omni.py",            "OMNIWeb combined solar wind + indices",          True),
    ("05_download_ssw_catalog.py",     "SSW Catalog (Butler 2015)",                      True),
    ("06_download_caic_avalanche.py",  "CAIC / EAWS / avalanche data",                  True),
    ("07_download_snotel.py",          "SNOTEL snowpack telemetry",                      True),
    ("08_download_dscovr_ace.py",      "DSCOVR / ACE solar wind",                       True),
    ("09_download_poes_meped.py",      "POES/MEPED particle precipitation",             True),
    ("12_setup_geomagnetic.py",        "SuperMAG / INTERMAGNET / Dst / AE setup",       True),
    ("14_download_norway_data.py",     "Norway avalanche (Varsom, Regobs, NVE)",        True),
    ("04_download_flare_catalog.py",   "NOAA flare catalog + HEK (large, last)",        True),
    ("03_download_goes_xrs.py",        "GOES-16/18 XRS NetCDF (large, last)",           True),
    ("13_setup_psp_sdo.py",            "Parker Solar Probe + SDO/JSOC setup",           True),
    # Setup-only (no credentials yet):
    ("10_setup_era5.py",               "ERA5 setup (credentials required)",             True),
    ("11_setup_nasa_earthdata.py",     "NASA Earthdata setup (credentials required)",   True),
]


def run_script(script_name: str, desc: str) -> dict:
    script_path = SCRIPTS / script_name
    if not script_path.exists():
        return {"script": script_name, "status": "MISSING", "elapsed": 0}

    start = time.time()
    print(f"\n{'='*70}")
    print(f"  Running: {script_name}")
    print(f"  Purpose: {desc}")
    print(f"  Started: {datetime.now().strftime('%H:%M:%S')}")
    print(f"{'='*70}")

    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=False,
            timeout=3600,  # 1 hour max per script
        )
        elapsed = time.time() - start
        status = "OK" if result.returncode == 0 else f"EXIT_{result.returncode}"
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        status = "TIMEOUT"
    except Exception as exc:
        elapsed = time.time() - start
        status = f"ERROR: {exc}"

    print(f"\n  → {status} in {elapsed:.1f}s")
    return {"script": script_name, "status": status, "elapsed": round(elapsed, 1), "desc": desc}


def main():
    print(f"\n{'#'*70}")
    print(f"  Solar-Magnetic-Analysis — Dataset Download")
    print(f"  Started: {datetime.now().isoformat()}")
    print(f"{'#'*70}\n")

    results = []
    for script, desc, enabled in DOWNLOAD_SCRIPTS:
        if not enabled:
            print(f"  SKIP  {script}")
            continue
        result = run_script(script, desc)
        results.append(result)

    # Save summary
    summary_path = LOG_DIR / f"download_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    summary_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    print(f"\n\n{'#'*70}")
    print(f"  DOWNLOAD SUMMARY")
    print(f"{'#'*70}")
    ok = [r for r in results if r["status"] == "OK"]
    fail = [r for r in results if r["status"] != "OK"]
    print(f"  ✓  {len(ok)} scripts succeeded")
    if fail:
        print(f"  ✗  {len(fail)} scripts had issues:")
        for r in fail:
            print(f"      {r['script']:45s}  {r['status']}")
    total = sum(r["elapsed"] for r in results)
    print(f"\n  Total time: {total/60:.1f} minutes")
    print(f"  Summary saved: {summary_path}")
    print(f"{'#'*70}\n")


if __name__ == "__main__":
    main()
