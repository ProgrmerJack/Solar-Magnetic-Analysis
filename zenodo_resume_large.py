#!/usr/bin/env python3
"""
Resume large-file uploads to Zenodo using curl.exe (avoids Windows TCP abort).
Run AFTER zenodo_upload.py to fill in any files that failed.
Usage: python zenodo_resume_large.py
"""
import requests
import subprocess
import os
import sys
import json
import time
from pathlib import Path

TOKEN = "v0vwEqX8u9dw6MUFZqAQJSGjwcqA3JImFA5zQbPJx4MIJrhlfQgVp77jJz7p"
DEPOSITS_JSON = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis\ZENODO_DEPOSITS.json")

def _auth():
    return {"Authorization": "Bearer " + TOKEN}

def get_deposit_info(dep_id):
    r = requests.get(f"https://zenodo.org/api/deposit/depositions/{dep_id}",
                     headers=_auth())
    r.raise_for_status()
    return r.json()

def list_uploaded(dep_id, bucket_override=None):
    """Return set of filenames already in the deposit."""
    d = get_deposit_info(dep_id)
    bucket = bucket_override or d["links"]["bucket"]
    return {f["filename"] for f in d.get("files", [])}, bucket

def upload_curl(bucket_url, filepath, filename=None):
    """Upload using curl.exe — handles large files reliably on Windows."""
    filename = filename or os.path.basename(filepath)
    size = os.path.getsize(filepath)
    size_str = f"{size/1e9:.2f} GB" if size >= 1e9 else f"{size/1e6:.1f} MB"
    url = f"{bucket_url}/{filename}"
    print(f"  ↑ curl: {filename} ({size_str}) …", flush=True)
    cmd = [
        "curl.exe", "-X", "PUT",
        "-H", f"Authorization: Bearer {TOKEN}",
        "-H", "Content-Type: application/octet-stream",
        "--upload-file", str(filepath),
        "--retry", "3",
        "--retry-delay", "10",
        "--retry-max-time", "7200",
        "--max-time", "14400",   # 4 hours max
        "--speed-limit", "1024", # abort if <1KB/s for 30s
        "--speed-time", "30",
        "--progress-bar",
        url,
    ]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode == 0:
        print(f"  ✓ {filename}")
        return True
    else:
        print(f"  ✗ curl exit code {result.returncode}")
        return False

def main():
    if not DEPOSITS_JSON.exists():
        print("ERROR: ZENODO_DEPOSITS.json not found. Run zenodo_upload.py first.")
        sys.exit(1)

    with open(DEPOSITS_JSON) as f:
        deps = json.load(f)

    dep_raw_id = deps["raw"]["id"]
    dep_proc_id = deps["processed"]["id"]
    bucket_raw_override = deps["raw"].get("bucket")
    bucket_proc_override = deps["processed"].get("bucket")

    print(f"Raw deposit:       {dep_raw_id}")
    print(f"Processed deposit: {dep_proc_id}")

    # ── Find temp dir with zip files ──────────────────────────────────────────
    import tempfile, glob as g
    tmp_dirs = sorted(g.glob(os.path.join(tempfile.gettempdir(), "zenodo_full_*")),
                      key=os.path.getmtime, reverse=True)
    if not tmp_dirs:
        print("ERROR: No zenodo_full_* temp directory found. Zips may have been deleted.")
        print("Re-run zenodo_upload.py to recreate them, or locate zip files manually.")
        sys.exit(1)
    tmp = Path(tmp_dirs[0])
    print(f"Using temp dir: {tmp}\n")

    # ══════════════════════════════════════════════════════════════════════════
    # DEPOSIT 1 — Raw data: upload missing files
    # ══════════════════════════════════════════════════════════════════════════
    print("=" * 60)
    print(f"DEPOSIT 1 (Raw, id={dep_raw_id})")
    print("=" * 60)
    uploaded_raw, bucket_raw = list_uploaded(dep_raw_id, bucket_raw_override)
    print(f"Already uploaded: {sorted(uploaded_raw)}\n")

    raw_files = [
        "raw_cryosphere.zip",
        "raw_atmospheric_indices.zip",
        "raw_aura_mls.zip",
        "raw_solar_other.zip",
        "raw_geomagnetic.zip",
        "raw_misc.zip",
        "raw_psp_mag_2018.zip",
        "raw_psp_mag_2019.zip",
        "raw_psp_mag_2020.zip",
        "raw_psp_mag_2021.zip",
        "raw_psp_mag_2022.zip",
        "raw_psp_mag_2023.zip",
        "raw_psp_mag_2024.zip",
        "raw_psp_mag_2025.zip",
        "DATA_README.md",
    ]

    for fname in raw_files:
        if fname in uploaded_raw:
            print(f"  ✓ already uploaded: {fname}")
            continue
        fpath = tmp / fname
        # README lives under a different name in the temp dir
        if fname == "DATA_README.md":
            fpath = tmp / "DATA_README_RAW.md"
        if not fpath.exists():
            print(f"  ✗ not found in temp: {fname}")
            continue
        if fpath.stat().st_size < 100:
            print(f"  - skipping empty: {fname}")
            continue
        upload_curl(bucket_raw, str(fpath), fname)
        time.sleep(2)

    # ══════════════════════════════════════════════════════════════════════════
    # DEPOSIT 2 — Processed data: upload missing files
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 60)
    print(f"DEPOSIT 2 (Processed, id={dep_proc_id})")
    print("=" * 60)
    uploaded_proc, bucket_proc = list_uploaded(dep_proc_id, bucket_proc_override)
    print(f"Already uploaded: {sorted(uploaded_proc)}\n")

    atm_dir = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric")
    sol_dir = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\solar")

    # Large individual files (uploaded directly, not via zip)
    large_atm = [f for f in sorted(atm_dir.iterdir())
                 if f.is_file() and (f.name.startswith("poes_") or
                                     f.name == "era5_polar_strat_gridded.nc")]
    large_sol_names = {
        "goes_r_particle.parquet", "goes_r_particle_goes16.parquet",
        "omni_1min.parquet", "goes_r_particle_goes17.parquet",
        "goes_xrs.parquet", "psp_mag.parquet",
        "goes_r_particle_goes18.parquet", "goes_legacy_particle.parquet",
    }
    large_sol = [sol_dir / n for n in sorted(large_sol_names)
                 if (sol_dir / n).exists()]

    proc_zips = [
        "processed_atmospheric_small.zip",
        "processed_cryosphere.zip",
        "processed_solar_small.zip",
        "processed_geomagnetic.zip",
        "analysis_panels.zip",
        "analysis_results.zip",
        "figures.zip",
    ]

    # Upload large atmospheric files
    print("  -- Large atmospheric files --")
    for fp in large_atm:
        zenodo_name = f"processed_atmospheric_{fp.name}"
        if zenodo_name in uploaded_proc:
            print(f"  ✓ already uploaded: {zenodo_name}")
            continue
        upload_curl(bucket_proc, str(fp), zenodo_name)
        time.sleep(2)

    # Upload large solar files
    print("\n  -- Large solar files --")
    for fp in large_sol:
        zenodo_name = f"processed_solar_{fp.name}"
        if zenodo_name in uploaded_proc:
            print(f"  ✓ already uploaded: {zenodo_name}")
            continue
        upload_curl(bucket_proc, str(fp), zenodo_name)
        time.sleep(2)

    # Upload zipped collections
    print("\n  -- Zipped collections --")
    for fname in proc_zips:
        if fname in uploaded_proc:
            print(f"  ✓ already uploaded: {fname}")
            continue
        fpath = tmp / fname
        if not fpath.exists():
            print(f"  ✗ not found: {fname}")
            continue
        upload_curl(bucket_proc, str(fpath), fname)
        time.sleep(2)

    # README
    readme_name = "DATA_README.md"
    if readme_name not in uploaded_proc:
        upload_curl(bucket_proc, str(tmp / "DATA_README_PROC.md"), readme_name)

    # ══════════════════════════════════════════════════════════════════════════
    # Final status
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 60)
    print("FINAL STATUS")
    print("=" * 60)
    for name, dep_id in [("RAW", dep_raw_id), ("PROCESSED", dep_proc_id)]:
        uploaded, _ = list_uploaded(dep_id)
        print(f"  {name} (id={dep_id}): {len(uploaded)} files")
        for f in sorted(uploaded):
            print(f"    ✓ {f}")
    print(f"\nView drafts: https://zenodo.org/deposit/")
    print("Re-run zenodo_upload.py --publish to mint DOIs.")

if __name__ == "__main__":
    main()
