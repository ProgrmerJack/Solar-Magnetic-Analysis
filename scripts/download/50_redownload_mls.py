"""
Re-download Aura/MLS L3 Daily Zonal files for all 4 EPP-relevant products.
Saves to data/atmospheric/aura_mls/{product}/ sub-directories.
Only downloads files not already on disk.
"""
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))

os.environ.setdefault("EARTHDATA_USERNAME", "andersonmark578")
os.environ.setdefault("EARTHDATA_PASSWORD", "Nasa.20080408@#")

import earthaccess

auth = earthaccess.login(strategy="environment", persist=False)
print(f"Authenticated: {auth.authenticated}")

BASE = Path(__file__).parents[2] / "data"
MLS_OUT = BASE / "atmospheric" / "aura_mls"
MLS_OUT.mkdir(parents=True, exist_ok=True)

TEMPORAL = ("2004-08-13", "2026-04-05")

# MLS products: (short_name, version, subdir, description)
MLS_PRODUCTS = [
    ("ML3DZHNO3",  "005", "HNO3",        "MLS L3 Daily Zonal HNO3"),
    ("ML3DZT",     "005", "Temperature", "MLS L3 Daily Zonal Temperature"),
    ("ML3DZN2O",   "005", "N2O",         "MLS L3 Daily Zonal N2O"),
    ("ML3DZO3",    "005", "O3",          "MLS L3 Daily Zonal Ozone"),
]


def download_mls_product(short_name, version, subdir, label):
    out_dir = MLS_OUT / subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Searching {label} ({short_name} v{version}) ...")
    granules = earthaccess.search_data(
        short_name=short_name,
        version=version,
        temporal=TEMPORAL,
        count=-1,
    )
    print(f"  Found {len(granules)} granules")

    pending = []
    for g in granules:
        links = g.data_links(access="onprem") or g.data_links()
        if not links:
            continue
        fname = Path(links[0]).name
        if not (out_dir / fname).exists():
            pending.append(g)

    print(f"  Already on disk: {len(granules) - len(pending)}")
    print(f"  To download: {len(pending)}")

    if not pending:
        print(f"  SKIP: all files present")
        return

    # Download in batches of 20 to avoid timeouts
    batch_size = 20
    total_dl = 0
    for i in range(0, len(pending), batch_size):
        batch = pending[i:i + batch_size]
        print(f"  Downloading batch {i//batch_size + 1}/{(len(pending)+batch_size-1)//batch_size} "
              f"({len(batch)} files) ...")
        try:
            earthaccess.download(batch, str(out_dir))
            total_dl += len(batch)
        except Exception as e:
            print(f"  WARN batch error: {e}")

    # Count final
    final = list(out_dir.glob("*.nc")) + list(out_dir.glob("*.he5")) + list(out_dir.glob("*.h5"))
    total_mb = sum(f.stat().st_size for f in final) / 1e6
    print(f"  Done: {len(final)} files on disk, {total_mb:.0f} MB total")


if __name__ == "__main__":
    for short_name, version, subdir, label in MLS_PRODUCTS:
        download_mls_product(short_name, version, subdir, label)
    print("\nMLS re-download complete.")
