"""
Download Dst from OMNI hourly low-resolution data (OMNI2 format).
Source: https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_YYYY.dat
OMNI2 column 41 (0-indexed: col 40) = Dst index (nT), fill=99999

Only saves Dst + Kp to avoid OOM on 60+ years of full OMNI2 data.
"""
import requests
import pandas as pd
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))

OUT = Path(__file__).parents[2] / "data" / "geomagnetic"
OUT.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni"
YEARS = list(range(1963, 2026))

# OMNI2 format: fixed-width, columns 0-indexed (55 total)
# Col 0: year, col 1: doy, col 2: hour
# Verified from 2003 Halloween storm (Oct 29 20:00 UTC, known Dst ≈ -253 nT)
DST_COL = 40  # 0-indexed → col 41 (1-indexed) = DST Index (nT), fill=99999
KP_COL  = 38  # 0-indexed → col 39 (1-indexed) = Kp*10, fill=99
AE_COL  = 41  # 0-indexed → col 42 (1-indexed) = AE Index (nT), fill=9999

DST_FILL = 99999
KP_FILL  = 99

all_rows = []
session = requests.Session()
session.headers["User-Agent"] = "SolarMagAvalanche/1.0"

print(f"Downloading OMNI2 hourly for {YEARS[0]}-{YEARS[-1]} ...")

for yr in YEARS:
    url = f"{BASE_URL}/omni2_{yr}.dat"
    try:
        r = session.get(url, timeout=60)
        r.raise_for_status()
        for line in r.text.splitlines():
            parts = line.split()
            if len(parts) < 45:
                continue
            year_val = int(parts[0])
            doy      = int(parts[1])
            hour     = int(parts[2])
            dst_raw  = float(parts[DST_COL])
            kp_raw   = float(parts[KP_COL])
            ae_raw   = float(parts[AE_COL]) if AE_COL < len(parts) else np.nan

            dst = np.nan if dst_raw >= DST_FILL else dst_raw
            kp  = np.nan if kp_raw  >= KP_FILL  else kp_raw
            ae  = np.nan if ae_raw  >= 9999      else ae_raw

            # Build timestamp
            ts = pd.Timestamp(year=year_val, month=1, day=1, tz="UTC") + \
                 pd.Timedelta(days=doy - 1, hours=hour)
            all_rows.append({"time": ts, "dst": dst, "kp_omni": kp, "ae": ae})

        print(f"  ✓ {yr}: {len(r.text.splitlines())} hours")
    except Exception as e:
        print(f"  WARN {yr}: {e}")

df = pd.DataFrame(all_rows).set_index("time").sort_index()
df.index = df.index.tz_convert("UTC")

out_path = OUT / "omni_hourly_dst.parquet"
df.to_parquet(out_path, engine="pyarrow", compression="snappy")

n_dst = df["dst"].notna().sum()
dst_min = df["dst"].min()
dst_max = df["dst"].max()
print(f"\nSaved: {len(df)} rows, Dst valid={n_dst} ({n_dst/len(df):.1%})")
print(f"Dst range: [{dst_min:.0f}, {dst_max:.0f}] nT")
print(f"Output: {out_path}")
