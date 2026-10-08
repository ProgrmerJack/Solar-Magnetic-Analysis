#!/usr/bin/env python3
"""
acquire_tc_data.py
==================
Inputs of 07_physical_decomposition/tc_class_law.py (registered 2026-10-08, committed
before any of these values was fetched).
  IBTrACS v04r01, all basins (NCEI, anonymous):
    .../ibtracs/v04r01/access/csv/ibtracs.ALL.list.v04r01.csv
  ERSST v5 monthly mean SST (NOAA PSL, anonymous):
    https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc
Raw files go to raw/ibtracs/ and raw/ersst/ (not committed). Per storm (SID), seasons
1982-2024, records at synoptic hours 00/06/12/18 UTC with a USA_WIND value:
  y_lmi    maximum USA_WIND (kt)
  y_ri     largest V(t+24h) - V(t) with both records over water (DIST2LAND > 0); NaN if
           the storm has no such pair
  genesis  first record with USA_WIND >= 34 kt (storms never reaching 34 kt are dropped)
  sst_g    ERSST v5 at the genesis month and the nearest grid cell with a value
           (great-circle distance), deg C
  basin, season, lat_g, lon_g, month_g kept for description.
Writes tc_storms.parquet next to this script.
"""
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
IBT = RAW / "ibtracs" / "ibtracs.ALL.list.v04r01.csv"
SST = RAW / "ersst" / "sst.mnmean.nc"
OUT = HERE / "tc_storms.parquet"
URL_IBT = ("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/"
           "v04r01/access/csv/ibtracs.ALL.list.v04r01.csv")
URL_SST = "https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"


def fetch(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    tmp.rename(dest)
    print(f"fetched {dest.name}: {dest.stat().st_size / 1e6:.1f} MB")


def storms():
    cols = ["SID", "SEASON", "BASIN", "ISO_TIME", "LAT", "LON", "DIST2LAND", "USA_WIND"]
    d = pd.read_csv(IBT, usecols=cols, skiprows=[1], keep_default_na=False, na_values=[" ", ""], low_memory=False)
    d["SEASON"] = pd.to_numeric(d.SEASON, errors="coerce")
    d = d[(d.SEASON >= 1982) & (d.SEASON <= 2024)].copy()
    for c in ("LAT", "LON", "DIST2LAND", "USA_WIND"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["t"] = pd.to_datetime(d.ISO_TIME)
    d = d[(d.t.dt.hour % 6 == 0) & (d.t.dt.minute == 0) & d.USA_WIND.notna()]
    out = []
    for sid, s in d.groupby("SID", sort=False):
        s = s.sort_values("t").drop_duplicates("t")
        if (s.USA_WIND >= 34).sum() == 0:
            continue
        g = s[s.USA_WIND >= 34].iloc[0]
        v = s.set_index("t")
        w = v.USA_WIND.where(v.DIST2LAND > 0)
        later = w.reindex(w.index + pd.Timedelta(hours=24)).values
        dv = later - w.values
        out.append({"sid": sid, "season": int(s.SEASON.iloc[0]), "basin": g.BASIN, "lat_g": g.LAT, "lon_g": g.LON,
                    "time_g": g.t, "y_lmi": float(s.USA_WIND.max()),
                    "y_ri": float(np.nanmax(dv)) if np.isfinite(dv).any() else np.nan})
    return pd.DataFrame(out)


def add_sst(df):
    ds = xr.open_dataset(SST)
    sst = ds.sst.load(); lat = np.radians(sst.lat.values); lon = np.radians(sst.lon.values)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    tidx = pd.DatetimeIndex(sst.time.values)
    vals = []
    for r in df.itertuples():
        k = tidx.get_indexer([pd.Timestamp(r.time_g.year, r.time_g.month, 1)])[0]
        if k < 0:
            vals.append(np.nan); continue
        f = sst.values[k]
        la, lo = np.radians(r.lat_g), np.radians(r.lon_g % 360)
        dist = np.arccos(np.clip(np.sin(la) * np.sin(LA) + np.cos(la) * np.cos(LA) * np.cos(LO - lo), -1, 1))
        dist = np.where(np.isfinite(f), dist, np.inf)
        vals.append(float(f.flat[np.argmin(dist)]))
    df["sst_g"] = vals
    df["month_g"] = df.time_g.dt.month
    return df


def main():
    fetch(URL_IBT, IBT); fetch(URL_SST, SST)
    df = add_sst(storms())
    df["time_g"] = df.time_g.astype("datetime64[ns]")
    df.to_parquet(OUT)
    print(f"{len(df)} storms, {df.season.nunique()} seasons; y_ri finite {df.y_ri.notna().sum()}; "
          f"sst_g finite {df.sst_g.notna().sum()} -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
