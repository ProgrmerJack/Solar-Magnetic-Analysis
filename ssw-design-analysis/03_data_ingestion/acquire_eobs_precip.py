#!/usr/bin/env python3
"""
acquire_eobs_precip.py
======================
Daily station-based precipitation (E-OBS v33.0e ensemble mean, 0.25 deg, 1950 to
present; Copernicus Climate Data Store "insitu-gridded-observations-europe") averaged
over three literature-defined boxes, for 08_literature_audit/criterion_precip_contrast.py
(registered before this script first ran; no E-OBS value had been examined by this
project).

BOXES (Dai et al. 2025, Weather Clim. Dynam. 6, 841, SNAPSI 2018 precipitation)
  IBERIA 36-44N, 10-1W;  SCANDINAVIA 55-70N, 4-20E;  BRITISH_ISLES 49-61N, 10W-1E
REDUCTION  cos-latitude mean over the land cells present that day; a box-day is kept
  only if at least 70% of the box's cells that ever report are present (fixed before
  any value was read); the share present is stored (COV_<box>).
The CDS key is read as in acquire_era5_u10_cds.client() and never printed.
Output: eobs_precip_regions_daily.parquet  date, IBERIA, SCANDINAVIA, BRITISH_ISLES, COV_*
        raw download kept in raw/eobs/ (read-only thereafter)
"""
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
RAW = HERE / "raw" / "eobs"
OUT = HERE / "eobs_precip_regions_daily.parquet"
BOXES = {"IBERIA": (36, 44, -10, -1), "SCANDINAVIA": (55, 70, 4, 20), "BRITISH_ISLES": (49, 61, -10, 1)}
MIN_COVER = 0.70
REQUEST = {"product_type": "ensemble_mean", "variable": ["precipitation_amount"], "grid_resolution": "0_25deg",
           "period": "full_period", "version": ["33_0e"]}


def download():
    RAW.mkdir(parents=True, exist_ok=True)
    target = RAW / "eobs_rr_v33.0e_0.25deg_full.zip"
    if not target.exists():
        import acquire_era5_u10_cds as U
        U.client().retrieve("insitu-gridded-observations-europe", REQUEST, str(target))
    ncs = sorted(RAW.glob("*.nc"))
    if not ncs:
        with zipfile.ZipFile(target) as z:
            z.extractall(RAW)
        ncs = sorted(RAW.glob("*.nc"))
    return ncs


def reduce(nc):
    import xarray as xr
    ds = xr.open_dataset(nc)
    rr = ds["rr"]
    cols = {"date": pd.to_datetime(rr["time"].values)}
    for b, (la0, la1, lo0, lo1) in BOXES.items():
        sub = rr.sel(latitude=slice(la0, la1), longitude=slice(lo0, lo1)).load()
        ok = np.isfinite(sub.values)
        ever = ok.any(axis=0)
        w = np.cos(np.deg2rad(sub.latitude.values))[:, None] * ever
        num = np.nansum(np.where(ok, sub.values, 0.0) * w[None], axis=(1, 2))
        den = (ok * w[None]).sum(axis=(1, 2))
        cover = ok[:, ever].mean(axis=1)
        val = np.where((cover >= MIN_COVER) & (den > 0), num / np.where(den > 0, den, 1), np.nan)
        if np.nanmin(val) < 0 or np.nanmax(val) > 200:
            raise ValueError(f"{b}: implausible box mean {np.nanmin(val)}..{np.nanmax(val)} mm/day")
        cols[b] = val; cols[f"COV_{b}"] = cover
        print(b, "cells", int(ever.sum()), "days kept", int(np.isfinite(val).sum()), "of", len(val), flush=True)
    return pd.DataFrame(cols)


def main():
    ncs = [p for p in download() if "rr" in p.name]
    if len(ncs) != 1:
        sys.exit(f"expected one rr file, found {[p.name for p in ncs]}")
    d = reduce(ncs[0])
    d.to_parquet(OUT)
    print(f"Saved -> {OUT.name} ({d.date.min().date()} .. {d.date.max().date()}, {len(d):,} days)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
