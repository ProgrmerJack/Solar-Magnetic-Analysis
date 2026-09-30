#!/usr/bin/env python3
"""
era5_nam_daily_mean_check.py
============================
Does the January 2019 Karpechko label at 1000 hPa survive a TRUE daily mean?

`03_data_ingestion/era5_nam_daily.parquet` builds NAM_1000 from the 00 UTC field
of each day. Under the Karpechko et al. (2017) criterion the January 2019 SSW
(central date 2019-01-02) is non-downward at 1000 hPa only because its mean NAM
over days +8..+52 is +0.019 sigma rather than negative (Fig. 4c), so a
construction choice of this size could flip it.

Rebuilding the whole record from daily means would read every chunk of the
WeatherBench2 store (~89 GB compressed). The label needs far less. The daily-mean
NAM differs from the 00 UTC one by

    NAM_dm - NAM_00 = -(d - clim_d(doy)) / sd(doy),   d = h_dailymean - h_00UTC

so only d in the 2019 window and its day-of-year climatology are needed. d is
measured over the same calendar window (10 Jan - 23 Feb) in 2019 and in
CLIM_YEARS; its climatology is their mean, with its standard error carried into
the result. sd(doy) is the 00 UTC product's own; the ratio of daily-mean to
00 UTC day-to-day spread over the sampled winters is reported to show that
substituting it cannot change a sign.

Output: results/current/9_literature/era5_nam_daily_mean_check.json
"""
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "9_literature"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "03_data_ingestion"))
import acquire_era5_nam as A                        # noqa: E402

NAM = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "era5_nam_daily.parquet"
ONSET = pd.Timestamp("2019-01-02")
WIN = (8, 52)                        # Karpechko et al. 2017 classification window
LEVEL = 1000
CLIM_YEARS = [y for y in range(2009, 2023) if y != 2019]
WORKERS = 6


def window(year):
    lo = pd.Timestamp(year=year, month=ONSET.month, day=ONSET.day) + pd.Timedelta(days=WIN[0])
    hi = pd.Timestamp(year=year, month=ONSET.month, day=ONSET.day) + pd.Timedelta(days=WIN[1])
    return lo, hi


def cap_6h(year):
    """6-hourly cos-lat polar-cap geopotential height (m) at LEVEL over the window,
    reduced exactly as acquire_era5_nam reduces it."""
    lo, hi = window(year)
    ds = xr.open_zarr(A.STORE, storage_options={"token": "anon"})
    z = ds["geopotential"].sel(level=LEVEL,
                               time=slice(lo, hi + pd.Timedelta(hours=18)))
    z = z.sel(latitude=z.latitude[z.latitude >= A.CAP_LAT])
    w = np.cos(np.deg2rad(z.latitude))
    s = (z.weighted(w).mean(dim=["latitude", "longitude"]) / A.G0).compute()
    return pd.Series(s.values, index=pd.to_datetime(s.time.values))


def daily(s6):
    day = s6.index.floor("D")
    n = s6.groupby(day).size()
    if not (n == 4).all():
        raise ValueError(f"days without four 6-hourly steps: {list(n[n != 4].index)}")
    dm = s6.groupby(day).mean()
    h00 = s6[s6.index.hour == 0]
    h00.index = h00.index.floor("D")
    return dm, h00


def main():
    nam = pd.read_parquet(NAM)
    h = nam[f"z{LEVEL}_m"]
    doy = h.index.dayofyear
    mu = h.groupby(doy).transform("mean")
    sd = h.groupby(doy).transform("std")
    rebuilt = -((h - mu) / sd)
    err = float((rebuilt - nam[f"nam_{LEVEL}"]).abs().max())
    print(f"00 UTC NAM_{LEVEL} rebuilt from its own heights: max|diff| = {err:.1e}")
    assert err < 1e-9

    years = [2019] + CLIM_YEARS
    with ThreadPoolExecutor(WORKERS) as ex:
        s6 = dict(zip(years, ex.map(cap_6h, years)))
    d, ratio = {}, []
    for y in years:
        dm, h00 = daily(s6[y])
        lo, hi = window(y)
        # the stored 00 UTC heights must be exactly what this reduction gives
        ref = h.reindex(h00.index)
        agree = float((ref - h00).abs().max())
        if agree > 1e-3:
            raise ValueError(f"{y}: 00 UTC cap height disagrees with the stored "
                             f"product by {agree:.4f} m")
        d[y] = (dm - h00).rename(y)
        ratio.append(float(dm.diff().std() / h00.diff().std()))
        print(f"  {y}: {len(dm)} days, mean d = {d[y].mean():+.3f} m, "
              f"00 UTC agreement {agree:.1e} m")

    # climatology of d by day of window (calendar day), from the other winters
    clim = pd.concat([d[y].reset_index(drop=True) for y in CLIM_YEARS], axis=1)
    dclim = clim.mean(axis=1)
    d19 = d[2019].reset_index(drop=True)
    idx = d[2019].index
    sd19 = sd.reindex(idx).values
    corr = -(d19.values - dclim.values) / sd19            # NAM_dm - NAM_00, per day
    nam00 = nam[f"nam_{LEVEL}"].reindex(idx).values
    namdm = nam00 + corr

    # standard error of the window-mean correction from the climatology estimate
    per_year_win = clim.div(sd19, axis=0).mean(axis=0)     # window-mean d/sd, per winter
    se = float(per_year_win.std(ddof=1) / np.sqrt(len(CLIM_YEARS)))

    m00, mdm = float(np.mean(nam00)), float(np.mean(namdm))
    f00, fdm = float(np.mean(nam00 < 0)), float(np.mean(namdm < 0))
    dw00 = bool(m00 < 0 and f00 > 0.5)
    dwdm = bool(mdm < 0 and fdm > 0.5)
    print(f"\nJanuary 2019, NAM_{LEVEL}, days +{WIN[0]}..+{WIN[1]} ({len(idx)} days):")
    print(f"  00 UTC      mean {m00:+.4f} sigma, fraction negative {f00:.3f} -> "
          f"{'DW' if dw00 else 'NDW'}")
    print(f"  daily mean  mean {mdm:+.4f} sigma (+-{se:.4f} from the climatology "
          f"estimate), fraction negative {fdm:.3f} -> {'DW' if dwdm else 'NDW'}")
    print(f"  spread ratio daily-mean / 00 UTC day-to-day change: "
          f"{np.mean(ratio):.3f} (range {min(ratio):.3f}-{max(ratio):.3f})")
    same = dw00 == dwdm
    print("label UNCHANGED" if same else "label CHANGES")

    res = {"question": "does the January 2019 Karpechko 1000 hPa label survive a "
                       "true daily-mean NAM?",
           "onset": str(ONSET.date()), "window_days": list(WIN), "level_hPa": LEVEL,
           "n_days": int(len(idx)),
           "inputs": {"nam_00utc": str(NAM.relative_to(ROOT)), "store": A.STORE,
                      "clim_years": CLIM_YEARS},
           "check_00utc_rebuild_max_abs": err,
           "nam_00utc": {"window_mean": round(m00, 4), "frac_negative": round(f00, 3),
                         "DW": dw00},
           "nam_daily_mean": {"window_mean": round(mdm, 4),
                              "window_mean_SE_from_climatology": round(se, 4),
                              "frac_negative": round(fdm, 3), "DW": dwdm},
           "window_mean_correction": round(mdm - m00, 4),
           "day_to_day_spread_ratio_dm_over_00": round(float(np.mean(ratio)), 3),
           "label_unchanged": same}
    (RESULTS / "era5_nam_daily_mean_check.json").write_text(
        json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> era5_nam_daily_mean_check.json")


if __name__ == "__main__":
    main()
