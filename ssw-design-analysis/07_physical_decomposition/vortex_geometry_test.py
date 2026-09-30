#!/usr/bin/env python3
"""
vortex_geometry_test.py
=======================
DO SPLIT AND DISPLACEMENT SSWs DIFFER IN THEIR SURFACE RESPONSE (OBSERVATIONS)?

Plan approved 2026-09-29 (referee: SNAPSI nudges only the zonal mean, so the most
often proposed event difference -- vortex geometry -- is untested; test it in
reanalysis, even at low power). Written before any classification was made.

GEOMETRY (Seviour, Mitchell & Gray 2013, GRL 40, 5268-5273, read 2026-09-29)
  NCEP-NCAR daily 10 hPa geopotential height, 20-90N (acquire_ncep_z10_fields.py).
  Vortex edge: one contour, the December-March zonal-mean height at 60N
  (climatology 1958-2024). Two-dimensional moments of (edge - Z) over the region
  Z < edge, on a polar-stereographic plane with spherical area weights:
  centroid latitude and aspect ratio (Matthewman et al. 2009 formulae).
  VALIDATION: Seviour et al.'s own event rule -- displaced if the centroid stays
  equatorward of 66N for >= 7 days, split if the aspect ratio stays above 2.4
  for >= 7 days, onsets >= 30 days apart, December-March -- applied to 1958-2009
  must give counts near theirs (17 displaced, 18 split in 52 winters with ERA).

CLASSIFICATION OF THE CATALOGUED SSWs (primary catalogue, 1958-2024)
  Window -10..+10 days around the central date. SPLIT if the aspect ratio exceeds
  2.4 on any day; otherwise DISPLACED if the centroid latitude falls below 66N on
  any day; otherwise UNCLASSIFIED. Sensitivity: the full 7-day persistence rule.

TESTS (observations; ERA5 1000 hPa NAM, era5_nam_daily.parquet, days +8..+25 and
+8..+52; ERA5 northern-Eurasian 2 m temperature anomaly, days +8..+24, as in
era5_regional_test.py; the downward label by conditions 1-2 on NAM 1000 hPa days
+8..+25)
  G1 difference of means, split minus displaced, for each outcome: 95% interval
     from 10,000 bootstrap resamples within class; two-sided permutation p
     (10,000 relabellings).
  G2 downward rate by class, Fisher exact test.
  Reading: a significant G1 (splits more negative) is a real, non-zonal event
  difference that SNAPSI cannot test; it is reported as such.

CORRECTIONS BEFORE THE FIRST RUN (2026-09-29, code review; the test had not been
run on the complete record): (i) area weights are PLANAR on the stereographic
plane, as the paper's eq. (1) integrates dx dy (the draft used spherical cos-lat
weights); (ii) the validation period is the paper's 52 winters, January 1958 to
March 2009 (the draft covered 51); (iii) the script refuses an incomplete Z10
cache (every year 1958-2024 must be present).

Output: results/current/2_event_study/vortex_geometry_test.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RESULTS = ROOT / "results" / "current" / "2_event_study"
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE))
from build_catalogue import load_catalogue          # noqa: E402
import era5_regional_test as ER                     # noqa: E402

NAME = "vortex_geometry_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N = 10000
CACHE = ING / "_ncep_z10"
LAT_THR, AR_THR, PERSIST, SEP = 66.0, 2.4, 7, 30
WIN = (-10, 10)


def load_z():
    have = sorted(int(f.stem.split("_")[1]) for f in CACHE.glob("z10_*.npz"))
    if have != list(range(1958, 2025)):
        raise SystemExit(f"Z10 cache incomplete: {len(have)} of 67 years")
    zs, ts = [], []
    for f in sorted(CACHE.glob("z10_*.npz")):
        with np.load(f) as d:
            zs.append(d["z"]); ts.append(pd.to_datetime(d["time"])); lat, lon = d["lat"], d["lon"]
    z = np.concatenate(zs); t = pd.DatetimeIndex(np.concatenate([x.values for x in ts]))
    o = np.argsort(t); return z[o], t[o], lat, lon


def moments(z, t, lat, lon, edge):
    """Daily centroid latitude (deg) and aspect ratio on a polar-stereographic plane."""
    LAT, LON = np.meshgrid(lat, lon, indexing="ij")
    colat = np.deg2rad(90.0 - LAT)
    r = np.tan(colat / 2.0)                                 # stereographic radius (unit sphere / 2)
    x, y = r * np.cos(np.deg2rad(LON)), r * np.sin(np.deg2rad(LON))
    # planar area on the stereographic plane (the paper integrates dx dy):
    # dA_plane / dA_sphere = 1 / (4 cos^4(colat/2)) for r = tan(colat/2)
    dA = np.cos(np.deg2rad(LAT)) / (4.0 * np.cos(colat / 2.0) ** 4)
    cl, ar = np.full(len(t), np.nan), np.full(len(t), np.nan)
    for k in range(len(t)):
        q = (edge - z[k]) * dA
        q = np.where(z[k] < edge, q, 0.0)
        m00 = q.sum()
        if m00 <= 0:
            continue
        xb, yb = (q * x).sum() / m00, (q * y).sum() / m00
        j20 = (q * (x - xb) ** 2).sum() / m00; j02 = (q * (y - yb) ** 2).sum() / m00
        j11 = (q * (x - xb) * (y - yb)).sum() / m00
        s = np.sqrt(4 * j11 ** 2 + (j20 - j02) ** 2)
        ar[k] = np.sqrt((j20 + j02 + s) / max(j20 + j02 - s, 1e-12))
        cl[k] = 90.0 - np.rad2deg(2.0 * np.arctan(np.hypot(xb, yb)))
    return pd.Series(cl, index=t), pd.Series(ar, index=t)


def seviour_events(cl, ar, years):
    djfm = cl.index.month.isin([12, 1, 2, 3])
    out = {"displaced": [], "split": []}
    for kind, cond in (("displaced", (cl < LAT_THR) & djfm), ("split", (ar > AR_THR) & djfm)):
        c = cond.values; last = None
        for i in range(len(c) - PERSIST + 1):
            if c[i] and (i == 0 or not c[i - 1]) and c[i:i + PERSIST].all():
                d = cond.index[i]
                if (d.month <= 3 and d.year in years) or (d.month == 12 and d.year + 1 in years):
                    if last is None or (d - last).days >= SEP:
                        out[kind].append(d); last = d
    return out


def main():
    rng = np.random.default_rng(SEED)
    z, t, lat, lon = load_z()
    i60 = int(np.argmin(np.abs(lat - 60.0)))
    djfm = t.month.isin([12, 1, 2, 3])
    edge = float(z[djfm][:, i60, :].mean())
    cl, ar = moments(z, t, lat, lon, edge)
    val = seviour_events(cl, ar, range(1958, 2010))       # Jan 1958 .. Mar 2009: 52 winters
    res = {"plan_approved": "2026-09-29", "seed": SEED, "n": N, "edge_m": round(edge, 1),
           "validation_1958_2009": {k: len(v) for k, v in val.items()},
           "validation_reference_seviour2013": {"displaced": 17, "split": 18, "winters": 52}}
    print(f"edge {edge:.0f} m; Seviour-rule events 1958-2009: {res['validation_1958_2009']} (paper: 17 D, 18 S)")
    cat = load_catalogue("primary")
    cls, cls7 = {}, {}
    for o in cat:
        w = (cl.index >= o + pd.Timedelta(days=WIN[0])) & (cl.index <= o + pd.Timedelta(days=WIN[1]))
        if w.sum() < 15:
            continue
        a, c = ar[w], cl[w]
        cls[o] = "split" if (a > AR_THR).any() else ("displaced" if (c < LAT_THR).any() else "unclassified")
        run_a = (a > AR_THR).astype(int).groupby(((a > AR_THR) != (a > AR_THR).shift()).cumsum()).transform("sum")
        run_c = (c < LAT_THR).astype(int).groupby(((c < LAT_THR) != (c < LAT_THR).shift()).cumsum()).transform("sum")
        cls7[o] = ("split" if ((a > AR_THR) & (run_a >= PERSIST)).any() else
                   "displaced" if ((c < LAT_THR) & (run_c >= PERSIST)).any() else "unclassified")
    res["classification"] = {str(o.date()): v for o, v in cls.items()}
    res["classification_7day"] = {str(o.date()): v for o, v in cls7.items()}
    # outcomes
    nam = pd.read_parquet(ER.NAM_FILE)["nam_1000"]
    tt = pd.read_parquet(ER.T_FILE).set_index("date")["NEURASIA"]; tt.index = pd.to_datetime(tt.index)
    base = tt[(tt.index.year >= 1959) & (tt.index.year <= 2022)]
    doy = base.groupby(base.index.dayofyear).mean().reindex(range(1, 367))
    doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    doy.index = range(1, 367)
    tan = tt - doy.reindex(tt.index.dayofyear).values

    def wm(s, o, a, b):
        v = s.reindex(pd.date_range(o + pd.Timedelta(days=a), o + pd.Timedelta(days=b)))
        return float(v.mean()) if v.notna().all() else np.nan

    rows = []
    for o, k in cls.items():
        n25 = nam.reindex(pd.date_range(o + pd.Timedelta(days=8), o + pd.Timedelta(days=25)))
        rows.append({"date": o, "cls": k, "cls7": cls7[o], "nam_8_25": wm(nam, o, 8, 25),
                     "nam_8_52": wm(nam, o, 8, 52), "t_neur": wm(tan, o, 8, 24),
                     "dw": float(n25.notna().all() and n25.mean() < 0 and (n25 < 0).mean() > 0.5)
                     if n25.notna().all() else np.nan})
    df = pd.DataFrame(rows)
    res["counts"] = df.cls.value_counts().to_dict(); res["counts_7day"] = df.cls7.value_counts().to_dict()
    res["tests"] = {}
    for scheme in ("cls", "cls7"):
        for oc in ("nam_8_25", "nam_8_52", "t_neur"):
            d = df[df[scheme].isin(["split", "displaced"]) & df[oc].notna()]
            s, dsp = d[d[scheme] == "split"][oc].values, d[d[scheme] == "displaced"][oc].values
            if len(s) < 3 or len(dsp) < 3:
                continue
            diff = float(s.mean() - dsp.mean())
            bs = [rng.choice(s, len(s)).mean() - rng.choice(dsp, len(dsp)).mean() for _ in range(N)]
            allv = np.concatenate([s, dsp])
            perm = []
            for _ in range(N):
                p = rng.permutation(allv); perm.append(p[:len(s)].mean() - p[len(s):].mean())
            res["tests"][f"{scheme}|{oc}"] = {
                "n_split": int(len(s)), "n_displaced": int(len(dsp)),
                "mean_split": round(float(s.mean()), 3), "mean_displaced": round(float(dsp.mean()), 3),
                "diff_split_minus_displaced": round(diff, 3),
                "ci95": [round(float(q), 3) for q in np.quantile(bs, [0.025, 0.975])],
                "p_perm_two_sided": round(float(np.mean(np.abs(perm) >= abs(diff))), 4)}
        d = df[df[scheme].isin(["split", "displaced"]) & df.dw.notna()]
        tab = [[int(((d[scheme] == c) & (d.dw == 1)).sum()), int(((d[scheme] == c) & (d.dw == 0)).sum())]
               for c in ("split", "displaced")]
        res["tests"][f"{scheme}|dw_rate"] = {"table_split_displaced_x_dw_ndw": tab,
                                             "dw_rate_split": round(tab[0][0] / max(sum(tab[0]), 1), 3),
                                             "dw_rate_displaced": round(tab[1][0] / max(sum(tab[1]), 1), 3),
                                             "fisher_p": round(float(fisher_exact(tab)[1]), 4)}
    for k, v in res["tests"].items():
        print(k, v)
    (RESULTS / "vortex_geometry_test.json").write_text(json.dumps(res, indent=2, default=str),
                                                       encoding="utf8", newline="\n")
    print("Saved -> vortex_geometry_test.json")


if __name__ == "__main__":
    main()
