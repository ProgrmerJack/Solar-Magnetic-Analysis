#!/usr/bin/env python3
"""
era5_regional_test.py
=====================
OBSERVED SSWs: IS THE REGIONAL COLD AFTER AN SSW WHAT ITS NAM SHIFT IMPLIES?

Plan approved 2026-09-26 (regional cold-risk plan). Written before the ERA5
regional temperature series was analysed. The observational counterpart of
snapsi_regional_test.py (T2, T3, T4 there).

DATA
  Events: load_catalogue("primary"), onsets whose days +8..+25 lie inside the
  ERA5 series (1959-01-01 .. 2023-01-10). Circulation: ERA5 NAM at 1000 hPa
  (era5_nam_daily.parquet, polar-cap form, standardised by day of year), mean
  over post-onset days +8..+25. Temperature: era5_t2m_regions_daily.parquet,
  anomaly from a day-of-year climatology (1959-2022, 31-day running mean of the
  daily means), mean over days +8..+24. Regions as in snapsi_regional_test.py;
  NEURASIA primary.

REFERENCE RELATION (SSW-free)
  Pseudo-onsets: every Nov-Mar day more than 135 days from every catalogued onset
  (the zone-free rule of the project), days +8..+25 covered. On them, OLS
  T = a + b N. It is the ordinary winter relation between the polar-cap
  circulation and regional temperature, with no SSW in the window.

TESTS
  S_T   mean T at events (the observed regional change after SSWs).
  R     mean over events of T - (a + b N): the part of the post-SSW regional
        anomaly NOT implied by the event's own NAM through the SSW-free relation.
        Null: 10,000 draws of one zone-free pseudo-onset per event within +-21
        calendar days of the event's date (any year), mean residual; two-sided
        p = P(|null| >= |R|). Falsifier (as in SNAPSI T2): p < 0.05 AND |R| > 0.2
        sigma of the SSW-free T in NEURASIA -> SSWs carry regional cold beyond
        their NAM shift in observations.
  LABEL residual regressed on the downward indicator (Karpechko conditions 1-2 on
        the daily NAM 1000 hPa over days +8..+25); coefficient with a 10,000-draw
        bootstrap over events. H0: 0.
  COLD  observed share of events with T below the SSW-free 10th percentile, and
        the share the relation predicts (Gaussian, mean a + b N, residual s.d.).

Output: results/current/2_event_study/era5_regional_test.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RESULTS = ROOT / "results" / "current" / "2_event_study"
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
from build_catalogue import load_catalogue          # noqa: E402

NAME = "era5_regional_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_NULL = 10000
N_BOOT = 10000
REGIONS = ["NEURASIA", "HI_EUROPE", "MID_EASIA", "MID_NAMER"]
PRIMARY = "NEURASIA"
NAM_WIN, T_WIN = (8, 25), (8, 24)
ZONE_SEP, CAL_WIN = 135, 21
NAM_FILE = ING / "era5_nam_daily.parquet"
T_FILE = ING / "era5_t2m_regions_daily.parquet"


def windows(nam, tan, anchors):
    """N (window mean), DW flag and regional T window means for each anchor."""
    rows = {}
    for a in anchors:
        nd = nam.reindex(pd.date_range(a + pd.Timedelta(days=NAM_WIN[0]),
                                       a + pd.Timedelta(days=NAM_WIN[1])))
        td = tan.reindex(pd.date_range(a + pd.Timedelta(days=T_WIN[0]),
                                       a + pd.Timedelta(days=T_WIN[1])))
        if nd.isna().any() or td.isna().any().any():
            continue
        rows[a] = {"N": float(nd.mean()), "DW": float((nd.mean() < 0) and ((nd < 0).mean() > 0.5)),
                   **{r: float(td[r].mean()) for r in REGIONS}}
    return pd.DataFrame.from_dict(rows, orient="index")


def main():
    rng = np.random.default_rng(SEED)
    nam = pd.read_parquet(NAM_FILE)["nam_1000"]
    t = pd.read_parquet(T_FILE).set_index("date")[REGIONS]
    t.index = pd.to_datetime(t.index)
    # day-of-year climatology, 1959-2022, 31-day running mean (circular)
    base = t[(t.index.year >= 1959) & (t.index.year <= 2022)]
    doy = base.groupby(base.index.dayofyear).mean().reindex(range(1, 367))
    doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True,
                                                                   min_periods=15).mean().iloc[15:-15]
    doy.index = range(1, 367)
    tan = t - doy.reindex(t.index.dayofyear).values

    cat = load_catalogue("primary")
    ev = windows(nam, tan, [pd.Timestamp(o) for o in cat])
    days = pd.date_range(nam.index.min(), nam.index.max())
    free = [d for d in days if d.month in (11, 12, 1, 2, 3)
            and np.all(np.abs((cat - d).days) > ZONE_SEP)]
    ps = windows(nam, tan, free)
    print(f"events used {len(ev)} of {len(cat)}; zone-free pseudo-onsets {len(ps)}", flush=True)

    res = {"plan_approved": "2026-09-26", "seed": SEED, "n_null": N_NULL, "n_boot": N_BOOT,
           "primary": PRIMARY, "nam_window": list(NAM_WIN), "t_window": list(T_WIN),
           "inputs": [str(NAM_FILE.relative_to(ROOT)), str(T_FILE.relative_to(ROOT))],
           "events": [str(d.date()) for d in ev.index], "n_events": int(len(ev)),
           "n_pseudo": int(len(ps)), "regions": {}}
    # candidate pseudo-onsets per event: +-21 calendar days, any year
    doy_ps = np.array([d.dayofyear for d in ps.index])
    cand = {}
    for e in ev.index:
        dd = np.abs(doy_ps - e.dayofyear)
        dd = np.minimum(dd, 365 - dd)
        cand[e] = np.where(dd <= CAL_WIN)[0]
    for r in REGIONS:
        b, a = np.polyfit(ps["N"], ps[r], 1)
        res_ps = ps[r] - (a + b * ps["N"])
        s_free, s_res = float(ps[r].std(ddof=1)), float(res_ps.std(ddof=2))
        resid = ev[r] - (a + b * ev["N"])
        R = float(resid.mean())
        null = np.array([np.mean([res_ps.iloc[rng.choice(cand[e])] for e in ev.index])
                         for _ in range(N_NULL)])
        p = float(np.mean(np.abs(null) >= abs(R)))
        # label: residual on DW indicator, bootstrap over events
        X = np.column_stack([np.ones(len(ev)), ev["DW"].values])
        coef = float(np.linalg.lstsq(X, resid.values, rcond=None)[0][1])
        bs = []
        for _ in range(N_BOOT):
            i = rng.integers(0, len(ev), len(ev))
            if 0 < ev["DW"].values[i].sum() < len(i):
                bs.append(np.linalg.lstsq(X[i], resid.values[i], rcond=None)[0][1])
        q10 = float(np.quantile(ps[r], 0.10))
        out = {"slope_T_on_N": round(float(b), 4), "intercept": round(float(a), 4),
               "sd_T_event_free": round(s_free, 4), "sd_residual": round(s_res, 4),
               "S_T_mean": round(float(ev[r].mean()), 4), "N_mean_events": round(float(ev["N"].mean()), 4),
               "R_mean": round(R, 4), "R_in_sd_units": round(R / s_free, 4),
               "R_null_q025_q975": [round(float(q), 4) for q in np.quantile(null, [0.025, 0.975])],
               "p_R_two_sided": round(p, 4),
               "share_of_S_T_explained": round(1 - R / float(ev[r].mean()), 3) if ev[r].mean() else None,
               "label_coef": round(coef, 4),
               "label_coef_ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
               "n_DW": int(ev["DW"].sum()),
               "cold_obs": round(float((ev[r] < q10).mean()), 4),
               "cold_pred": round(float(np.mean(norm.cdf((q10 - (a + b * ev["N"])) / s_res))), 4),
               "cold_event_free": round(float((ps[r] < q10).mean()), 4)}
        res["regions"][r] = out
        print(f"{r:10s} S_T {out['S_T_mean']:+.2f} K  R {R:+.2f} K ({out['R_in_sd_units']:+.2f} sd) "
              f"p {p:.3f}  explained {out['share_of_S_T_explained']}  label {coef:+.2f} "
              f"{out['label_coef_ci95']}  cold obs {out['cold_obs']:.2f} pred {out['cold_pred']:.2f}",
              flush=True)
    pr = res["regions"][PRIMARY]
    res["falsifier_primary"] = {"rule": "p_R < 0.05 AND |R| > 0.2 sd in NEURASIA",
                                "met": bool(pr["p_R_two_sided"] < 0.05 and abs(pr["R_in_sd_units"]) > 0.2)}
    print("falsifier met:", res["falsifier_primary"]["met"])
    (RESULTS / "era5_regional_test.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")


if __name__ == "__main__":
    main()
