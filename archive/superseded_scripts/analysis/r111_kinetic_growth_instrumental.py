#!/usr/bin/env python3
"""
r111_kinetic_growth_instrumental.py
===================================
A HUMAN-FREE test of the manuscript's structural mechanism.

Every observable tested so far is mediated by people:
  * avalanche accidents      -> driven by how many people were exposed
  * bulletins / danger / problems -> forecaster judgement
  * observed avalanche counts     -> observer effort
Each carries a confound that no statistical design can fully remove.

This script tests the mechanism with INSTRUMENTS ONLY.

PHYSICS
  The manuscript's structural limb is *kinetic-growth metamorphism*: faceted
  crystals and depth hoar form when the temperature gradient through the
  snowpack is steep. The classical threshold is

        TG = (T_base - T_surface) / HS  >  10 K/m       (Colbeck 1982; Akitaya 1974)

  with >20 K/m giving strong kinetic growth / depth hoar. Under a snowpack the
  base sits near 0 degC, and the snow surface tracks air temperature, so

        TG ~ (0 - T_air) / HS

  is computable from an automatic station reporting air temperature and snow
  depth. Nothing here depends on a human observing, reporting, or judging.

  The mechanism is physically coherent both ways: an SSW-driven cold, blocked
  regime lowers T_air (raising the numerator) AND suppresses snowfall (lowering
  HS), so both terms push TG up. The scientific content is that the snowpack
  crosses a NONLINEAR structural threshold, not merely that it gets colder --
  so both the continuous gradient and the threshold-crossing rate are reported.

DATA  SNOTEL daily, 945 automatic stations, 1980-2026 (10.59M station-days):
      snwd_mm (ultrasonic snow depth), tavg_c/tmin_c/tmax_c, prec_mm, wteq_mm.
      Covers ~28 canonical SSWs -- roughly five times the bulletin-era tests.

DESIGN  Station-days are aggregated to (state, date) counts and modelled as a
  RATE with offset log(number of reporting stations):

      n_in_regime ~ W + DOY harmonics + C(state x winter),  offset log(n_report)

  conditional Poisson (stratum effects profiled out), winter-block bootstrap
  resampling whole winters across all states together. Pre-onset placebo must
  be null. Snow-climate split (continental vs maritime) is pre-specified from
  R102 and from the manuscript's continental-specificity claim.

Output: data/results/r111_kinetic_growth_instrumental.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import condpois

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r111_kinetic_growth_instrumental.json"
SEASON = (11, 12, 1, 2, 3, 4)
N_BOOT = 1000
N_HARM = 3
MIN_HS_MM = 200          # >=20 cm snow: TG is meaningless on a thin/absent pack
MIN_STATIONS = 5         # per state-day

CONTINENTAL = {"CO", "UT", "WY", "MT", "ID", "NV", "NM", "AZ", "SD"}
MARITIME = {"WA", "OR", "CA", "AK"}

WINDOWS = {
    "placebo_pre": (-45, -16),
    "onset":       (-15, 15),
    "early_post":  (0, 14),
    "mid_post":    (15, 29),
    "late_post":   (30, 44),
    "post_15_44":  (15, 44),
}

# outcome -> (description, manuscript-predicted direction)
OUTCOMES = {
    "kinetic":      ("TG > 10 K/m: kinetic-growth regime (faceting)", "up"),
    "strong_kinetic": ("TG > 20 K/m: strong kinetic growth / depth hoar", "up"),
    "melt":         ("Tmax > 0 degC: melt-driven natural trigger", "down"),
    "new_snow":     ("prec > 10 mm with Tavg < 0: loading/new-snow trigger", "down"),
    "ros":          ("prec > 5 mm with Tmax > 1 degC: rain-on-snow trigger", "down"),
}


def state_of(sid):
    p = str(sid).split(":")
    return p[1] if len(p) > 1 else "??"


def build():
    sn = pd.read_parquet(ROOT / "data/processed/cryosphere/snotel_daily.parquet",
                         columns=["station_id", "snwd_mm", "tavg_c", "tmax_c", "prec_mm"])
    if not isinstance(sn.index, pd.DatetimeIndex):
        sn = sn.set_index("date")
    if sn.index.tz is None:
        sn.index = sn.index.tz_localize("UTC")
    sn = sn[np.isin(sn.index.month, SEASON)]
    sn["state"] = sn["station_id"].map(state_of)
    sn = sn[sn["state"].isin(CONTINENTAL | MARITIME)]

    hs_m = sn["snwd_mm"] / 1000.0
    ok = (sn["snwd_mm"] >= MIN_HS_MM) & sn["tavg_c"].notna()
    tg = np.where(ok, (0.0 - sn["tavg_c"]) / hs_m.replace(0, np.nan), np.nan)
    sn["tg"] = tg
    sn["valid_tg"] = ok.astype(int)
    sn["kinetic"] = ((sn["tg"] > 10) & ok).astype(int)
    sn["strong_kinetic"] = ((sn["tg"] > 20) & ok).astype(int)
    sn["melt"] = (sn["tmax_c"] > 0).astype(int)
    sn["new_snow"] = ((sn["prec_mm"] > 10) & (sn["tavg_c"] < 0)).astype(int)
    sn["ros"] = ((sn["prec_mm"] > 5) & (sn["tmax_c"] > 1)).astype(int)

    sn = sn.rename_axis("date").reset_index()
    agg = sn.groupby(["state", "date"]).agg(
        n_report=("station_id", "nunique"),
        n_valid_tg=("valid_tg", "sum"),
        tg_mean=("tg", "mean"),
        kinetic=("kinetic", "sum"),
        strong_kinetic=("strong_kinetic", "sum"),
        melt=("melt", "sum"),
        new_snow=("new_snow", "sum"),
        ros=("ros", "sum"),
    ).reset_index()
    agg = agg[agg["n_report"] >= MIN_STATIONS].copy()
    agg["climate"] = np.where(agg["state"].isin(CONTINENTAL), "continental", "maritime")
    agg["winter"] = np.where(agg["date"].dt.month >= 11,
                             agg["date"].dt.year + 1, agg["date"].dt.year)
    doy = agg["date"].dt.dayofyear.values
    for k in range(1, N_HARM + 1):
        agg[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        agg[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
    return agg


def tag(d, ev, a, b):
    f = np.zeros(len(d), int)
    dt = d["date"].values
    for o in ev:
        m = (dt >= np.datetime64(o + pd.Timedelta(days=a))) & \
            (dt <= np.datetime64(o + pd.Timedelta(days=b)))
        f[m.nonzero()[0]] = 1
    out = d.copy()
    out["W"] = f
    return out.groupby(["state", "winter"]).filter(lambda g: g["W"].nunique() == 2)


def run(p, col, label, seed=0):
    """Rate model: successes out of reporting stations."""
    expo = p["n_valid_tg"] if col in ("kinetic", "strong_kinetic") else p["n_report"]
    m = expo.values > 0
    p = p[m]
    expo = expo[m]
    y = p[col].values.astype(float)
    strata = pd.Categorical(
        p["state"].astype(str) + "|" + p["winter"].astype(str)).codes.astype(np.int64)
    H = [p[f"{t}{k}"].values for k in range(1, N_HARM + 1) for t in ("s", "c")]
    X = np.column_stack([p["W"].values.astype(float)] + H)
    off = np.log(expo.values.astype(float))

    b = condpois.fit(y, X, strata, offset=off)
    if b is None:
        return {"window": label, "status": "unidentified"}
    est = float(b[0])

    rng = np.random.default_rng(seed)
    st = pd.Categorical(p["state"]).codes.astype(np.int64)
    n_st = int(st.max()) + 1
    winters = np.array(sorted(p["winter"].unique()))
    idx_by_w = [np.flatnonzero((p["winter"] == w).values) for w in winters]
    bs = []
    for _ in range(N_BOOT):
        pick = rng.integers(0, len(winters), len(winters))
        rows = np.concatenate([idx_by_w[k] for k in pick])
        sb = np.concatenate([st[idx_by_w[k]] + j * n_st for j, k in enumerate(pick)])
        v = condpois.fit(y[rows], X[rows], sb, offset=off[rows])
        if v is not None and np.isfinite(v[0]):
            bs.append(float(v[0]))
    bs = np.array(bs)
    return {
        "window": label,
        "RR": round(float(np.exp(est)), 4),
        "CI95": [round(float(np.exp(np.percentile(bs, 2.5))), 4),
                 round(float(np.exp(np.percentile(bs, 97.5))), 4)] if len(bs) else None,
        "p_one_sided_up": float((bs <= 0).mean()) if len(bs) else None,
        "p_one_sided_down": float((bs >= 0).mean()) if len(bs) else None,
        "base_rate": round(float(y.sum() / expo.sum()), 4),
        "n_state_winters": int(p.groupby(["state", "winter"]).ngroups),
        "n_state_days": int(len(p)), "n_boot": int(len(bs)),
    }


def main():
    d = build()
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    ev = pd.to_datetime(can["date"]).dt.tz_localize("UTC")
    ev = pd.DatetimeIndex(ev[(ev >= d["date"].min()) & (ev <= d["date"].max())])
    print(f"SNOTEL aggregated: {len(d):,} state-days, {d['state'].nunique()} states, "
          f"{d['date'].min().date()}..{d['date'].max().date()}, "
          f"{d['winter'].nunique()} winters")
    print(f"Canonical SSWs covered: {len(ev)}")
    print(f"mean TG (K/m): {d['tg_mean'].mean():.1f}; "
          f"kinetic base rate: {d['kinetic'].sum()/max(d['n_valid_tg'].sum(),1):.3f}")

    res = {"n_state_days": int(len(d)), "n_ssw": len(ev),
           "ssw_dates": [str(x.date()) for x in ev],
           "outcomes": {k: v[0] for k, v in OUTCOMES.items()},
           "min_snow_depth_mm": MIN_HS_MM, "domains": {}}

    for dom, sub in (("continental", d[d["climate"] == "continental"]),
                     ("maritime", d[d["climate"] == "maritime"])):
        print(f"\n{'='*70}\n### {dom.upper()}  ({sub['state'].nunique()} states)")
        res["domains"][dom] = {}
        for col, (desc, pred) in OUTCOMES.items():
            print(f"\n--- {col}: {desc}  [predicted {pred}] ---")
            res["domains"][dom][col] = {}
            for lab, (a, b) in WINDOWS.items():
                p = tag(sub, ev, a, b)
                if p.empty or p["W"].nunique() < 2:
                    continue
                r = run(p, col, lab)
                res["domains"][dom][col][lab] = r
                if "RR" in r:
                    print(f"  {lab:12s} RR={r['RR']:5.2f}  "
                          f"CI[{r['CI95'][0]:.2f},{r['CI95'][1]:.2f}]  "
                          f"P(up)={r['p_one_sided_up']:.3f} "
                          f"P(down)={r['p_one_sided_down']:.3f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
