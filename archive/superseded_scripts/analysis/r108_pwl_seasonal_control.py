#!/usr/bin/env python3
"""
r108_pwl_seasonal_control.py
============================
Re-test the manuscript's SECOND decisive claim -- that persistent weak layers
DEEPEN during SSW windows (Pen_depth; 15/16 events, sign P = 0.0005) -- with
the seasonal and winter controls the original analysis lacks.

WHY THIS RE-TEST IS NECESSARY
  scripts/analysis/r82_pwl_depth_validation.py compares SNOWPACK station-days
  inside SSW windows against ALL other winter station-days pooled, with:
    * no day-of-year control -- but snowpack (and therefore the depth of any
      buried weak layer) deepens monotonically through the season, and SSWs
      cluster in mid-winter while the control pool is dominated by Nov/Dec and
      Mar/Apr days. Season alone predicts the sign of the reported effect.
    * no winter control    -- deep-snow winters differ systematically.
    * no station control   -- stations differ in altitude and climatology.
    * p-values from ~10^5 station-days that are not independent
      (130 stations x 30 days per event, autocorrelated).
  So the reported effect is confounded by season by construction.

DESIGN HERE
  Pen_depth ~ W + C(station x winter) + DOY harmonics
    * C(station x winter) removes each station's each-winter snowpack level.
    * DOY harmonics remove the seasonal deepening trend -- the key confounder.
    * Inference by winter-block bootstrap (station-days are not independent).
  Reported alongside the uncontrolled contrast so the confounding is visible.

  Secondary: event-level test in which each SSW window is compared with the
  DOY-matched expectation from the SAME winter (each winter its own control).

Catalog: canonical NOAA CSL / Butler et al. (2017) compendium.

Output: data/results/r108_pwl_seasonal_control.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r108_pwl_seasonal_control.json"
SEASON = (11, 12, 1, 2, 3, 4)
HALF = 15
N_BOOT = 1500
METRICS = {
    "Pen_depth":   ("deeper during SSW", +1),
    "min_ccl_pen": ("lower during SSW (easier propagation)", -1),
    "sk38_pwl":    ("lower during SSW (less stable)", -1),
    "HS_mod":      ("snow depth - CONFOUND CHECK", +1),
}


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def load():
    sp = pd.read_csv(ROOT / "data/cryosphere/envidat/weather_snowpack_danger.csv",
                     parse_dates=["datum"],
                     usecols=["datum", "station_code", "elevation_station",
                              "Pen_depth", "min_ccl_pen", "sk38_pwl",
                              "pwl_100", "HS_mod"])
    sp = sp[sp["datum"].dt.month.isin(SEASON)].copy()
    sp["winter"] = winter_of(sp["datum"].dt)
    sp["doy"] = sp["datum"].dt.dayofyear
    for k in (1, 2):
        sp[f"s{k}"] = np.sin(2 * np.pi * k * sp["doy"] / 365.25)
        sp[f"c{k}"] = np.cos(2 * np.pi * k * sp["doy"] / 365.25)
    return sp


def tag(sp, ev, half=HALF):
    f = np.zeros(len(sp), bool)
    eid = np.full(len(sp), -1)
    d = sp["datum"].values
    for i, o in enumerate(ev):
        m = (d >= np.datetime64(o - pd.Timedelta(days=half))) & \
            (d <= np.datetime64(o + pd.Timedelta(days=half)))
        f |= m
        eid[m] = i
    out = sp.copy()
    out["W"] = f.astype(int)
    out["event_id"] = eid
    return out


def _demean(A, codes, n_strata):
    """Subtract the stratum mean from every column (Frisch-Waugh-Lovell).

    ponytail: absorbing the station x winter effects by demeaning gives exactly
    the dummy-variable OLS coefficient on W, but avoids materialising a
    285,000 x ~3,000 design matrix once per bootstrap replicate.
    """
    cnt = np.bincount(codes, minlength=n_strata).astype(float)
    out = np.empty_like(A, dtype=float)
    for j in range(A.shape[1]):
        sm_ = np.bincount(codes, weights=A[:, j], minlength=n_strata)
        out[:, j] = A[:, j] - (sm_ / cnt)[codes]
    return out


def _prep(df, col):
    """Rows usable for the FE fit + stratum codes (both W levels present)."""
    d = df[[col, "W", "s1", "c1", "s2", "c2", "station_code", "winter"]].dropna()
    if len(d) < 200:
        return None
    key = d["station_code"].astype(str) + "|" + d["winter"].astype(str)
    codes = pd.Categorical(key).codes
    n = codes.max() + 1
    # a stratum informs W only if it holds both exposed and unexposed days
    wmin = np.full(n, 2.0)
    wmax = np.full(n, -1.0)
    np.minimum.at(wmin, codes, d["W"].values)
    np.maximum.at(wmax, codes, d["W"].values)
    keep = (wmax > wmin)[codes]
    d = d[keep]
    if len(d) < 200 or d["W"].nunique() < 2:
        return None
    _, codes = np.unique(pd.Categorical(
        d["station_code"].astype(str) + "|" + d["winter"].astype(str)).codes,
        return_inverse=True)
    return d, codes, codes.max() + 1


def _fit_prepped(d, codes, n_strata, col):
    A = np.column_stack([d[col].values.astype(float),
                         d["W"].values.astype(float),
                         d[["s1", "c1", "s2", "c2"]].values.astype(float)])
    Ad = _demean(A, codes, n_strata)
    y, X = Ad[:, 0], Ad[:, 1:]
    try:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(beta[0])
    except Exception:
        return None


def _fit(df, col):
    p = _prep(df, col)
    if p is None:
        return None, 0
    d, codes, n = p
    return _fit_prepped(d, codes, n, col), len(d)


def controlled(df, col, n_boot=N_BOOT, seed=0):
    p = _prep(df, col)
    if p is None:
        return {"status": "insufficient", "n": 0}
    d, codes, n_strata = p
    est = _fit_prepped(d, codes, n_strata, col)
    if est is None:
        return {"status": "fit_failed", "n": len(d)}

    # prebuilt arrays -> bootstrap without re-running the pandas prep
    A = np.column_stack([d[col].values.astype(float),
                         d["W"].values.astype(float),
                         d[["s1", "c1", "s2", "c2"]].values.astype(float)])
    # int64: Categorical.codes is int8/int16, and code + j*n_stn overflows it
    stn = pd.Categorical(d["station_code"]).codes.astype(np.int64)
    n_stn = int(stn.max()) + 1
    winters = np.array(sorted(d["winter"].unique()))
    idx_by_w = [np.flatnonzero((d["winter"] == w).values) for w in winters]

    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(winters), len(winters))
        rows = np.concatenate([idx_by_w[k] for k in pick])
        strata = np.concatenate([stn[idx_by_w[k]] + j * n_stn
                                 for j, k in enumerate(pick)])
        _, strata = np.unique(strata, return_inverse=True)
        Ad = _demean(A[rows], strata, strata.max() + 1)
        try:
            beta, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
            if np.isfinite(beta[0]):
                bs.append(float(beta[0]))
        except Exception:
            pass
    bs = np.array(bs)
    n = len(d)
    return {
        "effect_controlled": round(float(est), 4),
        "CI95": [round(float(np.percentile(bs, 2.5)), 4),
                 round(float(np.percentile(bs, 97.5)), 4)] if len(bs) else None,
        "p_two_sided": float(2 * min((bs >= 0).mean(), (bs <= 0).mean())) if len(bs) else None,
        "n_station_days": int(n), "n_boot": int(len(bs)),
    }


def uncontrolled(df, col):
    a = df.loc[df["W"] == 1, col].dropna()
    b = df.loc[df["W"] == 0, col].dropna()
    return {
        "ssw_mean": round(float(a.mean()), 3), "ctrl_mean": round(float(b.mean()), 3),
        "diff": round(float(a.mean() - b.mean()), 3),
        "mw_p_PSEUDOREPLICATED": float(
            stats.mannwhitneyu(a, b, alternative="two-sided").pvalue),
        "n_ssw_days": int(len(a)), "n_ctrl_days": int(len(b)),
        "note": "no season/winter/station control; p-value treats station-days as independent",
    }


def event_level(df, ev, col):
    """Each SSW window vs the DOY-matched expectation from its OWN winter."""
    rows = []
    for i, o in enumerate(ev):
        w = winter_of(pd.DatetimeIndex([o]))[0]
        g = df[df["winter"] == w]
        if len(g) < 200:
            continue
        inw = g[g["event_id"] == i][col].dropna()
        if len(inw) < 10:
            continue
        doys = set(g.loc[g["event_id"] == i, "doy"])
        # same winter, other days, but DOY-matched via a seasonal fit
        rest = g[g["event_id"] != i][[col, "doy"]].dropna()
        if len(rest) < 50:
            continue
        X = np.column_stack([np.ones(len(rest)),
                             np.sin(2 * np.pi * rest["doy"] / 365.25),
                             np.cos(2 * np.pi * rest["doy"] / 365.25),
                             np.sin(4 * np.pi * rest["doy"] / 365.25),
                             np.cos(4 * np.pi * rest["doy"] / 365.25)])
        try:
            beta = sm.OLS(rest[col].values.astype(float), X).fit().params
        except Exception:
            continue
        dd = np.array(sorted(doys))
        Xp = np.column_stack([np.ones(len(dd)),
                              np.sin(2 * np.pi * dd / 365.25), np.cos(2 * np.pi * dd / 365.25),
                              np.sin(4 * np.pi * dd / 365.25), np.cos(4 * np.pi * dd / 365.25)])
        exp = float((Xp @ beta).mean())
        rows.append({"onset": str(pd.Timestamp(o).date()),
                     "observed": round(float(inw.mean()), 3),
                     "expected_same_winter_doy_matched": round(exp, 3),
                     "diff": round(float(inw.mean()) - exp, 3)})
    d = np.array([r["diff"] for r in rows])
    n_pos = int((d > 0).sum())
    return {
        "n_events": len(rows), "events_positive": f"{n_pos}/{len(rows)}",
        "mean_diff": round(float(np.mean(d)), 4) if len(d) else None,
        "sign_p_two_sided": float(stats.binomtest(n_pos, len(d), 0.5).pvalue) if len(d) else None,
        "wilcoxon_p": float(stats.wilcoxon(d).pvalue) if len(d) > 5 else None,
        "per_event": rows,
    }


def main():
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    ev = pd.to_datetime(can["date"])
    sp = load()
    lo, hi = sp["datum"].min(), sp["datum"].max()
    ev = ev[(ev >= lo) & (ev <= hi)].reset_index(drop=True)
    print(f"SNOWPACK: {len(sp):,} station-days {lo.date()}..{hi.date()}, "
          f"{sp['station_code'].nunique()} stations")
    print(f"Canonical SSWs in range: {len(ev)}")

    df = tag(sp, ev)
    res = {"n_ssw": len(ev), "n_station_days": int(len(sp)), "metrics": {}}
    for col, (expect, sign) in METRICS.items():
        if col not in df.columns:
            continue
        print(f"\n--- {col}  (expected: {expect}) ---")
        u = uncontrolled(df, col)
        c = controlled(df, col)
        e = event_level(df, ev, col)
        print(f"  uncontrolled (as in r82): diff={u['diff']:+.3f}  "
              f"P={u['mw_p_PSEUDOREPLICATED']:.2e}  [confounded by season]")
        if "effect_controlled" in c:
            print(f"  season+winter+station controlled: {c['effect_controlled']:+.3f}  "
                  f"CI{c['CI95']}  P={c['p_two_sided']:.4f}")
        print(f"  event-level (own-winter DOY-matched): {e['events_positive']}  "
              f"mean {e['mean_diff']:+.3f}  sign P={e['sign_p_two_sided']}")
        res["metrics"][col] = {"expected": expect, "uncontrolled": u,
                               "controlled": c, "event_level": e}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
