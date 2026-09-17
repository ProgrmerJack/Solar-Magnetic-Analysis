#!/usr/bin/env python3
"""
r106_canonical_catalog_headline.py
==================================
Re-test the manuscript's HEADLINE result against the canonical SSW catalog.

The manuscript's event list (data/results/ssw_event_catalog.csv, 16 events) was
found to disagree with the authoritative NOAA CSL / Butler et al. (2017)
compendium (data/processed/atmospheric/ssw_canonical.csv):

  * 2012-01-11 is in the manuscript but is NOT a major SSW in ANY of the six
    reanalyses in the compendium.
  * 2000-03-20 (MAR 2000) and 2010-03-24 (MAR 2010) ARE major SSWs inside the
    Davos record era but are absent from the manuscript.

So the headline (gmRR = 0.32, 14/16 events, P = 0.004) rests on a wrong event
set. This script recomputes it under, for every catalog x design combination:

  CATALOGS   manuscript(16) | canonical(17) | canonical minus March events(15)
  DESIGNS    between-winter (the manuscript's: SSW windows vs day-of-year-matched
             days drawn from non-SSW winters)
             within-winter  (winter fixed effects + DOY harmonics; each winter is
             its own control -- the causally correct contrast)

Inference is by winter-block bootstrap throughout, because avalanche counts
cluster massively within winters (one cycle produces many events) and model
standard errors are anti-conservative by roughly an order of magnitude here.

Output: data/results/r106_canonical_catalog_headline.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r106_canonical_catalog_headline.json"
SEASON = (11, 12, 1, 2, 3, 4)
HALF = 15          # manuscript's +/-15 d window
N_BOOT = 4000


# ----------------------------------------------------------------- catalogs --
def catalogs(lo, hi):
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    cd = pd.to_datetime(can["date"]).dt.tz_localize("UTC")
    cd = cd[(cd >= lo) & (cd <= hi)]
    paper = pd.to_datetime(
        pd.read_csv(ROOT / "data/results/ssw_event_catalog.csv")["date"]).dt.tz_localize("UTC")
    return {
        "manuscript_16": pd.DatetimeIndex(sorted(paper)),
        "canonical_17": pd.DatetimeIndex(sorted(cd)),
        "canonical_no_march_15": pd.DatetimeIndex(sorted(cd[cd.dt.month != 3])),
    }


# ------------------------------------------------------------------- series --
def davos_counts():
    p = pd.read_parquet(ROOT / "data/processed/analysis_panel.parquet")
    if p.index.tz is not None:
        p.index = p.index.tz_convert("UTC").tz_localize(None)
    p.index = p.index.tz_localize("UTC")
    s = p["dry_natural_size_1234"].dropna().astype(float)
    return s[np.isin(s.index.month, SEASON)]


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


# ------------------------------------------------- design A: between-winter --
def between_winter(s, ev, half=HALF, n_boot=N_BOOT, seed=0):
    """Manuscript design: SSW windows vs DOY-matched days in NON-SSW winters."""
    sw = pd.Series(winter_of(s.index), index=s.index)
    ssw_w = set(winter_of(pd.DatetimeIndex(ev)))
    ctrl_mask = ~sw.isin(ssw_w)
    rows = []
    for o in ev:
        m = (s.index >= o - pd.Timedelta(days=half)) & (s.index <= o + pd.Timedelta(days=half))
        if m.sum() < 10:
            continue
        doys = set(s.index[m].dayofyear)
        ctrl = s[np.isin(s.index.dayofyear, list(doys)) & ctrl_mask.values]
        if len(ctrl) < 10:
            continue
        # +0.05 guard: a zero-count window otherwise sends log(RR) to -inf
        rows.append({"onset": str(o.date()),
                     "ssw_mean": float(s[m].mean()), "ctrl_mean": float(ctrl.mean()),
                     "rr": float((s[m].mean() + .05) / (ctrl.mean() + .05))})
    rr = np.array([r["rr"] for r in rows])
    n = len(rr)
    n_down = int((rr < 1).sum())
    # bootstrap over events (the manuscript's unit of inference)
    rng = np.random.default_rng(seed)
    bs = np.array([np.exp(np.mean(np.log(rr[rng.integers(0, n, n)]))) for _ in range(n_boot)])
    return {
        "design": "between-winter (manuscript)",
        "n_events": n,
        "gmRR": round(float(np.exp(np.mean(np.log(rr)))), 4),
        "gmRR_CI95": [round(float(np.percentile(bs, 2.5)), 4),
                      round(float(np.percentile(bs, 97.5)), 4)],
        "events_suppressed": f"{n_down}/{n}",
        "sign_test_p_one_sided": float(
            stats.binomtest(n_down, n, 0.5, alternative="greater").pvalue),
        "n_control_winters": int((~pd.Series(sorted(set(winter_of(s.index)))).isin(ssw_w)).sum()),
        "per_event": rows,
    }


# -------------------------------------------------- design B: within-winter --
def _frame(s, ev, a, b):
    df = pd.DataFrame({"count": s.values}, index=s.index)
    df["winter"] = winter_of(df.index)
    doy = df.index.dayofyear.values
    for k in (1, 2):
        df[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        df[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
    f = np.zeros(len(df), int)
    for o in ev:
        m = (df.index >= o + pd.Timedelta(days=a)) & (df.index <= o + pd.Timedelta(days=b))
        f[m.nonzero()[0]] = 1
    df["W"] = f
    return df.groupby("winter").filter(
        lambda g: g["W"].nunique() == 2 and g["count"].sum() > 0)


def _blocks(df):
    core = ["W", "s1", "c1", "s2", "c2"]
    return [(g["count"].values.astype(float), g[core].values.astype(float),
             g["W"].nunique() == 2) for _, g in df.groupby("winter")]


def _fit(blocks):
    keep = [b for b in blocks if b[2]]
    if len(keep) < 5:
        return None
    n = sum(len(b[0]) for b in keep)
    D = np.zeros((n, len(keep)))
    ys, Xs, i = [], [], 0
    for j, (y, X, _) in enumerate(keep):
        ys.append(y)
        Xs.append(X)
        D[i:i + len(y), j] = 1.0
        i += len(y)
    try:
        m = sm.GLM(np.concatenate(ys), np.hstack([np.vstack(Xs), D]),
                   family=sm.families.Poisson()).fit()
        return float(m.params[0])
    except Exception:
        return None


def within_winter(s, ev, a=-HALF, b=HALF, n_boot=N_BOOT, seed=0):
    df = _frame(s, ev, a, b)
    blocks = _blocks(df)
    est = _fit(blocks)
    if est is None:
        return {"design": "within-winter", "status": "insufficient"}
    rng = np.random.default_rng(seed)
    nb = len(blocks)
    bs = [_fit([blocks[k] for k in rng.integers(0, nb, nb)]) for _ in range(n_boot)]
    bs = np.array([v for v in bs if v is not None and np.isfinite(v)])
    return {
        "design": f"within-winter (winter FE), window [{a:+d},{b:+d}] d",
        "n_winters": int(df["winter"].nunique()),
        "IRR": round(float(np.exp(est)), 4),
        "IRR_CI95": [round(float(np.exp(np.percentile(bs, 2.5))), 4),
                     round(float(np.exp(np.percentile(bs, 97.5))), 4)],
        "p_one_sided_suppression": float((bs >= 0).mean()),
        "n_boot": int(len(bs)),
    }


def main():
    s = davos_counts()
    lo, hi = s.index.min(), s.index.max()
    print(f"Davos natural dry-slab counts: {len(s)} winter days "
          f"{lo.date()} .. {hi.date()}")
    cats = catalogs(lo, hi)
    res = {"record": {"start": str(lo.date()), "end": str(hi.date()), "n_days": len(s)},
           "catalogs": {k: [str(d.date()) for d in v] for k, v in cats.items()},
           "results": {}}

    for name, ev in cats.items():
        print(f"\n=== catalog: {name}  (n={len(ev)}) ===")
        A = between_winter(s, ev)
        B = within_winter(s, ev)
        print(f"  between-winter : gmRR={A['gmRR']:.3f} CI{A['gmRR_CI95']}  "
              f"{A['events_suppressed']}  sign P={A['sign_test_p_one_sided']:.4f}  "
              f"(control winters={A['n_control_winters']})")
        if "IRR" in B:
            print(f"  within-winter  : IRR ={B['IRR']:.3f} CI{B['IRR_CI95']}  "
                  f"P(suppression)={B['p_one_sided_suppression']:.3f}")
        res["results"][name] = {"between_winter": A, "within_winter": B}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
