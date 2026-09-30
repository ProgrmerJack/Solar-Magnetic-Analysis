#!/usr/bin/env python3
"""
r112_design_positive_control.py
===============================
POSITIVE CONTROL for the within-winter design used throughout R106-R111.

THE PROBLEM THIS SOLVES
  R106-R111 report that every SSW effect vanishes once each winter is used as
  its own control. There are two possible explanations and they have opposite
  consequences for the manuscript:

    (A) the avalanche effect really is a between-winter confound  -> nulls stand
    (B) the within-winter design is OVER-CONSERVATIVE. SSW surface impacts last
        ~60 days (Baldwin & Dunkerton 2001), which is a third of a snow season,
        so winter fixed effects could absorb the very signal being tested. If so
        the nulls are uninformative and the manuscript's original results stand.

  These cannot be distinguished from avalanche data alone. They CAN be
  distinguished with a positive control: apply the identical design to an
  SSW impact that is beyond dispute.

THE CONTROL
  The canonical, textbook consequence of a sudden stratospheric warming is a
  shift to the negative phase of the Arctic Oscillation persisting up to ~60
  days (Baldwin & Dunkerton, Science 294, 581, 2001). If the within-winter
  design cannot recover THAT, it cannot be trusted to detect anything, and
  every null in R106-R111 must be discarded.

  Data: CPC daily AO and NAO indices, 1950-2026 (all 47 canonical SSWs);
  NCEP daily 500 hPa height / SLP / 850 hPa zonal wind, 1979-2024.

  Both designs are run on the SAME series:
    between-winter : SSW window vs day-of-year-matched days in non-SSW winters
                     (the manuscript's design, and the field's standard composite)
    within-winter  : winter fixed effects + DOY harmonics, each winter its own
                     control, winter-block bootstrap (the design used in R106-R111)

INTERPRETATION
  within-winter recovers a strong negative AO  -> design validated; the avalanche
                                                  nulls are real (explanation A)
  within-winter kills the AO signal too        -> design over-conservative;
                                                  the nulls are artefacts and the
                                                  manuscript's results are reinstated
                                                  (explanation B)

Output: data/results/r112_design_positive_control.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r112_design_positive_control.json"
SEASON = (11, 12, 1, 2, 3, 4)
N_BOOT = 4000
N_HARM = 3

WINDOWS = {
    "placebo_pre": (-45, -16),
    "onset":       (-15, 15),
    "early_post":  (0, 14),
    "mid_post":    (15, 29),
    "late_post":   (30, 44),
    "post_0_60":   (0, 60),      # Baldwin-Dunkerton persistence window
    "post_15_44":  (15, 44),
}


def read_cpc(path):
    d = pd.read_csv(path, sep=r"\s+", header=None,
                    names=["year", "month", "day", "value"])
    d["date"] = pd.to_datetime(dict(year=d.year, month=d.month, day=d.day),
                               errors="coerce")
    d = d.dropna(subset=["date"])
    return d.set_index("date")["value"].astype(float).sort_index()


def load_series():
    s = {}
    s["AO"] = read_cpc(ROOT / "data/processed/atmospheric/ao_daily_cpc.txt")
    s["NAO"] = read_cpc(ROOT / "data/processed/atmospheric/nao_daily_cpc.txt")
    t = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_troposphere.parquet")
    if t.index.tz is not None:
        t.index = t.index.tz_convert("UTC").tz_localize(None)
    for c, name in (("hgt_500hPa_m", "Z500_NH"), ("slp_Pa", "SLP_NH"),
                    ("uwnd_850hPa_ms", "U850_NH")):
        if c in t.columns:
            s[name] = t[c].dropna()
    return s


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def frame(s, ev, a, b):
    d = pd.DataFrame({"y": s.values}, index=s.index)
    d = d[np.isin(d.index.month, SEASON)]
    d["winter"] = winter_of(d.index)
    doy = d.index.dayofyear.values
    for k in range(1, N_HARM + 1):
        d[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        d[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
    f = np.zeros(len(d), int)
    for o in ev:
        m = (d.index >= o + pd.Timedelta(days=a)) & (d.index <= o + pd.Timedelta(days=b))
        f[m.nonzero()[0]] = 1
    d["W"] = f
    return d.groupby("winter").filter(lambda g: g["W"].nunique() == 2)


def _demean(A, codes, n):
    cnt = np.bincount(codes, minlength=n).astype(float)
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        sm = np.bincount(codes, weights=A[:, j], minlength=n)
        out[:, j] = A[:, j] - (sm / cnt)[codes]
    return out


def within_winter(d, n_boot=N_BOOT, seed=0):
    """Winter fixed effects + DOY harmonics; effect of W in units of the series."""
    cols = ["y", "W"] + [f"{t}{k}" for k in range(1, N_HARM + 1) for t in ("s", "c")]
    A = d[cols].values.astype(float)
    codes = pd.Categorical(d["winter"]).codes.astype(np.int64)
    _, codes = np.unique(codes, return_inverse=True)
    n = codes.max() + 1

    def fit(A_, c_, n_):
        Ad = _demean(A_, c_, n_)
        try:
            beta, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
            return float(beta[0])
        except Exception:
            return None

    est = fit(A, codes, n)
    if est is None:
        return None
    rng = np.random.default_rng(seed)
    winters = np.arange(n)
    idx_by_w = [np.flatnonzero(codes == w) for w in winters]
    bs = []
    for _ in range(n_boot):
        pick = rng.integers(0, n, n)
        rows = np.concatenate([idx_by_w[k] for k in pick])
        cb = np.concatenate([np.full(len(idx_by_w[k]), j) for j, k in enumerate(pick)])
        v = fit(A[rows], cb, cb.max() + 1)
        if v is not None and np.isfinite(v):
            bs.append(v)
    bs = np.array(bs)
    return {
        "effect": round(float(est), 4),
        "CI95": [round(float(np.percentile(bs, 2.5)), 4),
                 round(float(np.percentile(bs, 97.5)), 4)] if len(bs) else None,
        "p_one_sided_negative": float((bs >= 0).mean()) if len(bs) else None,
        "p_two_sided": float(2 * min((bs >= 0).mean(), (bs <= 0).mean())) if len(bs) else None,
        "n_winters": int(n), "n_days": int(len(d)), "n_boot": int(len(bs)),
    }


def between_winter(s, ev, a, b):
    """Field-standard composite: SSW window vs DOY-matched days in non-SSW winters."""
    d = s[np.isin(s.index.month, SEASON)]
    w = pd.Series(winter_of(d.index), index=d.index)
    ssw_w = set(winter_of(pd.DatetimeIndex(ev)))
    ctrl_mask = (~w.isin(ssw_w)).values
    diffs = []
    for o in ev:
        m = (d.index >= o + pd.Timedelta(days=a)) & (d.index <= o + pd.Timedelta(days=b))
        if m.sum() < 5:
            continue
        doys = set(d.index[m].dayofyear)
        ctrl = d[np.isin(d.index.dayofyear, list(doys)) & ctrl_mask]
        if len(ctrl) < 5:
            continue
        diffs.append(float(d[m].mean() - ctrl.mean()))
    diffs = np.array(diffs)
    if len(diffs) < 3:
        return None
    return {
        "effect": round(float(np.mean(diffs)), 4),
        "n_events": len(diffs),
        "events_negative": f"{int((diffs < 0).sum())}/{len(diffs)}",
        "t_p_two_sided": float(stats.ttest_1samp(diffs, 0).pvalue),
        "sign_p_one_sided_neg": float(stats.binomtest(
            int((diffs < 0).sum()), len(diffs), 0.5, alternative="greater").pvalue),
    }


def main():
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    ev_all = pd.DatetimeIndex(pd.to_datetime(can["date"]).dropna())
    series = load_series()

    res = {"n_ssw_catalog": len(ev_all), "series": {}}
    print(f"Canonical SSWs: {len(ev_all)}")
    for name, s in series.items():
        ev = ev_all[(ev_all >= s.index.min()) & (ev_all <= s.index.max())]
        print(f"\n{'='*72}\n### {name}: {s.index.min().date()}..{s.index.max().date()}, "
              f"{len(ev)} SSWs covered")
        res["series"][name] = {"n_ssw": len(ev), "windows": {}}
        for lab, (a, b) in WINDOWS.items():
            d = frame(s, ev, a, b)
            if d.empty or d["W"].nunique() < 2:
                continue
            ww = within_winter(d)
            bw = between_winter(s, ev, a, b)
            res["series"][name]["windows"][lab] = {"within_winter": ww,
                                                   "between_winter": bw}
            if ww and bw:
                print(f"  {lab:12s} within {ww['effect']:+7.3f} "
                      f"[{ww['CI95'][0]:+.3f},{ww['CI95'][1]:+.3f}] P={ww['p_two_sided']:.4f}"
                      f"   |  between {bw['effect']:+7.3f} "
                      f"({bw['events_negative']} neg, P={bw['t_p_two_sided']:.4f})")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)

    ao = res["series"].get("AO", {}).get("windows", {}).get("post_0_60", {})
    if ao.get("within_winter"):
        e = ao["within_winter"]
        print(f"\n{'='*72}\nVERDICT (AO, days 0-60, Baldwin-Dunkerton window)")
        print(f"  within-winter effect = {e['effect']:+.3f} "
              f"[{e['CI95'][0]:+.3f}, {e['CI95'][1]:+.3f}], P = {e['p_two_sided']:.4g}")
        if e["effect"] < -0.15 and e["p_two_sided"] < 0.05:
            print("  -> DESIGN VALIDATED: it recovers the canonical negative-AO "
                  "response.\n     The avalanche nulls in R106-R111 are therefore "
                  "informative.")
        else:
            print("  -> DESIGN FAILS ITS POSITIVE CONTROL: it cannot recover a "
                  "signal that is\n     beyond dispute, so the avalanche nulls are "
                  "NOT evidence of absence.")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
