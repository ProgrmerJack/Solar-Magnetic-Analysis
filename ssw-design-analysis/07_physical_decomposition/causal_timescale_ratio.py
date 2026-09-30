#!/usr/bin/env python3
"""
causal_timescale_ratio.py
=========================
Turns this project's central NULL result into a causal diagnostic.

THE REFRAMING
  `FINDING_design_robustness.md` reports that the conventional between-winter
  composite (-0.954) and the within-winter estimator (-0.952) agree to 0.002
  sigma, and treats that as a dead end: no estimator gap left to explain.

  That reading was wrong. The two estimators contrast different things:

    between-winter   SSW-window days vs day-of-year-matched days from winters
                     with NO event. Includes anything that makes an SSW-hosting
                     winter anomalous as a whole.
    within-winter    the same window vs THAT winter's own baseline, >75 d from
                     onset. Isolates what changes AFTER onset.

  So their ratio is informative about mechanism:

      R = within-winter effect / between-winter effect

      R ~ 1   the response is EVENT-scale. The anomaly appears after onset
              relative to the same winter. Consistent with the stratosphere
              causing the surface anomaly.
      R ~ 0   the response is WINTER-scale. SSW-hosting winters are anomalous
              throughout, and nothing extra happens at onset. This is the
              signature of a COMMON CAUSE -- something makes a winter both
              SSW-prone and surface-anomalous.

  R is exactly the quantity a common-cause confounder attacks, and it needs no
  experiment: it is computable from any observational record or model run.

WHY THIS IS WORTH DOING
  MIROC6 forced the point. Its raw post-onset composite is -0.535 but its
  within-winter estimate is +0.055 -- R ~ 0. CanESM5 and the observations both
  have R ~ 1. The models are not disagreeing about the size of the surface
  response; they disagree about WHETHER THE STRATOSPHERE CAUSES IT. A model can
  reproduce the observed composite for entirely the wrong reason, and the
  standard composite cannot tell.

INFERENCE
  Winter-block bootstrap on both numerator and denominator jointly, so R carries
  a proper interval. Winters are the resampling unit throughout, as everywhere
  else in this project.

Output: causal_timescale_ratio.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "5_mechanism"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
import calibrate_seasonality as C                   # noqa: E402
import canonical_event_study as K                   # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402
import ensemble_precursor as EP                     # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
WINDOW = (15, 29)          # the bin where the observed AO response peaks
BASELINE_GAP = 75
N_BOOT = 1500
SEASON = (11, 12, 1, 2, 3, 4)


def per_event_pairs(d, onsets, ycol="y", wid_of=None):
    """Per-event (between-winter, within-winter, winter-id) triples.

    Working per event rather than returning two pooled means keeps pooling across
    ensemble members and the winter-block bootstrap straightforward, and avoids
    the winter-numbering collision that made the first version return NaN for
    every model: onsets were offset into per-member winter blocks but the winter
    of an onset was then recomputed from its calendar date, which has no offset.
    """
    wid_of = wid_of or (lambda o: int(o.year + 1 if o.month >= 11 else o.year))
    ssw_w = {wid_of(o) for o in onsets}
    non = d[~d["winter"].isin(ssw_w)]
    if len(non) < 200:
        return []
    clim = non.groupby("doy")[ycol].mean()
    out = []
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=WINDOW[0]), o + pd.Timedelta(days=WINDOW[1])
        w = d[(d.index >= lo) & (d.index <= hi)]
        if len(w) < 8:
            continue
        adj_w = w[ycol].mean() - clim.reindex(w["doy"]).mean()
        wtr = wid_of(o)
        own = d[d["winter"] == wtr]
        if not len(own):
            continue
        lag = (own.index - o).days.values
        base = own[np.abs(lag) > BASELINE_GAP]
        if len(base) < 15:
            continue
        adj_b = base[ycol].mean() - clim.reindex(base["doy"]).mean()
        if np.isfinite(adj_w) and np.isfinite(adj_b):
            out.append((float(adj_w), float(adj_w - adj_b), wtr))
    return out


def R_from_pairs(pairs, seed=20260731):
    if not pairs:
        return {"between": None, "within": None, "R": None,
                "R_CI95": None, "n_events": 0}
    b = np.array([p[0] for p in pairs])
    w = np.array([p[1] for p in pairs])
    wid = np.array([p[2] for p in pairs])
    R = float(w.mean() / b.mean())
    rng = np.random.default_rng(seed)
    uw = np.unique(wid)
    rs = []
    for _ in range(N_BOOT):
        pick = rng.choice(uw, len(uw), replace=True)
        sel = np.concatenate([np.flatnonzero(wid == x) for x in pick])
        bb, ww = b[sel].mean(), w[sel].mean()
        if abs(bb) > 1e-6:
            rs.append(ww / bb)
    rs = np.array(rs)
    return {"between": round(float(b.mean()), 4),
            "within": round(float(w.mean()), 4), "R": round(R, 4),
            "R_CI95": [round(float(np.percentile(rs, 2.5)), 3),
                       round(float(np.percentile(rs, 97.5)), 3)] if len(rs) else None,
            "n_events": int(len(b))}


def two_estimates(d, onsets, ycol="y"):
    """(between-winter, within-winter) effect over WINDOW, DOY-adjusted."""
    wi = d["winter"].values
    ssw_w = set(np.asarray([w for w in
                            d["winter"].values[np.isin(d.index, onsets)]]))
    # a winter counts as SSW-hosting if any onset falls in it
    ssw_w = set()
    for o in onsets:
        m = d.index == o
        if m.any():
            ssw_w.add(int(d["winter"].values[m][0]))
        else:
            ssw_w.add(int(o.year + 1 if o.month >= 11 else o.year))
    non = d[~d["winter"].isin(ssw_w)]
    if len(non) < 200:
        return np.nan, np.nan
    clim = non.groupby("doy")[ycol].mean()

    win_vals, btw_vals, wtn_vals = [], [], []
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=WINDOW[0]), o + pd.Timedelta(days=WINDOW[1])
        w = d[(d.index >= lo) & (d.index <= hi)]
        if len(w) < 8:
            continue
        adj_w = w[ycol].mean() - clim.reindex(w["doy"]).mean()
        # between-winter: window vs DOY-matched non-event winters (adj_w already is)
        btw_vals.append(adj_w)
        # within-winter: same window vs that winter's own distant baseline
        wtr = int(o.year + 1 if o.month >= 11 else o.year)
        own = d[d["winter"] == wtr]
        lag = (own.index - o).days.values
        base = own[np.abs(lag) > BASELINE_GAP]
        if len(base) < 15:
            continue
        adj_b = base[ycol].mean() - clim.reindex(base["doy"]).mean()
        wtn_vals.append(adj_w - adj_b)
    if not btw_vals or not wtn_vals:
        return np.nan, np.nan
    return float(np.mean(btw_vals)), float(np.mean(wtn_vals))


def bootstrap_R(d, onsets, ycol="y", seed=20260731):
    btw, wtn = two_estimates(d, onsets, ycol)
    rng = np.random.default_rng(seed)
    winters = np.array(sorted(d["winter"].unique()))
    idx = {w: np.flatnonzero((d["winter"] == w).values) for w in winters}
    on_by_w = {}
    for o in onsets:
        wtr = int(o.year + 1 if o.month >= 11 else o.year)
        on_by_w.setdefault(wtr, []).append(o)
    rs, bs, ws = [], [], []
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        rows, ons, wmap = [], [], []
        for j, w in enumerate(pick):
            rows.append(idx[w])
            for o in on_by_w.get(int(w), []):
                ons.append(o)
        sub = d.iloc[np.concatenate(rows)].copy()
        if not ons:
            continue
        try:
            b, wv = two_estimates(sub, pd.DatetimeIndex(sorted(set(ons))), ycol)
            if np.isfinite(b) and np.isfinite(wv) and abs(b) > 1e-6:
                rs.append(wv / b)
                bs.append(b)
                ws.append(wv)
        except Exception:
            pass
    rs = np.array(rs)
    return {"between": round(btw, 4), "within": round(wtn, 4),
            "R": round(wtn / btw, 4) if btw else None,
            "R_CI95": [round(float(np.percentile(rs, 2.5)), 3),
                       round(float(np.percentile(rs, 97.5)), 3)] if len(rs) else None,
            "n_boot_ok": int(len(rs))}


def obs():
    d = K.load_series("AO")
    on = load_catalogue("primary")
    on = on[(on >= d.index.min()) & (on <= d.index.max())]
    return d, on


def model(name):
    files = sorted(RAW.glob(f"{name}_*_zm.nc"))
    frames, ons = [], []
    for f in files:
        m = EP.load_member(f)
        full = m
        m = m[np.isin(m.index.month, SEASON)].dropna()
        o = EP.detect_ssw(full["u10"].values, full.index)   # full daily series: CP07 needs contiguous days
        mem = f.stem.split("_")[1]
        d = m.rename(columns={"am": "y"})[["y"]].copy()
        d["doy"] = d.index.dayofyear
        # offset each member into its own winter-numbering block so winters
        # never collide across members
        base = EP.winter_of(d.index) + 100000 * (len(frames) + 1)
        d["winter"] = base
        frames.append(d)
        ons.append(o)
    return pd.concat(frames), pd.DatetimeIndex(np.concatenate(
        [o.values for o in ons])) if ons else pd.DatetimeIndex([])


def main():
    res = {"window": list(WINDOW), "n_boot": N_BOOT, "datasets": {}}
    print("R = within-winter / between-winter effect at "
          f"+{WINDOW[0]}..+{WINDOW[1]} d\n")
    print("  R ~ 1  event-scale response -> consistent with CAUSAL stratospheric influence")
    print("  R ~ 0  winter-scale response -> signature of a COMMON CAUSE\n")
    print(f"{'dataset':14s} {'between':>9s} {'within':>9s} {'R':>7s} {'95% CI':>18s}")
    print("-" * 62)

    d, on = obs()
    r = R_from_pairs(per_event_pairs(d, on))
    res["datasets"]["observations"] = r
    print(f"{'observations':14s} {r['between']:+9.3f} {r['within']:+9.3f} "
          f"{r['R']:+7.2f} [{r['R_CI95'][0]:+7.2f},{r['R_CI95'][1]:+7.2f}]"
          f"  n={r['n_events']}")

    for name in ("CanESM5", "MIROC6"):
        files = sorted(RAW.glob(f"{name}_*_zm.nc"))
        if not files:
            continue
        pairs = []
        for k, f in enumerate(files):
            m = EP.load_member(f)
            full = m
            m = m[np.isin(m.index.month, SEASON)].dropna()
            o = EP.detect_ssw(full["u10"].values, full.index)   # full daily series: CP07 needs contiguous days
            dm = m.rename(columns={"am": "y"})[["y"]].copy()
            dm["doy"] = dm.index.dayofyear
            dm["winter"] = EP.winter_of(dm.index) + 100000 * (k + 1)
            off = 100000 * (k + 1)
            pairs += per_event_pairs(
                dm, o, wid_of=lambda x, off=off:
                    int((x.year + 1 if x.month >= 11 else x.year) + off))
        rm = R_from_pairs(pairs)
        res["datasets"][name] = rm
        ci = rm["R_CI95"] or [float("nan")] * 2
        print(f"{name:14s} {rm['between']:+9.3f} {rm['within']:+9.3f} "
              f"{rm['R']:+7.2f} [{ci[0]:+7.2f},{ci[1]:+7.2f}]  n={rm['n_events']}")

    (RESULTS / "causal_timescale_ratio.json").write_text(json.dumps(res, indent=2),
                                                       encoding="utf8")
    print("\n  A model can reproduce the observed COMPOSITE while getting R wrong,")
    print("  i.e. the right surface response for the wrong mechanism. The standard")
    print("  composite cannot detect that; R can.")
    print("\nSaved -> causal_timescale_ratio.json")


if __name__ == "__main__":
    main()
