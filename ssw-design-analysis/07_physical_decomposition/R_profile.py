#!/usr/bin/env python3
"""
R_profile.py
============
R as a function of lag -- and the test that REFUTED the causal reading of R.

THE HYPOTHESIS THIS SCRIPT WAS BUILT TO TEST (and which it killed)
  R = within-winter / between-winter was read as a causal fraction:
    pre-onset R ~ 0   the "precursor" is WINTER-scale, i.e. common cause
    pre-onset R ~ 1   a genuine coupled precursor mechanism
  A precursor cannot be caused by the event that follows it, so a low pre-onset
  R was expected to separate common cause from causal downward influence.

THE RESULT: R IS FLAT ACROSS LAG, SO IT IS NOT A CAUSAL FRACTION
                      observations (n=43)      CanESM5 (10 members)
    -30..-16 (pre)    +0.789 [0.34, 1.16]      +0.754 [0.65, 0.84]  n=794
    +15..+29 (post)   +0.824 [0.49, 1.09]      +0.825 [0.76, 0.89]  n=846
  Flat in BOTH systems, and in CanESM5 the intervals are only +-0.10 wide, so
  this is not an n=43 artefact. Flat at every baseline gap tested (60/75/90 d,
  see R_gap_sensitivity.py).

  MECHANISM: `within` subtracts a baseline that excludes |lag| <= BASELINE_GAP,
  so ANY anomaly inside that window -- before or after onset -- is absent from
  the baseline and returns R ~ 1. R measures TEMPORAL CONCENTRATION of the
  anomaly around the onset date, not causal fraction.

WHY THE VALIDATION DID NOT CATCH THIS
  validate_R_diagnostic.py simulates only two structures: whole-winter
  depression (f=0) and post-onset-only depression (f=1, `y[p:p+61] -= c`). It
  never simulates a precursor, so its "|bias| <= 0.047 across f in [0,1]" is
  true only within that two-structure family and does NOT license the causal
  reading. Add a precursor arm before any R-based claim is made again.

WHAT SURVIVES
  The association is temporally LOCALISED around onset rather than winter-scale.
  A short-timescale shared driver produces the same signature, so this is not
  evidence of downward causation.

Output: results/R_profile.json
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
import causal_timescale_ratio as CTR                # noqa: E402
import ensemble_precursor as EP                     # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
BINS = [(-60, -46), (-45, -31), (-30, -16), (-15, -1),
        (0, 14), (15, 29), (30, 44), (45, 60)]
LABELS = [f"{a:+d}..{b:+d}" for a, b in BINS]
BASELINE_GAP = 75
N_BOOT = 1500
SEASON = (11, 12, 1, 2, 3, 4)
MIN_DENOM = 0.15          # R is meaningless when the composite is ~0


def pairs_at(d, onsets, win, ycol="y", wid_of=None):
    """Per-event (between, within, winter) at an arbitrary lag window."""
    wid_of = wid_of or (lambda o: int(o.year + 1 if o.month >= 11 else o.year))
    ssw_w = {wid_of(o) for o in onsets}
    non = d[~d["winter"].isin(ssw_w)]
    if len(non) < 200:
        return []
    clim = non.groupby("doy")[ycol].mean()
    out = []
    for o in onsets:
        lo = o + pd.Timedelta(days=win[0])
        hi = o + pd.Timedelta(days=win[1])
        w = d[(d.index >= lo) & (d.index <= hi)]
        if len(w) < 6:
            continue
        adj_w = w[ycol].mean() - clim.reindex(w["doy"]).mean()
        own = d[d["winter"] == wid_of(o)]
        if not len(own):
            continue
        lag = (own.index - o).days.values
        base = own[np.abs(lag) > BASELINE_GAP]
        if len(base) < 15:
            continue
        adj_b = base[ycol].mean() - clim.reindex(base["doy"]).mean()
        if np.isfinite(adj_w) and np.isfinite(adj_b):
            out.append((float(adj_w), float(adj_w - adj_b), wid_of(o)))
    return out


def summarise(p, seed=20260801):
    """Winter-block bootstrap of R over a pooled list of per-event triples."""
    rng = np.random.default_rng(seed)
    B = np.array([x[0] for x in p]); W = np.array([x[1] for x in p])
    wid = np.array([x[2] for x in p]); uw = np.unique(wid)
    R = W.mean() / B.mean() if abs(B.mean()) > 1e-9 else np.nan
    rs = []
    for _ in range(N_BOOT):
        pick = rng.choice(uw, len(uw), replace=True)
        sel = np.concatenate([np.flatnonzero(wid == x) for x in pick])
        bb = B[sel].mean()
        if abs(bb) > 1e-9:
            rs.append(W[sel].mean() / bb)
    rs = np.array(rs)
    return {"between": round(float(B.mean()), 4),
            "within": round(float(W.mean()), 4),
            "R": round(float(R), 4),
            "R_CI95": [round(float(np.percentile(rs, 2.5)), 3),
                       round(float(np.percentile(rs, 97.5)), 3)]
            if len(rs) else None,
            "n": int(len(B)),
            "interpretable": bool(abs(B.mean()) > MIN_DENOM)}


def profile(d, onsets, ycol="y", wid_of=None, seed=20260801):
    rows = {}
    for (a, b), lab in zip(BINS, LABELS):
        p = pairs_at(d, onsets, (a, b), ycol, wid_of)
        rows[lab] = summarise(p, seed) if len(p) >= 15 else None
    return rows


def profile_members(members, ycol="y", seed=20260801):
    """Ensemble profile pooling PAIRS across members, never ROWS.

    Members of one model share the calendar 1850-2014, so concatenating their
    rows and then selecting an event window by date pulls in every member --
    only one of which had the SSW. That diluted the CanESM5 composite by ~3.5x
    (between -0.203 against the correct -0.716 at +15..+29) and is the route
    multimodel_R.py already avoids by building pairs one member at a time.
    """
    rows = {}
    for (a, b), lab in zip(BINS, LABELS):
        p = []
        for d, onsets, wid_of in members:
            p += pairs_at(d, onsets, (a, b), ycol, wid_of)
        rows[lab] = summarise(p, seed) if len(p) >= 15 else None
    return rows


def show(name, rows):
    print(f"\n=== {name} ===")
    print(f"{'bin':10s} {'between':>9s} {'within':>9s} {'R':>8s} {'95% CI':>18s}")
    print("-" * 60)
    for lab in LABELS:
        r = rows.get(lab)
        if not r:
            print(f"{lab:10s}  (too few events)")
            continue
        ci = r["R_CI95"] or [np.nan, np.nan]
        mark = "" if r["interpretable"] else "   <- composite ~0, R undefined"
        print(f"{lab:10s} {r['between']:+9.3f} {r['within']:+9.3f} {r['R']:+8.2f} "
              f"[{ci[0]:+7.2f},{ci[1]:+7.2f}]{mark}")


def main():
    out = {"bins": LABELS, "min_denominator": MIN_DENOM, "datasets": {}}
    d, on = CTR.obs()
    obs = profile(d, on)
    out["datasets"]["observations"] = obs
    show("OBSERVATIONS  (43 events)", obs)

    # the largest ensemble, as the high-precision check on the same profile
    files = sorted(RAW.glob("CanESM5_*_zm.nc"))
    if files:
        members, n_ev = [], 0
        for k, f in enumerate(files):
            m = EP.load_member(f)
            full = m
            m = m[np.isin(m.index.month, SEASON)].dropna()
            o = EP.detect_ssw(full["u10"].values, full.index)   # full daily series: CP07 needs contiguous days
            dm = m.rename(columns={"am": "y"})[["y"]].copy()
            dm["doy"] = dm.index.dayofyear
            off = 100000 * (k + 1)
            dm["winter"] = EP.winter_of(dm.index) + off
            members.append((dm, o, lambda x, off=off:
                            int((x.year + 1 if x.month >= 11 else x.year) + off)))
            n_ev += len(o)
        rows = profile_members(members)
        out["datasets"]["CanESM5"] = rows
        show(f"CanESM5  ({n_ev} events, {len(files)} members)", rows)

    print("\n=== THE DECOMPOSITION -- THE HYPOTHESIS ABOVE IS REFUTED ===")
    for name in ("observations", "CanESM5"):
        rows = out["datasets"].get(name)
        if not rows:
            continue
        pre, post = rows.get("-30..-16"), rows.get("+15..+29")
        if not (pre and post):
            continue
        print(f"  {name:14s} pre-onset R = {pre['R']:+.3f} {pre['R_CI95']}   "
              f"post-onset R = {post['R']:+.3f} {post['R_CI95']}")
    print("""
  R is FLAT across lag, in observations (n=43) AND in CanESM5 (n=794/846,
  CI +-0.10). A precursor cannot be caused by the event that follows it, so
  R at a pre-onset lag cannot be a causal fraction -- and it is the same 0.78
  as the post-onset value in both systems.

  MECHANISM: `within` subtracts a baseline that excludes |lag| <= 75 d, so ANY
  anomaly living inside that window -- before or after onset -- is absent from
  the baseline and returns R ~ 1. R measures TEMPORAL CONCENTRATION of the
  anomaly around the onset date, not causal fraction.

  WHY THE VALIDATION MISSED IT: validate_R_diagnostic.py simulates only two
  structures -- whole-winter depression (f=0) and post-onset-only depression
  (f=1, `y[p:p+61] -= c`). It never simulates a precursor, so it could not
  discover that a pre-onset anomaly also scores R ~ 0.8.

  The surviving claim is narrower: the association is temporally LOCALISED
  around onset rather than winter-scale. A short-timescale shared driver
  produces the same signature.""")

    (RESULTS / "R_profile.json").write_text(json.dumps(out, indent=2), encoding="utf8")
    print("\nSaved -> R_profile.json")


if __name__ == "__main__":
    main()
