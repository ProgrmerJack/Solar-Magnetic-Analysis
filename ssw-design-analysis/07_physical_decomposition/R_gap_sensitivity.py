#!/usr/bin/env python3
"""
R_gap_sensitivity.py
====================
THE TEST THAT DECIDES RESULT A.

R = within-winter / between-winter is flat across lag in observations (0.79
pre-onset, 0.82 post), which is why I concluded R measures temporal
CONCENTRATION rather than a causal fraction. That conclusion rests entirely on
the within-winter baseline, which is defined by an arbitrary analyst choice:
days more than BASELINE_GAP from any onset, with GAP fixed at 75.

The superseded sensitivity table shows the pre-onset estimate is far more
gap-sensitive than the post-onset one (-45..-31 moves -0.346 -> -0.455 -> -0.540
across gaps 60/75/90, while +15..+29 moves -1.000 -> -1.089 -> -1.138). If the
pre-onset numerator swings by a third under a choice I never tested, the
flat-lag finding is not established.

So: R(lag) at gaps 60, 75, 90.

  R stays flat at every gap        -> flat-lag is robust; R is concentration,
                                      Result A closes as a methods caution
  R becomes lag-dependent at some  -> the flat profile was an artefact of GAP=75
  gap                                 and Result A revives

Either outcome is decisive and neither depends on finding a new result.

Output: R_gap_sensitivity.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "5_mechanism"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
import causal_timescale_ratio as CTR                # noqa: E402

BINS = [(-45, -31), (-30, -16), (-15, -1), (0, 14), (15, 29), (30, 44)]
LABELS = [f"{a:+d}..{b:+d}" for a, b in BINS]
GAPS = [60, 75, 90]
N_BOOT = 1200


def pairs(d, onsets, win, gap):
    ssw_w = {int(o.year + 1 if o.month >= 11 else o.year) for o in onsets}
    non = d[~d["winter"].isin(ssw_w)]
    clim = non.groupby("doy")["y"].mean()
    out = []
    for o in onsets:
        lo, hi = o + pd.Timedelta(days=win[0]), o + pd.Timedelta(days=win[1])
        w = d[(d.index >= lo) & (d.index <= hi)]
        if len(w) < 6:
            continue
        aw = w["y"].mean() - clim.reindex(w["doy"]).mean()
        wtr = int(o.year + 1 if o.month >= 11 else o.year)
        own = d[d["winter"] == wtr]
        lag = (own.index - o).days.values
        base = own[np.abs(lag) > gap]
        if len(base) < 12:
            continue
        ab = base["y"].mean() - clim.reindex(base["doy"]).mean()
        if np.isfinite(aw) and np.isfinite(ab):
            out.append((float(aw), float(aw - ab), wtr))
    return out


def main():
    d, on = CTR.obs()
    res = {"gaps": GAPS, "bins": LABELS, "n_boot": N_BOOT, "table": {}}
    rng = np.random.default_rng(20260801)
    print("R(lag) across baseline gaps -- does the flat profile survive?\n")
    print(f"{'bin':10s}" + "".join(f"{'gap ' + str(g):>22s}" for g in GAPS))
    print("-" * 76)
    store = {}
    for lab, win in zip(LABELS, BINS):
        row = []
        for g in GAPS:
            p = pairs(d, on, win, g)
            if len(p) < 15:
                row.append(None); continue
            B = np.array([x[0] for x in p]); W = np.array([x[1] for x in p])
            wid = np.array([x[2] for x in p]); uw = np.unique(wid)
            R = W.mean() / B.mean() if abs(B.mean()) > 1e-9 else np.nan
            rs = []
            for _ in range(N_BOOT):
                pk = rng.choice(uw, len(uw), replace=True)
                sel = np.concatenate([np.flatnonzero(wid == q) for q in pk])
                bb = B[sel].mean()
                if abs(bb) > 1e-9:
                    rs.append(W[sel].mean() / bb)
            ci = [float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))] if rs else None
            row.append({"between": round(float(B.mean()), 4),
                        "within": round(float(W.mean()), 4),
                        "R": round(float(R), 4),
                        "R_CI95": [round(ci[0], 3), round(ci[1], 3)] if ci else None,
                        "n": int(len(B))})
        store[lab] = row
        cells = "".join(
            f"{(str(round(c['R'],2)) + ' [' + str(round(c['R_CI95'][0],2)) + ',' + str(round(c['R_CI95'][1],2)) + ']') if c else 'n/a':>22s}"
            for c in row)
        print(f"{lab:10s}{cells}")
    res["table"] = store

    print("\n=== VERDICT ===")
    pre = [store["-30..-16"][i]["R"] for i in range(len(GAPS))
           if store["-30..-16"][i]]
    post = [store["+15..+29"][i]["R"] for i in range(len(GAPS))
            if store["+15..+29"][i]]
    print(f"  pre-onset  -30..-16 R across gaps: "
          f"{', '.join(f'{v:+.2f}' for v in pre)}   spread {max(pre)-min(pre):.2f}")
    print(f"  post-onset +15..+29 R across gaps: "
          f"{', '.join(f'{v:+.2f}' for v in post)}   spread {max(post)-min(post):.2f}")
    flat = [abs(p - q) for p, q in zip(pre, post)]
    print(f"  |R_pre - R_post| at each gap: {', '.join(f'{v:.2f}' for v in flat)}")
    verdict = ("FLAT AT EVERY GAP -- R is concentration, Result A closes"
               if max(flat) < 0.25 else
               "LAG-DEPENDENT AT SOME GAP -- the flat profile was a GAP=75 artefact")
    print(f"  -> {verdict}")
    res["verdict"] = verdict
    res["pre_R_by_gap"] = pre
    res["post_R_by_gap"] = post
    (RESULTS / "R_gap_sensitivity.json").write_text(json.dumps(res, indent=2),
                                                  encoding="utf8")
    print("\nSaved -> R_gap_sensitivity.json")


if __name__ == "__main__":
    main()
