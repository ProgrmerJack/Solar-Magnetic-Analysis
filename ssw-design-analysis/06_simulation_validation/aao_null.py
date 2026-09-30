#!/usr/bin/env python3
"""
aao_null.py
===========
The AAO negative control FAILS, and it is not seed noise: -0.659 at +45..+60 d,
p=0.011..0.023 across 60 fixed seeds at 6,000 replicates (negcontrol_stability).
A Northern Hemisphere warming should not move the Southern annular mode.

Two readings, with very different consequences:

  A. THE ESTIMATOR MANUFACTURES SIGNAL.  Gate 6 fails, and every AO result in the
     project is suspect along with it.
  B. THE AAO SERIES BREAKS AN ASSUMPTION THE AO DOES NOT.  Gate 3 calibrated the
     pseudo-onset null on the AO only. The AAO is Southern-Hemisphere SUMMER
     during the NH SSW season, and it carries a strong ozone-driven positive
     trend; 3 annual harmonics plus winter fixed effects may not absorb either.
     Then the failure indicts the estimator FOR THE AAO, and says nothing about
     the AO -- whose own null is clean (Gate 3, bias +0.033, coverage 0.947).

This distinguishes them the only way that works: run the SAME clean pseudo-onset
null on the AAO series. Real-event influence is deleted from the data first
(-60..+75 d), then fake onsets are drawn from the surviving days at the observed
day-of-year set, exactly as gate3_clean_null.py does for the AO.

  If fake AAO onsets also produce ~-0.66 at +45..+60, the bin is picking up
  residual structure in the AAO, not a response -- reading B, and the AO results
  stand.
  If fake onsets centre on zero, the estimator is clean on the AAO too and the
  real-onset result is a genuine anomaly that needs explaining -- reading A.

Output: aao_null.json
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "3_calibration"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE))
import calibrate_seasonality as C
from build_catalogue import load_catalogue
import gate3_clean_null as G
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
import multi_index_event_study as M

N_DRAW = 2000
WATCH = "+45..+60"
ROOT = HERE.parents[1]


def main():
    d = M.load("aao")
    real = load_catalogue("primary")
    real = real[(real >= d.index.min()) & (real <= d.index.max())]
    print(f"AAO: {len(d)} winter days {d.index.min().date()}..{d.index.max().date()}, "
          f"{len(real)} real onsets in range")

    # observed profile, for reference
    prof, _ = M.run(d, real, seed=20260730)
    obs = prof[WATCH]["effect"]
    print(f"observed {WATCH}: {obs:+.3f} (p={prof[WATCH]['p_two_sided']:.4f})")

    # delete every day under the influence of a real onset, then draw fakes
    mask = G.real_influence_mask(d.index, real)
    clean = d[~mask]
    print(f"clean days after removing real-event zone {G.INFLUENCE}: "
          f"{len(clean)}/{len(d)} ({100*len(clean)/len(d):.0f}%)")

    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    rng = np.random.default_rng(20260730)
    vals, allbins = [], {l: [] for l in M.LABELS}
    for i in range(N_DRAW):
        fake = G.draw_clean(clean.index, doys, rng)
        if len(fake) < 10:
            continue
        dd = M.build(clean, fake)
        b = M.fit(dd)
        if b is None or not np.all(np.isfinite(b)):
            continue
        for j, l in enumerate(M.LABELS):
            allbins[l].append(b[j])
        vals.append(b[M.LABELS.index(WATCH)])
        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{N_DRAW} draws, running mean {np.mean(vals):+.4f}", flush=True)

    v = np.array(vals)
    res = {"n_draw": int(len(v)), "watch_bin": WATCH,
           "observed_effect": round(float(obs), 4),
           "null_mean": round(float(v.mean()), 4),
           "null_sd": round(float(v.std()), 4),
           "null_pct2.5": round(float(np.percentile(v, 2.5)), 4),
           "null_pct97.5": round(float(np.percentile(v, 97.5)), 4),
           "p_observed_vs_null": float((v <= obs).mean()),
           "per_bin_null_mean": {l: round(float(np.mean(allbins[l])), 4) for l in M.LABELS}}

    print(f"\n=== clean pseudo-onset null on the AAO, {len(v)} draws ===")
    print(f"  null mean at {WATCH}: {v.mean():+.4f}  (sd {v.std():.4f})")
    print(f"  null 95% range      : [{np.percentile(v,2.5):+.3f}, {np.percentile(v,97.5):+.3f}]")
    print(f"  observed            : {obs:+.3f}")
    print(f"  P(null <= observed) : {(v <= obs).mean():.4f}")
    print(f"\n  null mean by bin: "
          + "  ".join(f"{l} {np.mean(allbins[l]):+.3f}" for l in M.LABELS))
    print("\n  null far from zero at this bin -> reading B (AAO-specific residual "
          "structure; AO results stand)")
    print("  null centred on zero          -> reading A (estimator clean; the AAO "
          "anomaly is real and unexplained)")
    (RESULTS / "aao_null.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> aao_null.json")


if __name__ == "__main__":
    main()
