#!/usr/bin/env python3
"""
multimodel_R.py
===============
R across the CMIP6 archive, and the test that motivates it:

    DO MODELS THAT REPRODUCE THE OBSERVED SURFACE RESPONSE GET THE MECHANISM RIGHT?

  The field evaluates stratosphere-troposphere coupling with composites: does the
  model produce the observed post-SSW surface anomaly? That comparison cannot see
  mechanism. R can:

      R = within-winter effect / between-winter effect
      R ~ 1  event-scale -- the anomaly follows onset (causal downward influence)
      R ~ 0  winter-scale -- SSW-hosting winters are anomalous throughout
                            (common cause; the composite is right for the wrong reason)

  R is a validated estimator of the causal fraction: bias <= 0.047 across
  f in [0,1] under known truth, with f=0 and f=1 separable at n=43
  (validate_R_diagnostic.py). Observations give R = 0.824, and the pure
  common-cause null is rejected at P < 0.0007.

WHAT IS REPORTED PER MODEL
  between   the conventional composite -- what the field currently compares
  R         the mechanism -- what it cannot see
  and the correlation between the two ACROSS models. If models match the observed
  composite but scatter in R, composite-based evaluation is demonstrably blind to
  mechanism, which is the point.

Guards, because 40 models do not share conventions:
  - the 10 hPa level and 60N latitude are verified, not assumed (load_member raises)
  - SSW rate per winter is reported per model; an implausible rate flags a broken
    detection rather than a model result
  - any model that raises is reported as skipped WITH its reason, never dropped
    silently

Output: multimodel_R.json
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
import canonical_event_study as K                   # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402
import ensemble_precursor as EP                     # noqa: E402
import causal_timescale_ratio as CTR                # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
SEASON = (11, 12, 1, 2, 3, 4)
PLAUSIBLE_RATE = (0.15, 1.20)      # SSWs per winter; observations give ~0.6


def model_pairs(name):
    files = sorted(RAW.glob(f"{name}_*_zm.nc"))
    pairs, n_ssw, n_win = [], 0, 0
    for k, f in enumerate(files):
        m = EP.load_member(f)
        m = m[np.isin(m.index.month, SEASON)].dropna()
        if len(m) < 5000:
            continue
        o = EP.detect_ssw(m["u10"].values, m.index)
        d = m.rename(columns={"am": "y"})[["y"]].copy()
        d["doy"] = d.index.dayofyear
        off = 100000 * (k + 1)
        d["winter"] = EP.winter_of(d.index) + off
        n_ssw += len(o)
        n_win += d["winter"].nunique()
        pairs += CTR.per_event_pairs(
            d, o, wid_of=lambda x, off=off:
                int((x.year + 1 if x.month >= 11 else x.year) + off))
    return pairs, n_ssw, n_win, len(files)


def main():
    res = {"window": list(CTR.WINDOW), "models": {}, "skipped": {}}

    d, on = CTR.obs()
    o = CTR.R_from_pairs(CTR.per_event_pairs(d, on))
    res["observations"] = o
    print(f"{'dataset':18s} {'n_mem':>5s} {'SSW/win':>8s} {'between':>9s} "
          f"{'within':>8s} {'R':>7s} {'95% CI':>16s}")
    print("-" * 78)
    print(f"{'OBSERVATIONS':18s} {'-':>5s} {43/36:8.2f} {o['between']:+9.3f} "
          f"{o['within']:+8.3f} {o['R']:+7.2f} "
          f"[{o['R_CI95'][0]:+6.2f},{o['R_CI95'][1]:+6.2f}]")

    names = sorted({f.stem.split("_")[0] for f in RAW.glob("*_zm.nc")})
    for name in names:
        try:
            pairs, n_ssw, n_win, n_mem = model_pairs(name)
        except Exception as e:
            res["skipped"][name] = f"{type(e).__name__}: {str(e)[:110]}"
            print(f"{name:18s} SKIPPED -- {str(e)[:60]}")
            continue
        if len(pairs) < 20:
            res["skipped"][name] = f"only {len(pairs)} usable events"
            print(f"{name:18s} SKIPPED -- only {len(pairs)} usable events")
            continue
        rate = n_ssw / n_win if n_win else 0
        r = CTR.R_from_pairs(pairs)
        r.update({"n_members": n_mem, "ssw_per_winter": round(rate, 3),
                  "n_winters": int(n_win)})
        flag = "" if PLAUSIBLE_RATE[0] <= rate <= PLAUSIBLE_RATE[1] else "  <-- implausible SSW rate"
        r["rate_plausible"] = not bool(flag)
        res["models"][name] = r
        ci = r["R_CI95"] or [np.nan, np.nan]
        print(f"{name:18s} {n_mem:5d} {rate:8.2f} {r['between']:+9.3f} "
              f"{r['within']:+8.3f} {r['R']:+7.2f} "
              f"[{ci[0]:+6.2f},{ci[1]:+6.2f}]{flag}")

    # ---- the headline test ----
    use = {k: v for k, v in res["models"].items()
           if v.get("rate_plausible") and v["R_CI95"]
           and abs(v["between"]) > 0.15}          # R needs a non-trivial denominator
    if len(use) >= 4:
        b = np.array([v["between"] for v in use.values()])
        R = np.array([v["R"] for v in use.values()])
        from scipy import stats
        rr, pp = stats.pearsonr(b, R)
        obs_b, obs_R = o["between"], o["R"]
        print(f"\n=== DOES MATCHING THE COMPOSITE IMPLY MATCHING THE MECHANISM? ===")
        print(f"  {len(use)} models with a usable composite and R")
        print(f"  composite (between-winter): observed {obs_b:+.3f}, "
              f"models {b.min():+.3f}..{b.max():+.3f}")
        print(f"  mechanism (R):              observed {obs_R:+.3f}, "
              f"models {R.min():+.3f}..{R.max():+.3f}")
        print(f"  corr(composite, R) across models = {rr:+.3f}  (P = {pp:.3f})")
        near = {k: v for k, v in use.items() if abs(v["between"] - obs_b) < 0.25}
        if len(near) >= 3:
            rs = [v["R"] for v in near.values()]
            print(f"\n  Of the {len(near)} models within 0.25 sigma of the observed "
                  f"composite,\n  R ranges {min(rs):+.2f} to {max(rs):+.2f} "
                  f"(observed {obs_R:+.2f}).")
            print("  -> matching the composite does NOT pin down the mechanism"
                  if (max(rs) - min(rs)) > 0.3 else
                  "  -> models matching the composite also agree on mechanism")
        if 0 < len(near) < 3:
            print(f"\n  only {len(near)} model(s) within 0.25 sigma of the observed"
                  f" composite -- too few to say anything about whether matching"
                  f" the composite implies matching the mechanism")
        res["headline"] = {"n_models": len(use), "corr_composite_R": round(float(rr), 4),
                           "p": float(pp),
                           "R_range": [round(float(R.min()), 3), round(float(R.max()), 3)],
                           "composite_range": [round(float(b.min()), 3),
                                               round(float(b.max()), 3)],
                           "n_near_observed_composite": len(near),
                           "R_range_among_near": (
                               [round(float(min(v["R"] for v in near.values())), 3),
                                round(float(max(v["R"] for v in near.values())), 3)]
                               if len(near) >= 3 else None)}
    else:
        print(f"\n  only {len(use)} usable models -- headline test needs >=4")

    (RESULTS / "multimodel_R.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print(f"\nSaved -> multimodel_R.json")


if __name__ == "__main__":
    main()
