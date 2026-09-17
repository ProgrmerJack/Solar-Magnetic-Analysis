#!/usr/bin/env python3
"""
era_split.py
============
common_period.py showed the pre-onset AO anomaly is -0.82 (p=0.005) on the full
1950-2026 record but -0.37 (p=0.31) on 1980-2019, in every catalogue. That is
consistent with the "precursor" being a pre-satellite-era artefact, but two
separate windows is not a test of a difference.

This fits the difference directly: the same event-study design, with every bin
interacted with a PRE1980 indicator, so the pre/post-1980 contrast is a single
estimated coefficient with its own winter-block bootstrap interval.

  H0: the -30..-16 effect is the same in both eras.

Two candidate explanations if the eras do differ, distinguished at the end:
  DATING   pre-satellite onset dates are less accurate, so a mis-dated event
           smears its post-onset response backwards into the pre-onset bins.
           Predicts: the POST-onset effect should also be attenuated pre-1980
           (the response is spread over more bins, not lost).
  REAL     the atmosphere behaved differently. Predicts: post-onset effect
           comparable across eras, only the precursor differs.

Output: era_split.json
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "2_event_study"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
from build_catalogue import load_catalogue
import canonical_event_study as K
import calibrate_seasonality as C

SPLIT = 1980
N_BOOT = 4000
WATCH = ["-30..-16", "-15..-1", "+0..+14", "+15..+29"]


def fit_interacted(dd):
    """Bin effects, each interacted with pre-1980. Returns (early, late) blocks."""
    y = dd["y"].values.astype(float)
    S = C.harmonics(dd["doy"].values, 3)
    early = (dd["winter"].values < SPLIT).astype(float)
    D = np.column_stack([(dd["bin"].values == i).astype(float) for i in range(len(K.BINS))])
    X = np.column_stack([D * early[:, None], D * (1 - early)[:, None], S])
    A = K.C_demean(np.column_stack([y, X]), dd["winter"].values)
    b, *_ = np.linalg.lstsq(A[:, 1:], A[:, 0], rcond=None)
    n = len(K.BINS)
    return b[:n], b[n:2 * n]


def main():
    d = K.load_series("AO")
    on = load_catalogue("primary")
    on = on[(on >= d.index.min()) & (on <= d.index.max())]
    dd = K.build(d, on)
    n_early = int((pd.DatetimeIndex(on).year < SPLIT).sum())
    print(f"primary events: {n_early} pre-{SPLIT}, {len(on)-n_early} from {SPLIT}")

    e0, l0 = fit_interacted(dd)
    rng = np.random.default_rng(20260730)
    winters = np.array(sorted(dd["winter"].unique()))
    idx = {w: np.flatnonzero((dd["winter"] == w).values) for w in winters}
    boot = []
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        rows = np.concatenate([idx[w] for w in pick])
        sub = dd.iloc[rows].copy()
        sub["winter_orig"] = sub["winter"].values
        sub["winter"] = np.concatenate([np.full(len(idx[w]), j) for j, w in enumerate(pick)])
        try:
            yy = sub["y"].values.astype(float)
            S = C.harmonics(sub["doy"].values, 3)
            early = (sub["winter_orig"].values < SPLIT).astype(float)
            D = np.column_stack([(sub["bin"].values == i).astype(float)
                                 for i in range(len(K.BINS))])
            X = np.column_stack([D * early[:, None], D * (1 - early)[:, None], S])
            A = K.C_demean(np.column_stack([yy, X]), sub["winter"].values)
            b, *_ = np.linalg.lstsq(A[:, 1:], A[:, 0], rcond=None)
            if np.all(np.isfinite(b)):
                boot.append(b)
        except Exception:
            pass
    boot = np.array(boot)
    n = len(K.BINS)
    print(f"{len(boot)}/{N_BOOT} bootstrap fits converged\n")

    res = {"split_year": SPLIT, "n_boot": N_BOOT, "n_events_early": n_early,
           "n_events_late": int(len(on) - n_early), "bins": {}}
    print(f"{'bin':10s} {'pre-1980':>18s} {'1980+':>18s} {'difference':>22s}")
    print("-" * 72)
    for i, lab in enumerate(K.LABELS):
        de = boot[:, i] - boot[:, n + i]
        p = float(2 * min((de >= 0).mean(), (de <= 0).mean()))
        res["bins"][lab] = {
            "effect_pre1980": round(float(e0[i]), 4),
            "effect_1980on": round(float(l0[i]), 4),
            "difference": round(float(e0[i] - l0[i]), 4),
            "difference_CI95": [round(float(np.percentile(de, 2.5)), 4),
                                round(float(np.percentile(de, 97.5)), 4)],
            "difference_p": p}
        star = " *" if p < 0.05 else ""
        print(f"{lab:10s} {e0[i]:+18.3f} {l0[i]:+18.3f} "
              f"{e0[i]-l0[i]:+10.3f} (p={p:.3f}){star}")

    b = res["bins"]
    print("\n=== which explanation ===")
    pre_d = b["-30..-16"]["difference"]
    post_d = b["+15..+29"]["difference"]
    print(f"  pre-onset  era difference {pre_d:+.3f} (p={b['-30..-16']['difference_p']:.3f})")
    print(f"  post-onset era difference {post_d:+.3f} (p={b['+15..+29']['difference_p']:.3f})")
    print("  DATING predicts post-onset ALSO attenuated pre-1980 "
          "(response smeared, not lost);")
    print("  REAL   predicts post-onset comparable across eras.")
    res["interpretation_note"] = (
        "dating-error explanation requires the post-onset effect to be attenuated "
        "pre-1980 as well; a real change does not")
    (RESULTS / "era_split.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> era_split.json")


if __name__ == "__main__":
    main()
