#!/usr/bin/env python3
"""
isolation_interaction.py
========================
HONESTY CHECK on isolated_events.py.

That script showed the pre-onset anomaly at -0.55 (p=0.112) for isolated events
and -2.75 for clustered ones. It is tempting to read that as "clustering causes
the precursor". It does not establish that: the isolated set is SMALLER (29 vs
43), so losing significance there is partly lost power, and the two subsets were
never compared against each other.

This fits both in one model, with every bin interacted with CLUSTERED, so the
difference is an estimated coefficient with its own interval -- the same
treatment era_split.py gives the era question (which came back p=0.115, i.e. NOT
significant, and is reported as such).

  H0: the -30..-16 effect is the same for isolated and clustered events.

If this also fails to reach significance, then the correct claim is the weaker
one: the pre-onset anomaly is NOT ROBUST across design choices, not that either
design choice has been shown to cause it.

Output: isolation_interaction.json
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

N_BOOT = 4000
ISOLATION_DAYS = 120


def main():
    d = K.load_series("AO")
    on = load_catalogue("primary")
    on = on[(on >= d.index.min()) & (on <= d.index.max())]
    t = np.sort(on.values.astype("datetime64[D]").astype(int))
    gap = np.full(len(t), 10**6)
    gap[:-1] = np.minimum(gap[:-1], np.diff(t))
    gap[1:] = np.minimum(gap[1:], np.diff(t))
    iso_flag = gap > ISOLATION_DAYS

    # tag each retained day by whether its NEAREST onset is a clustered event
    dd = K.build(d, pd.DatetimeIndex(pd.to_datetime(t, unit="D")))
    dts = dd.index.values.astype("datetime64[D]").astype(int)
    nearest = np.abs(dts[:, None] - t[None, :]).argmin(axis=1)
    dd = dd.assign(clustered=(~iso_flag)[nearest].astype(float))
    print(f"{len(t)} events: {int(iso_flag.sum())} isolated, {int((~iso_flag).sum())} clustered")
    print(f"retained days: {int((dd['clustered']==0).sum())} nearest-isolated, "
          f"{int((dd['clustered']==1).sum())} nearest-clustered")

    def fit(df):
        y = df["y"].values.astype(float)
        S = C.harmonics(df["doy"].values, 3)
        cl = df["clustered"].values
        D = np.column_stack([(df["bin"].values == i).astype(float)
                             for i in range(len(K.BINS))])
        X = np.column_stack([D * (1 - cl)[:, None], D * cl[:, None], S])
        A = K.C_demean(np.column_stack([y, X]), df["winter"].values)
        b, *_ = np.linalg.lstsq(A[:, 1:], A[:, 0], rcond=None)
        return b

    b0 = fit(dd)
    n = len(K.BINS)
    rng = np.random.default_rng(20260730)
    winters = np.array(sorted(dd["winter"].unique()))
    idx = {w: np.flatnonzero((dd["winter"] == w).values) for w in winters}
    boot = []
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        sub = dd.iloc[np.concatenate([idx[w] for w in pick])].copy()
        sub["winter"] = np.concatenate([np.full(len(idx[w]), j) for j, w in enumerate(pick)])
        try:
            b = fit(sub)
            if np.all(np.isfinite(b)):
                boot.append(b)
        except Exception:
            pass
    boot = np.array(boot)
    print(f"{len(boot)}/{N_BOOT} bootstrap fits converged\n")

    res = {"n_boot": N_BOOT, "isolation_days": ISOLATION_DAYS, "bins": {}}
    print(f"{'bin':10s} {'isolated':>12s} {'clustered':>12s} {'difference':>24s}")
    print("-" * 64)
    for i, lab in enumerate(K.LABELS):
        diff = boot[:, i] - boot[:, n + i]
        p = float(2 * min((diff >= 0).mean(), (diff <= 0).mean()))
        res["bins"][lab] = {
            "effect_isolated": round(float(b0[i]), 4),
            "effect_clustered": round(float(b0[n + i]), 4),
            "difference": round(float(b0[i] - b0[n + i]), 4),
            "difference_CI95": [round(float(np.percentile(diff, 2.5)), 4),
                                round(float(np.percentile(diff, 97.5)), 4)],
            "difference_p": p}
        print(f"{lab:10s} {b0[i]:+12.3f} {b0[n+i]:+12.3f} "
              f"{b0[i]-b0[n+i]:+10.3f} (p={p:.3f}){' *' if p < 0.05 else ''}")

    pre = res["bins"]["-30..-16"]
    print(f"\n=== verdict on the clustering explanation ===")
    print(f"  pre-onset difference {pre['difference']:+.3f} "
          f"95% CI [{pre['difference_CI95'][0]:+.3f}, {pre['difference_CI95'][1]:+.3f}] "
          f"p={pre['difference_p']:.3f}")
    print("  significant -> clustering demonstrably inflates the precursor")
    print("  not significant -> the honest claim is only that the pre-onset")
    print("                     anomaly is NOT ROBUST, with cause unproven")
    (RESULTS / "isolation_interaction.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> isolation_interaction.json")


if __name__ == "__main__":
    main()
