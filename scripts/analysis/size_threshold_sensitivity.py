#!/usr/bin/env python3
"""
Size-threshold sensitivity analysis: compute event-level RR for size>=1 vs size>=2.
Addresses reviewer concern that size-1 events may include inconsistently reported small sloughs.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load raw observations
obs = pd.read_csv(os.path.join(ROOT, "data", "cryosphere", "davos_avalanches", "avalanche_observations.csv"), sep=";")
obs["date"] = pd.to_datetime(obs["date_release"])

# Load SSW catalog
ssw = pd.read_csv(os.path.join(ROOT, "data", "results", "ssw_event_catalog.csv"))
ssw["onset"] = pd.to_datetime(ssw["date"])

print(f"Total observations: {len(obs)}")
print(f"SSW events: {len(ssw)}")

results = {}
for size_min in [1, 2, 3]:
    subset = obs[
        (obs["snow_type"] == "dry")
        & (obs["trigger_type"] == "NATURAL")
        & (obs["aval_size_class"] >= size_min)
    ]
    daily = subset.groupby("date").size().reset_index(name="count")
    print(f"\nSize >= {size_min}: {len(subset)} avalanches, {len(daily)} days with observations")

    # Build date frame spanning study period
    all_dates = pd.date_range("1998-11-01", "2019-04-30", freq="D")
    date_df = pd.DataFrame({"date": all_dates})
    date_df = date_df.merge(daily, on="date", how="left").fillna({"count": 0})
    date_df["doy"] = date_df["date"].dt.dayofyear
    date_df["month"] = date_df["date"].dt.month
    date_df["is_winter"] = date_df["month"].isin([11, 12, 1, 2, 3, 4])

    # Flag SSW windows
    date_df["in_ssw"] = False
    for _, ev in ssw.iterrows():
        onset = ev["onset"]
        mask = (date_df["date"] >= onset - pd.Timedelta(days=15)) & (
            date_df["date"] <= onset + pd.Timedelta(days=15)
        )
        date_df.loc[mask, "in_ssw"] = True

    # Compute event-level RR
    rrs = []
    for _, ev in ssw.iterrows():
        onset = ev["onset"]
        window = date_df[
            (date_df["date"] >= onset - pd.Timedelta(days=15))
            & (date_df["date"] <= onset + pd.Timedelta(days=15))
        ]
        observed = window["count"].sum()

        # DOY-matched expected from non-SSW winters
        expected = 0
        for _, day in window.iterrows():
            d = day["doy"]
            ref = date_df[
                (~date_df["in_ssw"])
                & (date_df["is_winter"])
                & (date_df["doy"].between(d - 3, d + 3))
            ]
            expected += ref["count"].mean()

        rr = observed / expected if expected > 0 else np.nan
        rrs.append(rr)

    rrs = np.array(rrs)
    log_rrs = np.log(rrs[rrs > 0])

    n_below = int(np.sum(rrs < 1))
    n_total = len(rrs)
    sign_p = stats.binomtest(n_below, n_total, 0.5, alternative="greater").pvalue
    geom_mean = float(np.exp(np.mean(log_rrs)))
    t_stat, t_p = stats.ttest_1samp(log_rrs, 0)
    wilcox_stat, wilcox_p = stats.wilcoxon(log_rrs, alternative="less")

    # Bootstrap 95% CI for geometric mean
    boot_means = []
    rng = np.random.default_rng(42)
    for _ in range(10000):
        sample = rng.choice(log_rrs, size=len(log_rrs), replace=True)
        boot_means.append(np.exp(np.mean(sample)))
    ci_lo, ci_hi = np.percentile(boot_means, [2.5, 97.5])

    print(f"  Events below 1: {n_below}/{n_total}")
    print(f"  Sign test P: {sign_p:.4f}")
    print(f"  Geometric mean RR: {geom_mean:.3f} [{ci_lo:.3f}, {ci_hi:.3f}]")
    print(f"  Wilcoxon P: {wilcox_p:.4f}")
    print(f"  t-test P: {t_p:.4f}")
    print(f"  RRs: {[round(r, 3) for r in rrs]}")

    results[f"size_ge_{size_min}"] = {
        "n_avalanches": int(len(subset)),
        "n_below_1": n_below,
        "n_total": n_total,
        "sign_p": float(sign_p),
        "geom_mean_rr": float(geom_mean),
        "ci_95_lo": float(ci_lo),
        "ci_95_hi": float(ci_hi),
        "wilcoxon_p": float(wilcox_p),
        "t_p": float(t_p),
        "rrs": [float(r) for r in rrs],
    }

# Save results
outpath = os.path.join(ROOT, "data", "results", "size_threshold_sensitivity.json")
with open(outpath, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {outpath}")
