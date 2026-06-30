#!/usr/bin/env python3
"""
r103_monitoring_blindspot.py
============================
Addresses the two NG "killers" by relocating the decisive + novel claim away from
the unconfirmable compound-event TYPE toward the MONITORING BLIND SPOT, which is
(a) decisively testable and (b) broadly significant (a generalizable failure mode
of activity-based hazard monitoring, not "SSW -> cold").

PART A - Count-danger decoupling (the blind spot):
  Normally forecaster danger tracks observable natural avalanche activity. We fit
  that mapping danger = f(natural_count) on NON-SSW winter days (isotonic), then ask
  whether, during SSW windows, observed danger EXCEEDS what the observable activity
  predicts. A positive residual = hazard that count-based monitoring cannot see.
  Tested per event (n=16) and as a magnitude.

PART B - Natural/human trigger-fraction shift (Davos individual avalanches):
  Across 13,918 Davos avalanches (1998-2019) with trigger type, does the NATURAL
  fraction fall and the HUMAN fraction rise during SSW windows vs DOY-matched controls?

Outputs: data/results/r103_monitoring_blindspot.json
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.isotonic import IsotonicRegression

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r103_monitoring_blindspot.json"
WIN = 15


def load_panel():
    p = pd.read_parquet(ROOT / "data/processed/analysis_panel.parquet")
    if p.index.tz is not None:
        p.index = p.index.tz_convert("UTC").tz_localize(None)
    m = (p["is_winter"].astype(bool) & p["dry_natural_size_1234"].notna()
         & p["max_danger_corr"].notna())
    return p[m]


def part_a_blindspot(p):
    """Does danger exceed count-predicted danger during SSW windows?"""
    nonssw = p[p["ssw_within_15d"] == 0]
    ssw = p[p["ssw_within_15d"] == 1]
    x_nz = nonssw["dry_natural_size_1234"].values.astype(float)
    y_nz = nonssw["max_danger_corr"].values.astype(float)
    # monotone non-SSW mapping count -> expected danger
    iso = IsotonicRegression(out_of_bounds="clip").fit(x_nz, y_nz)

    p = p.copy()
    p["danger_pred"] = iso.predict(p["dry_natural_size_1234"].values.astype(float))
    p["residual"] = p["max_danger_corr"] - p["danger_pred"]

    # per-event residual (event = SSW window), winter id for clustering
    cat = pd.read_csv(ROOT / "data/results/ssw_event_catalog.csv")
    onsets = pd.to_datetime(cat["date"]).dt.tz_localize(None)
    per_event = []
    for d in onsets:
        w = p[(p.index >= d - pd.Timedelta(days=WIN)) & (p.index <= d + pd.Timedelta(days=WIN))]
        if len(w) == 0:
            continue
        per_event.append({"onset": str(d.date()),
                          "mean_residual": float(w["residual"].mean()),
                          "mean_danger": float(w["max_danger_corr"].mean()),
                          "mean_count": float(w["dry_natural_size_1234"].mean()),
                          "n_days": int(len(w))})
    res = np.array([e["mean_residual"] for e in per_event])
    n = len(res)
    n_pos = int((res > 0).sum())
    sign_p = stats.binomtest(n_pos, n, 0.5, alternative="greater").pvalue
    t, t_p = stats.ttest_1samp(res, 0)
    # descriptive high-N magnitude (clustering acknowledged)
    ssw_resid = ssw["max_danger_corr"].values - iso.predict(ssw["dry_natural_size_1234"].values.astype(float))
    return {
        "design": "danger residual vs count-predicted danger (non-SSW isotonic mapping)",
        "n_events": n,
        "events_with_hidden_hazard": f"{n_pos}/{n}",
        "sign_test_p_one_sided": float(sign_p),
        "mean_event_residual_danger_levels": float(np.mean(res)),
        "t_stat": float(t), "t_p_two_sided": float(t_p),
        "per_event": per_event,
        "highN_mean_residual_all_ssw_days": float(np.mean(ssw_resid)),
        "highN_n_ssw_days": int(len(ssw_resid)),
        "interpretation": (
            "Positive residual = forecaster-assessed danger exceeds what observable "
            "natural-avalanche counts predict -> hazard invisible to count-based "
            "monitoring. This is the predictable blind spot."),
    }


def part_b_davos_triggers(p):
    """Natural/human trigger fraction in SSW windows vs DOY-matched controls."""
    dv = pd.read_csv(ROOT / "data/cryosphere/davos_avalanches/avalanche_observations.csv",
                     sep=None, engine="python")
    dv["date"] = pd.to_datetime(dv["date_release"], errors="coerce")
    dv = dv.dropna(subset=["date"])
    dv["nat"] = (dv["trigger_type"] == "NATURAL").astype(int)
    dv["hum"] = (dv["trigger_type"] == "HUMAN").astype(int)
    dv = dv.set_index("date").sort_index()

    cat = pd.read_csv(ROOT / "data/results/ssw_event_catalog.csv")
    onsets = pd.to_datetime(cat["date"]).dt.tz_localize(None)
    ssw_years = set()
    for d in onsets:
        ssw_years |= {d.year, (d + pd.Timedelta(days=WIN)).year, (d - pd.Timedelta(days=WIN)).year}

    rows = []
    for d in onsets:
        win = dv[(dv.index >= d - pd.Timedelta(days=WIN)) & (dv.index <= d + pd.Timedelta(days=WIN))]
        if len(win) < 5:
            continue
        # DOY-matched control: same calendar window across non-SSW years
        doys = [(d.dayofyear + k - 1) % 366 + 1 for k in range(-WIN, WIN + 1)]
        lo, hi = min(doys), max(doys)
        idx = dv.index
        if hi - lo <= 2 * WIN + 2:
            cm = (idx.dayofyear >= lo) & (idx.dayofyear <= hi)
        else:
            cm = (idx.dayofyear >= lo) | (idx.dayofyear <= hi)
        ctrl = dv[cm & (~idx.year.isin(ssw_years))]
        if len(ctrl) < 5:
            continue
        rows.append({
            "onset": str(d.date()),
            "ssw_nat_frac": float(win["nat"].mean()), "ctrl_nat_frac": float(ctrl["nat"].mean()),
            "ssw_hum_frac": float(win["hum"].mean()), "ctrl_hum_frac": float(ctrl["hum"].mean()),
            "n_ssw": int(len(win)), "n_ctrl": int(len(ctrl)),
        })
    dnat = np.array([r["ssw_nat_frac"] - r["ctrl_nat_frac"] for r in rows])
    dhum = np.array([r["ssw_hum_frac"] - r["ctrl_hum_frac"] for r in rows])
    n = len(rows)
    return {
        "n_events": n, "n_davos_avalanches": int(len(dv)),
        "natural_fraction_down": {
            "events_lower": f"{int((dnat < 0).sum())}/{n}",
            "mean_delta": float(np.mean(dnat)),
            "sign_p_one_sided": float(stats.binomtest(int((dnat < 0).sum()), n, 0.5, alternative="greater").pvalue),
            "t_p": float(stats.ttest_1samp(dnat, 0).pvalue)},
        "human_fraction_up": {
            "events_higher": f"{int((dhum > 0).sum())}/{n}",
            "mean_delta": float(np.mean(dhum)),
            "sign_p_one_sided": float(stats.binomtest(int((dhum > 0).sum()), n, 0.5, alternative="greater").pvalue),
            "t_p": float(stats.ttest_1samp(dhum, 0).pvalue)},
        "per_event": rows,
    }


def main():
    p = load_panel()
    print(f"Panel winter-days with count+danger: {len(p)}")
    A = part_a_blindspot(p)
    print("\n=== PART A: monitoring blind spot ===")
    print(f"  events with hidden hazard (danger > count-predicted): {A['events_with_hidden_hazard']}")
    print(f"  mean residual: {A['mean_event_residual_danger_levels']:+.3f} danger levels; "
          f"sign P={A['sign_test_p_one_sided']:.4f}; t-P={A['t_p_two_sided']:.4f}")
    B = part_b_davos_triggers(p)
    print("\n=== PART B: Davos trigger-fraction dissociation ===")
    print(f"  natural fraction down: {B['natural_fraction_down']['events_lower']}, "
          f"mean Δ={B['natural_fraction_down']['mean_delta']:+.3f}, "
          f"sign P={B['natural_fraction_down']['sign_p_one_sided']:.4f}")
    print(f"  human fraction up:    {B['human_fraction_up']['events_higher']}, "
          f"mean Δ={B['human_fraction_up']['mean_delta']:+.3f}, "
          f"sign P={B['human_fraction_up']['sign_p_one_sided']:.4f}")

    out = {"part_a_blindspot": A, "part_b_davos_triggers": B}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
