"""
03_chemistry_response.py — Superposed Epoch Analysis of MLS Chemistry
======================================================================
Tests whether geomagnetic disturbance events produce measurable
perturbations in polar stratospheric chemistry (ozone, temperature).
This is the upstream manipulation check — if chemistry doesn't respond,
the causal chain cannot be established.
"""
import sys
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, FIGURES, LOG, load_panel

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def compute_anomalies(series: pd.Series, window: int = 61) -> pd.Series:
    """Compute anomalies relative to centered running mean (±30 day climatology)."""
    clim = series.rolling(window, center=True, min_periods=20).mean()
    return series - clim


def superposed_epoch(panel: pd.DataFrame, event_dates, variable: str,
                     before: int = 15, after: int = 30) -> dict:
    """
    Superposed epoch analysis centered on event dates.
    Returns epoch matrix (events × lag), mean profile, CI, and permutation p-value.
    """
    if variable not in panel.columns:
        return None

    # Compute anomalies
    anom = compute_anomalies(panel[variable].dropna())

    lags = np.arange(-before, after + 1)
    epoch_matrix = []

    for ed in event_dates:
        row = []
        for lag in lags:
            target = ed + pd.Timedelta(days=int(lag))
            if target in anom.index:
                row.append(anom.loc[target])
            else:
                row.append(np.nan)
        epoch_matrix.append(row)

    epoch = np.array(epoch_matrix)

    # Mean and SEM across events
    with np.errstate(all="ignore"):
        mean_profile = np.nanmean(epoch, axis=0)
        sem_profile = np.nanstd(epoch, axis=0, ddof=1) / np.sqrt(
            np.sum(~np.isnan(epoch), axis=0)
        )

    # Block permutation test: 2000 iterations
    n_perms = 2000
    n_events = len(event_dates)

    # Generate sham event sets (random winter NDJFM days, not near real events)
    winter_days = panel.index[
        (panel["is_winter"] == 1) & (panel["geo_event"] == 0) &
        (panel["post_event_0_30d"] == 0) & panel[variable].notna()
    ]

    perm_means = []
    rng = np.random.default_rng(42)
    for _ in range(n_perms):
        sham_dates = rng.choice(winter_days, size=min(n_events, len(winter_days)),
                                replace=False)
        sham_epoch = []
        for sd in sham_dates:
            row = []
            for lag in lags:
                target = sd + pd.Timedelta(days=int(lag))
                if target in anom.index:
                    row.append(anom.loc[target])
                else:
                    row.append(np.nan)
            sham_epoch.append(row)
        sham_arr = np.array(sham_epoch)
        perm_means.append(np.nanmean(sham_arr, axis=0))

    perm_means = np.array(perm_means)

    # P-values for each lag: fraction of permutations with more extreme mean
    p_values = np.zeros(len(lags))
    for j in range(len(lags)):
        observed = mean_profile[j]
        null_dist = perm_means[:, j]
        null_dist = null_dist[~np.isnan(null_dist)]
        if len(null_dist) > 0:
            # Two-sided: fraction more extreme than observed
            p_values[j] = np.mean(np.abs(null_dist) >= np.abs(observed))

    # Aggregate p-values for pre-specified windows
    window_results = {}
    for w_start, w_end in [(0, 7), (8, 14), (15, 21)]:
        idx = (lags >= w_start) & (lags <= w_end)
        window_mean = np.nanmean(mean_profile[idx])

        # Permutation test for window mean
        perm_window = np.nanmean(perm_means[:, idx], axis=1)
        perm_window = perm_window[~np.isnan(perm_window)]
        if len(perm_window) > 0:
            p_window = np.mean(np.abs(perm_window) >= np.abs(window_mean))
        else:
            p_window = 1.0

        window_results[f"d{w_start}_{w_end}"] = {
            "mean_anomaly": float(window_mean),
            "p_value": float(p_window),
            "significant": p_window < 0.05,
        }

    return {
        "variable": variable,
        "n_events": n_events,
        "lags": lags.tolist(),
        "mean_profile": [float(x) if not np.isnan(x) else None for x in mean_profile],
        "sem_profile": [float(x) if not np.isnan(x) else None for x in sem_profile],
        "p_values": [float(x) for x in p_values],
        "window_results": window_results,
        "n_permutations": n_perms,
    }


def main():
    panel = load_panel(winter_only=False)

    # Get winter event dates where MLS data exists
    events = panel[(panel["geo_event"] == 1) & (panel["is_winter"] == 1)]
    event_dates = events.index

    # Chemistry variables to test
    chem_vars = [
        "mls_o3_lev_1p0hpa", "mls_o3_lev_2p0hpa", "mls_o3_lev_4p6hpa",
        "mls_o3_lev_10p0hpa",
        "mls_t_lev_1p0hpa", "mls_t_lev_2p0hpa", "mls_t_lev_10p0hpa",
    ]

    # Filter to events with MLS coverage (2004+)
    mls_start = pd.Timestamp("2004-08-01", tz="UTC")
    event_dates_mls = event_dates[event_dates >= mls_start]
    LOG.info("Events with MLS coverage: %d / %d", len(event_dates_mls), len(event_dates))

    results = {}
    for var in chem_vars:
        LOG.info("SEA for %s ...", var)
        res = superposed_epoch(panel, event_dates_mls, var, before=15, after=30)
        if res is not None:
            results[var] = res
            wr = res["window_results"]
            LOG.info("  d0-7:  mean=%.4f p=%.3f %s",
                     wr["d0_7"]["mean_anomaly"], wr["d0_7"]["p_value"],
                     "***" if wr["d0_7"]["significant"] else "")
            LOG.info("  d8-14: mean=%.4f p=%.3f %s",
                     wr["d8_14"]["mean_anomaly"], wr["d8_14"]["p_value"],
                     "***" if wr["d8_14"]["significant"] else "")

    # Save results
    out = RESULTS / "chemistry_response.json"
    out.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    LOG.info("Results saved to %s", out)

    # Summary table
    print(f"\n{'='*70}")
    print("CHEMISTRY MANIPULATION CHECK — Superposed Epoch Analysis")
    print(f"{'='*70}")
    print(f"Events: {len(event_dates_mls)} winter geomagnetic disturbances (2004+)")
    print(f"Permutations: 2000")
    print(f"\n{'Variable':<25} {'Window':<10} {'Mean Anom':>10} {'p-value':>10} {'Sig':>5}")
    print("-" * 65)
    for var, res in results.items():
        short = var.replace("mls_o3_lev_", "O3 ").replace("mls_t_lev_", "T ")
        for wname, wr in res["window_results"].items():
            sig = "***" if wr["p_value"] < 0.01 else "**" if wr["p_value"] < 0.05 else "*" if wr["p_value"] < 0.1 else ""
            print(f"{short:<25} {wname:<10} {wr['mean_anomaly']:>10.4f} {wr['p_value']:>10.3f} {sig:>5}")

    # Decision gate
    any_sig = any(
        res["window_results"]["d0_7"]["significant"] or
        res["window_results"]["d8_14"]["significant"]
        for res in results.values()
    )
    print(f"\n{'='*70}")
    if any_sig:
        print("DECISION GATE: PASSED — Chemistry responds to geomagnetic disturbances")
    else:
        print("DECISION GATE: NOTE — Weak chemistry response. Consider:")
        print("  - Adjusting threshold definitions")
        print("  - Testing individual event subsets (strong storms only)")
        print("  - Examining seasonal stratification")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
