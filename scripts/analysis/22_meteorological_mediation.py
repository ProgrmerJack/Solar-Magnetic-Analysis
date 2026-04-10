"""
22_meteorological_mediation.py — Address Opus 4.5 Critic's Key Demands
=======================================================================
1. Meteorological intermediate analysis: NAO + precipitation around SSW events
2. Formal pre vs post SSW difference test (bootstrap the difference-of-differences)
3. Norwegian data classification check
4. Absolute rate reporting
"""

import json, logging, sys, gc, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis._analysis_utils import RESULTS, PROCESSED, FIGURES

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger("met_mediation")

def main():
    results = {}

    # ================================================================
    # Load panel
    # ================================================================
    panel = pd.read_parquet(PROCESSED / "analysis_panel_v2.parquet")
    winter = panel[panel["is_winter"] == 1].copy()
    LOG.info("Panel: %d winter days", len(winter))

    # ================================================================
    # PART 1: Absolute rate reporting
    # ================================================================
    LOG.info("=" * 60)
    LOG.info("PART 1: Absolute Avalanche Rates")
    LOG.info("=" * 60)

    geo_exposed = winter[winter["post_event_1_3d"] == 1]
    geo_unexposed = winter[winter["post_event_1_3d"] == 0]

    rates = {}
    for col, label in [("dry_natural_size_1234", "dry_natural"),
                        ("wet_natural_size_1234", "wet_natural"),
                        ("aai_all_natural", "all_natural"),
                        ("norway_aval_count", "norway")]:
        if col not in winter.columns:
            continue
        base = winter[col].mean()
        exp_rate = geo_exposed[col].mean() if len(geo_exposed) > 0 else np.nan
        unexp_rate = geo_unexposed[col].mean() if len(geo_unexposed) > 0 else np.nan
        rates[label] = {
            "overall_mean": round(base, 3),
            "overall_median": round(winter[col].median(), 1),
            "exposed_mean": round(exp_rate, 3),
            "unexposed_mean": round(unexp_rate, 3),
            "n_exposed_days": int(len(geo_exposed)),
            "n_unexposed_days": int(len(geo_unexposed)),
            "absolute_diff": round(exp_rate - unexp_rate, 3),
        }
        LOG.info("  %s: base=%.2f, exposed=%.2f, unexposed=%.2f, diff=%.3f",
                 label, base, exp_rate, unexp_rate, exp_rate - unexp_rate)

    results["part1_absolute_rates"] = rates
    _save(results)
    gc.collect()

    # ================================================================
    # PART 2: Meteorological intermediates around SSW events
    # ================================================================
    LOG.info("=" * 60)
    LOG.info("PART 2: Meteorological Intermediates Around SSW Events")
    LOG.info("=" * 60)

    ssw_cat = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw_cat.index = ssw_cat.index.tz_localize(None)

    met_cols = []
    for c in ["nao_daily", "era5_t2m_alps", "era5_precip_alps", "era5_u10_65N",
              "f107", "modis_snow_frac", "qbo_u50"]:
        if c in winter.columns:
            met_cols.append(c)
    LOG.info("  Available met columns: %s", met_cols)

    ssw_met = {}
    for _, row in ssw_cat.iterrows():
        ssw_date = row.name if hasattr(row.name, 'date') else pd.Timestamp(row.name)
        if ssw_date < winter.index.min() or ssw_date > winter.index.max():
            continue

        pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
        post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))

        pre_data = winter.loc[pre_mask]
        post_data = winter.loc[post_mask]

        if len(pre_data) == 0 or len(post_data) == 0:
            continue

        event_met = {"date": str(ssw_date.date()), "n_pre": len(pre_data), "n_post": len(post_data)}
        for c in met_cols:
            event_met[f"{c}_pre"] = round(pre_data[c].mean(), 4) if c in pre_data else np.nan
            event_met[f"{c}_post"] = round(post_data[c].mean(), 4) if c in post_data else np.nan
            event_met[f"{c}_diff"] = round(
                (post_data[c].mean() - pre_data[c].mean()) if c in post_data else np.nan, 4
            )

        ssw_met[str(ssw_date.date())] = event_met

    # Aggregate across SSW events
    met_summary = {}
    for c in met_cols:
        diffs = [v[f"{c}_diff"] for v in ssw_met.values()
                 if f"{c}_diff" in v and not np.isnan(v[f"{c}_diff"])]
        if len(diffs) >= 5:
            from scipy import stats
            t_stat, t_p = stats.ttest_1samp(diffs, 0)
            w_stat, w_p = stats.wilcoxon(diffs)
            n_neg = sum(1 for d in diffs if d < 0)
            met_summary[c] = {
                "n_events": len(diffs),
                "mean_diff": round(np.mean(diffs), 4),
                "median_diff": round(np.median(diffs), 4),
                "n_negative": n_neg,
                "t_stat": round(t_stat, 3),
                "t_p": round(t_p, 6),
                "wilcoxon_p": round(w_p, 6),
            }
            LOG.info("  %s: mean_diff=%.4f, %d/%d neg, t_p=%.4f, W_p=%.4f",
                     c, np.mean(diffs), n_neg, len(diffs), t_p, w_p)

    results["part2_met_intermediates"] = {
        "per_event": ssw_met,
        "summary": met_summary,
    }
    _save(results)
    gc.collect()

    # ================================================================
    # PART 3: Formal pre vs post SSW difference test
    # ================================================================
    LOG.info("=" * 60)
    LOG.info("PART 3: Formal Pre vs Post SSW Difference Test")
    LOG.info("=" * 60)

    from scipy import stats

    upg = json.load(open(RESULTS / "tier2_upgrade.json"))
    p2 = upg["part2_ssw_battery"]

    formal_tests = {}
    for outcome in ["dry_natural", "norway", "all_natural"]:
        d = p2[outcome]
        diffs = d.get("diffs", [])
        if not diffs or len(diffs) < 5:
            continue

        # Need pre-SSW diffs too — compute from panel
        pre_diffs = []
        post_diffs = []
        for _, row in ssw_cat.iterrows():
            ssw_date = row.name if hasattr(row.name, 'date') else pd.Timestamp(row.name)
            if ssw_date < winter.index.min() or ssw_date > winter.index.max():
                continue

            col_map = {
                "dry_natural": "dry_natural_size_1234",
                "all_natural": "aai_all_natural",
                "norway": "norway_aval_count",
            }
            col = col_map[outcome]
            if col not in winter.columns:
                continue

            # Post-SSW window (0-15d)
            post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
            post_val = winter.loc[post_mask, col].mean() if post_mask.sum() > 0 else np.nan

            # Pre-SSW window (15-0d before)
            pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
            pre_val = winter.loc[pre_mask, col].mean() if pre_mask.sum() > 0 else np.nan

            # Control: same calendar days in adjacent winters
            controls = []
            for offset in [-1, 1, -2, 2]:
                for window_type, dt_start, dt_end in [
                    ("post", ssw_date, ssw_date + pd.Timedelta(days=15)),
                    ("pre", ssw_date - pd.Timedelta(days=15), ssw_date),
                ]:
                    ctrl_start = dt_start + pd.DateOffset(years=offset)
                    ctrl_end = dt_end + pd.DateOffset(years=offset)
                    ctrl_mask = (winter.index >= ctrl_start) & (winter.index < ctrl_end)
                    if ctrl_mask.sum() > 5:
                        if window_type == "post":
                            controls.append(("post", winter.loc[ctrl_mask, col].mean()))
                        else:
                            controls.append(("pre", winter.loc[ctrl_mask, col].mean()))

            post_ctrls = [v for t, v in controls if t == "post"]
            pre_ctrls = [v for t, v in controls if t == "pre"]

            if not np.isnan(post_val) and post_ctrls:
                post_diffs.append(post_val - np.mean(post_ctrls))
            if not np.isnan(pre_val) and pre_ctrls:
                pre_diffs.append(pre_val - np.mean(pre_ctrls))

        if len(post_diffs) >= 5 and len(pre_diffs) >= 5:
            # Paired difference: post_diff - pre_diff
            n = min(len(post_diffs), len(pre_diffs))
            paired_diff = [post_diffs[i] - pre_diffs[i] for i in range(n)]

            t_stat, t_p = stats.ttest_1samp(paired_diff, 0)
            w_stat, w_p = stats.wilcoxon(paired_diff)

            # Bootstrap the difference
            np.random.seed(42)
            boot_diffs = []
            for _ in range(10000):
                idx = np.random.choice(n, n, replace=True)
                boot_diffs.append(np.mean([paired_diff[i] for i in idx]))
            ci_lo = np.percentile(boot_diffs, 2.5)
            ci_hi = np.percentile(boot_diffs, 97.5)

            formal_tests[outcome] = {
                "n_events": n,
                "mean_post_diff": round(np.mean(post_diffs[:n]), 3),
                "mean_pre_diff": round(np.mean(pre_diffs[:n]), 3),
                "mean_paired_diff": round(np.mean(paired_diff), 3),
                "t_stat": round(t_stat, 3),
                "t_p": round(t_p, 6),
                "wilcoxon_p": round(w_p, 6),
                "bootstrap_ci": [round(ci_lo, 3), round(ci_hi, 3)],
                "post_larger_than_pre": abs(np.mean(post_diffs[:n])) > abs(np.mean(pre_diffs[:n])),
            }
            LOG.info("  %s: post=%.3f, pre=%.3f, diff=%.3f, t_p=%.4f, CI=[%.3f, %.3f]",
                     outcome, np.mean(post_diffs[:n]), np.mean(pre_diffs[:n]),
                     np.mean(paired_diff), t_p, ci_lo, ci_hi)

    results["part3_pre_vs_post_formal"] = formal_tests
    _save(results)
    gc.collect()

    # ================================================================
    # PART 4: Norwegian data classification check
    # ================================================================
    LOG.info("=" * 60)
    LOG.info("PART 4: Norwegian Data Classification")
    LOG.info("=" * 60)

    norway_info = {}
    norway_cols = [c for c in winter.columns if "norway" in c.lower()]
    LOG.info("  Norway columns: %s", norway_cols)
    for c in norway_cols:
        norway_info[c] = {
            "mean": round(winter[c].mean(), 3),
            "sum": round(winter[c].sum(), 1),
            "n_nonzero": int((winter[c] > 0).sum()),
        }
        LOG.info("  %s: mean=%.3f, sum=%.0f, n_nonzero=%d",
                 c, winter[c].mean(), winter[c].sum(), (winter[c] > 0).sum())

    results["part4_norway_classification"] = norway_info
    _save(results)

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 60)
    print("METEOROLOGICAL MEDIATION ANALYSIS COMPLETE")
    print("=" * 60)
    for k in results:
        print(f"  {k}: OK")


def _save(results):
    out = RESULTS / "met_mediation.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2, default=str)


if __name__ == "__main__":
    main()
