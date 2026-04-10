"""
24_fresh_analysis.py — Comprehensive fresh statistical analysis for Nature Geoscience paper.
Run parts sequentially to avoid OOM:
    python 24_fresh_analysis.py --part 1
    python 24_fresh_analysis.py --part 2
    ...through --part 7

Part 1: Primary SSW matched comparison (Swiss dry + Norwegian)
Part 2: Type specificity (MH case-crossover) + geomagnetic storms
Part 3: Temporal structure (multi-window + event-study distributed lag)
Part 4: Meteorological chain (4-step) + NAO + U850 direct test
Part 5: Norwegian cold-regime + SNOTEL + prediction model
Part 6: Robustness battery (specification curve + LOOCV)
Part 7: Publication-quality figures for LaTeX
"""
import sys, json, warnings, pathlib, gc
import numpy as np
import pandas as pd
from scipy import stats
warnings.filterwarnings('ignore')

REPO = pathlib.Path(__file__).resolve().parents[2]
PROCESSED = REPO / 'data' / 'processed'
RESULTS   = REPO / 'data' / 'results'
FIGURES   = REPO / 'data' / 'figures'
RESULTS.mkdir(exist_ok=True, parents=True)
FIGURES.mkdir(exist_ok=True, parents=True)

def get_part():
    for i, a in enumerate(sys.argv[1:], 1):
        if a == '--part' and i < len(sys.argv):
            return int(sys.argv[i+1])
    return 1

def load_data():
    panel = pd.read_parquet(PROCESSED / 'analysis_panel_v2.parquet')
    ssw = pd.read_parquet(PROCESSED / 'atmospheric' / 'ssw_catalog.parquet')
    ssw.index = pd.to_datetime(ssw.index).tz_localize(None)
    winter = panel[panel['is_winter'] == 1].copy()
    ssw_dates = ssw.index[(ssw.index >= winter.index.min()) & (ssw.index <= winter.index.max())]
    print(f"Panel: {len(panel)} rows | Winter: {len(winter)} rows | SSW events: {len(ssw_dates)}")
    return panel, winter, ssw, ssw_dates

def get_winter_year(dt):
    return dt.year if dt.month < 7 else dt.year + 1

def get_ssw_winter_ids(ssw_dates):
    """Get the winter_id strings for each SSW date."""
    ids = set()
    for sd in ssw_dates:
        if sd.month >= 11:
            ids.add(f"{sd.year}/{sd.year+1}")
        elif sd.month <= 6:
            ids.add(f"{sd.year-1}/{sd.year}")
    return ids

def get_matched_control(ssw_date, ssw_dates, winter, col, window_start_offset, window_days):
    """
    Match using ALL non-SSW winters with day-of-season ±3 matching.
    This is the robust matching used in the original analysis (script 16).
    """
    w_start = ssw_date + pd.Timedelta(days=window_start_offset)
    w_end = w_start + pd.Timedelta(days=window_days)
    mask = (winter.index >= w_start) & (winter.index < w_end)
    event_val = winter.loc[mask, col].dropna()
    if len(event_val) == 0:
        return np.nan, np.nan
    event_mean = event_val.mean()

    # Get day_of_season for SSW event window
    if len(event_val) > 0:
        dos_vals = winter.loc[mask, 'day_of_season'].values
        dos_lo = dos_vals.min() - 3
        dos_hi = dos_vals.max() + 3
    else:
        return event_mean, np.nan

    ssw_wids = get_ssw_winter_ids(ssw_dates)
    all_wids = winter['winter_id'].dropna().unique()
    non_ssw_wids = [w for w in all_wids if w not in ssw_wids]

    ctrl_means = []
    for wid in non_ssw_wids:
        ctrl_winter = winter[winter['winter_id'] == wid]
        matched = ctrl_winter[(ctrl_winter['day_of_season'] >= dos_lo) &
                              (ctrl_winter['day_of_season'] <= dos_hi)]
        if len(matched) >= 5:
            cmean = matched[col].dropna().mean()
            if not np.isnan(cmean):
                ctrl_means.append(cmean)

    if len(ctrl_means) < 3:
        return event_mean, np.nan
    return event_mean, np.mean(ctrl_means)

def exact_signflip_pvalue(diffs):
    """Exact sign-flip permutation test (2^n permutations)."""
    diffs = np.array([d for d in diffs if not np.isnan(d)])
    n = len(diffs)
    if n > 20:
        # Monte Carlo for large n
        obs = np.abs(np.mean(diffs))
        count = 0
        nperms = 100000
        for _ in range(nperms):
            signs = np.random.choice([-1, 1], size=n)
            if np.abs(np.mean(diffs * signs)) >= obs:
                count += 1
        return count / nperms
    obs = np.abs(np.mean(diffs))
    count = 0
    total = 2**n
    for i in range(total):
        signs = np.array([1 if (i >> j) & 1 else -1 for j in range(n)])
        if np.abs(np.mean(diffs * signs)) >= obs:
            count += 1
    return count / total

def save_json(data, name):
    path = RESULTS / name
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, default=str)
    print(f"  Saved: {path}")

part = get_part()
print(f"\n{'='*70}")
print(f"  PART {part}: COMPREHENSIVE FRESH ANALYSIS")
print(f"{'='*70}\n")

# ================================================================
# PART 1: PRIMARY SSW MATCHED COMPARISON
# ================================================================
if part == 1:
    panel, winter, ssw, ssw_dates = load_data()
    results = {"part": 1, "description": "Primary SSW matched comparison"}

    for col, label in [('dry_natural_size_1234', 'swiss_dry_slab'),
                       ('norway_aval_count', 'norway_total')]:
        print(f"\n--- {label} ---")
        diffs = []
        event_details = []

        for ssw_date in ssw_dates:
            ev, ctrl = get_matched_control(ssw_date, ssw_dates, winter, col, 0, 15)
            if np.isnan(ev) or np.isnan(ctrl):
                continue
            diff = ev - ctrl
            diffs.append(diff)
            event_details.append({
                "date": str(ssw_date.date()),
                "event_mean": round(ev, 3),
                "control_mean": round(ctrl, 3),
                "difference": round(diff, 3)
            })

        diffs_arr = np.array(diffs)
        n = len(diffs_arr)
        n_neg = int(np.sum(diffs_arr < 0))
        mean_d = float(np.mean(diffs_arr))
        median_d = float(np.median(diffs_arr))

        # Four statistical tests
        t_stat, t_p = stats.ttest_1samp(diffs_arr, 0)
        sign_p = float(stats.binomtest(n_neg, n, 0.5).pvalue)
        wilcox_stat, wilcox_p = stats.wilcoxon(diffs_arr, alternative='two-sided')
        perm_p = exact_signflip_pvalue(diffs_arr)

        # Bootstrap 95% CI
        np.random.seed(42)
        boot_means = [np.mean(np.random.choice(diffs_arr, size=n, replace=True)) for _ in range(10000)]
        ci_lo, ci_hi = float(np.percentile(boot_means, 2.5)), float(np.percentile(boot_means, 97.5))

        # Effect size (Cohen's d)
        cohens_d = float(mean_d / np.std(diffs_arr, ddof=1))

        print(f"  N events: {n}")
        print(f"  N negative: {n_neg}/{n}")
        print(f"  Mean diff: {mean_d:.3f} (median: {median_d:.3f})")
        print(f"  Bootstrap 95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]")
        print(f"  Paired t: P={t_p:.6f}")
        print(f"  Sign test: P={sign_p:.6f}")
        print(f"  Wilcoxon: P={wilcox_p:.6f}")
        print(f"  Exact perm: P={perm_p:.6f}")
        print(f"  Cohen's d: {cohens_d:.3f}")

        # Leave-one-out
        loo = []
        for i in range(n):
            sub = np.delete(diffs_arr, i)
            _, p_loo = stats.ttest_1samp(sub, 0)
            loo.append({"dropped": event_details[i]["date"], "p": round(float(p_loo), 6)})

        results[label] = {
            "n_events": n,
            "n_negative": n_neg,
            "mean_diff": round(mean_d, 4),
            "median_diff": round(median_d, 4),
            "bootstrap_ci_lo": round(ci_lo, 4),
            "bootstrap_ci_hi": round(ci_hi, 4),
            "ttest_p": round(float(t_p), 6),
            "sign_p": round(float(sign_p), 6),
            "wilcoxon_p": round(float(wilcox_p), 6),
            "perm_p": round(float(perm_p), 6),
            "cohens_d": round(cohens_d, 4),
            "events": event_details,
            "leave_one_out": loo
        }

    save_json(results, "fresh_part1_ssw_primary.json")
    print("\n[PART 1 COMPLETE]")

# ================================================================
# PART 2: TYPE SPECIFICITY + GEOMAGNETIC CASE-CROSSOVER
# ================================================================
elif part == 2:
    panel, winter, ssw, ssw_dates = load_data()
    results = {"part": 2, "description": "Type specificity and geomagnetic storm case-crossover"}

    # Mantel-Haenszel rate ratio computation
    def mantel_haenszel_rr(exposed_events, exposed_days, unexposed_events, unexposed_days):
        """Compute MH rate ratio across strata with RBG variance."""
        num = 0; den = 0
        var_num = 0
        for e1, d1, e0, d0 in zip(exposed_events, exposed_days, unexposed_events, unexposed_days):
            T = d1 + d0
            if T == 0: continue
            num += e1 * d0 / T
            den += e0 * d1 / T
            var_num += (e1 * d0**2 + e0 * d1**2) / T**2  # RBG variance
        if den == 0:
            return np.nan, (np.nan, np.nan), np.nan, np.nan
        rr = num / den
        se_ln_rr = np.sqrt(var_num / (num * den))
        ci_lo = np.exp(np.log(rr) - 1.96 * se_ln_rr)
        ci_hi = np.exp(np.log(rr) + 1.96 * se_ln_rr)
        z = np.log(rr) / se_ln_rr
        p = 2 * stats.norm.sf(np.abs(z))
        return float(rr), (float(ci_lo), float(ci_hi)), float(z), float(p)

    # Stratify by winter_id x 15-day period
    winter['period'] = winter['day_of_season'] // 15
    winter['stratum'] = winter['winter_id'].astype(str) + '_' + winter['period'].astype(str)
    strata = winter.groupby('stratum')

    for col, label in [('dry_natural_size_1234', 'dry_natural'),
                       ('wet_natural_size_1234', 'wet_natural'),
                       ('natural_size_1234', 'all_natural')]:
        print(f"\n--- MH: {label} vs geomag 1-3d ---")
        exp_e, exp_d, unexp_e, unexp_d = [], [], [], []
        for _, g in strata:
            exposed = g['post_event_1_3d'] == 1
            n_exp = exposed.sum()
            n_unexp = (~exposed).sum()
            if n_exp == 0 or n_unexp == 0:
                continue
            exp_e.append(g.loc[exposed, col].sum())
            exp_d.append(n_exp)
            unexp_e.append(g.loc[~exposed, col].sum())
            unexp_d.append(n_unexp)

        rr, ci, z, p = mantel_haenszel_rr(exp_e, exp_d, unexp_e, unexp_d)
        n_strata = len(exp_e)

        # Absolute rates
        total_exp_events = sum(exp_e)
        total_exp_days = sum(exp_d)
        total_unexp_events = sum(unexp_e)
        total_unexp_days = sum(unexp_d)
        rate_exp = total_exp_events / total_exp_days if total_exp_days > 0 else 0
        rate_unexp = total_unexp_events / total_unexp_days if total_unexp_days > 0 else 0

        print(f"  RR={rr:.3f}, 95% CI [{ci[0]:.3f}, {ci[1]:.3f}], z={z:.2f}, P={p:.6f}")
        print(f"  Rate exposed: {rate_exp:.3f}/day, unexposed: {rate_unexp:.3f}/day")
        print(f"  N strata: {n_strata}")

        results[f"mh_{label}_geomag"] = {
            "rate_ratio": round(rr, 4),
            "ci_lo": round(ci[0], 4),
            "ci_hi": round(ci[1], 4),
            "z_stat": round(z, 3),
            "p_value": round(p, 6),
            "n_strata": n_strata,
            "rate_exposed": round(rate_exp, 4),
            "rate_unexposed": round(rate_unexp, 4),
            "total_exposed_events": float(total_exp_events),
            "total_exposed_days": int(total_exp_days),
            "total_unexposed_events": float(total_unexp_events),
            "total_unexposed_days": int(total_unexp_days)
        }

    # Isolated events stratum-width sensitivity
    print("\n--- Isolated event stratum-width sensitivity ---")
    geo_events = winter[winter['geo_event'] == 1].index
    # Decluster: keep events >14d apart
    sorted_events = sorted(geo_events)
    isolated = []
    for i, ev in enumerate(sorted_events):
        is_iso = True
        for j, other in enumerate(sorted_events):
            if i != j and abs((ev - other).days) <= 14:
                is_iso = False
                break
        if is_iso:
            isolated.append(ev)
    print(f"  Isolated events: {len(isolated)} of {len(sorted_events)}")

    sensitivity = []
    for sw in [7, 10, 15, 21]:
        winter_copy = winter.copy()
        winter_copy['period_sw'] = winter_copy['day_of_season'] // sw
        winter_copy['stratum_sw'] = winter_copy['winter_id'].astype(str) + '_' + winter_copy['period_sw'].astype(str)

        for direction, day_range in [('post', (1, 3)), ('pre', (-7, -3))]:
            # Mark exposure
            winter_copy['exposed_dir'] = 0
            for ev in isolated:
                lo = ev + pd.Timedelta(days=day_range[0])
                hi = ev + pd.Timedelta(days=day_range[1])
                mask = (winter_copy.index >= lo) & (winter_copy.index <= hi)
                winter_copy.loc[mask, 'exposed_dir'] = 1

            strata_sw = winter_copy.groupby('stratum_sw')
            exp_e, exp_d, unexp_e, unexp_d = [], [], [], []
            for _, g in strata_sw:
                exposed = g['exposed_dir'] == 1
                n_exp = exposed.sum()
                n_unexp = (~exposed).sum()
                if n_exp == 0 or n_unexp == 0:
                    continue
                exp_e.append(g.loc[exposed, 'dry_natural_size_1234'].sum())
                exp_d.append(n_exp)
                unexp_e.append(g.loc[~exposed, 'dry_natural_size_1234'].sum())
                unexp_d.append(n_unexp)

            rr, ci, z, p = mantel_haenszel_rr(exp_e, exp_d, unexp_e, unexp_d)
            print(f"  SW={sw}d {direction}: RR={rr:.3f}, P={p:.6f}")
            sensitivity.append({
                "stratum_width": sw,
                "direction": direction,
                "rr": round(rr, 4),
                "ci_lo": round(ci[0], 4),
                "ci_hi": round(ci[1], 4),
                "p": round(p, 6)
            })

    results["isolated_sensitivity"] = sensitivity
    results["n_isolated_events"] = len(isolated)

    # Falsification: summer months
    summer = panel[panel['is_summer'] == 1].copy()
    summer['period'] = (summer['day_of_year'] - 120) // 15
    summer['stratum'] = summer.index.year.astype(str) + '_' + summer['period'].astype(str)
    summer_strata = summer.groupby('stratum')
    exp_e, exp_d, unexp_e, unexp_d = [], [], [], []
    for _, g in summer_strata:
        exposed = g['post_event_1_3d'] == 1
        n_exp = exposed.sum()
        n_unexp = (~exposed).sum()
        if n_exp == 0 or n_unexp == 0:
            continue
        exp_e.append(g.loc[exposed, 'dry_natural_size_1234'].sum())
        exp_d.append(n_exp)
        unexp_e.append(g.loc[~exposed, 'dry_natural_size_1234'].sum())
        unexp_d.append(n_unexp)
    rr_s, ci_s, z_s, p_s = mantel_haenszel_rr(exp_e, exp_d, unexp_e, unexp_d)
    print(f"\n  Summer falsification: RR={rr_s:.3f}, P={p_s:.4f}")
    results["summer_falsification"] = {
        "rr": round(rr_s, 4), "p": round(p_s, 6)
    }

    # Placebo tests (10 random catalogs)
    print("\n--- Placebo test (10 random catalogs) ---")
    np.random.seed(123)
    placebo_results = []
    for trial in range(10):
        fake_events = np.random.choice(winter.index, size=len(geo_events), replace=False)
        winter_tmp = winter.copy()
        winter_tmp['fake_exp'] = 0
        for fe in fake_events:
            lo = fe + pd.Timedelta(days=1)
            hi = fe + pd.Timedelta(days=3)
            mask = (winter_tmp.index >= lo) & (winter_tmp.index <= hi)
            winter_tmp.loc[mask, 'fake_exp'] = 1
        winter_tmp['stratum_p'] = winter_tmp['winter_id'].astype(str) + '_' + (winter_tmp['day_of_season'] // 15).astype(str)
        strata_p = winter_tmp.groupby('stratum_p')
        exp_e, exp_d, unexp_e, unexp_d = [], [], [], []
        for _, g in strata_p:
            exposed = g['fake_exp'] == 1
            n_exp = exposed.sum()
            n_unexp = (~exposed).sum()
            if n_exp == 0 or n_unexp == 0:
                continue
            exp_e.append(g.loc[exposed, 'dry_natural_size_1234'].sum())
            exp_d.append(n_exp)
            unexp_e.append(g.loc[~exposed, 'dry_natural_size_1234'].sum())
            unexp_d.append(n_unexp)
        rr_pl, _, _, p_pl = mantel_haenszel_rr(exp_e, exp_d, unexp_e, unexp_d)
        print(f"  Placebo {trial+1}: RR={rr_pl:.3f}, P={p_pl:.4f}")
        placebo_results.append({"trial": trial+1, "rr": round(rr_pl, 4), "p": round(p_pl, 6)})
    results["placebo_tests"] = placebo_results
    n_sig_placebo = sum(1 for pr in placebo_results if pr["p"] < 0.05)
    print(f"  Placebo significant: {n_sig_placebo}/10")
    results["placebo_n_significant"] = n_sig_placebo

    save_json(results, "fresh_part2_type_specificity.json")
    print("\n[PART 2 COMPLETE]")

# ================================================================
# PART 3: TEMPORAL STRUCTURE + EVENT-STUDY
# ================================================================
elif part == 3:
    panel, winter, ssw, ssw_dates = load_data()
    results = {"part": 3, "description": "Temporal structure and event-study"}

    # Multi-window temporal structure
    windows = [
        ("pre_15_0d", -15, 15),
        ("post_0_7d", 0, 7),
        ("post_0_15d", 0, 15),
        ("post_0_30d", 0, 30),
        ("post_15_30d", 15, 15),
    ]

    for col, label in [('dry_natural_size_1234', 'swiss_dry'),
                       ('norway_aval_count', 'norway')]:
        print(f"\n--- Temporal windows: {label} ---")
        window_results = []
        for wname, offset, days in windows:
            diffs = []
            for ssw_date in ssw_dates:
                ev, ctrl = get_matched_control(ssw_date, ssw_dates, winter, col, offset, days)
                if not np.isnan(ev) and not np.isnan(ctrl):
                    diffs.append(ev - ctrl)
            diffs = np.array(diffs)
            n_neg = int(np.sum(diffs < 0))
            mean_d = float(np.mean(diffs))
            t_stat, t_p = stats.ttest_1samp(diffs, 0)
            _, w_p = stats.wilcoxon(diffs) if len(diffs) > 5 else (0, 1)
            sign_p = float(stats.binomtest(n_neg, len(diffs), 0.5).pvalue)
            print(f"  {wname}: mean={mean_d:.3f}, {n_neg}/{len(diffs)} neg, t P={t_p:.4f}, W P={w_p:.4f}, sign P={sign_p:.4f}")
            window_results.append({
                "window": wname, "mean_diff": round(mean_d, 4),
                "n_negative": f"{n_neg}/{len(diffs)}",
                "ttest_p": round(float(t_p), 6),
                "wilcoxon_p": round(float(w_p), 6),
                "sign_p": round(float(sign_p), 6)
            })
        results[f"temporal_{label}"] = window_results

    # Formal pre-vs-post test
    print("\n--- Formal pre-vs-post test ---")
    for col, label in [('dry_natural_size_1234', 'swiss_dry'),
                       ('norway_aval_count', 'norway')]:
        pre_diffs = []
        post_diffs = []
        for ssw_date in ssw_dates:
            ev_pre, ctrl_pre = get_matched_control(ssw_date, ssw_dates, winter, col, -15, 15)
            ev_post, ctrl_post = get_matched_control(ssw_date, ssw_dates, winter, col, 0, 15)
            if not any(np.isnan([ev_pre, ctrl_pre, ev_post, ctrl_post])):
                pre_diffs.append(ev_pre - ctrl_pre)
                post_diffs.append(ev_post - ctrl_post)
        diff_of_diff = np.array(post_diffs) - np.array(pre_diffs)
        valid = diff_of_diff[~np.isnan(diff_of_diff)]
        perm_p = exact_signflip_pvalue(valid)
        boot_means = [np.mean(np.random.choice(valid, len(valid), replace=True)) for _ in range(10000)]
        ci = (float(np.percentile(boot_means, 2.5)), float(np.percentile(boot_means, 97.5)))
        print(f"  {label}: mean post-pre diff={np.mean(valid):.3f}, perm P={perm_p:.4f}, CI [{ci[0]:.3f}, {ci[1]:.3f}]")
        results[f"pre_vs_post_{label}"] = {
            "mean_diff_of_diff": round(float(np.mean(valid)), 4),
            "perm_p": round(perm_p, 4),
            "ci_lo": round(ci[0], 4),
            "ci_hi": round(ci[1], 4),
            "n_events": len(valid)
        }

    # Event-study distributed lag
    print("\n--- Event-study distributed lag ---")
    lags = range(-30, 61)
    for col, label in [('dry_natural_size_1234', 'swiss_dry'),
                       ('norway_aval_count', 'norway')]:
        # Seasonal expectation
        seasonal = winter.groupby('day_of_season')[col].mean()
        lag_data = []
        for lag in lags:
            anomalies = []
            for ssw_date in ssw_dates:
                target = ssw_date + pd.Timedelta(days=lag)
                if target in winter.index:
                    val = winter.loc[target, col]
                    dos = winter.loc[target, 'day_of_season']
                    if dos in seasonal.index and not np.isnan(val):
                        anomalies.append(val - seasonal[dos])
            if len(anomalies) >= 5:
                arr = np.array(anomalies)
                t_stat, t_p = stats.ttest_1samp(arr, 0)
                lag_data.append({
                    "lag": int(lag),
                    "mean_anomaly": round(float(np.mean(arr)), 4),
                    "se": round(float(np.std(arr, ddof=1) / np.sqrt(len(arr))), 4),
                    "p": round(float(t_p), 4),
                    "n": len(arr)
                })
        # Summary stats for pre/post windows
        pre_anoms = [d["mean_anomaly"] for d in lag_data if -15 <= d["lag"] < 0]
        post_anoms = [d["mean_anomaly"] for d in lag_data if 0 <= d["lag"] < 15]
        late_anoms = [d["mean_anomaly"] for d in lag_data if 15 <= d["lag"] < 30]
        print(f"  {label}: pre mean={np.mean(pre_anoms):.3f}, post mean={np.mean(post_anoms):.3f}, late mean={np.mean(late_anoms):.3f}")
        results[f"event_study_{label}"] = {
            "lag_data": lag_data,
            "pre_mean": round(float(np.mean(pre_anoms)), 4),
            "post_mean": round(float(np.mean(post_anoms)), 4),
            "late_mean": round(float(np.mean(late_anoms)), 4),
            "n_sig_pre": sum(1 for d in lag_data if -15 <= d["lag"] < 0 and d["p"] < 0.05),
            "n_sig_post": sum(1 for d in lag_data if 0 <= d["lag"] < 15 and d["p"] < 0.05)
        }

    save_json(results, "fresh_part3_temporal.json")
    print("\n[PART 3 COMPLETE]")

# ================================================================
# PART 4: METEOROLOGICAL CHAIN + NAO + U850
# ================================================================
elif part == 4:
    panel, winter, ssw, ssw_dates = load_data()
    results = {"part": 4, "description": "Meteorological chain, NAO, and U850 tests"}

    # 4-step pre-specified met chain
    met_vars = [
        ("ncep_t_10hpa", "T 10hPa (K)"),
        ("ncep_u_10hpa", "U 10hPa (m/s)"),
        ("ncep_u850_nh", "U 850hPa (m/s)"),
        ("ncep_z500_nh", "Z500 (m)"),
        ("ncep_slp_nh", "SLP (hPa)")
    ]

    print("--- 4-step meteorological chain ---")
    chain_results = []
    for col, label in met_vars:
        pre_vals = []
        post_vals = []
        for ssw_date in ssw_dates:
            pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
            post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
            pre_v = winter.loc[pre_mask, col].dropna()
            post_v = winter.loc[post_mask, col].dropna()
            if len(pre_v) > 0 and len(post_v) > 0:
                pre_vals.append(pre_v.mean())
                post_vals.append(post_v.mean())
        pre_arr = np.array(pre_vals)
        post_arr = np.array(post_vals)
        change = post_arr - pre_arr
        t_stat, t_p = stats.ttest_rel(post_arr, pre_arr)
        try:
            w_stat, w_p = stats.wilcoxon(change)
        except:
            w_p = np.nan
        mean_change = float(np.mean(change))
        print(f"  {label}: change={mean_change:.2f}, t P={t_p:.6f}, W P={w_p:.6f}")
        chain_results.append({
            "variable": col,
            "label": label,
            "mean_pre": round(float(np.mean(pre_arr)), 2),
            "mean_post": round(float(np.mean(post_arr)), 2),
            "mean_change": round(mean_change, 2),
            "ttest_p": round(float(t_p), 6),
            "wilcoxon_p": round(float(w_p) if not np.isnan(w_p) else 1.0, 6),
            "bonferroni_p": round(min(1.0, float(t_p) * 4), 6)  # 4 primary tests
        })
    results["met_chain"] = chain_results

    # Event-by-event: strat warming magnitude vs avalanche decrease
    print("\n--- Event-by-event correlations ---")
    strat_changes = []
    aval_diffs = []
    u850_changes = []
    for ssw_date in ssw_dates:
        pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
        post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
        
        t_pre = winter.loc[pre_mask, 'ncep_t_10hpa'].dropna()
        t_post = winter.loc[post_mask, 'ncep_t_10hpa'].dropna()
        u850_pre = winter.loc[pre_mask, 'ncep_u850_nh'].dropna()
        u850_post = winter.loc[post_mask, 'ncep_u850_nh'].dropna()
        
        ev_aval, ctrl_aval = get_matched_control(ssw_date, ssw_dates, winter, 'dry_natural_size_1234', 0, 15)
        
        if len(t_pre) > 0 and len(t_post) > 0 and not np.isnan(ev_aval) and not np.isnan(ctrl_aval):
            strat_changes.append(t_post.mean() - t_pre.mean())
            aval_diffs.append(ev_aval - ctrl_aval)
            if len(u850_pre) > 0 and len(u850_post) > 0:
                u850_changes.append(u850_post.mean() - u850_pre.mean())

    # T10 vs avalanche
    r_t, p_t = stats.pearsonr(strat_changes, aval_diffs)
    rho_t, prho_t = stats.spearmanr(strat_changes, aval_diffs)
    print(f"  ΔT10hPa vs Δdry_aval: Pearson r={r_t:.3f} P={p_t:.4f}, Spearman ρ={rho_t:.3f} P={prho_t:.4f}")
    results["event_corr_t10_aval"] = {
        "pearson_r": round(r_t, 4), "pearson_p": round(float(p_t), 4),
        "spearman_rho": round(rho_t, 4), "spearman_p": round(float(prho_t), 4),
        "n": len(strat_changes)
    }

    # U850 vs avalanche (event-level)
    if len(u850_changes) == len(aval_diffs):
        r_u, p_u = stats.pearsonr(u850_changes, aval_diffs)
        print(f"  ΔU850 vs Δdry_aval: r={r_u:.3f} P={p_u:.4f}")
        results["event_corr_u850_aval"] = {
            "pearson_r": round(r_u, 4), "pearson_p": round(float(p_u), 4),
            "n": len(u850_changes)
        }

    # U850 vs avalanche (daily)
    valid_mask = winter['ncep_u850_nh'].notna() & winter['dry_natural_size_1234'].notna()
    r_daily, p_daily = stats.pearsonr(
        winter.loc[valid_mask, 'ncep_u850_nh'],
        winter.loc[valid_mask, 'dry_natural_size_1234']
    )
    print(f"  Daily U850 vs dry_aval: r={r_daily:.4f} P={p_daily:.4f} (n={valid_mask.sum()})")
    results["daily_u850_aval"] = {
        "pearson_r": round(r_daily, 4), "pearson_p": round(float(p_daily), 4),
        "n": int(valid_mask.sum())
    }

    # NAO analysis
    print("\n--- NAO-avalanche relationship ---")
    nao_valid = winter['nao_daily'].notna() & winter['dry_natural_size_1234'].notna()
    r_nao, p_nao = stats.pearsonr(
        winter.loc[nao_valid, 'nao_daily'],
        winter.loc[nao_valid, 'dry_natural_size_1234']
    )
    # NAO terciles
    nao_vals = winter.loc[nao_valid, 'nao_daily']
    aval_vals = winter.loc[nao_valid, 'dry_natural_size_1234']
    q33, q67 = nao_vals.quantile([0.333, 0.667])
    nao_neg = aval_vals[nao_vals <= q33].mean()
    nao_mid = aval_vals[(nao_vals > q33) & (nao_vals <= q67)].mean()
    nao_pos = aval_vals[nao_vals > q67].mean()
    print(f"  NAO-dry_aval: r={r_nao:.4f}, P={p_nao:.4f} (n={nao_valid.sum()})")
    print(f"  NAO terciles: negative={nao_neg:.3f}, middle={nao_mid:.3f}, positive={nao_pos:.3f}/day")
    results["nao_avalanche"] = {
        "pearson_r": round(r_nao, 4), "pearson_p": round(float(p_nao), 4),
        "n": int(nao_valid.sum()),
        "nao_negative_rate": round(float(nao_neg), 4),
        "nao_middle_rate": round(float(nao_mid), 4),
        "nao_positive_rate": round(float(nao_pos), 4)
    }

    # SSW -> NAO (mean shift test)
    nao_ssw = []
    nao_ctrl = []
    for ssw_date in ssw_dates:
        post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
        nao_post = winter.loc[post_mask, 'nao_daily'].dropna()
        if len(nao_post) > 0:
            nao_ssw.append(nao_post.mean())
    overall_nao_mean = winter['nao_daily'].dropna().mean()
    t_ssw_nao, p_ssw_nao = stats.ttest_1samp(nao_ssw, overall_nao_mean)
    print(f"  SSW→NAO: post-SSW NAO mean={np.mean(nao_ssw):.3f}, climatology={overall_nao_mean:.3f}, P={p_ssw_nao:.4f}")
    results["ssw_nao"] = {
        "post_ssw_nao_mean": round(float(np.mean(nao_ssw)), 4),
        "climatology_mean": round(float(overall_nao_mean), 4),
        "ttest_p": round(float(p_ssw_nao), 4)
    }

    save_json(results, "fresh_part4_met_chain.json")
    print("\n[PART 4 COMPLETE]")

# ================================================================
# PART 5: NORWEGIAN COLD-REGIME + SNOTEL + PREDICTION MODEL
# ================================================================
elif part == 5:
    panel, winter, ssw, ssw_dates = load_data()
    results = {"part": 5, "description": "Norwegian cold-regime, SNOTEL, prediction model"}

    # Norwegian cold-regime stratification
    print("--- Norwegian cold-regime stratification ---")
    t100 = winter['ncep_t_100hpa']
    t100_median = t100.median()
    print(f"  100hPa T median: {t100_median:.1f} K")

    ssw_mask = winter['ssw_within_15d'] == 1 if 'ssw_within_15d' in winter.columns else pd.Series(False, index=winter.index)
    if ssw_mask.sum() == 0:
        for ssw_date in ssw_dates:
            mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
            ssw_mask = ssw_mask | mask

    for regime, rmask in [("cold", t100 < t100_median), ("warm", t100 >= t100_median)]:
        ssw_in_regime = winter.loc[ssw_mask & rmask, 'norway_aval_count'].dropna()
        non_ssw_in_regime = winter.loc[~ssw_mask & rmask, 'norway_aval_count'].dropna()
        if len(ssw_in_regime) > 5 and len(non_ssw_in_regime) > 5:
            diff = ssw_in_regime.mean() - non_ssw_in_regime.mean()
            t_stat, t_p = stats.ttest_ind(ssw_in_regime, non_ssw_in_regime)
            u_stat, u_p = stats.mannwhitneyu(ssw_in_regime, non_ssw_in_regime, alternative='two-sided')
            print(f"  {regime}: SSW mean={ssw_in_regime.mean():.2f}, non-SSW mean={non_ssw_in_regime.mean():.2f}, Δ={diff:.2f}, t P={t_p:.6f}, U P={u_p:.6f}")
            results[f"norway_{regime}"] = {
                "ssw_mean": round(ssw_in_regime.mean(), 4),
                "non_ssw_mean": round(non_ssw_in_regime.mean(), 4),
                "difference": round(diff, 4),
                "ttest_p": round(float(t_p), 6),
                "mannwhitney_p": round(float(u_p), 6),
                "n_ssw_days": len(ssw_in_regime),
                "n_non_ssw_days": len(non_ssw_in_regime)
            }

    # SNOTEL western US
    print("\n--- SNOTEL western US ---")
    if 'snotel_swe_mean' in winter.columns:
        winter['dswe'] = winter['snotel_swe_mean'].diff()
        for col, label in [('dswe', 'dSWE/dt'), ('snotel_temp_mean', 'Temperature')]:
            ssw_vals = winter.loc[ssw_mask, col].dropna()
            non_ssw_vals = winter.loc[~ssw_mask, col].dropna()
            if len(ssw_vals) > 5 and len(non_ssw_vals) > 5:
                diff = ssw_vals.mean() - non_ssw_vals.mean()
                t_stat, t_p = stats.ttest_ind(ssw_vals, non_ssw_vals)
                print(f"  {label}: SSW mean={ssw_vals.mean():.3f}, non-SSW mean={non_ssw_vals.mean():.3f}, Δ={diff:.3f}, P={t_p:.4f}")
                results[f"snotel_{label}"] = {
                    "ssw_mean": round(float(ssw_vals.mean()), 4),
                    "non_ssw_mean": round(float(non_ssw_vals.mean()), 4),
                    "difference": round(diff, 4),
                    "ttest_p": round(float(t_p), 6)
                }
    else:
        print("  SNOTEL columns not available")

    # Prediction model comparison (leave-one-winter-out) using statsmodels Poisson
    print("\n--- Prediction model comparison ---")
    import statsmodels.api as sm

    y_col = 'dry_natural_size_1234'
    winter_pred = winter[[y_col, 'day_of_season', 'day_of_season_sq', 'winter_id']].copy()
    for c in ['ssw_within_15d', 'ncep_t_10hpa', 'post_event_1_3d']:
        if c in winter.columns:
            winter_pred[c] = winter[c]
    winter_pred = winter_pred.fillna(0)
    winter_pred['y'] = winter_pred[y_col].clip(lower=0).astype(int)

    base_cols = ['day_of_season', 'day_of_season_sq']
    enh_cols = base_cols.copy()
    for c in ['ssw_within_15d', 'ncep_t_10hpa', 'post_event_1_3d']:
        if c in winter_pred.columns:
            enh_cols.append(c)

    winters_list = sorted(winter_pred['winter_id'].unique())
    rmse_base, rmse_enh, ll_base, ll_enh = [], [], [], []

    for w in winters_list:
        test_mask = winter_pred['winter_id'] == w
        train = winter_pred[~test_mask]
        test = winter_pred[test_mask]
        if len(test) < 10 or len(train) < 50:
            continue
        try:
            X_b_tr = sm.add_constant(train[base_cols].astype(float))
            X_b_te = sm.add_constant(test[base_cols].astype(float))
            m_b = sm.GLM(train['y'], X_b_tr, family=sm.families.Poisson()).fit(disp=0)
            pred_b = m_b.predict(X_b_te)

            X_e_tr = sm.add_constant(train[enh_cols].astype(float))
            X_e_te = sm.add_constant(test[enh_cols].astype(float))
            m_e = sm.GLM(train['y'], X_e_tr, family=sm.families.Poisson()).fit(disp=0)
            pred_e = m_e.predict(X_e_te)

            y_test = test['y'].values.astype(float)
            rmse_base.append(np.sqrt(np.mean((y_test - pred_b)**2)))
            rmse_enh.append(np.sqrt(np.mean((y_test - pred_e)**2)))
            # Poisson log-likelihood
            ll_b = np.sum(y_test * np.log(pred_b + 1e-10) - pred_b)
            ll_e = np.sum(y_test * np.log(pred_e + 1e-10) - pred_e)
            ll_base.append(ll_b)
            ll_enh.append(ll_e)
        except Exception as e:
            print(f"  Winter {w}: {e}")
            continue

    rmse_base_arr = np.array(rmse_base)
    rmse_enh_arr = np.array(rmse_enh)
    ll_base_arr = np.array(ll_base)
    ll_enh_arr = np.array(ll_enh)
    t_rmse, p_rmse = stats.ttest_rel(rmse_base_arr, rmse_enh_arr)
    t_ll, p_ll = stats.ttest_rel(ll_base_arr, ll_enh_arr)
    n_improved = int(np.sum(rmse_enh_arr < rmse_base_arr))
    n_ll_improved = int(np.sum(ll_enh_arr > ll_base_arr))
    print(f"  RMSE: base={np.mean(rmse_base_arr):.3f}, enhanced={np.mean(rmse_enh_arr):.3f}, P={p_rmse:.4f}")
    print(f"  LogLik: base={np.mean(ll_base_arr):.1f}, enhanced={np.mean(ll_enh_arr):.1f}, P={p_ll:.4f}")
    print(f"  Folds RMSE improved: {n_improved}/{len(rmse_base_arr)}, LL improved: {n_ll_improved}/{len(ll_base_arr)}")
    results["prediction"] = {
        "rmse_base": round(float(np.mean(rmse_base_arr)), 4),
        "rmse_enhanced": round(float(np.mean(rmse_enh_arr)), 4),
        "rmse_paired_t_p": round(float(p_rmse), 4),
        "ll_base": round(float(np.mean(ll_base_arr)), 1),
        "ll_enhanced": round(float(np.mean(ll_enh_arr)), 1),
        "ll_paired_t_p": round(float(p_ll), 4),
        "n_folds": len(rmse_base_arr),
        "n_improved_rmse": n_improved,
        "n_improved_ll": n_ll_improved
    }

    # In-sample full model coefficients
    X_full = sm.add_constant(winter_pred[enh_cols].astype(float))
    m_full = sm.GLM(winter_pred['y'], X_full, family=sm.families.Poisson()).fit(disp=0)
    coefs = {}
    for i, c in enumerate(['const'] + enh_cols):
        coefs[c] = {"beta": round(float(m_full.params.iloc[i]), 4),
                     "p": round(float(m_full.pvalues.iloc[i]), 6)}
    print(f"  In-sample coefs:")
    for k, v in coefs.items():
        print(f"    {k}: beta={v['beta']}, p={v['p']}")
    results["insample_coefs"] = coefs

    save_json(results, "fresh_part5_secondary.json")
    print("\n[PART 5 COMPLETE]")

# ================================================================
# PART 6: ROBUSTNESS BATTERY
# ================================================================
elif part == 6:
    panel, winter, ssw, ssw_dates = load_data()
    results = {"part": 6, "description": "Robustness: specification curve and LOOCV"}

    def mantel_haenszel_rr(exposed_events, exposed_days, unexposed_events, unexposed_days):
        num = 0; den = 0; var_num = 0
        for e1, d1, e0, d0 in zip(exposed_events, exposed_days, unexposed_events, unexposed_days):
            T = d1 + d0
            if T == 0: continue
            num += e1 * d0 / T
            den += e0 * d1 / T
            var_num += (e1 * d0**2 + e0 * d1**2) / T**2
        if den == 0:
            return np.nan, np.nan
        rr = num / den
        se = np.sqrt(var_num / (num * den))
        p = 2 * stats.norm.sf(np.abs(np.log(rr) / se))
        return float(rr), float(p)

    # Specification curve: 60 variants
    print("--- Specification curve (60 variants) ---")
    outcomes = [('dry_natural_size_1234', 'dry'), ('wet_natural_size_1234', 'wet'), ('natural_size_1234', 'all_nat')]
    exposures = [
        ('post_event_1_3d', 'geomag_1_3d'),
        ('post_event_5_21d', 'geomag_5_21d'),
        ('ssw_within_15d', 'ssw_0_15d'),
        ('post_event_15_30d', 'geomag_15_30d'),
    ]
    geographies = [('swiss', ['dry_natural_size_1234', 'wet_natural_size_1234', 'natural_size_1234']),
                   ('norway', ['norway_aval_count'])]

    spec_results = []
    for out_col, out_label in outcomes:
        for exp_col, exp_label in exposures:
            if exp_col not in winter.columns:
                continue
            winter['period_sc'] = winter['day_of_season'] // 15
            winter['stratum_sc'] = winter['winter_id'].astype(str) + '_' + winter['period_sc'].astype(str)
            strata = winter.groupby('stratum_sc')
            exp_e, exp_d, unexp_e, unexp_d = [], [], [], []
            for _, g in strata:
                exposed = g[exp_col] == 1
                n_exp = exposed.sum()
                n_unexp = (~exposed).sum()
                if n_exp == 0 or n_unexp == 0:
                    continue
                exp_e.append(g.loc[exposed, out_col].sum())
                exp_d.append(n_exp)
                unexp_e.append(g.loc[~exposed, out_col].sum())
                unexp_d.append(n_unexp)
            rr, p = mantel_haenszel_rr(exp_e, exp_d, unexp_e, unexp_d)
            spec_results.append({
                "outcome": out_label, "exposure": exp_label, "geography": "swiss",
                "rr": round(rr, 4) if not np.isnan(rr) else None,
                "p": round(p, 6) if not np.isnan(p) else None,
                "decrease": rr < 1 if not np.isnan(rr) else False
            })

    # Norwegian specs
    if 'norway_aval_count' in winter.columns:
        for exp_col, exp_label in exposures:
            if exp_col not in winter.columns:
                continue
            winter['period_sc'] = winter['day_of_season'] // 15
            winter['stratum_sc'] = winter['winter_id'].astype(str) + '_' + winter['period_sc'].astype(str)
            strata = winter.groupby('stratum_sc')
            exp_e, exp_d, unexp_e, unexp_d = [], [], [], []
            for _, g in strata:
                exposed = g[exp_col] == 1
                n_exp = exposed.sum()
                n_unexp = (~exposed).sum()
                if n_exp == 0 or n_unexp == 0:
                    continue
                exp_e.append(g.loc[exposed, 'norway_aval_count'].sum())
                exp_d.append(n_exp)
                unexp_e.append(g.loc[~exposed, 'norway_aval_count'].sum())
                unexp_d.append(n_unexp)
            rr, p = mantel_haenszel_rr(exp_e, exp_d, unexp_e, unexp_d)
            spec_results.append({
                "outcome": "norway", "exposure": exp_label, "geography": "norway",
                "rr": round(rr, 4) if not np.isnan(rr) else None,
                "p": round(p, 6) if not np.isnan(p) else None,
                "decrease": rr < 1 if not np.isnan(rr) else False
            })

    valid_specs = [s for s in spec_results if s["rr"] is not None]
    n_decrease = sum(1 for s in valid_specs if s["decrease"])
    pct_decrease = n_decrease / len(valid_specs) * 100 if valid_specs else 0
    median_rr = float(np.median([s["rr"] for s in valid_specs])) if valid_specs else np.nan
    print(f"  Total specs: {len(valid_specs)}")
    print(f"  Decrease: {n_decrease}/{len(valid_specs)} ({pct_decrease:.0f}%)")
    print(f"  Median RR: {median_rr:.3f}")
    results["spec_curve"] = {
        "n_specs": len(valid_specs),
        "n_decrease": n_decrease,
        "pct_decrease": round(pct_decrease, 1),
        "median_rr": round(median_rr, 4),
        "specs": valid_specs
    }

    # LOOCV (leave-one-winter-out matched SSW)
    print("\n--- LOOCV (leave-one-winter-out) ---")
    all_winters = sorted(winter['winter_id'].unique())
    loo_rrs = []
    for held_out in all_winters:
        subset = winter[winter['winter_id'] != held_out]
        ssw_in_sub = ssw_dates[ssw_dates.isin(subset.index) |
                               ((ssw_dates >= subset.index.min()) & (ssw_dates <= subset.index.max()))]
        diffs = []
        for ssw_date in ssw_in_sub:
            ev, ctrl = get_matched_control(ssw_date, ssw_in_sub, subset, 'dry_natural_size_1234', 0, 15)
            if not np.isnan(ev) and not np.isnan(ctrl) and ctrl > 0:
                diffs.append(ev / ctrl)
        if diffs:
            rr = float(np.mean(diffs))
            loo_rrs.append({"winter": held_out, "mean_rr": round(rr, 4)})

    loo_arr = np.array([l["mean_rr"] for l in loo_rrs])
    n_below_1 = int(np.sum(loo_arr < 1))
    print(f"  All folds RR<1: {n_below_1}/{len(loo_arr)}")
    print(f"  Mean RR: {np.mean(loo_arr):.3f}, SD: {np.std(loo_arr):.3f}")
    results["loocv"] = {
        "folds": loo_rrs,
        "n_below_1": n_below_1,
        "n_folds": len(loo_rrs),
        "mean_rr": round(float(np.mean(loo_arr)), 4),
        "sd_rr": round(float(np.std(loo_arr)), 4)
    }

    # Winter fixed-effects
    print("\n--- Winter fixed-effects ---")
    winter['period_wfe'] = winter['day_of_season'] // 15
    winter['stratum_wfe'] = winter['winter_id'].astype(str) + '_' + winter['period_wfe'].astype(str)
    strata = winter.groupby('stratum_wfe')
    for col, label in [('dry_natural_size_1234', 'dry_slab'), ('natural_size_1234', 'all_natural')]:
        for exp_col, elabel in [('post_event_1_3d', 'geomag'), ('ssw_within_15d', 'ssw')]:
            if exp_col not in winter.columns:
                continue
            exp_e, exp_d, unexp_e, unexp_d = [], [], [], []
            for _, g in strata:
                exposed = g[exp_col] == 1
                n_exp = exposed.sum()
                n_unexp = (~exposed).sum()
                if n_exp == 0 or n_unexp == 0:
                    continue
                exp_e.append(g.loc[exposed, col].sum())
                exp_d.append(n_exp)
                unexp_e.append(g.loc[~exposed, col].sum())
                unexp_d.append(n_unexp)
            rr, p = mantel_haenszel_rr(exp_e, exp_d, unexp_e, unexp_d)
            print(f"  WFE {label} x {elabel}: RR={rr:.3f}, P={p:.6f}")
            results[f"wfe_{label}_{elabel}"] = {"rr": round(rr, 4), "p": round(p, 6)}

    save_json(results, "fresh_part6_robustness.json")
    print("\n[PART 6 COMPLETE]")

# ================================================================
# PART 7: PUBLICATION-QUALITY FIGURES
# ================================================================
elif part == 7:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    panel, winter, ssw, ssw_dates = load_data()
    
    # Load fresh results
    with open(RESULTS / 'fresh_part1_ssw_primary.json') as f:
        r1 = json.load(f)
    with open(RESULTS / 'fresh_part3_temporal.json') as f:
        r3 = json.load(f)

    plt.rcParams.update({
        'font.size': 9, 'axes.labelsize': 10, 'axes.titlesize': 11,
        'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'font.family': 'sans-serif', 'figure.dpi': 300
    })

    # FIGURE 1: SSW matched comparison bar plots
    print("--- Figure 1: SSW matched comparison ---")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0))
    
    for ax_idx, (key, title, ylabel) in enumerate([
        ('swiss_dry_slab', 'Swiss dry slab', 'Δ avalanches/day'),
        ('norway_total', 'Norway total', 'Δ avalanches/day'),
    ]):
        ax = axes[ax_idx]
        events = r1[key]['events']
        dates = [e['date'] for e in events]
        diffs = [e['difference'] for e in events]
        colors = ['#2ecc71' if d < 0 else '#e74c3c' for d in diffs]
        ax.bar(range(len(diffs)), diffs, color=colors, edgecolor='k', linewidth=0.3)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_title(title, fontweight='bold')
        ax.set_ylabel(ylabel)
        ax.set_xlabel('SSW event')
        ax.set_xticks(range(len(diffs)))
        ax.set_xticklabels([d[:4] for d in dates], rotation=45, fontsize=6)
        mean_d = r1[key]['mean_diff']
        ci_lo = r1[key]['bootstrap_ci_lo']
        ci_hi = r1[key]['bootstrap_ci_hi']
        ax.text(0.02, 0.02, f'Δ={mean_d:.2f}\n[{ci_lo:.2f}, {ci_hi:.2f}]',
                transform=ax.transAxes, fontsize=7, va='bottom',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Panel c: all-natural (compute on the fly)
    ax = axes[2]
    diffs_all = []
    dates_all = []
    for ssw_date in ssw_dates:
        ev, ctrl = get_matched_control(ssw_date, ssw_dates, winter, 'natural_size_1234', 0, 15)
        if not np.isnan(ev) and not np.isnan(ctrl):
            diffs_all.append(ev - ctrl)
            dates_all.append(str(ssw_date.date()))
    colors = ['#2ecc71' if d < 0 else '#e74c3c' for d in diffs_all]
    ax.bar(range(len(diffs_all)), diffs_all, color=colors, edgecolor='k', linewidth=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.set_title('Swiss all-natural', fontweight='bold')
    ax.set_ylabel('Δ avalanches/day')
    ax.set_xlabel('SSW event')
    ax.set_xticks(range(len(diffs_all)))
    ax.set_xticklabels([d[:4] for d in dates_all], rotation=45, fontsize=6)
    
    for i, ax in enumerate(axes):
        ax.text(-0.05, 1.05, chr(97+i), transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')
    
    fig.tight_layout()
    fig.savefig(FIGURES / 'fig1_ssw_matched_fresh.pdf', bbox_inches='tight')
    fig.savefig(FIGURES / 'fig1_ssw_matched_fresh.png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved fig1_ssw_matched_fresh.pdf/png")

    # FIGURE 2: Type specificity (forest plot)
    print("--- Figure 2: Type specificity ---")
    with open(RESULTS / 'fresh_part2_type_specificity.json') as f:
        r2 = json.load(f)
    
    fig, ax = plt.subplots(figsize=(4.5, 3.0))
    labels = ['Dry slab\n(geomag 1-3d)', 'All natural\n(geomag 1-3d)', 'Wet natural\n(geomag 1-3d)']
    keys = ['mh_dry_natural_geomag', 'mh_all_natural_geomag', 'mh_wet_natural_geomag']
    y_pos = [2, 1, 0]
    for i, (lab, key) in enumerate(zip(labels, keys)):
        if key in r2:
            rr = r2[key]['rate_ratio']
            ci_lo = r2[key]['ci_lo']
            ci_hi = r2[key]['ci_hi']
            p = r2[key]['p_value']
            color = '#e74c3c' if p < 0.05 else '#95a5a6'
            ax.errorbar(rr, y_pos[i], xerr=[[rr-ci_lo], [ci_hi-rr]], fmt='o',
                       color=color, markersize=8, capsize=4, linewidth=1.5)
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            ax.text(ci_hi + 0.02, y_pos[i], f'RR={rr:.2f} {sig}', va='center', fontsize=8)
    ax.axvline(1, color='k', linestyle='--', linewidth=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Mantel-Haenszel Rate Ratio')
    ax.set_title('Type Specificity: Geomagnetic Storm Exposure', fontweight='bold')
    ax.set_xlim(0.3, 1.5)
    fig.tight_layout()
    fig.savefig(FIGURES / 'fig2_type_specificity_fresh.pdf', bbox_inches='tight')
    fig.savefig(FIGURES / 'fig2_type_specificity_fresh.png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved fig2_type_specificity_fresh.pdf/png")

    # FIGURE 3: Stratum-width sensitivity
    print("--- Figure 3: Stratum-width sensitivity ---")
    fig, ax = plt.subplots(figsize=(5, 3.5))
    sens = r2['isolated_sensitivity']
    widths = sorted(set(s['stratum_width'] for s in sens))
    post_rrs = [next(s['rr'] for s in sens if s['stratum_width'] == w and s['direction'] == 'post') for w in widths]
    pre_rrs = [next(s['rr'] for s in sens if s['stratum_width'] == w and s['direction'] == 'pre') for w in widths]
    x = np.arange(len(widths))
    ax.bar(x - 0.15, post_rrs, 0.3, label='Post-event', color='#3498db', edgecolor='k', linewidth=0.3)
    ax.bar(x + 0.15, pre_rrs, 0.3, label='Pre-event', color='#e67e22', edgecolor='k', linewidth=0.3)
    ax.axhline(1, color='k', linestyle='--', linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{w}d' for w in widths])
    ax.set_xlabel('Stratum width')
    ax.set_ylabel('Mantel-Haenszel Rate Ratio')
    ax.set_title('Isolated Event Stratum-Width Sensitivity', fontweight='bold')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / 'fig3_stratum_sensitivity_fresh.pdf', bbox_inches='tight')
    fig.savefig(FIGURES / 'fig3_stratum_sensitivity_fresh.png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved fig3_stratum_sensitivity_fresh.pdf/png")

    # FIGURE 4: Event-study distributed lag
    print("--- Figure 4: Event-study ---")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for ax_idx, (key, title) in enumerate([
        ('event_study_swiss_dry', 'Swiss dry slab'),
        ('event_study_norway', 'Norway total')
    ]):
        ax = axes[ax_idx]
        lag_data = r3[key]['lag_data']
        lags_arr = [d['lag'] for d in lag_data]
        means = [d['mean_anomaly'] for d in lag_data]
        ses = [d['se'] for d in lag_data]
        upper = [m + 1.96*s for m, s in zip(means, ses)]
        lower = [m - 1.96*s for m, s in zip(means, ses)]
        
        ax.fill_between(lags_arr, lower, upper, alpha=0.2, color='steelblue')
        ax.plot(lags_arr, means, color='navy', linewidth=0.5, alpha=0.3)
        # Rolling mean
        means_arr = np.array(means)
        kernel = np.ones(7) / 7
        smooth = np.convolve(means_arr, kernel, mode='same')
        ax.plot(lags_arr, smooth, color='navy', linewidth=1.5, label='7-day mean')
        ax.axhline(0, color='k', linewidth=0.5)
        ax.axvline(0, color='red', linewidth=1, linestyle='--', label='SSW onset')
        ax.set_xlabel('Days relative to SSW onset')
        ax.set_ylabel('Anomaly (events/day)')
        ax.set_title(title, fontweight='bold')
        ax.legend(fontsize=7)
        ax.text(-0.05, 1.05, chr(97+ax_idx), transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')
    
    fig.tight_layout()
    fig.savefig(FIGURES / 'fig4_event_study_fresh.pdf', bbox_inches='tight')
    fig.savefig(FIGURES / 'fig4_event_study_fresh.png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  Saved fig4_event_study_fresh.pdf/png")

    print("\n[PART 7 COMPLETE - ALL FIGURES SAVED]")
    gc.collect()

print(f"\n{'='*70}")
print(f"  ALL DONE")
print(f"{'='*70}")
