"""
19_tier2_upgrade.py — Upgrade Tier 2 findings with rigorous methodology
=========================================================================
Based on critic feedback:
PART 1: Permutation-based specification curve test (not binomial)
PART 2: SSW all-natural — Wilcoxon + sign + influence + bootstrap CI
PART 3: LOOCV — held-out deviance improvement (paired test)
PART 4: Isolated events + stratum-width sensitivity for pre/post
PART 5: Stratum FE Poisson + placebo retest
PART 6: Expanded specification curve (design-based specs only)
"""
import sys, gc, json, logging, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, LOG

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

RESULTS_FILE = RESULTS / "tier2_upgrade.json"


def load_panel():
    panel = pd.read_parquet(PROCESSED / "analysis_panel_v2.parquet")
    return panel[panel["is_winter"] == 1].copy()


def load_ssw():
    ssw_cat = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    if hasattr(ssw_cat.index, 'tz') and ssw_cat.index.tz:
        ssw_cat.index = ssw_cat.index.tz_localize(None)
    return ssw_cat


def save(results):
    with open(RESULTS_FILE, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=str)


def mh_rate_ratio(valid, outcome_col, exposure_col, stratum_col):
    """Compute Mantel-Haenszel stratified rate ratio."""
    exposed_sum = unexposed_sum = 0.0
    n_strata = 0
    for _, group in valid.groupby(stratum_col):
        exp = group[group[exposure_col] == 1]
        unexp = group[group[exposure_col] == 0]
        if len(exp) == 0 or len(unexp) == 0:
            continue
        n_total = len(exp) + len(unexp)
        exposed_sum += exp[outcome_col].sum() * len(unexp) / n_total
        unexposed_sum += unexp[outcome_col].sum() * len(exp) / n_total
        n_strata += 1
    if exposed_sum <= 0 or unexposed_sum <= 0:
        return None
    rr = exposed_sum / unexposed_sum
    ln_rr = np.log(rr)
    se = np.sqrt(1.0 / exposed_sum + 1.0 / unexposed_sum)
    z = ln_rr / se
    p = 2 * stats.norm.sf(abs(z))
    return {"rr": float(rr), "z": float(z), "p": float(p),
            "ci_lo": float(np.exp(ln_rr - 1.96*se)),
            "ci_hi": float(np.exp(ln_rr + 1.96*se)),
            "n_strata": n_strata}


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: Permutation-Based Specification Curve Test
# ═══════════════════════════════════════════════════════════════════════════════
def part1_perm_spec_curve(results):
    """
    Shuffle event labels within each winter (preserving temporal structure),
    rerun all specs, compare observed % RR<1 / median RR to null distribution.
    Uses MH case-crossover (design-based) not NB GLM.
    """
    LOG.info("=" * 70)
    LOG.info("PART 1: Permutation Spec Curve (500 iterations, MH-based)")
    LOG.info("=" * 70)

    w = load_panel()
    rng = np.random.RandomState(42)
    N_PERM = 500

    outcomes = [
        ("dry_natural", "dry_natural_size_1234"),
        ("all_natural", "aai_all_natural"),
    ]
    exposures = [("1_3d", "post_event_1_3d"), ("5_21d", "post_event_5_21d")]
    strata_widths = [7, 10, 15, 21]

    def run_all_specs(w_df):
        """Run all MH specs and return list of RRs."""
        rrs = []
        for out_name, out_col in outcomes:
            if out_col not in w_df.columns:
                continue
            y = w_df[out_col].dropna()
            valid = w_df.loc[y.index].copy()
            for exp_name, exp_col in exposures:
                if exp_col not in valid.columns:
                    continue
                for sw in strata_widths:
                    valid["_strat"] = valid["winter_id"].astype(str) + "_" + \
                                      (valid["day_of_season"] // sw).astype(int).astype(str)
                    res = mh_rate_ratio(valid, out_col, exp_col, "_strat")
                    if res:
                        rrs.append(res["rr"])
        return rrs

    # Observed
    obs_rrs = run_all_specs(w)
    obs_n = len(obs_rrs)
    obs_pct_decrease = 100 * np.mean([r < 1 for r in obs_rrs])
    obs_median_rr = np.median(obs_rrs)
    LOG.info("  Observed: %d specs, %.1f%% decrease, median RR=%.3f",
             obs_n, obs_pct_decrease, obs_median_rr)

    # Permutation null: circular shift events within each winter
    null_pct = []
    null_median = []
    for perm_i in range(N_PERM):
        w_perm = w.copy()
        for exp_name, exp_col in exposures:
            for wid, group in w_perm.groupby("winter_id"):
                if wid is None:
                    continue
                vals = group[exp_col].values
                shift = rng.randint(0, len(group))
                w_perm.loc[group.index, exp_col] = np.roll(vals, shift)

        perm_rrs = run_all_specs(w_perm)
        if perm_rrs:
            null_pct.append(np.mean([r < 1 for r in perm_rrs]) * 100)
            null_median.append(np.median(perm_rrs))

        if (perm_i + 1) % 100 == 0:
            LOG.info("    Permutation %d/%d done", perm_i + 1, N_PERM)

    null_pct = np.array(null_pct)
    null_median = np.array(null_median)

    p_pct = np.mean(null_pct >= obs_pct_decrease)
    p_median = np.mean(null_median <= obs_median_rr)

    section = {
        "n_specs": obs_n,
        "observed_pct_decrease": float(obs_pct_decrease),
        "observed_median_rr": float(obs_median_rr),
        "null_mean_pct": float(null_pct.mean()),
        "null_std_pct": float(null_pct.std()),
        "p_value_pct_decrease": float(p_pct),
        "p_value_median_rr": float(p_median),
        "n_permutations": N_PERM,
        "method": "Circular shift within winter, MH case-crossover specs",
        "interpretation": (
            "p_pct tests whether observed fraction of specs with RR<1 exceeds chance. "
            "p_median tests whether observed median RR is lower than chance."
        ),
    }
    LOG.info("  Perm p(pct_decrease): %.4f, Perm p(median_rr): %.4f", p_pct, p_median)

    results["part1_perm_spec_curve"] = section
    save(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: SSW All-Natural — Complete Statistical Battery
# ═══════════════════════════════════════════════════════════════════════════════
def part2_ssw_all_natural(results):
    """
    Report t-test, sign test, Wilcoxon signed-rank, exact permutation,
    bootstrap CI, and leave-one-out influence for SSW all-natural.
    """
    LOG.info("=" * 70)
    LOG.info("PART 2: SSW All-Natural — Full Statistical Battery")
    LOG.info("=" * 70)

    w = load_panel()
    ssw_cat = load_ssw()

    ssw_winter_ids = set()
    for sd in ssw_cat.index:
        if sd.month >= 11:
            ssw_winter_ids.add(f"{sd.year}/{sd.year+1}")
        elif sd.month <= 3:
            ssw_winter_ids.add(f"{sd.year-1}/{sd.year}")
    non_ssw_winters = [wid for wid in w["winter_id"].dropna().unique()
                       if wid not in ssw_winter_ids]

    section = {}

    for out_name, out_col in [("all_natural", "aai_all_natural"),
                               ("dry_natural", "dry_natural_size_1234"),
                               ("norway", "norway_aval_count")]:
        if out_col not in w.columns:
            continue

        diffs = []
        event_detail = []
        for sd in ssw_cat.index:
            if sd.month >= 11:
                sw = f"{sd.year}/{sd.year+1}"
                dos = (sd - pd.Timestamp(sd.year, 11, 1)).days
            elif sd.month <= 3:
                sw = f"{sd.year-1}/{sd.year}"
                dos = (sd - pd.Timestamp(sd.year-1, 11, 1)).days
            else:
                continue
            if sw not in w["winter_id"].values:
                continue

            post = w[(w.index >= sd) & (w.index < sd + pd.Timedelta(days=15))]
            ssw_mean = post[out_col].mean()
            if np.isnan(ssw_mean):
                continue

            ctrl_means = []
            for ctrl_wid in non_ssw_winters:
                ctrl = w[w["winter_id"] == ctrl_wid]
                if len(ctrl) == 0:
                    continue
                matched = ctrl[(ctrl["day_of_season"] >= dos - 3) &
                               (ctrl["day_of_season"] <= dos + 15 + 3)]
                if len(matched) >= 5:
                    cm = matched[out_col].mean()
                    if not np.isnan(cm):
                        ctrl_means.append(cm)

            if len(ctrl_means) >= 3:
                diff = ssw_mean - np.mean(ctrl_means)
                diffs.append(diff)
                event_detail.append({
                    "date": str(sd.date()), "ssw_mean": float(ssw_mean),
                    "ctrl_mean": float(np.mean(ctrl_means)), "diff": float(diff)
                })

        if len(diffs) < 5:
            continue

        d = np.array(diffs)

        # 1. Paired t-test
        t_stat, t_p = stats.ttest_1samp(d, 0)

        # 2. Sign test
        n_neg = np.sum(d < 0)
        n_total = len(d)
        # Exact binomial
        sign_p = 2 * min(stats.binom.cdf(n_neg, n_total, 0.5),
                          1 - stats.binom.cdf(n_neg - 1, n_total, 0.5))

        # 3. Wilcoxon signed-rank (robust to outliers)
        try:
            w_stat, w_p = stats.wilcoxon(d, alternative="two-sided")
        except Exception:
            w_stat, w_p = np.nan, np.nan

        # 4. Exact permutation test (2^15 = 32768 possible sign flips)
        n_perm = min(2**len(d), 50000)
        rng = np.random.RandomState(123)
        obs_mean = d.mean()
        if 2**len(d) <= 50000:
            # Exact enumeration
            count_extreme = 0
            for i in range(2**len(d)):
                signs = np.array([(i >> j) & 1 for j in range(len(d))]) * 2 - 1
                perm_mean = (d * signs).mean()
                if abs(perm_mean) >= abs(obs_mean):
                    count_extreme += 1
            exact_p = count_extreme / (2**len(d))
        else:
            perm_means = []
            for _ in range(50000):
                signs = rng.choice([-1, 1], size=len(d))
                perm_means.append((d * signs).mean())
            exact_p = np.mean(np.abs(perm_means) >= abs(obs_mean))

        # 5. Bootstrap CI (10000 resamples)
        boot_means = []
        for _ in range(10000):
            resample = rng.choice(d, size=len(d), replace=True)
            boot_means.append(resample.mean())
        boot_means = np.array(boot_means)
        ci_lo = np.percentile(boot_means, 2.5)
        ci_hi = np.percentile(boot_means, 97.5)

        # 6. Leave-one-out influence
        loo_results = []
        for i in range(len(d)):
            d_loo = np.delete(d, i)
            t_loo, p_loo = stats.ttest_1samp(d_loo, 0)
            n_neg_loo = np.sum(d_loo < 0)
            loo_results.append({
                "dropped": event_detail[i]["date"],
                "dropped_diff": float(d[i]),
                "remaining_mean": float(d_loo.mean()),
                "t_p": float(p_loo),
                "n_neg": int(n_neg_loo),
            })

        # Identify the influential outlier
        loo_df = pd.DataFrame(loo_results)
        most_influential = loo_df.loc[loo_df["t_p"].idxmin()]

        section[out_name] = {
            "n_events": n_total,
            "mean_diff": float(d.mean()),
            "median_diff": float(np.median(d)),
            "t_test": {"stat": float(t_stat), "p": float(t_p)},
            "sign_test": {"n_negative": int(n_neg), "n_total": n_total,
                          "p": float(sign_p)},
            "wilcoxon": {"stat": float(w_stat) if not np.isnan(w_stat) else None,
                         "p": float(w_p) if not np.isnan(w_p) else None},
            "exact_permutation": {"p": float(exact_p),
                                   "method": "exact sign-flip" if 2**len(d) <= 50000 else "50000 random sign-flips"},
            "bootstrap_ci": {"ci_2.5": float(ci_lo), "ci_97.5": float(ci_hi)},
            "leave_one_out": {
                "most_influential": {
                    "date": most_influential["dropped"],
                    "diff": float(most_influential["dropped_diff"]),
                    "p_without": float(most_influential["t_p"]),
                },
                "n_folds_p_lt_0.05_without": int((loo_df["t_p"] < 0.05).sum()),
                "all_loo_n_neg_gte_12": bool((loo_df["n_neg"] >= 12).all()),
            },
        }

        LOG.info("  %s (n=%d): mean=%.3f, median=%.3f", out_name, n_total, d.mean(), np.median(d))
        LOG.info("    t-test p=%.4f | sign p=%.4f | Wilcoxon p=%.4f | perm p=%.4f",
                 t_p, sign_p, w_p if not np.isnan(w_p) else -1, exact_p)
        LOG.info("    Bootstrap CI: [%.3f, %.3f]", ci_lo, ci_hi)
        LOG.info("    Most influential: %s (diff=%.2f, p_without=%.4f)",
                 most_influential["dropped"], most_influential["dropped_diff"],
                 most_influential["t_p"])

    results["part2_ssw_battery"] = section
    save(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: LOOCV — Held-Out Deviance Improvement
# ═══════════════════════════════════════════════════════════════════════════════
def part3_loocv_deviance(results):
    """
    For each winter: fit baseline (season only) and event model (season + exposure),
    compute held-out log-likelihood. Paired test across winters.
    """
    LOG.info("=" * 70)
    LOG.info("PART 3: LOOCV — Held-Out Deviance Improvement")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    for hyp_name, out_col, exp_col in [
        ("dry_geomag_1_3d", "dry_natural_size_1234", "post_event_1_3d"),
        ("dry_geomag_5_21d", "dry_natural_size_1234", "post_event_5_21d"),
        ("all_geomag_1_3d", "aai_all_natural", "post_event_1_3d"),
    ]:
        if out_col not in w.columns:
            continue

        winters = sorted(w["winter_id"].dropna().unique())
        base_logliks = []
        event_logliks = []
        fold_rrs = []
        improvements = []

        for wid in winters:
            train = w[w["winter_id"] != wid]
            test = w[w["winter_id"] == wid]
            y_train = train[out_col].dropna()
            y_test = test[out_col].dropna()
            if len(y_train) < 100 or len(y_test) < 10:
                continue

            base_cols = ["day_of_season", "day_of_season_sq"]
            event_cols = base_cols + [exp_col]

            try:
                # Baseline
                X_tr_b = sm.add_constant(train.loc[y_train.index, base_cols])
                X_te_b = sm.add_constant(test.loc[y_test.index, base_cols])
                mask_tr = y_train.notna() & X_tr_b.notna().all(axis=1)
                mask_te = y_test.notna() & X_te_b.notna().all(axis=1)
                m_base = sm.GLM(y_train[mask_tr], X_tr_b[mask_tr],
                                family=sm.families.Poisson()).fit(maxiter=50)
                pred_base = m_base.predict(X_te_b[mask_te])
                y_te = y_test[mask_te].values
                # Poisson log-likelihood
                ll_base = np.sum(y_te * np.log(np.maximum(pred_base, 1e-10)) - pred_base)

                # Event model
                X_tr_e = sm.add_constant(train.loc[y_train.index, event_cols])
                X_te_e = sm.add_constant(test.loc[y_test.index, event_cols])
                mask_tr_e = y_train.notna() & X_tr_e.notna().all(axis=1)
                mask_te_e = y_test.notna() & X_te_e.notna().all(axis=1)
                m_event = sm.GLM(y_train[mask_tr_e], X_tr_e[mask_tr_e],
                                 family=sm.families.Poisson()).fit(maxiter=50)
                pred_event = m_event.predict(X_te_e[mask_te_e])
                y_te_e = y_test[mask_te_e].values
                ll_event = np.sum(y_te_e * np.log(np.maximum(pred_event, 1e-10)) - pred_event)

                base_logliks.append(ll_base)
                event_logliks.append(ll_event)
                improvements.append(ll_event - ll_base)

                if exp_col in m_event.params.index:
                    fold_rrs.append(float(np.exp(m_event.params[exp_col])))
            except Exception:
                pass

        if improvements:
            imp = np.array(improvements)
            t_imp, p_imp = stats.ttest_1samp(imp, 0)
            try:
                w_imp, wp_imp = stats.wilcoxon(imp, alternative="greater")
            except Exception:
                w_imp, wp_imp = np.nan, np.nan

            n_improved = np.sum(imp > 0)
            rr_arr = np.array(fold_rrs) if fold_rrs else np.array([])

            section[hyp_name] = {
                "n_folds": len(imp),
                "mean_ll_improvement": float(imp.mean()),
                "median_ll_improvement": float(np.median(imp)),
                "n_folds_improved": int(n_improved),
                "pct_improved": float(100 * n_improved / len(imp)),
                "paired_t_p": float(p_imp),
                "wilcoxon_p": float(wp_imp) if not np.isnan(wp_imp) else None,
                "fold_rrs": {
                    "mean": float(rr_arr.mean()) if len(rr_arr) > 0 else None,
                    "std": float(rr_arr.std()) if len(rr_arr) > 0 else None,
                    "all_lt_1": bool(np.all(rr_arr < 1)) if len(rr_arr) > 0 else None,
                    "sign_test_p": float(2 * min(
                        stats.binom.cdf(np.sum(rr_arr < 1), len(rr_arr), 0.5),
                        1 - stats.binom.cdf(np.sum(rr_arr < 1) - 1, len(rr_arr), 0.5)
                    )) if len(rr_arr) > 0 else None,
                },
            }
            LOG.info("  %s: mean_LL_imp=%.3f, %d/%d improved, paired_t p=%.4f, Wilcoxon p=%s",
                     hyp_name, imp.mean(), n_improved, len(imp), p_imp,
                     "%.4f" % wp_imp if not np.isnan(wp_imp) else "N/A")
            if len(rr_arr) > 0:
                LOG.info("    Fold RRs: mean=%.3f, std=%.3f, all<1=%s",
                         rr_arr.mean(), rr_arr.std(), np.all(rr_arr < 1))

    results["part3_loocv_deviance"] = section
    save(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 4: Isolated Events + Stratum-Width Sensitivity
# ═══════════════════════════════════════════════════════════════════════════════
def part4_isolated_events(results):
    """
    1. Identify isolated events (>14 days from nearest other event).
    2. Run MH case-crossover on isolated events only.
    3. Stratum-width sensitivity: 7, 10, 15, 21 days.
    4. Test pre vs post asymmetry with isolated events.
    """
    LOG.info("=" * 70)
    LOG.info("PART 4: Isolated Events + Stratum-Width Sensitivity")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    event_dates = w.index[w["geo_event"] == 1].values
    # Find isolated events (>14 days from any other event)
    isolated = []
    for i, ed in enumerate(event_dates):
        is_isolated = True
        for j, other in enumerate(event_dates):
            if i == j:
                continue
            gap = abs((ed - other) / np.timedelta64(1, 'D'))
            if gap <= 14:
                is_isolated = False
                break
        if is_isolated:
            isolated.append(ed)

    LOG.info("  %d/%d events are isolated (>14d gap)", len(isolated), len(event_dates))

    # Create isolated event indicators
    w["isolated_post_1_3d"] = 0
    w["isolated_pre_7_3d"] = 0
    for ed in isolated:
        ed_ts = pd.Timestamp(ed)
        mask_post = (w.index > ed_ts) & (w.index <= ed_ts + pd.Timedelta(days=3))
        w.loc[mask_post, "isolated_post_1_3d"] = 1
        mask_pre = (w.index >= ed_ts - pd.Timedelta(days=7)) & \
                   (w.index <= ed_ts - pd.Timedelta(days=3))
        w.loc[mask_pre, "isolated_pre_7_3d"] = 1

    # Stratum-width sensitivity for isolated events
    for sw in [7, 10, 15, 21]:
        w["_strat"] = w["winter_id"].astype(str) + "_" + \
                      (w["day_of_season"] // sw).astype(int).astype(str)

        for out_name, out_col in [("dry_natural", "dry_natural_size_1234"),
                                   ("all_natural", "aai_all_natural")]:
            if out_col not in w.columns:
                continue
            y = w[out_col].dropna()
            valid = w.loc[y.index]

            # Post-event
            res_post = mh_rate_ratio(valid, out_col, "isolated_post_1_3d", "_strat")
            # Pre-event
            res_pre = mh_rate_ratio(valid, out_col, "isolated_pre_7_3d", "_strat")

            key = f"{out_name}_sw{sw}"
            result = {"stratum_width": sw}
            if res_post:
                result["post_1_3d"] = res_post
            if res_pre:
                result["pre_7_3d"] = res_pre

            # Asymmetry: is post stronger than pre?
            if res_post and res_pre:
                result["post_rr"] = res_post["rr"]
                result["pre_rr"] = res_pre["rr"]
                result["asymmetry_ratio"] = res_post["rr"] / res_pre["rr"] if res_pre["rr"] > 0 else None

            section[key] = result

            if res_post:
                LOG.info("  %s sw=%d: POST RR=%.3f p=%.4f | PRE RR=%.3f p=%.4f",
                         out_name, sw,
                         res_post["rr"], res_post["p"],
                         res_pre["rr"] if res_pre else float('nan'),
                         res_pre["p"] if res_pre else float('nan'))

    # Also run all events (not just isolated) across strata widths for comparison
    section["n_isolated"] = len(isolated)
    section["n_total_events"] = len(event_dates)

    # Summary across stratum widths
    dry_post_rrs = []
    dry_pre_rrs = []
    for sw in [7, 10, 15, 21]:
        key = f"dry_natural_sw{sw}"
        if key in section:
            if "post_1_3d" in section[key]:
                dry_post_rrs.append(section[key]["post_1_3d"]["rr"])
            if "pre_7_3d" in section[key]:
                dry_pre_rrs.append(section[key]["pre_7_3d"]["rr"])

    if dry_post_rrs and dry_pre_rrs:
        section["dry_summary"] = {
            "post_rrs_across_strata": dry_post_rrs,
            "pre_rrs_across_strata": dry_pre_rrs,
            "post_mean": float(np.mean(dry_post_rrs)),
            "pre_mean": float(np.mean(dry_pre_rrs)),
            "post_stronger_in_n_widths": int(sum(1 for po, pr in zip(dry_post_rrs, dry_pre_rrs)
                                                  if po < pr)),
        }
        LOG.info("  Dry isolated: post mean RR=%.3f, pre mean RR=%.3f, post<pre in %d/4 widths",
                 np.mean(dry_post_rrs), np.mean(dry_pre_rrs),
                 sum(1 for po, pr in zip(dry_post_rrs, dry_pre_rrs) if po < pr))

    results["part4_isolated_events"] = section
    save(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 5: Stratum FE Poisson + Placebo Retest
# ═══════════════════════════════════════════════════════════════════════════════
def part5_poisson_fe_placebo(results):
    """
    Replace NB GLM with Poisson + stratum FE + robust SE.
    Then rerun 10 placebo event sets to verify false-positive rate.
    """
    LOG.info("=" * 70)
    LOG.info("PART 5: Poisson Stratum FE + Placebo Retest")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}
    rng = np.random.RandomState(456)

    w["period_15"] = (w["day_of_season"] // 15).astype(int)
    w["stratum"] = w["winter_id"].astype(str) + "_p" + w["period_15"].astype(str)

    # Create stratum dummies
    strat_dummies = pd.get_dummies(w["stratum"], prefix="s", drop_first=True)

    for out_name, out_col in [("dry_natural", "dry_natural_size_1234"),
                               ("all_natural", "aai_all_natural")]:
        if out_col not in w.columns:
            continue

        y = w[out_col].dropna()
        idx = y.index

        for exp_name, exp_col in [("post_1_3d", "post_event_1_3d")]:
            # Fit Poisson with stratum FE
            X = pd.concat([
                w.loc[idx, [exp_col]],
                strat_dummies.loc[idx]
            ], axis=1)
            X = sm.add_constant(X)
            mask = y.notna() & X.notna().all(axis=1)
            y_c = y[mask].astype(float)
            X_c = X[mask].astype(float)

            if len(y_c) < 100:
                continue

            try:
                model = sm.GLM(y_c, X_c, family=sm.families.Poisson())
                result = model.fit(maxiter=100, cov_type="HC1")  # robust SE
                coef = float(result.params[exp_col])
                rr = float(np.exp(coef))
                p = float(result.pvalues[exp_col])
                ci = result.conf_int().loc[exp_col]

                key = f"{out_name}_{exp_name}"
                section[key] = {
                    "rr": rr, "p": p,
                    "ci_lo": float(np.exp(ci[0])),
                    "ci_hi": float(np.exp(ci[1])),
                    "method": "Poisson GLM + stratum FE + HC1 robust SE",
                    "n_strata_dummies": int(strat_dummies.shape[1]),
                }
                LOG.info("  PoissonFE %s: RR=%.3f [%.3f-%.3f] p=%.4f",
                         key, rr, np.exp(ci[0]), np.exp(ci[1]), p)

                # Placebo test: 10 random event sets
                placebo_ps = []
                placebo_rrs = []
                n_events = int(w[exp_col].sum())
                winter_dates = w.index.values

                for trial in range(10):
                    placebo_dates = rng.choice(winter_dates, size=int(w["geo_event"].sum()),
                                               replace=False)
                    w["_placebo"] = 0
                    for pd_date in placebo_dates:
                        pm = (w.index > pd_date) & (w.index <= pd_date + pd.Timedelta(days=3))
                        w.loc[pm, "_placebo"] = 1

                    X_p = pd.concat([
                        w.loc[idx, ["_placebo"]],
                        strat_dummies.loc[idx]
                    ], axis=1)
                    X_p = sm.add_constant(X_p)
                    mask_p = y.notna() & X_p.notna().all(axis=1)

                    try:
                        m_p = sm.GLM(y[mask_p], X_p[mask_p].astype(float),
                                     family=sm.families.Poisson())
                        r_p = m_p.fit(maxiter=50, cov_type="HC1")
                        placebo_rrs.append(float(np.exp(r_p.params["_placebo"])))
                        placebo_ps.append(float(r_p.pvalues["_placebo"]))
                    except Exception:
                        pass

                n_sig = sum(1 for p in placebo_ps if p < 0.05)
                section[key + "_placebo"] = {
                    "n_trials": len(placebo_ps),
                    "n_significant": n_sig,
                    "mean_rr": float(np.mean(placebo_rrs)) if placebo_rrs else None,
                    "mean_p": float(np.mean(placebo_ps)) if placebo_ps else None,
                    "verdict": "PASS" if n_sig <= 1 else "FAIL",
                }
                LOG.info("    Placebo: %d/%d significant (PASS if ≤1)",
                         n_sig, len(placebo_ps))

            except Exception as e:
                LOG.warning("  PoissonFE %s_%s FAILED: %s", out_name, exp_name, e)

    results["part5_poisson_fe"] = section
    save(results)
    del w, strat_dummies; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 6: SSW Expanded — Longer Windows + Weather-Matched
# ═══════════════════════════════════════════════════════════════════════════════
def part6_ssw_expanded(results):
    """
    Expand SSW analysis:
    1. Multiple post-SSW windows (0-7, 0-15, 0-30, 15-30)
    2. Pre-SSW comparison for causal direction
    3. Bootstrap CI on matched diffs
    """
    LOG.info("=" * 70)
    LOG.info("PART 6: SSW Expanded Windows + Bootstrap CI")
    LOG.info("=" * 70)

    w = load_panel()
    ssw_cat = load_ssw()
    section = {}
    rng = np.random.RandomState(789)

    ssw_winter_ids = set()
    for sd in ssw_cat.index:
        if sd.month >= 11:
            ssw_winter_ids.add(f"{sd.year}/{sd.year+1}")
        elif sd.month <= 3:
            ssw_winter_ids.add(f"{sd.year-1}/{sd.year}")
    non_ssw = [wid for wid in w["winter_id"].dropna().unique()
               if wid not in ssw_winter_ids]

    windows = [("pre_15_0", -15, 0), ("post_0_7", 0, 7),
               ("post_0_15", 0, 15), ("post_0_30", 0, 30), ("post_15_30", 15, 30)]

    for out_name, out_col in [("all_natural", "aai_all_natural"),
                               ("dry_natural", "dry_natural_size_1234"),
                               ("norway", "norway_aval_count")]:
        if out_col not in w.columns:
            continue

        for win_name, d_lo, d_hi in windows:
            diffs = []
            for sd in ssw_cat.index:
                if sd.month >= 11:
                    sw = f"{sd.year}/{sd.year+1}"
                    dos = (sd - pd.Timestamp(sd.year, 11, 1)).days
                elif sd.month <= 3:
                    sw = f"{sd.year-1}/{sd.year}"
                    dos = (sd - pd.Timestamp(sd.year-1, 11, 1)).days
                else:
                    continue
                if sw not in w["winter_id"].values:
                    continue

                ssw_win = w[(w.index >= sd + pd.Timedelta(days=d_lo)) &
                            (w.index < sd + pd.Timedelta(days=d_hi))]
                ssw_mean = ssw_win[out_col].mean()
                if np.isnan(ssw_mean):
                    continue

                ctrl_means = []
                for ctrl_wid in non_ssw:
                    ctrl = w[w["winter_id"] == ctrl_wid]
                    if len(ctrl) == 0:
                        continue
                    matched = ctrl[(ctrl["day_of_season"] >= dos + d_lo - 3) &
                                   (ctrl["day_of_season"] <= dos + d_hi + 3)]
                    if len(matched) >= max(3, (d_hi - d_lo) // 2):
                        cm = matched[out_col].mean()
                        if not np.isnan(cm):
                            ctrl_means.append(cm)
                if len(ctrl_means) >= 3:
                    diffs.append(ssw_mean - np.mean(ctrl_means))

            if len(diffs) >= 5:
                d = np.array(diffs)
                t_stat, t_p = stats.ttest_1samp(d, 0)
                try:
                    w_stat, w_p = stats.wilcoxon(d, alternative="two-sided")
                except Exception:
                    w_stat, w_p = np.nan, np.nan
                n_neg = np.sum(d < 0)
                sign_p = 2 * min(stats.binom.cdf(n_neg, len(d), 0.5),
                                  1 - stats.binom.cdf(n_neg - 1, len(d), 0.5))

                # Bootstrap CI
                boot = [rng.choice(d, len(d), replace=True).mean() for _ in range(10000)]
                boot = np.array(boot)

                key = f"{out_name}_{win_name}"
                section[key] = {
                    "n": len(d), "mean": float(d.mean()), "median": float(np.median(d)),
                    "t_p": float(t_p),
                    "wilcoxon_p": float(w_p) if not np.isnan(w_p) else None,
                    "sign_test": {"n_neg": int(n_neg), "n_total": len(d), "p": float(sign_p)},
                    "boot_ci": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                }
                LOG.info("  %s: mean=%.3f, %d/%d neg, t_p=%.4f, Wilcoxon=%.4f, sign=%.4f",
                         key, d.mean(), n_neg, len(d), t_p,
                         w_p if not np.isnan(w_p) else -1, sign_p)

    results["part6_ssw_expanded"] = section
    save(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    results = {}
    if RESULTS_FILE.exists():
        try:
            results = json.loads(open(RESULTS_FILE, encoding="utf-8").read())
        except Exception:
            results = {}

    parts = [
        ("part1_perm_spec_curve", part1_perm_spec_curve),
        ("part2_ssw_battery", part2_ssw_all_natural),
        ("part3_loocv_deviance", part3_loocv_deviance),
        ("part4_isolated_events", part4_isolated_events),
        ("part5_poisson_fe", part5_poisson_fe_placebo),
        ("part6_ssw_expanded", part6_ssw_expanded),
    ]

    for name, func in parts:
        if name in results:
            LOG.info("Skipping %s (already done)", name)
            continue
        try:
            results = func(results)
            gc.collect()
        except Exception as e:
            import traceback
            LOG.error("FAILED %s: %s", name, e)
            results[name] = {"error": str(e), "traceback": traceback.format_exc()}
            save(results)

    # Final summary
    print("\n" + "=" * 70)
    print("TIER 2 UPGRADE COMPLETE")
    print("=" * 70)
    for name, _ in parts:
        status = "OK" if name in results and "error" not in results.get(name, {}) else "FAILED"
        print(f"  {name}: {status}")


if __name__ == "__main__":
    main()
