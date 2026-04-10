"""
16_robustness_resolution.py — Address all robustness limitations
================================================================
Implements reviewer-proof analysis hierarchy:

PART 1: Case-Crossover Design (within-winter seasonal matching)
PART 2: Within-Winter Anomaly Analysis
PART 3: Winter Fixed-Effects Model
PART 4: Block Bootstrap Permutation (for SSW + dry avalanche)
PART 5: Specification Curve (all reasonable model variants)
PART 6: Falsification Battery (wet, placebo, pre-event, shifted)
PART 7: Prediction-Based LOOCV (deviance gain from event indicators)
PART 8: SSW Matched Comparison (tight within-stratum matching)
"""
import sys, gc, json, logging
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, LOG

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

RESULTS_FILE = RESULTS / "robustness_resolution.json"


def load_panel():
    panel = pd.read_parquet(PROCESSED / "analysis_panel_v2.parquet")
    return panel[panel["is_winter"] == 1].copy()


def save_results(results):
    with open(RESULTS_FILE, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=str)


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: Case-Crossover Design
# ═══════════════════════════════════════════════════════════════════════════════
def part1_case_crossover(results):
    """
    Time-stratified case-crossover: for each event day, controls are
    non-event days in the same 15-day stratum of the same winter.
    Uses Mantel-Haenszel rate ratio — non-parametric, no model assumptions.
    """
    LOG.info("=" * 70)
    LOG.info("PART 1: Case-Crossover Design")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    # Create strata: (winter_id, 15-day period)
    w["period"] = (w["day_of_season"] // 15).astype(int)
    w["stratum"] = w["winter_id"].astype(str) + "_p" + w["period"].astype(str)

    for outcome_name, outcome_col in [
        ("all_natural", "aai_all_natural"),
        ("dry_natural", "dry_natural_size_1234"),
        ("wet_natural", "wet_natural_size_1234"),
    ]:
        if outcome_col not in w.columns:
            continue

        for exposure_name, exposure_col in [
            ("geomag_1_3d", "post_event_1_3d"),
            ("geomag_5_21d", "post_event_5_21d"),
        ]:
            if exposure_col not in w.columns:
                continue

            y = w[outcome_col].dropna()
            valid = w.loc[y.index].copy()

            # Mantel-Haenszel stratified rate ratio
            exposed_sum = 0.0
            unexposed_sum = 0.0
            exposed_weight = 0.0
            unexposed_weight = 0.0
            n_strata_used = 0

            for stratum_id, group in valid.groupby("stratum"):
                exp = group[group[exposure_col] == 1]
                unexp = group[group[exposure_col] == 0]
                n_exp = len(exp)
                n_unexp = len(unexp)
                if n_exp == 0 or n_unexp == 0:
                    continue
                n_total = n_exp + n_unexp
                y_exp = exp[outcome_col].sum()
                y_unexp = unexp[outcome_col].sum()

                # MH components
                exposed_sum += y_exp * n_unexp / n_total
                unexposed_sum += y_unexp * n_exp / n_total
                n_strata_used += 1

            if unexposed_sum > 0 and exposed_sum > 0:
                mh_rr = exposed_sum / unexposed_sum
                # Variance via Greenland-Robins
                ln_rr = np.log(mh_rr)
                # Approximate SE using Rothman formula
                se_ln_rr = np.sqrt(1.0 / exposed_sum + 1.0 / unexposed_sum)
                z = ln_rr / se_ln_rr
                p = 2 * stats.norm.sf(abs(z))

                key = "%s_%s" % (outcome_name, exposure_name)
                section[key] = {
                    "mh_rate_ratio": float(mh_rr),
                    "ci_low": float(np.exp(ln_rr - 1.96 * se_ln_rr)),
                    "ci_high": float(np.exp(ln_rr + 1.96 * se_ln_rr)),
                    "z": float(z),
                    "p_value": float(p),
                    "n_strata": n_strata_used,
                    "method": "Mantel-Haenszel stratified by winter x 15-day period",
                }
                LOG.info("  MH %s: RR=%.3f [%.3f-%.3f] z=%.2f p=%.4f (%d strata)",
                         key, mh_rr, np.exp(ln_rr - 1.96*se_ln_rr),
                         np.exp(ln_rr + 1.96*se_ln_rr), z, p, n_strata_used)

    # SSW case-crossover: SSW within-stratum matching
    LOG.info("  SSW case-crossover...")
    ssw_cat = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw_cat.index = ssw_cat.index.tz_localize(None) if hasattr(ssw_cat.index, 'tz') and ssw_cat.index.tz else ssw_cat.index

    # Create SSW exposure windows
    for window, d_lo, d_hi in [("ssw_0_15", 0, 15), ("ssw_15_30", 15, 30)]:
        w[window] = 0
        for sd in ssw_cat.index:
            mask = (w.index >= sd + pd.Timedelta(days=d_lo)) & \
                   (w.index < sd + pd.Timedelta(days=d_hi))
            w.loc[mask, window] = 1

        for outcome_name, outcome_col in [
            ("all_natural", "aai_all_natural"),
            ("dry_natural", "dry_natural_size_1234"),
            ("norway", "norway_aval_count"),
        ]:
            if outcome_col not in w.columns:
                continue
            y = w[outcome_col].dropna()
            valid = w.loc[y.index].copy()

            exposed_sum = 0.0
            unexposed_sum = 0.0
            n_strata_used = 0

            for stratum_id, group in valid.groupby("stratum"):
                exp = group[group[window] == 1]
                unexp = group[group[window] == 0]
                if len(exp) == 0 or len(unexp) == 0:
                    continue
                n_total = len(exp) + len(unexp)
                y_exp = exp[outcome_col].sum()
                y_unexp = unexp[outcome_col].sum()
                exposed_sum += y_exp * len(unexp) / n_total
                unexposed_sum += y_unexp * len(exp) / n_total
                n_strata_used += 1

            if unexposed_sum > 0 and exposed_sum > 0:
                mh_rr = exposed_sum / unexposed_sum
                ln_rr = np.log(mh_rr)
                se_ln_rr = np.sqrt(1.0 / exposed_sum + 1.0 / unexposed_sum)
                z = ln_rr / se_ln_rr
                p = 2 * stats.norm.sf(abs(z))
                key = "%s_%s" % (outcome_name, window)
                section[key] = {
                    "mh_rate_ratio": float(mh_rr),
                    "ci_low": float(np.exp(ln_rr - 1.96 * se_ln_rr)),
                    "ci_high": float(np.exp(ln_rr + 1.96 * se_ln_rr)),
                    "z": float(z),
                    "p_value": float(p),
                    "n_strata": n_strata_used,
                }
                LOG.info("  MH %s: RR=%.3f [%.3f-%.3f] p=%.4f (%d strata)",
                         key, mh_rr, np.exp(ln_rr - 1.96*se_ln_rr),
                         np.exp(ln_rr + 1.96*se_ln_rr), p, n_strata_used)

    results["part1_case_crossover"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: Within-Winter Anomaly Analysis
# ═══════════════════════════════════════════════════════════════════════════════
def part2_anomaly_analysis(results):
    """
    Compute within-winter anomalies: daily value minus that winter's
    running mean (21-day centered, excluding event windows).
    Then test if event-day anomalies are negative.
    """
    LOG.info("=" * 70)
    LOG.info("PART 2: Within-Winter Anomaly Analysis")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    for outcome_name, outcome_col in [
        ("all_natural", "aai_all_natural"),
        ("dry_natural", "dry_natural_size_1234"),
        ("wet_natural", "wet_natural_size_1234"),
        ("norway", "norway_aval_count"),
    ]:
        if outcome_col not in w.columns:
            continue

        anomalies = pd.Series(np.nan, index=w.index)

        for wid, group in w.groupby("winter_id"):
            if wid is None:
                continue
            y_winter = group[outcome_col]
            if y_winter.notna().sum() < 30:
                continue

            # Compute running mean excluding event windows
            non_event = y_winter.copy()
            non_event[group["post_event_0_30d"] == 1] = np.nan
            baseline = non_event.rolling(21, center=True, min_periods=7).mean()
            # Fill gaps in baseline with winter mean
            baseline = baseline.fillna(y_winter.mean())
            if (baseline > 0).any():
                anom = y_winter - baseline
                anomalies.loc[group.index] = anom

        w["anomaly_" + outcome_name] = anomalies

        # Test: are post-event anomalies negative?
        for exposure_name, exposure_col in [
            ("geomag_1_3d", "post_event_1_3d"),
            ("geomag_5_21d", "post_event_5_21d"),
        ]:
            anom_col = "anomaly_" + outcome_name
            valid = w[[anom_col, exposure_col]].dropna()
            exp = valid[valid[exposure_col] == 1][anom_col]
            unexp = valid[valid[exposure_col] == 0][anom_col]

            if len(exp) > 10 and len(unexp) > 10:
                t, p = stats.ttest_ind(exp, unexp, equal_var=False)
                # Also Wilcoxon rank-sum (non-parametric)
                u_stat, p_mann = stats.mannwhitneyu(exp, unexp, alternative="two-sided")

                key = "%s_%s" % (outcome_name, exposure_name)
                section[key] = {
                    "exposed_mean_anomaly": float(exp.mean()),
                    "unexposed_mean_anomaly": float(unexp.mean()),
                    "diff": float(exp.mean() - unexp.mean()),
                    "t_stat": float(t),
                    "t_p_value": float(p),
                    "mann_whitney_U": float(u_stat),
                    "mann_whitney_p": float(p_mann),
                    "n_exposed": int(len(exp)),
                    "n_unexposed": int(len(unexp)),
                }
                LOG.info("  Anomaly %s: exp=%.4f unexp=%.4f t=%.2f p=%.4f MW_p=%.4f",
                         key, exp.mean(), unexp.mean(), t, p, p_mann)

    # SSW anomalies
    ssw_cat = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw_cat.index = ssw_cat.index.tz_localize(None) if hasattr(ssw_cat.index, 'tz') and ssw_cat.index.tz else ssw_cat.index
    for window, d_lo, d_hi in [("ssw_0_15", 0, 15), ("ssw_15_30", 15, 30)]:
        w[window] = 0
        for sd in ssw_cat.index:
            mask = (w.index >= sd + pd.Timedelta(days=d_lo)) & \
                   (w.index < sd + pd.Timedelta(days=d_hi))
            w.loc[mask, window] = 1

        for outcome_name in ["all_natural", "dry_natural", "norway"]:
            anom_col = "anomaly_" + outcome_name
            if anom_col not in w.columns:
                continue
            valid = w[[anom_col, window]].dropna()
            exp = valid[valid[window] == 1][anom_col]
            unexp = valid[valid[window] == 0][anom_col]
            if len(exp) > 5 and len(unexp) > 5:
                t, p = stats.ttest_ind(exp, unexp, equal_var=False)
                u_stat, p_mann = stats.mannwhitneyu(exp, unexp, alternative="two-sided")
                key = "%s_%s" % (outcome_name, window)
                section[key] = {
                    "exposed_mean_anomaly": float(exp.mean()),
                    "unexposed_mean_anomaly": float(unexp.mean()),
                    "diff": float(exp.mean() - unexp.mean()),
                    "t_stat": float(t),
                    "t_p_value": float(p),
                    "mann_whitney_p": float(p_mann),
                    "n_exposed": int(len(exp)),
                }
                LOG.info("  SSW Anomaly %s: exp=%.4f unexp=%.4f t=%.2f p=%.4f",
                         key, exp.mean(), unexp.mean(), t, p)

    results["part2_anomaly"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: Winter Fixed-Effects Model
# ═══════════════════════════════════════════════════════════════════════════════
def part3_winter_fixed_effects(results):
    """
    NB GLM with winter dummy variables (fixed effects).
    Effect identified PURELY from within-winter variation.
    """
    LOG.info("=" * 70)
    LOG.info("PART 3: Winter Fixed-Effects Model")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    # Create winter dummies
    winter_dummies = pd.get_dummies(w["winter_id"], prefix="w", drop_first=True)

    for outcome_name, outcome_col in [
        ("all_natural", "aai_all_natural"),
        ("dry_natural", "dry_natural_size_1234"),
        ("norway", "norway_aval_count"),
    ]:
        if outcome_col not in w.columns:
            continue
        y = w[outcome_col].dropna()
        idx = y.index

        for exposure_name, exposure_col in [
            ("geomag_1_3d", "post_event_1_3d"),
            ("geomag_5_21d", "post_event_5_21d"),
        ]:
            # Base covariates + winter FE
            X = pd.concat([
                w.loc[idx, [exposure_col, "day_of_season", "day_of_season_sq"]],
                winter_dummies.loc[idx]
            ], axis=1)
            X = sm.add_constant(X)
            mask = y.notna() & X.notna().all(axis=1)
            y_c = y[mask].astype(float)
            X_c = X[mask].astype(float)

            if len(y_c) < 100:
                continue

            try:
                model = sm.GLM(y_c, X_c, family=sm.families.NegativeBinomial())
                result = model.fit(maxiter=100, method="IRLS")
                coef = float(result.params[exposure_col])
                se = float(result.bse[exposure_col])
                z = float(result.tvalues[exposure_col])
                p = float(result.pvalues[exposure_col])
                rr = float(np.exp(coef))
                ci = result.conf_int().loc[exposure_col]

                key = "%s_%s" % (outcome_name, exposure_name)
                section[key] = {
                    "rr": rr,
                    "ci_low": float(np.exp(ci[0])),
                    "ci_high": float(np.exp(ci[1])),
                    "coef": coef,
                    "se": se,
                    "z": z,
                    "p_value": p,
                    "n": int(len(y_c)),
                    "n_winter_dummies": int(winter_dummies.shape[1]),
                    "method": "NB GLM with winter fixed effects",
                }
                LOG.info("  WinterFE %s: RR=%.3f [%.3f-%.3f] p=%.4f",
                         key, rr, np.exp(ci[0]), np.exp(ci[1]), p)
            except Exception as e:
                LOG.warning("  WinterFE %s: FAILED %s",
                            "%s_%s" % (outcome_name, exposure_name), e)

    # SSW with winter FE
    ssw_cat = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw_cat.index = ssw_cat.index.tz_localize(None) if hasattr(ssw_cat.index, 'tz') and ssw_cat.index.tz else ssw_cat.index
    for window, d_lo, d_hi in [("ssw_0_15", 0, 15), ("ssw_15_30", 15, 30)]:
        w[window] = 0
        for sd in ssw_cat.index:
            mask = (w.index >= sd + pd.Timedelta(days=d_lo)) & \
                   (w.index < sd + pd.Timedelta(days=d_hi))
            w.loc[mask, window] = 1

        for outcome_name, outcome_col in [
            ("all_natural", "aai_all_natural"),
            ("dry_natural", "dry_natural_size_1234"),
            ("norway", "norway_aval_count"),
        ]:
            if outcome_col not in w.columns:
                continue
            y = w[outcome_col].dropna()
            idx = y.index
            X = pd.concat([
                w.loc[idx, [window, "day_of_season", "day_of_season_sq"]],
                winter_dummies.loc[idx]
            ], axis=1)
            X = sm.add_constant(X)
            mask = y.notna() & X.notna().all(axis=1)

            try:
                model = sm.GLM(y[mask].astype(float), X[mask].astype(float),
                               family=sm.families.NegativeBinomial())
                result = model.fit(maxiter=100)
                coef = float(result.params[window])
                rr = float(np.exp(coef))
                p = float(result.pvalues[window])
                ci = result.conf_int().loc[window]
                key = "%s_%s" % (outcome_name, window)
                section[key] = {
                    "rr": rr,
                    "ci_low": float(np.exp(ci[0])),
                    "ci_high": float(np.exp(ci[1])),
                    "p_value": p,
                    "n": int(mask.sum()),
                }
                LOG.info("  WinterFE SSW %s: RR=%.3f [%.3f-%.3f] p=%.4f",
                         key, rr, np.exp(ci[0]), np.exp(ci[1]), p)
            except Exception as e:
                LOG.warning("  WinterFE SSW %s: FAILED %s", key, e)

    results["part3_winter_fe"] = section
    save_results(results)
    del w, winter_dummies; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 4: Block Bootstrap Permutation
# ═══════════════════════════════════════════════════════════════════════════════
def part4_block_bootstrap(results):
    """
    Circular block bootstrap: shift event dates by random offset
    within each winter, preserving temporal structure.
    Tests SSW, dry avalanche, and stratospheric pathway.
    """
    LOG.info("=" * 70)
    LOG.info("PART 4: Block Bootstrap Permutation (2000 iterations)")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}
    rng = np.random.RandomState(42)
    N_BOOT = 2000

    # For each hypothesis, compute observed statistic, then bootstrap null
    hypotheses = [
        ("dry_geomag_1_3d", "dry_natural_size_1234", "post_event_1_3d"),
        ("all_geomag_5_21d", "aai_all_natural", "post_event_5_21d"),
    ]

    for hyp_name, outcome_col, exposure_col in hypotheses:
        if outcome_col not in w.columns or exposure_col not in w.columns:
            continue

        y = w[outcome_col].dropna()
        valid = w.loc[y.index]

        # Observed: mean difference (exposed - unexposed)
        exp_mean = y[valid[exposure_col] == 1].mean()
        unexp_mean = y[valid[exposure_col] == 0].mean()
        obs_diff = exp_mean - unexp_mean

        # Bootstrap: circular shift within each winter
        null_diffs = []
        for b in range(N_BOOT):
            shifted_exposure = pd.Series(0, index=valid.index)
            for wid, group in valid.groupby("winter_id"):
                if wid is None:
                    continue
                exp_vals = group[exposure_col].values
                shift = rng.randint(0, len(group))
                shifted = np.roll(exp_vals, shift)
                shifted_exposure.loc[group.index] = shifted

            null_exp = y[shifted_exposure == 1].mean()
            null_unexp = y[shifted_exposure == 0].mean()
            if not np.isnan(null_exp) and not np.isnan(null_unexp):
                null_diffs.append(null_exp - null_unexp)

        if null_diffs:
            null_arr = np.array(null_diffs)
            # Two-sided p-value
            p_two = np.mean(np.abs(null_arr) >= abs(obs_diff))
            # One-sided (testing for decrease)
            p_one = np.mean(null_arr <= obs_diff)

            section[hyp_name] = {
                "observed_diff": float(obs_diff),
                "null_mean": float(null_arr.mean()),
                "null_std": float(null_arr.std()),
                "p_two_sided": float(p_two),
                "p_one_sided_decrease": float(p_one),
                "n_bootstrap": len(null_diffs),
                "method": "Circular block bootstrap within winter",
            }
            LOG.info("  Block boot %s: obs=%.4f, null_mean=%.4f, p_two=%.4f, p_one=%.4f",
                     hyp_name, obs_diff, null_arr.mean(), p_two, p_one)

    # SSW bootstrap: shuffle which winters get SSW
    LOG.info("  SSW bootstrap (reassigning SSW to random winters)...")
    ssw_cat = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw_cat.index = ssw_cat.index.tz_localize(None) if hasattr(ssw_cat.index, 'tz') and ssw_cat.index.tz else ssw_cat.index

    # Create observed SSW 0-15d indicator
    w["ssw_0_15"] = 0
    for sd in ssw_cat.index:
        mask = (w.index >= sd) & (w.index < sd + pd.Timedelta(days=15))
        w.loc[mask, "ssw_0_15"] = 1

    for outcome_name, outcome_col in [
        ("all_natural", "aai_all_natural"),
        ("dry_natural", "dry_natural_size_1234"),
        ("norway", "norway_aval_count"),
    ]:
        if outcome_col not in w.columns:
            continue
        y = w[outcome_col].dropna()
        valid_idx = y.index

        # Observed MH-like statistic
        obs_exp = y[w.loc[valid_idx, "ssw_0_15"] == 1].mean()
        obs_unexp = y[w.loc[valid_idx, "ssw_0_15"] == 0].mean()
        obs_ratio = obs_exp / obs_unexp if obs_unexp > 0 else np.nan

        # Compute within-stratum observed statistic
        w["period"] = (w["day_of_season"] // 15).astype(int)
        obs_strat_diffs = []
        for wid, group in w.loc[valid_idx].groupby("winter_id"):
            if wid is None:
                continue
            for pid, subg in group.groupby("period"):
                exp = subg[subg["ssw_0_15"] == 1][outcome_col]
                unexp = subg[subg["ssw_0_15"] == 0][outcome_col]
                if len(exp) > 0 and len(unexp) > 0:
                    obs_strat_diffs.append(exp.mean() - unexp.mean())

        obs_strat_mean = np.mean(obs_strat_diffs) if obs_strat_diffs else 0

        # Bootstrap: randomly assign SSW-like windows to different positions
        winters = w["winter_id"].dropna().unique()
        n_ssw_winters = w.loc[w["ssw_0_15"] == 1, "winter_id"].nunique()

        null_strat_means = []
        for b in range(N_BOOT):
            w["ssw_boot"] = 0
            # Pick random winters and random day_of_season for SSW placement
            boot_winters = rng.choice(winters, size=n_ssw_winters, replace=False)
            for bw in boot_winters:
                wgroup = w[w["winter_id"] == bw]
                if len(wgroup) < 15:
                    continue
                start_idx = rng.randint(0, max(1, len(wgroup) - 15))
                boot_dates = wgroup.index[start_idx:start_idx + 15]
                w.loc[boot_dates, "ssw_boot"] = 1

            strat_diffs = []
            for wid, group in w.loc[valid_idx].groupby("winter_id"):
                if wid is None:
                    continue
                for pid, subg in group.groupby("period"):
                    exp = subg[subg["ssw_boot"] == 1][outcome_col]
                    unexp = subg[subg["ssw_boot"] == 0][outcome_col]
                    if len(exp) > 0 and len(unexp) > 0:
                        strat_diffs.append(exp.mean() - unexp.mean())
            if strat_diffs:
                null_strat_means.append(np.mean(strat_diffs))

        if null_strat_means:
            null_arr = np.array(null_strat_means)
            p_two = np.mean(np.abs(null_arr) >= abs(obs_strat_mean))
            p_one = np.mean(null_arr <= obs_strat_mean)

            key = "%s_ssw_0_15" % outcome_name
            section[key] = {
                "observed_strat_diff": float(obs_strat_mean),
                "null_mean": float(null_arr.mean()),
                "null_std": float(null_arr.std()),
                "p_two_sided": float(p_two),
                "p_one_sided_decrease": float(p_one),
                "n_bootstrap": len(null_strat_means),
                "method": "SSW placement bootstrap with stratified statistic",
            }
            LOG.info("  SSW boot %s: obs=%.4f, null=%.4f, p_two=%.4f",
                     key, obs_strat_mean, null_arr.mean(), p_two)

    results["part4_block_bootstrap"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 5: Specification Curve
# ═══════════════════════════════════════════════════════════════════════════════
def part5_specification_curve(results):
    """
    Run ALL reasonable model variants. Report distribution of effect sizes.
    Follows Simonsohn et al. (2020) specification curve analysis.
    """
    LOG.info("=" * 70)
    LOG.info("PART 5: Specification Curve Analysis")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    specifications = []

    # Specification dimensions:
    # Outcome: [all_natural, dry_natural]
    # Exposure: [1_3d, 5_21d]
    # Season control: [quadratic, cubic, winter_FE]
    # Confounders: [none, nao, nao+qbo, nao+qbo+f107, full]
    # Model: [NB_GLM]

    outcomes = [
        ("all_natural", "aai_all_natural"),
        ("dry_natural", "dry_natural_size_1234"),
    ]
    exposures = [
        ("1_3d", "post_event_1_3d"),
        ("5_21d", "post_event_5_21d"),
    ]
    season_controls = {
        "quadratic": ["day_of_season", "day_of_season_sq"],
        "cubic": None,  # will create
    }
    confounder_sets = {
        "none": [],
        "nao": ["nao_daily"],
        "nao_qbo": ["nao_daily", "qbo_u50"],
        "nao_qbo_f107": ["nao_daily", "qbo_u50", "f107"],
        "full": ["nao_daily", "qbo_u50", "f107", "snotel_swe_mean", "modis_snow_frac"],
    }

    # Create cubic term
    w["day_of_season_cu"] = w["day_of_season"] ** 3

    # Winter dummies for FE specs
    winter_dummies = pd.get_dummies(w["winter_id"], prefix="w", drop_first=True)

    for out_name, out_col in outcomes:
        if out_col not in w.columns:
            continue
        y = w[out_col].dropna()
        idx = y.index

        for exp_name, exp_col in exposures:
            for season_name in ["quadratic", "cubic", "winter_FE"]:
                for conf_name, conf_cols in confounder_sets.items():
                    # Build X
                    cols = [exp_col]
                    if season_name == "quadratic":
                        cols += ["day_of_season", "day_of_season_sq"]
                    elif season_name == "cubic":
                        cols += ["day_of_season", "day_of_season_sq", "day_of_season_cu"]

                    cols += [c for c in conf_cols if c in w.columns]
                    cols = list(dict.fromkeys(cols))  # unique

                    if season_name == "winter_FE":
                        X = pd.concat([w.loc[idx, cols], winter_dummies.loc[idx]], axis=1)
                    else:
                        X = w.loc[idx, cols]

                    X = sm.add_constant(X)
                    mask = y.notna() & X.notna().all(axis=1)
                    y_c = y[mask].astype(float)
                    X_c = X[mask].astype(float)

                    if len(y_c) < 100:
                        continue

                    try:
                        model = sm.GLM(y_c, X_c, family=sm.families.NegativeBinomial())
                        result = model.fit(maxiter=50, method="IRLS")
                        coef = float(result.params[exp_col])
                        rr = float(np.exp(coef))
                        p = float(result.pvalues[exp_col])

                        specifications.append({
                            "outcome": out_name,
                            "exposure": exp_name,
                            "season": season_name,
                            "confounders": conf_name,
                            "rr": rr,
                            "coef": coef,
                            "p_value": p,
                            "n": int(len(y_c)),
                            "significant": p < 0.05,
                            "direction_decrease": rr < 1.0,
                        })
                    except Exception:
                        pass

    if specifications:
        df_spec = pd.DataFrame(specifications)
        section["n_specifications"] = len(df_spec)
        section["n_significant"] = int(df_spec["significant"].sum())
        section["n_decrease"] = int(df_spec["direction_decrease"].sum())
        section["pct_significant"] = float(100 * df_spec["significant"].mean())
        section["pct_decrease"] = float(100 * df_spec["direction_decrease"].mean())
        section["median_rr"] = float(df_spec["rr"].median())
        section["mean_rr"] = float(df_spec["rr"].mean())
        section["rr_range"] = [float(df_spec["rr"].min()), float(df_spec["rr"].max())]

        # By outcome
        for out in df_spec["outcome"].unique():
            sub = df_spec[df_spec["outcome"] == out]
            section["spec_" + out] = {
                "n": int(len(sub)),
                "n_sig": int(sub["significant"].sum()),
                "n_decrease": int(sub["direction_decrease"].sum()),
                "median_rr": float(sub["rr"].median()),
                "rr_range": [float(sub["rr"].min()), float(sub["rr"].max())],
            }

        # Store individual specs
        section["specifications"] = specifications

        LOG.info("  %d specifications: %d significant (%.1f%%), %d decrease (%.1f%%)",
                 len(df_spec), df_spec["significant"].sum(),
                 100 * df_spec["significant"].mean(),
                 df_spec["direction_decrease"].sum(),
                 100 * df_spec["direction_decrease"].mean())
        LOG.info("  Median RR=%.3f, range [%.3f - %.3f]",
                 df_spec["rr"].median(), df_spec["rr"].min(), df_spec["rr"].max())

    results["part5_spec_curve"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 6: Falsification Battery
# ═══════════════════════════════════════════════════════════════════════════════
def part6_falsification(results):
    """
    Tests that SHOULD be null:
    1. Wet avalanches (if dry-specific mechanism is correct)
    2. Pre-event window (7-3 days BEFORE events)
    3. Placebo event dates (random dates)
    4. Summer season
    5. Reversed-sign test (exposure = 0 days, control = post-event)
    """
    LOG.info("=" * 70)
    LOG.info("PART 6: Falsification Battery")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    y_dry = w["dry_natural_size_1234"].dropna() if "dry_natural_size_1234" in w.columns else None
    y_wet = w["wet_natural_size_1234"].dropna() if "wet_natural_size_1234" in w.columns else None
    y_all = w["aai_all_natural"].dropna()

    base_cols = ["day_of_season", "day_of_season_sq"]

    # 1. Wet avalanche control
    if y_wet is not None and len(y_wet) > 100:
        X = w.loc[y_wet.index, ["post_event_1_3d"] + base_cols]
        X = sm.add_constant(X)
        mask = y_wet.notna() & X.notna().all(axis=1)
        try:
            model = sm.GLM(y_wet[mask], X[mask], family=sm.families.NegativeBinomial())
            result = model.fit(maxiter=50)
            p = float(result.pvalues["post_event_1_3d"])
            rr = float(np.exp(result.params["post_event_1_3d"]))
            section["wet_control"] = {"rr": rr, "p_value": p, "expected": "null"}
            LOG.info("  Wet control: RR=%.3f p=%.4f (expected null)", rr, p)
        except Exception as e:
            section["wet_control"] = {"error": str(e)}

    # 2. Pre-event window (should be null — no causal effect before event)
    event_dates = w.index[w["geo_event"] == 1]
    w["pre_event_7_3d"] = 0
    for ed in event_dates:
        mask = (w.index >= ed - pd.Timedelta(days=7)) & \
               (w.index <= ed - pd.Timedelta(days=3))
        w.loc[mask, "pre_event_7_3d"] = 1

    X = w.loc[y_all.index, ["pre_event_7_3d"] + base_cols]
    X = sm.add_constant(X)
    mask = y_all.notna() & X.notna().all(axis=1)
    try:
        model = sm.GLM(y_all[mask], X[mask], family=sm.families.NegativeBinomial())
        result = model.fit(maxiter=50)
        p = float(result.pvalues["pre_event_7_3d"])
        rr = float(np.exp(result.params["pre_event_7_3d"]))
        section["pre_event_control"] = {"rr": rr, "p_value": p, "expected": "null"}
        LOG.info("  Pre-event control: RR=%.3f p=%.4f (expected null)", rr, p)
    except Exception as e:
        section["pre_event_control"] = {"error": str(e)}

    # 3. Placebo events (10 random date sets, average p-value)
    LOG.info("  Running 10 placebo event sets...")
    rng = np.random.RandomState(123)
    n_real_events = w["geo_event"].sum()
    winter_dates = w.index.values
    placebo_rrs = []
    placebo_ps = []

    for trial in range(10):
        placebo_dates = rng.choice(winter_dates, size=n_real_events, replace=False)
        w["placebo_post"] = 0
        for pd_date in placebo_dates:
            mask = (w.index > pd_date) & (w.index <= pd_date + pd.Timedelta(days=3))
            w.loc[mask, "placebo_post"] = 1

        X = w.loc[y_all.index, ["placebo_post"] + base_cols]
        X = sm.add_constant(X)
        mask_v = y_all.notna() & X.notna().all(axis=1)
        try:
            model = sm.GLM(y_all[mask_v], X[mask_v], family=sm.families.NegativeBinomial())
            result = model.fit(maxiter=30)
            placebo_rrs.append(float(np.exp(result.params["placebo_post"])))
            placebo_ps.append(float(result.pvalues["placebo_post"]))
        except Exception:
            pass

    if placebo_rrs:
        section["placebo_events"] = {
            "n_trials": len(placebo_rrs),
            "mean_rr": float(np.mean(placebo_rrs)),
            "std_rr": float(np.std(placebo_rrs)),
            "mean_p": float(np.mean(placebo_ps)),
            "n_significant": int(sum(1 for p in placebo_ps if p < 0.05)),
            "expected": "null (mean RR~1, few significant)",
        }
        LOG.info("  Placebo: mean RR=%.3f, mean p=%.3f, %d/10 significant",
                 np.mean(placebo_rrs), np.mean(placebo_ps),
                 sum(1 for p in placebo_ps if p < 0.05))

    # 4. Summer season control
    w_full = pd.read_parquet(PROCESSED / "analysis_panel_v2.parquet")
    summer = w_full[w_full["is_summer"] == 1].copy()
    y_sum = summer["aai_all_natural"].dropna()
    if len(y_sum) > 100:
        X = summer.loc[y_sum.index, ["post_event_1_3d", "day_of_year"]]
        X = sm.add_constant(X)
        mask = y_sum.notna() & X.notna().all(axis=1)
        try:
            model = sm.GLM(y_sum[mask], X[mask], family=sm.families.NegativeBinomial())
            result = model.fit(maxiter=50)
            p = float(result.pvalues["post_event_1_3d"])
            rr = float(np.exp(result.params["post_event_1_3d"]))
            section["summer_control"] = {"rr": rr, "p_value": p, "expected": "null"}
            LOG.info("  Summer control: RR=%.3f p=%.4f (expected null)", rr, p)
        except Exception as e:
            section["summer_control"] = {"error": str(e)}

    # Also dry-specific placebo test
    if y_dry is not None:
        # Dry avalanche pre-event (should be null)
        X = w.loc[y_dry.index, ["pre_event_7_3d"] + base_cols]
        X = sm.add_constant(X)
        mask = y_dry.notna() & X.notna().all(axis=1)
        try:
            model = sm.GLM(y_dry[mask], X[mask], family=sm.families.NegativeBinomial())
            result = model.fit(maxiter=50)
            p = float(result.pvalues["pre_event_7_3d"])
            rr = float(np.exp(result.params["pre_event_7_3d"]))
            section["dry_pre_event"] = {"rr": rr, "p_value": p, "expected": "null"}
            LOG.info("  Dry pre-event: RR=%.3f p=%.4f (expected null)", rr, p)
        except Exception as e:
            section["dry_pre_event"] = {"error": str(e)}

    # Score
    n_pass = 0
    n_total = 0
    for key, val in section.items():
        if isinstance(val, dict) and "expected" in val and "p_value" in val:
            n_total += 1
            if val["p_value"] > 0.05:
                n_pass += 1
    section["score"] = "%d/%d null as expected" % (n_pass, n_total)
    LOG.info("  Falsification score: %d/%d PASS", n_pass, n_total)

    results["part6_falsification"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 7: Prediction-Based LOOCV
# ═══════════════════════════════════════════════════════════════════════════════
def part7_prediction_loocv(results):
    """
    Leave-one-winter-out: compare predictive deviance
    with and without event indicators.
    """
    LOG.info("=" * 70)
    LOG.info("PART 7: Prediction-Based LOOCV")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    for outcome_name, outcome_col, exposure_col in [
        ("dry_geomag", "dry_natural_size_1234", "post_event_1_3d"),
        ("all_geomag", "aai_all_natural", "post_event_1_3d"),
        ("all_strat", "aai_all_natural", "post_event_5_21d"),
    ]:
        if outcome_col not in w.columns:
            continue

        winters = w["winter_id"].dropna().unique()
        base_devs = []
        event_devs = []
        rrs = []

        for wid in winters:
            train = w[w["winter_id"] != wid]
            test = w[w["winter_id"] == wid]

            y_train = train[outcome_col].dropna()
            y_test = test[outcome_col].dropna()
            if len(y_train) < 100 or len(y_test) < 10:
                continue

            base_cols = ["day_of_season", "day_of_season_sq"]

            # Baseline model (season only)
            X_train_base = sm.add_constant(train.loc[y_train.index, base_cols])
            X_test_base = sm.add_constant(test.loc[y_test.index, base_cols])

            # Event model (season + exposure)
            event_cols = base_cols + [exposure_col]
            X_train_event = sm.add_constant(train.loc[y_train.index, event_cols])
            X_test_event = sm.add_constant(test.loc[y_test.index, event_cols])

            mask_train_b = y_train.notna() & X_train_base.notna().all(axis=1)
            mask_train_e = y_train.notna() & X_train_event.notna().all(axis=1)
            mask_test_b = y_test.notna() & X_test_base.notna().all(axis=1)
            mask_test_e = y_test.notna() & X_test_event.notna().all(axis=1)

            try:
                # Baseline
                m_base = sm.GLM(y_train[mask_train_b],
                                X_train_base[mask_train_b],
                                family=sm.families.NegativeBinomial()).fit(maxiter=50)
                pred_base = m_base.predict(X_test_base[mask_test_b])
                dev_base = np.mean((y_test[mask_test_b] - pred_base) ** 2)

                # Event
                m_event = sm.GLM(y_train[mask_train_e],
                                 X_train_event[mask_train_e],
                                 family=sm.families.NegativeBinomial()).fit(maxiter=50)
                pred_event = m_event.predict(X_test_event[mask_test_e])
                dev_event = np.mean((y_test[mask_test_e] - pred_event) ** 2)

                base_devs.append(dev_base)
                event_devs.append(dev_event)

                if exposure_col in m_event.params.index:
                    rrs.append(float(np.exp(m_event.params[exposure_col])))

            except Exception:
                pass

        if base_devs:
            base_arr = np.array(base_devs)
            event_arr = np.array(event_devs)
            improvement = base_arr - event_arr
            pct_improved = np.mean(improvement > 0) * 100

            section[outcome_name] = {
                "n_folds": len(base_devs),
                "mean_mse_base": float(base_arr.mean()),
                "mean_mse_event": float(event_arr.mean()),
                "mean_improvement": float(improvement.mean()),
                "pct_folds_improved": float(pct_improved),
                "mean_rr_across_folds": float(np.mean(rrs)) if rrs else None,
                "std_rr_across_folds": float(np.std(rrs)) if rrs else None,
                "all_rr_same_direction": all(r < 1 for r in rrs) if rrs else None,
            }
            LOG.info("  LOOCV %s: base_MSE=%.2f event_MSE=%.2f imp=%.2f (%.0f%% folds better)",
                     outcome_name, base_arr.mean(), event_arr.mean(),
                     improvement.mean(), pct_improved)
            if rrs:
                LOG.info("    RR across folds: mean=%.3f std=%.3f all<1=%s",
                         np.mean(rrs), np.std(rrs), all(r < 1 for r in rrs))

    results["part7_loocv"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# PART 8: SSW Matched Comparison
# ═══════════════════════════════════════════════════════════════════════════════
def part8_ssw_matched(results):
    """
    For each SSW event, match to the same day-of-season window in non-SSW
    winters. Compare avalanche activity. This is design-based, not model-based.
    """
    LOG.info("=" * 70)
    LOG.info("PART 8: SSW Matched Comparison")
    LOG.info("=" * 70)

    w = load_panel()
    section = {}

    ssw_cat = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw_cat.index = ssw_cat.index.tz_localize(None) if hasattr(ssw_cat.index, 'tz') and ssw_cat.index.tz else ssw_cat.index

    # Get SSW winters
    ssw_winter_ids = set()
    for sd in ssw_cat.index:
        if sd.month >= 11:
            ssw_winter_ids.add("%d/%d" % (sd.year, sd.year + 1))
        elif sd.month <= 3:
            ssw_winter_ids.add("%d/%d" % (sd.year - 1, sd.year))

    non_ssw_winters = [wid for wid in w["winter_id"].dropna().unique()
                       if wid not in ssw_winter_ids]

    for outcome_name, outcome_col in [
        ("all_natural", "aai_all_natural"),
        ("dry_natural", "dry_natural_size_1234"),
        ("norway", "norway_aval_count"),
    ]:
        if outcome_col not in w.columns:
            continue

        matched_diffs = []
        ssw_event_stats = []

        for sd in ssw_cat.index:
            ssw_dos = None
            ssw_winter = None
            if sd.month >= 11:
                ssw_winter = "%d/%d" % (sd.year, sd.year + 1)
                ssw_dos = (sd - pd.Timestamp(sd.year, 11, 1)).days
            elif sd.month <= 3:
                ssw_winter = "%d/%d" % (sd.year - 1, sd.year)
                ssw_dos = (sd - pd.Timestamp(sd.year - 1, 11, 1)).days

            if ssw_winter not in w["winter_id"].values:
                continue

            # SSW window: 0-15 days post-SSW
            ssw_window = w[(w.index >= sd) & (w.index < sd + pd.Timedelta(days=15))]
            ssw_mean = ssw_window[outcome_col].mean()
            if np.isnan(ssw_mean):
                continue

            # Matched controls: same day_of_season ±3 in non-SSW winters
            control_means = []
            for ctrl_wid in non_ssw_winters:
                ctrl = w[w["winter_id"] == ctrl_wid]
                if len(ctrl) == 0:
                    continue
                # Find matching days
                matched = ctrl[(ctrl["day_of_season"] >= ssw_dos - 3) &
                               (ctrl["day_of_season"] <= ssw_dos + 15 + 3)]
                if len(matched) >= 5:
                    ctrl_mean = matched[outcome_col].mean()
                    if not np.isnan(ctrl_mean):
                        control_means.append(ctrl_mean)

            if len(control_means) >= 3:
                mean_ctrl = np.mean(control_means)
                diff = ssw_mean - mean_ctrl
                matched_diffs.append(diff)
                ssw_event_stats.append({
                    "ssw_date": str(sd.date()),
                    "ssw_mean": float(ssw_mean),
                    "control_mean": float(mean_ctrl),
                    "diff": float(diff),
                    "n_control_winters": len(control_means),
                })

        if matched_diffs:
            diffs = np.array(matched_diffs)
            t, p = stats.ttest_1samp(diffs, 0)
            # Also sign test
            n_neg = np.sum(diffs < 0)
            n_total = len(diffs)
            sign_p = stats.binom_test(n_neg, n_total, 0.5) if hasattr(stats, 'binom_test') else \
                     2 * min(stats.binom.cdf(n_neg, n_total, 0.5),
                             1 - stats.binom.cdf(n_neg - 1, n_total, 0.5))

            section[outcome_name] = {
                "n_ssw_events": len(matched_diffs),
                "mean_diff": float(diffs.mean()),
                "std_diff": float(diffs.std()),
                "t_stat": float(t),
                "t_p_value": float(p),
                "n_negative": int(n_neg),
                "n_total": n_total,
                "sign_test_p": float(sign_p),
                "pct_decrease": float(100 * n_neg / n_total),
                "events": ssw_event_stats,
                "method": "Matched SSW vs same day-of-season in non-SSW winters",
            }
            LOG.info("  SSW matched %s: mean_diff=%.3f t=%.2f p=%.4f (%d/%d negative, sign_p=%.4f)",
                     outcome_name, diffs.mean(), t, p, n_neg, n_total, sign_p)

    results["part8_ssw_matched"] = section
    save_results(results)
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
        ("part1_case_crossover", part1_case_crossover),
        ("part2_anomaly", part2_anomaly_analysis),
        ("part3_winter_fe", part3_winter_fixed_effects),
        ("part4_block_bootstrap", part4_block_bootstrap),
        ("part5_spec_curve", part5_specification_curve),
        ("part6_falsification", part6_falsification),
        ("part7_loocv", part7_prediction_loocv),
        ("part8_ssw_matched", part8_ssw_matched),
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
            save_results(results)

    # Summary
    print("\n" + "=" * 70)
    print("ROBUSTNESS RESOLUTION COMPLETE")
    print("=" * 70)
    for name, _ in parts:
        if name in results:
            if isinstance(results[name], dict) and "error" in results[name] and \
               isinstance(results[name].get("error"), str) and "traceback" in results[name]:
                print("  %s: FAILED" % name)
            else:
                print("  %s: OK" % name)
        else:
            print("  %s: MISSING" % name)


if __name__ == "__main__":
    main()
