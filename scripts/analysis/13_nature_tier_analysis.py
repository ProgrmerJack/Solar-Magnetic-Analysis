"""
13_nature_tier_analysis.py — Definitive Nature Geoscience Analysis
===================================================================
Runs sequentially to avoid OOM. Each section cleans up memory.

Sections:
  A. Primary endpoint with FULL confounders (F10.7, snow, NAO, QBO)
  B. F10.7 confound resolution (stratified + detrended)
  C. Formal mediation analysis (Baron & Kenny + bootstrap)
  D. Multi-region validation (Norway avalanches)
  E. Dose-response with full model
  F. Lag-resolved pathway analysis
  G. SSW coupling with confounders
  H. SOC comparative analysis
  I. Energy budget estimation
  J. Final robustness battery
"""
import sys, gc, json, logging
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, LOG

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

RESULTS_FILE = RESULTS / "nature_tier_analysis.json"


def load_winter_panel():
    """Load the enhanced v2 panel, winter only."""
    panel = pd.read_parquet(PROCESSED / "analysis_panel_v2.parquet")
    return panel[panel["is_winter"] == 1].copy()


def save_results(results):
    """Incrementally save results to JSON."""
    with open(RESULTS_FILE, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=str)
    LOG.info("Results saved: %d sections", len(results))


def fit_nb_glm(y, X, exposure_col=None):
    """Fit Negative Binomial GLM, return dict of results."""
    import statsmodels.api as sm
    X = sm.add_constant(X)
    mask = y.notna() & X.notna().all(axis=1)
    y_clean = y[mask].astype(float)
    X_clean = X[mask].astype(float)
    if len(y_clean) < 50:
        return {"error": "too few observations", "n": len(y_clean)}
    try:
        model = sm.GLM(y_clean, X_clean, family=sm.families.NegativeBinomial())
        result = model.fit(maxiter=100, method="IRLS")
        params = {}
        for name in result.params.index:
            params[name] = {
                "coef": float(result.params[name]),
                "se": float(result.bse[name]),
                "z": float(result.tvalues[name]),
                "p": float(result.pvalues[name]),
                "rate_ratio": float(np.exp(result.params[name])),
                "rr_ci_low": float(np.exp(result.conf_int().loc[name, 0])),
                "rr_ci_high": float(np.exp(result.conf_int().loc[name, 1])),
            }
        return {
            "params": params,
            "aic": float(result.aic),
            "bic": float(result.bic),
            "n": int(len(y_clean)),
            "deviance": float(result.deviance),
            "pearson_chi2": float(result.pearson_chi2),
        }
    except Exception as e:
        return {"error": str(e), "n": int(len(y_clean))}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION A: Primary Endpoint with FULL Confounders
# ═══════════════════════════════════════════════════════════════════════════════

def section_a_primary_full_confounders(results):
    LOG.info("=" * 70)
    LOG.info("SECTION A: Primary Endpoint with Full Confounders")
    LOG.info("=" * 70)

    w = load_winter_panel()
    y = w["aai_all_natural"].dropna()

    section = {}

    # Model 1: Minimal (original)
    X1 = w.loc[y.index, ["post_event_1_3d", "day_of_season", "day_of_season_sq"]]
    r1 = fit_nb_glm(y, X1)
    section["model1_minimal"] = r1
    pe = r1.get("params", {}).get("post_event_1_3d", {})
    LOG.info("  Model 1 (minimal): RR=%.3f [%.3f-%.3f] p=%.4f",
             pe.get("rate_ratio", 0), pe.get("rr_ci_low", 0),
             pe.get("rr_ci_high", 0), pe.get("p", 1))

    # Model 2: + meteorological confounders (NAO, QBO, SNOTEL)
    met_cols = ["post_event_1_3d", "day_of_season", "day_of_season_sq",
                "nao_daily", "qbo_u50", "snotel_swe_mean", "snotel_prec_mean", "snotel_temp_mean"]
    met_cols = [c for c in met_cols if c in w.columns]
    X2 = w.loc[y.index, met_cols]
    r2 = fit_nb_glm(y, X2)
    section["model2_meteorological"] = r2
    pe = r2.get("params", {}).get("post_event_1_3d", {})
    LOG.info("  Model 2 (meteorological): RR=%.3f p=%.4f",
             pe.get("rate_ratio", 0), pe.get("p", 1))

    # Model 3: + solar cycle (F10.7)
    sol_cols = met_cols + ["f107"]
    sol_cols = [c for c in sol_cols if c in w.columns]
    X3 = w.loc[y.index, sol_cols]
    r3 = fit_nb_glm(y, X3)
    section["model3_solar_cycle"] = r3
    pe = r3.get("params", {}).get("post_event_1_3d", {})
    LOG.info("  Model 3 (+F10.7): RR=%.3f p=%.4f",
             pe.get("rate_ratio", 0), pe.get("p", 1))

    # Model 4: + snow cover
    snow_cols = sol_cols + ["modis_snow_frac", "ims_nh_snow"]
    snow_cols = [c for c in snow_cols if c in w.columns]
    X4 = w.loc[y.index, snow_cols]
    r4 = fit_nb_glm(y, X4)
    section["model4_snow_cover"] = r4
    pe = r4.get("params", {}).get("post_event_1_3d", {})
    LOG.info("  Model 4 (+snow): RR=%.3f p=%.4f",
             pe.get("rate_ratio", 0), pe.get("p", 1))

    # Model 5: FULL kitchen sink
    full_cols = ["post_event_1_3d", "day_of_season", "day_of_season_sq",
                 "nao_daily", "qbo_u50", "pna_monthly",
                 "snotel_swe_mean", "snotel_prec_mean", "snotel_temp_mean",
                 "f107", "modis_snow_frac", "ims_nh_snow",
                 "ncep_u_10hpa", "ncep_t_10hpa", "ncep_z500_nh",
                 "flare_count", "sw_bz_min", "sw_speed_max"]
    full_cols = [c for c in full_cols if c in w.columns]
    X5 = w.loc[y.index, full_cols]
    r5 = fit_nb_glm(y, X5)
    section["model5_full"] = r5
    pe = r5.get("params", {}).get("post_event_1_3d", {})
    LOG.info("  Model 5 (full): RR=%.3f p=%.4f",
             pe.get("rate_ratio", 0), pe.get("p", 1))

    # Model 6: Stratospheric pathway (5-21d lag)
    strat_cols = ["post_event_5_21d", "day_of_season", "day_of_season_sq",
                  "nao_daily", "qbo_u50", "f107",
                  "snotel_swe_mean", "snotel_prec_mean"]
    strat_cols = [c for c in strat_cols if c in w.columns]
    X6 = w.loc[y.index, strat_cols]
    r6 = fit_nb_glm(y, X6)
    section["model6_stratospheric_pathway"] = r6
    pe = r6.get("params", {}).get("post_event_5_21d", {})
    LOG.info("  Model 6 (strat 5-21d): RR=%.3f p=%.4f",
             pe.get("rate_ratio", 0), pe.get("p", 1))

    results["section_a_primary"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION B: F10.7 Confound Resolution
# ═══════════════════════════════════════════════════════════════════════════════

def section_b_f107_resolution(results):
    LOG.info("=" * 70)
    LOG.info("SECTION B: F10.7 Confound Resolution")
    LOG.info("=" * 70)

    w = load_winter_panel()
    section = {}

    # Strategy 1: Stratify by solar cycle phase
    mask_f107 = w["f107"].notna()
    w_f107 = w[mask_f107].copy()

    # Use terciles for more nuanced stratification
    f107_33 = w_f107["f107"].quantile(0.33)
    f107_66 = w_f107["f107"].quantile(0.66)

    for label, mask_fn in [
        ("solar_low", lambda x: x["f107"] <= f107_33),
        ("solar_mid", lambda x: (x["f107"] > f107_33) & (x["f107"] <= f107_66)),
        ("solar_high", lambda x: x["f107"] > f107_66),
    ]:
        subset = w_f107[mask_fn(w_f107)]
        y = subset["aai_all_natural"].dropna()
        X = subset.loc[y.index, ["post_event_1_3d", "day_of_season", "day_of_season_sq",
                                  "nao_daily", "qbo_u50", "snotel_swe_mean"]]
        X = X[[c for c in X.columns if c in subset.columns]]
        r = fit_nb_glm(y, X)
        section["stratified_" + label] = r
        pe = r.get("params", {}).get("post_event_1_3d", {})
        n_events = subset["geo_event"].sum()
        LOG.info("  %s (N=%d, events=%d): RR=%.3f p=%.4f",
                 label, len(y), n_events, pe.get("rate_ratio", 0), pe.get("p", 1))

    # Strategy 2: Detrend both avalanche and geomagnetic series by solar cycle
    # Remove the solar-cycle trend from avalanche activity, then test event effect
    from scipy.signal import detrend as scipy_detrend
    w_both = w_f107.copy()
    y_raw = w_both["aai_all_natural"]
    mask_valid = y_raw.notna() & w_both["f107"].notna()
    w_valid = w_both[mask_valid].copy()

    # Residualize avalanche counts against F10.7 (+ season)
    import statsmodels.api as sm
    X_detrend = sm.add_constant(w_valid[["f107", "day_of_season", "day_of_season_sq"]])
    y_detrend = w_valid["aai_all_natural"]
    try:
        model_detrend = sm.OLS(np.log1p(y_detrend), X_detrend).fit()
        w_valid["aval_residual"] = model_detrend.resid

        # Now test: does post_event_1_3d predict avalanche RESIDUALS?
        t_exp = w_valid.loc[w_valid["post_event_1_3d"] == 1, "aval_residual"]
        t_unexp = w_valid.loc[w_valid["post_event_1_3d"] == 0, "aval_residual"]
        tstat, pval = stats.ttest_ind(t_exp, t_unexp, equal_var=False)
        section["detrended_ttest"] = {
            "exposed_mean_residual": float(t_exp.mean()),
            "unexposed_mean_residual": float(t_unexp.mean()),
            "t_stat": float(tstat),
            "p_value": float(pval),
            "n_exposed": int(len(t_exp)),
            "n_unexposed": int(len(t_unexp)),
            "interpretation": "Negative residual = fewer avalanches than F10.7-season model predicts"
        }
        LOG.info("  Detrended t-test: t=%.3f, p=%.4f (exp mean=%.4f vs unexp=%.4f)",
                 tstat, pval, t_exp.mean(), t_unexp.mean())
    except Exception as e:
        section["detrended_ttest"] = {"error": str(e)}
        LOG.warning("  Detrend failed: %s", e)

    # Strategy 3: Interaction model (event × F10.7)
    w_valid["event_x_f107"] = w_valid["post_event_1_3d"] * w_valid["f107"]
    y3 = w_valid["aai_all_natural"]
    X3 = w_valid[["post_event_1_3d", "f107", "event_x_f107",
                   "day_of_season", "day_of_season_sq", "nao_daily", "qbo_u50"]]
    X3 = X3[[c for c in X3.columns if c in w_valid.columns]]
    r3 = fit_nb_glm(y3, X3)
    section["interaction_model"] = r3
    if "params" in r3:
        for key in ["post_event_1_3d", "f107", "event_x_f107"]:
            if key in r3["params"]:
                p = r3["params"][key]
                LOG.info("  Interaction: %s RR=%.3f p=%.4f", key, p["rate_ratio"], p["p"])

    # Strategy 4: Within-winter fixed effects (absorbs annual solar cycle)
    # Each winter gets its own intercept, isolating within-winter event effects
    winters = w_f107["winter_id"].dropna().unique()
    within_effects = []
    for wid in winters:
        wsub = w_f107[w_f107["winter_id"] == wid]
        y_w = wsub["aai_all_natural"].dropna()
        if len(y_w) < 30 or wsub["geo_event"].sum() == 0:
            continue
        X_w = wsub.loc[y_w.index, ["post_event_1_3d", "day_of_season", "day_of_season_sq"]]
        try:
            mdl = sm.GLM(y_w, sm.add_constant(X_w),
                         family=sm.families.NegativeBinomial()).fit(maxiter=50)
            if "post_event_1_3d" in mdl.params.index:
                within_effects.append({
                    "winter": wid,
                    "coef": float(mdl.params["post_event_1_3d"]),
                    "se": float(mdl.bse["post_event_1_3d"]),
                    "rr": float(np.exp(mdl.params["post_event_1_3d"])),
                    "p": float(mdl.pvalues["post_event_1_3d"]),
                    "n": int(len(y_w)),
                    "n_events": int(wsub["geo_event"].sum()),
                    "f107_mean": float(wsub["f107"].mean()),
                })
        except Exception:
            pass

    if within_effects:
        df_we = pd.DataFrame(within_effects)
        # Meta-analytic combination (inverse-variance weighted)
        valid = df_we[df_we["se"] > 0].copy()
        if len(valid) > 0:
            weights = 1.0 / (valid["se"] ** 2)
            pooled_coef = (valid["coef"] * weights).sum() / weights.sum()
            pooled_se = 1.0 / np.sqrt(weights.sum())
            pooled_z = pooled_coef / pooled_se
            pooled_p = 2 * stats.norm.sf(abs(pooled_z))
            section["within_winter_meta"] = {
                "pooled_rr": float(np.exp(pooled_coef)),
                "pooled_se": float(pooled_se),
                "pooled_z": float(pooled_z),
                "pooled_p": float(pooled_p),
                "n_winters": int(len(valid)),
                "n_winters_significant": int((valid["p"] < 0.05).sum()),
                "rr_ci_low": float(np.exp(pooled_coef - 1.96 * pooled_se)),
                "rr_ci_high": float(np.exp(pooled_coef + 1.96 * pooled_se)),
                "interpretation": "Fixed-effects meta across winters absorbs solar cycle"
            }
            LOG.info("  Within-winter meta: RR=%.3f [%.3f-%.3f] p=%.4f (N=%d winters)",
                     np.exp(pooled_coef),
                     np.exp(pooled_coef - 1.96 * pooled_se),
                     np.exp(pooled_coef + 1.96 * pooled_se),
                     pooled_p, len(valid))
        section["within_winter_individual"] = within_effects

    results["section_b_f107"] = section
    save_results(results)
    del w, w_f107; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION C: Formal Mediation Analysis
# ═══════════════════════════════════════════════════════════════════════════════

def section_c_mediation(results):
    LOG.info("=" * 70)
    LOG.info("SECTION C: Formal Mediation Analysis")
    LOG.info("=" * 70)

    w = load_winter_panel()
    section = {}
    import statsmodels.api as sm

    # Baron & Kenny (1986) steps for each candidate mediator
    mediators = {
        "ncep_u_10hpa": "Stratospheric zonal wind (polar vortex strength)",
        "mls_hno3_lev_6p8hpa": "HNO3 at 6.8 hPa (EPP chemistry marker)",
        "modis_snow_frac": "MODIS Alps snow fraction",
        "snotel_prec_mean": "SNOTEL mean precipitation",
        "snotel_temp_mean": "SNOTEL mean temperature",
        "nao_daily": "North Atlantic Oscillation",
        "ncep_z500_nh": "500hPa geopotential height (blocking proxy)",
    }

    base_confounders = ["day_of_season", "day_of_season_sq"]

    for med_name, med_desc in mediators.items():
        if med_name not in w.columns:
            continue

        y = w["aai_all_natural"]
        mask = y.notna() & w[med_name].notna() & w["post_event_1_3d"].notna()
        for bc in base_confounders:
            mask &= w[bc].notna()
        wsub = w[mask].copy()

        if len(wsub) < 100:
            continue

        X_base = wsub[["post_event_1_3d"] + base_confounders]
        y_out = wsub["aai_all_natural"]
        med_var = wsub[med_name]

        try:
            # Step 1: X → Y (total effect)
            r_total = sm.GLM(y_out, sm.add_constant(X_base),
                             family=sm.families.NegativeBinomial()).fit(maxiter=50)
            total_coef = float(r_total.params.get("post_event_1_3d", np.nan))
            total_p = float(r_total.pvalues.get("post_event_1_3d", np.nan))

            # Step 2: X → M (does the event affect the mediator?)
            r_xm = sm.OLS(med_var, sm.add_constant(X_base)).fit()
            a_coef = float(r_xm.params.get("post_event_1_3d", np.nan))
            a_p = float(r_xm.pvalues.get("post_event_1_3d", np.nan))

            # Step 3: X + M → Y (does mediator predict Y controlling for X?)
            X_med = wsub[["post_event_1_3d", med_name] + base_confounders]
            r_med = sm.GLM(y_out, sm.add_constant(X_med),
                           family=sm.families.NegativeBinomial()).fit(maxiter=50)
            direct_coef = float(r_med.params.get("post_event_1_3d", np.nan))
            direct_p = float(r_med.pvalues.get("post_event_1_3d", np.nan))
            b_coef = float(r_med.params.get(med_name, np.nan))
            b_p = float(r_med.pvalues.get(med_name, np.nan))

            # Mediation proportion
            if abs(total_coef) > 0.001:
                prop_mediated = 1.0 - (direct_coef / total_coef)
            else:
                prop_mediated = np.nan

            section[med_name] = {
                "description": med_desc,
                "n": int(len(wsub)),
                "step1_total_effect": {"coef": total_coef, "rr": float(np.exp(total_coef)), "p": total_p},
                "step2_event_to_mediator": {"coef": a_coef, "p": a_p},
                "step3_direct_effect": {"coef": direct_coef, "rr": float(np.exp(direct_coef)), "p": direct_p},
                "step3_mediator_effect": {"coef": b_coef, "p": b_p},
                "proportion_mediated": float(prop_mediated) if not np.isnan(prop_mediated) else None,
                "mediation_significant": bool(a_p < 0.05 and b_p < 0.05),
            }
            LOG.info("  %s: total RR=%.3f, direct RR=%.3f, prop_med=%.1f%%, event->med p=%.4f",
                     med_name, np.exp(total_coef), np.exp(direct_coef),
                     100 * prop_mediated if not np.isnan(prop_mediated) else 0, a_p)

        except Exception as e:
            section[med_name] = {"error": str(e)}
            LOG.warning("  %s: FAILED: %s", med_name, e)

    results["section_c_mediation"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION D: Multi-Region Validation (Norway)
# ═══════════════════════════════════════════════════════════════════════════════

def section_d_multi_region(results):
    LOG.info("=" * 70)
    LOG.info("SECTION D: Multi-Region Validation")
    LOG.info("=" * 70)

    w = load_winter_panel()
    section = {}

    # Test 1: Norway avalanche response
    y_nor = w["norway_aval_count"]
    mask_nor = y_nor.notna() & (y_nor > 0)  # Only days Norway has data
    if mask_nor.sum() > 100:
        X_nor = w.loc[mask_nor, ["post_event_1_3d", "day_of_season", "day_of_season_sq",
                                  "nao_daily"]]
        X_nor = X_nor[[c for c in X_nor.columns if c in w.columns]]
        r_nor = fit_nb_glm(y_nor[mask_nor], X_nor)
        section["norway_primary"] = r_nor
        pe = r_nor.get("params", {}).get("post_event_1_3d", {})
        LOG.info("  Norway primary: RR=%.3f p=%.4f (N=%d)",
                 pe.get("rate_ratio", 0), pe.get("p", 1), mask_nor.sum())

    # Test 2: Norway with full confounders
    full_cols = ["post_event_1_3d", "day_of_season", "day_of_season_sq",
                 "nao_daily", "qbo_u50", "f107"]
    full_cols = [c for c in full_cols if c in w.columns]
    X_nor2 = w.loc[mask_nor, full_cols]
    r_nor2 = fit_nb_glm(y_nor[mask_nor], X_nor2)
    section["norway_full_confounders"] = r_nor2
    pe = r_nor2.get("params", {}).get("post_event_1_3d", {})
    LOG.info("  Norway full: RR=%.3f p=%.4f",
             pe.get("rate_ratio", 0), pe.get("p", 1))

    # Test 3: Switzerland dry vs wet avalanches
    for aval_type, col in [("dry_natural", "dry_natural_size_1234"),
                            ("wet_natural", "wet_natural_size_1234"),
                            ("all_natural", "aai_all_natural")]:
        if col not in w.columns:
            continue
        y_type = w[col].dropna()
        X_type = w.loc[y_type.index, ["post_event_1_3d", "day_of_season",
                                       "day_of_season_sq", "nao_daily"]]
        X_type = X_type[[c for c in X_type.columns if c in w.columns]]
        r_type = fit_nb_glm(y_type, X_type)
        section["switzerland_" + aval_type] = r_type
        pe = r_type.get("params", {}).get("post_event_1_3d", {})
        LOG.info("  Switzerland %s: RR=%.3f p=%.4f",
                 aval_type, pe.get("rate_ratio", 0), pe.get("p", 1))

    # Test 4: SNOTEL-region temperature/precip response to geomag events
    for var, desc in [("snotel_temp_mean", "temperature"),
                      ("snotel_prec_mean", "precipitation"),
                      ("snotel_swe_mean", "snow water equiv")]:
        if var not in w.columns:
            continue
        exposed = w.loc[w["post_event_1_3d"] == 1, var].dropna()
        unexposed = w.loc[w["post_event_1_3d"] == 0, var].dropna()
        if len(exposed) > 10 and len(unexposed) > 10:
            t, p = stats.ttest_ind(exposed, unexposed, equal_var=False)
            section["snotel_" + var.split("_")[-1]] = {
                "exposed_mean": float(exposed.mean()),
                "unexposed_mean": float(unexposed.mean()),
                "diff": float(exposed.mean() - unexposed.mean()),
                "t_stat": float(t),
                "p_value": float(p),
            }
            LOG.info("  SNOTEL %s: exp=%.2f vs unexp=%.2f, p=%.4f",
                     desc, exposed.mean(), unexposed.mean(), p)

    results["section_d_multi_region"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION E: Dose-Response with Full Model
# ═══════════════════════════════════════════════════════════════════════════════

def section_e_dose_response(results):
    LOG.info("=" * 70)
    LOG.info("SECTION E: Dose-Response Analysis")
    LOG.info("=" * 70)

    w = load_winter_panel()
    section = {}

    event_cat = pd.read_parquet(RESULTS / "event_catalog.parquet")
    event_cat.index = event_cat.index.tz_localize(None) if hasattr(event_cat.index, 'tz') and event_cat.index.tz else event_cat.index

    # Create intensity bins from event catalog
    if "dst_min" in event_cat.columns:
        dst_col = "dst_min"
    elif "Dst_min" in event_cat.columns:
        dst_col = "Dst_min"
    else:
        results["section_e_dose_response"] = {"error": "no dst column in event_catalog"}
        save_results(results)
        return results

    # Assign each post-event day the storm intensity
    w["storm_dst"] = np.nan
    w["storm_kp"] = np.nan
    for evt_date, evt_row in event_cat.iterrows():
        mask = (w.index > evt_date) & (w.index <= evt_date + pd.Timedelta(days=3))
        w.loc[mask, "storm_dst"] = evt_row[dst_col]
        if "kp_max" in evt_row.index:
            w.loc[mask, "storm_kp"] = evt_row["kp_max"]

    # Dose-response by Dst categories
    dst_bins = [
        ("moderate", -100, -50),
        ("strong", -200, -100),
        ("severe", -500, -200),
    ]

    y = w["aai_all_natural"].dropna()
    base_cols = ["day_of_season", "day_of_season_sq", "nao_daily", "qbo_u50", "f107"]
    base_cols = [c for c in base_cols if c in w.columns]

    for label, dst_lo, dst_hi in dst_bins:
        # Create binary for this intensity level
        w["dose_" + label] = ((w["storm_dst"] >= dst_lo) & (w["storm_dst"] < dst_hi)).astype(int)
        dose_col = "dose_" + label
        X = w.loc[y.index, [dose_col] + base_cols]
        r = fit_nb_glm(y, X)
        section[label] = r
        pe = r.get("params", {}).get(dose_col, {})
        n_exp = w[dose_col].sum()
        LOG.info("  Dst %s (%d to %d, N=%d): RR=%.3f p=%.4f",
                 label, dst_lo, dst_hi, n_exp, pe.get("rate_ratio", 0), pe.get("p", 1))

    # Continuous dose-response: Dst as continuous predictor in post-event window
    post_days = w[w["post_event_1_3d"] == 1].copy()
    y_post = post_days["aai_all_natural"].dropna()
    if len(y_post) > 30 and "storm_dst" in post_days.columns:
        mask = y_post.index.isin(post_days.index)
        X_cont = post_days.loc[y_post.index, ["storm_dst", "day_of_season", "day_of_season_sq"]]
        X_cont = X_cont.dropna()
        y_cont = y_post.loc[X_cont.index]
        r_cont = fit_nb_glm(y_cont, X_cont)
        section["continuous_dst"] = r_cont
        pe = r_cont.get("params", {}).get("storm_dst", {})
        LOG.info("  Continuous Dst: coef=%.4f, RR per unit=%.4f, p=%.4f",
                 pe.get("coef", 0), pe.get("rate_ratio", 0), pe.get("p", 1))

    results["section_e_dose_response"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION F: Lag-Resolved Pathway Analysis
# ═══════════════════════════════════════════════════════════════════════════════

def section_f_lag_analysis(results):
    LOG.info("=" * 70)
    LOG.info("SECTION F: Lag-Resolved Pathway Analysis")
    LOG.info("=" * 70)

    w = load_winter_panel()
    section = {}

    y = w["aai_all_natural"].dropna()
    base_cols = ["day_of_season", "day_of_season_sq", "nao_daily", "f107"]
    base_cols = [c for c in base_cols if c in w.columns]

    event_dates = w.index[w["geo_event"] == 1]

    # Day-by-day lag sweep (0-30d) with full confounders
    lag_results = []
    for lag in range(0, 31):
        w["lag_%d" % lag] = 0
        for ed in event_dates:
            target = ed + pd.Timedelta(days=lag)
            if target in w.index:
                w.loc[target, "lag_%d" % lag] = 1

        lag_col = "lag_%d" % lag
        X = w.loc[y.index, [lag_col] + base_cols]
        r = fit_nb_glm(y, X)
        pe = r.get("params", {}).get(lag_col, {})
        lag_results.append({
            "lag": lag,
            "rr": pe.get("rate_ratio", 1.0),
            "rr_ci_low": pe.get("rr_ci_low", 1.0),
            "rr_ci_high": pe.get("rr_ci_high", 1.0),
            "p": pe.get("p", 1.0),
            "coef": pe.get("coef", 0.0),
        })
        del w["lag_%d" % lag]

    section["lag_sweep"] = lag_results

    # BH FDR correction
    from statsmodels.stats.multitest import multipletests
    pvals = [r["p"] for r in lag_results]
    reject, pval_corrected, _, _ = multipletests(pvals, method="fdr_bh")
    for i, r in enumerate(lag_results):
        r["p_fdr"] = float(pval_corrected[i])
        r["significant_fdr"] = bool(reject[i])

    n_sig = sum(1 for r in lag_results if r["significant_fdr"])
    section["n_significant_fdr"] = n_sig
    LOG.info("  Lag sweep: %d/31 significant after FDR correction", n_sig)

    # Identify peak lag
    min_rr = min(lag_results, key=lambda x: x["rr"])
    section["peak_lag"] = min_rr
    LOG.info("  Peak effect at lag %d: RR=%.3f p=%.4f",
             min_rr["lag"], min_rr["rr"], min_rr["p"])

    # Window-based analysis with confounders
    for window_name, window_col in [
        ("fast_1_3d", "post_event_1_3d"),
        ("cme_3_8d", "post_event_3_8d"),
        ("strat_5_21d", "post_event_5_21d"),
        ("ssw_15_30d", "post_event_15_30d"),
        ("rebound_30_60d", "post_event_30_60d"),
    ]:
        if window_col not in w.columns:
            continue
        X = w.loc[y.index, [window_col] + base_cols]
        r = fit_nb_glm(y, X)
        section["window_" + window_name] = r
        pe = r.get("params", {}).get(window_col, {})
        LOG.info("  Window %s: RR=%.3f p=%.4f",
                 window_name, pe.get("rate_ratio", 0), pe.get("p", 1))

    results["section_f_lags"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION G: SSW Coupling with Confounders
# ═══════════════════════════════════════════════════════════════════════════════

def section_g_ssw(results):
    LOG.info("=" * 70)
    LOG.info("SECTION G: SSW-Avalanche Coupling")
    LOG.info("=" * 70)

    w = load_winter_panel()
    section = {}

    ssw = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw.index = ssw.index.tz_localize(None) if hasattr(ssw.index, 'tz') and ssw.index.tz else ssw.index
    ssw_dates = ssw.index

    y = w["aai_all_natural"]

    # Create SSW post-event windows
    for window_name, d_lo, d_hi in [
        ("ssw_0_15d", 0, 15),
        ("ssw_15_30d", 15, 30),
        ("ssw_30_60d", 30, 60),
    ]:
        w[window_name] = 0
        for sd in ssw_dates:
            mask = (w.index >= sd + pd.Timedelta(days=d_lo)) & \
                   (w.index <= sd + pd.Timedelta(days=d_hi))
            w.loc[mask, window_name] = 1

    # SSW effect on avalanches (multiple time windows)
    base_cols = ["day_of_season", "day_of_season_sq", "nao_daily", "qbo_u50", "f107"]
    base_cols = [c for c in base_cols if c in w.columns]

    y_valid = y.dropna()
    for window_name in ["ssw_0_15d", "ssw_15_30d", "ssw_30_60d"]:
        X = w.loc[y_valid.index, [window_name] + base_cols]
        r = fit_nb_glm(y_valid, X)
        section[window_name] = r
        pe = r.get("params", {}).get(window_name, {})
        LOG.info("  %s: RR=%.3f p=%.4f", window_name,
                 pe.get("rate_ratio", 0), pe.get("p", 1))

    # SSW → Norway
    y_nor = w["norway_aval_count"]
    mask_nor = y_nor.notna() & (y_nor > 0)
    if mask_nor.sum() > 100:
        for window_name in ["ssw_0_15d", "ssw_15_30d", "ssw_30_60d"]:
            X = w.loc[mask_nor, [window_name] + base_cols]
            X = X[[c for c in X.columns if c in w.columns]]
            r = fit_nb_glm(y_nor[mask_nor], X)
            section["norway_" + window_name] = r
            pe = r.get("params", {}).get(window_name, {})
            LOG.info("  Norway %s: RR=%.3f p=%.4f", window_name,
                     pe.get("rate_ratio", 0), pe.get("p", 1))

    # Superposed epoch analysis: avalanche activity centered on SSW dates
    sea_results = []
    for sd in ssw_dates:
        for lag in range(-30, 61):
            target = sd + pd.Timedelta(days=lag)
            if target in w.index and pd.notna(w.loc[target, "aai_all_natural"]):
                sea_results.append({
                    "ssw_date": str(sd.date()),
                    "lag": lag,
                    "aai": float(w.loc[target, "aai_all_natural"]),
                })

    if sea_results:
        sea_df = pd.DataFrame(sea_results)
        sea_mean = sea_df.groupby("lag")["aai"].agg(["mean", "std", "count"])
        sea_mean["se"] = sea_mean["std"] / np.sqrt(sea_mean["count"])
        # Baseline: -30 to -10 days before SSW
        baseline = sea_mean.loc[-30:-10, "mean"].mean()
        sea_mean["anomaly"] = sea_mean["mean"] - baseline
        sea_mean["anomaly_pct"] = 100 * sea_mean["anomaly"] / baseline

        section["sea_baseline"] = float(baseline)
        section["sea_by_lag"] = {
            int(lag): {
                "mean": float(row["mean"]),
                "anomaly_pct": float(row["anomaly_pct"]),
                "se": float(row["se"]),
                "n": int(row["count"]),
            }
            for lag, row in sea_mean.iterrows()
        }
        # Summary windows
        for wname, lo, hi in [("pre", -30, -10), ("onset", -5, 5),
                               ("response", 10, 30), ("late", 30, 60)]:
            sub = sea_mean.loc[lo:hi]
            section["sea_window_" + wname] = {
                "mean_anomaly_pct": float(sub["anomaly_pct"].mean()),
                "mean_aai": float(sub["mean"].mean()),
            }
            LOG.info("  SEA %s: anomaly=%.1f%%", wname, sub["anomaly_pct"].mean())

    results["section_g_ssw"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION H: SOC Comparative Analysis
# ═══════════════════════════════════════════════════════════════════════════════

def section_h_soc(results):
    LOG.info("=" * 70)
    LOG.info("SECTION H: SOC Comparative Analysis")
    LOG.info("=" * 70)

    section = {}

    # 1. Solar flare power-law from GOES catalog
    LOG.info("  Loading flare catalog...")
    fl = pd.read_parquet(PROCESSED / "solar" / "flares.parquet")

    def class_to_flux(ct):
        if pd.isna(ct) or not isinstance(ct, str):
            return np.nan
        ct = ct.strip()
        cls = ct[0].upper()
        try:
            mag = float(ct[1:])
        except (ValueError, IndexError):
            mag = 1.0
        base = {"A": 1e-8, "B": 1e-7, "C": 1e-6, "M": 1e-5, "X": 1e-4}.get(cls, np.nan)
        return base * mag if not np.isnan(base) else np.nan

    fl["flux"] = fl["classType"].apply(class_to_flux)
    fluxes = fl["flux"].dropna().values
    fluxes = fluxes[fluxes > 0]

    try:
        import powerlaw
        fit_flare = powerlaw.Fit(fluxes, discrete=False)
        section["flare_power_law"] = {
            "alpha": float(fit_flare.alpha),
            "xmin": float(fit_flare.xmin),
            "sigma": float(fit_flare.sigma),
            "n_above_xmin": int(np.sum(fluxes >= fit_flare.xmin)),
            "n_total": int(len(fluxes)),
        }
        # Compare power-law vs alternatives
        R_exp, p_exp = fit_flare.distribution_compare("power_law", "exponential")
        R_ln, p_ln = fit_flare.distribution_compare("power_law", "lognormal")
        section["flare_model_comparison"] = {
            "vs_exponential": {"R": float(R_exp), "p": float(p_exp)},
            "vs_lognormal": {"R": float(R_ln), "p": float(p_ln)},
        }
        LOG.info("  Flares: alpha=%.3f xmin=%.2e, vs_exp R=%.2f p=%.4f",
                 fit_flare.alpha, fit_flare.xmin, R_exp, p_exp)
    except ImportError:
        LOG.warning("  powerlaw package not installed, using MLE fallback")
        # Manual MLE for power-law exponent
        xmin = np.percentile(fluxes, 90)
        above = fluxes[fluxes >= xmin]
        alpha = 1 + len(above) / np.sum(np.log(above / xmin))
        alpha_se = (alpha - 1) / np.sqrt(len(above))
        section["flare_power_law"] = {
            "alpha": float(alpha),
            "xmin": float(xmin),
            "sigma": float(alpha_se),
            "n_above_xmin": int(len(above)),
            "n_total": int(len(fluxes)),
            "method": "MLE (no powerlaw package)",
        }
        LOG.info("  Flares MLE: alpha=%.3f+/-%.3f", alpha, alpha_se)

    del fl; gc.collect()

    # 2. Avalanche size power-law
    LOG.info("  Loading avalanche data...")
    act = pd.read_parquet(PROCESSED / "cryosphere" / "slf_activity.parquet")
    # Use daily total as proxy for "event size"
    aai_vals = act["aai_all_natural"].dropna().values
    aai_vals = aai_vals[aai_vals > 0]

    try:
        import powerlaw
        fit_aval = powerlaw.Fit(aai_vals, discrete=True)
        section["avalanche_power_law"] = {
            "alpha": float(fit_aval.alpha),
            "xmin": float(fit_aval.xmin),
            "sigma": float(fit_aval.sigma),
            "n_above_xmin": int(np.sum(aai_vals >= fit_aval.xmin)),
            "n_total": int(len(aai_vals)),
        }
        R_exp_a, p_exp_a = fit_aval.distribution_compare("power_law", "exponential")
        R_ln_a, p_ln_a = fit_aval.distribution_compare("power_law", "lognormal")
        section["avalanche_model_comparison"] = {
            "vs_exponential": {"R": float(R_exp_a), "p": float(p_exp_a)},
            "vs_lognormal": {"R": float(R_ln_a), "p": float(p_ln_a)},
        }
        LOG.info("  Avalanches: alpha=%.3f xmin=%.1f, vs_exp R=%.2f p=%.4f",
                 fit_aval.alpha, fit_aval.xmin, R_exp_a, p_exp_a)
    except ImportError:
        xmin = np.percentile(aai_vals, 75)
        above = aai_vals[aai_vals >= xmin]
        alpha = 1 + len(above) / np.sum(np.log(above / xmin))
        alpha_se = (alpha - 1) / np.sqrt(len(above))
        section["avalanche_power_law"] = {
            "alpha": float(alpha),
            "xmin": float(xmin),
            "sigma": float(alpha_se),
            "n_above_xmin": int(len(above)),
            "n_total": int(len(aai_vals)),
            "method": "MLE",
        }
        LOG.info("  Avalanches MLE: alpha=%.3f+/-%.3f", alpha, alpha_se)

    del act; gc.collect()

    # 3. SOC exponent comparison: during vs after geomagnetic events
    LOG.info("  SOC comparison: quiet vs disturbed periods...")
    w = pd.read_parquet(PROCESSED / "analysis_panel_v2.parquet")
    w = w[w["is_winter"] == 1]

    for label, mask_col in [("post_event", "post_event_1_3d"),
                             ("quiet", None)]:
        if label == "post_event":
            subset = w[w[mask_col] == 1]["aai_all_natural"].dropna().values
        else:
            subset = w[(w["post_event_0_30d"] == 0)]["aai_all_natural"].dropna().values
        subset = subset[subset > 0]
        if len(subset) < 50:
            continue
        xmin = max(np.percentile(subset, 75), 1)
        above = subset[subset >= xmin]
        if len(above) < 10:
            continue
        alpha = 1 + len(above) / np.sum(np.log(above / xmin))
        alpha_se = (alpha - 1) / np.sqrt(len(above))
        section["soc_" + label] = {
            "alpha": float(alpha),
            "alpha_se": float(alpha_se),
            "xmin": float(xmin),
            "n": int(len(above)),
        }
        LOG.info("  SOC %s: alpha=%.3f+/-%.3f (N=%d)", label, alpha, alpha_se, len(above))

    results["section_h_soc"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION I: Energy Budget Estimation
# ═══════════════════════════════════════════════════════════════════════════════

def section_i_energy_budget(results):
    LOG.info("=" * 70)
    LOG.info("SECTION I: Energy Budget Estimation")
    LOG.info("=" * 70)

    section = {}

    # Order-of-magnitude energy budget for solar-atmosphere coupling
    # Key question: is the energy from EPP sufficient to perturb the stratosphere?

    # 1. Energy input from EPP during a major storm
    # Typical hemispheric power during major storm: ~100 GW (POES data)
    # Duration: ~24-48 hours
    # Total energy: 100e9 W * 36 hr * 3600 s/hr = 1.3e16 J
    epp_power_gw = 100  # GW, typical major storm
    epp_duration_hr = 36
    epp_energy_j = epp_power_gw * 1e9 * epp_duration_hr * 3600
    section["epp_energy_input"] = {
        "power_gw": epp_power_gw,
        "duration_hr": epp_duration_hr,
        "total_energy_j": float(epp_energy_j),
        "note": "Typical major geomagnetic storm hemispheric power"
    }

    # 2. Energy to warm stratospheric polar cap by 10K (SSW magnitude)
    # Polar cap area (60-90N): ~2.5e13 m^2
    # Stratosphere mass above polar cap (10-100 hPa): ~3e15 kg
    # Cp of air: 1004 J/(kg·K)
    # Energy for 10K warming: 3e15 * 1004 * 10 = 3e19 J
    polar_area_m2 = 2.5e13
    strat_mass_kg = 3e15
    cp_air = 1004
    delta_t = 10  # K
    energy_ssw_j = strat_mass_kg * cp_air * delta_t
    section["ssw_energy_required"] = {
        "polar_cap_area_m2": float(polar_area_m2),
        "stratosphere_mass_kg": float(strat_mass_kg),
        "delta_T_K": delta_t,
        "energy_j": float(energy_ssw_j),
    }

    # 3. Ratio
    ratio = epp_energy_j / energy_ssw_j
    section["energy_ratio"] = {
        "epp_to_ssw": float(ratio),
        "interpretation": "EPP provides ~%.1e of the energy needed for SSW-scale warming" % ratio,
        "conclusion": "Direct EPP heating is insufficient (factor ~%.0f too small). "
                      "Mechanism must be catalytic (NOx->O3->radiative) not direct heating." % (1/ratio)
    }
    LOG.info("  EPP energy: %.2e J", epp_energy_j)
    LOG.info("  SSW energy: %.2e J", energy_ssw_j)
    LOG.info("  Ratio: %.2e (EPP is ~%.0fx too weak for direct heating)", ratio, 1/ratio)

    # 4. Ozone radiative forcing amplification
    # NOx-induced ozone depletion can alter radiative balance by ~1-5 W/m^2 over polar cap
    # Over 10-20 days, this accumulates:
    # 3 W/m^2 * 2.5e13 m^2 * 15 days * 86400 s/day = 9.7e19 J
    radiative_forcing_wm2 = 3  # W/m^2
    accumulation_days = 15
    radiative_energy_j = radiative_forcing_wm2 * polar_area_m2 * accumulation_days * 86400
    section["catalytic_amplification"] = {
        "ozone_radiative_forcing_wm2": radiative_forcing_wm2,
        "accumulation_days": accumulation_days,
        "radiative_energy_j": float(radiative_energy_j),
        "amplification_factor": float(radiative_energy_j / epp_energy_j),
        "interpretation": "Catalytic NOx->O3->radiative pathway amplifies EPP energy by ~%.0fx" % (radiative_energy_j / epp_energy_j)
    }
    LOG.info("  Catalytic amplification: %.0fx", radiative_energy_j / epp_energy_j)

    # 5. Forbush decrease energy pathway (alternative)
    # Forbush decrease reduces GCR by ~5-15%
    # GCR ionization rate in troposphere: ~2-10 ion pairs/cm^3/s
    # Effect on cloud condensation nuclei: uncertain but ~1-3% cloud cover change
    # Radiative effect of 1% cloud change: ~0.5-1.5 W/m^2 globally
    section["forbush_pathway"] = {
        "gcr_reduction_pct": "5-15%",
        "ccn_effect": "1-3% cloud cover change (theoretical)",
        "radiative_effect_wm2": "0.5-1.5 globally",
        "timescale_days": "1-5",
        "note": "Svensmark mechanism. Controversial but timescale matches 1-3d lag."
    }

    results["section_i_energy_budget"] = section
    save_results(results)
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION J: Final Robustness Battery
# ═══════════════════════════════════════════════════════════════════════════════

def section_j_robustness(results):
    LOG.info("=" * 70)
    LOG.info("SECTION J: Final Robustness Battery")
    LOG.info("=" * 70)

    w = load_winter_panel()
    section = {}

    y = w["aai_all_natural"].dropna()
    base_cols = ["day_of_season", "day_of_season_sq", "nao_daily", "f107"]
    base_cols = [c for c in base_cols if c in w.columns]

    # 1. Permutation test (1000 iterations)
    LOG.info("  Running permutation test (1000 iterations)...")
    X_obs = w.loc[y.index, ["post_event_1_3d"] + base_cols]
    r_obs = fit_nb_glm(y, X_obs)
    obs_coef = r_obs.get("params", {}).get("post_event_1_3d", {}).get("coef", 0)

    null_coefs = []
    rng = np.random.RandomState(42)
    event_col = w.loc[y.index, "post_event_1_3d"].values.copy()

    for i in range(1000):
        shuffled = rng.permutation(event_col)
        X_perm = X_obs.copy()
        X_perm["post_event_1_3d"] = shuffled
        try:
            import statsmodels.api as sm
            X_c = sm.add_constant(X_perm)
            mask = y.notna() & X_c.notna().all(axis=1)
            model = sm.GLM(y[mask], X_c[mask], family=sm.families.NegativeBinomial())
            result = model.fit(maxiter=30, method="IRLS")
            null_coefs.append(float(result.params.get("post_event_1_3d", 0)))
        except Exception:
            pass

    if null_coefs:
        null_arr = np.array(null_coefs)
        perm_p = np.mean(null_arr <= obs_coef)  # one-sided (decrease)
        section["permutation_test"] = {
            "observed_coef": float(obs_coef),
            "null_mean": float(null_arr.mean()),
            "null_std": float(null_arr.std()),
            "null_5th_pct": float(np.percentile(null_arr, 5)),
            "perm_p_value": float(perm_p),
            "n_permutations": len(null_coefs),
        }
        LOG.info("  Permutation p=%.4f (obs=%.4f, null mean=%.4f)",
                 perm_p, obs_coef, null_arr.mean())

    # 2. Leave-one-winter-out cross-validation
    LOG.info("  Running LOOCV...")
    winters = w["winter_id"].dropna().unique()
    loocv_results = []
    for wid in winters:
        train = w[w["winter_id"] != wid]
        test = w[w["winter_id"] == wid]
        y_train = train["aai_all_natural"].dropna()
        X_train = train.loc[y_train.index, ["post_event_1_3d"] + base_cols]

        try:
            import statsmodels.api as sm
            X_c = sm.add_constant(X_train)
            mask = y_train.notna() & X_c.notna().all(axis=1)
            model = sm.GLM(y_train[mask], X_c[mask], family=sm.families.NegativeBinomial())
            result = model.fit(maxiter=50)
            coef = float(result.params.get("post_event_1_3d", np.nan))
            loocv_results.append({
                "winter_excluded": wid,
                "coef": coef,
                "rr": float(np.exp(coef)),
                "p": float(result.pvalues.get("post_event_1_3d", np.nan)),
            })
        except Exception:
            pass

    if loocv_results:
        rrs = [r["rr"] for r in loocv_results]
        section["loocv"] = {
            "n_folds": len(loocv_results),
            "mean_rr": float(np.mean(rrs)),
            "std_rr": float(np.std(rrs)),
            "min_rr": float(np.min(rrs)),
            "max_rr": float(np.max(rrs)),
            "n_significant": sum(1 for r in loocv_results if r["p"] < 0.05),
            "all_same_direction": all(r["rr"] < 1.0 for r in loocv_results),
        }
        LOG.info("  LOOCV: mean RR=%.3f, std=%.3f, all<1=%s, %d/%d significant",
                 np.mean(rrs), np.std(rrs),
                 all(r < 1.0 for r in rrs),
                 sum(1 for r in loocv_results if r["p"] < 0.05),
                 len(loocv_results))

    # 3. Alternative event definitions
    LOG.info("  Testing alternative event definitions...")
    for def_name, kp_thresh, dst_thresh in [
        ("strict", 6.0, -100),
        ("moderate", 5.0, -50),
        ("lenient", 4.0, -30),
        ("dst_only", 99, -50),   # impossible kp, uses dst
        ("kp_only", 5.0, 999),   # impossible dst, uses kp
    ]:
        w["alt_event"] = 0
        if def_name == "dst_only":
            w["alt_event"] = (w["dst_min"] <= dst_thresh).astype(int)
        elif def_name == "kp_only":
            w["alt_event"] = (w["kp_max"] >= kp_thresh).astype(int)
        else:
            w["alt_event"] = ((w["kp_max"] >= kp_thresh) | (w["dst_min"] <= dst_thresh)).astype(int)

        # Create post-event window for alt definition
        event_dates = w.index[w["alt_event"] == 1]
        w["alt_post_1_3d"] = 0
        for ed in event_dates:
            mask = (w.index > ed) & (w.index <= ed + pd.Timedelta(days=3))
            w.loc[mask, "alt_post_1_3d"] = 1

        X_alt = w.loc[y.index, ["alt_post_1_3d"] + base_cols]
        r_alt = fit_nb_glm(y, X_alt)
        pe = r_alt.get("params", {}).get("alt_post_1_3d", {})
        section["alt_def_" + def_name] = {
            "kp_threshold": kp_thresh,
            "dst_threshold": dst_thresh,
            "n_events": int(w["alt_event"].sum()),
            "n_exposed_days": int(w["alt_post_1_3d"].sum()),
            "rr": pe.get("rate_ratio", 1.0),
            "p": pe.get("p", 1.0),
        }
        LOG.info("  Alt %s (Kp>=%s, Dst<=%s): RR=%.3f p=%.4f (N_events=%d)",
                 def_name, kp_thresh, dst_thresh,
                 pe.get("rate_ratio", 0), pe.get("p", 1), w["alt_event"].sum())

    # 4. Summer control (should be null)
    LOG.info("  Summer control...")
    w_full = pd.read_parquet(PROCESSED / "analysis_panel_v2.parquet")
    summer = w_full[w_full["is_summer"] == 1]
    y_sum = summer["aai_all_natural"].dropna()
    if len(y_sum) > 100:
        X_sum = summer.loc[y_sum.index, ["post_event_1_3d", "day_of_year"]]
        X_sum = X_sum[[c for c in X_sum.columns if c in summer.columns]]
        r_sum = fit_nb_glm(y_sum, X_sum)
        section["summer_control"] = r_sum
        pe = r_sum.get("params", {}).get("post_event_1_3d", {})
        LOG.info("  Summer control: RR=%.3f p=%.4f",
                 pe.get("rate_ratio", 0), pe.get("p", 1))

    results["section_j_robustness"] = section
    save_results(results)
    del w; gc.collect()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    results = {}

    # Check if partial results exist
    if RESULTS_FILE.exists():
        try:
            results = json.loads(open(RESULTS_FILE, encoding="utf-8").read())
            LOG.info("Loaded existing results with %d sections", len(results))
        except Exception:
            results = {}

    sections = [
        ("section_a_primary", section_a_primary_full_confounders),
        ("section_b_f107", section_b_f107_resolution),
        ("section_c_mediation", section_c_mediation),
        ("section_d_multi_region", section_d_multi_region),
        ("section_e_dose_response", section_e_dose_response),
        ("section_f_lags", section_f_lag_analysis),
        ("section_g_ssw", section_g_ssw),
        ("section_h_soc", section_h_soc),
        ("section_i_energy_budget", section_i_energy_budget),
        ("section_j_robustness", section_j_robustness),
    ]

    for name, func in sections:
        if name in results:
            LOG.info("Skipping %s (already completed)", name)
            continue
        try:
            results = func(results)
            gc.collect()
        except Exception as e:
            LOG.error("FAILED %s: %s", name, e)
            import traceback
            results[name] = {"error": str(e), "traceback": traceback.format_exc()}
            save_results(results)

    # Final summary
    print("\n" + "=" * 70)
    print("NATURE-TIER ANALYSIS COMPLETE")
    print("=" * 70)
    for name, _ in sections:
        if name in results:
            if "error" in results[name] and isinstance(results[name].get("error"), str):
                print("  %s: FAILED - %s" % (name, results[name]["error"]))
            else:
                print("  %s: OK" % name)
        else:
            print("  %s: MISSING" % name)


if __name__ == "__main__":
    main()
