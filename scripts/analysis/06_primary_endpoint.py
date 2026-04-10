"""
06_primary_endpoint.py — Primary Endpoint Test: Avalanche Activity
===================================================================
Tests the central hypothesis: are geomagnetic disturbance events
followed by elevated Swiss natural avalanche activity?

Model: Negative Binomial regression for daily count data.
"""
import sys
from pathlib import Path
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import NegativeBinomial
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import NegativeBinomial as NB_family
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, LOG, load_panel

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def run_primary_model(panel: pd.DataFrame, outcome_var: str = "aai_all_natural",
                      exposure_var: str = "post_event_5_21d",
                      label: str = "primary") -> dict:
    """
    Run negative binomial / Poisson regression for daily avalanche counts.
    Returns rate ratio, CI, and model diagnostics.
    """
    df = panel.copy()

    # Ensure outcome is non-negative integer-like
    df[outcome_var] = df[outcome_var].fillna(0)
    # Round for count model (AAI can be fractional due to weighting)
    df["y"] = np.round(df[outcome_var]).astype(int).clip(lower=0)

    # Build design matrix
    covariates = [exposure_var]

    # Add confounders
    confound_cols = []
    for c in ["day_of_season", "day_of_season_sq", "nao_daily", "qbo_u50",
              "ncep_z500_nh", "ncep_slp_nh"]:
        if c in df.columns and df[c].notna().sum() > len(df) * 0.5:
            confound_cols.append(c)

    covariates.extend(confound_cols)

    # Winter fixed effects (dummies)
    winters = df["winter_id"].dropna().unique()
    if len(winters) > 1:
        winter_dummies = pd.get_dummies(df["winter_id"], prefix="w", drop_first=True, dtype=float)
        for col in winter_dummies.columns:
            df[col] = winter_dummies[col].astype(float)
            covariates.append(col)

    # Month fixed effects
    month_dummies = pd.get_dummies(df["month"], prefix="m", drop_first=True, dtype=float)
    for col in month_dummies.columns:
        df[col] = month_dummies[col].astype(float)
        covariates.append(col)

    # Drop NaN
    all_cols = ["y"] + covariates
    df_clean = df[all_cols].dropna()

    if len(df_clean) < 50:
        return {"error": f"Too few observations: {len(df_clean)}"}

    Y = df_clean["y"].astype(float)
    X = df_clean[covariates].astype(float)
    X = sm.add_constant(X)

    # Fit negative binomial model
    try:
        model = sm.GLM(Y, X, family=sm.families.NegativeBinomial(alpha=1.0))
        result = model.fit(maxiter=200)
    except Exception:
        LOG.info("NB GLM failed, trying Poisson...")
        model = sm.GLM(Y, X, family=sm.families.Poisson())
        result = model.fit(maxiter=200)

    # Extract exposure coefficient
    beta = float(result.params[exposure_var])
    se = float(result.bse[exposure_var])
    pval = float(result.pvalues[exposure_var])
    ci = result.conf_int(alpha=0.05).loc[exposure_var]

    rate_ratio = np.exp(beta)
    rr_ci_lo = np.exp(float(ci.iloc[0]))
    rr_ci_hi = np.exp(float(ci.iloc[1]))

    # Summary stats
    n_exposed = int(df_clean[exposure_var].sum())
    n_unexposed = int(len(df_clean) - n_exposed)
    mean_exposed = float(Y[df_clean[exposure_var] == 1].mean())
    mean_unexposed = float(Y[df_clean[exposure_var] == 0].mean())
    raw_ratio = mean_exposed / mean_unexposed if mean_unexposed > 0 else np.nan

    res = {
        "label": label,
        "outcome": outcome_var,
        "exposure": exposure_var,
        "n_obs": int(result.nobs),
        "n_exposed_days": n_exposed,
        "n_unexposed_days": n_unexposed,
        "mean_count_exposed": mean_exposed,
        "mean_count_unexposed": mean_unexposed,
        "raw_ratio": float(raw_ratio),
        "beta": beta,
        "se": se,
        "p_value": pval,
        "rate_ratio": float(rate_ratio),
        "rr_ci_lower": float(rr_ci_lo),
        "rr_ci_upper": float(rr_ci_hi),
        "aic": float(result.aic),
        "bic": float(result.bic_llf) if hasattr(result, "bic_llf") else None,
        "deviance": float(result.deviance),
        "n_confounders": len(confound_cols),
        "n_winter_fe": len(winters) - 1 if len(winters) > 1 else 0,
    }
    return res


def main():
    panel = load_panel(winter_only=True)

    # Filter to period with SLF data
    panel = panel[panel["aai_all_natural"].notna()].copy()

    LOG.info("Analysis panel: %d winter days with SLF data", len(panel))

    all_results = {}

    # ─── PRIMARY MODEL: Stratospheric pathway (5-21d) ─────────────────
    LOG.info("Running primary model: stratospheric pathway (5-21d)...")
    res = run_primary_model(panel, "aai_all_natural", "post_event_5_21d", "primary_strat")
    all_results["primary_strat_5_21d"] = res
    LOG.info("  Rate ratio: %.3f [%.3f, %.3f], p=%.4f",
             res["rate_ratio"], res["rr_ci_lower"], res["rr_ci_upper"], res["p_value"])

    # ─── FAST PATHWAY: Direct effect (1-3d) ───────────────────────────
    LOG.info("Running fast pathway model (1-3d)...")
    res = run_primary_model(panel, "aai_all_natural", "post_event_1_3d", "fast_1_3d")
    all_results["fast_pathway_1_3d"] = res
    LOG.info("  Rate ratio: %.3f [%.3f, %.3f], p=%.4f",
             res["rate_ratio"], res["rr_ci_lower"], res["rr_ci_upper"], res["p_value"])

    # ─── SECONDARY OUTCOMES ───────────────────────────────────────────
    secondary_outcomes = [
        ("dry_natural_size_1234", "Dry natural avalanches"),
        ("wet_natural_size_1234", "Wet natural avalanches"),
        ("natural_size_234", "Natural avalanches size ≥ 2"),
        ("natural_size_1234", "Natural avalanches all sizes"),
        ("max_size", "Max avalanche size"),
    ]

    for var, desc in secondary_outcomes:
        if var in panel.columns:
            LOG.info("Secondary: %s (%s)", var, desc)
            res = run_primary_model(panel, var, "post_event_5_21d", f"secondary_{var}")
            all_results[f"secondary_{var}"] = res
            LOG.info("  RR=%.3f, p=%.4f", res["rate_ratio"], res["p_value"])

    # ─── DOSE-RESPONSE: Split events by intensity ─────────────────────
    # Create strong event indicator (Kp >= 6 or Dst <= -100)
    panel_dose = panel.copy()
    panel_dose["strong_event"] = 0
    event_days = panel_dose[panel_dose["geo_event"] == 1]
    strong_days = event_days[(event_days["kp_max"] >= 6) | (event_days["dst_min"] <= -100)]
    for ed in strong_days.index:
        m = (panel_dose.index >= ed + pd.Timedelta(days=5)) & \
            (panel_dose.index <= ed + pd.Timedelta(days=21))
        panel_dose.loc[m, "strong_event"] = 1

    LOG.info("Dose-response: strong events...")
    res = run_primary_model(panel_dose, "aai_all_natural", "strong_event", "dose_strong")
    all_results["dose_strong"] = res
    LOG.info("  Strong event RR=%.3f, p=%.4f", res["rate_ratio"], res["p_value"])

    # ─── EXTENDED PERIOD (no MLS requirement) ─────────────────────────
    panel_full = load_panel(winter_only=True)
    panel_full = panel_full[panel_full["aai_all_natural"].notna()].copy()
    LOG.info("Extended period: %d days", len(panel_full))
    res = run_primary_model(panel_full, "aai_all_natural", "post_event_5_21d", "extended")
    all_results["extended_full_period"] = res
    LOG.info("  Extended RR=%.3f, p=%.4f", res["rate_ratio"], res["p_value"])

    # Save results
    out = RESULTS / "primary_endpoint.json"
    out.write_text(json.dumps(all_results, indent=2, default=str), encoding="utf-8")

    # Summary table
    print(f"\n{'='*80}")
    print("PRIMARY ENDPOINT RESULTS")
    print(f"{'='*80}")
    print(f"\n{'Model':<30} {'RR':>6} {'95% CI':>16} {'p-value':>10} {'N':>6}")
    print("-" * 75)
    for key, r in all_results.items():
        if "error" in r:
            print(f"{key:<30} ERROR: {r['error']}")
            continue
        sig = "***" if r["p_value"] < 0.001 else "**" if r["p_value"] < 0.01 else "*" if r["p_value"] < 0.05 else "†" if r["p_value"] < 0.1 else ""
        print(f"{key:<30} {r['rate_ratio']:>6.3f} [{r['rr_ci_lower']:.3f}, {r['rr_ci_upper']:.3f}]"
              f" {r['p_value']:>10.4f} {r['n_obs']:>6} {sig}")

    # Key interpretation
    pr = all_results.get("primary_strat_5_21d", {})
    if "rate_ratio" in pr:
        pct_change = (pr["rate_ratio"] - 1) * 100
        print(f"\n{'='*80}")
        print(f"PRIMARY RESULT: {pct_change:+.1f}% change in natural avalanche activity")
        print(f"  Rate ratio: {pr['rate_ratio']:.3f} [{pr['rr_ci_lower']:.3f}, {pr['rr_ci_upper']:.3f}]")
        print(f"  Raw mean exposed: {pr['mean_count_exposed']:.2f} vs unexposed: {pr['mean_count_unexposed']:.2f}")
        print(f"  p = {pr['p_value']:.4f}")


if __name__ == "__main__":
    main()
