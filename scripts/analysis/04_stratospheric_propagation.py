"""
04_stratospheric_propagation.py — Distributed-Lag Model for Vortex Response
=============================================================================
Tests whether geomagnetic disturbances weaken the polar vortex.
Uses NCEP 10 hPa zonal wind and temperature as response variables.
"""
import sys
from pathlib import Path
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS, WLS

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, LOG, load_panel

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def distributed_lag_model(panel: pd.DataFrame, response_var: str,
                          max_lag: int = 30) -> dict:
    """
    Fit distributed-lag model:
      Y_t = α + Σ(l=0..max_lag) β_l × Event_{t-l} + controls + ε_t
    with Newey-West HAC standard errors.
    """
    df = panel.copy()

    # Create lagged event indicators
    for lag in range(max_lag + 1):
        df[f"event_lag_{lag}"] = df["geo_event"].shift(lag).fillna(0)

    # Control variables
    controls = []
    for c in ["day_of_season", "day_of_season_sq", "qbo_u50", "nao_daily"]:
        if c in df.columns and df[c].notna().sum() > 100:
            controls.append(c)

    # Response variable
    if response_var not in df.columns:
        LOG.warning("Variable %s not in panel", response_var)
        return None

    # Drop NaN
    lag_cols = [f"event_lag_{l}" for l in range(max_lag + 1)]
    all_cols = [response_var] + lag_cols + controls
    df_clean = df[all_cols].dropna()

    if len(df_clean) < 100:
        LOG.warning("Too few observations for %s: %d", response_var, len(df_clean))
        return None

    Y = df_clean[response_var]
    X = df_clean[lag_cols + controls]
    X = sm.add_constant(X)

    # OLS with Newey-West HAC SEs (bandwidth = max_lag + 1)
    model = OLS(Y, X).fit(cov_type="HAC",
                           cov_kwds={"maxlags": max_lag + 1})

    # Extract lag coefficients
    lag_coefs = []
    lag_pvals = []
    lag_ci_lo = []
    lag_ci_hi = []
    ci = model.conf_int(alpha=0.05)

    for lag in range(max_lag + 1):
        col = f"event_lag_{lag}"
        lag_coefs.append(float(model.params[col]))
        lag_pvals.append(float(model.pvalues[col]))
        lag_ci_lo.append(float(ci.loc[col, 0]))
        lag_ci_hi.append(float(ci.loc[col, 1]))

    # Cumulative effect in pre-specified windows
    def window_effect(start, end):
        cols = [f"event_lag_{l}" for l in range(start, end + 1)]
        coefs = [float(model.params[c]) for c in cols]
        cumulative = sum(coefs)
        # Wald test for joint significance
        R = np.zeros((end - start + 1, len(model.params)))
        for i, c in enumerate(cols):
            R[i, list(model.params.index).index(c)] = 1
        try:
            wald = model.wald_test(R)
            p_joint = float(wald.pvalue)
        except Exception:
            p_joint = np.nan
        return {
            "cumulative_effect": cumulative,
            "mean_effect": cumulative / (end - start + 1),
            "p_joint": p_joint,
        }

    windows = {
        "d0_7": window_effect(0, 7),
        "d5_21": window_effect(5, min(21, max_lag)),
        "d0_30": window_effect(0, max_lag),
    }

    return {
        "response_variable": response_var,
        "n_obs": int(model.nobs),
        "r_squared": float(model.rsquared),
        "lag_coefficients": lag_coefs,
        "lag_pvalues": lag_pvals,
        "lag_ci_lower": lag_ci_lo,
        "lag_ci_upper": lag_ci_hi,
        "window_effects": windows,
        "controls_used": controls,
    }


def main():
    panel = load_panel(winter_only=True)

    # Identify NCEP stratospheric variables
    strat_vars = [c for c in panel.columns if c.startswith("ncep_") and
                  ("10hpa" in c or "30hpa" in c or "50hpa" in c)]

    LOG.info("Stratospheric variables: %s", strat_vars)

    results = {}
    for var in strat_vars:
        LOG.info("Distributed-lag model for %s ...", var)
        res = distributed_lag_model(panel, var, max_lag=30)
        if res is not None:
            results[var] = res
            w = res["window_effects"]
            LOG.info("  d5-21 effect: %.4f (p_joint=%.3f)",
                     w["d5_21"]["mean_effect"], w["d5_21"]["p_joint"])

    # Save
    out = RESULTS / "stratospheric_propagation.json"
    out.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    LOG.info("Results saved to %s", out)

    # Summary
    print(f"\n{'='*70}")
    print("STRATOSPHERIC PROPAGATION — Distributed-Lag Model")
    print(f"{'='*70}")
    for var, res in results.items():
        print(f"\n{var}:")
        print(f"  N={res['n_obs']}, R²={res['r_squared']:.4f}")
        for wname, w in res["window_effects"].items():
            sig = "***" if w["p_joint"] < 0.01 else "**" if w["p_joint"] < 0.05 else "*" if w["p_joint"] < 0.1 else ""
            print(f"  {wname}: cumulative={w['cumulative_effect']:.4f}, "
                  f"mean={w['mean_effect']:.4f}, p_joint={w['p_joint']:.3f} {sig}")


if __name__ == "__main__":
    main()
