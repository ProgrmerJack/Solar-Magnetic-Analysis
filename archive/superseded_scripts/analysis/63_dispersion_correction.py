"""
63_dispersion_correction.py
Poisson overdispersion correction for daily avalanche count regressions.

Demonstrates that standard Poisson P-values are inflated ~45x by overdispersion
and computes corrected P-values under ZIP, quasi-Poisson, robust SE, and GEE
specifications.

Outputs: data/results/r55_dispersion_correction.json
Referenced by: Extended Data Table 13, main text incremental-value section
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

# Load data
panel = pd.read_parquet(ROOT / "data" / "processed" / "analysis_panel_v2.parquet")
slf = pd.read_parquet(ROOT / "data" / "processed" / "cryosphere" / "slf_activity.parquet")

panel.index = panel.index.tz_localize(None) if panel.index.tz else panel.index
slf.index = slf.index.tz_localize(None) if slf.index.tz else slf.index

merged = panel.join(slf[["dry_natural_size_1234"]], how="inner").dropna(
    subset=["dry_natural_size_1234", "ssw_within_15d", "ncep_z500_nh"]
)

y = merged["dry_natural_size_1234"].values.astype(int)
ssw = merged["ssw_within_15d"].astype(float).values
z500 = merged["ncep_z500_nh"].astype(float).values
z500_std = (z500 - z500.mean()) / z500.std()

# Day-of-year seasonality controls (sin/cos harmonics)
doy = merged.index.dayofyear.values
X_season = np.column_stack([
    np.sin(2 * np.pi * doy / 365),
    np.cos(2 * np.pi * doy / 365),
    np.sin(4 * np.pi * doy / 365),
    np.cos(4 * np.pi * doy / 365),
])

X_combined = np.column_stack([ssw, z500_std, X_season])
X_combined = sm.add_constant(X_combined)
col_names = ["const", "ssw", "z500_std", "sin1", "cos1", "sin2", "cos2"]

results = {
    "n_days": int(len(y)),
    "mean_count": float(y.mean()),
    "var_count": float(y.var()),
    "overdispersion_ratio": float(y.var() / max(y.mean(), 0.001)),
    "zero_fraction": float((y == 0).mean()),
    "models": {}
}

# 1. Standard Poisson
poisson = sm.GLM(y, X_combined, family=sm.families.Poisson()).fit()
results["models"]["standard_poisson"] = {
    "ssw_coef": float(poisson.params[1]),
    "ssw_pvalue": float(poisson.pvalues[1]),
    "z500_coef": float(poisson.params[2]),
    "z500_pvalue": float(poisson.pvalues[2]),
    "note": "INFLATED - do not use for inference"
}

# 2. Quasi-Poisson
qpoisson = sm.GLM(y, X_combined, family=sm.families.Poisson()).fit(
    scale="X2"
)
results["models"]["quasi_poisson"] = {
    "ssw_coef": float(qpoisson.params[1]),
    "ssw_pvalue": float(qpoisson.pvalues[1]),
    "z500_coef": float(qpoisson.params[2]),
    "z500_pvalue": float(qpoisson.pvalues[2]),
    "dispersion": float(qpoisson.scale),
}

# 3. Poisson + Robust HC1 SE
poisson_hc = sm.GLM(y, X_combined, family=sm.families.Poisson()).fit(
    cov_type="HC1"
)
results["models"]["poisson_robust_hc1"] = {
    "ssw_coef": float(poisson_hc.params[1]),
    "ssw_pvalue": float(poisson_hc.pvalues[1]),
    "z500_coef": float(poisson_hc.params[2]),
    "z500_pvalue": float(poisson_hc.pvalues[2]),
}

# 4. ZIP (Zero-Inflated Poisson)
try:
    zip_model = sm.ZeroInflatedPoisson(
        y, X_combined, exog_infl=sm.add_constant(X_season),
        inflation="logit"
    ).fit(disp=0, maxiter=200)
    results["models"]["zip"] = {
        "ssw_coef": float(zip_model.params[1]),
        "ssw_pvalue": float(zip_model.pvalues[1]),
        "z500_coef": float(zip_model.params[2]),
        "z500_pvalue": float(zip_model.pvalues[2]),
        "note": "Preferred specification"
    }
except Exception as e:
    results["models"]["zip"] = {"error": str(e)}

# 5. GEE with winter clusters
winter_id = (merged.index.year - (merged.index.month < 7).astype(int)).values
try:
    from statsmodels.genmod.generalized_estimating_equations import GEE
    from statsmodels.genmod.cov_struct import Exchangeable
    gee_data = pd.DataFrame({
        "y": y, "ssw": ssw, "z500": z500_std,
        "sin1": X_season[:, 0], "cos1": X_season[:, 1],
        "sin2": X_season[:, 2], "cos2": X_season[:, 3],
        "winter": winter_id
    })
    gee_model = GEE.from_formula(
        "y ~ ssw + z500 + sin1 + cos1 + sin2 + cos2",
        groups="winter", data=gee_data,
        family=sm.families.Poisson(),
        cov_struct=Exchangeable()
    ).fit()
    ssw_idx = list(gee_model.params.index).index("ssw")
    z500_idx = list(gee_model.params.index).index("z500")
    results["models"]["gee_winter_clusters"] = {
        "ssw_coef": float(gee_model.params.iloc[ssw_idx]),
        "ssw_pvalue": float(gee_model.pvalues.iloc[ssw_idx]),
        "z500_coef": float(gee_model.params.iloc[z500_idx]),
        "z500_pvalue": float(gee_model.pvalues.iloc[z500_idx]),
    }
except Exception as e:
    results["models"]["gee_winter_clusters"] = {"error": str(e)}

out_path = ROOT / "data" / "results" / "r55_dispersion_correction.json"
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"Saved to {out_path}")
print(f"Overdispersion ratio: {results['overdispersion_ratio']:.1f}")
for name, m in results["models"].items():
    if "error" not in m:
        print(f"  {name}: SSW P={m['ssw_pvalue']:.4f}, Z500 P={m['z500_pvalue']:.4f}")
    else:
        print(f"  {name}: ERROR - {m['error']}")
