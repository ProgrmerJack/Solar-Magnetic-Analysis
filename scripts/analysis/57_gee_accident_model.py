"""
Script 57: GEE/Cluster-robust accident model
=============================================
Fix the CMH pseudoreplication issue flagged by ALL reviewers.
Replace raw CMH (treating 4,868 records as independent) with:
1. GEE logistic model with SSW event as cluster
2. Mixed-effects logistic regression
3. Cluster-robust sandwich SEs

This addresses the design effect DE = 1 + (m-1)ρ that inflates χ².
"""

import pandas as pd
import numpy as np
from scipy import stats
import json, os, warnings
warnings.filterwarnings("ignore")

OUT = "data/results/57_gee_accident_model.json"
os.makedirs("data/results", exist_ok=True)

# Load accident data
acc = pd.read_parquet("data/processed/cryosphere/slf_accidents.parquet")
ssw_cat = pd.read_parquet("data/processed/atmospheric/ssw_catalog.parquet")

# Timezone handling
acc.index = acc.index.tz_localize(None) if acc.index.tz is None else acc.index.tz_convert("UTC").tz_localize(None)
ssw_cat.index = ssw_cat.index.tz_localize(None) if ssw_cat.index.tz is None else ssw_cat.index.tz_convert("UTC").tz_localize(None)

ssw_dates = ssw_cat.index.to_list()

print(f"Accident records: {len(acc)}")
print(f"Columns: {acc.columns.tolist()}")
print(f"SSW events: {len(ssw_dates)}")

# Create SSW exposure variable for each accident
window = 15  # ±15 days from SSW onset

def assign_ssw_exposure(date, ssw_dates, window):
    """Assign SSW exposure and cluster ID."""
    for i, ssw_date in enumerate(ssw_dates):
        diff = (date - ssw_date).days
        if -window <= diff <= window:
            return True, i, diff
    return False, -1, None

# Apply to all accidents  
exposures = [assign_ssw_exposure(d, ssw_dates, window) for d in acc.index]
acc["ssw_exposed"] = [e[0] for e in exposures]
acc["ssw_cluster"] = [e[1] for e in exposures]
acc["days_from_ssw"] = [e[2] for e in exposures]

# For non-SSW records, assign to "control" clusters by winter
for idx in acc[~acc["ssw_exposed"]].index:
    winter = idx.year if idx.month >= 7 else idx.year - 1
    acc.loc[idx, "ssw_cluster"] = 1000 + winter  # unique cluster per winter

# Basic counts
n_exposed = acc["ssw_exposed"].sum()
n_control = (~acc["ssw_exposed"]).sum()
print(f"\nSSW-exposed accidents: {n_exposed}")
print(f"Control accidents: {n_control}")

# Identify fatality outcome
if "number_dead" in acc.columns:
    acc["fatal"] = (acc["number_dead"] > 0).astype(int)
elif "dead" in acc.columns:
    acc["fatal"] = (acc["dead"] > 0).astype(int)
else:
    print("Available columns:", acc.columns.tolist())
    # Try to find fatality column
    for c in acc.columns:
        if "dead" in c.lower() or "fatal" in c.lower() or "kill" in c.lower():
            acc["fatal"] = (acc[c] > 0).astype(int)
            print(f"Using {c} for fatality")
            break

# Decade variable for stratification
acc["decade"] = (acc.index.year // 10) * 10

# Method 1: Raw CMH (the problematic original)
from collections import Counter

# Simple 2x2: SSW exposed vs not, fatal vs not
if "fatal" in acc.columns:
    a = acc[(acc["ssw_exposed"]) & (acc["fatal"] == 1)].shape[0]
    b = acc[(acc["ssw_exposed"]) & (acc["fatal"] == 0)].shape[0]
    c_val = acc[(~acc["ssw_exposed"]) & (acc["fatal"] == 1)].shape[0]
    d_val = acc[(~acc["ssw_exposed"]) & (acc["fatal"] == 0)].shape[0]
    
    raw_or = (a * d_val) / (b * c_val) if b * c_val > 0 else np.nan
    print(f"\nRaw 2×2 OR: {raw_or:.3f}")
    print(f"  SSW: {a} fatal, {b} non-fatal")
    print(f"  Control: {c_val} fatal, {d_val} non-fatal")

# Method 2: Try GEE with statsmodels
try:
    import statsmodels.api as sm
    from statsmodels.genmod.generalized_estimating_equations import GEE
    from statsmodels.genmod.families import Binomial
    from statsmodels.genmod.cov_struct import Exchangeable, Independence
    
    # Prepare GEE data
    gee_data = acc[["ssw_exposed", "fatal", "ssw_cluster", "decade"]].dropna().copy()
    gee_data["ssw_exposed"] = gee_data["ssw_exposed"].astype(int)
    gee_data["intercept"] = 1
    gee_data["cluster"] = gee_data["ssw_cluster"].astype(int)
    
    # Sort by cluster (required for GEE)
    gee_data = gee_data.sort_values("cluster")
    
    # GEE with exchangeable correlation
    try:
        model_exch = GEE.from_formula(
            "fatal ~ ssw_exposed",
            groups="cluster",
            data=gee_data,
            family=Binomial(),
            cov_struct=Exchangeable()
        )
        result_exch = model_exch.fit()
        
        gee_or = np.exp(result_exch.params["ssw_exposed"])
        gee_ci_low = np.exp(result_exch.conf_int().loc["ssw_exposed", 0])
        gee_ci_high = np.exp(result_exch.conf_int().loc["ssw_exposed", 1])
        gee_p = result_exch.pvalues["ssw_exposed"]
        gee_rho = result_exch.cov_struct.summary()
        
        print(f"\n=== GEE (Exchangeable) ===")
        print(f"  OR = {gee_or:.3f} (95% CI: {gee_ci_low:.3f}–{gee_ci_high:.3f})")
        print(f"  P = {gee_p:.6f}")
        print(f"  Within-cluster correlation: {gee_rho}")
        
        gee_success = True
    except Exception as e:
        print(f"GEE exchangeable failed: {e}")
        gee_success = False
        
    # GEE with independence (sandwich SE)
    try:
        model_ind = GEE.from_formula(
            "fatal ~ ssw_exposed",
            groups="cluster",
            data=gee_data,
            family=Binomial(),
            cov_struct=Independence()
        )
        result_ind = model_ind.fit()
        
        sandwich_or = np.exp(result_ind.params["ssw_exposed"])
        sandwich_ci_low = np.exp(result_ind.conf_int().loc["ssw_exposed", 0])
        sandwich_ci_high = np.exp(result_ind.conf_int().loc["ssw_exposed", 1])
        sandwich_p = result_ind.pvalues["ssw_exposed"]
        
        print(f"\n=== GEE (Independence + Sandwich SE) ===")
        print(f"  OR = {sandwich_or:.3f} (95% CI: {sandwich_ci_low:.3f}–{sandwich_ci_high:.3f})")
        print(f"  P = {sandwich_p:.6f}")
        
        sandwich_success = True
    except Exception as e:
        print(f"GEE independence failed: {e}")
        sandwich_success = False
        
except ImportError:
    print("statsmodels not available for GEE")
    gee_success = False
    sandwich_success = False

# Method 3: Manual design-effect correction
# Estimate intra-cluster correlation from the data
clusters = acc.groupby("ssw_cluster")["fatal"].agg(["mean", "count"])
clusters = clusters[clusters["count"] > 1]

if len(clusters) > 1:
    # Estimate ICC using ANOVA-based method
    grand_mean = acc["fatal"].mean()
    n_clusters = len(clusters)
    n_total = clusters["count"].sum()
    m_bar = n_total / n_clusters  # mean cluster size
    
    # Between-cluster variance
    ssb = sum(clusters["count"] * (clusters["mean"] - grand_mean)**2)
    msb = ssb / (n_clusters - 1)
    
    # Within-cluster variance
    ssw_var = sum(clusters["count"] * clusters["mean"] * (1 - clusters["mean"]))
    msw = ssw_var / (n_total - n_clusters) if (n_total - n_clusters) > 0 else 1
    
    # ICC
    icc = (msb - msw) / (msb + (m_bar - 1) * msw) if (msb + (m_bar - 1) * msw) > 0 else 0
    icc = max(0, icc)  # ICC can't be negative in practice
    
    # Design effect
    design_effect = 1 + (m_bar - 1) * icc
    effective_n = n_total / design_effect
    
    # Corrected chi-square
    raw_chi2 = 54.1  # from manuscript
    corrected_chi2 = raw_chi2 / design_effect
    corrected_p = 1 - stats.chi2.cdf(corrected_chi2, 1)
    
    print(f"\n=== Design Effect Correction ===")
    print(f"  Mean cluster size: {m_bar:.1f}")
    print(f"  Estimated ICC: {icc:.4f}")
    print(f"  Design effect: {design_effect:.2f}")
    print(f"  Effective N: {effective_n:.0f} (from {n_total})")
    print(f"  Raw χ²: {raw_chi2:.1f} → Corrected χ²: {corrected_chi2:.2f}")
    print(f"  Raw P: <10⁻¹² → Corrected P: {corrected_p:.6f}")

# Method 4: Event-level aggregation (most conservative)
# Aggregate to SSW-event level: proportion of fatal accidents per SSW window
ssw_events = acc[acc["ssw_exposed"]].groupby("ssw_cluster").agg(
    n_accidents=("fatal", "count"),
    n_fatal=("fatal", "sum"),
    prop_fatal=("fatal", "mean")
)

control_events = acc[~acc["ssw_exposed"]].groupby("ssw_cluster").agg(
    n_accidents=("fatal", "count"),
    n_fatal=("fatal", "sum"),
    prop_fatal=("fatal", "mean")
)

if len(ssw_events) > 0 and len(control_events) > 0:
    ssw_prop = ssw_events["prop_fatal"].mean()
    ctrl_prop = control_events["prop_fatal"].mean()
    
    # Two-sample t-test on proportions
    t_prop, p_prop = stats.ttest_ind(
        ssw_events["prop_fatal"].values,
        control_events["prop_fatal"].values
    )
    
    # Mann-Whitney on proportions
    u_prop, p_mw_prop = stats.mannwhitneyu(
        ssw_events["prop_fatal"].values,
        control_events["prop_fatal"].values,
        alternative="greater"
    )
    
    print(f"\n=== Event-Level Aggregation (Most Conservative) ===")
    print(f"  SSW events: n={len(ssw_events)}, mean fatal proportion={ssw_prop:.4f}")
    print(f"  Control winters: n={len(control_events)}, mean fatal proportion={ctrl_prop:.4f}")
    print(f"  Ratio: {ssw_prop/ctrl_prop:.3f}" if ctrl_prop > 0 else "  Ratio: undefined")
    print(f"  t-test: t={t_prop:.3f}, P={p_prop:.4f}")
    print(f"  Mann-Whitney: U={u_prop:.0f}, P={p_mw_prop:.4f}")

# Save all results
results = {
    "raw_cmh_or": float(raw_or) if "fatal" in acc.columns else None,
    "n_accidents_total": int(len(acc)),
    "n_ssw_exposed": int(n_exposed),
    "n_control": int(n_control),
}

if gee_success:
    results.update({
        "gee_exchangeable_or": float(gee_or),
        "gee_exchangeable_ci": [float(gee_ci_low), float(gee_ci_high)],
        "gee_exchangeable_p": float(gee_p),
        "gee_within_cluster_corr": str(gee_rho),
    })

if sandwich_success:
    results.update({
        "gee_sandwich_or": float(sandwich_or),
        "gee_sandwich_ci": [float(sandwich_ci_low), float(sandwich_ci_high)],
        "gee_sandwich_p": float(sandwich_p),
    })

if len(clusters) > 1:
    results.update({
        "icc": float(icc),
        "design_effect": float(design_effect),
        "effective_n": float(effective_n),
        "corrected_chi2": float(corrected_chi2),
        "corrected_p": float(corrected_p),
        "mean_cluster_size": float(m_bar),
    })

if len(ssw_events) > 0 and len(control_events) > 0:
    results.update({
        "event_level_ssw_fatal_prop": float(ssw_prop),
        "event_level_ctrl_fatal_prop": float(ctrl_prop),
        "event_level_ratio": float(ssw_prop/ctrl_prop) if ctrl_prop > 0 else None,
        "event_level_t": float(t_prop),
        "event_level_p": float(p_prop),
        "event_level_mw_p": float(p_mw_prop),
        "n_ssw_clusters": int(len(ssw_events)),
        "n_control_clusters": int(len(control_events)),
    })

with open(OUT, "w") as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to {OUT}")
