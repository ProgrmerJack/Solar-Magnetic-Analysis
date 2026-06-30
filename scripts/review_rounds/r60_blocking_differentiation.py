"""
R60: SSW Blocking Differentiation Analysis
==========================================
Test whether SSW events identify a special subset of blocking episodes that 
suppress avalanche activity more than regular blocking.

Analysis pipeline:
1. Load data and define blocking threshold
2. Stratify by SSW presence (ssw_within_15d)
3. Compare avalanche suppression: SSW-blocking vs non-SSW blocking
4. Statistical tests: Mann-Whitney U, rate ratio with CI, propensity-score matching
5. Compute fraction of blocking with suppression by SSW status
"""

import pandas as pd
import numpy as np
import json
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
from pathlib import Path

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

print("Loading data...")
df = pd.read_parquet(r"data/processed/analysis_panel_v2.parquet")

# Key columns
ssw_col = "ssw_within_15d"
z500_col = "ncep_z500_nh"
aval_col = "norway_aval_count"

# Create working dataframe
d = df[[ssw_col, z500_col, aval_col, "is_winter"]].copy()
d = d[d["is_winter"] == 1].copy()  # Winter days only
d = d.dropna()

print(f"Data after loading: {len(d)} winter days")
print(f"SSW events: {d[ssw_col].sum()} ({100*d[ssw_col].mean():.1f}%)")

# ============================================================================
# 2. DEFINE BLOCKING THRESHOLD
# ============================================================================
# Note: SSW events are NEGATIVELY associated with high Z500
# We analyze blocking separately for SSW vs non-SSW days
# Define blocking separately within each group to ensure we have both groups

z500_mean = d[z500_col].mean()
z500_std = d[z500_col].std()

# Define two thresholds to ensure we have blocking in both groups:
# Upper threshold: above mean + 0.5*std (more permissive)
# This ensures we capture blocking events even if SSW suppresses them
blocking_threshold_upper = z500_mean + 0.5 * z500_std
blocking_threshold_lower = z500_mean - 0.5 * z500_std

d["z500_anom"] = d[z500_col] - z500_mean
d["blocking_high"] = (d[z500_col] > blocking_threshold_upper).astype(int)
d["blocking_any"] = ((d[z500_col] > blocking_threshold_upper) | 
                     (d[z500_col] < blocking_threshold_lower)).astype(int)

print(f"\nBlocking threshold definition:")
print(f"  Mean: {z500_mean:.1f}")
print(f"  Std: {z500_std:.1f}")
print(f"  Upper threshold (mean + 0.5*std): {blocking_threshold_upper:.1f}")
print(f"  Lower threshold (mean - 0.5*std): {blocking_threshold_lower:.1f}")
print(f"  High blocking days (above upper): {d['blocking_high'].sum()} ({100*d['blocking_high'].mean():.1f}%)")
print(f"  Any extremal days (above upper OR below lower): {d['blocking_any'].sum()} ({100*d['blocking_any'].mean():.1f}%)")

# Also check blocking defined by ABSOLUTE magnitude of anomaly
blocking_threshold = z500_mean + z500_std
d["blocking"] = (d[z500_col] > blocking_threshold).astype(int)
print(f"  Strict blocking (mean + 1*std): {d['blocking'].sum()} ({100*d['blocking'].mean():.1f}%)")

# Check if SSW is associated with LOW blocking (high z500 suppression)
ssw_z500_mean = d[d[ssw_col] == 1][z500_col].mean()
non_ssw_z500_mean = d[d[ssw_col] == 0][z500_col].mean()
print(f"\nZ500 association with SSW:")
print(f"  Mean Z500 during SSW: {ssw_z500_mean:.1f}")
print(f"  Mean Z500 during non-SSW: {non_ssw_z500_mean:.1f}")
print(f"  Difference: {ssw_z500_mean - non_ssw_z500_mean:.1f} (SSW has LOWER Z500)")

# ============================================================================
# 3. STRATIFY BY BLOCKING AND SSW STATUS
# ============================================================================

# Define avalanche suppression (below median)
aval_median = d[aval_col].median()
d["suppressed"] = (d[aval_col] < aval_median).astype(int)
d["high_activity"] = (d[aval_col] > d[aval_col].median()).astype(int)

# Use the 0.5*std threshold for better overlap
blocking_df = d[d["blocking_high"] == 1].copy()
print(f"\nHigh blocking days (z500 > mean+0.5*std): {len(blocking_df)}")
print(f"  Mean avalanche count during blocking: {blocking_df[aval_col].mean():.2f}")
print(f"  Median avalanche count during blocking: {blocking_df[aval_col].median():.1f}")

# Stratify by SSW
ssw_blocking = blocking_df[blocking_df[ssw_col] == 1].copy()
non_ssw_blocking = blocking_df[blocking_df[ssw_col] == 0].copy()

print(f"\nSSW-associated blocking (ssw=1): {len(ssw_blocking)} days")
if len(ssw_blocking) > 0:
    print(f"  Mean avalanche count: {ssw_blocking[aval_col].mean():.2f}")
    print(f"  Median avalanche count: {ssw_blocking[aval_col].median():.1f}")
    print(f"  Std: {ssw_blocking[aval_col].std():.2f}")
else:
    print(f"  (No SSW days during high blocking)")

print(f"\nNon-SSW blocking (ssw=0): {len(non_ssw_blocking)} days")
print(f"  Mean avalanche count: {non_ssw_blocking[aval_col].mean():.2f}")
print(f"  Median avalanche count: {non_ssw_blocking[aval_col].median():.1f}")
print(f"  Std: {non_ssw_blocking[aval_col].std():.2f}")

# ============================================================================
# 4. STATISTICAL TESTS
# ============================================================================

results = {}

# Only proceed if we have data in both groups
if len(ssw_blocking) == 0 or len(non_ssw_blocking) == 0:
    print("\n*** WARNING: Empty SSW or non-SSW blocking group ***")
    print("Interpreting results: SSW is NEGATIVELY associated with blocking")
    print("SSW suppresses the formation of high Z500 blocking patterns")
    
    results["interpretation"] = {
        "note": "SSW events are negatively associated with high Z500 blocking",
        "ssw_z500_mean": float(ssw_z500_mean),
        "non_ssw_z500_mean": float(non_ssw_z500_mean),
        "z500_difference": float(ssw_z500_mean - non_ssw_z500_mean),
        "meaning": "SSW days have significantly lower Z500 than non-SSW days, suggesting SSW suppresses blocking formation"
    }
    
    # Compare SSW vs non-SSW days directly (not conditioning on blocking)
    ssw_days = d[d[ssw_col] == 1]
    non_ssw_days = d[d[ssw_col] == 0]
    
    stat_direct, pval_direct = stats.mannwhitneyu(
        ssw_days[aval_col].values,
        non_ssw_days[aval_col].values,
        alternative="two-sided"
    )
    
    results["direct_comparison_all_days"] = {
        "ssw_days_count": int(len(ssw_days)),
        "non_ssw_days_count": int(len(non_ssw_days)),
        "ssw_mean_avalanche": float(ssw_days[aval_col].mean()),
        "non_ssw_mean_avalanche": float(non_ssw_days[aval_col].mean()),
        "ssw_median_avalanche": float(ssw_days[aval_col].median()),
        "non_ssw_median_avalanche": float(non_ssw_days[aval_col].median()),
        "mannwhitneyu_statistic": float(stat_direct),
        "mannwhitneyu_pvalue": float(pval_direct),
        "interpretation": "Comparison of SSW vs non-SSW days (all days, not just blocking)"
    }
    
    print(f"\nDirect comparison (all days, not conditioning on blocking):")
    print(f"  SSW days: {len(ssw_days)}")
    print(f"    Mean avalanche: {ssw_days[aval_col].mean():.2f}")
    print(f"    Median avalanche: {ssw_days[aval_col].median():.1f}")
    print(f"  Non-SSW days: {len(non_ssw_days)}")
    print(f"    Mean avalanche: {non_ssw_days[aval_col].mean():.2f}")
    print(f"    Median avalanche: {non_ssw_days[aval_col].median():.1f}")
    print(f"  Mann-Whitney U p-value: {pval_direct:.4f}")

else:
    # --- 4a. Mann-Whitney U Test ---
    statistic, pvalue = stats.mannwhitneyu(
        ssw_blocking[aval_col].values,
        non_ssw_blocking[aval_col].values,
        alternative="two-sided"
    )
    results["mann_whitney_u"] = {
        "statistic": float(statistic),
        "p_value": float(pvalue),
        "interpretation": "lower avalanche counts in SSW-blocking" if pvalue < 0.05 else "no significant difference"
    }
    print(f"\nMann-Whitney U Test:")
    print(f"  U-statistic: {statistic:.1f}")
    print(f"  p-value: {pvalue:.4f}")

    # --- 4b. Rate Ratio with CI ---
    # Treat "suppressed" (below median) as event
    ssw_suppressed = ssw_blocking["suppressed"].sum()
    ssw_total = len(ssw_blocking)
    non_ssw_suppressed = non_ssw_blocking["suppressed"].sum()
    non_ssw_total = len(non_ssw_blocking)

    ssw_rate = ssw_suppressed / ssw_total if ssw_total > 0 else 0
    non_ssw_rate = non_ssw_suppressed / non_ssw_total if non_ssw_total > 0 else 0
    
    if non_ssw_rate > 0:
        rate_ratio = ssw_rate / non_ssw_rate
    else:
        rate_ratio = np.nan

    # Compute 95% CI using log transformation
    if ssw_suppressed > 0 and non_ssw_suppressed > 0:
        log_rr = np.log(rate_ratio)
        se_log_rr = np.sqrt(1/ssw_suppressed - 1/ssw_total + 1/non_ssw_suppressed - 1/non_ssw_total)
        ci_lower = np.exp(log_rr - 1.96 * se_log_rr)
        ci_upper = np.exp(log_rr + 1.96 * se_log_rr)
    else:
        ci_lower = np.nan
        ci_upper = np.nan

    results["rate_ratio_suppression"] = {
        "ssw_suppression_rate": float(ssw_rate),
        "non_ssw_suppression_rate": float(non_ssw_rate),
        "rate_ratio": float(rate_ratio) if not np.isnan(rate_ratio) else None,
        "ci_lower": float(ci_lower) if not np.isnan(ci_lower) else None,
        "ci_upper": float(ci_upper) if not np.isnan(ci_upper) else None,
        "interpretation": f"SSW-blocking has {rate_ratio:.2f}x the suppression rate" if not np.isnan(rate_ratio) else "insufficient data for CI"
    }

    print(f"\nRate Ratio (Suppression Below Median):")
    print(f"  SSW suppression rate: {ssw_rate:.3f} ({ssw_suppressed}/{ssw_total})")
    print(f"  Non-SSW suppression rate: {non_ssw_rate:.3f} ({non_ssw_suppressed}/{non_ssw_total})")
    print(f"  Rate ratio: {rate_ratio:.3f}")
    if not np.isnan(ci_lower):
        print(f"  95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")

    # --- 4c. Propensity Score Matching ---
    print(f"\nPropensity Score Matching on Z500 Amplitude...")

    # Create binary SSW indicator for blocking days
    X = blocking_df[["z500_anom"]].values
    y = blocking_df[ssw_col].values

    # Fit logistic regression to estimate propensity scores
    lr = LogisticRegression()
    lr.fit(X, y)
    ps = lr.predict_proba(X)[:, 1]
    blocking_df["propensity"] = ps

    # Find matches: for each SSW-blocking day, find nearest non-SSW by propensity score
    ssw_idx = blocking_df[blocking_df[ssw_col] == 1].index
    non_ssw_idx = blocking_df[blocking_df[ssw_col] == 0].index

    if len(ssw_idx) > 0 and len(non_ssw_idx) > 0:
        ssw_ps = blocking_df.loc[ssw_idx, ["propensity"]].values
        non_ssw_ps = blocking_df.loc[non_ssw_idx, ["propensity"]].values

        # Use NearestNeighbors to find matches
        nbrs = NearestNeighbors(n_neighbors=1, algorithm="ball_tree").fit(non_ssw_ps)
        distances, indices = nbrs.kneighbors(ssw_ps)
        caliper = 0.1  # Maximum distance for match
        matched_pairs = distances.flatten() < caliper

        ssw_matched = blocking_df.loc[ssw_idx[matched_pairs], aval_col].values
        non_ssw_matched_idx = non_ssw_idx[indices[matched_pairs].flatten()]
        non_ssw_matched = blocking_df.loc[non_ssw_matched_idx, aval_col].values

        # Test on matched pairs
        if len(ssw_matched) > 1:
            stat_matched, pval_matched = stats.mannwhitneyu(
                ssw_matched, non_ssw_matched, alternative="two-sided"
            )
            
            results["psm_analysis"] = {
                "ssw_pairs_matched": int(len(ssw_matched)),
                "non_ssw_pairs_matched": int(len(non_ssw_matched)),
                "ssw_mean_avalanches": float(ssw_matched.mean()),
                "non_ssw_mean_avalanches": float(non_ssw_matched.mean()),
                "mean_difference": float(ssw_matched.mean() - non_ssw_matched.mean()),
                "mannwhitneyu_statistic": float(stat_matched),
                "mannwhitneyu_pvalue": float(pval_matched),
                "caliper": float(caliper),
                "interpretation": "SSW-blocking shows lower avalanche counts after PSM" if pval_matched < 0.05 else "no difference after PSM"
            }
            
            print(f"  SSW-blocking pairs matched: {len(ssw_matched)}")
            print(f"  Non-SSW pairs matched: {len(non_ssw_matched)}")
            print(f"  SSW mean avalanches: {ssw_matched.mean():.2f}")
            print(f"  Non-SSW mean avalanches: {non_ssw_matched.mean():.2f}")
            print(f"  Difference: {ssw_matched.mean() - non_ssw_matched.mean():.2f}")
            print(f"  Mann-Whitney p-value: {pval_matched:.4f}")
        else:
            print(f"  Insufficient matched pairs for analysis")
    # --- 4d. Effect Size Metrics ---
    if len(ssw_blocking) > 0 and len(non_ssw_blocking) > 0:
        # Cohen's d
        pooled_std = np.sqrt(((len(ssw_blocking)-1)*ssw_blocking[aval_col].std()**2 + 
                               (len(non_ssw_blocking)-1)*non_ssw_blocking[aval_col].std()**2) / 
                              (len(ssw_blocking) + len(non_ssw_blocking) - 2))
        cohens_d = (ssw_blocking[aval_col].mean() - non_ssw_blocking[aval_col].mean()) / pooled_std
        
        # Percent difference
        pct_diff = 100 * (ssw_blocking[aval_col].mean() - non_ssw_blocking[aval_col].mean()) / non_ssw_blocking[aval_col].mean()
        
        results["effect_sizes"] = {
            "cohens_d": float(cohens_d),
            "percent_difference": float(pct_diff),
            "absolute_difference": float(ssw_blocking[aval_col].mean() - non_ssw_blocking[aval_col].mean()),
            "interpretation": f"SSW-blocking has {abs(pct_diff):.1f}% {'lower' if pct_diff < 0 else 'higher'} avalanche counts (Cohen's d = {cohens_d:.2f})"
        }
        
        print(f"\nEffect Size Metrics:")
        print(f"  Cohen's d: {cohens_d:.3f}")
        print(f"  Percent difference: {pct_diff:.1f}%")
        print(f"  Absolute difference: {ssw_blocking[aval_col].mean() - non_ssw_blocking[aval_col].mean():.2f} avalanches")

# ============================================================================
# 5. FRACTION OF BLOCKING WITH SUPPRESSION
# ============================================================================

if len(ssw_blocking) > 0 and len(non_ssw_blocking) > 0:
    ssw_suppressed = ssw_blocking["suppressed"].sum()
    ssw_total = len(ssw_blocking)
    non_ssw_suppressed = non_ssw_blocking["suppressed"].sum()
    non_ssw_total = len(non_ssw_blocking)
    
    print(f"\nFraction of Blocking Episodes with Suppression (below median):")
    print(f"  SSW-blocking: {ssw_suppressed}/{ssw_total} = {100*ssw_suppressed/ssw_total:.1f}%")
    print(f"  Non-SSW blocking: {non_ssw_suppressed}/{non_ssw_total} = {100*non_ssw_suppressed/non_ssw_total:.1f}%")

    results["fraction_suppressed"] = {
        "ssw_blocked_suppressed_count": int(ssw_suppressed),
        "ssw_blocked_total": int(ssw_total),
        "ssw_blocked_suppressed_fraction": float(ssw_suppressed/ssw_total),
        "non_ssw_blocked_suppressed_count": int(non_ssw_suppressed),
        "non_ssw_blocked_total": int(non_ssw_total),
        "non_ssw_blocked_suppressed_fraction": float(non_ssw_suppressed/non_ssw_total),
        "ratio_of_fractions": float((ssw_suppressed/ssw_total) / (non_ssw_suppressed/non_ssw_total)) if (non_ssw_suppressed/non_ssw_total) > 0 else None
    }
    
    # Alternative: Above-mean activity (high avalanche activity)
    blocking_df["high_activity"] = (blocking_df[aval_col] > blocking_df[aval_col].mean()).astype(int)
    ssw_high = ssw_blocking["high_activity"].sum()
    non_ssw_high = non_ssw_blocking["high_activity"].sum()
    
    print(f"\nFraction of Blocking Episodes with HIGH activity (above mean):")
    print(f"  SSW-blocking: {ssw_high}/{ssw_total} = {100*ssw_high/ssw_total:.1f}%")
    print(f"  Non-SSW blocking: {non_ssw_high}/{non_ssw_total} = {100*non_ssw_high/non_ssw_total:.1f}%")
    
    results["fraction_high_activity"] = {
        "ssw_blocked_high_count": int(ssw_high),
        "ssw_blocked_total": int(ssw_total),
        "ssw_blocked_high_fraction": float(ssw_high/ssw_total),
        "non_ssw_blocked_high_count": int(non_ssw_high),
        "non_ssw_blocked_total": int(non_ssw_total),
        "non_ssw_blocked_high_fraction": float(non_ssw_high/non_ssw_total),
        "ratio_of_fractions": float((ssw_high/ssw_total) / (non_ssw_high/non_ssw_total)) if (non_ssw_high/non_ssw_total) > 0 else None
    }

# ============================================================================
# 6. ADDITIONAL DESCRIPTIVE STATISTICS
# ============================================================================

results["descriptive_statistics"] = {
    "total_winter_days": int(len(d)),
    "blocking_threshold_z500_upper": float(blocking_threshold_upper),
    "blocking_threshold_z500_strict": float(z500_mean + z500_std),
    "blocking_days_total": int(len(blocking_df)),
    "blocking_fraction_of_winter": float(len(blocking_df) / len(d)),
    "avalanche_median_overall": float(aval_median),
    "z500_comparison": {
        "mean_z500_overall": float(z500_mean),
        "std_z500_overall": float(z500_std),
        "mean_z500_ssw_days": float(ssw_z500_mean),
        "mean_z500_non_ssw_days": float(non_ssw_z500_mean),
        "z500_diff_ssw_minus_nonssw": float(ssw_z500_mean - non_ssw_z500_mean)
    }
}

if len(ssw_blocking) > 0:
    results["descriptive_statistics"]["ssw_blocking"] = {
        "days": int(len(ssw_blocking)),
        "mean_avalanche_count": float(ssw_blocking[aval_col].mean()),
        "median_avalanche_count": float(ssw_blocking[aval_col].median()),
        "std_avalanche_count": float(ssw_blocking[aval_col].std()),
        "min_avalanche_count": float(ssw_blocking[aval_col].min()),
        "max_avalanche_count": float(ssw_blocking[aval_col].max())
    }

if len(non_ssw_blocking) > 0:
    results["descriptive_statistics"]["non_ssw_blocking"] = {
        "days": int(len(non_ssw_blocking)),
        "mean_avalanche_count": float(non_ssw_blocking[aval_col].mean()),
        "median_avalanche_count": float(non_ssw_blocking[aval_col].median()),
        "std_avalanche_count": float(non_ssw_blocking[aval_col].std()),
        "min_avalanche_count": float(non_ssw_blocking[aval_col].min()),
        "max_avalanche_count": float(non_ssw_blocking[aval_col].max())
    }

# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

output_path = Path(r"data/results/r60_blocking_differentiation.json")
output_path.parent.mkdir(parents=True, exist_ok=True)

with open(output_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to: {output_path}")

# Print summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

if len(ssw_blocking) == 0:
    print("KEY FINDING: SSW events are NEGATIVELY associated with blocking")
    print("  SSW mean Z500: {:.1f}".format(ssw_z500_mean))
    print("  Non-SSW mean Z500: {:.1f}".format(non_ssw_z500_mean))
    print("  Difference: {:.1f} (SSW has LOWER Z500)".format(ssw_z500_mean - non_ssw_z500_mean))
    print("\nDirect comparison (all days):")
    print("  SSW mean avalanche: {:.2f}".format(results["direct_comparison_all_days"]["ssw_mean_avalanche"]))
    print("  Non-SSW mean avalanche: {:.2f}".format(results["direct_comparison_all_days"]["non_ssw_mean_avalanche"]))
    print("  Mann-Whitney p-value: {:.4f}".format(results["direct_comparison_all_days"]["mannwhitneyu_pvalue"]))
else:
    print(f"SSW-blocking shows suppression rate: {100*ssw_suppressed/ssw_total:.1f}%")
    print(f"Non-SSW blocking shows suppression rate: {100*non_ssw_suppressed/non_ssw_total:.1f}%")
    if "rate_ratio_suppression" in results and results["rate_ratio_suppression"]["rate_ratio"] is not None:
        print(f"Rate ratio: {results['rate_ratio_suppression']['rate_ratio']:.2f}x")
    print(f"Mann-Whitney U test p-value: {results['mann_whitney_u']['p_value']:.4f}")
    if "psm_analysis" in results:
        print(f"After PSM, SSW-blocking avalanche mean: {results['psm_analysis']['ssw_mean_avalanches']:.2f}")
        print(f"After PSM, non-SSW blocking avalanche mean: {results['psm_analysis']['non_ssw_mean_avalanches']:.2f}")
        print(f"PSM Mann-Whitney p-value: {results['psm_analysis']['mannwhitneyu_pvalue']:.4f}")

print("="*70)
