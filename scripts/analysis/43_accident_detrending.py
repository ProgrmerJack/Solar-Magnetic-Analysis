"""
Script 43: Formal accident detrending (Cochran-Mantel-Haenszel stratified by decade)

Addresses R5's blocking weakness: the pre-2000 RR<1 confound.
Tests whether SSW-accident association holds after controlling for secular recreation trend.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
from pathlib import Path

DATA = Path("data/processed")
OUT = Path("data/results")

# Load data
accidents = pd.read_parquet(DATA / "cryosphere/slf_accidents.parquet")
ssw_cat = pd.read_parquet(DATA / "atmospheric/ssw_catalog.parquet")

# Reset index to get date column
acc = accidents.reset_index()
acc['date'] = pd.to_datetime(acc['date']).dt.tz_localize(None)

# Winter season assignment (Nov-Apr)
def get_winter(dt):
    if dt.month >= 11:
        return f"{dt.year}/{dt.year+1}"
    elif dt.month <= 4:
        return f"{dt.year-1}/{dt.year}"
    return None

acc['winter'] = acc['date'].apply(get_winter)
acc = acc.dropna(subset=['winter'])

# SSW windows
ssw_windows = {}
for onset in ssw_cat.index:
    onset_dt = pd.Timestamp(onset).tz_localize(None)
    ssw_windows[onset_dt] = (onset_dt - pd.Timedelta(days=15), onset_dt + pd.Timedelta(days=15))

def is_ssw_day(dt):
    for onset, (start, end) in ssw_windows.items():
        if start <= dt <= end:
            return True
    return False

acc['is_ssw'] = acc['date'].apply(is_ssw_day)

# Assign decades
acc['decade'] = (acc['date'].dt.year // 10) * 10

# Count accidents per day in winter months
# Create daily grid
date_range = pd.date_range(start='1970-01-01', end='2025-04-30', freq='D')
daily = pd.DataFrame({'date': date_range})
daily = daily[daily['date'].dt.month.isin([11, 12, 1, 2, 3, 4])]
daily['is_ssw'] = daily['date'].apply(is_ssw_day)
daily['decade'] = (daily['date'].dt.year // 10) * 10

# Count accidents per day
acc_counts = acc.groupby('date').size().reset_index(name='n_acc')
daily = daily.merge(acc_counts, on='date', how='left')
daily['n_acc'] = daily['n_acc'].fillna(0).astype(int)

# Has accident (binary for CMH)
daily['has_acc'] = (daily['n_acc'] > 0).astype(int)

print("=== DECADE-STRATIFIED ANALYSIS ===\n")

# Cochran-Mantel-Haenszel-like analysis
# For each decade, compute 2x2 table: SSW/non-SSW × accident/no-accident
cmh_tables = []
decade_results = {}

for decade in sorted(daily['decade'].unique()):
    d = daily[daily['decade'] == decade]
    
    ssw_d = d[d['is_ssw']]
    non_d = d[~d['is_ssw']]
    
    if len(ssw_d) < 10 or len(non_d) < 10:
        continue
    
    # Binary: day with ≥1 accident
    a = ssw_d['has_acc'].sum()  # SSW days with accidents
    b = len(ssw_d) - a          # SSW days without
    c = non_d['has_acc'].sum()  # non-SSW days with accidents
    d_val = len(non_d) - c      # non-SSW days without
    
    # Rate-based comparison
    ssw_rate = ssw_d['n_acc'].mean()
    non_rate = non_d['n_acc'].mean()
    rr = ssw_rate / non_rate if non_rate > 0 else np.nan
    
    # Store 2x2 table for CMH
    table = np.array([[a, b], [c, d_val]])
    cmh_tables.append(table)
    
    decade_results[str(decade)] = {
        'ssw_days': int(len(ssw_d)),
        'non_ssw_days': int(len(non_d)),
        'ssw_acc_rate': float(ssw_rate),
        'non_ssw_acc_rate': float(non_rate),
        'rr': float(rr),
        'table_2x2': [[int(a), int(b)], [int(c), int(d_val)]]
    }
    
    print(f"{decade}s: SSW rate={ssw_rate:.3f}, non-SSW rate={non_rate:.3f}, RR={rr:.3f} (n_ssw={len(ssw_d)}, n_non={len(non_d)})")

# Manual CMH computation
# CMH odds ratio = Σ(a_i * d_i / n_i) / Σ(b_i * c_i / n_i)
# CMH test statistic
numerator_or = 0
denominator_or = 0
numerator_chi2 = 0
denominator_chi2 = 0

for table in cmh_tables:
    a, b = table[0]
    c, d_val = table[1]
    n = a + b + c + d_val
    n1 = a + b  # SSW row total
    n0 = c + d_val  # non-SSW row total
    m1 = a + c  # accident column total
    m0 = b + d_val  # no-accident column total
    
    numerator_or += (a * d_val) / n
    denominator_or += (b * c) / n
    
    # CMH test
    E_a = n1 * m1 / n
    V_a = n1 * n0 * m1 * m0 / (n**2 * (n - 1))
    numerator_chi2 += (a - E_a)
    denominator_chi2 += V_a

cmh_or = numerator_or / denominator_or if denominator_or > 0 else np.nan
cmh_chi2 = numerator_chi2**2 / denominator_chi2 if denominator_chi2 > 0 else np.nan
cmh_p = 1 - stats.chi2.cdf(cmh_chi2, df=1)

print(f"\n=== COCHRAN-MANTEL-HAENSZEL (DECADE-STRATIFIED) ===")
print(f"CMH Common Odds Ratio: {cmh_or:.3f}")
print(f"CMH chi-square: {cmh_chi2:.3f}")
print(f"CMH P-value: {cmh_p:.4f}")

# Also do rate-ratio detrending: regress accident rate on year, compute SSW residuals
print(f"\n=== DETRENDED RATE-RATIO ===")
# Compute annual winter rates
daily['year'] = daily['date'].dt.year
annual = daily.groupby(['year', 'is_ssw']).agg(
    total_days=('n_acc', 'count'),
    total_acc=('n_acc', 'sum')
).reset_index()
annual['rate'] = annual['total_acc'] / annual['total_days']

# Overall annual rate (for detrending)
annual_all = daily.groupby('year').agg(
    total_days=('n_acc', 'count'),
    total_acc=('n_acc', 'sum')
).reset_index()
annual_all['rate'] = annual_all['total_acc'] / annual_all['total_days']

# Linear trend
slope, intercept, r, p_trend, se = stats.linregress(annual_all['year'], annual_all['rate'])
print(f"Secular trend: slope={slope:.5f} acc/day/yr, R²={r**2:.3f}, P={p_trend:.4f}")

# Detrended rates: subtract predicted trend
daily['predicted_rate'] = intercept + slope * daily['date'].dt.year
daily['detrended_acc'] = daily['n_acc'] - daily['predicted_rate']

# SSW vs non-SSW on detrended residuals
ssw_detrended = daily[daily['is_ssw']]['detrended_acc']
non_detrended = daily[~daily['is_ssw']]['detrended_acc']

t_stat, t_p = stats.ttest_ind(ssw_detrended, non_detrended)
u_stat, u_p = stats.mannwhitneyu(ssw_detrended, non_detrended, alternative='two-sided')

print(f"Detrended SSW residual: {ssw_detrended.mean():.4f}")
print(f"Detrended non-SSW residual: {non_detrended.mean():.4f}")
print(f"Difference: {ssw_detrended.mean() - non_detrended.mean():.4f}")
print(f"T-test: t={t_stat:.3f}, P={t_p:.4f}")
print(f"Mann-Whitney: U={u_stat:.0f}, P={u_p:.4f}")

# Save results
output = {
    'decade_results': decade_results,
    'cmh_common_odds_ratio': float(cmh_or),
    'cmh_chi_square': float(cmh_chi2),
    'cmh_p_value': float(cmh_p),
    'secular_trend': {
        'slope': float(slope),
        'r_squared': float(r**2),
        'p': float(p_trend)
    },
    'detrended': {
        'ssw_mean_residual': float(ssw_detrended.mean()),
        'non_ssw_mean_residual': float(non_detrended.mean()),
        'difference': float(ssw_detrended.mean() - non_detrended.mean()),
        't_test_p': float(t_p),
        'mann_whitney_p': float(u_p)
    }
}

with open(OUT / '43_accident_detrending.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {OUT / '43_accident_detrending.json'}")
