#!/usr/bin/env python3
"""
Script 31: Stratospheric Downward Propagation & Avalanche Response
===================================================================
Computes per-event stratospheric metrics and tests whether SSWs with
stronger downward propagation produce greater avalanche suppression.

Key analyses:
1. Per-event vortex deceleration rate (Deltau10 from pre-SSW to post-SSW)
2. Downward propagation index (T anomaly propagation from 10->100 hPa)
3. Cumulative polar cap warming (integral of T10 anomaly)
4. Group comparison: propagating vs non-propagating SSWs
5. Multi-predictor model: does stratospheric info add value beyond Z500?
6. Aura MLS independent validation of stratospheric state

Output: data/results/downward_propagation.json
        figures/downward_propagation_composite.pdf (if matplotlib available)
"""

import pandas as pd
import numpy as np
from scipy import stats
import json
import os
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(BASE, 'data', 'results')

# ── Load data ──────────────────────────────────────────────────────
panel = pd.read_parquet(os.path.join(BASE, 'data', 'processed', 'analysis_panel_v2.parquet'))
dose = pd.read_csv(os.path.join(RESULTS, 'event_level_dose_response.csv'))
catalog = pd.read_csv(os.path.join(RESULTS, 'ssw_event_catalog.csv'))

ssw_dates = pd.to_datetime(dose['date'])
log_rr = dose['log_rr'].values
rr = dose['rr'].values

# Pressure levels available
pressure_levels = [10, 20, 30, 50, 70, 100]  # hPa
t_cols = [f'ncep_t_{p}hpa' for p in pressure_levels]
u_cols = [f'ncep_u_{p}hpa' for p in pressure_levels]
z_cols = [f'ncep_z_{p}hpa' for p in pressure_levels]

# ── 1. Per-event stratospheric metrics ─────────────────────────────
print("=" * 70)
print("1. COMPUTING PER-EVENT STRATOSPHERIC METRICS")
print("=" * 70)

event_metrics = []
panel.index = pd.to_datetime(panel.index)

for i, onset in enumerate(ssw_dates):
    # Windows
    pre_window = panel.loc[(panel.index >= onset - pd.Timedelta(days=30)) &
                           (panel.index < onset - pd.Timedelta(days=5))]
    onset_window = panel.loc[(panel.index >= onset - pd.Timedelta(days=5)) &
                             (panel.index <= onset + pd.Timedelta(days=5))]
    post_window = panel.loc[(panel.index > onset + pd.Timedelta(days=5)) &
                            (panel.index <= onset + pd.Timedelta(days=30))]
    late_window = panel.loc[(panel.index > onset + pd.Timedelta(days=15)) &
                            (panel.index <= onset + pd.Timedelta(days=45))]

    # Climatological baseline (DJF means)
    doy = onset.day_of_year
    clim_mask = (panel.index.month.isin([12, 1, 2]))
    clim = panel.loc[clim_mask]

    m = {}
    m['date'] = str(onset.date())
    m['log_rr'] = log_rr[i]
    m['rr'] = rr[i]

    # a) Vortex deceleration: change in u10 from pre to post
    if len(pre_window) > 0 and len(post_window) > 0:
        m['u10_pre'] = pre_window['ncep_u_10hpa'].mean()
        m['u10_post'] = post_window['ncep_u_10hpa'].mean()
        m['u10_decel'] = m['u10_pre'] - m['u10_post']  # positive = deceleration
        m['u10_onset'] = onset_window['ncep_u_10hpa'].mean() if len(onset_window) > 0 else np.nan
    else:
        m['u10_pre'] = m['u10_post'] = m['u10_decel'] = m['u10_onset'] = np.nan

    # b) Temperature anomalies at each level (relative to DJF mean)
    for p in pressure_levels:
        col = f'ncep_t_{p}hpa'
        clim_mean = clim[col].mean()
        if len(post_window) > 0:
            m[f't{p}_anom_post'] = post_window[col].mean() - clim_mean
        else:
            m[f't{p}_anom_post'] = np.nan
        if len(onset_window) > 0:
            m[f't{p}_anom_onset'] = onset_window[col].mean() - clim_mean
        else:
            m[f't{p}_anom_onset'] = np.nan

    # c) Downward propagation index: T anomaly at 100hPa relative to 10hPa
    if not np.isnan(m.get('t10_anom_post', np.nan)) and not np.isnan(m.get('t100_anom_post', np.nan)):
        # If T warming propagates downward, both should be positive
        m['propagation_ratio'] = m['t100_anom_post'] / m['t10_anom_post'] if m['t10_anom_post'] != 0 else np.nan
    else:
        m['propagation_ratio'] = np.nan

    # d) Cumulative warming at 10hPa (integral over 0 to +30 days)
    event_30d = panel.loc[(panel.index >= onset) & (panel.index <= onset + pd.Timedelta(days=30))]
    if len(event_30d) > 0:
        m['cumulative_t10'] = (event_30d['ncep_t_10hpa'] - clim['ncep_t_10hpa'].mean()).sum()
    else:
        m['cumulative_t10'] = np.nan

    # e) Wind reversal duration at 10hPa (number of days with u10 < 0)
    event_60d = panel.loc[(panel.index >= onset) & (panel.index <= onset + pd.Timedelta(days=60))]
    if len(event_60d) > 0:
        m['u10_reversal_days'] = (event_60d['ncep_u_10hpa'] < 0).sum()
    else:
        m['u10_reversal_days'] = np.nan

    # f) Lower-stratosphere coupling: u at 100hPa change
    if len(pre_window) > 0 and len(post_window) > 0:
        m['u100_decel'] = pre_window['ncep_u_100hpa'].mean() - post_window['ncep_u_100hpa'].mean()
        m['u100_post'] = post_window['ncep_u_100hpa'].mean()
    else:
        m['u100_decel'] = m['u100_post'] = np.nan

    # g) Peak temperature anomaly time at each level (lag of max T)
    event_full = panel.loc[(panel.index >= onset - pd.Timedelta(days=5)) &
                           (panel.index <= onset + pd.Timedelta(days=45))]
    if len(event_full) > 5:
        for p in pressure_levels:
            col = f'ncep_t_{p}hpa'
            t_series = event_full[col] - clim[col].mean()
            peak_idx = t_series.idxmax()
            m[f't{p}_peak_lag'] = (peak_idx - onset).days
    else:
        for p in pressure_levels:
            m[f't{p}_peak_lag'] = np.nan

    # h) Surface coupling metrics
    if len(late_window) > 0:
        m['z500_late'] = late_window['ncep_z500_nh'].mean() if 'ncep_z500_nh' in late_window.columns else np.nan
        m['slp_late'] = late_window['ncep_slp_nh'].mean() if 'ncep_slp_nh' in late_window.columns else np.nan
        m['u850_late'] = late_window['ncep_u850_nh'].mean() if 'ncep_u850_nh' in late_window.columns else np.nan
    else:
        m['z500_late'] = m['slp_late'] = m['u850_late'] = np.nan

    # i) Propagation depth: lowest level where T anomaly > 1 K
    propagation_depth = 10  # starts at 10 hPa
    for p in [20, 30, 50, 70, 100]:
        if m.get(f't{p}_anom_post', 0) > 1.0:
            propagation_depth = p
    m['propagation_depth_hpa'] = propagation_depth

    # j) Downward propagation speed (lag difference between 10hPa and deepest level peak)
    t10_lag = m.get('t10_peak_lag', np.nan)
    deepest_lag = m.get(f't{propagation_depth}_peak_lag', np.nan)
    if not np.isnan(t10_lag) and not np.isnan(deepest_lag) and propagation_depth > 10:
        m['propagation_speed_days'] = deepest_lag - t10_lag
    else:
        m['propagation_speed_days'] = np.nan

    event_metrics.append(m)

df_events = pd.DataFrame(event_metrics)
print(f"\nComputed metrics for {len(df_events)} SSW events")
print(f"Columns: {list(df_events.columns)}")

# ── 2. Correlation analysis ────────────────────────────────────────
print("\n" + "=" * 70)
print("2. STRATOSPHERIC PREDICTORS vs AVALANCHE RESPONSE (log RR)")
print("=" * 70)

strat_predictors = [
    'u10_decel', 'u10_reversal_days', 'cumulative_t10',
    't10_anom_post', 't50_anom_post', 't100_anom_post',
    'u100_decel', 'propagation_depth_hpa',
    'propagation_ratio', 'propagation_speed_days',
    'z500_late', 'u850_late'
]

correlation_results = {}
for pred in strat_predictors:
    vals = df_events[pred].values
    mask = ~np.isnan(vals) & ~np.isnan(log_rr)
    if mask.sum() >= 5:
        r_p, p_p = stats.pearsonr(vals[mask], log_rr[mask])
        r_s, p_s = stats.spearmanr(vals[mask], log_rr[mask])
        correlation_results[pred] = {
            'n': int(mask.sum()),
            'pearson_r': round(r_p, 4),
            'pearson_p': round(p_p, 4),
            'spearman_rho': round(r_s, 4),
            'spearman_p': round(p_s, 4),
        }
        sig = '*' if p_s < 0.10 else ''
        print(f"  {pred:30s}: rho={r_s:+.3f} P={p_s:.4f} {sig}  (n={mask.sum()})")

# ── 3. Group comparison: propagating vs non-propagating SSWs ───────
print("\n" + "=" * 70)
print("3. GROUP COMPARISON: PROPAGATING vs NON-PROPAGATING SSWs")
print("=" * 70)

# Define propagating events: those where T anomaly reaches at least 50 hPa
# (i.e., propagation_depth >= 50)
propagating = df_events['propagation_depth_hpa'] >= 50
prop_rr = log_rr[propagating.values]
nonprop_rr = log_rr[~propagating.values]

print(f"\nPropagating SSWs (depth >= 50 hPa): n = {propagating.sum()}")
print(f"  Mean log(RR) = {prop_rr.mean():.3f} (RR = {np.exp(prop_rr.mean()):.3f})")
print(f"Non-propagating SSWs (depth < 50 hPa): n = (~propagating).sum()")
print(f"  Mean log(RR) = {nonprop_rr.mean():.3f} (RR = {np.exp(nonprop_rr.mean()):.3f})")

# Mann-Whitney U test
if len(prop_rr) >= 3 and len(nonprop_rr) >= 3:
    u_stat, u_p = stats.mannwhitneyu(prop_rr, nonprop_rr, alternative='less')
    print(f"\nMann-Whitney U: U={u_stat:.1f}, P={u_p:.4f}")
    group_mw_p = u_p
else:
    group_mw_p = np.nan
    print("\nInsufficient group sizes for Mann-Whitney test")

# Effect size (Cohen's d)
pooled_sd = np.sqrt((prop_rr.var() * (len(prop_rr)-1) + nonprop_rr.var() * (len(nonprop_rr)-1)) / (len(prop_rr) + len(nonprop_rr) - 2))
cohens_d = (prop_rr.mean() - nonprop_rr.mean()) / pooled_sd if pooled_sd > 0 else np.nan
print(f"Cohen's d = {cohens_d:.3f}")

# Also try tertile split on u10_decel
print("\n--- Tertile split on vortex deceleration (u10_decel) ---")
u10_decel = df_events['u10_decel'].values
terciles = np.nanpercentile(u10_decel, [33, 67])
strong = log_rr[u10_decel > terciles[1]]
weak = log_rr[u10_decel < terciles[0]]
print(f"Strong deceleration (top tercile): n={len(strong)}, mean log(RR)={strong.mean():.3f}")
print(f"Weak deceleration (bottom tercile): n={len(weak)}, mean log(RR)={weak.mean():.3f}")
if len(strong) >= 3 and len(weak) >= 3:
    u2, p2 = stats.mannwhitneyu(strong, weak, alternative='two-sided')
    print(f"Mann-Whitney U: U={u2:.1f}, P={p2:.4f}")

# Tertile split on cumulative T10
print("\n--- Tertile split on cumulative T10 warming ---")
cum_t10 = df_events['cumulative_t10'].values
terc_t = np.nanpercentile(cum_t10, [33, 67])
strong_t = log_rr[cum_t10 > terc_t[1]]
weak_t = log_rr[cum_t10 < terc_t[0]]
print(f"Strong warming (top tercile): n={len(strong_t)}, mean log(RR)={strong_t.mean():.3f}")
print(f"Weak warming (bottom tercile): n={len(weak_t)}, mean log(RR)={weak_t.mean():.3f}")
if len(strong_t) >= 3 and len(weak_t) >= 3:
    u3, p3 = stats.mannwhitneyu(strong_t, weak_t, alternative='two-sided')
    print(f"Mann-Whitney U: U={u3:.1f}, P={p3:.4f}")

# ── 4. Downward propagation composite ─────────────────────────────
print("\n" + "=" * 70)
print("4. DOWNWARD PROPAGATION COMPOSITE (T anomaly by level and lag)")
print("=" * 70)

# Compute composite T anomalies for each pressure level and lag
lags = range(-30, 46)
composite = {}
for p in pressure_levels:
    col = f'ncep_t_{p}hpa'
    clim_mean = panel.loc[panel.index.month.isin([12, 1, 2]), col].mean()
    lag_means = []
    for lag in lags:
        vals = []
        for onset in ssw_dates:
            target_date = onset + pd.Timedelta(days=lag)
            if target_date in panel.index:
                vals.append(panel.loc[target_date, col] - clim_mean)
        lag_means.append(np.mean(vals) if vals else np.nan)
    composite[f't{p}hpa'] = lag_means

print("Composite computed for 6 levels x 76 lag days")

# Peak lag at each level
for p in pressure_levels:
    key = f't{p}hpa'
    series = np.array(composite[key])
    if not np.all(np.isnan(series)):
        peak_lag = list(lags)[np.nanargmax(series)]
        peak_val = np.nanmax(series)
        print(f"  {p:4d} hPa: peak at lag {peak_lag:+3d} days, anomaly = {peak_val:+.2f} K")

# ── 5. Multi-predictor analysis ───────────────────────────────────
print("\n" + "=" * 70)
print("5. INCREMENTAL VALUE OF STRATOSPHERIC PREDICTORS")
print("=" * 70)

from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_squared_error

z500 = dose['z500_nh_m'].values
y = log_rr

# Model 1: Z500 alone
lr1 = LinearRegression().fit(z500.reshape(-1, 1), y)
r2_z500 = lr1.score(z500.reshape(-1, 1), y)

# Model 2: Z500 + u10_decel
X2 = np.column_stack([z500, df_events['u10_decel'].fillna(0).values])
lr2 = LinearRegression().fit(X2, y)
r2_z500_u10 = lr2.score(X2, y)

# Model 3: Z500 + propagation_depth
X3 = np.column_stack([z500, df_events['propagation_depth_hpa'].fillna(0).values])
lr3 = LinearRegression().fit(X3, y)
r2_z500_prop = lr3.score(X3, y)

# Model 4: Z500 + cumulative_t10
X4 = np.column_stack([z500, df_events['cumulative_t10'].fillna(0).values])
lr4 = LinearRegression().fit(X4, y)
r2_z500_cumt10 = lr4.score(X4, y)

# Model 5: Z500 + u10_decel + propagation_depth
X5 = np.column_stack([z500, df_events['u10_decel'].fillna(0).values,
                       df_events['propagation_depth_hpa'].fillna(0).values])
lr5 = LinearRegression().fit(X5, y)
r2_full = lr5.score(X5, y)

print(f"  Z500 alone:                    R2 = {r2_z500:.4f}")
print(f"  Z500 + u10 deceleration:       R2 = {r2_z500_u10:.4f}  (Delta = {r2_z500_u10-r2_z500:+.4f})")
print(f"  Z500 + propagation depth:      R2 = {r2_z500_prop:.4f}  (Delta = {r2_z500_prop-r2_z500:+.4f})")
print(f"  Z500 + cumulative T10:         R2 = {r2_z500_cumt10:.4f}  (Delta = {r2_z500_cumt10-r2_z500:+.4f})")
print(f"  Z500 + u10_decel + prop_depth: R2 = {r2_full:.4f}  (Delta = {r2_full-r2_z500:+.4f})")

# LOO cross-validated R2
loo = LeaveOneOut()
for name, X in [('Z500', z500.reshape(-1, 1)),
                ('Z500+u10_decel', X2),
                ('Z500+u10_decel+prop', X5)]:
    y_pred = np.zeros(len(y))
    for train_idx, test_idx in loo.split(X):
        lr = LinearRegression().fit(X[train_idx], y[train_idx])
        y_pred[test_idx] = lr.predict(X[test_idx])
    loo_r2 = 1 - np.sum((y - y_pred)**2) / np.sum((y - y.mean())**2)
    print(f"  LOO R2 ({name}): {loo_r2:.4f}")

# ── 6. Aura MLS independent validation ────────────────────────────
print("\n" + "=" * 70)
print("6. AURA MLS INDEPENDENT VALIDATION (2004-2019)")
print("=" * 70)

# MLS temperature columns in the panel
mls_t_cols = [c for c in panel.columns if c.startswith('mls_t_')]
print(f"MLS temperature columns: {mls_t_cols}")

mls_events = []
for i, onset in enumerate(ssw_dates):
    year = onset.year
    if year < 2004:
        continue

    pre = panel.loc[(panel.index >= onset - pd.Timedelta(days=30)) &
                    (panel.index < onset - pd.Timedelta(days=5))]
    post = panel.loc[(panel.index > onset + pd.Timedelta(days=5)) &
                     (panel.index <= onset + pd.Timedelta(days=30))]

    if len(pre) > 0 and len(post) > 0:
        m = {'date': str(onset.date()), 'log_rr': log_rr[i]}
        for col in mls_t_cols:
            pre_mean = pre[col].dropna().mean()
            post_mean = post[col].dropna().mean()
            if not np.isnan(pre_mean) and not np.isnan(post_mean):
                m[f'{col}_change'] = post_mean - pre_mean
        mls_events.append(m)

print(f"Events with MLS data: {len(mls_events)}")

if mls_events:
    mls_df = pd.DataFrame(mls_events)
    # Cross-validate: MLS T at 10hPa vs NCEP T at 10hPa
    ncep_vals = []
    mls_vals = []
    for _, row in mls_df.iterrows():
        onset = pd.Timestamp(row['date'])
        idx = list(ssw_dates).index(onset)
        ncep_t10_change = df_events.loc[idx, 't10_anom_post'] if idx < len(df_events) else np.nan
        mls_t10_change = row.get('mls_t_lev_10p0hpa_change', np.nan)
        if not np.isnan(ncep_t10_change) and not np.isnan(mls_t10_change):
            ncep_vals.append(ncep_t10_change)
            mls_vals.append(mls_t10_change)

    if len(ncep_vals) >= 3:
        r_val, p_val = stats.pearsonr(ncep_vals, mls_vals)
        print(f"\nNCEP vs MLS temperature change at 10 hPa:")
        print(f"  Pearson r = {r_val:.3f}, P = {p_val:.4f}, n = {len(ncep_vals)}")
        print(f"  -> {'STRONG' if r_val > 0.7 else 'Moderate' if r_val > 0.4 else 'Weak'} independent validation")

    # MLS T change vs avalanche response
    for col in mls_t_cols:
        change_col = f'{col}_change'
        if change_col in mls_df.columns:
            vals = mls_df[change_col].values
            lr_vals = mls_df['log_rr'].values
            mask = ~np.isnan(vals) & ~np.isnan(lr_vals)
            if mask.sum() >= 5:
                r, p = stats.spearmanr(vals[mask], lr_vals[mask])
                if p < 0.15:
                    print(f"  {col} change vs log(RR): rho={r:+.3f}, P={p:.4f}")

# ── 7. E-P Flux Proxy ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("7. EDDY HEAT FLUX PROXY (using zonal wind deceleration)")
print("=" * 70)

# Since we don't have v'T' directly, we use the fundamental relationship:
# Wave-driven vortex deceleration (du/dt) is proportional to EP flux divergence
# (Charney & Drazin 1961, Andrews et al. 1987)
#
# Specifically: du/dt = (1/a cos(φ)) * ∂F/∂z + residual
# So the zonal wind deceleration IS a measure of the wave forcing

print("Using vortex deceleration as EP flux divergence proxy")
print("(du/dt at 10 hPa ~ ∇·F by the transformed Eulerian mean equation)")
print()

# Compute multi-level deceleration as integrated wave forcing
for plev in pressure_levels:
    u_col = f'ncep_u_{plev}hpa'
    decel_vals = []
    for onset in ssw_dates:
        pre = panel.loc[(panel.index >= onset - pd.Timedelta(days=20)) &
                        (panel.index < onset)]
        post = panel.loc[(panel.index > onset) &
                         (panel.index <= onset + pd.Timedelta(days=20))]
        if len(pre) > 0 and len(post) > 0:
            decel_vals.append(pre[u_col].mean() - post[u_col].mean())
        else:
            decel_vals.append(np.nan)

    decel_arr = np.array(decel_vals)
    mask = ~np.isnan(decel_arr) & ~np.isnan(log_rr)
    if mask.sum() >= 5:
        r_val, p_val = stats.spearmanr(decel_arr[mask], log_rr[mask])
        print(f"  {plev:4d} hPa wind deceleration vs log(RR): rho={r_val:+.3f}, P={p_val:.4f}")

# Integrated deceleration across levels (proxy for column-integrated EP flux)
integrated_decel = np.zeros(len(ssw_dates))
for p_idx, plev2 in enumerate(pressure_levels):
    u_col = f'ncep_u_{plev2}hpa'
    for i, onset in enumerate(ssw_dates):
        pre = panel.loc[(panel.index >= onset - pd.Timedelta(days=20)) &
                        (panel.index < onset)]
        post = panel.loc[(panel.index > onset) &
                         (panel.index <= onset + pd.Timedelta(days=20))]
        if len(pre) > 0 and len(post) > 0:
            integrated_decel[i] += pre[u_col].mean() - post[u_col].mean()

r_int, p_int = stats.spearmanr(integrated_decel, log_rr)
print(f"\n  Column-integrated deceleration vs log(RR): rho={r_int:+.3f}, P={p_int:.4f}")

# ── 8. Conditional independence test ──────────────────────────────
print("\n" + "=" * 70)
print("8. DOES SSW ADD VALUE BEYOND TROPOSPHERIC STATE?")
print("=" * 70)

# Partial correlation: log(RR) ~ u10_decel | Z500
from functools import partial

def partial_corr(x, y, z):
    """Partial Spearman correlation of x and y controlling for z."""
    n = len(x)
    # Rank-transform
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    # Residualize
    resid_x = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    resid_y = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r, p = stats.pearsonr(resid_x, resid_y)
    return r, p

# Test: does stratospheric deceleration predict avalanche response
# AFTER controlling for Z500?
mask = ~np.isnan(df_events['u10_decel'].values)
u10_d = df_events['u10_decel'].values[mask]
z500_v = z500[mask]
lr_v = log_rr[mask]

r_partial, p_partial = partial_corr(u10_d, lr_v, z500_v)
print(f"Partial correlation (u10_decel -> log(RR) | Z500):")
print(f"  rho_partial = {r_partial:+.3f}, P = {p_partial:.4f}")

# Also test Z500 controlling for u10_decel
r_z500_partial, p_z500_partial = partial_corr(z500_v, lr_v, u10_d)
print(f"\nPartial correlation (Z500 -> log(RR) | u10_decel):")
print(f"  rho_partial = {r_z500_partial:+.3f}, P = {p_z500_partial:.4f}")

# ── 9. Pre-conditioning analysis ──────────────────────────────────
print("\n" + "=" * 70)
print("9. PRE-SSW TROPOSPHERIC CONDITIONING")
print("=" * 70)

# Test: does the tropospheric state BEFORE the SSW predict the response?
# This tests the "common cause" hypothesis
for lag_start, lag_end, label in [(-30, -10, "Pre-SSW (-30 to -10d)"),
                                    (-10, 0, "Immediate pre (-10 to 0d)")]:
    pre_z500 = []
    pre_t2m = []
    for onset in ssw_dates:
        window = panel.loc[(panel.index >= onset + pd.Timedelta(days=lag_start)) &
                           (panel.index < onset + pd.Timedelta(days=lag_end))]
        if len(window) > 0 and 'ncep_z500_nh' in window.columns:
            pre_z500.append(window['ncep_z500_nh'].mean())
        else:
            pre_z500.append(np.nan)

    pre_z500 = np.array(pre_z500)
    mask = ~np.isnan(pre_z500)
    if mask.sum() >= 5:
        r, p = stats.spearmanr(pre_z500[mask], log_rr[mask])
        print(f"  {label} Z500 vs log(RR): rho={r:+.3f}, P={p:.4f}")

# ── SAVE RESULTS ──────────────────────────────────────────────────
results = {
    'n_events': len(df_events),
    'per_event_metrics': df_events.to_dict(orient='records'),
    'correlations': correlation_results,
    'group_comparison': {
        'propagating_n': int(propagating.sum()),
        'non_propagating_n': int((~propagating).sum()),
        'propagating_mean_logRR': round(float(prop_rr.mean()), 4),
        'non_propagating_mean_logRR': round(float(nonprop_rr.mean()), 4),
        'propagating_RR': round(float(np.exp(prop_rr.mean())), 4),
        'non_propagating_RR': round(float(np.exp(nonprop_rr.mean())), 4),
        'mann_whitney_p': round(float(group_mw_p), 4) if not np.isnan(group_mw_p) else None,
        'cohens_d': round(float(cohens_d), 4),
    },
    'multi_predictor': {
        'r2_z500': round(r2_z500, 4),
        'r2_z500_u10': round(r2_z500_u10, 4),
        'r2_z500_prop': round(r2_z500_prop, 4),
        'r2_z500_cumt10': round(r2_z500_cumt10, 4),
        'r2_full': round(r2_full, 4),
    },
    'partial_correlations': {
        'u10_decel_given_z500': {'rho': round(r_partial, 4), 'p': round(p_partial, 4)},
        'z500_given_u10_decel': {'rho': round(r_z500_partial, 4), 'p': round(p_z500_partial, 4)},
    },
    'composite_propagation': {
        'lags': list(lags),
        'levels_hPa': pressure_levels,
    },
    'ep_flux_proxy': {
        'integrated_decel_rho': round(r_int, 4),
        'integrated_decel_p': round(p_int, 4),
    },
}

with open(os.path.join(RESULTS, 'downward_propagation.json'), 'w') as f:
    json.dump(results, f, indent=2, default=str)

print("\n" + "=" * 70)
print("RESULTS SAVED to data/results/downward_propagation.json")
print("=" * 70)

