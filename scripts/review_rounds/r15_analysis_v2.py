"""
R15 Comprehensive Analysis v2 — Fixed column names and parsers
==============================================================
"""
import pandas as pd
import numpy as np
from scipy import stats
from datetime import timedelta
import json, warnings, os
warnings.filterwarnings('ignore')

results = {}

# =============================================================================
# PHASE 1: Load ALL SSW catalogs
# =============================================================================
print("=" * 80)
print("PHASE 1: SSW Catalogs")
print("=" * 80)

ssw_butler = pd.read_csv('data/processed/atmospheric/butler_ssw_compendium_era5.csv', parse_dates=['date'])
print("Butler SSW: %d events (1979-2023)" % len(ssw_butler))

ssw_existing = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_existing.index = ssw_existing.index.tz_localize(None)
print("Existing SSW catalog: %d events" % len(ssw_existing))
print("Existing SSW dates:\n%s" % ssw_existing.index.tolist()[:20])

# Check overlap with European Alps range (2011-12 to 2015-04)
alps_range = pd.Timestamp('2011-12-14'), pd.Timestamp('2015-04-16')
ssw_in_alps_existing = ssw_existing[(ssw_existing.index >= alps_range[0]) & (ssw_existing.index <= alps_range[1])]
print("\nSSW events in Alps range (existing catalog): %d" % len(ssw_in_alps_existing))
for d in ssw_in_alps_existing.index:
    print("  %s" % d)

ssw_in_alps_butler = ssw_butler[(ssw_butler['date'] >= alps_range[0]) & (ssw_butler['date'] <= alps_range[1])]
print("SSW events in Alps range (Butler): %d" % len(ssw_in_alps_butler))
for _, r in ssw_in_alps_butler.iterrows():
    print("  %s (%s)" % (r['date'], r['event_name']))

# Merge: use all unique SSW dates from both catalogs
all_ssw_dates = sorted(set(list(ssw_existing.index) + list(ssw_butler['date'])))
print("\nTotal unique SSW dates: %d" % len(all_ssw_dates))

# =============================================================================
# PHASE 2: Multi-Country Replication (European Alps)
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 2: Multi-Country Replication")
print("=" * 80)

alps = pd.read_csv('data/cryosphere/european_alps/data_dmax.csv', sep=';', parse_dates=['date'])
# Clean NaN countries
alps = alps.dropna(subset=['country'])
print("Alps data: %d rows, %s" % (len(alps), alps['country'].unique().tolist()))

# Use ALL SSW dates in range
ssw_alps_dates = [d for d in all_ssw_dates if alps_range[0] <= d <= alps_range[1]]
print("SSW events in Alps range (merged): %d" % len(ssw_alps_dates))
for d in ssw_alps_dates:
    print("  %s" % d)

# If only 1-2 events, also do a broader "matched control" approach
# For each SSW, compare 30d post with the SAME calendar window across all non-SSW years

country_results = {}
for country in ['AT', 'FR', 'IT', 'DE', 'CH']:
    cdata = alps[alps['country'] == country].copy()
    n_regions = cdata['warningRegion'].nunique()
    daily_mean = cdata.groupby('date')['dangerLevelMax'].mean()
    
    ssw_danger = []
    ctrl_danger = []
    event_details = []
    
    for ssw_date in ssw_alps_dates:
        ssw_date = pd.Timestamp(ssw_date)
        post = daily_mean[(daily_mean.index >= ssw_date) & (daily_mean.index < ssw_date + timedelta(days=30))]
        pre = daily_mean[(daily_mean.index >= ssw_date - timedelta(days=30)) & (daily_mean.index < ssw_date)]
        
        if len(post) > 5 and len(pre) > 5:
            ssw_danger.append(post.mean())
            ctrl_danger.append(pre.mean())
            event_details.append({
                'event': str(ssw_date.date()),
                'pre_mean': float(pre.mean()),
                'post_mean': float(post.mean()),
                'change_pct': float((post.mean() - pre.mean()) / pre.mean() * 100)
            })
    
    if len(ssw_danger) >= 1:
        country_results[country] = {
            'n_regions': n_regions,
            'n_events': len(ssw_danger),
            'events': event_details
        }
        
        print("\n%s: %d regions, %d SSW events" % (country, n_regions, len(ssw_danger)))
        for ed in event_details:
            print("  %s: pre=%.3f, post=%.3f (%+.1f%%)" % (ed['event'], ed['pre_mean'], ed['post_mean'], ed['change_pct']))

# Also: use calendar matching — for each SSW event, compare with same dates in other years
print("\n--- Calendar-Matched Control Analysis ---")
for country in ['AT', 'FR', 'IT', 'DE', 'CH']:
    cdata = alps[alps['country'] == country].copy()
    daily_mean = cdata.groupby('date')['dangerLevelMax'].mean()
    
    for ssw_date in ssw_alps_dates:
        ssw_date = pd.Timestamp(ssw_date)
        # SSW year window
        post_ssw = daily_mean[(daily_mean.index >= ssw_date) & 
                               (daily_mean.index < ssw_date + timedelta(days=30))]
        
        # Same calendar window in non-SSW years
        control_means = []
        for year_offset in [-1, 1, 2]:  # Try adjacent years
            ctrl_date = ssw_date + pd.DateOffset(years=year_offset)
            ctrl_vals = daily_mean[(daily_mean.index >= ctrl_date) & 
                                    (daily_mean.index < ctrl_date + timedelta(days=30))]
            if len(ctrl_vals) > 10:
                control_means.append(ctrl_vals.mean())
        
        if len(post_ssw) > 5 and len(control_means) >= 1:
            ctrl_mean = np.mean(control_means)
            ssw_mean = post_ssw.mean()
            change = (ssw_mean - ctrl_mean) / ctrl_mean * 100
            print("  %s %s: SSW=%.3f, calendar-ctrl=%.3f (%+.1f%%)" % (
                country, ssw_date.date(), ssw_mean, ctrl_mean, change))

results['european_alps_replication'] = country_results

# =============================================================================
# PHASE 3: SNOWPACK Stability Index Analysis
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 3: SNOWPACK Stability Index Analysis")
print("=" * 80)

sp = pd.read_csv('data/cryosphere/swiss_snowpack/data_rf2_tidy.csv', 
                  parse_dates=['datum'], low_memory=False)
print("SNOWPACK data: %d rows × %d cols" % sp.shape)
print("Date range: %s to %s" % (sp['datum'].min(), sp['datum'].max()))
print("Stations: %d" % sp['station_code'].nunique())

# Key stability columns
stability_cols = ['ssi_pwl', 'sk38_pwl', 'sn38_pwl', 'ccl_pwl', 
                  'ssi_pwl_100', 'sk38_pwl_100', 'sn38_pwl_100', 'ccl_pwl_100',
                  'pwl_100', 'pwl_100_15', 'base_pwl']
# Key snow property columns
snow_cols = ['HS_mod', 'HS_meas', 'HN24', 'SWE', 'hoar_size', 'Pen_depth', 'min_ccl_pen']
# Temperature / radiation columns
temp_cols = ['TA', 'TSS_mod', 'TS0', 'TS1', 'TS2', 'ISWR', 'ILWR', 'LWR_net']
# Danger level
has_danger = 'dangerLevel' in sp.columns

# SSW events in snowpack range
ssw_in_sp = [d for d in all_ssw_dates if sp['datum'].min() <= d <= sp['datum'].max()]
# Filter to winter only
ssw_in_sp_winter = [d for d in ssw_in_sp if d.month in [11, 12, 1, 2, 3]]
print("SSW events in SNOWPACK range (winter): %d" % len(ssw_in_sp_winter))

# Compute daily station-averaged values
all_test_cols = stability_cols + snow_cols + temp_cols + (['dangerLevel'] if has_danger else [])
existing_cols = [c for c in all_test_cols if c in sp.columns]
print("Testing %d variables" % len(existing_cols))

daily_sp = sp.groupby('datum')[existing_cols].mean()

# SSW composite analysis
print("\n--- SSW Composite Response (30-day windows) ---")
sp_results = {}

for col in existing_cols:
    series = daily_sp[col].dropna()
    if len(series) < 100:
        continue
    
    ssw_vals = []
    ctrl_vals = []
    event_details = []
    
    for ssw_date in ssw_in_sp_winter:
        ssw_date = pd.Timestamp(ssw_date)
        post = series[(series.index >= ssw_date) & (series.index < ssw_date + timedelta(days=30))]
        pre = series[(series.index >= ssw_date - timedelta(days=30)) & (series.index < ssw_date)]
        
        if len(post) > 10 and len(pre) > 10:
            ssw_vals.append(post.mean())
            ctrl_vals.append(pre.mean())
            event_details.append({
                'date': str(ssw_date.date()),
                'pre': float(pre.mean()),
                'post': float(post.mean())
            })
    
    if len(ssw_vals) >= 5:
        ssw_arr = np.array(ssw_vals)
        ctrl_arr = np.array(ctrl_vals)
        diff = ssw_arr - ctrl_arr
        
        _, p_wilcox = stats.wilcoxon(ctrl_arr, ssw_arr)
        _, p_ttest = stats.ttest_rel(ctrl_arr, ssw_arr)
        
        pct_change = (ssw_arr.mean() - ctrl_arr.mean()) / abs(ctrl_arr.mean()) * 100 if ctrl_arr.mean() != 0 else 0
        n_increase = int(np.sum(diff > 0))
        n_decrease = int(np.sum(diff < 0))
        cohens_d = diff.mean() / diff.std() if diff.std() > 0 else 0
        
        sp_results[col] = {
            'n_events': len(ssw_vals),
            'pre_mean': float(ctrl_arr.mean()),
            'post_mean': float(ssw_arr.mean()),
            'pct_change': float(pct_change),
            'n_increase': n_increase,
            'n_decrease': n_decrease,
            'p_wilcoxon': float(p_wilcox),
            'p_ttest': float(p_ttest),
            'cohens_d': float(cohens_d)
        }
        
        sig = "***" if min(p_wilcox, p_ttest) < 0.001 else "**" if min(p_wilcox, p_ttest) < 0.01 else "*" if min(p_wilcox, p_ttest) < 0.05 else "†" if min(p_wilcox, p_ttest) < 0.1 else ""
        direction = "↑" if pct_change > 0 else "↓"
        print("  %s %s: pre=%.3f, post=%.3f, %s%.1f%%, %d↑/%d↓, Pw=%.4f, Pt=%.4f, d=%.2f %s" % (
            direction, col, ctrl_arr.mean(), ssw_arr.mean(), direction, abs(pct_change), 
            n_increase, n_decrease, p_wilcox, p_ttest, cohens_d, sig))

results['snowpack_stability'] = sp_results

# Highlight key stability interpretation
print("\n--- KEY STABILITY FINDINGS ---")
for col in ['ssi_pwl', 'sk38_pwl', 'sn38_pwl', 'ccl_pwl']:
    if col in sp_results:
        r = sp_results[col]
        # For stability indices: HIGHER = MORE STABLE (less avalanche prone)
        direction = "MORE STABLE" if r['pct_change'] > 0 else "LESS STABLE"
        print("  %s: %+.1f%% change → %s after SSW (P_w=%.4f)" % (col, r['pct_change'], direction, r['p_wilcoxon']))

# =============================================================================
# PHASE 4: Davos Individual Avalanche Analysis  
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 4: Davos Avalanche Analysis")
print("=" * 80)

davos_obs = pd.read_csv('data/cryosphere/davos_avalanches/avalanche_observations.csv', 
                         sep=';', parse_dates=['date_release'])
davos_daily = pd.read_csv('data/cryosphere/davos_avalanches/daily_activity.csv',
                           sep=';', parse_dates=['date'])
print("Observations: %d | Daily: %d days" % (len(davos_obs), len(davos_daily)))

# Use daily AAI (Avalanche Activity Index) for dry natural — more robust than raw counts
aai_dry = 'AAI_all.dry' if 'AAI_all.dry' in davos_daily.columns else None
aai_dry_nat = None
for c in davos_daily.columns:
    if 'dry' in c.lower() and 'natural' in c.lower():
        aai_dry_nat = c
        break

print("AAI dry column: %s" % aai_dry)
print("AAI dry natural column: %s" % aai_dry_nat)

# Use AAI_all.dry as primary metric
target_col = aai_dry_nat if aai_dry_nat else aai_dry
if target_col is None:
    target_col = 'AAI_all.dry'

# SSW events in Davos range
davos_dates = davos_daily.set_index('date')
ssw_in_davos = [d for d in all_ssw_dates 
                if davos_dates.index.min() <= d <= davos_dates.index.max() 
                and d.month in [11, 12, 1, 2, 3]]
print("SSW events in Davos range (winter): %d" % len(ssw_in_davos))

# AAI-based superposed epoch
print("\n--- AAI-based SSW Response ---")
for col_name in ['AAI_all.dry', 'AAI_all', 'AAI_all.wet']:
    if col_name not in davos_dates.columns:
        continue
    
    series = davos_dates[col_name].dropna()
    ssw_vals = []
    ctrl_vals = []
    events = []
    
    for ssw_date in ssw_in_davos:
        ssw_date = pd.Timestamp(ssw_date)
        post = series[(series.index >= ssw_date) & (series.index < ssw_date + timedelta(days=30))]
        pre = series[(series.index >= ssw_date - timedelta(days=30)) & (series.index < ssw_date)]
        
        if len(post) > 10 and len(pre) > 10:
            ssw_vals.append(post.mean())
            ctrl_vals.append(pre.mean())
            events.append({'date': str(ssw_date.date()), 'pre': float(pre.mean()), 'post': float(post.mean())})
    
    if len(ssw_vals) >= 5:
        ssw_a = np.array(ssw_vals)
        ctrl_a = np.array(ctrl_vals)
        _, p_w = stats.wilcoxon(ctrl_a, ssw_a)
        _, p_t = stats.ttest_rel(ctrl_a, ssw_a)
        rr = ssw_a.mean() / ctrl_a.mean() if ctrl_a.mean() > 0 else 0
        n_dec = np.sum(ssw_a < ctrl_a)
        
        sig = "*" if min(p_w, p_t) < 0.05 else ""
        print("  %s: n=%d, pre=%.2f, post=%.2f, RR=%.3f, %d/%d↓, Pw=%.4f %s" % (
            col_name, len(ssw_vals), ctrl_a.mean(), ssw_a.mean(), rr, n_dec, len(ssw_vals), p_w, sig))
        
        results['davos_%s' % col_name.replace('.', '_')] = {
            'n_events': len(ssw_vals),
            'pre_mean': float(ctrl_a.mean()),
            'post_mean': float(ssw_a.mean()),
            'rate_ratio': float(rr),
            'n_decrease': int(n_dec),
            'p_wilcoxon': float(p_w),
            'p_ttest': float(p_t),
            'events': events
        }

# Size-class stratified analysis
print("\n--- Size-Class Stratified (observations) ---")
for size_min in [1, 2, 3]:
    subset = davos_obs[(davos_obs['snow_type'] == 'dry') & 
                        (davos_obs['trigger_type'] == 'NATURAL') &
                        (davos_obs['aval_size_class'] >= size_min)]
    daily = subset.groupby('date_release').size()
    full_range = pd.date_range(davos_obs['date_release'].min(), davos_obs['date_release'].max())
    daily = daily.reindex(full_range, fill_value=0)
    
    ssw_vals = []
    ctrl_vals = []
    for ssw_date in ssw_in_davos:
        ssw_date = pd.Timestamp(ssw_date)
        post = daily[(daily.index >= ssw_date) & (daily.index < ssw_date + timedelta(days=30))]
        pre = daily[(daily.index >= ssw_date - timedelta(days=30)) & (daily.index < ssw_date)]
        if len(post) > 10 and len(pre) > 10:
            ssw_vals.append(post.mean())
            ctrl_vals.append(pre.mean())
    
    if len(ssw_vals) >= 5:
        ssw_a = np.array(ssw_vals)
        ctrl_a = np.array(ctrl_vals)
        _, p_w = stats.wilcoxon(ctrl_a, ssw_a)
        rr = ssw_a.mean() / ctrl_a.mean() if ctrl_a.mean() > 0 else 0
        n_dec = np.sum(ssw_a < ctrl_a)
        print("  Size≥%d: n=%d, pre=%.3f, post=%.3f, RR=%.3f, %d/%d↓, P=%.4f" % (
            size_min, len(ssw_vals), ctrl_a.mean(), ssw_a.mean(), rr, n_dec, len(ssw_vals), p_w))

# =============================================================================
# PHASE 5: Swiss Panel — Proper Dry Avalanche Columns
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 5: Swiss Panel — dry_natural_size_1234")
print("=" * 80)

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
panel_dates = panel.index if isinstance(panel.index, pd.DatetimeIndex) else pd.to_datetime(panel.index)

for target in ['dry_natural_size_1234', 'aai_all_dry']:
    if target not in panel.columns:
        continue
    
    series = panel[target]
    ssw_in_panel = [d for d in all_ssw_dates 
                    if panel_dates.min() <= d <= panel_dates.max()
                    and d.month in [11, 12, 1, 2, 3]]
    
    ssw_vals = []
    ctrl_vals = []
    event_details = []
    
    for ssw_date in ssw_in_panel:
        ssw_date = pd.Timestamp(ssw_date)
        post_mask = (panel_dates >= ssw_date) & (panel_dates < ssw_date + timedelta(days=30))
        pre_mask = (panel_dates >= ssw_date - timedelta(days=30)) & (panel_dates < ssw_date)
        post_data = series[post_mask].dropna()
        pre_data = series[pre_mask].dropna()
        
        if len(post_data) > 10 and len(pre_data) > 10:
            ssw_vals.append(post_data.mean())
            ctrl_vals.append(pre_data.mean())
            event_details.append({
                'date': str(ssw_date.date()),
                'pre': float(pre_data.mean()),
                'post': float(post_data.mean()),
                'change_pct': float((post_data.mean() - pre_data.mean()) / pre_data.mean() * 100) if pre_data.mean() > 0 else 0
            })
    
    if len(ssw_vals) >= 5:
        ssw_a = np.array(ssw_vals)
        ctrl_a = np.array(ctrl_vals)
        diff = ssw_a - ctrl_a
        _, p_w = stats.wilcoxon(ctrl_a, ssw_a)
        _, p_t = stats.ttest_rel(ctrl_a, ssw_a)
        rr = ssw_a.mean() / ctrl_a.mean() if ctrl_a.mean() > 0 else 0
        n_dec = np.sum(ssw_a < ctrl_a)
        cohens_d = diff.mean() / diff.std() if diff.std() > 0 else 0
        
        # Sign test
        n_neg = int(np.sum(diff < 0))
        p_sign = stats.binomtest(n_neg, np.sum(diff != 0), 0.5).pvalue
        
        # Permutation test
        n_perm = 10000
        obs_diff = diff.mean()
        perm_diffs = np.zeros(n_perm)
        combined = np.concatenate([ctrl_a, ssw_a])
        n = len(ctrl_a)
        for i in range(n_perm):
            np.random.shuffle(combined)
            perm_diffs[i] = combined[n:].mean() - combined[:n].mean()
        p_perm = np.mean(np.abs(perm_diffs) >= np.abs(obs_diff))
        
        # Bootstrap CI
        n_boot = 10000
        boot_diffs = np.zeros(n_boot)
        for i in range(n_boot):
            idx = np.random.choice(len(diff), size=len(diff), replace=True)
            boot_diffs[i] = diff[idx].mean()
        ci_lo, ci_hi = np.percentile(boot_diffs, [2.5, 97.5])
        
        print("\n  %s (n=%d events):" % (target, len(ssw_vals)))
        print("    Pre-SSW: %.3f ± %.3f" % (ctrl_a.mean(), ctrl_a.std()))
        print("    Post-SSW: %.3f ± %.3f" % (ssw_a.mean(), ssw_a.std()))
        print("    Rate ratio: %.3f" % rr)
        print("    Events with decrease: %d/%d (%.0f%%)" % (n_dec, len(ssw_vals), n_dec/len(ssw_vals)*100))
        print("    Cohen's d: %.3f" % cohens_d)
        print("    Wilcoxon P: %.4f" % p_w)
        print("    Paired t P: %.4f" % p_t)
        print("    Sign test P: %.4f" % p_sign)
        print("    Permutation P: %.4f" % p_perm)
        print("    Bootstrap 95%% CI: [%.3f, %.3f]" % (ci_lo, ci_hi))
        
        results['panel_%s' % target] = {
            'n_events': len(ssw_vals),
            'pre_mean': float(ctrl_a.mean()),
            'post_mean': float(ssw_a.mean()),
            'rate_ratio': float(rr),
            'n_decrease': int(n_dec),
            'cohens_d': float(cohens_d),
            'p_wilcoxon': float(p_w),
            'p_ttest': float(p_t),
            'p_sign': float(p_sign),
            'p_permutation': float(p_perm),
            'bootstrap_ci': [float(ci_lo), float(ci_hi)],
            'events': event_details
        }
        
        for ed in event_details:
            d = "↓" if ed['change_pct'] < 0 else "↑"
            print("      %s: %.2f → %.2f (%s%.1f%%)" % (ed['date'], ed['pre'], ed['post'], d, abs(ed['change_pct'])))

# =============================================================================
# PHASE 6: ENSO/QBO Conditioning with Proper Columns
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 6: ENSO/QBO Conditioning")
print("=" * 80)

target = 'dry_natural_size_1234' if 'dry_natural_size_1234' in panel.columns else 'aai_all_dry'
series = panel[target]

for phase_var in ['enso_phase', 'qbo_phase']:
    print("\n--- %s ---" % phase_var)
    for phase in ssw_butler[phase_var].unique():
        phase_ssws = ssw_butler[ssw_butler[phase_var] == phase]
        phase_dates = [d for d in phase_ssws['date'] 
                       if panel_dates.min() <= d <= panel_dates.max() and d.month in [11,12,1,2,3]]
        
        if len(phase_dates) < 2:
            continue
        
        ssw_v, ctrl_v = [], []
        for d in phase_dates:
            d = pd.Timestamp(d)
            post = series[(panel_dates >= d) & (panel_dates < d + timedelta(days=30))].dropna()
            pre = series[(panel_dates >= d - timedelta(days=30)) & (panel_dates < d)].dropna()
            if len(post) > 10 and len(pre) > 10:
                ssw_v.append(post.mean())
                ctrl_v.append(pre.mean())
        
        if len(ssw_v) >= 2:
            ssw_a = np.array(ssw_v)
            ctrl_a = np.array(ctrl_v)
            n_dec = np.sum(ssw_a < ctrl_a)
            pct = (ssw_a.mean() - ctrl_a.mean()) / ctrl_a.mean() * 100 if ctrl_a.mean() > 0 else 0
            p_val = stats.wilcoxon(ctrl_a, ssw_a)[1] if len(ssw_v) >= 6 else stats.ttest_rel(ctrl_a, ssw_a)[1]
            print("  %s=%s (n=%d): pre=%.2f, post=%.2f, %+.1f%%, %d/%d↓, P=%.4f" % (
                phase_var, phase, len(ssw_v), ctrl_a.mean(), ssw_a.mean(), pct, n_dec, len(ssw_v), p_val))

# =============================================================================
# PHASE 7: SSW Type (Displacement vs Split) with Proper Columns
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 7: SSW Type Stratification")
print("=" * 80)

split_events = ['FEB 2018', 'JAN 2009', 'FEB 1979', 'JAN 1985', 'FEB 1989']
ssw_typed = ssw_butler.copy()
ssw_typed['ssw_type'] = ssw_typed['event_name'].apply(lambda x: 'split' if x in split_events else 'displacement')

for ssw_type in ['displacement', 'split']:
    typed = ssw_typed[ssw_typed['ssw_type'] == ssw_type]
    typed_dates = [d for d in typed['date'] 
                   if panel_dates.min() <= d <= panel_dates.max() and d.month in [11,12,1,2,3]]
    
    if len(typed_dates) < 2:
        print("  %s: too few events (%d)" % (ssw_type, len(typed_dates)))
        continue
    
    ssw_v, ctrl_v = [], []
    for d in typed_dates:
        d = pd.Timestamp(d)
        post = series[(panel_dates >= d) & (panel_dates < d + timedelta(days=30))].dropna()
        pre = series[(panel_dates >= d - timedelta(days=30)) & (panel_dates < d)].dropna()
        if len(post) > 10 and len(pre) > 10:
            ssw_v.append(post.mean())
            ctrl_v.append(pre.mean())
    
    if len(ssw_v) >= 2:
        ssw_a = np.array(ssw_v)
        ctrl_a = np.array(ctrl_v)
        n_dec = np.sum(ssw_a < ctrl_a)
        pct = (ssw_a.mean() - ctrl_a.mean()) / ctrl_a.mean() * 100 if ctrl_a.mean() > 0 else 0
        print("  %s (n=%d): pre=%.2f, post=%.2f, %+.1f%%, %d/%d↓" % (
            ssw_type, len(ssw_v), ctrl_a.mean(), ssw_a.mean(), pct, n_dec, len(ssw_v)))
        
        results['ssw_type_%s' % ssw_type] = {
            'n_events': len(ssw_v),
            'pre_mean': float(ctrl_a.mean()),
            'post_mean': float(ssw_a.mean()),
            'pct_change': float(pct),
            'n_decrease': int(n_dec)
        }

# =============================================================================
# PHASE 8: Swiss Danger Descriptions (fixed parsing)
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 8: Swiss Danger Descriptions (dry vs wet)")
print("=" * 80)

try:
    desc = pd.read_csv('data/cryosphere/swiss_snowpack/danger_descriptions_2012_2020.csv',
                        encoding='latin-1')  # comma delimiter, not semicolon
    print("Loaded: %d rows, columns: %s" % (len(desc), desc.columns.tolist()))
    desc['date'] = pd.to_datetime(desc['validFromDate'])
    
    ssw_in_desc = [d for d in all_ssw_dates 
                   if desc['date'].min() <= d <= desc['date'].max() and d.month in [11,12,1,2,3]]
    print("SSW events in range: %d" % len(ssw_in_desc))
    
    for problem_type in ['dry', 'wet']:
        subset = desc[desc['problem'] == problem_type].copy()
        daily_dl = subset.groupby('date')['dangerlevel'].mean()
        
        ssw_v, ctrl_v = [], []
        for d in ssw_in_desc:
            d = pd.Timestamp(d)
            post = daily_dl[(daily_dl.index >= d) & (daily_dl.index < d + timedelta(days=30))]
            pre = daily_dl[(daily_dl.index >= d - timedelta(days=30)) & (daily_dl.index < d)]
            if len(post) > 5 and len(pre) > 5:
                ssw_v.append(post.mean())
                ctrl_v.append(pre.mean())
        
        if len(ssw_v) >= 2:
            ssw_a = np.array(ssw_v)
            ctrl_a = np.array(ctrl_v)
            n_dec = np.sum(ssw_a < ctrl_a)
            pct = (ssw_a.mean() - ctrl_a.mean()) / ctrl_a.mean() * 100 if ctrl_a.mean() > 0 else 0
            p = stats.wilcoxon(ctrl_a, ssw_a)[1] if len(ssw_v) >= 6 else stats.ttest_rel(ctrl_a, ssw_a)[1]
            print("  %s problem (n=%d): pre=%.2f, post=%.2f, %+.1f%%, %d/%d↓, P=%.4f" % (
                problem_type, len(ssw_v), ctrl_a.mean(), ssw_a.mean(), pct, n_dec, len(ssw_v), p))
            
            results['danger_desc_%s' % problem_type] = {
                'n_events': len(ssw_v),
                'pre_mean': float(ctrl_a.mean()),
                'post_mean': float(ssw_a.mean()),
                'pct_change': float(pct),
                'n_decrease': int(n_dec),
                'p_value': float(p)
            }
except Exception as e:
    print("Error: %s" % str(e))
    import traceback; traceback.print_exc()

# =============================================================================
# PHASE 9: ERA5 Stratospheric Analysis (existing processed data)
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 9: ERA5 Stratospheric Response")
print("=" * 80)

try:
    era5_strat = pd.read_parquet('data/processed/atmospheric/era5_polar_strat_means.parquet')
    print("ERA5 strat means: %d rows × %d cols" % era5_strat.shape)
    print("Columns: %s" % era5_strat.columns.tolist())
    
    era5_dates = era5_strat.index if isinstance(era5_strat.index, pd.DatetimeIndex) else pd.to_datetime(era5_strat.index)
    if era5_dates.tz is not None:
        era5_dates = era5_dates.tz_localize(None)
    
    ssw_in_era5 = [d for d in all_ssw_dates if era5_dates.min() <= d <= era5_dates.max()]
    print("SSW events in ERA5 range: %d" % len(ssw_in_era5))
    
    for col in era5_strat.columns[:10]:
        series_e = era5_strat[col].dropna()
        if isinstance(series_e.index, pd.DatetimeIndex) and series_e.index.tz is not None:
            series_e.index = series_e.index.tz_localize(None)
        
        ssw_v, ctrl_v = [], []
        for d in ssw_in_era5:
            d = pd.Timestamp(d)
            post = series_e[(series_e.index >= d) & (series_e.index < d + timedelta(days=30))]
            pre = series_e[(series_e.index >= d - timedelta(days=30)) & (series_e.index < d)]
            if len(post) > 5 and len(pre) > 5:
                ssw_v.append(post.mean())
                ctrl_v.append(pre.mean())
        
        if len(ssw_v) >= 5:
            ssw_a = np.array(ssw_v)
            ctrl_a = np.array(ctrl_v)
            _, p_w = stats.wilcoxon(ctrl_a, ssw_a)
            pct = (ssw_a.mean() - ctrl_a.mean()) / abs(ctrl_a.mean()) * 100 if ctrl_a.mean() != 0 else 0
            sig = "***" if p_w < 0.001 else "**" if p_w < 0.01 else "*" if p_w < 0.05 else ""
            print("  %s: %+.2f%%, P=%.4f %s" % (col, pct, p_w, sig))
except Exception as e:
    print("Error: %s" % str(e))

# =============================================================================
# PHASE 10: Norwegian Replication (existing data)
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 10: Norwegian Replication Check")
print("=" * 80)

try:
    if 'norway_aval_count' in panel.columns:
        norway = panel['norway_aval_count'].dropna()
        print("Norway data: %d non-null days" % len(norway))
        
        ssw_in_norway = [d for d in all_ssw_dates 
                         if panel_dates.min() <= d <= panel_dates.max() and d.month in [11,12,1,2,3]]
        
        ssw_v, ctrl_v = [], []
        for d in ssw_in_norway:
            d = pd.Timestamp(d)
            post = norway[(panel_dates >= d) & (panel_dates < d + timedelta(days=30))]
            pre = norway[(panel_dates >= d - timedelta(days=30)) & (panel_dates < d)]
            if len(post) > 5 and len(pre) > 5:
                ssw_v.append(post.mean())
                ctrl_v.append(pre.mean())
        
        if len(ssw_v) >= 2:
            ssw_a = np.array(ssw_v)
            ctrl_a = np.array(ctrl_v)
            n_dec = np.sum(ssw_a < ctrl_a)
            print("  Norway: n=%d, pre=%.2f, post=%.2f, %d/%d↓" % (
                len(ssw_v), ctrl_a.mean(), ssw_a.mean(), n_dec, len(ssw_v)))
except Exception as e:
    print("Error: %s" % str(e))

# =============================================================================
# PHASE 11: Comprehensive Evidence Summary
# =============================================================================
print("\n" + "=" * 80)
print("COMPREHENSIVE EVIDENCE SUMMARY")
print("=" * 80)

print("\n1. MULTI-COUNTRY (European Alps 2011-2015)")
for c, r in sorted(results.get('european_alps_replication', {}).items()):
    print("   %s: %d regions, %d events — %s" % (c, r['n_regions'], r['n_events'],
          ', '.join(['%+.1f%%' % e['change_pct'] for e in r['events']])))

print("\n2. SNOWPACK STABILITY (%d variables)" % len(results.get('snowpack_stability', {})))
for col in ['ssi_pwl', 'sk38_pwl', 'sn38_pwl', 'ccl_pwl', 'HS_mod', 'HN24', 'TA', 'dangerLevel']:
    if col in results.get('snowpack_stability', {}):
        r = results['snowpack_stability'][col]
        sig = "✓" if r['p_wilcoxon'] < 0.05 else "○" if r['p_wilcoxon'] < 0.1 else ""
        print("   %s %s: %+.1f%%, d=%.2f, P=%.4f" % (sig, col, r['pct_change'], r['cohens_d'], r['p_wilcoxon']))

print("\n3. DAVOS INDIVIDUAL AVALANCHES")
for key in sorted([k for k in results if k.startswith('davos_')]):
    r = results[key]
    print("   %s: RR=%.3f, %d/%d↓, P=%.4f" % (key, r['rate_ratio'], r['n_decrease'], r['n_events'], r['p_wilcoxon']))

print("\n4. SWISS PANEL (dry_natural_size_1234)")
if 'panel_dry_natural_size_1234' in results:
    r = results['panel_dry_natural_size_1234']
    print("   RR=%.3f, %d/%d↓, Pw=%.4f, Pt=%.4f, Pperm=%.4f" % (
        r['rate_ratio'], r['n_decrease'], r['n_events'], r['p_wilcoxon'], r['p_ttest'], r['p_permutation']))
    print("   Bootstrap CI: [%.3f, %.3f]" % tuple(r['bootstrap_ci']))
    print("   Cohen's d: %.3f" % r['cohens_d'])

print("\n5. SSW TYPE STRATIFICATION")
for t in ['displacement', 'split']:
    key = 'ssw_type_%s' % t
    if key in results:
        r = results[key]
        print("   %s (n=%d): %+.1f%%, %d/%d↓" % (t.title(), r['n_events'], r['pct_change'], r['n_decrease'], r['n_events']))

print("\n6. DANGER DESCRIPTIONS (dry vs wet)")
for t in ['dry', 'wet']:
    key = 'danger_desc_%s' % t
    if key in results:
        r = results[key]
        print("   %s: %+.1f%%, %d/%d↓, P=%.4f" % (t, r['pct_change'], r['n_decrease'], r['n_events'], r['p_value']))

# Save results
def make_serializable(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.Timestamp):
        return obj.isoformat()
    return obj

def deep_convert(d):
    if isinstance(d, dict):
        return {k: deep_convert(v) for k, v in d.items()}
    elif isinstance(d, list):
        return [deep_convert(i) for i in d]
    return make_serializable(d)

with open('data/results/r15_comprehensive_v2.json', 'w') as f:
    json.dump(deep_convert(results), f, indent=2, default=str)
print("\nResults saved to data/results/r15_comprehensive_v2.json")
