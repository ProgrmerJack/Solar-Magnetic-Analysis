"""
Comprehensive Norwegian Varsom SSW analysis with 10 regions × 4 SSW events.
Proper statistical testing: Mann-Whitney, permutation, event-level concordance.
Also: DOY-matched controls, per-region analysis, latitudinal gradient.
"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats
from datetime import datetime, timedelta

# Load data
with open('data/cryosphere/norway_expanded/varsom_targeted.json') as f:
    records = json.load(f)

df = pd.DataFrame(records)
df['date'] = pd.to_datetime(df['valid_from']).dt.date
df = df[df['danger_level'] > 0].copy()
df['doy'] = pd.to_datetime(df['valid_from']).dt.dayofyear

print(f"Total valid records: {len(df)}")
print(f"Regions: {df['region_name'].nunique()}: {sorted(df['region_name'].unique())}")
print(f"Date range: {df['date'].min()} to {df['date'].max()}")

# SSW events in our data
ssw_dates = {
    '2018-02-12': datetime(2018, 2, 12).date(),
    '2019-01-02': datetime(2019, 1, 2).date(),
    '2021-01-05': datetime(2021, 1, 5).date(),
    '2023-02-16': datetime(2023, 2, 16).date(),
}

# Separate SSW and control data
ssw_df = df[df['ssw_date'].notna()].copy()
ctrl_df = df[df['ssw_date'].isna()].copy()

# ============================================================
# 1. SSW ±15d WINDOW vs DOY-MATCHED CONTROL
# ============================================================
print("\n" + "="*60)
print("1. SSW ±15d WINDOW ANALYSIS")
print("="*60)

ssw_window_records = []
for _, row in ssw_df.iterrows():
    if row['ssw_date']:
        ssw_dt = ssw_dates.get(row['ssw_date'])
        if ssw_dt:
            rec_dt = row['date']
            delta = (rec_dt - ssw_dt).days
            if -15 <= delta <= 15:
                ssw_window_records.append(row)

ssw_win = pd.DataFrame(ssw_window_records)
print(f"SSW ±15d records: {len(ssw_win)}")

# DOY-matched control
ssw_doys = set(ssw_win['doy'].unique())
ctrl_matched = ctrl_df[ctrl_df['doy'].isin(ssw_doys)]
print(f"DOY-matched control records: {len(ctrl_matched)}")

# Mann-Whitney U test
u_stat, mw_p = stats.mannwhitneyu(ssw_win['danger_level'], ctrl_matched['danger_level'], alternative='less')
print(f"\nSSW mean danger: {ssw_win['danger_level'].mean():.3f} (n={len(ssw_win)})")
print(f"Control mean danger: {ctrl_matched['danger_level'].mean():.3f} (n={len(ctrl_matched)})")
print(f"Mann-Whitney U P (one-sided): {mw_p:.6f}")

# Cohen's d
d_ssw = ssw_win['danger_level'].values
d_ctrl = ctrl_matched['danger_level'].values
pooled_std = np.sqrt((d_ssw.var() + d_ctrl.var()) / 2)
cohens_d = (d_ssw.mean() - d_ctrl.mean()) / pooled_std
print(f"Cohen's d: {cohens_d:.3f}")

# Permutation test
n_perm = 10000
obs_diff = d_ssw.mean() - d_ctrl.mean()
combined = np.concatenate([d_ssw, d_ctrl])
n_ssw = len(d_ssw)
perm_diffs = np.zeros(n_perm)
rng = np.random.default_rng(42)
for i in range(n_perm):
    perm = rng.permutation(combined)
    perm_diffs[i] = perm[:n_ssw].mean() - perm[n_ssw:].mean()
perm_p = np.mean(perm_diffs <= obs_diff)
print(f"Permutation P (10,000 iter): {perm_p:.4f}")

# ============================================================
# 2. PER-EVENT ANALYSIS
# ============================================================
print("\n" + "="*60)
print("2. PER-EVENT ANALYSIS")
print("="*60)

event_results = []
for ssw_name, ssw_dt in ssw_dates.items():
    evt = ssw_win[ssw_win['ssw_date'] == ssw_name]
    if len(evt) == 0:
        continue
    
    # DOY-matched control for this event
    evt_doys = set(evt['doy'].unique())
    ctrl_evt = ctrl_matched[ctrl_matched['doy'].isin(evt_doys)]
    
    evt_mean = evt['danger_level'].mean()
    ctrl_mean = ctrl_evt['danger_level'].mean() if len(ctrl_evt) > 0 else np.nan
    decrease = evt_mean < ctrl_mean if not np.isnan(ctrl_mean) else None
    
    # Mann-Whitney for this event
    if len(ctrl_evt) > 0 and len(evt) > 5:
        _, evt_p = stats.mannwhitneyu(evt['danger_level'], ctrl_evt['danger_level'], alternative='less')
    else:
        evt_p = np.nan
    
    event_results.append({
        'ssw_date': ssw_name,
        'n_region_days': len(evt),
        'n_regions': evt['region_id'].nunique(),
        'ssw_mean': evt_mean,
        'ctrl_mean': ctrl_mean,
        'decrease': decrease,
        'p_value': evt_p,
    })
    
    print(f"SSW {ssw_name}: n={len(evt)}, {evt['region_id'].nunique()} regions, "
          f"danger={evt_mean:.3f} vs ctrl={ctrl_mean:.3f}, "
          f"decrease={decrease}, P={evt_p:.4f}" if not np.isnan(evt_p) else 
          f"SSW {ssw_name}: n={len(evt)}, danger={evt_mean:.3f}")

n_decrease = sum(1 for e in event_results if e['decrease'])
n_events = len(event_results)
print(f"\nEvents with decrease: {n_decrease}/{n_events}")
sign_p = stats.binomtest(n_decrease, n_events, 0.5, alternative='greater').pvalue if n_events > 0 else 1.0
print(f"Sign test P: {sign_p:.4f}")

# ============================================================
# 3. PER-REGION ANALYSIS
# ============================================================
print("\n" + "="*60)
print("3. PER-REGION ANALYSIS")
print("="*60)

# Region approximate latitudes (from north to south)
region_lats = {
    3009: ('Nord-Troms', 69.5),
    3010: ('Lyngen', 69.6),
    3011: ('Tromsø', 69.3),
    3013: ('Indre Troms', 68.8),
    3022: ('Trollheimen', 62.8),
    3023: ('Romsdal', 62.5),
    3024: ('Sunnmøre', 62.2),
    3027: ('Indre Fjordane', 61.5),
    3028: ('Jotunheimen', 61.5),
    3029: ('Indre Sogn', 61.2),
}

region_results = []
for region_id in sorted(ssw_win['region_id'].unique()):
    reg_ssw = ssw_win[ssw_win['region_id'] == region_id]
    reg_ctrl = ctrl_matched[ctrl_matched['region_id'] == region_id]
    
    if len(reg_ssw) == 0:
        continue
    
    ssw_mean = reg_ssw['danger_level'].mean()
    ctrl_mean = reg_ctrl['danger_level'].mean() if len(reg_ctrl) > 0 else np.nan
    
    if len(reg_ctrl) > 5:
        _, reg_p = stats.mannwhitneyu(reg_ssw['danger_level'], reg_ctrl['danger_level'], alternative='less')
    else:
        reg_p = np.nan
    
    name, lat = region_lats.get(region_id, ('?', None))
    region_results.append({
        'region_id': region_id,
        'name': name,
        'lat': lat,
        'ssw_mean': ssw_mean,
        'ctrl_mean': ctrl_mean,
        'ratio': ssw_mean / ctrl_mean if ctrl_mean > 0 else np.nan,
        'decrease': ssw_mean < ctrl_mean if not np.isnan(ctrl_mean) else None,
        'p_value': reg_p,
    })
    
    print(f"  {name} ({lat}°N): SSW={ssw_mean:.2f}, Ctrl={ctrl_mean:.2f}, "
          f"ratio={ssw_mean/ctrl_mean:.3f}, P={reg_p:.4f}" if not np.isnan(reg_p) else
          f"  {name}: SSW={ssw_mean:.2f}")

n_reg_decrease = sum(1 for r in region_results if r['decrease'])
print(f"\nRegions with decrease: {n_reg_decrease}/{len(region_results)}")

# ============================================================
# 4. PHASE-RESOLVED ANALYSIS (Pre/Post SSW onset)
# ============================================================
print("\n" + "="*60)
print("4. PHASE-RESOLVED ANALYSIS")
print("="*60)

phases = {
    'Pre (-15 to -1d)': (-15, -1),
    'Post (0 to +15d)': (0, 15),
    'Early (0 to +7d)': (0, 7),
    'Late (+8 to +15d)': (8, 15),
}

for phase_name, (d_start, d_end) in phases.items():
    phase_records = []
    for _, row in ssw_df.iterrows():
        if row['ssw_date']:
            ssw_dt = ssw_dates.get(row['ssw_date'])
            if ssw_dt:
                delta = (row['date'] - ssw_dt).days
                if d_start <= delta <= d_end:
                    phase_records.append(row)
    
    if not phase_records:
        continue
    
    phase_df = pd.DataFrame(phase_records)
    phase_doys = set(phase_df['doy'].unique())
    phase_ctrl = ctrl_df[ctrl_df['doy'].isin(phase_doys)]
    
    if len(phase_ctrl) > 0:
        _, phase_p = stats.mannwhitneyu(phase_df['danger_level'], phase_ctrl['danger_level'], alternative='less')
        print(f"{phase_name}: SSW={phase_df['danger_level'].mean():.3f} (n={len(phase_df)}), "
              f"Ctrl={phase_ctrl['danger_level'].mean():.3f} (n={len(phase_ctrl)}), P={phase_p:.6f}")

# ============================================================
# 5. BOOTSTRAP CI for effect size
# ============================================================
print("\n" + "="*60)
print("5. BOOTSTRAP CI")
print("="*60)

n_boot = 10000
boot_diffs = np.zeros(n_boot)
for i in range(n_boot):
    ssw_sample = rng.choice(d_ssw, size=len(d_ssw), replace=True)
    ctrl_sample = rng.choice(d_ctrl, size=len(d_ctrl), replace=True)
    boot_diffs[i] = ssw_sample.mean() - ctrl_sample.mean()

ci_lo, ci_hi = np.percentile(boot_diffs, [2.5, 97.5])
print(f"Bootstrap 95% CI for danger difference: [{ci_lo:.3f}, {ci_hi:.3f}]")
print(f"Mean difference: {boot_diffs.mean():.3f}")
print(f"CI excludes zero: {ci_hi < 0}")

# ============================================================
# 6. EXPANDED META-ANALYSIS: Norwegian concordance
# ============================================================
print("\n" + "="*60)
print("6. META-ANALYSIS CONCORDANCE")
print("="*60)

# Count event-region pairs
concordance = []
for ssw_name, ssw_dt in ssw_dates.items():
    for region_id in ssw_win['region_id'].unique():
        evt_reg = ssw_win[(ssw_win['ssw_date'] == ssw_name) & (ssw_win['region_id'] == region_id)]
        ctrl_reg = ctrl_matched[ctrl_matched['region_id'] == region_id]
        
        if len(evt_reg) > 0 and len(ctrl_reg) > 0:
            evt_mean = evt_reg['danger_level'].mean()
            ctrl_mean = ctrl_reg['danger_level'].mean()
            concordance.append({
                'ssw': ssw_name,
                'region': region_id,
                'decrease': evt_mean < ctrl_mean,
            })

n_conc = sum(1 for c in concordance if c['decrease'])
n_total = len(concordance)
print(f"Event-region pairs: {n_total}")
print(f"Concordant (decrease): {n_conc}/{n_total} = {n_conc/n_total*100:.1f}%")
conc_p = stats.binomtest(n_conc, n_total, 0.5, alternative='greater').pvalue
print(f"Binomial P: {conc_p:.6f}")

# ============================================================
# SAVE RESULTS
# ============================================================
results = {
    'overview': {
        'n_regions': int(ssw_win['region_id'].nunique()),
        'n_ssw_events': len(ssw_dates),
        'n_region_days_ssw': len(ssw_win),
        'n_region_days_ctrl': len(ctrl_matched),
        'ssw_mean_danger': float(ssw_win['danger_level'].mean()),
        'ctrl_mean_danger': float(ctrl_matched['danger_level'].mean()),
        'mann_whitney_p': float(mw_p),
        'cohens_d': float(cohens_d),
        'permutation_p': float(perm_p),
        'bootstrap_ci': [float(ci_lo), float(ci_hi)],
    },
    'per_event': event_results,
    'per_region': region_results,
    'concordance': {
        'n_pairs': n_total,
        'n_decrease': n_conc,
        'pct_decrease': float(n_conc/n_total*100),
        'binomial_p': float(conc_p),
    },
}

os.makedirs('data/results', exist_ok=True)
with open('data/results/r38_norway_expanded.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nSaved to data/results/r38_norway_expanded.json")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*60)
print("SUMMARY: EXPANDED NORWEGIAN REPLICATION")
print("="*60)
print(f"• 10 mountain regions across 61-70°N latitude")
print(f"• 4 SSW events (2018, 2019, 2021, 2023)")
print(f"• SSW ±15d danger: {ssw_win['danger_level'].mean():.3f} vs control {ctrl_matched['danger_level'].mean():.3f}")
print(f"• Mann-Whitney P = {mw_p:.6f}")
print(f"• Permutation P = {perm_p:.4f}")
print(f"• Cohen's d = {cohens_d:.3f}")
print(f"• Bootstrap 95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]")
print(f"• Events with decrease: {n_decrease}/{n_events}")
print(f"• Region-event concordance: {n_conc}/{n_total} ({n_conc/n_total*100:.1f}%)")
print(f"• Regions with decrease: {n_reg_decrease}/{len(region_results)}")
