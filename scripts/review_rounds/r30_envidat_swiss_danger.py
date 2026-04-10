"""
R30: EnviDat Swiss Regional Danger Level Analysis
Analyzes SSW impact on re-analyzed danger levels across 146 Swiss warning sectors.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json, os, warnings
warnings.filterwarnings('ignore')

OUT = 'data/results'
os.makedirs(OUT, exist_ok=True)

# ── 1. Load and harmonize all 3 danger-level files ──────────────────────────
print("=" * 70)
print("R30: EnviDat Swiss Regional Danger Level Analysis")
print("=" * 70)

df1 = pd.read_csv('data/cryosphere/envidat/swiss_danger_2001_2020.csv')
df1 = df1.rename(columns={'validDate': 'date', 'sectorId': 'sector'})
df1['date'] = pd.to_datetime(df1['date'])
# 2001-2020 has dryWet column; keep only dry
if 'dryWet' in df1.columns:
    df1 = df1[df1['dryWet'] == 'dry'].copy()
df1 = df1[['date', 'sector', 'dangerLevel']].copy()

df2 = pd.read_csv('data/cryosphere/envidat/swiss_danger_2020_2023.csv')
df2 = df2.rename(columns={'validDate': 'date', 'sectorId': 'sector'})
df2['date'] = pd.to_datetime(df2['date'])
df2 = df2[['date', 'sector', 'dangerLevel']].copy()

df3 = pd.read_csv('data/cryosphere/envidat/swiss_danger_2023_2024.csv')
df3 = df3.rename(columns={'date': 'date', 'sector_id': 'sector', 'dangerlevel_tidy': 'dangerLevel'})
df3['date'] = pd.to_datetime(df3['date'])
df3 = df3[['date', 'sector', 'dangerLevel']].copy()

df = pd.concat([df1, df2, df3], ignore_index=True)
df = df.dropna(subset=['dangerLevel'])
df['dangerLevel'] = df['dangerLevel'].astype(int)

print(f"\nCombined dataset: {len(df):,} region-days")
print(f"Date range: {df['date'].min().date()} to {df['date'].max().date()}")
print(f"Unique sectors: {df['sector'].nunique()}")
print(f"Danger level distribution:")
for lv in sorted(df['dangerLevel'].unique()):
    n = (df['dangerLevel'] == lv).sum()
    print(f"  Level {lv}: {n:,} ({100*n/len(df):.1f}%)")

# ── 2. Load SSW catalog ─────────────────────────────────────────────────────
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_cat = ssw_cat.reset_index()
if 'onset_date' in ssw_cat.columns:
    ssw_dates = pd.to_datetime(ssw_cat['onset_date']).dt.tz_localize(None)
else:
    ssw_dates = pd.to_datetime(ssw_cat.iloc[:, 0]).dt.tz_localize(None)

# Filter to data range
data_start, data_end = df['date'].min(), df['date'].max()
ssw_dates = ssw_dates[(ssw_dates >= data_start - pd.Timedelta(days=30)) & 
                       (ssw_dates <= data_end + pd.Timedelta(days=30))]
# Keep only winter SSWs (Nov-Mar)
ssw_dates = ssw_dates[ssw_dates.dt.month.isin([11, 12, 1, 2, 3])]

print(f"\nSSW events in data range: {len(ssw_dates)}")
for d in sorted(ssw_dates):
    print(f"  {d.date()}")

# ── 3. Define SSW windows and DOY-matched controls ──────────────────────────
WINDOW = 30  # days after SSW onset

# For each sector-day, mark if within SSW window
df['doy'] = df['date'].dt.dayofyear
df['season'] = df['date'].apply(lambda x: x.year if x.month >= 9 else x.year - 1)
df['is_winter'] = df['date'].dt.month.isin([11, 12, 1, 2, 3, 4])
df_winter = df[df['is_winter']].copy()

# Mark SSW windows
df_winter['ssw_window'] = False
df_winter['ssw_event'] = None
for ssw in ssw_dates:
    mask = (df_winter['date'] >= ssw) & (df_winter['date'] < ssw + pd.Timedelta(days=WINDOW))
    df_winter.loc[mask, 'ssw_window'] = True
    df_winter.loc[mask, 'ssw_event'] = ssw.strftime('%Y-%m-%d')

# DOY-matched control: same DOY range from non-SSW winters
df_winter['doy_adj'] = df_winter['date'].apply(
    lambda x: x.timetuple().tm_yday if x.month <= 6 else x.timetuple().tm_yday - 365
)

print(f"\nWinter region-days: {len(df_winter):,}")
print(f"SSW window region-days: {df_winter['ssw_window'].sum():,}")
print(f"Control region-days: {(~df_winter['ssw_window']).sum():,}")

# ── 4. Aggregate analysis: SSW vs DOY-matched control ───────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 1: Overall SSW vs Control (DOY-matched)")
print("=" * 70)

results = {}

# For each SSW event, get the DOY range, then find control days with same DOY range
# from OTHER winters
event_results = []
for ssw in sorted(ssw_dates):
    ssw_str = ssw.strftime('%Y-%m-%d')
    ssw_season = ssw.year if ssw.month >= 9 else ssw.year - 1
    
    # SSW window data
    ssw_mask = df_winter['ssw_event'] == ssw_str
    ssw_data = df_winter[ssw_mask]
    
    if len(ssw_data) == 0:
        continue
    
    # DOY range for this SSW window
    doy_min = ssw_data['doy_adj'].min()
    doy_max = ssw_data['doy_adj'].max()
    
    # Control: same DOY range, different seasons, NOT in any SSW window
    ctrl_mask = (
        (df_winter['doy_adj'] >= doy_min) & 
        (df_winter['doy_adj'] <= doy_max) & 
        (df_winter['season'] != ssw_season) & 
        (~df_winter['ssw_window'])
    )
    ctrl_data = df_winter[ctrl_mask]
    
    if len(ctrl_data) == 0:
        continue
    
    ssw_mean = ssw_data['dangerLevel'].mean()
    ctrl_mean = ctrl_data['dangerLevel'].mean()
    diff = ssw_mean - ctrl_mean
    
    # Mann-Whitney test
    u_stat, u_p = stats.mannwhitneyu(
        ssw_data['dangerLevel'], ctrl_data['dangerLevel'], alternative='two-sided'
    )
    
    # Effect size (rank-biserial)
    n1, n2 = len(ssw_data), len(ctrl_data)
    rbc = 1 - (2 * u_stat) / (n1 * n2)
    
    event_results.append({
        'ssw_date': ssw_str,
        'ssw_mean': round(ssw_mean, 3),
        'ctrl_mean': round(ctrl_mean, 3),
        'diff': round(diff, 3),
        'n_ssw': n1,
        'n_ctrl': n2,
        'mw_p': u_p,
        'rank_biserial': round(rbc, 3),
        'direction': 'decrease' if diff < 0 else 'increase'
    })
    
    sig = '*' if u_p < 0.05 else ''
    print(f"  {ssw_str}: SSW={ssw_mean:.2f} Ctrl={ctrl_mean:.2f} "
          f"Δ={diff:+.3f} P={u_p:.4f}{sig} n={n1}/{n2} rbc={rbc:.3f}")

n_decrease = sum(1 for e in event_results if e['diff'] < 0)
n_total = len(event_results)
print(f"\n  Events with decrease: {n_decrease}/{n_total}")

# Sign test
sign_p = stats.binomtest(n_decrease, n_total, 0.5).pvalue
print(f"  Sign test P = {sign_p:.6f}")

# Aggregate means
ssw_means = [e['ssw_mean'] for e in event_results]
ctrl_means = [e['ctrl_mean'] for e in event_results]
t_stat, t_p = stats.ttest_rel(ssw_means, ctrl_means)
wilcox_stat, wilcox_p = stats.wilcoxon([e['diff'] for e in event_results])
print(f"  Paired t-test: t={t_stat:.3f}, P={t_p:.6f}")
print(f"  Wilcoxon signed-rank: P={wilcox_p:.6f}")

# Cohen's d
diffs = np.array([e['diff'] for e in event_results])
cohen_d = diffs.mean() / diffs.std()
print(f"  Mean Δ = {diffs.mean():.3f}, Cohen's d = {cohen_d:.3f}")

results['overall'] = {
    'n_events': n_total,
    'n_decrease': n_decrease,
    'sign_p': sign_p,
    't_p': t_p,
    'wilcoxon_p': wilcox_p,
    'mean_diff': round(float(diffs.mean()), 4),
    'cohen_d': round(float(cohen_d), 3),
    'events': event_results
}

# ── 5. Sector-level analysis ────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 2: Sector-Level Concordance")
print("=" * 70)

sector_results = []
for sector in sorted(df_winter['sector'].unique()):
    sec_data = df_winter[df_winter['sector'] == sector]
    
    sec_ssw = sec_data[sec_data['ssw_window']]
    sec_ctrl = sec_data[~sec_data['ssw_window']]
    
    if len(sec_ssw) < 10 or len(sec_ctrl) < 50:
        continue
    
    ssw_m = sec_ssw['dangerLevel'].mean()
    ctrl_m = sec_ctrl['dangerLevel'].mean()
    diff = ssw_m - ctrl_m
    
    try:
        u, p = stats.mannwhitneyu(sec_ssw['dangerLevel'], sec_ctrl['dangerLevel'], alternative='two-sided')
    except:
        p = 1.0
    
    sector_results.append({
        'sector': sector,
        'ssw_mean': round(ssw_m, 3),
        'ctrl_mean': round(ctrl_m, 3),
        'diff': round(diff, 3),
        'p': p,
        'n_ssw': len(sec_ssw),
        'n_ctrl': len(sec_ctrl)
    })

n_sec_decrease = sum(1 for s in sector_results if s['diff'] < 0)
n_sec_sig_decrease = sum(1 for s in sector_results if s['diff'] < 0 and s['p'] < 0.05)
n_sec_total = len(sector_results)
sec_sign_p = stats.binomtest(n_sec_decrease, n_sec_total, 0.5).pvalue

print(f"  Sectors analyzed: {n_sec_total}")
print(f"  Sectors with decrease: {n_sec_decrease}/{n_sec_total} ({100*n_sec_decrease/n_sec_total:.1f}%)")
print(f"  Sectors with sig decrease (P<0.05): {n_sec_sig_decrease}")
print(f"  Sector sign test P = {sec_sign_p:.6f}")

# Mean effect across sectors
sec_diffs = [s['diff'] for s in sector_results]
print(f"  Mean sector Δ = {np.mean(sec_diffs):.4f}")
print(f"  Median sector Δ = {np.median(sec_diffs):.4f}")

results['sector_level'] = {
    'n_sectors': n_sec_total,
    'n_decrease': n_sec_decrease,
    'n_sig_decrease': n_sec_sig_decrease,
    'sign_p': sec_sign_p,
    'mean_diff': round(float(np.mean(sec_diffs)), 4),
    'median_diff': round(float(np.median(sec_diffs)), 4)
}

# ── 6. Phase-resolved analysis (Pre/During/Post SSW) ────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 3: Phase-Resolved (Pre / During / Post SSW)")
print("=" * 70)

phases = {
    'pre_15d': (-15, 0),
    'during_15d': (0, 15),
    'post_15d': (15, 30),
    'post_30d': (30, 45)
}

for phase_name, (d_start, d_end) in phases.items():
    phase_diffs = []
    for ssw in sorted(ssw_dates):
        ssw_season = ssw.year if ssw.month >= 9 else ssw.year - 1
        
        # Phase window
        phase_mask = (
            (df_winter['date'] >= ssw + pd.Timedelta(days=d_start)) &
            (df_winter['date'] < ssw + pd.Timedelta(days=d_end))
        )
        phase_data = df_winter[phase_mask]
        
        if len(phase_data) == 0:
            continue
        
        # DOY-matched control
        doy_min = phase_data['doy_adj'].min()
        doy_max = phase_data['doy_adj'].max()
        ctrl_mask = (
            (df_winter['doy_adj'] >= doy_min) & 
            (df_winter['doy_adj'] <= doy_max) & 
            (df_winter['season'] != ssw_season) & 
            (~df_winter['ssw_window'])
        )
        ctrl_data = df_winter[ctrl_mask]
        
        if len(ctrl_data) == 0:
            continue
        
        phase_diffs.append(phase_data['dangerLevel'].mean() - ctrl_data['dangerLevel'].mean())
    
    if len(phase_diffs) > 2:
        arr = np.array(phase_diffs)
        n_neg = sum(1 for d in arr if d < 0)
        t, tp = stats.ttest_1samp(arr, 0)
        print(f"  {phase_name}: mean Δ={arr.mean():+.4f}, "
              f"{n_neg}/{len(arr)} decrease, t={t:.3f}, P={tp:.4f}")
        results[phase_name] = {
            'mean_diff': round(float(arr.mean()), 4),
            'n_decrease': n_neg,
            'n_total': len(arr),
            't_p': tp
        }

# ── 7. Danger level ≥3 (Considerable+) frequency analysis ───────────────────
print("\n" + "=" * 70)
print("ANALYSIS 4: High Danger (≥3) Frequency During SSW")
print("=" * 70)

for threshold in [3, 4]:
    label = {3: 'Considerable+', 4: 'High+'}[threshold]
    event_freqs = []
    
    for ssw in sorted(ssw_dates):
        ssw_str = ssw.strftime('%Y-%m-%d')
        ssw_season = ssw.year if ssw.month >= 9 else ssw.year - 1
        
        ssw_mask = df_winter['ssw_event'] == ssw_str
        ssw_data = df_winter[ssw_mask]
        
        if len(ssw_data) == 0:
            continue
        
        doy_min = ssw_data['doy_adj'].min()
        doy_max = ssw_data['doy_adj'].max()
        ctrl_mask = (
            (df_winter['doy_adj'] >= doy_min) & 
            (df_winter['doy_adj'] <= doy_max) & 
            (df_winter['season'] != ssw_season) & 
            (~df_winter['ssw_window'])
        )
        ctrl_data = df_winter[ctrl_mask]
        
        if len(ctrl_data) == 0:
            continue
        
        ssw_freq = (ssw_data['dangerLevel'] >= threshold).mean()
        ctrl_freq = (ctrl_data['dangerLevel'] >= threshold).mean()
        
        event_freqs.append({
            'ssw': ssw_str,
            'ssw_freq': ssw_freq,
            'ctrl_freq': ctrl_freq,
            'rr': ssw_freq / ctrl_freq if ctrl_freq > 0 else np.nan
        })
    
    rrs = [e['rr'] for e in event_freqs if not np.isnan(e['rr'])]
    n_rr_below1 = sum(1 for r in rrs if r < 1)
    
    if len(rrs) > 2:
        geo_mean_rr = np.exp(np.mean(np.log(np.array(rrs))))
        t, tp = stats.ttest_1samp(np.log(np.array(rrs)), 0)
        print(f"  {label}: geo-mean RR={geo_mean_rr:.3f}, "
              f"{n_rr_below1}/{len(rrs)} events RR<1, log-RR t-test P={tp:.4f}")
        results[f'high_danger_{threshold}'] = {
            'geo_mean_rr': round(float(geo_mean_rr), 3),
            'n_below_1': n_rr_below1,
            'n_total': len(rrs),
            'log_rr_p': tp
        }

# ── 8. Cross-validation with existing SLF count data ────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 5: Count vs Rating Discrepancy (Existing Panel)")
print("=" * 70)

try:
    panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
    panel = panel.reset_index()
    panel['date'] = pd.to_datetime(panel['time'])
    
    # For each SSW, compare count-based and rating-based measures
    for ssw in sorted(ssw_dates)[:5]:  # First 5 for display
        ssw_str = ssw.strftime('%Y-%m-%d')
        ssw_season = ssw.year if ssw.month >= 9 else ssw.year - 1
        
        # Count data from panel
        count_mask = (panel['date'] >= ssw) & (panel['date'] < ssw + pd.Timedelta(days=30))
        if 'dry_natural_size_1234' in panel.columns:
            ssw_counts = panel.loc[count_mask, 'dry_natural_size_1234']
            
            doy_min = panel.loc[count_mask, 'date'].dt.dayofyear.min() if count_mask.sum() > 0 else 0
            doy_max = panel.loc[count_mask, 'date'].dt.dayofyear.max() if count_mask.sum() > 0 else 0
            
            ctrl_mask = (
                (panel['date'].dt.dayofyear >= doy_min) & 
                (panel['date'].dt.dayofyear <= doy_max) & 
                (panel['date'].apply(lambda x: x.year if x.month >= 9 else x.year - 1) != ssw_season)
            )
            ctrl_counts = panel.loc[ctrl_mask, 'dry_natural_size_1234']
            
            # Rating data from EnviDat
            rat_mask = df_winter['ssw_event'] == ssw_str
            ssw_rating = df_winter.loc[rat_mask, 'dangerLevel'].mean() if rat_mask.sum() > 0 else np.nan
            
            if count_mask.sum() > 0 and ctrl_mask.sum() > 0 and not np.isnan(ssw_rating):
                count_rr = ssw_counts.mean() / ctrl_counts.mean() if ctrl_counts.mean() > 0 else np.nan
                print(f"  {ssw_str}: Count RR={count_rr:.3f}, Rating SSW mean={ssw_rating:.2f}")
except Exception as e:
    print(f"  Panel comparison skipped: {e}")

# ── 9. Human-triggered avalanche analysis ────────────────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 6: Human-Triggered Avalanche Rates")
print("=" * 70)

try:
    ht = pd.read_csv('data/cryosphere/envidat/human-triggered_avalanches.csv')
    ht['date'] = pd.to_datetime(ht['date'])
    ht['season'] = ht['date'].apply(lambda x: x.year if x.month >= 9 else x.year - 1)
    ht['doy_adj'] = ht['date'].apply(
        lambda x: x.timetuple().tm_yday if x.month <= 6 else x.timetuple().tm_yday - 365
    )
    
    ht_results = []
    for ssw in sorted(ssw_dates):
        ssw_str = ssw.strftime('%Y-%m-%d')
        ssw_season = ssw.year if ssw.month >= 9 else ssw.year - 1
        
        ssw_mask = (ht['date'] >= ssw) & (ht['date'] < ssw + pd.Timedelta(days=WINDOW))
        ssw_data = ht[ssw_mask]
        
        if len(ssw_data) == 0:
            continue
        
        doy_min = ssw_data['doy_adj'].min()
        doy_max = ssw_data['doy_adj'].max()
        ctrl_mask = (
            (ht['doy_adj'] >= doy_min) & 
            (ht['doy_adj'] <= doy_max) & 
            (ht['season'] != ssw_season) & 
            (~((ht['date'] >= ssw - pd.Timedelta(days=15)) & (ht['date'] < ssw + pd.Timedelta(days=45))))
        )
        ctrl_data = ht[ctrl_mask]
        
        if len(ctrl_data) == 0 or 'rhoSz2' not in ht.columns:
            continue
        
        ssw_rate = ssw_data['rhoSz2'].mean()
        ctrl_rate = ctrl_data['rhoSz2'].mean()
        rr = ssw_rate / ctrl_rate if ctrl_rate > 0 else np.nan
        
        ht_results.append({
            'ssw': ssw_str,
            'ssw_rate': ssw_rate,
            'ctrl_rate': ctrl_rate,
            'rr': rr
        })
        
        print(f"  {ssw_str}: SSW rate={ssw_rate:.5f}, Ctrl rate={ctrl_rate:.5f}, RR={rr:.3f}")
    
    if ht_results:
        rrs = [e['rr'] for e in ht_results if not np.isnan(e['rr'])]
        n_below = sum(1 for r in rrs if r < 1)
        geo_rr = np.exp(np.mean(np.log(np.array(rrs)))) if rrs else np.nan
        print(f"\n  Human-triggered: {n_below}/{len(rrs)} events RR<1, geo-mean RR={geo_rr:.3f}")
        results['human_triggered'] = {
            'n_events': len(rrs),
            'n_decrease': n_below,
            'geo_mean_rr': round(float(geo_rr), 3)
        }
except Exception as e:
    print(f"  Human-triggered analysis skipped: {e}")

# ── 10. Rutschblock stability test analysis ──────────────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS 7: Rutschblock Test Results")
print("=" * 70)

try:
    rb = pd.read_csv('data/cryosphere/envidat/rutschblock_test_results.csv')
    rb['date'] = pd.to_datetime(rb['date'])
    rb['season'] = rb['date'].apply(lambda x: x.year if x.month >= 9 else x.year - 1)
    rb['doy_adj'] = rb['date'].apply(
        lambda x: x.timetuple().tm_yday if x.month <= 6 else x.timetuple().tm_yday - 365
    )
    
    print(f"  Total Rutschblock tests: {len(rb):,}")
    print(f"  Date range: {rb['date'].min().date()} to {rb['date'].max().date()}")
    print(f"  rbStab values: {sorted(rb['rbStab'].unique())}")
    
    rb_results = []
    for ssw in sorted(ssw_dates):
        ssw_str = ssw.strftime('%Y-%m-%d')
        ssw_season = ssw.year if ssw.month >= 9 else ssw.year - 1
        
        ssw_mask = (rb['date'] >= ssw) & (rb['date'] < ssw + pd.Timedelta(days=WINDOW))
        ssw_data = rb[ssw_mask]
        
        if len(ssw_data) < 3:
            continue
        
        doy_min = ssw_data['doy_adj'].min()
        doy_max = ssw_data['doy_adj'].max()
        ctrl_mask = (
            (rb['doy_adj'] >= doy_min) & 
            (rb['doy_adj'] <= doy_max) & 
            (rb['season'] != ssw_season)
        )
        ctrl_data = rb[ctrl_mask]
        
        if len(ctrl_data) < 10:
            continue
        
        ssw_stab = ssw_data['rbStab'].mean()
        ctrl_stab = ctrl_data['rbStab'].mean()
        diff = ssw_stab - ctrl_stab
        
        u, p = stats.mannwhitneyu(ssw_data['rbStab'], ctrl_data['rbStab'], alternative='two-sided')
        
        rb_results.append({
            'ssw': ssw_str,
            'ssw_stab': round(ssw_stab, 2),
            'ctrl_stab': round(ctrl_stab, 2),
            'diff': round(diff, 2),
            'p': p,
            'n_ssw': len(ssw_data),
            'n_ctrl': len(ctrl_data)
        })
        
        sig = '*' if p < 0.05 else ''
        print(f"  {ssw_str}: SSW={ssw_stab:.2f} Ctrl={ctrl_stab:.2f} "
              f"Δ={diff:+.2f} P={p:.4f}{sig} n={len(ssw_data)}/{len(ctrl_data)}")
    
    if rb_results:
        rb_diffs = [r['diff'] for r in rb_results]
        n_pos = sum(1 for d in rb_diffs if d > 0)  # Higher = more stable
        print(f"\n  Events with HIGHER stability during SSW: {n_pos}/{len(rb_results)}")
        print(f"  Mean Δ stability = {np.mean(rb_diffs):+.3f}")
        results['rutschblock'] = {
            'n_events': len(rb_results),
            'n_more_stable': n_pos,
            'mean_diff': round(float(np.mean(rb_diffs)), 3)
        }
except Exception as e:
    print(f"  Rutschblock analysis skipped: {e}")

# ── 11. Combined multi-source concordance ────────────────────────────────────
print("\n" + "=" * 70)
print("GRAND SYNTHESIS: Multi-Source Swiss Concordance")
print("=" * 70)

print(f"""
Source                  | Direction | N events | Key statistic
------------------------|-----------|----------|------------------
EnviDat danger ratings  | {'↓' if results['overall']['mean_diff'] < 0 else '↑'} {abs(results['overall']['mean_diff']):.3f}   | {results['overall']['n_events']}       | Sign {results['overall']['n_decrease']}/{results['overall']['n_events']}, P={results['overall']['sign_p']:.4f}
Sector concordance      | {results['sector_level']['n_decrease']}/{results['sector_level']['n_sectors']} ↓   | {results['sector_level']['n_sectors']} sectors | Sign P={results['sector_level']['sign_p']:.4f}
SLF occurrence counts   | ↓ 68%     | 15       | t P=0.003, sign P=0.001
Norway NVE ratings      | ↓ 0.27    | 4        | MW P<0.0001
""")

# Save results
with open(f'{OUT}/r30_envidat_analysis.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to {OUT}/r30_envidat_analysis.json")
print("=" * 70)
