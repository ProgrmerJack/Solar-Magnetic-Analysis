"""
R16: SNOWPACK Stability Index Analysis — Mechanism Validation
Analyzes SNOWPACK-simulated stability indices around SSW events
to directly validate the sintering/stabilization mechanism.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json, warnings
warnings.filterwarnings('ignore')

results = {}

# ── Load SNOWPACK data ──
sp = pd.read_csv('data/cryosphere/swiss_snowpack/data_rf2_tidy.csv')
sp['date'] = pd.to_datetime(sp['datum'])
print(f"SNOWPACK: {len(sp)} rows, {sp['date'].min()} to {sp['date'].max()}")
print(f"Stations: {sp['station_code'].nunique() if 'station_code' in sp.columns else 'N/A'}")
print(f"Columns (stability-related): {[c for c in sp.columns if any(x in c.lower() for x in ['ssi','sk3','sn3','ccl','pwl','hs_','swe','ta_','ts0','ts1','ts2','iswr','hn24','danger'])]}")

# ── Load SSW catalog ──
ssw_exist = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_exist = ssw_exist.reset_index()
ssw_exist['onset_date'] = pd.to_datetime(ssw_exist['onset_date']).dt.tz_localize(None)
print(f"\nSSW catalog: {len(ssw_exist)} events")

# Also load Butler catalog if available
try:
    butler = pd.read_csv('data/processed/atmospheric/butler_ssw_compendium_era5.csv')
    butler['onset_date'] = pd.to_datetime(butler['date'])
    print(f"Butler SSW: {len(butler)} events")
except:
    butler = pd.DataFrame()

# Merge catalogs
all_dates = set(ssw_exist['onset_date'].dt.date)
if len(butler) > 0:
    all_dates |= set(butler['onset_date'].dt.date)
ssw_dates = sorted(all_dates)
print(f"Unique SSW dates: {len(ssw_dates)}")

# Filter to SNOWPACK date range
sp_min, sp_max = sp['date'].min(), sp['date'].max()
ssw_in_range = [d for d in ssw_dates if sp_min.date() <= d <= sp_max.date()]
print(f"SSW events in SNOWPACK range ({sp_min.date()} to {sp_max.date()}): {len(ssw_in_range)}")
for d in ssw_in_range:
    print(f"  {d}")

# ── Stability index columns ──
stability_cols = [c for c in sp.columns if any(x in c.lower() for x in ['ssi_pwl','sk38_pwl','sn38_pwl','ccl_pwl'])]
weather_cols = [c for c in sp.columns if c in ['HS_mod','HS_meas','SWE','TA','TS0','TS1','TS2','ISWR','HN24','HN24_7d','dangerLevel','min_ccl_pen']]
pwl_cols = [c for c in sp.columns if c.startswith('pwl_')]

print(f"\nStability columns: {stability_cols}")
print(f"Weather columns: {weather_cols[:10]}")
print(f"PWL columns: {pwl_cols[:10]}")

# ── Phase 1: SSW Composite Analysis ──
print("\n" + "="*70)
print("PHASE 1: SNOWPACK STABILITY AROUND SSW EVENTS")
print("="*70)

window = 30  # days post-SSW
pre_window = 30  # days pre-SSW for control

all_cols = stability_cols + weather_cols[:8] + pwl_cols[:5]
# Remove duplicates and missing
all_cols = [c for c in all_cols if c in sp.columns]

phase1_results = {}
for col in all_cols:
    pre_vals = []
    post_vals = []
    event_results = []
    
    for ssw_date in ssw_in_range:
        ssw_dt = pd.Timestamp(ssw_date)
        
        # Pre-SSW window (control)
        pre_mask = (sp['date'] >= ssw_dt - pd.Timedelta(days=pre_window)) & \
                   (sp['date'] < ssw_dt)
        pre_data = sp.loc[pre_mask, col].dropna()
        
        # Post-SSW window
        post_mask = (sp['date'] >= ssw_dt) & \
                    (sp['date'] < ssw_dt + pd.Timedelta(days=window))
        post_data = sp.loc[post_mask, col].dropna()
        
        if len(pre_data) > 10 and len(post_data) > 10:
            pre_mean = pre_data.mean()
            post_mean = post_data.mean()
            pre_vals.append(pre_mean)
            post_vals.append(post_mean)
            pct = ((post_mean - pre_mean) / abs(pre_mean) * 100) if pre_mean != 0 else 0
            event_results.append({
                'date': str(ssw_date),
                'pre_mean': round(pre_mean, 4),
                'post_mean': round(post_mean, 4),
                'pct_change': round(pct, 1),
                'direction': 'decrease' if post_mean < pre_mean else 'increase'
            })
    
    if len(pre_vals) >= 3:
        pre_arr = np.array(pre_vals)
        post_arr = np.array(post_vals)
        diffs = post_arr - pre_arr
        n_decrease = sum(1 for d in diffs if d < 0)
        n_total = len(diffs)
        
        # Statistical tests
        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(pre_arr, post_arr)
        ttest_stat, ttest_p = stats.ttest_rel(pre_arr, post_arr)
        sign_p = stats.binomtest(n_decrease, n_total, 0.5).pvalue
        
        mean_pct = np.mean([(p - q) / abs(q) * 100 if q != 0 else 0 for p, q in zip(post_arr, pre_arr)])
        
        r = {
            'n_events': n_total,
            'n_decrease': n_decrease,
            'mean_pct_change': round(mean_pct, 1),
            'wilcoxon_p': round(wilcoxon_p, 6),
            'ttest_p': round(ttest_p, 6),
            'sign_p': round(sign_p, 6),
            'cohens_d': round(np.mean(diffs) / np.std(diffs, ddof=1), 3) if np.std(diffs) > 0 else 0,
            'events': event_results
        }
        phase1_results[col] = r
        
        sig = "***" if wilcoxon_p < 0.001 else "**" if wilcoxon_p < 0.01 else "*" if wilcoxon_p < 0.05 else "†" if wilcoxon_p < 0.1 else ""
        print(f"  {col:20s}: {mean_pct:+6.1f}%, {n_decrease}/{n_total} decrease, "
              f"P_wilcoxon={wilcoxon_p:.4f} {sig}, d={r['cohens_d']:+.3f}")

results['phase1_stability'] = phase1_results

# ── Phase 2: Calendar-Matched Control ──
print("\n" + "="*70)
print("PHASE 2: CALENDAR-MATCHED SNOWPACK ANALYSIS")
print("="*70)

# For each SSW event, find the same calendar dates in non-SSW years
ssw_years = set(d.year for d in ssw_in_range)
# Also include years with nearby SSWs (within 60 days)
all_years = set(sp['date'].dt.year.unique())
# Non-SSW winters: years where no SSW falls in the SNOWPACK range
winter_years = sorted(all_years - ssw_years)
print(f"SSW years: {sorted(ssw_years)}")
print(f"Non-SSW years (potential controls): {winter_years}")

phase2_results = {}
key_cols = [c for c in stability_cols + weather_cols[:6] if c in sp.columns]

for col in key_cols:
    ssw_means = []
    ctrl_means = []
    
    for ssw_date in ssw_in_range:
        ssw_dt = pd.Timestamp(ssw_date)
        month, day = ssw_dt.month, ssw_dt.day
        
        # SSW post-window mean
        post_mask = (sp['date'] >= ssw_dt) & \
                    (sp['date'] < ssw_dt + pd.Timedelta(days=window))
        ssw_data = sp.loc[post_mask, col].dropna()
        
        if len(ssw_data) < 10:
            continue
        
        ssw_mean = ssw_data.mean()
        
        # Same calendar dates in control years
        ctrl_year_means = []
        for yr in winter_years:
            try:
                ctrl_start = pd.Timestamp(year=yr, month=month, day=day)
                ctrl_mask = (sp['date'] >= ctrl_start) & \
                            (sp['date'] < ctrl_start + pd.Timedelta(days=window))
                ctrl_data = sp.loc[ctrl_mask, col].dropna()
                if len(ctrl_data) > 10:
                    ctrl_year_means.append(ctrl_data.mean())
            except:
                pass
        
        if len(ctrl_year_means) >= 1:
            ctrl_mean = np.mean(ctrl_year_means)
            ssw_means.append(ssw_mean)
            ctrl_means.append(ctrl_mean)
    
    if len(ssw_means) >= 3:
        ssw_arr = np.array(ssw_means)
        ctrl_arr = np.array(ctrl_means)
        diffs = ssw_arr - ctrl_arr
        n_decrease = sum(1 for d in diffs if d < 0)
        
        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(ssw_arr, ctrl_arr)
        mean_pct = np.mean([(s - c) / abs(c) * 100 if c != 0 else 0 for s, c in zip(ssw_arr, ctrl_arr)])
        
        phase2_results[col] = {
            'n_events': len(ssw_means),
            'n_decrease': n_decrease,
            'mean_pct_change': round(mean_pct, 1),
            'wilcoxon_p': round(wilcoxon_p, 6),
            'n_control_years': len(winter_years)
        }
        
        sig = "***" if wilcoxon_p < 0.001 else "**" if wilcoxon_p < 0.01 else "*" if wilcoxon_p < 0.05 else "†" if wilcoxon_p < 0.1 else ""
        print(f"  {col:20s}: {mean_pct:+6.1f}% (cal-matched), {n_decrease}/{len(ssw_means)} decrease, "
              f"P={wilcoxon_p:.4f} {sig}")

results['phase2_calendar_matched'] = phase2_results

# ── Phase 3: Phase-Resolved (Early/Mid/Late) ──
print("\n" + "="*70)
print("PHASE 3: PHASE-RESOLVED SNOWPACK STABILITY")
print("="*70)

phases = [
    ('early', 0, 14),
    ('mid', 15, 29),
    ('late', 30, 44)
]

phase3_results = {}
for col in key_cols:
    phase3_results[col] = {}
    for phase_name, start, end in phases:
        ssw_vals = []
        ctrl_vals = []
        
        for ssw_date in ssw_in_range:
            ssw_dt = pd.Timestamp(ssw_date)
            month, day = ssw_dt.month, ssw_dt.day
            
            # SSW phase window
            phase_mask = (sp['date'] >= ssw_dt + pd.Timedelta(days=start)) & \
                         (sp['date'] < ssw_dt + pd.Timedelta(days=end+1))
            ssw_data = sp.loc[phase_mask, col].dropna()
            
            if len(ssw_data) < 5:
                continue
            
            # Calendar-matched controls
            ctrl_year_vals = []
            for yr in winter_years:
                try:
                    ctrl_start = pd.Timestamp(year=yr, month=month, day=day) + pd.Timedelta(days=start)
                    ctrl_mask = (sp['date'] >= ctrl_start) & \
                                (sp['date'] < ctrl_start + pd.Timedelta(days=end-start+1))
                    ctrl_data = sp.loc[ctrl_mask, col].dropna()
                    if len(ctrl_data) > 5:
                        ctrl_year_vals.append(ctrl_data.mean())
                except:
                    pass
            
            if ctrl_year_vals:
                ssw_vals.append(ssw_data.mean())
                ctrl_vals.append(np.mean(ctrl_year_vals))
        
        if len(ssw_vals) >= 3:
            ssw_arr = np.array(ssw_vals)
            ctrl_arr = np.array(ctrl_vals)
            diffs = ssw_arr - ctrl_arr
            n_dec = sum(1 for d in diffs if d < 0)
            
            try:
                w_stat, w_p = stats.wilcoxon(ssw_arr, ctrl_arr)
            except:
                w_p = 1.0
            
            mean_pct = np.mean([(s-c)/abs(c)*100 if c!=0 else 0 for s,c in zip(ssw_arr, ctrl_arr)])
            
            phase3_results[col][phase_name] = {
                'n': len(ssw_vals),
                'mean_pct': round(mean_pct, 1),
                'n_decrease': n_dec,
                'wilcoxon_p': round(w_p, 4)
            }
            
            sig = "*" if w_p < 0.05 else "†" if w_p < 0.1 else ""
            if any(x in col for x in ['ssi','sk3','sn3','ccl']):
                print(f"  {col:15s} {phase_name:5s}: {mean_pct:+6.1f}%, {n_dec}/{len(ssw_vals)} dec, P={w_p:.4f} {sig}")

results['phase3_phase_resolved'] = phase3_results

# ── Phase 4: European Alps Multi-Country ──
print("\n" + "="*70)
print("PHASE 4: EUROPEAN ALPS MULTI-COUNTRY REPLICATION")
print("="*70)

try:
    alps = pd.read_csv('data/cryosphere/european_alps/data_dmax.csv', sep=';')
    alps['date'] = pd.to_datetime(alps['date'])
    alps = alps.dropna(subset=['country', 'dangerLevelMax'])
    print(f"Alps data: {len(alps)} rows, {alps['date'].min().date()} to {alps['date'].max().date()}")
    print(f"Countries: {alps['country'].value_counts().to_dict()}")
    
    alps_range = [d for d in ssw_dates if alps['date'].min().date() <= d <= alps['date'].max().date()]
    print(f"SSW events in Alps range: {len(alps_range)}")
    
    phase4_results = {}
    for country in sorted(alps['country'].unique()):
        if pd.isna(country):
            continue
        cdata = alps[alps['country'] == country]
        
        ssw_dangers = []
        ctrl_dangers = []
        events_detail = []
        
        for ssw_date in alps_range:
            ssw_dt = pd.Timestamp(ssw_date)
            
            # Post-SSW 15 days
            post = cdata[(cdata['date'] >= ssw_dt) & 
                         (cdata['date'] < ssw_dt + pd.Timedelta(days=15))]
            
            # Control: same calendar dates, other years
            ctrl_vals = []
            for yr_offset in [-2, -1, 1, 2]:
                try:
                    ctrl_dt = ssw_dt.replace(year=ssw_dt.year + yr_offset)
                    ctrl = cdata[(cdata['date'] >= ctrl_dt) & 
                                 (cdata['date'] < ctrl_dt + pd.Timedelta(days=15))]
                    if len(ctrl) > 3:
                        ctrl_vals.append(ctrl['dangerLevelMax'].mean())
                except:
                    pass
            
            if len(post) > 3 and ctrl_vals:
                ssw_mean = post['dangerLevelMax'].mean()
                ctrl_mean = np.mean(ctrl_vals)
                ssw_dangers.append(ssw_mean)
                ctrl_dangers.append(ctrl_mean)
                pct = (ssw_mean - ctrl_mean) / ctrl_mean * 100 if ctrl_mean > 0 else 0
                events_detail.append({
                    'date': str(ssw_date),
                    'ssw_danger': round(ssw_mean, 2),
                    'ctrl_danger': round(ctrl_mean, 2),
                    'pct_change': round(pct, 1)
                })
        
        if len(ssw_dangers) >= 2:
            ssw_arr = np.array(ssw_dangers)
            ctrl_arr = np.array(ctrl_dangers)
            diffs = ssw_arr - ctrl_arr
            n_dec = sum(1 for d in diffs if d < 0)
            mean_pct = np.mean([(s-c)/c*100 if c>0 else 0 for s,c in zip(ssw_arr, ctrl_arr)])
            
            try:
                mw_stat, mw_p = stats.mannwhitneyu(ssw_arr, ctrl_arr, alternative='two-sided')
            except:
                mw_p = 1.0
            
            phase4_results[country] = {
                'n_events': len(ssw_dangers),
                'n_decrease': n_dec,
                'mean_pct_change': round(mean_pct, 1),
                'mw_p': round(mw_p, 4),
                'events': events_detail
            }
            
            dir_str = "↓" if n_dec > len(ssw_dangers)/2 else "↑"
            print(f"  {country}: {mean_pct:+6.1f}%, {n_dec}/{len(ssw_dangers)} decrease {dir_str}, P={mw_p:.4f}")
    
    results['phase4_european_alps'] = phase4_results

except Exception as e:
    print(f"  Alps analysis error: {e}")
    results['phase4_european_alps'] = {'error': str(e)}

# ── Phase 5: ENSO/QBO Conditioning ──
print("\n" + "="*70)
print("PHASE 5: ENSO/QBO CONDITIONING OF SSW-AVALANCHE LINK")
print("="*70)

try:
    butler = pd.read_csv('data/processed/atmospheric/butler_ssw_compendium_era5.csv')
    butler['onset_date'] = pd.to_datetime(butler['onset_date'])
    
    # Load panel data for Swiss avalanches
    panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
    panel['date'] = pd.to_datetime(panel['date'])
    
    phase5_results = {}
    
    # For each ENSO phase
    for enso_phase in ['L', 'E', 'N']:
        phase_events = butler[butler['enso_phase'] == enso_phase]
        if len(phase_events) == 0:
            continue
        
        pre_rates = []
        post_rates = []
        
        for _, row in phase_events.iterrows():
            ssw_dt = row['onset_date']
            
            # Find in panel
            pre = panel[(panel['date'] >= ssw_dt - pd.Timedelta(days=30)) & 
                        (panel['date'] < ssw_dt)]
            post = panel[(panel['date'] >= ssw_dt) & 
                         (panel['date'] < ssw_dt + pd.Timedelta(days=15))]
            
            aval_col = 'dry_natural_size_1234' if 'dry_natural_size_1234' in panel.columns else 'aai_all_dry'
            
            if len(pre) > 5 and len(post) > 5:
                pre_rates.append(pre[aval_col].mean())
                post_rates.append(post[aval_col].mean())
        
        if pre_rates:
            pre_arr = np.array(pre_rates)
            post_arr = np.array(post_rates)
            mean_change = np.mean((post_arr - pre_arr) / np.where(pre_arr != 0, pre_arr, 1) * 100)
            n_dec = sum(1 for p, q in zip(post_arr, pre_arr) if p < q)
            
            phase5_results[f'ENSO_{enso_phase}'] = {
                'n_events': len(pre_rates),
                'mean_pct_change': round(mean_change, 1),
                'n_decrease': n_dec,
                'label': {'L': 'La Niña', 'E': 'El Niño', 'N': 'Neutral'}[enso_phase]
            }
            print(f"  ENSO {enso_phase} ({phase5_results[f'ENSO_{enso_phase}']['label']}): "
                  f"{mean_change:+.1f}%, {n_dec}/{len(pre_rates)} decrease")
    
    # For each QBO phase
    for qbo_phase in ['W', 'E']:
        phase_events = butler[butler['qbo_phase'] == qbo_phase]
        if len(phase_events) == 0:
            continue
        
        pre_rates = []
        post_rates = []
        
        for _, row in phase_events.iterrows():
            ssw_dt = row['onset_date']
            pre = panel[(panel['date'] >= ssw_dt - pd.Timedelta(days=30)) & 
                        (panel['date'] < ssw_dt)]
            post = panel[(panel['date'] >= ssw_dt) & 
                         (panel['date'] < ssw_dt + pd.Timedelta(days=15))]
            
            aval_col = 'dry_natural_size_1234' if 'dry_natural_size_1234' in panel.columns else 'aai_all_dry'
            
            if len(pre) > 5 and len(post) > 5:
                pre_rates.append(pre[aval_col].mean())
                post_rates.append(post[aval_col].mean())
        
        if pre_rates:
            pre_arr = np.array(pre_rates)
            post_arr = np.array(post_rates)
            mean_change = np.mean((post_arr - pre_arr) / np.where(pre_arr != 0, pre_arr, 1) * 100)
            n_dec = sum(1 for p, q in zip(post_arr, pre_arr) if p < q)
            
            phase5_results[f'QBO_{qbo_phase}'] = {
                'n_events': len(pre_rates),
                'mean_pct_change': round(mean_change, 1),
                'n_decrease': n_dec,
                'label': {'W': 'Westerly', 'E': 'Easterly'}[qbo_phase]
            }
            print(f"  QBO {qbo_phase} ({phase5_results[f'QBO_{qbo_phase}']['label']}): "
                  f"{mean_change:+.1f}%, {n_dec}/{len(pre_rates)} decrease")
    
    results['phase5_enso_qbo'] = phase5_results

except Exception as e:
    print(f"  ENSO/QBO error: {e}")
    import traceback; traceback.print_exc()
    results['phase5_enso_qbo'] = {'error': str(e)}

# ── Phase 6: Danger Level Integration ──
print("\n" + "="*70)
print("PHASE 6: SWISS DANGER DESCRIPTIONS ANALYSIS")
print("="*70)

try:
    danger = pd.read_csv('data/cryosphere/swiss_snowpack/danger_descriptions_2012_2020.csv', 
                          encoding='latin-1')
    danger['date'] = pd.to_datetime(danger['validFromDate'])
    print(f"Danger descriptions: {len(danger)} rows")
    print(f"Problems: {danger['problem'].value_counts().to_dict() if 'problem' in danger.columns else 'N/A'}")
    
    # Analyze dry problem danger levels around SSW
    dry = danger[danger['problem'].str.lower().str.contains('dry', na=False)] if 'problem' in danger.columns else danger
    
    ssw_in_danger = [d for d in ssw_dates if danger['date'].min().date() <= d <= danger['date'].max().date()]
    print(f"SSW events in danger range: {len(ssw_in_danger)}")
    
    ssw_levels = []
    ctrl_levels = []
    
    for ssw_date in ssw_in_danger:
        ssw_dt = pd.Timestamp(ssw_date)
        post = dry[(dry['date'] >= ssw_dt) & (dry['date'] < ssw_dt + pd.Timedelta(days=15))]
        
        ctrl_vals = []
        for yr_offset in [-2, -1, 1, 2]:
            try:
                ctrl_dt = ssw_dt.replace(year=ssw_dt.year + yr_offset)
                ctrl = dry[(dry['date'] >= ctrl_dt) & (dry['date'] < ctrl_dt + pd.Timedelta(days=15))]
                if len(ctrl) > 0 and 'dangerlevel' in ctrl.columns:
                    ctrl_vals.append(ctrl['dangerlevel'].mean())
            except:
                pass
        
        if len(post) > 0 and ctrl_vals and 'dangerlevel' in post.columns:
            ssw_levels.append(post['dangerlevel'].mean())
            ctrl_levels.append(np.mean(ctrl_vals))
    
    if ssw_levels:
        ssw_arr = np.array(ssw_levels)
        ctrl_arr = np.array(ctrl_levels)
        mean_pct = np.mean((ssw_arr - ctrl_arr) / ctrl_arr * 100)
        n_dec = sum(1 for s, c in zip(ssw_arr, ctrl_arr) if s < c)
        
        try:
            w_stat, w_p = stats.wilcoxon(ssw_arr, ctrl_arr)
        except:
            w_p = 1.0
        
        results['phase6_danger'] = {
            'n_events': len(ssw_levels),
            'mean_danger_ssw': round(np.mean(ssw_arr), 2),
            'mean_danger_ctrl': round(np.mean(ctrl_arr), 2),
            'pct_change': round(mean_pct, 1),
            'n_decrease': n_dec,
            'wilcoxon_p': round(w_p, 4)
        }
        print(f"  SSW danger: {np.mean(ssw_arr):.2f} vs ctrl: {np.mean(ctrl_arr):.2f} "
              f"({mean_pct:+.1f}%), {n_dec}/{len(ssw_levels)} decrease, P={w_p:.4f}")

except Exception as e:
    print(f"  Danger analysis error: {e}")
    import traceback; traceback.print_exc()

# ── Phase 7: Forecast Data Analysis ──
print("\n" + "="*70)
print("PHASE 7: SNOWPACK FORECAST DATA ANALYSIS")
print("="*70)

try:
    forecast = pd.read_csv('data/cryosphere/swiss_snowpack/data_rf1_forecast.csv', nrows=100)
    print(f"Forecast columns: {list(forecast.columns)[:20]}")
    print(f"Forecast shape (sample): {forecast.shape}")
except Exception as e:
    print(f"  Forecast data error: {e}")

# ── Phase 8: Comprehensive Summary ──
print("\n" + "="*70)
print("PHASE 8: KEY FINDINGS SUMMARY")
print("="*70)

print("\n--- SNOWPACK Stability (Pre vs Post SSW) ---")
for col, r in sorted(phase1_results.items(), key=lambda x: x[1]['wilcoxon_p']):
    if r['wilcoxon_p'] < 0.1:
        sig = "***" if r['wilcoxon_p'] < 0.001 else "**" if r['wilcoxon_p'] < 0.01 else "*" if r['wilcoxon_p'] < 0.05 else "†"
        print(f"  {col:20s}: {r['mean_pct_change']:+6.1f}%, "
              f"{r['n_decrease']}/{r['n_events']} dec, "
              f"P={r['wilcoxon_p']:.4f} {sig}, d={r['cohens_d']:+.3f}")

# Save results
with open('data/results/r16_snowpack_mechanism.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to data/results/r16_snowpack_mechanism.json")
print("DONE")
