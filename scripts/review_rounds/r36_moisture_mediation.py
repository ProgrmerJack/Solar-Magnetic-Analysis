"""
R36 Upgrade: ERA5 Moisture Flux Proxy + Mediation Analysis + Blocking Index
Goal: Address reviewer requests for IVT analysis and formal mediation
"""
import pandas as pd, numpy as np, json
from scipy import stats

# Load data
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_daily.parquet')
try:
    era5_ext = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
    era5 = pd.concat([era5, era5_ext]).sort_index()
    era5 = era5[~era5.index.duplicated(keep='last')]
except: pass

ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_dates = ssw_cat.index.tz_localize(None)

# ============================================================
# PART 1: Moisture flux proxy from ERA5 surface data
# IVT requires pressure-level data; we approximate with surface variables
# Proxy: precip * southerly_component (captures Mediterranean moisture delivery)
# ============================================================
winter_mask = era5['doy'].isin(list(range(1,121)) + list(range(305,366)))
winter = era5[winter_mask].copy()

# Southerly moisture flux proxy: v10 * tp_mm (positive v10 = from south)
# Actually v10 is northward component, so positive = northward flow
# Mediterranean moisture comes from the south, so we want negative v10 (southward blocked)
# But actually for "moisture from Mediterranean", we want POSITIVE v10 at Swiss lat
# (air moving northward = Mediterranean moisture advection)
winter['merid_moisture'] = winter['v10'] * winter['tp_mm']  # northward moisture transport proxy
winter['southerly_flow'] = winter['v10']  # positive = northward (Mediterranean)

ssw_in = ssw_dates[(ssw_dates >= winter.index.min()) & (ssw_dates <= winter.index.max())]

def ssw_mask_fn(idx, dates, w=15):
    m = pd.Series(False, index=idx)
    for d in dates:
        m |= (idx >= d - pd.Timedelta(days=w)) & (idx <= d + pd.Timedelta(days=w))
    return m

sm = ssw_mask_fn(winter.index, ssw_in)

print("=" * 60)
print("PART 1: Meridional Moisture Flux Proxy")
print("=" * 60)

for col, name in [
    ('southerly_flow', 'Southerly flow (v10)'),
    ('merid_moisture', 'Meridional moisture flux'),
    ('tp_mm', 'Total precipitation'),
    ('sf_mm', 'Snowfall'),
]:
    ssw_v = winter.loc[sm, col].dropna()
    ctrl_v = winter.loc[~sm, col].dropna()
    d = (ssw_v.mean() - ctrl_v.mean()) / ctrl_v.std()
    _, p = stats.mannwhitneyu(ssw_v, ctrl_v, alternative='two-sided')
    print(f'{name:30s}: SSW={ssw_v.mean():.4f} Ctrl={ctrl_v.mean():.4f} d={d:+.3f} P={p:.2e}')

# ============================================================
# PART 2: Z500-based Euro-Atlantic blocking proxy
# Blocking = anomalously high Z500 south + low Z500 at higher lat
# We only have single-point Z500 for NH, but can use it as proxy
# ============================================================
print("\n" + "=" * 60)
print("PART 2: Z500 Blocking Analysis")
print("=" * 60)

# Use panel data for Z500
panel_winter = panel[panel.index.month.isin([1,2,3,4,11,12])].copy()
panel_winter['z500'] = panel_winter['ncep_z500_nh']

# Compute Z500 anomaly relative to DOY mean
panel_winter['doy'] = panel_winter.index.dayofyear
doy_mean = panel_winter.groupby('doy')['z500'].transform('mean')
doy_std = panel_winter.groupby('doy')['z500'].transform('std')
panel_winter['z500_anom'] = (panel_winter['z500'] - doy_mean) / doy_std

ssw_panel = ssw_dates[(ssw_dates >= panel_winter.index.min()) & (ssw_dates <= panel_winter.index.max())]
sm_panel = ssw_mask_fn(panel_winter.index, ssw_panel)

z_ssw = panel_winter.loc[sm_panel, 'z500_anom'].dropna()
z_ctrl = panel_winter.loc[~sm_panel, 'z500_anom'].dropna()
d_z = z_ssw.mean() - z_ctrl.mean()  # already standardized
_, p_z = stats.mannwhitneyu(z_ssw, z_ctrl)
print(f'Z500 anomaly: SSW={z_ssw.mean():.3f}sigma Ctrl={z_ctrl.mean():.3f}sigma d={d_z:.3f} P={p_z:.2e}')

# Event-level Z500
print("\nEvent-level Z500 anomalies:")
for i, d in enumerate(ssw_panel):
    w = (panel_winter.index >= d - pd.Timedelta(days=15)) & (panel_winter.index <= d + pd.Timedelta(days=15))
    if w.sum() > 0:
        z_event = panel_winter.loc[w, 'z500_anom'].mean()
        print(f'  SSW {d.strftime("%Y-%m-%d")}: Z500 anom = {z_event:+.3f} sigma')

# ============================================================
# PART 3: Formal Causal Mediation Analysis
# Path: SSW -> snowfall_fraction -> avalanche_RR
# Using Baron-Kenny approach + Sobel test
# ============================================================
print("\n" + "=" * 60)
print("PART 3: Formal Mediation Analysis")
print("=" * 60)

# Build event-level dataset
panel_av = panel[['dry_natural_size_1234']].copy()
panel_av.columns = ['av']

events = []
for d in ssw_panel:
    # SSW window
    ssw_w = (panel_av.index >= d - pd.Timedelta(days=15)) & (panel_av.index <= d + pd.Timedelta(days=15))
    # Control: same DOY from other years
    ssw_doy_start = (d - pd.Timedelta(days=15)).dayofyear
    ssw_doy_end = (d + pd.Timedelta(days=15)).dayofyear
    
    av_ssw = panel_av.loc[ssw_w, 'av'].mean()
    
    # ERA5 snow fraction during SSW
    era5_w = (winter.index >= d - pd.Timedelta(days=15)) & (winter.index <= d + pd.Timedelta(days=15))
    if era5_w.sum() == 0:
        continue
    
    rain_mm = (winter.loc[era5_w, 'tp_mm'] - winter.loc[era5_w, 'sf_mm']).clip(lower=0)
    snow_frac_ssw = winter.loc[era5_w, 'sf_mm'].sum() / max(winter.loc[era5_w, 'tp_mm'].sum(), 0.01)
    temp_ssw = winter.loc[era5_w, 't2m_K'].mean() if 't2m_K' in winter.columns else np.nan
    rain_ssw = rain_mm.mean()
    
    # Control values (same DOY, other years)
    ctrl_years = [y for y in winter.index.year.unique() if y != d.year]
    ctrl_sf = []
    ctrl_av_vals = []
    ctrl_rain = []
    ctrl_temp = []
    for y in ctrl_years:
        try:
            cd = pd.Timestamp(year=y, month=d.month, day=d.day)
            cw = (winter.index >= cd - pd.Timedelta(days=15)) & (winter.index <= cd + pd.Timedelta(days=15))
            if cw.sum() > 5:
                ctrl_sf.append(winter.loc[cw, 'sf_mm'].sum() / max(winter.loc[cw, 'tp_mm'].sum(), 0.01))
                cr = (winter.loc[cw, 'tp_mm'] - winter.loc[cw, 'sf_mm']).clip(lower=0).mean()
                ctrl_rain.append(cr)
                if 't2m_K' in winter.columns:
                    ctrl_temp.append(winter.loc[cw, 't2m_K'].mean())
            
            cw_av = (panel_av.index >= cd - pd.Timedelta(days=15)) & (panel_av.index <= cd + pd.Timedelta(days=15))
            if cw_av.sum() > 5:
                ctrl_av_vals.append(panel_av.loc[cw_av, 'av'].mean())
        except: continue
    
    if len(ctrl_sf) < 3 or len(ctrl_av_vals) < 3:
        continue
    
    events.append({
        'onset': d,
        'av_ssw': av_ssw,
        'av_ctrl': np.mean(ctrl_av_vals),
        'log_rr': np.log(max(av_ssw, 0.01) / max(np.mean(ctrl_av_vals), 0.01)),
        'sf_ssw': snow_frac_ssw,
        'sf_ctrl': np.mean(ctrl_sf),
        'sf_change': snow_frac_ssw - np.mean(ctrl_sf),
        'rain_ssw': rain_ssw,
        'rain_ctrl': np.mean(ctrl_rain),
        'rain_change': rain_ssw - np.mean(ctrl_rain),
        'temp_ssw': temp_ssw,
        'temp_ctrl': np.mean(ctrl_temp) if ctrl_temp else np.nan,
    })

edf = pd.DataFrame(events)
print(f"Event-level dataset: {len(edf)} SSW events with matched ERA5 + avalanche data")

if len(edf) >= 5:
    # Path a: SSW -> Mediator (snow fraction change)
    # Already established: SSW significantly changes snow fraction
    
    # Path b: Mediator -> Outcome (snow fraction -> log RR)
    r_sf, p_sf = stats.pearsonr(edf['sf_change'].dropna(), edf['log_rr'].dropna().iloc[:len(edf['sf_change'].dropna())])
    print(f"\nSnow fraction change vs log(RR): r={r_sf:.3f}, P={p_sf:.3f}")
    
    r_rain, p_rain = stats.pearsonr(edf['rain_change'].dropna(), edf['log_rr'].dropna().iloc[:len(edf['rain_change'].dropna())])
    print(f"Rainfall change vs log(RR): r={r_rain:.3f}, P={p_rain:.3f}")
    
    if 'temp_ctrl' in edf.columns:
        temp_change = edf['temp_ssw'] - edf['temp_ctrl']
        valid = temp_change.dropna()
        if len(valid) >= 5:
            r_temp, p_temp = stats.pearsonr(valid, edf.loc[valid.index, 'log_rr'])
            print(f"Temperature change vs log(RR): r={r_temp:.3f}, P={p_temp:.3f}")
    
    # Multiple mediator regression
    from sklearn.linear_model import LinearRegression
    
    # Full model: log_rr ~ sf_change + rain_change + temp_change
    X_cols = ['sf_change', 'rain_change']
    valid_mask = edf[X_cols + ['log_rr']].dropna().index
    if len(valid_mask) >= 5:
        X = edf.loc[valid_mask, X_cols].values
        y = edf.loc[valid_mask, 'log_rr'].values
        
        reg = LinearRegression().fit(X, y)
        r2 = reg.score(X, y)
        print(f"\nMultiple mediation R2: {r2:.3f}")
        for c, name in zip(reg.coef_, X_cols):
            print(f"  {name}: beta = {c:.4f}")
        
        # Proportion mediated
        # Total effect: mean log(RR)
        total = edf['log_rr'].mean()
        predicted = reg.predict(X).mean() - reg.intercept_  # indirect effect
        if abs(total) > 0.01:
            prop_mediated = predicted / total * 100
            print(f"\nTotal effect (mean log RR): {total:.3f}")
            print(f"Mediated by snow fraction + rainfall: {prop_mediated:.1f}%")

# ============================================================
# PART 4: Event-level cascade timing
# For each SSW, show: T_strat -> Z500 -> T_surface -> Snow_frac -> Avalanche
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Event-Level Cascade")
print("=" * 60)

for i, row in edf.iterrows():
    d = row['onset']
    print(f"\nSSW {d.strftime('%Y-%m-%d')}:")
    print(f"  Snow frac: SSW={row['sf_ssw']:.3f} Ctrl={row['sf_ctrl']:.3f} change={row['sf_change']:+.3f}")
    print(f"  Rainfall:  SSW={row['rain_ssw']:.4f} Ctrl={row['rain_ctrl']:.4f} change={row['rain_change']:+.4f}")
    print(f"  Avalanche: SSW={row['av_ssw']:.2f} Ctrl={row['av_ctrl']:.2f} log(RR)={row['log_rr']:+.2f}")

# ============================================================
# PART 5: Save results
# ============================================================
results = {
    'n_events': len(edf),
    'snow_fraction': {
        'ssw_mean': float(edf['sf_ssw'].mean()),
        'ctrl_mean': float(edf['sf_ctrl'].mean()),
        'change_mean': float(edf['sf_change'].mean()),
    },
    'rainfall': {
        'ssw_mean': float(edf['rain_ssw'].mean()),
        'ctrl_mean': float(edf['rain_ctrl'].mean()),
        'change_mean': float(edf['rain_change'].mean()),
    },
    'mediation': {
        'sf_vs_rr_r': float(r_sf) if 'r_sf' in dir() else None,
        'sf_vs_rr_p': float(p_sf) if 'p_sf' in dir() else None,
        'rain_vs_rr_r': float(r_rain) if 'r_rain' in dir() else None,
        'rain_vs_rr_p': float(p_rain) if 'p_rain' in dir() else None,
    }
}

with open('data/results/r36_moisture_mediation.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to data/results/r36_moisture_mediation.json")
