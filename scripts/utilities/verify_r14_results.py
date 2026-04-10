"""Verify key R14 results referenced in the paper."""
import pandas as pd
import numpy as np
import json

# Load datasets
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_extended.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_cat.index = ssw_cat.index.tz_localize(None)

print('Panel shape:', panel.shape)
print('ERA5 shape:', era5.shape, '  Range:', era5.index.min(), 'to', era5.index.max())
print('SSW events:', len(ssw_cat))
print('SSW types:', ssw_cat['type'].value_counts().to_dict())
print()

# Load R14 mechanism breakthrough results
with open('data/results/r14_mechanism_breakthrough.json') as f:
    r14 = json.load(f)

# Phase 4: Propagation
print('=== PHASE 4: PROPAGATION (displacement) ===')
if 'phase4_propagation' in r14:
    p4 = r14['phase4_propagation']
    if 'displacement' in p4:
        disp = p4['displacement']
        if 'levels' in disp:
            for lev in disp['levels']:
                print("  %s: peak lag=%s, anom=%.1f K, P=%.6f" % (
                    lev.get('level','?'), lev.get('peak_lag','?'),
                    lev.get('peak_anom_K', 0), lev.get('peak_P', 1)))
        elif isinstance(disp, dict):
            for k, v in disp.items():
                if isinstance(v, dict):
                    print("  %s: %s" % (k, v))

# Phase 7: NH circulation
print()
print('=== PHASE 7: NH CIRCULATION ===')
if 'phase7_nh_circulation' in r14:
    p7 = r14['phase7_nh_circulation']
    for key, val in p7.items():
        if isinstance(val, dict):
            print("  %s: anom=%s, P=%s" % (key, val.get('anomaly','?'), val.get('P','?')))

# Phase 1: Composites
print()
print('=== PHASE 1: ERA5 COMPOSITES ===')
if 'phase1_composites' in r14:
    p1 = r14['phase1_composites']
    for ssw_type, data in p1.items():
        if isinstance(data, dict):
            print("  %s:" % ssw_type)
            for var, stats in data.items():
                if isinstance(stats, dict):
                    print("    %s: anom=%s, P=%s" % (var, stats.get('anomaly','?'), stats.get('P','?')))

# Check mediation results from console output (re-derive key numbers)
print()
print('=== MEDIATION VERIFICATION ===')
# Wind mediation: SSW -> wind anomaly -> avalanches
if 'wind_speed_anom' in era5.columns:
    ssw_mask = panel.get('ssw_active', pd.Series(False, index=panel.index))
    # Merge ERA5 with panel on date
    merged = panel.join(era5[['wind_speed_anom']], how='inner')
    if 'ssw_active' in merged.columns and 'wind_speed_anom' in merged.columns:
        ssw_days = merged[merged['ssw_active'] == True]
        ctrl_days = merged[merged['ssw_active'] == False]
        wind_diff = ssw_days['wind_speed_anom'].mean() - ctrl_days['wind_speed_anom'].mean()
        from scipy import stats
        t_stat, p_val = stats.ttest_ind(ssw_days['wind_speed_anom'].dropna(), 
                                         ctrl_days['wind_speed_anom'].dropna())
        print("  Wind speed anomaly SSW vs ctrl: diff=%.4f, P=%.4f" % (wind_diff, p_val))

print()
print('=== PAPER CLAIMS VERIFICATION ===')

# 1. Sign test: 14/16 decrease
if 'ssw_event_diffs' in panel.columns or True:
    # Load sintering results to check
    try:
        with open('data/results/sintering_extended.json') as f:
            sint = json.load(f)
        print("  Sintering: n=%d events" % sint.get('n_events', '?'))
        print("  Displacement: mean=%.1f%%, P=%s" % (
            sint.get('displacement',{}).get('mean_pct',0), 
            sint.get('displacement',{}).get('P','?')))
    except Exception as e:
        print("  Could not load sintering results:", e)

# 2. Check Norwegian data exists
try:
    nve = pd.read_parquet('data/processed/avalanche/norway_nve_danger.parquet')
    print("  Norway NVE: %d rows" % len(nve))
except:
    try:
        nve = pd.read_csv('data/processed/avalanche/norway_nve_danger.csv')
        print("  Norway NVE: %d rows" % len(nve))
    except Exception as e:
        print("  Norway NVE: could not load -", e)

# 3. Check Utah data
try:
    utah = pd.read_parquet('data/processed/avalanche/utah_daily_counts.parquet')
    print("  Utah: %d rows" % len(utah))
except:
    try:
        utah = pd.read_csv('data/processed/avalanche/utah_daily_counts.csv')
        print("  Utah: %d rows" % len(utah))
    except Exception as e:
        print("  Utah: could not load -", e)

print()
print("Verification complete.")
