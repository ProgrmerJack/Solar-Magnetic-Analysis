"""
Bootstrap confidence intervals for Blinder-Oaxaca decomposition components.

Resamples the 16 SSW events (with replacement) 10,000 times and
re-computes the composition/within/interaction shares of the B-O
decomposition, yielding 95% percentile CIs on each component.

Output: data/results/bo_bootstrap_cis.json
"""
import pandas as pd
import numpy as np
import json, os, sys

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load panel
panel = pd.read_parquet(f'{BASE}/data/processed/analysis_panel_v2.parquet')
panel.index = pd.DatetimeIndex(panel.index)

# Load SSW catalog
ssw_cat = pd.read_parquet(f'{BASE}/data/processed/atmospheric/ssw_catalog.parquet')
if 'date' in ssw_cat.columns:
    ssw_dates = pd.DatetimeIndex(ssw_cat['date'])
elif 'onset' in ssw_cat.columns:
    ssw_dates = pd.DatetimeIndex(ssw_cat['onset'])
else:
    ssw_dates = pd.DatetimeIndex(ssw_cat.index)
if ssw_dates.tz is not None:
    ssw_dates = ssw_dates.tz_localize(None)
ssw_dates = ssw_dates[(ssw_dates >= '1998-11-01') & (ssw_dates <= '2019-06-01')]
ssw_dates = ssw_dates.sort_values()

outcome_col = 'dry_natural_size_1234'
ssw_mask = panel['ssw_within_15d'].fillna(False).astype(bool)

# Define weather regimes using ERA5 T2m and precip median splits
t2m_col = 'era5_t2m' if 'era5_t2m' in panel.columns else None
precip_col = 'era5_tp' if 'era5_tp' in panel.columns else None

# Find the right column names - merge ERA5 data into panel
era5_file = f'{BASE}/data/processed/era5_swiss_alps_extended.parquet'
if not os.path.exists(era5_file):
    era5_file = f'{BASE}/data/processed/era5_swiss_alps_daily.parquet'
era5 = pd.read_parquet(era5_file)
era5.index = pd.DatetimeIndex(era5.index)

# Merge ERA5 into panel
panel = panel.join(era5[['t2m_K', 'tp_mm']], how='left')
t2m_col = 't2m_K'
precip_col = 'tp_mm'

print(f"Using T2m: {t2m_col}, Precip: {precip_col}")

# Winter days only (Nov-Apr)
winter = panel[panel.index.month.isin([11, 12, 1, 2, 3, 4])].copy()
winter['ssw'] = ssw_mask.reindex(winter.index, fill_value=False)

# Define regimes by median split
t2m_med = winter[t2m_col].median()
precip_med = winter[precip_col].median()
winter['regime'] = 'warm-wet'
winter.loc[(winter[t2m_col] < t2m_med) & (winter[precip_col] < precip_med), 'regime'] = 'cold-dry'
winter.loc[(winter[t2m_col] < t2m_med) & (winter[precip_col] >= precip_med), 'regime'] = 'cold-wet'
winter.loc[(winter[t2m_col] >= t2m_med) & (winter[precip_col] < precip_med), 'regime'] = 'warm-dry'

def bo_decomposition(data, ssw_col='ssw', outcome='dry_natural_size_1234', regime_col='regime'):
    """Compute Blinder-Oaxaca decomposition of rate gap."""
    ssw_data = data[data[ssw_col] == True]
    ctrl_data = data[data[ssw_col] == False]
    
    if len(ssw_data) == 0 or len(ctrl_data) == 0:
        return None
    
    regimes = data[regime_col].unique()
    
    # Regime frequencies
    f_ssw = ssw_data[regime_col].value_counts(normalize=True)
    f_ctrl = ctrl_data[regime_col].value_counts(normalize=True)
    
    # Within-regime rates
    r_ssw = ssw_data.groupby(regime_col)[outcome].mean()
    r_ctrl = ctrl_data.groupby(regime_col)[outcome].mean()
    
    # Overall rates
    R_ssw = ssw_data[outcome].mean()
    R_ctrl = ctrl_data[outcome].mean()
    gap = R_ssw - R_ctrl
    
    if abs(gap) < 1e-10:
        return None
    
    # Three-fold decomposition
    composition = 0
    within = 0
    interaction = 0
    
    for reg in regimes:
        df = f_ssw.get(reg, 0) - f_ctrl.get(reg, 0)
        dr = r_ssw.get(reg, 0) - r_ctrl.get(reg, 0)
        composition += df * r_ctrl.get(reg, 0)
        within += f_ctrl.get(reg, 0) * dr
        interaction += df * dr
    
    return {
        'gap': gap,
        'composition': composition,
        'within': within,
        'interaction': interaction,
        'composition_pct': 100 * composition / gap if gap != 0 else 0,
        'within_pct': 100 * within / gap if gap != 0 else 0,
        'interaction_pct': 100 * interaction / gap if gap != 0 else 0,
        'R_ssw': R_ssw,
        'R_ctrl': R_ctrl
    }

# Point estimate
point = bo_decomposition(winter)
print(f"\n=== Point Estimate ===")
print(f"Gap: {point['gap']:.4f}")
print(f"Composition: {point['composition_pct']:.1f}%")
print(f"Within: {point['within_pct']:.1f}%")
print(f"Interaction: {point['interaction_pct']:.1f}%")

# Bootstrap: resample SSW event windows
np.random.seed(42)
n_boot = 10000
boot_comp = []
boot_within = []
boot_inter = []
boot_gap = []

# Get event-level windows for resampling
event_windows = []
for onset in ssw_dates:
    start = onset - pd.Timedelta(days=15)
    end = onset + pd.Timedelta(days=15)
    mask = (winter.index >= start) & (winter.index <= end)
    if mask.sum() > 0:
        event_windows.append(winter.index[mask])

n_events = len(event_windows)
ctrl_days = winter[winter['ssw'] == False]

for b in range(n_boot):
    # Resample events with replacement
    idx = np.random.choice(n_events, size=n_events, replace=True)
    boot_ssw_idx = pd.DatetimeIndex([])
    for i in idx:
        boot_ssw_idx = boot_ssw_idx.append(event_windows[i])
    boot_ssw_idx = boot_ssw_idx.unique()
    
    # Build bootstrap dataset
    boot_data = winter.copy()
    boot_data['ssw'] = boot_data.index.isin(boot_ssw_idx)
    
    result = bo_decomposition(boot_data)
    if result is not None:
        boot_comp.append(result['composition_pct'])
        boot_within.append(result['within_pct'])
        boot_inter.append(result['interaction_pct'])
        boot_gap.append(result['gap'])

boot_comp = np.array(boot_comp)
boot_within = np.array(boot_within)
boot_inter = np.array(boot_inter)
boot_gap = np.array(boot_gap)

def ci(arr, alpha=0.05):
    return [float(np.percentile(arr, 100*alpha/2)), float(np.percentile(arr, 100*(1-alpha/2)))]

output = {
    'point_estimate': {
        'gap': round(point['gap'], 4),
        'composition_pct': round(point['composition_pct'], 1),
        'within_pct': round(point['within_pct'], 1),
        'interaction_pct': round(point['interaction_pct'], 1),
    },
    'bootstrap': {
        'n_resamples': n_boot,
        'composition_pct_CI95': [round(x, 1) for x in ci(boot_comp)],
        'within_pct_CI95': [round(x, 1) for x in ci(boot_within)],
        'interaction_pct_CI95': [round(x, 1) for x in ci(boot_inter)],
        'gap_CI95': [round(x, 4) for x in ci(boot_gap)],
        'composition_pct_median': round(float(np.median(boot_comp)), 1),
        'within_pct_median': round(float(np.median(boot_within)), 1),
        'interaction_pct_median': round(float(np.median(boot_inter)), 1),
    }
}

os.makedirs(f'{BASE}/data/results', exist_ok=True)
with open(f'{BASE}/data/results/bo_bootstrap_cis.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n=== Bootstrap 95% CIs (n={n_boot}) ===")
print(f"Composition: {point['composition_pct']:.1f}% [{ci(boot_comp)[0]:.1f}, {ci(boot_comp)[1]:.1f}]")
print(f"Within:      {point['within_pct']:.1f}% [{ci(boot_within)[0]:.1f}, {ci(boot_within)[1]:.1f}]")
print(f"Interaction: {point['interaction_pct']:.1f}% [{ci(boot_inter)[0]:.1f}, {ci(boot_inter)[1]:.1f}]")
print(f"Gap:         {point['gap']:.4f} [{ci(boot_gap)[0]:.4f}, {ci(boot_gap)[1]:.4f}]")
