"""
Window-sensitivity event-flip analysis.

Identifies which specific SSW events change sign (RR crosses 1.0)
at different window widths (±5, ±7, ±10, ±15 days), explaining the
non-monotonic P-value pattern in the sign test.

Output: data/results/window_sensitivity_events.json
"""
import pandas as pd
import numpy as np
import json, os, sys

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, BASE)

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
windows = [5, 7, 10, 15]

def compute_rr(onset, half_window, panel_df, outcome):
    """Compute RR for a single event at a given window width."""
    start = onset - pd.Timedelta(days=half_window)
    end = onset + pd.Timedelta(days=half_window)
    event_days = panel_df.loc[start:end]
    obs = event_days[outcome].sum()

    # DOY-matched control from non-SSW winters
    doys = event_days.index.dayofyear.unique()
    ssw_mask = panel_df['ssw_within_15d'].fillna(False).astype(bool)
    ctrl_mask = ~ssw_mask & panel_df.index.dayofyear.isin(doys)
    ctrl_days = panel_df.loc[ctrl_mask]

    if len(ctrl_days) == 0 or len(event_days) == 0:
        return np.nan
    # Scale control to same number of days
    n_ssw_winters = 1
    n_ctrl_winters = ctrl_days.groupby(ctrl_days.index.year).ngroups
    if n_ctrl_winters == 0:
        return np.nan
    exp = ctrl_days[outcome].sum() * (len(event_days) / len(ctrl_days))
    if exp == 0:
        return np.nan
    return obs / exp

results = {}
for w in windows:
    rrs = []
    for onset in ssw_dates:
        rr = compute_rr(onset, w, panel, outcome_col)
        rrs.append({
            'onset': str(onset.date()),
            'window': f'±{w}d',
            'RR': round(rr, 4) if not np.isnan(rr) else None,
            'direction': 'decrease' if rr < 1 else 'increase' if rr > 1 else 'neutral'
        })
    results[f'pm{w}d'] = rrs

# Compare across windows: find events that flip
flips = []
for i, onset in enumerate(ssw_dates):
    event_rrs = {}
    for w in windows:
        key = f'pm{w}d'
        rr = results[key][i]['RR']
        event_rrs[f'±{w}d'] = rr
    
    # Check for sign changes
    signs = {k: ('decrease' if v < 1 else 'increase') for k, v in event_rrs.items() if v is not None}
    unique_signs = set(signs.values())
    if len(unique_signs) > 1:
        flips.append({
            'onset': str(onset.date()),
            'rr_by_window': event_rrs,
            'signs': signs,
            'flip_type': 'decrease→increase' if signs.get('±7d') == 'decrease' else 'increase→decrease'
        })

# Summary statistics
summary = {}
for w in windows:
    key = f'pm{w}d'
    n_decrease = sum(1 for r in results[key] if r['RR'] is not None and r['RR'] < 1)
    n_increase = sum(1 for r in results[key] if r['RR'] is not None and r['RR'] >= 1)
    n_total = sum(1 for r in results[key] if r['RR'] is not None)
    from scipy import stats
    p_sign = stats.binomtest(n_decrease, n_total, 0.5).pvalue
    summary[f'±{w}d'] = {
        'n_decrease': n_decrease,
        'n_increase': n_increase,
        'n_total': n_total,
        'P_sign': round(p_sign, 4)
    }

output = {
    'window_results': results,
    'flipping_events': flips,
    'summary': summary,
    'n_flipping': len(flips),
    'flipping_onsets': [f['onset'] for f in flips]
}

os.makedirs(f'{BASE}/data/results', exist_ok=True)
with open(f'{BASE}/data/results/window_sensitivity_events.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n=== Window Sensitivity Event-Flip Analysis ===")
print(f"\nSign-test summary by window width:")
for w, s in summary.items():
    print(f"  {w}: {s['n_decrease']}/{s['n_total']} decrease, P = {s['P_sign']}")

print(f"\n{len(flips)} events flip sign across windows:")
for flip in flips:
    print(f"\n  {flip['onset']}:")
    for w, rr in flip['rr_by_window'].items():
        sign = flip['signs'].get(w, '?')
        print(f"    {w}: RR = {rr:.3f} ({sign})")
