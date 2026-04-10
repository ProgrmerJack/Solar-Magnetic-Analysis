"""
Per-event dose-response: does sintering enhancement predict avalanche reduction?
This tests whether SSW events that produce more surface warming (and thus more sintering)
also produce larger avalanche reductions.
"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.join(os.path.dirname(__file__), '..')

# Load sintering results
with open(os.path.join(BASE, 'data', 'results', 'sintering_extended.json')) as f:
    sinter = json.load(f)

# Load Swiss panel data
panel = pd.read_parquet(os.path.join(BASE, 'data', 'processed', 'analysis_panel_v2.parquet'))
panel.index = pd.to_datetime(panel.index).tz_localize(None)
aval_col = [c for c in panel.columns if 'dry' in c.lower() and 'nat' in c.lower()][0]
panel_clean = panel[[aval_col]].dropna()

# Winter baseline rate per season
panel_clean['winter'] = panel_clean.index.year + (panel_clean.index.month >= 10).astype(int)
winter_baseline = panel_clean.groupby('winter')[aval_col].mean()

# For each SSW event, compute avalanche rate ratio
results = []
for event in sinter['per_event']:
    ssw_date = pd.Timestamp(event['ssw_date'])
    winter = ssw_date.year + (ssw_date.month >= 10)
    
    # SSW window: -5 to +15 days
    window_start = ssw_date - pd.Timedelta(days=5)
    window_end = ssw_date + pd.Timedelta(days=15)
    
    ssw_aval = panel_clean.loc[window_start:window_end, aval_col].dropna()
    
    if len(ssw_aval) < 5:
        print(f"  {event['ssw_date']}: insufficient avalanche data ({len(ssw_aval)} days)")
        continue
    
    # Baseline: same winter, outside SSW windows
    ssw_mask = (panel_clean.index >= window_start) & (panel_clean.index <= window_end)
    same_winter = panel_clean['winter'] == winter
    ctrl_aval = panel_clean.loc[same_winter & ~ssw_mask, aval_col].dropna()
    
    if len(ctrl_aval) < 10:
        print(f"  {event['ssw_date']}: insufficient control data ({len(ctrl_aval)} days)")
        continue
    
    rr = ssw_aval.mean() / ctrl_aval.mean() if ctrl_aval.mean() > 0 else np.nan
    
    results.append({
        'ssw_date': event['ssw_date'],
        'sintering_pct': event['sintering_enhancement_pct'],
        'delta_T': event['delta_T_K'],
        'avalanche_rr': round(rr, 3),
        'ssw_rate': round(ssw_aval.mean(), 2),
        'ctrl_rate': round(ctrl_aval.mean(), 2),
        'n_ssw_days': len(ssw_aval),
        'n_ctrl_days': len(ctrl_aval)
    })

print(f"\n{'='*60}")
print(f"PER-EVENT DOSE-RESPONSE: Sintering vs Avalanche Reduction")
print(f"{'='*60}")
print(f"{'SSW Date':<14} {'ΔT':>6} {'Sinter%':>8} {'Aval RR':>8} {'SSW rate':>9} {'Ctrl rate':>10}")
print("-" * 60)
for r in sorted(results, key=lambda x: x['ssw_date']):
    print(f"{r['ssw_date']:<14} {r['delta_T']:>+6.2f} {r['sintering_pct']:>+7.1f}% {r['avalanche_rr']:>8.3f} {r['ssw_rate']:>9.2f} {r['ctrl_rate']:>10.2f}")

# Correlations
sinter_vals = np.array([r['sintering_pct'] for r in results])
rr_vals = np.array([r['avalanche_rr'] for r in results])
dt_vals = np.array([r['delta_T'] for r in results])

# Higher sintering → lower RR (expect negative correlation)
r_sinter_rr, p_sinter_rr = stats.spearmanr(sinter_vals, rr_vals)
r_dt_rr, p_dt_rr = stats.spearmanr(dt_vals, rr_vals)

print(f"\nCorrelation: sintering_enhancement vs avalanche_RR")
print(f"  Spearman r = {r_sinter_rr:.3f}, P = {p_sinter_rr:.4f}")
print(f"Correlation: delta_T vs avalanche_RR")
print(f"  Spearman r = {r_dt_rr:.3f}, P = {p_dt_rr:.4f}")

# The key question: do warming events suppress avalanches more?
warming = [(r['sintering_pct'], r['avalanche_rr']) for r in results if r['sintering_pct'] > 0]
cooling = [(r['sintering_pct'], r['avalanche_rr']) for r in results if r['sintering_pct'] <= 0]

if warming and cooling:
    warm_rr = [w[1] for w in warming]
    cool_rr = [c[1] for c in cooling]
    mw_stat, mw_p = stats.mannwhitneyu(warm_rr, cool_rr, alternative='less')
    print(f"\nWarming SSWs (n={len(warming)}): mean RR = {np.mean(warm_rr):.3f}")
    print(f"Cooling SSWs (n={len(cooling)}): mean RR = {np.mean(cool_rr):.3f}")
    print(f"Mann-Whitney (warming < cooling): U={mw_stat}, P={mw_p:.4f}")

# Summary
n_decrease = sum(1 for r in results if r['avalanche_rr'] < 1)
print(f"\nAvalanche decrease: {n_decrease}/{len(results)} events")
print(f"Mean RR: {np.mean(rr_vals):.3f}")

# Save
output = {
    'n_events': len(results),
    'dose_response_sinter_rr': {'r': round(r_sinter_rr, 3), 'p': round(p_sinter_rr, 4)},
    'dose_response_deltaT_rr': {'r': round(r_dt_rr, 3), 'p': round(p_dt_rr, 4)},
    'per_event': results
}
out_path = os.path.join(BASE, 'data', 'results', 'sintering_dose_response.json')
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {out_path}")
