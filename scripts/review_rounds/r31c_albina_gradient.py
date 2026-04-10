"""R31c: ALBINA out-of-sample geographic gradient test.

Uses ALBINA data (Austria/Italy) from 2 SSW events (2023, 2024) to test the
EAWS geographic gradient prediction: Central-to-Northern Alps should show
neutral-to-increased danger.

Also decomposes by sub-region (AT-07=Tyrol, IT-32-BZ=South Tyrol, IT-32-TN=Trentino).
"""
import json, sys, os
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu, binomtest

sys.stdout.reconfigure(encoding='utf-8')

with open('data/cryosphere/albina_bulletins/albina_ssw_bulletins.json') as f:
    data = json.load(f)

print(f"Total records: {len(data)}")

# ── Overall analysis ──────────────────────────────────────────────────────────
print("\n" + "="*80)
print("ALBINA SSW ANALYSIS (Austria/Italy)")
print("="*80)

results = {'overall': {}, 'by_event': {}, 'by_region': {}, 'phase_resolved': {}}

# By event
for ssw_event in sorted(set(d['ssw_event'] for d in data)):
    records = [d for d in data if d['ssw_event'] == ssw_event]
    in_window = [d for d in records if d['in_ssw_window']]
    out_window = [d for d in records if not d['in_ssw_window']]
    
    in_danger = [d['danger_max'] for d in in_window if d['danger_max'] is not None]
    out_danger = [d['danger_max'] for d in out_window if d['danger_max'] is not None]
    
    if in_danger and out_danger:
        diff = np.mean(in_danger) - np.mean(out_danger)
        stat, p = mannwhitneyu(in_danger, out_danger, alternative='two-sided')
        d_cohen = diff / np.sqrt((np.std(in_danger)**2 + np.std(out_danger)**2) / 2)
        
        print(f"\nSSW {ssw_event}:")
        print(f"  In-window: mean={np.mean(in_danger):.3f} ± {np.std(in_danger):.3f} (n={len(in_danger)})")
        print(f"  Out-window: mean={np.mean(out_danger):.3f} ± {np.std(out_danger):.3f} (n={len(out_danger)})")
        print(f"  Δ = {diff:+.3f}, Cohen's d = {d_cohen:+.3f}, MW P = {p:.4f}")
        
        results['by_event'][ssw_event] = {
            'in_mean': float(np.mean(in_danger)),
            'out_mean': float(np.mean(out_danger)),
            'diff': float(diff),
            'cohen_d': float(d_cohen),
            'mw_p': float(p),
            'n_in': len(in_danger),
            'n_out': len(out_danger),
        }

# By parent region
print("\n" + "="*80)
print("BY REGION")
print("="*80)

for parent in sorted(set(d['parent_region'] for d in data)):
    region_data = [d for d in data if d['parent_region'] == parent]
    in_danger = [d['danger_max'] for d in region_data if d['in_ssw_window'] and d['danger_max'] is not None]
    out_danger = [d['danger_max'] for d in region_data if not d['in_ssw_window'] and d['danger_max'] is not None]
    
    if in_danger and out_danger:
        diff = np.mean(in_danger) - np.mean(out_danger)
        stat, p = mannwhitneyu(in_danger, out_danger, alternative='two-sided')
        d_cohen = diff / np.sqrt((np.std(in_danger)**2 + np.std(out_danger)**2) / 2)
        
        region_names = {'AT-07': 'Tyrol (Austria)', 'IT-32-BZ': 'South Tyrol (Italy)', 'IT-32-TN': 'Trentino (Italy)'}
        rname = region_names.get(parent, parent)
        
        print(f"\n{rname} ({parent}):")
        print(f"  SSW: {np.mean(in_danger):.3f} ± {np.std(in_danger):.3f} (n={len(in_danger)})")
        print(f"  Ctrl: {np.mean(out_danger):.3f} ± {np.std(out_danger):.3f} (n={len(out_danger)})")
        print(f"  Δ = {diff:+.3f}, d = {d_cohen:+.3f}, P = {p:.4f}")
        
        results['by_region'][parent] = {
            'name': rname,
            'in_mean': float(np.mean(in_danger)),
            'out_mean': float(np.mean(out_danger)),
            'diff': float(diff),
            'cohen_d': float(d_cohen),
            'mw_p': float(p),
        }

# By event x region
print("\n" + "="*80)
print("EVENT x REGION MATRIX")
print("="*80)

for ssw_event in sorted(set(d['ssw_event'] for d in data)):
    print(f"\nSSW {ssw_event}:")
    for parent in sorted(set(d['parent_region'] for d in data)):
        rdata = [d for d in data if d['ssw_event'] == ssw_event and d['parent_region'] == parent]
        in_d = [d['danger_max'] for d in rdata if d['in_ssw_window'] and d['danger_max'] is not None]
        out_d = [d['danger_max'] for d in rdata if not d['in_ssw_window'] and d['danger_max'] is not None]
        
        if in_d and out_d:
            diff = np.mean(in_d) - np.mean(out_d)
            region_names = {'AT-07': 'Tyrol', 'IT-32-BZ': 'S.Tyrol', 'IT-32-TN': 'Trentino'}
            print(f"  {region_names.get(parent, parent):10s}: SSW={np.mean(in_d):.2f}, Ctrl={np.mean(out_d):.2f}, Δ={diff:+.2f}")

# Phase-resolved
print("\n" + "="*80)
print("PHASE-RESOLVED (all regions combined)")
print("="*80)

phases = {
    'pre': (-15, -1),
    'during': (0, 7),
    'post': (8, 15),
    'late': (16, 30),
}

for phase_name, (d_start, d_end) in phases.items():
    phase_data = [d for d in data if d_start <= d['day_offset'] <= d_end]
    ctrl_data = [d for d in data if d['day_offset'] < -15 or d['day_offset'] > 15]
    
    in_d = [d['danger_max'] for d in phase_data if d['danger_max'] is not None]
    out_d = [d['danger_max'] for d in ctrl_data if d['danger_max'] is not None]
    
    if in_d and out_d:
        diff = np.mean(in_d) - np.mean(out_d)
        stat, p = mannwhitneyu(in_d, out_d, alternative='two-sided')
        d_cohen = diff / np.sqrt((np.std(in_d)**2 + np.std(out_d)**2) / 2)
        
        print(f"  {phase_name:6s} [{d_start:+3d} to {d_end:+3d}d]: "
              f"mean={np.mean(in_d):.3f}, ctrl={np.mean(out_d):.3f}, "
              f"Δ={diff:+.3f} (d={d_cohen:+.3f}) P={p:.4f}")
        
        results['phase_resolved'][phase_name] = {
            'mean': float(np.mean(in_d)),
            'ctrl': float(np.mean(out_d)),
            'diff': float(diff),
            'cohen_d': float(d_cohen),
            'mw_p': float(p),
        }

# ── Geographic gradient comparison ────────────────────────────────────────────
print("\n" + "="*80)
print("GEOGRAPHIC GRADIENT COMPARISON")
print("="*80)
print("EAWS prediction: Western Alps decrease, Central neutral, Northern increase")
print("ALBINA covers Central Alps → predict NEUTRAL to slight increase")
print()

# Overall ALBINA result
all_in = [d['danger_max'] for d in data if d['in_ssw_window'] and d['danger_max'] is not None]
all_out = [d['danger_max'] for d in data if not d['in_ssw_window'] and d['danger_max'] is not None]
overall_diff = np.mean(all_in) - np.mean(all_out)
_, overall_p = mannwhitneyu(all_in, all_out, alternative='two-sided')

print(f"Overall ALBINA: Δ = {overall_diff:+.3f}, P = {overall_p:.4f}")
print(f"  2023 SSW: Δ = {results['by_event']['2023-02-16']['diff']:+.3f}")
print(f"  2024 SSW: Δ = {results['by_event']['2024-03-04']['diff']:+.3f}")
print()
print("Comparison with EAWS (2012-2013 events):")
print("  France:       Δ = -0.40 (decrease) ← Western Alps")
print("  Switzerland:  Δ ≈  0.00 (neutral)  ← Central Alps")
print("  Austria:      Δ ≈  0.00 (neutral)  ← Central Alps")
print("  Germany:      Δ = +0.89 (increase) ← Northern Alps")
print(f"  ALBINA 2023:  Δ = {results['by_event']['2023-02-16']['diff']:+.3f}       ← Central Alps (validation)")
print(f"  ALBINA 2024:  Δ = {results['by_event']['2024-03-04']['diff']:+.3f}       ← Central Alps")

results['overall'] = {
    'in_mean': float(np.mean(all_in)),
    'out_mean': float(np.mean(all_out)),
    'diff': float(overall_diff),
    'mw_p': float(overall_p),
}

# Save
with open('data/results/r31c_albina_gradient.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print("\nResults saved to data/results/r31c_albina_gradient.json")
