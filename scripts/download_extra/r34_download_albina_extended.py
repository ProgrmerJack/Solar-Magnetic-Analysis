"""
Download ALBINA data for ALL SSW events within coverage (2018+).
Also downloads control winters for comparison.
Extends from 2 events to 5 events for Austrian/Italian replication.
"""
import requests
import json
import time
import os
from datetime import datetime, timedelta

SSW_DATES = [
    '2018-02-12',  # Beast from the East - MAJOR SSW
    '2019-01-02',  # January 2019 SSW
    '2021-01-05',  # January 2021 SSW  
    '2023-02-16',  # February 2023 SSW
    '2024-03-04',  # March 2024 SSW
]

REGIONS = ['AT-07', 'IT-32-BZ', 'IT-32-TN', 'AT-05', 'AT-08']

DANGER_MAP = {
    'low': 1, 'moderate': 2, 'considerable': 3,
    'high': 4, 'very_high': 5, 'no_snow': 0, 'no_rating': 0,
}

WINDOW = 30  # ±30 days for SSW + control buffer
output_dir = r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\albina_bulletins'
os.makedirs(output_dir, exist_ok=True)

all_records = []
failed = []
session = requests.Session()
session.headers['User-Agent'] = 'AcademicResearch/1.0 (SSW-avalanche study)'

for ssw_str in SSW_DATES:
    ssw_date = datetime.strptime(ssw_str, '%Y-%m-%d')
    print(f"\n{'='*60}")
    print(f"SSW Event: {ssw_str}")
    
    for offset in range(-WINDOW, WINDOW + 1):
        date = ssw_date + timedelta(days=offset)
        date_str = date.strftime('%Y-%m-%d')
        
        for region in REGIONS:
            url = f'https://static.avalanche.report/bulletins/{date_str}/{region}.json'
            try:
                r = session.get(url, timeout=10)
                if r.status_code == 200:
                    data = r.json()
                    for bulletin in data:
                        regions = bulletin.get('regions', [])
                        forenoon = bulletin.get('forenoon', {})
                        afternoon = bulletin.get('afternoon', {})
                        
                        danger_above = forenoon.get('dangerRatingAbove', 'no_rating')
                        danger_below = forenoon.get('dangerRatingBelow', 'no_rating')
                        danger_above_num = DANGER_MAP.get(danger_above, 0)
                        danger_below_num = DANGER_MAP.get(danger_below, 0)
                        danger_max = max(danger_above_num, danger_below_num)
                        
                        # Also get afternoon if available
                        pm_above = afternoon.get('dangerRatingAbove', danger_above) if afternoon else danger_above
                        pm_below = afternoon.get('dangerRatingBelow', danger_below) if afternoon else danger_below
                        
                        for reg in regions:
                            all_records.append({
                                'date': date_str,
                                'ssw_event': ssw_str,
                                'day_offset': offset,
                                'in_ssw_window': abs(offset) <= 15,
                                'region': reg,
                                'parent_region': region,
                                'danger_above': danger_above,
                                'danger_below': danger_below,
                                'danger_above_num': danger_above_num,
                                'danger_below_num': danger_below_num,
                                'danger_max': danger_max,
                            })
                elif r.status_code != 404:
                    failed.append((date_str, region, r.status_code))
            except Exception as e:
                failed.append((date_str, region, str(e)))
            
            time.sleep(0.05)
    
    print(f"  Records so far: {len(all_records)}")

# Save
outpath = os.path.join(output_dir, 'albina_all_ssw_extended.json')
with open(outpath, 'w') as f:
    json.dump(all_records, f, indent=2)

print(f"\n{'='*60}")
print(f"Total records: {len(all_records)}")
print(f"Failed: {len(failed)}")
print(f"Saved to {outpath}")

from collections import Counter
by_ssw = Counter(r['ssw_event'] for r in all_records)
print(f"\nBy SSW event:")
for s, n in sorted(by_ssw.items()):
    print(f"  {s}: {n}")

by_region = Counter(r['parent_region'] for r in all_records)
print(f"\nBy parent region:")
for reg, n in sorted(by_region.items()):
    print(f"  {reg}: {n}")

# Quick analysis: danger during SSW vs control
print(f"\n{'='*60}")
print("QUICK SSW vs CONTROL ANALYSIS")
ssw_records = [r for r in all_records if r['in_ssw_window'] and r['danger_max'] > 0]
ctrl_records = [r for r in all_records if not r['in_ssw_window'] and r['danger_max'] > 0]

import numpy as np
ssw_danger = np.array([r['danger_max'] for r in ssw_records])
ctrl_danger = np.array([r['danger_max'] for r in ctrl_records])
print(f"SSW mean danger: {ssw_danger.mean():.3f} (n={len(ssw_danger)})")
print(f"Ctrl mean danger: {ctrl_danger.mean():.3f} (n={len(ctrl_danger)})")
print(f"Difference: {ssw_danger.mean() - ctrl_danger.mean():.3f}")

# Per-event analysis
print("\nPer-event danger levels:")
for ssw_str in SSW_DATES:
    ssw_d = [r['danger_max'] for r in all_records if r['ssw_event'] == ssw_str and r['in_ssw_window'] and r['danger_max'] > 0]
    ctrl_d = [r['danger_max'] for r in all_records if r['ssw_event'] == ssw_str and not r['in_ssw_window'] and r['danger_max'] > 0]
    if ssw_d and ctrl_d:
        print(f"  {ssw_str}: SSW={np.mean(ssw_d):.2f} (n={len(ssw_d)}), Ctrl={np.mean(ctrl_d):.2f} (n={len(ctrl_d)}), diff={np.mean(ssw_d)-np.mean(ctrl_d):+.2f}")

# Per-region per-event
print("\nPer-region analysis (Tyrol AT-07 vs Trentino IT-32-TN):")
for ssw_str in SSW_DATES:
    for reg in ['AT-07', 'IT-32-TN']:
        ssw_d = [r['danger_max'] for r in all_records if r['ssw_event'] == ssw_str and r['in_ssw_window'] and r['parent_region'] == reg and r['danger_max'] > 0]
        ctrl_d = [r['danger_max'] for r in all_records if r['ssw_event'] == ssw_str and not r['in_ssw_window'] and r['parent_region'] == reg and r['danger_max'] > 0]
        if ssw_d and ctrl_d:
            print(f"  {ssw_str} {reg}: SSW={np.mean(ssw_d):.2f}, Ctrl={np.mean(ctrl_d):.2f}, diff={np.mean(ssw_d)-np.mean(ctrl_d):+.2f}")
