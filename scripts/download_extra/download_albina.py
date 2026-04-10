"""
Download ALBINA avalanche bulletins for SSW event analysis.
Covers: Tyrol (AT-07), South Tyrol (IT-32-BZ), Trentino (IT-32-TN)
Also tries: Salzburg (AT-05), Vorarlberg (AT-08), Carinthia (AT-02), Styria (AT-06)
Date range: 2018-2025 (ALBINA launched ~2018)
"""
import requests
import json
import time
import os
from datetime import datetime, timedelta

# SSW events (Butler catalog + recent)
SSW_DATES = [
    '2018-02-12',  # Beast from the East
    '2019-01-02',
    '2021-01-05',
    '2023-02-16',
    '2024-03-04',  # latest SSW
]

# Regions to try
REGIONS = [
    'AT-07',    # Tyrol
    'IT-32-BZ', # South Tyrol  
    'IT-32-TN', # Trentino
    'AT-05',    # Salzburg
    'AT-08',    # Vorarlberg
    'AT-02',    # Carinthia
    'AT-06',    # Styria
]

DANGER_MAP = {
    'low': 1,
    'moderate': 2,
    'considerable': 3,
    'high': 4,
    'very_high': 5,
    'no_snow': 0,
    'no_rating': 0,
}

WINDOW = 15  # ±15 days

output_dir = r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\albina_bulletins'
os.makedirs(output_dir, exist_ok=True)

all_records = []
failed_dates = []

for ssw_date_str in SSW_DATES:
    ssw_date = datetime.strptime(ssw_date_str, '%Y-%m-%d')
    print(f"\n{'='*60}")
    print(f"SSW Event: {ssw_date_str}")
    print(f"{'='*60}")
    
    # Download ±30 days (SSW window + control buffer)
    for day_offset in range(-30, 31):
        date = ssw_date + timedelta(days=day_offset)
        date_str = date.strftime('%Y-%m-%d')
        
        for region in REGIONS:
            url = f'https://static.avalanche.report/bulletins/{date_str}/{region}.json'
            
            try:
                r = requests.get(url, timeout=10, headers={'User-Agent': 'ResearchBot/1.0'})
                
                if r.status_code == 200:
                    data = r.json()
                    
                    for bulletin in data:
                        regions = bulletin.get('regions', [])
                        forenoon = bulletin.get('forenoon', {})
                        
                        danger_above = forenoon.get('dangerRatingAbove', 'no_rating')
                        danger_below = forenoon.get('dangerRatingBelow', 'no_rating')
                        
                        danger_above_num = DANGER_MAP.get(danger_above, 0)
                        danger_below_num = DANGER_MAP.get(danger_below, 0)
                        danger_max = max(danger_above_num, danger_below_num)
                        
                        for reg in regions:
                            all_records.append({
                                'date': date_str,
                                'ssw_event': ssw_date_str,
                                'day_offset': day_offset,
                                'in_ssw_window': abs(day_offset) <= WINDOW,
                                'region': reg,
                                'parent_region': region,
                                'danger_above': danger_above,
                                'danger_below': danger_below,
                                'danger_above_num': danger_above_num,
                                'danger_below_num': danger_below_num,
                                'danger_max': danger_max,
                            })
                elif r.status_code == 404:
                    pass  # Normal for dates before region was added
                else:
                    failed_dates.append((date_str, region, r.status_code))
                    
            except Exception as e:
                failed_dates.append((date_str, region, str(e)))
            
            time.sleep(0.1)  # Be polite
    
    print(f"  Records so far: {len(all_records)}")

# Save raw data
with open(os.path.join(output_dir, 'albina_ssw_bulletins.json'), 'w') as f:
    json.dump(all_records, f, indent=2)

print(f"\n{'='*60}")
print(f"Total records: {len(all_records)}")
print(f"Failed requests: {len(failed_dates)}")
print(f"Saved to {output_dir}/albina_ssw_bulletins.json")

# Quick summary
import collections
by_region = collections.Counter(r['parent_region'] for r in all_records)
print(f"\nRecords by parent region:")
for reg, count in sorted(by_region.items()):
    print(f"  {reg}: {count}")

by_ssw = collections.Counter(r['ssw_event'] for r in all_records)
print(f"\nRecords by SSW event:")
for ssw, count in sorted(by_ssw.items()):
    print(f"  {ssw}: {count}")
