"""
Download Norwegian Varsom Simple API data for SSW analysis.
Simple endpoint is faster - just need DangerLevel by date/region.
Focus on SSW windows and matched control periods.
"""
import urllib.request, json, os, time, sys
import pandas as pd

os.makedirs('data/cryosphere/norway_expanded', exist_ok=True)

# Key inland mountain regions covering latitudinal gradient
regions = [3009, 3010, 3011, 3013, 3022, 3023, 3024, 3027, 3028, 3029]

# SSW events (date, winter)
ssw_events = [
    ('2013-01-06', '2012-2013'),
    ('2018-02-12', '2017-2018'),
    ('2019-01-02', '2018-2019'),
    ('2021-01-05', '2020-2021'),
    ('2023-02-16', '2022-2023'),
]

# Download 3-month windows per winter (Nov-Apr) using quarterly chunks
# The API seems faster for shorter date ranges
all_records = []
success = 0
fail = 0

for region_id in regions:
    for ssw_date, winter in ssw_events:
        # Download 3-month period centered on SSW (but within the same winter)
        ssw_year = int(ssw_date[:4])
        ssw_month = int(ssw_date[5:7])
        
        # Get Dec-Mar window
        start = f'{ssw_year-1 if ssw_month <= 3 else ssw_year}-12-01'
        end = f'{ssw_year if ssw_month <= 3 else ssw_year+1}-03-31'
        
        url = f'https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/AvalancheWarningByRegion/Simple/{region_id}/1/{start}/{end}'
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 Research'})
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.load(resp)
            
            for rec in data:
                all_records.append({
                    'region_id': region_id,
                    'region_name': rec.get('RegionName', ''),
                    'danger_level': int(rec.get('DangerLevel', 0)),
                    'valid_from': rec.get('ValidFrom', ''),
                    'ssw_date': ssw_date,
                    'winter': winter,
                })
            success += 1
            sys.stdout.write(f"  {region_id} {winter}: {len(data)} recs\n")
            sys.stdout.flush()
        except Exception as e:
            fail += 1
            sys.stdout.write(f"  {region_id} {winter}: FAIL\n")
            sys.stdout.flush()
        time.sleep(0.5)
    
    # Also get control winters (no SSW)
    control_periods = [
        ('2014-12-01', '2015-03-31', '2014-2015'),
        ('2015-12-01', '2016-03-31', '2015-2016'),
        ('2016-12-01', '2017-03-31', '2016-2017'),
        ('2019-12-01', '2020-03-31', '2019-2020'),
        ('2021-12-01', '2022-03-31', '2021-2022'),
    ]
    
    for start, end, winter in control_periods:
        url = f'https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/AvalancheWarningByRegion/Simple/{region_id}/1/{start}/{end}'
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 Research'})
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.load(resp)
            
            for rec in data:
                all_records.append({
                    'region_id': region_id,
                    'region_name': rec.get('RegionName', ''),
                    'danger_level': int(rec.get('DangerLevel', 0)),
                    'valid_from': rec.get('ValidFrom', ''),
                    'ssw_date': None,
                    'winter': winter,
                })
            success += 1
            sys.stdout.write(f"  {region_id} {winter} (ctrl): {len(data)} recs\n")
            sys.stdout.flush()
        except Exception as e:
            fail += 1
            sys.stdout.write(f"  {region_id} {winter} (ctrl): FAIL\n")
            sys.stdout.flush()
        time.sleep(0.5)

print(f"\nDone: {success} success, {fail} fail, {len(all_records)} total records")

with open('data/cryosphere/norway_expanded/varsom_targeted.json', 'w') as f:
    json.dump(all_records, f)

# Analyze
if all_records:
    df = pd.DataFrame(all_records)
    df['date'] = pd.to_datetime(df['valid_from']).dt.date
    df = df[df['danger_level'] > 0]
    
    ssw_df = df[df['ssw_date'].notna()]
    ctrl_df = df[df['ssw_date'].isna()]
    
    print(f"\nSSW periods: {len(ssw_df)} region-days, mean danger={ssw_df['danger_level'].mean():.3f}")
    print(f"Control:     {len(ctrl_df)} region-days, mean danger={ctrl_df['danger_level'].mean():.3f}")
    
    # For SSW periods, compute ±15d window around each SSW date
    from datetime import datetime, timedelta
    ssw_window = []
    matched_ctrl = []
    
    for _, row in ssw_df.iterrows():
        if row['ssw_date']:
            ssw_dt = datetime.strptime(row['ssw_date'], '%Y-%m-%d').date()
            rec_dt = row['date']
            delta = (rec_dt - ssw_dt).days
            if -15 <= delta <= 15:
                ssw_window.append(row)
    
    ssw_win_df = pd.DataFrame(ssw_window)
    if len(ssw_win_df) > 0:
        print(f"\nSSW ±15d window: {len(ssw_win_df)} region-days, mean danger={ssw_win_df['danger_level'].mean():.3f}")
        print(f"Full-winter control: mean danger={ctrl_df['danger_level'].mean():.3f}")
        
        # Per-event analysis
        for ssw_date in ssw_df['ssw_date'].unique():
            if ssw_date:
                evt = ssw_win_df[ssw_win_df['ssw_date'] == ssw_date]
                if len(evt) > 0:
                    print(f"  SSW {ssw_date}: n={len(evt)}, mean danger={evt['danger_level'].mean():.3f}")
