"""
Focused EAWS download: test which regions have data for recent dates,
then download systematically for post-2018 SSW events only.
"""
import json, os, csv, sys
from datetime import datetime, timedelta

os.makedirs('data/cryosphere/eaws_danger', exist_ok=True)

from avacore.pyAvaCore import get_bulletins

DANGER_MAP = {
    'low': 1, 'moderate': 2, 'considerable': 3, 'high': 4, 'very_high': 5,
    'no_snow': 0, 'no_rating': 0, '': 0, None: 0,
}

REGIONS = [
    'AT-02', 'AT-03', 'AT-04', 'AT-05', 'AT-06', 'AT-07', 'AT-08',
    'IT-32-BZ', 'IT-32-TN', 'DE-BY', 'SI', 'FR', 'ES-CT-L',
    'NO', 'CZ', 'PL', 'PL-12', 'SE',
]

# Step 1: Test which regions have data for a recent date
print("=== Testing region availability (2024-01-15) ===")
working_regions = []
for region in REGIONS:
    try:
        bulletins = get_bulletins(region_id=region, date='2024-01-15', lang='en')
        n = len(bulletins.bulletins)
        if n > 0:
            # Get max danger
            max_d = 0
            for b in bulletins.bulletins:
                for dr in b.dangerRatings:
                    val = dr.mainValue
                    d = DANGER_MAP.get(val, 0) if isinstance(val, str) else int(val or 0)
                    max_d = max(max_d, d)
            print("  {} -> {} bulletins, max danger {}".format(region, n, max_d))
            working_regions.append(region)
        else:
            print("  {} -> 0 bulletins".format(region))
    except Exception as e:
        print("  {} -> ERROR {}".format(region, str(e)[:60]))

print("\nWorking regions: {}".format(working_regions))

# Step 2: Test earliest available date for AT-07
print("\n=== Testing earliest data for AT-07 ===")
for year in [2017, 2018, 2019, 2020]:
    test_date = '{}-01-15'.format(year)
    try:
        b = get_bulletins(region_id='AT-07', date=test_date, lang='en')
        print("  {} -> {} bulletins".format(test_date, len(b.bulletins)))
    except Exception as e:
        print("  {} -> ERROR".format(test_date))

# Step 3: Download data for SSW events that fall within data availability
# SSW events: 2019-01-01, 2021-01-05, 2023-02-16 (skip 2018-02-12)
SSW_EVENTS = ['2019-01-01', '2021-01-05', '2023-02-16']
CONTROLS = ['2020-01-15', '2022-01-15', '2024-01-15']

def download_day(region, date_str):
    try:
        bulletins = get_bulletins(region_id=region, date=date_str, lang='en')
        records = []
        for b in bulletins.bulletins:
            region_ids = [r.regionID for r in (b.regions or []) if hasattr(r, 'regionID')]
            
            for dr in b.dangerRatings:
                val = dr.mainValue
                danger_num = DANGER_MAP.get(val, 0) if isinstance(val, str) else int(val or 0)
                elev = ''
                if dr.elevation:
                    if dr.elevation.lowerBound:
                        elev = 'hi'
                    elif dr.elevation.upperBound:
                        elev = 'lo'
                records.append({
                    'date': date_str,
                    'provider_region': region,
                    'n_subregions': len(region_ids),
                    'danger_text': val,
                    'danger_num': danger_num,
                    'elevation': elev,
                })
        return records
    except:
        return []

print("\n=== Downloading SSW + Control Windows ===")
all_records = []
events = [(d, 'ssw') for d in SSW_EVENTS] + [(d, 'ctrl') for d in CONTROLS]

for center_str, event_type in events:
    center = datetime.strptime(center_str, '%Y-%m-%d')
    
    for region in working_regions:
        sys.stdout.write("{} {} {} ".format(event_type, center_str, region))
        sys.stdout.flush()
        
        region_count = 0
        for offset in range(-20, 21):
            d = center + timedelta(days=offset)
            date_str = d.strftime('%Y-%m-%d')
            records = download_day(region, date_str)
            for r in records:
                r['event_type'] = event_type
                r['event_center'] = center_str
                r['offset_days'] = offset
            all_records.extend(records)
            region_count += len(records)
        
        sys.stdout.write("-> {} records\n".format(region_count))
        sys.stdout.flush()

# Save
outfile = 'data/cryosphere/eaws_danger/eaws_danger_all.csv'
print("\n=== Saving {} records ===".format(len(all_records)))

if all_records:
    import pandas as pd
    df = pd.DataFrame(all_records)
    df.to_csv(outfile, index=False)
    print("Saved to {}".format(outfile))
    
    print("\nBy region:")
    print(df.groupby('provider_region').size().to_string())
    print("\nBy event:")
    print(df.groupby(['event_type', 'event_center']).size().to_string())
