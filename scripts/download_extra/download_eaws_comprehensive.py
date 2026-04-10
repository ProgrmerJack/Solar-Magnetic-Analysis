"""
Comprehensive EAWS multi-country avalanche danger download using pyAvaCore.
Downloads Austrian, French, Italian, Swiss, German, Slovenian data for SSW events.
"""
import json, os, csv, sys, traceback
from datetime import datetime, timedelta

# Danger level text to number mapping
DANGER_MAP = {
    'low': 1, 'moderate': 2, 'considerable': 3, 'high': 4, 'very_high': 5,
    'no_snow': 0, 'no_rating': 0, '': 0, None: 0,
    '1': 1, '2': 2, '3': 3, '4': 4, '5': 5,
}

# Regions to download
REGIONS = [
    # Austria
    'AT-02', 'AT-03', 'AT-04', 'AT-05', 'AT-06', 'AT-07', 'AT-08',
    # Italy (Alto Adige/South Tyrol + Trentino)
    'IT-32-BZ', 'IT-32-TN',
    # Germany (Bavaria)
    'DE-BY',
    # Slovenia 
    'SI',
    # France
    'FR',
    # Spain (Aran)
    'ES-CT-L',
    # Norway
    'NO',
    # Czech Republic
    'CZ',
    # Poland
    'PL', 'PL-12',
    # Sweden
    'SE',
]

# SSW events in EAWS era (2018+)
SSW_EVENTS = ['2018-02-12', '2019-01-01', '2021-01-05', '2023-02-16']
CONTROL_DATES = ['2020-01-15', '2022-01-15', '2024-01-15']

os.makedirs('data/cryosphere/eaws_danger', exist_ok=True)

from avacore.pyAvaCore import get_bulletins

def download_day(region, date_str):
    """Download and parse danger ratings for a single region/date."""
    try:
        bulletins = get_bulletins(region_id=region, date=date_str, lang='en')
        records = []
        for b in bulletins.bulletins:
            # Get regions this bulletin covers
            region_names = []
            region_ids = []
            if hasattr(b, 'regions') and b.regions:
                for r in b.regions:
                    region_ids.append(r.regionID if hasattr(r, 'regionID') else '')
                    region_names.append(r.name if hasattr(r, 'name') else '')
            
            # Get danger ratings
            for dr in b.dangerRatings:
                val = dr.mainValue
                if isinstance(val, str):
                    danger_num = DANGER_MAP.get(val, 0)
                else:
                    danger_num = int(val) if val else 0
                
                elev = ''
                if dr.elevation:
                    if dr.elevation.lowerBound:
                        elev = 'above_' + str(dr.elevation.lowerBound)
                    elif dr.elevation.upperBound:
                        elev = 'below_' + str(dr.elevation.upperBound)
                
                records.append({
                    'date': date_str,
                    'provider_region': region,
                    'bulletin_regions': ';'.join(region_ids[:5]),
                    'danger_text': val,
                    'danger_num': danger_num,
                    'elevation': elev,
                    'time_period': dr.validTimePeriod or 'all_day',
                })
        return records
    except Exception as e:
        return []

# Download all data
all_records = []
n_success = 0
n_fail = 0

# For each SSW event and control, download ±20 day window
events = [(d, 'ssw') for d in SSW_EVENTS] + [(d, 'ctrl') for d in CONTROL_DATES]

for center_date_str, event_type in events:
    center_date = datetime.strptime(center_date_str, '%Y-%m-%d')
    
    for region in REGIONS:
        sys.stdout.write("{} {} {} ... ".format(event_type, center_date_str, region))
        sys.stdout.flush()
        
        region_records = []
        for offset in range(-20, 21):
            d = center_date + timedelta(days=offset)
            date_str = d.strftime('%Y-%m-%d')
            
            records = download_day(region, date_str)
            for r in records:
                r['event_type'] = event_type
                r['event_center'] = center_date_str
                r['offset_days'] = offset
            region_records.extend(records)
        
        if region_records:
            all_records.extend(region_records)
            n_success += 1
            sys.stdout.write("{} records\n".format(len(region_records)))
        else:
            n_fail += 1
            sys.stdout.write("no data\n")
        sys.stdout.flush()

print("\n=== Download Summary ===")
print("Total records: {}".format(len(all_records)))
print("Successful region-events: {}".format(n_success))
print("Failed/empty: {}".format(n_fail))

# Save to CSV
outfile = 'data/cryosphere/eaws_danger/eaws_danger_all.csv'
if all_records:
    fieldnames = ['date', 'provider_region', 'bulletin_regions', 'danger_text', 
                  'danger_num', 'elevation', 'time_period', 'event_type', 
                  'event_center', 'offset_days']
    with open(outfile, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_records)
    print("Saved to {}".format(outfile))
    
    # Summary by region
    import pandas as pd
    df = pd.DataFrame(all_records)
    print("\nRecords by provider region:")
    for region, count in df.groupby('provider_region').size().items():
        print("  {}: {} records".format(region, count))
    
    print("\nRecords by event:")
    for (etype, ecenter), count in df.groupby(['event_type', 'event_center']).size().items():
        print("  {} {}: {} records".format(etype, ecenter, count))
else:
    print("No records downloaded!")
