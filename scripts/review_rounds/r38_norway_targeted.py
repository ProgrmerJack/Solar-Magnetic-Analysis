"""
Targeted Norwegian Varsom download - only key mountain regions and SSW-relevant winters.
Use shorter timeout, fewer regions, focused on SSW event dates.
"""
import urllib.request, json, os, time, sys

os.makedirs('data/cryosphere/norway_expanded', exist_ok=True)

# SSW events that fall within Varsom data availability (2013+)
# Major SSWs: 2013-01-06, 2018-02-12, 2019-01-02, 2021-01-05, 2023-02-16
ssw_winters = {
    '2012-2013': ('2012-11-01', '2013-04-30'),
    '2017-2018': ('2017-11-01', '2018-04-30'),
    '2018-2019': ('2018-11-01', '2019-04-30'),
    '2020-2021': ('2020-11-01', '2021-04-30'),
    '2022-2023': ('2022-11-01', '2023-04-30'),
}

# Non-SSW control winters
control_winters = {
    '2013-2014': ('2013-11-01', '2014-04-30'),
    '2014-2015': ('2014-11-01', '2015-04-30'),
    '2015-2016': ('2015-11-01', '2016-04-30'),
    '2016-2017': ('2016-11-01', '2017-04-30'),
    '2019-2020': ('2019-11-01', '2020-04-30'),
    '2021-2022': ('2021-11-01', '2022-04-30'),
    '2023-2024': ('2023-11-01', '2024-04-30'),
}

# Key inland mountain regions (A-varsling = TypeId 10)
# Focus on 6 key inland regions covering different latitudes
key_regions = [3003, 3009, 3013, 3022, 3028, 3029]  # Spread from Svalbard to southern Norway

all_records = []
all_winters = {**ssw_winters, **control_winters}

for region_id in key_regions:
    for winter_name, (start, end) in sorted(all_winters.items()):
        url = f'https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/AvalancheWarningByRegion/Detail/{region_id}/1/{start}/{end}'
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 Research'})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.load(resp)
            
            for rec in data:
                all_records.append({
                    'region_id': region_id,
                    'region_name': rec.get('RegionName', ''),
                    'danger_level': rec.get('DangerLevel', None),
                    'valid_from': rec.get('ValidFrom', ''),
                    'valid_to': rec.get('ValidTo', ''),
                    'winter': winter_name,
                    'is_ssw_winter': winter_name in ssw_winters,
                })
            
            n = len(data)
            sys.stdout.write(f"  {region_id} {winter_name}: {n} recs\n")
            sys.stdout.flush()
        except Exception as e:
            sys.stdout.write(f"  {region_id} {winter_name}: FAIL ({str(e)[:50]})\n")
            sys.stdout.flush()
        
        time.sleep(0.3)

print(f"\nTotal records: {len(all_records)}")

with open('data/cryosphere/norway_expanded/varsom_targeted.json', 'w') as f:
    json.dump(all_records, f)

# Quick summary
import pandas as pd
df = pd.DataFrame(all_records)
df['date'] = pd.to_datetime(df['valid_from']).dt.date
df = df[df['danger_level'].notna() & (df['danger_level'] > 0)]
print(f"Valid records: {len(df)}")
print(f"Regions: {df['region_name'].nunique()}")
print(f"Date range: {df['date'].min()} to {df['date'].max()}")
print(f"SSW winters: {df[df['is_ssw_winter']].shape[0]}")
print(f"Control winters: {df[~df['is_ssw_winter']].shape[0]}")
print(f"\nDanger level by SSW status:")
print(df.groupby('is_ssw_winter')['danger_level'].describe())
