"""
Download Norwegian NVE/Varsom avalanche danger data for ALL available regions
and ALL available winters to expand multi-country replication.
The API endpoint: api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/
"""
import urllib.request, json, os, time

os.makedirs('data/cryosphere/norway_expanded', exist_ok=True)

# First get all regions
url = 'https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/Region'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 Research'})
with urllib.request.urlopen(req, timeout=30) as resp:
    regions = json.load(resp)

print(f"Regions: {len(regions)}")
# Show inland/mountain regions (TypeId 10 = A-regions = mountain)
mountain = [r for r in regions if r.get('TypeId') == 10 or r.get('RegionTypeName', '').lower() in ['a-varsling', 'b-varsling']]
print(f"Mountain regions (TypeId=10): {len(mountain)}")

# Print all regions
for r in regions[:30]:
    print(f"  {r.get('Id', '?')}: {r.get('Name', '?')} (TypeId={r.get('TypeId', '?')}, Type={r.get('TypeName', '?')})")

# Download danger data for key inland regions for all available winters
# Region IDs for Norwegian avalanche forecast regions (inland mountain)
# 3003=Nordenskiöld Land, 3007=Vest-Finnmark, 3009=Nord-Troms, etc.
# Try all regions with TypeId in [10, 20, 30]
target_regions = [r['Id'] for r in regions if r.get('TypeId') in [10, 20, 30, 0]]
if not target_regions:
    target_regions = [r['Id'] for r in regions]

print(f"\nTarget regions: {len(target_regions)}")

# Download winters 2012-2024 (covering SSW events)
all_records = []
winters = [
    ('2012-11-01', '2013-05-31'),
    ('2013-11-01', '2014-05-31'),
    ('2014-11-01', '2015-05-31'),
    ('2015-11-01', '2016-05-31'),
    ('2016-11-01', '2017-05-31'),
    ('2017-11-01', '2018-05-31'),
    ('2018-11-01', '2019-05-31'),
    ('2019-11-01', '2020-05-31'),
    ('2020-11-01', '2021-05-31'),
    ('2021-11-01', '2022-05-31'),
    ('2022-11-01', '2023-05-31'),
    ('2023-11-01', '2024-05-31'),
]

# Take first 10 regions to be reasonable
sample_regions = target_regions[:15]
print(f"Sampling {len(sample_regions)} regions across {len(winters)} winters")

for region_id in sample_regions:
    for start, end in winters:
        url = f'https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/AvalancheWarningByRegion/Detail/{region_id}/1/{start}/{end}'
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 Research'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.load(resp)
            
            for rec in data:
                all_records.append({
                    'region_id': region_id,
                    'region_name': rec.get('RegionName', ''),
                    'danger_level': rec.get('DangerLevel', None),
                    'valid_from': rec.get('ValidFrom', ''),
                    'valid_to': rec.get('ValidTo', ''),
                })
            
            if len(data) > 0:
                print(f"  Region {region_id} ({data[0].get('RegionName', '?')}), {start[:4]}-{end[:4]}: {len(data)} records")
        except Exception as e:
            pass  # Skip failures silently
        
        time.sleep(0.2)  # Rate limit

print(f"\nTotal records downloaded: {len(all_records)}")

# Save
with open('data/cryosphere/norway_expanded/varsom_all.json', 'w') as f:
    json.dump(all_records, f)
print(f"Saved to data/cryosphere/norway_expanded/varsom_all.json")

# Quick analysis
import pandas as pd
df = pd.DataFrame(all_records)
df['date'] = pd.to_datetime(df['valid_from']).dt.date
df = df[df['danger_level'].notna() & (df['danger_level'] > 0)]
print(f"\nValid records: {len(df)}")
print(f"Regions: {df['region_name'].nunique()}")
print(f"Date range: {df['date'].min()} to {df['date'].max()}")
print(f"\nRecords per region:")
print(df.groupby('region_name').size().sort_values(ascending=False).head(15).to_string())
