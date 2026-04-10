"""
Download ALBINA (Euregio: Austria/South Tyrol/Trentino) avalanche bulletin data
from public files at avalanche.report/albina_files/
Available from 2018-12-04 onwards.
"""
import urllib.request, json, os, time, sys
from datetime import datetime, timedelta

os.makedirs('data/cryosphere/albina', exist_ok=True)

# First check structure of one bulletin
url = 'https://avalanche.report/albina_files/2019-01-15/avalanche_report.json'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 Research'})
with urllib.request.urlopen(req, timeout=15) as resp:
    sample = json.load(resp)

print(f"Sample type: {type(sample)}")
if isinstance(sample, list):
    print(f"N items: {len(sample)}")
    if sample:
        print(f"Keys: {list(sample[0].keys())}")
        # Find danger level
        for item in sample[:3]:
            dl = item.get('dangerRatings', item.get('dangerLevel', item.get('maxDangerRating', '?')))
            region = item.get('regions', item.get('region', item.get('regionID', '?')))
            print(f"  Region: {region}, Danger: {dl}")
elif isinstance(sample, dict):
    print(f"Keys: {list(sample.keys())}")
    # Check nested structure
    for k, v in sample.items():
        if isinstance(v, list):
            print(f"  {k}: list of {len(v)}")
            if v and isinstance(v[0], dict):
                print(f"    First keys: {list(v[0].keys())}")
        elif isinstance(v, dict):
            print(f"  {k}: dict with {list(v.keys())[:5]}")
        else:
            print(f"  {k}: {type(v).__name__} = {str(v)[:100]}")

print("\n--- Full first item ---")
if isinstance(sample, list) and sample:
    print(json.dumps(sample[0], indent=2, default=str)[:3000])
elif isinstance(sample, dict):
    # Try to find bulletins key
    for k in ['bulletins', 'avalancheReport', 'reports']:
        if k in sample:
            items = sample[k]
            if isinstance(items, list) and items:
                print(f"Found in {k}:")
                print(json.dumps(items[0], indent=2, default=str)[:3000])
                break
    else:
        print(json.dumps(sample, indent=2, default=str)[:3000])
