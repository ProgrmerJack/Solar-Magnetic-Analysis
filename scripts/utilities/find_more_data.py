"""Search for additional avalanche data sources - ALBINA, Japanese, Scottish, Icelandic"""
import urllib.request, json, sys, os
from datetime import datetime, timedelta

sys.stdout.reconfigure(encoding='utf-8')

print("="*80)
print("1. ALBINA (Austria/South Tyrol/Trentino) - Archived bulletins")
print("="*80)

# ALBINA has an API for archived bulletins
# Format: https://api.avalanche.report/albina/api/bulletins?date=YYYY-MM-DD&lang=en
# SSW events in ALBINA range (2017+): 2018-02-12, 2019-01-02, 2021-01-05, 2023-02-16
ssw_dates = ['2018-02-12', '2019-01-02', '2021-01-05', '2023-02-16']

for ssw_date in ssw_dates:
    try:
        url = f'https://api.avalanche.report/albina/api/bulletins?date={ssw_date}&lang=en'
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=15)
        data = json.loads(resp.read())
        
        if isinstance(data, dict):
            bulletins = data.get('bulletins', data.get('data', []))
            if not bulletins and 'dangerRatings' in str(data):
                bulletins = [data]
        elif isinstance(data, list):
            bulletins = data
        else:
            bulletins = []
        
        print(f"\n  SSW {ssw_date}: Got {len(bulletins)} bulletins")
        if bulletins:
            b = bulletins[0]
            if isinstance(b, dict):
                keys = list(b.keys())[:10]
                print(f"    Keys: {keys}")
                # Look for danger ratings
                if 'dangerRatings' in b:
                    print(f"    Danger ratings: {b['dangerRatings']}")
                if 'validTime' in b:
                    print(f"    Valid time: {b['validTime']}")
    except Exception as e:
        print(f"  SSW {ssw_date}: Error - {e}")

print("\n" + "="*80)
print("2. ALBINA - Try different API format")  
print("="*80)

# Try older API format
try:
    url = 'https://api.avalanche.report/albina/api/bulletins/2023-02-16'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=15)
    data = json.loads(resp.read())
    print(f"  Type: {type(data).__name__}")
    if isinstance(data, list) and len(data) > 0:
        b = data[0]
        print(f"  Keys: {list(b.keys())[:15]}")
        if 'dangerRatings' in b:
            for dr in b['dangerRatings'][:3]:
                print(f"    Rating: {dr}")
        if 'regions' in b:
            print(f"  Regions: {len(b['regions'])} regions")
            for r in b['regions'][:3]:
                print(f"    {r}")
    elif isinstance(data, dict):
        print(f"  Keys: {list(data.keys())[:15]}")
except Exception as e:
    print(f"  Error: {e}")

print("\n" + "="*80)
print("3. Japanese NIED Snow Disaster data")
print("="*80)

try:
    url = 'https://www.bosai.go.jp/seppyo/avalanche_db/en/'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=15)
    html = resp.read().decode('utf-8', errors='replace')
    print(f"  Got {len(html)} chars of HTML")
    # Check for data links
    import re
    links = re.findall(r'href="([^"]*(?:csv|json|download)[^"]*)"', html, re.I)
    print(f"  Data links found: {len(links)}")
    for l in links[:5]:
        print(f"    {l}")
except Exception as e:
    print(f"  Error: {e}")

print("\n" + "="*80)
print("4. Scottish SAIS avalanche data") 
print("="*80)

try:
    url = 'https://www.sais.gov.uk/api'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=10)
    data = json.loads(resp.read())
    print(f"  Got: {type(data).__name__}")
except Exception as e:
    print(f"  Error: {e}")

try:
    url = 'https://www.sais.gov.uk/'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=10)
    html = resp.read().decode('utf-8', errors='replace')
    import re
    api_links = re.findall(r'href="([^"]*(?:api|data|bulletin)[^"]*)"', html, re.I)
    print(f"  Links: {api_links[:5]}")
except Exception as e:
    print(f"  Error: {e}")

print("\n" + "="*80)
print("5. Icelandic Met Office avalanche data")
print("="*80)

try:
    url = 'https://en.vedur.is/avalanches/forecast/'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=10)
    html = resp.read().decode('utf-8', errors='replace')
    print(f"  Got {len(html)} chars")
except Exception as e:
    print(f"  Error: {e}")

print("\n" + "="*80)
print("6. Check existing ALBINA data in repo")
print("="*80)

albina_dir = 'data/cryosphere/albina_bulletins'
if os.path.exists(albina_dir):
    files = os.listdir(albina_dir)
    print(f"  Found {len(files)} files in {albina_dir}")
    for f in sorted(files)[:5]:
        fpath = os.path.join(albina_dir, f)
        size = os.path.getsize(fpath)
        print(f"    {f} ({size:,} bytes)")
    if len(files) > 5:
        print(f"    ... and {len(files)-5} more")
    for f in sorted(files)[-5:]:
        fpath = os.path.join(albina_dir, f)
        size = os.path.getsize(fpath)
        print(f"    {f} ({size:,} bytes)")
else:
    print(f"  No ALBINA directory found")

# Check for existing EAWS data
eaws_dir = 'data/cryosphere/eaws'
if os.path.exists(eaws_dir):
    files = os.listdir(eaws_dir)
    print(f"\n  Found {len(files)} files in {eaws_dir}")
    for f in sorted(files)[:10]:
        fpath = os.path.join(eaws_dir, f)
        size = os.path.getsize(fpath)
        print(f"    {f} ({size:,} bytes)")
else:
    print(f"\n  No EAWS directory found")
