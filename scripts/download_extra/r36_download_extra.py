"""
Try to download additional avalanche databases for multi-country replication:
1. Canadian Avalanche Centre (Avalanche Canada)
2. Japanese NIED avalanche data
3. Scottish Avalanche Information Service (SAIS)
4. New Zealand avalanche data
5. Icelandic Meteorological Office avalanche data
"""
import urllib.request, json, os

os.makedirs('data/cryosphere/canada', exist_ok=True)
os.makedirs('data/cryosphere/japan', exist_ok=True)
os.makedirs('data/cryosphere/scotland', exist_ok=True)

results = {}

# 1. Canadian Avalanche Centre - MIN database
urls_to_try = {
    'cac_min_api': 'https://www.avalanche.ca/api/min/en/submissions?page=1&per_page=100',
    'cac_incidents': 'https://www.avalanche.ca/api/incidents',
    'cac_forecasts_api': 'https://www.avalanche.ca/api/forecasts',
    'cac_bulletin': 'https://www.avalanche.ca/api/bulletin-archive/latest.json',
}

for name, url in urls_to_try.items():
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Research)', 'Accept': 'application/json'})
        resp = urllib.request.urlopen(req, timeout=15)
        data = resp.read()
        with open(f'data/cryosphere/canada/{name}.json', 'wb') as f:
            f.write(data)
        results[name] = f'OK ({len(data)} bytes)'
        print(f'{name}: OK ({len(data)} bytes)')
    except Exception as e:
        results[name] = str(e)
        print(f'{name}: FAIL - {e}')

# 2. Japanese NIED Snow and Ice Research Center
urls_jp = {
    'nied_snow': 'https://www.bosai.go.jp/snow/snow/index.html',
    'nied_avalanche_db': 'https://www.bosai.go.jp/snow/avalanche/avalanche_data.html',
}
for name, url in urls_jp.items():
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=10)
        data = resp.read()
        results[name] = f'Page loaded ({len(data)} bytes)'
        print(f'{name}: Page loaded ({len(data)} bytes)')
    except Exception as e:
        results[name] = str(e)
        print(f'{name}: FAIL - {e}')

# 3. Scottish Avalanche Information Service (SAIS)
urls_sais = {
    'sais_incidents': 'https://www.sais.gov.uk/api/incidents',
    'sais_blog': 'https://www.sais.gov.uk/api/blog-posts?page=1',
}
for name, url in urls_sais.items():
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
        resp = urllib.request.urlopen(req, timeout=10)
        data = resp.read()
        with open(f'data/cryosphere/scotland/{name}.json', 'wb') as f:
            f.write(data)
        results[name] = f'OK ({len(data)} bytes)'
        print(f'{name}: OK ({len(data)} bytes)')
    except Exception as e:
        results[name] = str(e)
        print(f'{name}: FAIL - {e}')

# 4. Lawinen.report (ALBINA) - try all available dates
urls_albina = {}
for year in [2020, 2021, 2022, 2023, 2024, 2025]:
    for month in [1, 2, 3, 12]:
        if year == 2025 and month == 12:
            continue
        if year == 2020 and month == 12:
            continue
        date = f'{year}-{month:02d}-15'
        urls_albina[f'albina_{date}'] = f'https://avalanche.report/albina_files/latest/{date}_en.json'

for name, url in list(urls_albina.items())[:5]:  # Just test a few
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=10)
        data = resp.read()
        results[name] = f'OK ({len(data)} bytes)'
        print(f'{name}: OK ({len(data)} bytes)')
    except Exception as e:
        results[name] = str(e)
        print(f'{name}: FAIL - {e}')

# 5. Try EAWS data API
urls_eaws = {
    'eaws_regions': 'https://avalanches.org/wp-admin/admin-ajax.php?action=get_regions',
    'eaws_api_v2': 'https://api.avalanches.org/v2/regions',
}
for name, url in urls_eaws.items():
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=10)
        data = resp.read()
        results[name] = f'OK ({len(data)} bytes)'
        print(f'{name}: OK ({len(data)} bytes)')
    except Exception as e:
        results[name] = str(e)
        print(f'{name}: FAIL - {e}')

# 6. Try Avalanche Canada detailed forecast archive  
urls_cac2 = {
    'cac_forecasts_regions': 'https://www.avalanche.ca/api/forecasts/en',
    'cac_archive': 'https://www.avalanche.ca/api/bulletin-archive/en/2024-01-15',
}
for name, url in urls_cac2.items():
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
        resp = urllib.request.urlopen(req, timeout=15)
        data = resp.read()
        with open(f'data/cryosphere/canada/{name}.json', 'wb') as f:
            f.write(data)
        results[name] = f'OK ({len(data)} bytes)'
        print(f'{name}: OK ({len(data)} bytes)')
    except Exception as e:
        results[name] = str(e)
        print(f'{name}: FAIL - {e}')

print('\n=== Summary ===')
successes = sum(1 for v in results.values() if v.startswith('OK'))
print(f'Successes: {successes}/{len(results)}')
for k, v in results.items():
    status = 'OK' if v.startswith('OK') or 'loaded' in v else 'FAIL'
    print(f'  {k:30s}: {status}')
