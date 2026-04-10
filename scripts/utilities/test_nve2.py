"""Quick check which seasons have data for Jotunheimen (3028)."""
import json, urllib.request

BASE = 'https://api01.nve.no/hydrology/forecast/avalanche/v6.2.1/api'

for y1 in range(2017, 2025):
    y2 = y1 + 1
    url = f'{BASE}/AvalancheWarningByRegion/Simple/3028/1/{y1}-12-01/{y2}-02-28'
    try:
        req = urllib.request.Request(url)
        req.add_header('Accept', 'application/json')
        resp = urllib.request.urlopen(req, timeout=30)
        data = json.loads(resp.read().decode('utf-8'))
        print(f'{y1}-{y2}: {len(data)} records')
    except Exception as e:
        print(f'{y1}-{y2}: ERROR {e}')
