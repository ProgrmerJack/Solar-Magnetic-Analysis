"""Quick test of NVE API for key regions."""
import json, urllib.request, time

BASE = 'https://api01.nve.no/hydrology/forecast/avalanche/v6.2.1/api'
regions = [3010, 3022, 3028]  # Lyngen, Trollheimen, Jotunheimen

for rid in regions:
    url = f'{BASE}/AvalancheWarningByRegion/Simple/{rid}/1/2018-12-01/2019-02-28'
    try:
        req = urllib.request.Request(url)
        req.add_header('Accept', 'application/json')
        resp = urllib.request.urlopen(req, timeout=30)
        data = json.loads(resp.read().decode('utf-8'))
        if data:
            levels = [d['DangerLevel'] for d in data]
            print(f'Region {rid}: {len(data)} days')
            print(f'  Danger levels: {levels[:10]}...')
        else:
            print(f'Region {rid}: empty')
    except Exception as e:
        print(f'Region {rid}: ERROR {e}')
    time.sleep(1)
