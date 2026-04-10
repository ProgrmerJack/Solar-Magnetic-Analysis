import urllib.request
import json

# Try CAIC directly and other data sources
urls = [
    ('CAIC data', 'https://avalanche.state.co.us/caic/acc/acc_repo.php?date1=01/01/2019&date2=12/31/2019'),
    ('CAIC api', 'https://api.cotrip.org/v1/avalanche'),
    ('Avalanche Canada', 'https://www.avalanche.ca/api/acs/products?date=2024-01-05'),
    ('Avalanche Canada 2', 'https://avcan-api.avalanche.ca/products?date=2024-01-05'),
]

for name, url in urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
        resp = urllib.request.urlopen(req, timeout=10)
        ct = resp.headers.get('Content-Type', '')
        data = resp.read(1000).decode('utf-8', errors='replace')
        print(f'{name}: Status {resp.status}, CT: {ct[:40]}')
        if 'json' in ct:
            print(f'  {data[:300]}')
        else:
            print(f'  HTML response (not API)')
    except Exception as e:
        print(f'{name}: {str(e)[:150]}')
    print()

# Try NSIDC for avalanche data
url = 'https://nsidc.org/api/search?q=avalanche&format=json'
try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=10)
    data = json.loads(resp.read())
    print(f'NSIDC: {type(data)}')
    if isinstance(data, dict):
        print(f'Keys: {list(data.keys())[:5]}')
except Exception as e:
    print(f'NSIDC: {e}')
