"""Search for French EPA avalanche data on data.gouv.fr"""
import urllib.request, json, os, sys

sys.stdout.reconfigure(encoding='utf-8')

# Search for avalanche datasets on French open data portal
urls = [
    'https://www.data.gouv.fr/api/1/datasets/?q=avalanche&page_size=20',
    'https://www.data.gouv.fr/api/1/datasets/?q=EPA+neige&page_size=10',
]

for url in urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=15)
        data = json.loads(resp.read())
        print(f"Query: {url.split('q=')[1].split('&')[0]}")
        print(f"Total results: {data.get('total', 0)}")
        for ds in data.get('data', []):
            title = ds.get('title', 'N/A')
            org = ds.get('organization', {})
            org_name = org.get('name', 'N/A') if org else 'N/A'
            desc = ds.get('description', '')[:200]
            resources = ds.get('resources', [])
            csv_res = [r for r in resources if r.get('format','').lower() in ('csv', 'json', 'xlsx')]
            print(f"\n  Dataset: {title}")
            print(f"  Org: {org_name}")
            print(f"  Desc: {desc}")
            print(f"  Data resources: {len(csv_res)}")
            for r in csv_res[:3]:
                rtitle = r.get('title', '?')
                rurl = r.get('url', '?')[:150]
                rfmt = r.get('format', '?')
                rsize = r.get('filesize', 0)
                print(f"    [{rfmt}] {rtitle} ({rsize} bytes)")
                print(f"    URL: {rurl}")
        print("\n" + "="*80)
    except Exception as e:
        print(f"Error: {e}")

# Also try CLPA (Carte de Localisation des Phénomènes d'Avalanche)
print("\nTrying CLPA/RTM avalanche data...")
try:
    url = 'https://www.data.gouv.fr/api/1/datasets/?q=CLPA+avalanche&page_size=5'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=15)
    data = json.loads(resp.read())
    for ds in data.get('data', []):
        print(f"  {ds.get('title', 'N/A')}")
        for r in ds.get('resources', [])[:3]:
            print(f"    [{r.get('format','?')}] {r.get('url','?')[:120]}")
except Exception as e:
    print(f"Error: {e}")

# Try Avalanche.org API for US data
print("\nTrying US avalanche.org API...")
try:
    url = 'https://api.avalanche.org/v2/public/products?type=accident&page=1&per_page=5'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=15)
    data = json.loads(resp.read())
    print(f"  Got {len(data) if isinstance(data, list) else 'dict'} results")
    if isinstance(data, list):
        for item in data[:3]:
            print(f"  {item.get('date', '?')} - {item.get('title', '?')[:80]}")
    elif isinstance(data, dict):
        for k, v in list(data.items())[:5]:
            print(f"  {k}: {str(v)[:100]}")
except Exception as e:
    print(f"Error: {e}")

# Try Canadian Avalanche Centre
print("\nTrying Canadian avalanche data...")
try:
    url = 'https://www.avalanche.ca/api/min/en/submissions?period=2023-01-01:2023-02-28&rows=5'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=15)
    data = json.loads(resp.read())
    print(f"  Got {type(data).__name__} with {len(data)} items")
    if isinstance(data, list) and len(data) > 0:
        print(f"  Keys: {list(data[0].keys())[:10] if isinstance(data[0], dict) else '?'}")
except Exception as e:
    print(f"Error: {e}")
