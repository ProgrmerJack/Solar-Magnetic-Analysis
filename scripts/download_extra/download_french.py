"""Download French avalanche data from data.gouv.fr"""
import urllib.request, json, os

os.makedirs('data/cryosphere/french_avalanche', exist_ok=True)

# Search data.gouv.fr for avalanche datasets
url = 'https://www.data.gouv.fr/api/1/datasets/?q=avalanche+neige&page_size=20'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req, timeout=15) as resp:
    data = json.loads(resp.read())
    total = data.get('total', 0)
    print("Broader search: {} results".format(total))
    for item in data.get('data', [])[:10]:
        title = item.get('title', '')
        org = item.get('organization', {})
        org_name = org.get('name', 'N/A') if org else 'N/A'
        print("  - {} (org: {})".format(title, org_name))
        for res in item.get('resources', [])[:3]:
            rtitle = res.get('title', '')
            fmt = res.get('format', '')
            rurl = res.get('url', '')[:120]
            print("    -> {} [{}] {}".format(rtitle, fmt, rurl))

# Also search specifically for avalanche bulletins / BRA
print("\n--- Search for BRA ---")
url2 = 'https://www.data.gouv.fr/api/1/datasets/?q=bulletin+risque+avalanche&page_size=10'
req2 = urllib.request.Request(url2, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req2, timeout=15) as resp2:
    data2 = json.loads(resp2.read())
    for item in data2.get('data', [])[:5]:
        title = item.get('title', '')
        did = item.get('id', '')
        print("  {} ({})".format(title, did))
        for res in item.get('resources', [])[:3]:
            print("    {} [{}]".format(res.get('title', ''), res.get('format', '')))
            print("    URL: {}".format(res.get('url', '')[:150]))

# Try Meteofrance open data  
print("\n--- Meteofrance Open Data ---")
url3 = 'https://donneespubliques.meteofrance.fr/inspire/CLPA/'
try:
    req3 = urllib.request.Request(url3, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req3, timeout=10) as resp3:
        content = resp3.read().decode('utf-8', errors='replace')[:2000]
        print("CLPA response ({} chars): {}".format(len(content), content[:500]))
except Exception as e:
    print("CLPA error: {}".format(e))
