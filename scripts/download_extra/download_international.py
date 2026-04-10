"""Download additional international avalanche data sources"""
import urllib.request, json, os, csv

# ============================================================
# 1. AVALANCHE CANADA - try all known endpoints
# ============================================================
print("=== Avalanche Canada ===")
ca_urls = [
    'https://www.avalanche.ca/api/forecasts',
    'https://www.avalanche.ca/api/bulletin-archive',
    'https://api.avalanche.ca/forecasts',
    'https://avalanche.ca/api/forecasts/archive',
    'https://www.avalanche.ca/api/min/observations',
]
for url in ca_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
            print("  {} -> {} bytes".format(url.split('.ca')[-1], len(data)))
            try:
                j = json.loads(data)
                if isinstance(j, list):
                    print("    List of {} items".format(len(j)))
                    if j:
                        print("    Keys: {}".format(list(j[0].keys()) if isinstance(j[0], dict) else type(j[0])))
                elif isinstance(j, dict):
                    print("    Dict keys: {}".format(list(j.keys())[:10]))
            except:
                print("    Not JSON, first 200 chars: {}".format(data.decode('utf-8', errors='replace')[:200]))
    except Exception as e:
        print("  {} -> {}".format(url.split('.ca')[-1], str(e)[:80]))

# ============================================================
# 2. AUSTRIAN LAWINEN.REPORT (EAWS/ALBINA) 
# ============================================================
print("\n=== Austrian ALBINA ===")
# Try the ALBINA API that serves lawinen.report
au_urls = [
    'https://api.avalanche.report/albina/api/bulletins?date=2024-01-15&lang=en',
    'https://avalanche.report/albina/api/bulletins?date=2024-01-15&lang=en', 
    'https://api.avalanche.report/v1/bulletins?date=2024-01-15',
    'https://static.avalanche.report/bulletins/2024/en/2024-01-15_en.json',
    'https://static.avalanche.report/bulletins/2024/2024-01-15.json',
]
for url in au_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
            print("  {} -> {} bytes".format(url.split('report')[-1], len(data)))
            try:
                j = json.loads(data)
                if isinstance(j, list):
                    print("    List of {} items".format(len(j)))
                elif isinstance(j, dict):
                    print("    Keys: {}".format(list(j.keys())[:10]))
            except:
                print("    Not JSON")
    except Exception as e:
        print("  {} -> {}".format(url.split('report')[-1], str(e)[:80]))

# ============================================================
# 3. JAPANESE SNOW DISASTER DATABASE
# ============================================================
print("\n=== Japan NIED ===")
jp_urls = [
    'https://www.bosai.go.jp/seppyo/',
    'https://snow.bosai.go.jp/',
]
for url in jp_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read().decode('utf-8', errors='replace')
            print("  {} -> {} bytes".format(url, len(data)))
            # Look for API or data links
            for line in data.split('\n'):
                lower = line.lower()
                if 'api' in lower or 'download' in lower or 'data' in lower or 'csv' in lower:
                    stripped = line.strip()[:120]
                    if stripped:
                        print("    Link: {}".format(stripped))
    except Exception as e:
        print("  {} -> {}".format(url, str(e)[:80]))

# ============================================================
# 4. SCOTTISH AVALANCHE (SAIS)
# ============================================================
print("\n=== Scottish SAIS ===")
sais_urls = [
    'https://www.sais.gov.uk/api/forecasts',
    'https://www.sais.gov.uk/api/',
    'https://data.gov.uk/dataset/?q=avalanche',
]
for url in sais_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
            print("  {} -> {} bytes".format(url, len(data)))
    except Exception as e:
        print("  {} -> {}".format(url, str(e)[:80]))

# ============================================================
# 5. ICELAND AVALANCHE DATA
# ============================================================
print("\n=== Iceland IMO ===")
ice_urls = [
    'https://en.vedur.is/avalanches/imo-avalanche-database/',
    'https://gatt.vedur.is/api/',
]
for url in ice_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
            print("  {} -> {} bytes".format(url, len(data)))
    except Exception as e:
        print("  {} -> {}".format(url, str(e)[:80]))
