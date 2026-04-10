"""
Try to find ALBINA static JSON bulletins and avalanche.org observation data.
"""
import urllib.request, json, os, re

# ============================================================
# 1. ALBINA - Try various static file URL patterns
# ============================================================
print("=== ALBINA Static Files ===")
# The ALBINA project stores bulletins as static JSON/CAAML files
# Try various URL patterns for a known date
date = '2024-01-15'
patterns = [
    'https://static.avalanche.report/bulletins/{}/{}.json'.format(date[:4], date),
    'https://static.avalanche.report/bulletins/{}/{}_en.json'.format(date[:4], date),
    'https://static.avalanche.report/bulletins/{}/{}.xml'.format(date[:4], date),
    'https://avalanche.report/content/bulletins/{}/{}.json'.format(date[:4], date),
    'https://avalanche.report/content/{}.json'.format(date),
    'https://avalanche.report/bulletin/{}'.format(date),  
    'https://avalanche.report/content/bulletins/{}.json'.format(date),
    # CAAMLv6 format
    'https://avalanche.report/caaml/en/{}'.format(date),
    # Try with PM suffix (bulletins are typically PM)
    'https://static.avalanche.report/bulletins/{}/{}_{}.json'.format(date[:4], date, '1700'),
]

for url in patterns:
    try:
        req = urllib.request.Request(url, headers={
            'User-Agent': 'Mozilla/5.0',
            'Accept': 'application/json,application/xml,*/*'
        })
        with urllib.request.urlopen(req, timeout=8) as resp:
            data = resp.read()
            ct = resp.headers.get('Content-Type', '')
            print("OK  {} -> {} bytes [{}]".format(url[35:], len(data), ct[:30]))
            text = data.decode('utf-8', errors='replace')[:500]
            if 'json' in ct:
                j = json.loads(data)
                print("    {}".format(str(j)[:300]))
            elif 'xml' in ct:
                print("    {}".format(text[:300]))
            else:
                print("    {}".format(text[:200]))
    except urllib.error.HTTPError as e:
        if e.code != 404:
            print("ERR {} -> {}".format(url[35:], e))
    except Exception as e:
        print("ERR {} -> {}".format(url[35:], str(e)[:60]))

# ============================================================
# 2. ALBINA GitLab API - check for data repository
# ============================================================
print("\n=== ALBINA GitLab ===")
gl_urls = [
    'https://gitlab.com/api/v4/projects/albina-euregio',
    'https://gitlab.com/api/v4/groups/albina-euregio/projects?per_page=20',
]
for url in gl_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
            if isinstance(data, list):
                for p in data:
                    print("  Project: {} ({})".format(p.get('name', ''), p.get('web_url', '')))
            elif isinstance(data, dict):
                print("  {}".format(str(data)[:300]))
    except Exception as e:
        print("  {} -> {}".format(url[35:], str(e)[:60]))

# ============================================================
# 3. Avalanche.org - download actual avalanche OBSERVATIONS
# ============================================================
print("\n=== Avalanche.org Observations ===")
# We know the forecast endpoint works. Let's try to get observation reports
obs_urls = [
    'https://api.avalanche.org/v2/public/products?type=observation&date_start=2024-01-01&date_end=2024-01-15',
    'https://api.avalanche.org/v2/public/obs?date_start=2024-01-01&date_end=2024-01-15',
    'https://api.avalanche.org/v2/public/products?type=avalanche&date_start=2024-01-01&date_end=2024-01-15',
    'https://api.avalanche.org/v2/public/products?type=accident&date_start=2024-01-01&date_end=2024-01-15',
    'https://api.avalanche.org/v2/public/products?type=incident&date_start=2024-01-01&date_end=2024-01-15',
]
for url in obs_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
            if isinstance(data, list):
                print("  {} -> {} items".format(url.split('?')[0].split('/')[-1], len(data)))
                if data and isinstance(data[0], dict):
                    print("    Keys: {}".format(list(data[0].keys())[:8]))
                    # Check product types
                    types = set(d.get('product_type', '') for d in data[:20])
                    print("    Types: {}".format(types))
            elif isinstance(data, dict):
                print("  {} -> dict keys {}".format(url.split('?')[0].split('/')[-1], list(data.keys())[:5]))
    except Exception as e:
        print("  {} -> {}".format(url.split('?')[0].split('/')[-1], str(e)[:60]))

# ============================================================  
# 4. Try Meteo-France avalanche danger API
# ============================================================
print("\n=== Meteo-France BRA API ===")
mf_urls = [
    'https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA/',
    'https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA/BRA.CHABLAIS.20240115140029.pdf',
    'https://donneespubliques.meteofrance.fr/donnees_libres/Txt/BRA/',
    'https://donneespubliques.meteofrance.fr/donnees_libres/Xml/BRA/',
]
for url in mf_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
            ct = resp.headers.get('Content-Type', '')
            print("  {} -> {} bytes [{}]".format(url.split('BRA')[-1] or '/BRA/', len(data), ct[:30]))
            if 'html' in ct or 'text' in ct:
                text = data.decode('utf-8', errors='replace')
                # Look for file links
                links = re.findall(r'href="([^"]*\.(pdf|xml|json|csv|txt))"', text, re.I)
                print("    Files: {}".format(links[:10]))
    except Exception as e:
        print("  {} -> {}".format(url.split('BRA')[-1] or '/BRA/', str(e)[:60]))
