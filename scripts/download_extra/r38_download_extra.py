"""
Try alternative data sources to expand multi-country evidence:
1. MeteoSwiss/IDAWEB station data (open data portal)
2. Austrian ZAMG avalanche warnings (Lawinen.report API)
3. French Météo-France/ANENA (data.gouv.fr)
4. Japanese NIED avalanche (open data)
5. Icelandic Met Office
6. Canadian avalanche (alternative URLs)
7. Scottish SAIS (alternative URLs)
8. Catalan ICGC avalanche database
"""
import urllib.request, urllib.error, json, os

os.makedirs('data/cryosphere/multi_country', exist_ok=True)

# Track results
results = []

urls = {
    # Austrian Lawinen.report (successor to ALBINA)
    'austria_lawinen_api': 'https://api.avalanche.report/albina/api/bulletins?date=2024-02-01&lang=en',
    'austria_lawinen_v2': 'https://avalanche.report/albina/api/bulletins?date=2024-01-15&lang=en',
    'austria_lawinen_caaml': 'https://avalanche.report/caaml/en/latest',
    
    # French data.gouv.fr (ANENA)
    'france_anena_datagouv': 'https://www.data.gouv.fr/api/1/datasets/?q=avalanche+france',
    'france_meteo_montagne': 'https://donneespubliques.meteofrance.fr/donnees_libres/Txt/BRA/',
    
    # Norwegian (alternative API endpoints)
    'norway_varsom_v2': 'https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/AvalancheWarningByRegion/Detail/3003/1/2024-01-01/2024-03-31',
    'norway_varsom_regions': 'https://api01.nve.no/hydrology/forecast/avalanche/v6.3.0/api/Region',
    
    # Japanese NIED snow disaster
    'japan_nied_snow': 'https://www.bosai.go.jp/snow/snow_disaster/download.html',
    
    # Canadian Avalanche Association
    'canada_avalanche_api': 'https://www.avalanche.ca/api/min/en/submissions?period=2024',
    'canada_cac_forecasts': 'https://www.avalanche.ca/api/forecasts/en/archive',
    
    # Icelandic Met Office
    'iceland_vedur': 'https://api.vedur.is/avalanches/v1/warnings',
    
    # Catalan ICGC
    'catalan_icgc': 'https://www.icgc.cat/en/Public-Administration-and-Enterprises/Downloads/Avalanches',
    
    # Swiss additional (WSL/SLF)
    'slf_opendata': 'https://www.slf.ch/en/avalanche-bulletin-and-snow-situation/measured-values.html',
    'envidat_datasets': 'https://www.envidat.ch/api/3/action/package_search?q=avalanche&rows=20',
    
    # EAWS API v2
    'eaws_bulletins_api': 'https://api.avalanche.report/v2/bulletins/latest',
    'eaws_regions': 'https://avalanches.org/wp-json/eaws/v1/regions',
    
    # MeteoSwiss IDAWEB
    'meteoswiss_opendata': 'https://data.geo.admin.ch/ch.meteoschweiz.klima/nbcn-daily/list.json',
}

for name, url in urls.items():
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 Research'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = resp.read()
            status = resp.status
            ct = resp.headers.get('Content-Type', '')
            size = len(data)
            
            # Save if substantial
            if size > 100:
                ext = '.json' if 'json' in ct else '.html' if 'html' in ct else '.txt'
                fpath = f'data/cryosphere/multi_country/{name}{ext}'
                with open(fpath, 'wb') as f:
                    f.write(data)
                result = f"SUCCESS: {status}, {size:,} bytes -> {fpath}"
            else:
                result = f"SUCCESS but tiny: {status}, {size} bytes"
            
    except urllib.error.HTTPError as e:
        result = f"HTTP {e.code}: {e.reason}"
    except urllib.error.URLError as e:
        result = f"URL Error: {e.reason}"
    except Exception as e:
        result = f"Error: {str(e)[:80]}"
    
    results.append((name, result))
    print(f"  {name}: {result}")

# Summary
print(f"\n{'='*60}")
successes = [r for r in results if r[1].startswith('SUCCESS')]
print(f"Successful: {len(successes)}/{len(results)}")
for name, result in successes:
    print(f"  ✓ {name}: {result}")
