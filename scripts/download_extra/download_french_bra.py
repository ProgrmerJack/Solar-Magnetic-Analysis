"""
Download French Météo-France BRA (Bulletin d'estimation du Risque d'Avalanche) data.
Try multiple public endpoints to get daily danger levels for Alpine massifs.
"""
import requests
import json
import os
import time
from datetime import datetime, timedelta

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'cryosphere', 'france', 'bra_data')
os.makedirs(OUT_DIR, exist_ok=True)

# Alpine massif IDs from Météo-France
ALPINE_MASSIFS = {
    1: 'Chablais', 2: 'Mont-Blanc', 3: 'Aravis', 4: 'Chartreuse',
    5: 'Belledonne', 6: 'Grandes-Rousses', 7: 'Vercors', 8: 'Oisans',
    9: 'Haute-Tarentaise', 10: 'Beaufortain', 11: 'Bauges',
    12: 'Vanoise', 13: 'Haute-Maurienne', 14: 'Maurienne',
    15: 'Ubaye', 16: 'Devoluy', 17: 'Champsaur',
    18: 'Embrunais-Parpaillon', 19: 'Queyras', 20: 'Thabor',
    21: 'Pelvoux', 22: 'Mercantour', 23: 'Haut-Var/Haut-Verdon'
}

# SSW events that overlap with French BRA data (2016+)
SSW_DATES = ['2018-02-12', '2019-01-02', '2021-01-05', '2023-02-16']

def try_donneespubliques():
    """Try Météo-France public data endpoint."""
    print("=== Trying donneespubliques.meteofrance.fr ===")
    
    # Try the BRA archive endpoint
    urls_to_try = [
        'https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA/',
        'https://donneespubliques.meteofrance.fr/client/gfx/BRA/',
        'https://donneespubliques.meteofrance.fr/api/BRA/',
    ]
    
    for url in urls_to_try:
        try:
            r = requests.get(url, timeout=15, allow_redirects=True)
            print(f"  {url}: status={r.status_code}, size={len(r.content)}")
            if r.status_code == 200:
                print(f"  Content preview: {r.text[:300]}")
                return True
        except Exception as e:
            print(f"  {url}: {e}")
    return False

def try_data_gouv():
    """Search data.gouv.fr for BRA datasets."""
    print("\n=== Searching data.gouv.fr ===")
    
    url = "https://www.data.gouv.fr/api/1/datasets/"
    params = {'q': 'bulletin estimation risque avalanche', 'page_size': 10}
    
    try:
        r = requests.get(url, params=params, timeout=15)
        if r.status_code == 200:
            data = r.json()
            print(f"  Found {data['total']} datasets")
            for ds in data.get('data', [])[:5]:
                print(f"  - {ds['title']}")
                print(f"    URL: {ds['page']}")
                for res in ds.get('resources', [])[:3]:
                    print(f"    Resource: {res['title']} ({res.get('format', '?')}) -> {res['url']}")
            return data
        else:
            print(f"  Status: {r.status_code}")
    except Exception as e:
        print(f"  Error: {e}")
    return None

def try_meteofrance_open_data():
    """Try Météo-France open data portal."""
    print("\n=== Trying Météo-France Open Data ===")
    
    # Météo-France has been transitioning to a new open data portal
    urls = [
        'https://public-api.meteofrance.fr/public/DPBRA/v1/liste-massif',
        'https://public-api.meteofrance.fr/public/DPBRA/v1/',
        'https://rpcache-aa.meteofrance.com/internet2018client/2.0/report?domain=BRA&report_type=BULLETINBRA&report_subtype=BRA&massif=MONT-BLANC',
    ]
    
    headers = {'Accept': 'application/json'}
    
    for url in urls:
        try:
            r = requests.get(url, headers=headers, timeout=15, allow_redirects=True)
            print(f"  {url}: status={r.status_code}")
            if r.status_code == 200:
                ct = r.headers.get('content-type', '')
                if 'json' in ct:
                    print(f"  JSON: {json.dumps(r.json(), indent=2)[:500]}")
                else:
                    print(f"  Content ({ct}): {r.text[:300]}")
                return r
        except Exception as e:
            print(f"  {url}: {e}")
    return None

def try_meteofrance_bra_xml(massif_id, date_str):
    """Try to get BRA XML for a specific massif and date."""
    
    # Various known URL patterns for BRA data
    urls = [
        f'https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA/BRA.{massif_id}.{date_str.replace("-","")}.xml',
        f'https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA/BRA_{massif_id}_{date_str.replace("-","")}.xml',
        f'https://rpcache-aa.meteofrance.com/internet2018client/2.0/report?domain=BRA&report_type=BULLETINBRA&report_subtype=BRA&massif={ALPINE_MASSIFS.get(massif_id, "")}&date={date_str}',
    ]
    
    for url in urls:
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200 and len(r.content) > 100:
                return r.content
        except:
            pass
    return None

def try_meteofrance_rpcache():
    """Try the rpcache endpoint used by the Météo-France website."""
    print("\n=== Trying rpcache endpoint ===")
    
    # This is the endpoint used by the actual Météo-France website
    base = 'https://rpcache-aa.meteofrance.com/internet2018client/2.0'
    
    # Try getting massif list
    try:
        r = requests.get(f'{base}/report?domain=BRA&report_type=BULLETINBRA&report_subtype=BRA_LIST',
                        timeout=15,
                        headers={'User-Agent': 'Mozilla/5.0'})
        print(f"  BRA_LIST: status={r.status_code}")
        if r.status_code == 200:
            print(f"  Content: {r.text[:500]}")
    except Exception as e:
        print(f"  BRA_LIST error: {e}")
    
    # Try specific massif BRA
    massifs_to_try = ['MONT-BLANC', 'BELLEDONNE', 'VANOISE', 'OISANS']
    for massif in massifs_to_try:
        try:
            url = f'{base}/report?domain=BRA&report_type=BULLETINBRA&report_subtype=BRA&massif={massif}'
            r = requests.get(url, timeout=10, headers={'User-Agent': 'Mozilla/5.0'})
            print(f"  {massif}: status={r.status_code}, size={len(r.content)}")
            if r.status_code == 200 and len(r.content) > 200:
                # Save this
                fname = os.path.join(OUT_DIR, f'bra_latest_{massif}.xml')
                with open(fname, 'wb') as f:
                    f.write(r.content)
                print(f"  Saved: {fname}")
                # Print first few lines
                print(f"  Preview: {r.text[:300]}")
        except Exception as e:
            print(f"  {massif}: {e}")

def try_meteofrance_api_nonauthenticated():
    """Try the public non-authenticated Météo-France API endpoints."""
    print("\n=== Trying public MF API endpoints ===")
    
    urls = [
        # New public API (launched ~2024)
        'https://public-api.meteofrance.fr/public/DPBRA/v1/massif/liste',
        'https://public-api.meteofrance.fr/public/DPBRA/v1/massif/2',  # Mont-Blanc
        # Try without auth
        'https://webservice.meteofrance.com/bra?massif=MONT-BLANC',
    ]
    
    for url in urls:
        try:
            r = requests.get(url, timeout=10, 
                           headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'})
            print(f"  {url}: status={r.status_code}")
            if r.status_code == 200:
                print(f"  Content: {r.text[:500]}")
        except Exception as e:
            print(f"  {url}: {e}")

def download_from_data_gouv_bra():
    """Check if data.gouv.fr has downloadable BRA archives."""
    print("\n=== Checking data.gouv.fr for BRA archives ===")
    
    # Known dataset IDs related to avalanche bulletins
    dataset_urls = [
        'https://www.data.gouv.fr/api/1/datasets/?q=BRA+avalanche&page_size=20',
        'https://www.data.gouv.fr/api/1/datasets/?q=risque+avalanche+bulletin&page_size=20',
        'https://www.data.gouv.fr/api/1/datasets/?q=neige+avalanche+meteo+france&page_size=20',
    ]
    
    all_resources = []
    seen_ids = set()
    
    for api_url in dataset_urls:
        try:
            r = requests.get(api_url, timeout=15)
            if r.status_code == 200:
                for ds in r.json().get('data', []):
                    if ds['id'] not in seen_ids:
                        seen_ids.add(ds['id'])
                        print(f"\n  Dataset: {ds['title']}")
                        print(f"  ID: {ds['id']}")
                        for res in ds.get('resources', [])[:5]:
                            fmt = res.get('format', '?')
                            size = res.get('filesize', 0)
                            print(f"    [{fmt}] {res['title']} ({size} bytes)")
                            print(f"    URL: {res['url']}")
                            if fmt in ('xml', 'json', 'csv', 'zip', 'XML', 'JSON', 'CSV'):
                                all_resources.append({
                                    'dataset': ds['title'],
                                    'resource': res['title'],
                                    'format': fmt,
                                    'url': res['url'],
                                    'size': size
                                })
        except Exception as e:
            print(f"  Error querying: {e}")
    
    # Try downloading the most promising resources
    if all_resources:
        print(f"\n  Found {len(all_resources)} downloadable resources")
        for res in all_resources[:3]:
            try:
                print(f"\n  Downloading: {res['resource']} ({res['format']})")
                r = requests.get(res['url'], timeout=30, stream=True)
                if r.status_code == 200:
                    ext = res['format'].lower()
                    safe_name = res['resource'].replace('/', '_').replace('\\', '_')[:50]
                    fname = os.path.join(OUT_DIR, f'{safe_name}.{ext}')
                    with open(fname, 'wb') as f:
                        for chunk in r.iter_content(8192):
                            f.write(chunk)
                    fsize = os.path.getsize(fname)
                    print(f"  Saved: {fname} ({fsize} bytes)")
                else:
                    print(f"  Status: {r.status_code}")
            except Exception as e:
                print(f"  Download error: {e}")
    
    return all_resources

if __name__ == '__main__':
    print("French BRA Data Download")
    print("="*60)
    
    # Try all endpoints
    try_donneespubliques()
    try_data_gouv()
    try_meteofrance_open_data()
    try_meteofrance_rpcache()
    try_meteofrance_api_nonauthenticated()
    download_from_data_gouv_bra()
    
    print("\n" + "="*60)
    print("Download attempts complete. Check data/cryosphere/france/bra_data/")
