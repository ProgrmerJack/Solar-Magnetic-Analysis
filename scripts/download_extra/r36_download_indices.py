"""
Download NOAA CPC blocking index + NAO/AO indices for mechanism analysis.
Also try downloading Japanese avalanche data and NOAA blocking patterns.
"""
import urllib.request, os, json

os.makedirs('data/atmospheric', exist_ok=True)

downloads = {}

# 1. NOAA CPC Daily NAO Index
try:
    url = 'https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii'
    urllib.request.urlretrieve(url, 'data/atmospheric/nao_monthly.txt')
    downloads['nao_monthly'] = 'success'
    print('NAO monthly: OK')
except Exception as e:
    downloads['nao_monthly'] = str(e)
    print(f'NAO monthly: FAILED - {e}')

# 2. Daily NAO index
try:
    url = 'https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/nao.shtml'
    # Actually the daily index is at:
    url = 'https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.nao.index.b500101.current.ascii'
    urllib.request.urlretrieve(url, 'data/atmospheric/nao_daily.txt')
    downloads['nao_daily'] = 'success'
    print('NAO daily: OK')
except Exception as e:
    downloads['nao_daily'] = str(e)
    print(f'NAO daily: FAILED - {e}')

# 3. Daily AO index  
try:
    url = 'https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.ao.index.b500101.current.ascii'
    urllib.request.urlretrieve(url, 'data/atmospheric/ao_daily.txt')
    downloads['ao_daily'] = 'success'
    print('AO daily: OK')
except Exception as e:
    downloads['ao_daily'] = str(e)
    print(f'AO daily: FAILED - {e}')

# 4. NOAA Blocking Index (Tibaldi-Molteni style)
try:
    url = 'https://ftp.cpc.ncep.noaa.gov/data/indices/teleconnections/blocking.index.dat'
    urllib.request.urlretrieve(url, 'data/atmospheric/blocking_index.dat')
    downloads['blocking'] = 'success'
    print('Blocking index: OK')
except Exception as e:
    downloads['blocking'] = str(e)
    print(f'Blocking index: FAILED - {e}')

# 5. Euro-Atlantic blocking from NOAA
try:
    url = 'https://ftp.cpc.ncep.noaa.gov/data/indices/teleconnections/blocks.dat'
    urllib.request.urlretrieve(url, 'data/atmospheric/blocks.dat')
    downloads['blocks'] = 'success'
    print('Blocks: OK')
except Exception as e:
    downloads['blocks'] = str(e)
    print(f'Blocks: FAILED - {e}')

# 6. PNA index (Pacific-North American)
try:
    url = 'https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.pna.index.b500101.current.ascii'
    urllib.request.urlretrieve(url, 'data/atmospheric/pna_daily.txt')
    downloads['pna_daily'] = 'success'
    print('PNA daily: OK')
except Exception as e:
    downloads['pna_daily'] = str(e)
    print(f'PNA daily: FAILED - {e}')

# 7. Euro-Atlantic teleconnections (EA, EAWR, SCAND patterns)
for pattern, fname in [('ea', 'ea'), ('eawr', 'eawr'), ('scand', 'scand')]:
    try:
        url = f'https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.{pattern}.index.b500101.current.ascii'
        urllib.request.urlretrieve(url, f'data/atmospheric/{fname}_daily.txt')
        downloads[f'{fname}_daily'] = 'success'
        print(f'{fname.upper()} daily: OK')
    except Exception as e:
        downloads[f'{fname}_daily'] = str(e)
        print(f'{fname.upper()} daily: FAILED - {e}')

# 8. QBO index from FU Berlin
try:
    url = 'https://www.geo.fu-berlin.de/met/ag/strat/produkte/qbo/singapore.dat'
    urllib.request.urlretrieve(url, 'data/atmospheric/qbo_singapore.dat')
    downloads['qbo'] = 'success'
    print('QBO: OK')
except Exception as e:
    downloads['qbo'] = str(e)
    print(f'QBO: FAILED - {e}')

# 9. Try downloading EnviDat additional data (Swiss avalanche occurrence count)
try:
    url = 'https://www.envidat.ch/dataset/avalanche-bulletin-danger-levels-2012-2022/resource/0cad1b4e-edc7-4aec-bb7d-d0da4eeafe90/download/count_data.csv'
    urllib.request.urlretrieve(url, 'data/cryosphere/envidat_count_data.csv')
    downloads['envidat_counts'] = 'success'
    print('EnviDat count data: OK')
except Exception as e:
    downloads['envidat_counts'] = str(e)
    print(f'EnviDat count data: FAILED - {e}')

# 10. Try EAWS API for recent data
try:
    url = 'https://api.avalanche.report/v2/avalanche-report?lang=en&date=2024-01-15'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=15)
    data = resp.read().decode()
    with open('data/cryosphere/eaws_api_test.json', 'w') as f:
        f.write(data)
    downloads['eaws_api'] = 'success'
    print(f'EAWS API test: OK ({len(data)} bytes)')
except Exception as e:
    downloads['eaws_api'] = str(e)
    print(f'EAWS API: FAILED - {e}')

print('\n=== Summary ===')
for k, v in downloads.items():
    status = 'OK' if v == 'success' else 'FAIL'
    print(f'  {k:25s}: {status}')
