"""
Download public SLF/WSL avalanche datasets from EnviDat (Swiss Environmental Data)
No authentication required — all datasets are CC-BY licensed.
"""
import requests
import json
from pathlib import Path

OUT = Path('C:/Users/Jack0/Solar-Magnetic-Analysis/data/cryosphere/slf_avalanche')
OUT.mkdir(parents=True, exist_ok=True)

BASE = 'https://www.envidat.ch'
sess = requests.Session()
sess.headers['User-Agent'] = 'Solar-Magnetic-Analysis/1.0 (academic research)'

# High-value public SLF datasets on EnviDat
DATASETS = [
    'avalanche-accidents-in-switzerland-since-1970-71',   # PRIMARY: accidents since 1970
    'snow-avalanche-data-davos',                           # Davos 1999-2019 event data
    'data_wet_aval_model',                                 # weather+snowpack+avalanche for ML model
    'stability-tests-avalanche-observations-switzerland-norway',
    'simulated-avalanche-problem-types-at-weissfluhjoch-1999-2017',
]

for ds_id in DATASETS:
    r = sess.get(BASE + '/api/3/action/package_show?id=' + ds_id, timeout=15)
    info = r.json().get('result', {})
    title = info.get('title', ds_id)
    resources = info.get('resources', [])
    print('Dataset:', title)
    print('  Resources:', len(resources))

    for res in resources:
        url = res.get('url', '')
        fmt = res.get('format', 'dat').lower()
        name = res.get('name', res.get('id', 'file'))[:50].replace('/', '_').replace(' ', '_')
        if not url:
            continue
        fname = OUT / (ds_id[:25] + '_' + name + '.' + fmt)
        if fname.exists() and fname.stat().st_size > 100:
            print('    skip (exists):', fname.name)
            continue
        rdata = sess.get(url, timeout=120)
        if rdata.status_code == 200:
            fname.write_bytes(rdata.content)
            print('    saved:', fname.name, '(' + str(len(rdata.content)//1024) + ' KB)')
        else:
            print('    HTTP', rdata.status_code, ':', url[:70])
    print()

print('EnviDat SLF download complete. Output:', OUT)
