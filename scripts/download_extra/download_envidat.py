"""
Download expanded avalanche datasets from EnviDat:
1. Re-analyzed Swiss danger levels (2001-2023) — 22 seasons!
2. Weather/snowpack/danger data for predictions
3. Observational data (stability tests, danger signs)
"""
import requests
import os

datasets = {
    'swiss_danger_2020_2023': {
        'url': 'https://www.envidat.ch/dataset/6fa1a7f9-dea0-4bf7-aa98-a9d5d1e41625/resource/1f7aa436-aa12-488d-9897-217b8f467670/download/data_2020-2021_2022-2023.csv',
        'file': 'swiss_danger_2020_2023.csv',
        'desc': 'Re-analyzed Swiss danger levels 2020-2023'
    },
    'swiss_danger_2001_2020': {
        'url': 'https://www.envidat.ch/dataset/6fa1a7f9-dea0-4bf7-aa98-a9d5d1e41625/resource/b4197f26-0e2a-4831-85dc-f3e66bfffe92/download/data_2001-2002_2019-2020.csv',
        'file': 'swiss_danger_2001_2020.csv',
        'desc': 'Re-analyzed Swiss danger levels 2001-2020'
    },
    'weather_snowpack_danger': {
        'url': 'https://www.envidat.ch/dataset/54d7dc19-e70c-4998-bfb6-caffe41c83e6/resource/11c8a22e-3e76-4f1f-8880-17ebf6525e9a/download/Data_weather_snowpack_danger_forecast.csv',
        'file': 'weather_snowpack_danger_forecast.csv',
        'desc': 'Weather/snowpack/danger for ML predictions'
    },
    'human_triggered_avalanches': {
        'url': 'https://www.envidat.ch/dataset/dfacb44d-49d2-4841-b3ee-a4e9eee8d22f/resource/5fee4967-1a72-4ef1-be4a-cfdbe9184802/download/human_triggered_avalanches.csv',
        'file': 'human_triggered_avalanches.csv',
        'desc': 'Human-triggered avalanche observations'
    },
    'danger_signs': {
        'url': 'https://www.envidat.ch/dataset/dfacb44d-49d2-4841-b3ee-a4e9eee8d22f/resource/576b0dd4-fd5a-4c25-998c-e50ef630428d/download/danger_signs.csv',
        'file': 'danger_signs.csv',
        'desc': 'Observed danger signs'
    },
}

outdir = r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\envidat'
os.makedirs(outdir, exist_ok=True)

for key, ds in datasets.items():
    outpath = os.path.join(outdir, ds['file'])
    if os.path.exists(outpath):
        print(f"[SKIP] {ds['desc']} already exists")
        continue
    
    print(f"[DOWNLOAD] {ds['desc']}...")
    try:
        r = requests.get(ds['url'], timeout=60, stream=True)
        if r.status_code == 200:
            with open(outpath, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            size = os.path.getsize(outpath)
            print(f"  -> {outpath} ({size:,} bytes)")
        else:
            print(f"  ERROR: HTTP {r.status_code}")
    except Exception as e:
        print(f"  ERROR: {e}")

print("\nDone! Checking files:")
for f in os.listdir(outdir):
    path = os.path.join(outdir, f)
    size = os.path.getsize(path)
    print(f"  {f}: {size:,} bytes")
