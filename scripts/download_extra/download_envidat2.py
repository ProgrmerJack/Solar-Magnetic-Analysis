"""Download Swiss re-analyzed danger levels from EnviDat."""
import requests, os

outdir = r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\envidat'
os.makedirs(outdir, exist_ok=True)

datasets = [
    ('https://www.envidat.ch/dataset/6fa1a7f9-dea0-4bf7-aa98-a9d5d1e41625/resource/b4197f26-0e2a-4831-85dc-f3e66bfffe92/download/dangerleveldataswitzerland_tidy_ft.csv',
     'swiss_danger_2001_2020.csv'),
    ('https://www.envidat.ch/dataset/6fa1a7f9-dea0-4bf7-aa98-a9d5d1e41625/resource/1f7aa436-aa12-488d-9897-217b8f467670/download/dangerlevel_tidy_2020-2021_2022-2023.csv',
     'swiss_danger_2020_2023.csv'),
    ('https://www.envidat.ch/dataset/6fa1a7f9-dea0-4bf7-aa98-a9d5d1e41625/resource/a6564128-022a-421e-bb43-7caefd92bcee/download/data_2023-2024.csv',
     'swiss_danger_2023_2024.csv'),
]

for url, fname in datasets:
    outpath = os.path.join(outdir, fname)
    print(f"Downloading {fname}...")
    try:
        r = requests.get(url, timeout=60, stream=True)
        print(f"  HTTP {r.status_code}")
        if r.status_code == 200:
            with open(outpath, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            size = os.path.getsize(outpath)
            print(f"  -> {size:,} bytes")
        else:
            print(f"  Failed: {r.text[:200]}")
    except Exception as e:
        print(f"  ERROR: {e}")

# Also download the weather/snowpack data
weather_url = 'https://www.envidat.ch/dataset/54d7dc19-e70c-4998-bfb6-caffe41c83e6/resource/11c8a22e-3e76-4f1f-8880-17ebf6525e9a/download/Data_weather_snowpack_danger_forecast.csv'
outpath = os.path.join(outdir, 'weather_snowpack_danger.csv')
print(f"\nDownloading weather/snowpack data...")

# First get actual filename
r2 = requests.get('https://www.envidat.ch/api/3/action/package_show?id=54d7dc19-e70c-4998-bfb6-caffe41c83e6', timeout=15)
data = r2.json()
for res in data.get('result', {}).get('resources', []):
    name = res.get('name', '')
    url = res.get('url', '')
    fmt = res.get('format', '')
    print(f"  Resource: {name} ({fmt})")
    print(f"  URL: {url}")
    if fmt == 'CSV' and 'weather' in name.lower():
        print(f"  Downloading...")
        r3 = requests.get(url, timeout=120, stream=True)
        if r3.status_code == 200:
            with open(outpath, 'wb') as f:
                for chunk in r3.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"  -> {os.path.getsize(outpath):,} bytes")

# Also get observational data
r4 = requests.get('https://www.envidat.ch/api/3/action/package_show?id=dfacb44d-49d2-4841-b3ee-a4e9eee8d22f', timeout=15)
data4 = r4.json()
for res in data4.get('result', {}).get('resources', []):
    name = res.get('name', '')
    url = res.get('url', '')
    fmt = res.get('format', '')
    if fmt == 'CSV':
        fname = name.replace(' ', '_').lower()[:40] + '.csv'
        outpath = os.path.join(outdir, fname)
        print(f"\nDownloading {name}...")
        r5 = requests.get(url, timeout=60, stream=True)
        if r5.status_code == 200:
            with open(outpath, 'wb') as f:
                for chunk in r5.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"  -> {os.path.getsize(outpath):,} bytes")

print("\nAll files:")
for f in os.listdir(outdir):
    size = os.path.getsize(os.path.join(outdir, f))
    print(f"  {f}: {size:,} bytes")
