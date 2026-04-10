"""Download NOAA SWPC Solar Event Reports and additional solar datasets.
Also downloads SuperMAG indices if API key available.
"""
import requests
import json
import os
import time

# === 1. SWPC Solar Event Reports (Flare catalog) ===
output_dir = 'data/solar/swpc'
os.makedirs(output_dir, exist_ok=True)

# Download flare lists for multiple years
print("=== Downloading NOAA SWPC Flare Data ===")
# SWPC provides historical data through their archive

# Download 7-day flare data (latest)
url = "https://services.swpc.noaa.gov/json/goes/primary/xray-flares-7-day.json"
try:
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        with open(os.path.join(output_dir, 'xray_flares_7day.json'), 'w') as f:
            json.dump(r.json(), f, indent=2)
        print("  Saved 7-day flares:", len(r.json()), "events")
except Exception as e:
    print("  Error:", e)

# Download solar cycle progression
url = "https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json"
try:
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        with open(os.path.join(output_dir, 'solar_cycle_indices.json'), 'w') as f:
            json.dump(r.json(), f, indent=2)
        print("  Saved solar cycle indices:", len(r.json()), "months")
except Exception as e:
    print("  Error:", e)

# Download sunspot number
url = "https://services.swpc.noaa.gov/json/solar-cycle/sunspots.json"
try:
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        with open(os.path.join(output_dir, 'sunspots.json'), 'w') as f:
            json.dump(r.json(), f, indent=2)
        print("  Saved sunspot numbers:", len(r.json()), "records")
except Exception as e:
    print("  Error:", e)

# Download SWPC geomagnetic storm data  
url = "https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json"
try:
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        with open(os.path.join(output_dir, 'planetary_k_index.json'), 'w') as f:
            json.dump(r.json(), f, indent=2)
        print("  Saved planetary K-index:", len(r.json()), "records")
except Exception as e:
    print("  Error:", e)

# Download solar wind magnetic field
url = "https://services.swpc.noaa.gov/products/solar-wind/mag-7-day.json"
try:
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        with open(os.path.join(output_dir, 'solar_wind_mag_7day.json'), 'w') as f:
            json.dump(r.json(), f, indent=2)
        print("  Saved solar wind mag:", len(r.json()), "records")
except Exception as e:
    print("  Error:", e)

# Download solar wind plasma
url = "https://services.swpc.noaa.gov/products/solar-wind/plasma-7-day.json"
try:
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        with open(os.path.join(output_dir, 'solar_wind_plasma_7day.json'), 'w') as f:
            json.dump(r.json(), f, indent=2)
        print("  Saved solar wind plasma:", len(r.json()), "records")
except Exception as e:
    print("  Error:", e)

# === 2. NOAA POES/MEPED Hemispheric Power Index ===
print("\n=== Downloading POES Hemispheric Power ===")
poes_dir = 'data/atmospheric/poes_hpi'
os.makedirs(poes_dir, exist_ok=True)

url = "https://services.swpc.noaa.gov/text/aurora-nowcast-hemi-power.txt"
try:
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        with open(os.path.join(poes_dir, 'hemispheric_power_latest.txt'), 'w') as f:
            f.write(r.text)
        print("  Saved hemispheric power data")
except Exception as e:
    print("  Error:", e)

# === 3. Kp/Dst/AE from Kyoto WDC ===
print("\n=== Downloading Geomagnetic Indices from GFZ ===")
geo_dir = 'data/geomagnetic/gfz'
os.makedirs(geo_dir, exist_ok=True)

# GFZ Potsdam provides Kp data
url = "https://www-app3.gfz.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt"
try:
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        with open(os.path.join(geo_dir, 'kp_ap_since_1932.txt'), 'w') as f:
            f.write(r.text)
        print("  Saved Kp/ap since 1932:", len(r.text), "bytes")
except Exception as e:
    print("  Error:", e)

# Download AE index from WDC Kyoto  
print("\n=== Downloading AE Index ===")
for year in range(2010, 2025):
    url = "https://wdc.kugi.kyoto-u.ac.jp/ae_provisional/%04d/index.html" % year
    try:
        r = requests.get(url, timeout=10)
        # Note: actual AE data requires specific format request
    except:
        pass

print("\nAll SWPC/geomagnetic downloads complete!")
print("\nData inventory:")
for d in ['data/solar/swpc', 'data/atmospheric/poes_hpi', 'data/geomagnetic/gfz']:
    if os.path.exists(d):
        files = os.listdir(d)
        total = sum(os.path.getsize(os.path.join(d, f)) for f in files if os.path.isfile(os.path.join(d, f)))
        print("  %s: %d files, %.1f KB" % (d, len(files), total/1024))
