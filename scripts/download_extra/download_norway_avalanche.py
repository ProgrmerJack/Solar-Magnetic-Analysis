"""
Download Norwegian avalanche danger ratings from NVE API.
Mainland Norway mountain regions, winters 2012/13 - 2024/25.
"""
import json, os, time, urllib.request
from datetime import datetime

OUT_DIR = r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\norway_nve"
os.makedirs(OUT_DIR, exist_ok=True)

# Mainland Norwegian avalanche forecast regions (IDs discovered from API)
# We'll scan a range and keep those that return data
REGIONS_TO_TRY = list(range(3003, 3050))

# Winter seasons to download (Nov-Apr)
SEASONS = []
for y in range(2012, 2025):
    SEASONS.append((f"{y}-11-01", f"{y+1}-04-30", f"{y}-{y+1}"))

BASE = "https://api01.nve.no/hydrology/forecast/avalanche/v6.2.1/api"

def fetch_json(url, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header('Accept', 'application/json')
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode('utf-8'))
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2)
            else:
                return None

# Step 1: Discover active regions by probing winter 2019-2020
print("=== Discovering active Norwegian avalanche regions ===")
active_regions = {}
for rid in REGIONS_TO_TRY:
    data = fetch_json(f"{BASE}/AvalancheWarningByRegion/Simple/{rid}/1/2019-12-01/2019-12-31")
    if data and len(data) > 0:
        name = data[0].get('RegionName', f'Region_{rid}')
        active_regions[rid] = name
        print(f"  Found: {rid} = {name} ({len(data)} days)")
    time.sleep(0.3)  # Be polite to API

print(f"\nActive regions: {len(active_regions)}")

# Save region list
with open(os.path.join(OUT_DIR, "regions.json"), 'w') as f:
    json.dump(active_regions, f, indent=2)

# Step 2: Download all seasons for all active regions
all_data = []
for rid, rname in active_regions.items():
    for start, end, season_label in SEASONS:
        data = fetch_json(f"{BASE}/AvalancheWarningByRegion/Simple/{rid}/1/{start}/{end}")
        if data and len(data) > 0:
            for rec in data:
                rec['Season'] = season_label
            all_data.extend(data)
            print(f"  {rname} ({rid}) {season_label}: {len(data)} days")
        time.sleep(0.3)

print(f"\nTotal records downloaded: {len(all_data)}")

# Save raw data
with open(os.path.join(OUT_DIR, "nve_danger_all.json"), 'w', encoding='utf-8') as f:
    json.dump(all_data, f, indent=1, ensure_ascii=False)

# Convert to simplified CSV
import csv
csv_path = os.path.join(OUT_DIR, "nve_danger_all.csv")
with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['date', 'region_id', 'region_name', 'danger_level', 'season'])
    for rec in all_data:
        date = rec['ValidFrom'][:10]
        writer.writerow([date, rec['RegionId'], rec['RegionName'], 
                        rec['DangerLevel'], rec['Season']])

print(f"Saved CSV: {csv_path}")
print("Done!")
