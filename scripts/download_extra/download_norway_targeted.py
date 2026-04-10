"""
Targeted download of Norwegian avalanche danger ratings from NVE API.
Focus on key mountain regions for SSW analysis windows.
"""
import json, os, time, urllib.request, csv
from datetime import datetime, timedelta

OUT_DIR = r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\norway_nve"
os.makedirs(OUT_DIR, exist_ok=True)

BASE = "https://api01.nve.no/hydrology/forecast/avalanche/v6.2.1/api"

# Key mainland mountain regions (skip Svalbard, coastal, and low-altitude areas)
MOUNTAIN_REGIONS = {
    3010: "Lyngen",
    3011: "Tromsø",  
    3012: "Sør-Troms",
    3013: "Indre Troms",
    3015: "Ofoten",
    3022: "Trollheimen",
    3023: "Romsdal",
    3024: "Sunnmøre",
    3025: "Nord-Gudbrandsdalen",
    3027: "Indre Fjordane",
    3028: "Jotunheimen",
    3029: "Indre Sogn",
    3031: "Voss",
    3032: "Hallingdal",
    3034: "Hardanger",
}

# SSW events from catalog (winter months, 2013-2024)
SSW_DATES = [
    "2013-01-06", "2018-02-12", "2019-01-02", "2021-01-05",
    "2023-02-16", "2024-03-04",
]

def fetch_json(url, retries=2):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header('Accept', 'application/json')
            with urllib.request.urlopen(req, timeout=20) as resp:
                return json.loads(resp.read().decode('utf-8'))
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(1)
            else:
                return None

# Download full winter seasons that contain SSW events
all_data = []
seasons_done = set()

for ssw_str in SSW_DATES:
    ssw_dt = datetime.strptime(ssw_str, "%Y-%m-%d")
    # Define winter season
    if ssw_dt.month >= 9:
        season_start = f"{ssw_dt.year}-11-01"
        season_end = f"{ssw_dt.year+1}-04-30"
        season_label = f"{ssw_dt.year}-{ssw_dt.year+1}"
    else:
        season_start = f"{ssw_dt.year-1}-11-01"
        season_end = f"{ssw_dt.year}-04-30"
        season_label = f"{ssw_dt.year-1}-{ssw_dt.year}"
    
    if season_label in seasons_done:
        continue
    seasons_done.add(season_label)
    
    print(f"\n=== Season {season_label} (SSW: {ssw_str}) ===")
    
    for rid, rname in MOUNTAIN_REGIONS.items():
        url = f"{BASE}/AvalancheWarningByRegion/Simple/{rid}/1/{season_start}/{season_end}"
        data = fetch_json(url)
        if data and len(data) > 0:
            for rec in data:
                rec['Season'] = season_label
                rec['SSW_date'] = ssw_str
            all_data.extend(data)
            print(f"  {rname}: {len(data)} days")
        time.sleep(0.2)

# Also download 2 control seasons (no SSW)
print("\n=== Control seasons ===")
for ctrl_season, ctrl_start, ctrl_end in [
    ("2015-2016", "2015-11-01", "2016-04-30"),
    ("2016-2017", "2016-11-01", "2017-04-30"),
    ("2020-2021_ctrl", "2019-11-01", "2020-04-30"),
]:
    print(f"\n  Season {ctrl_season}:")
    for rid, rname in MOUNTAIN_REGIONS.items():
        url = f"{BASE}/AvalancheWarningByRegion/Simple/{rid}/1/{ctrl_start}/{ctrl_end}"
        data = fetch_json(url)
        if data and len(data) > 0:
            for rec in data:
                rec['Season'] = ctrl_season
                rec['SSW_date'] = 'control'
            all_data.extend(data)
            print(f"    {rname}: {len(data)} days")
        time.sleep(0.2)

print(f"\n\nTotal records: {len(all_data)}")

# Save as JSON and CSV
with open(os.path.join(OUT_DIR, "nve_danger_mountain.json"), 'w', encoding='utf-8') as f:
    json.dump(all_data, f, indent=1, ensure_ascii=False)

csv_path = os.path.join(OUT_DIR, "nve_danger_mountain.csv")
with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['date','region_id','region_name','danger_level','season','ssw_date'])
    for rec in all_data:
        writer.writerow([
            rec['ValidFrom'][:10], rec['RegionId'], rec['RegionName'],
            rec['DangerLevel'], rec['Season'], rec.get('SSW_date','')
        ])

print(f"Saved: {csv_path}")
print("Done!")
