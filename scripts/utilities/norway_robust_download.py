"""
Robust NVE download with retries. Downloads SSW windows + controls 
for 5 key mountain regions with 60s timeout and 5 retries.
"""
import json, os, urllib.request, time, csv
import numpy as np
from datetime import datetime, timedelta
from collections import defaultdict

OUT = r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\norway_nve"
os.makedirs(OUT, exist_ok=True)

BASE = 'https://api01.nve.no/hydrology/forecast/avalanche/v6.2.1/api'

# 5 key inland mountain regions
REGIONS = {3010: "Lyngen", 3022: "Trollheimen", 3025: "Nord-Gudbrandsdalen",
           3028: "Jotunheimen", 3034: "Hardanger"}

def fetch(url, retries=5, timeout=60):
    for i in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header('Accept', 'application/json')
            resp = urllib.request.urlopen(req, timeout=timeout)
            return json.loads(resp.read().decode('utf-8'))
        except Exception as e:
            print(f"    Retry {i+1}/{retries}: {e}")
            time.sleep(3 * (i+1))
    return None

# Windows to download: 2-month blocks around SSW events + controls
WINDOWS = [
    # SSW events
    ("2017-12-15", "2018-03-15", "2018-02-12", "ssw"),
    ("2018-12-01", "2019-02-28", "2019-01-02", "ssw"),
    ("2020-11-15", "2021-02-15", "2021-01-05", "ssw"),
    ("2022-12-15", "2023-03-31", "2023-02-16", "ssw"),
    # Control periods (no SSW)
    ("2019-12-01", "2020-02-28", "control", "ctrl"),
    ("2021-12-01", "2022-02-28", "control", "ctrl"),
    ("2023-12-01", "2024-02-28", "control", "ctrl"),
]

all_records = []
for start, end, ssw_date, wtype in WINDOWS:
    print(f"\n{'SSW '+ssw_date if wtype=='ssw' else 'Control '+start[:7]}: {start} to {end}")
    for rid, rname in REGIONS.items():
        url = f"{BASE}/AvalancheWarningByRegion/Simple/{rid}/1/{start}/{end}"
        data = fetch(url)
        if data and len(data) > 0:
            for rec in data:
                all_records.append({
                    'date': rec['ValidFrom'][:10],
                    'region_id': rec['RegionId'],
                    'region_name': rec['RegionName'],
                    'danger_level': int(rec['DangerLevel']),
                    'ssw_date': ssw_date,
                    'window_type': wtype,
                })
            print(f"  {rname}: {len(data)} days")
        else:
            print(f"  {rname}: NO DATA")
        time.sleep(1)

print(f"\nTotal records: {len(all_records)}")

# Save
csv_path = os.path.join(OUT, "nve_ssw_analysis.csv")
with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['date','region_id','region_name','danger_level','ssw_date','window_type'])
    w.writeheader()
    w.writerows(all_records)
print(f"Saved: {csv_path}")

# ===== ANALYSIS =====
if len(all_records) < 50:
    print("Insufficient data for analysis")
    exit()

print("\n\n========== SSW vs CONTROL ANALYSIS ==========")

# Build daily mean danger
daily = defaultdict(list)
for r in all_records:
    daily[(r['date'], r['window_type'], r['ssw_date'])].append(r['danger_level'])

# Per-SSW event analysis
ssw_events = [w[2] for w in WINDOWS if w[3] == 'ssw']
for ssw_date in ssw_events:
    ssw_dt = datetime.strptime(ssw_date, "%Y-%m-%d")
    
    ssw_vals = []
    for d in range(-15, 16):
        dt = ssw_dt + timedelta(days=d)
        key = dt.strftime("%Y-%m-%d")
        for rec in all_records:
            if rec['date'] == key and rec['ssw_date'] == ssw_date:
                ssw_vals.append(rec['danger_level'])
    
    if ssw_vals:
        print(f"SSW {ssw_date}: mean danger in ±15d window = {np.mean(ssw_vals):.3f} (n={len(ssw_vals)})")

# Aggregate SSW vs control
ssw_all = [r['danger_level'] for r in all_records if r['window_type'] == 'ssw']
ctrl_all = [r['danger_level'] for r in all_records if r['window_type'] == 'ctrl']

if ssw_all and ctrl_all:
    from scipy import stats
    print(f"\nSSW periods: mean={np.mean(ssw_all):.3f} (n={len(ssw_all)})")
    print(f"Control periods: mean={np.mean(ctrl_all):.3f} (n={len(ctrl_all)})")
    print(f"Difference: {np.mean(ssw_all)-np.mean(ctrl_all):+.3f}")
    
    u, p = stats.mannwhitneyu(ssw_all, ctrl_all, alternative='two-sided')
    print(f"Mann-Whitney: U={u:.0f}, P={p:.4f}")
    
    # Chi-square on danger level distribution
    from collections import Counter
    ssw_counts = Counter(ssw_all)
    ctrl_counts = Counter(ctrl_all)
    print(f"\nDanger level distribution:")
    print(f"  Level | SSW (%)  | Control (%)")
    for lvl in range(1, 6):
        s_pct = ssw_counts.get(lvl, 0) / len(ssw_all) * 100
        c_pct = ctrl_counts.get(lvl, 0) / len(ctrl_all) * 100
        print(f"    {lvl}   | {s_pct:6.1f}  | {c_pct:6.1f}")

# Phase-resolved
print("\n=== Phase-resolved (relative to SSW onset) ===")
for phase, d0, d1 in [("Pre d-15..d-1", -15, -1), ("Post d0..d+15", 0, 15), ("Late d+16..d+30", 16, 30)]:
    phase_vals = []
    for ssw_date in ssw_events:
        ssw_dt = datetime.strptime(ssw_date, "%Y-%m-%d")
        for d in range(d0, d1+1):
            dt = ssw_dt + timedelta(days=d)
            key = dt.strftime("%Y-%m-%d")
            for rec in all_records:
                if rec['date'] == key and rec['ssw_date'] == ssw_date:
                    phase_vals.append(rec['danger_level'])
    
    if phase_vals and ctrl_all:
        u, p = stats.mannwhitneyu(phase_vals, ctrl_all, alternative='two-sided')
        print(f"  {phase}: mean={np.mean(phase_vals):.3f} (n={len(phase_vals)}) vs ctrl={np.mean(ctrl_all):.3f}, MW P={p:.4f}")

print("\nDone!")
