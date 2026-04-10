"""
Download Norwegian avalanche danger data and analyze SSW response.
Efficient: downloads only needed windows, analyzes immediately.
"""
import json, os, urllib.request, time, csv
import numpy as np
from datetime import datetime, timedelta
from collections import defaultdict

OUT_DIR = r"C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\norway_nve"
os.makedirs(OUT_DIR, exist_ok=True)

BASE = 'https://api01.nve.no/hydrology/forecast/avalanche/v6.2.1/api'

# Key mountain regions (inland/alpine, not coastal)
REGIONS = {
    3010: "Lyngen", 3011: "Tromsø", 3013: "Indre Troms",
    3015: "Ofoten", 3022: "Trollheimen", 3023: "Romsdal",
    3024: "Sunnmøre", 3025: "Nord-Gudbrandsdalen",
    3028: "Jotunheimen", 3029: "Indre Sogn",
    3031: "Voss", 3032: "Hallingdal", 3034: "Hardanger",
}

# SSW events (winter only, after NVE system started ~2013)
SSW_EVENTS = [
    "2013-01-06", "2018-02-12", "2019-01-02",
    "2021-01-05", "2023-02-16",
]

def fetch_json(url):
    try:
        req = urllib.request.Request(url)
        req.add_header('Accept', 'application/json')
        resp = urllib.request.urlopen(req, timeout=30)
        return json.loads(resp.read().decode('utf-8'))
    except:
        return None

# Step 1: Download full winter for each SSW season + controls
print("=== Downloading Norwegian avalanche danger data ===")
all_records = []

# Determine unique winter seasons
seasons = set()
for ssw in SSW_EVENTS:
    dt = datetime.strptime(ssw, "%Y-%m-%d")
    yr = dt.year if dt.month <= 6 else dt.year + 1
    seasons.add((yr-1, yr))

# Add control seasons
seasons.add((2015, 2016))
seasons.add((2016, 2017))
seasons.add((2014, 2015))

for y1, y2 in sorted(seasons):
    start = f"{y1}-11-01"
    end = f"{y2}-04-30"
    label = f"{y1}-{y2}"
    n_total = 0
    for rid, rname in REGIONS.items():
        url = f"{BASE}/AvalancheWarningByRegion/Simple/{rid}/1/{start}/{end}"
        data = fetch_json(url)
        if data and len(data) > 0:
            for rec in data:
                all_records.append({
                    'date': rec['ValidFrom'][:10],
                    'region_id': rec['RegionId'],
                    'region_name': rec['RegionName'],
                    'danger_level': int(rec['DangerLevel']),
                    'season': label,
                })
            n_total += len(data)
        time.sleep(0.3)
    print(f"  {label}: {n_total} records")

print(f"\nTotal records: {len(all_records)}")

# Save CSV
csv_path = os.path.join(OUT_DIR, "nve_mountain_danger.csv")
with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['date','region_id','region_name','danger_level','season'])
    writer.writeheader()
    writer.writerows(all_records)
print(f"Saved: {csv_path}")

# Step 2: SSW Analysis
print("\n\n========== SSW ANALYSIS ON NORWEGIAN DANGER LEVELS ==========")

# Build daily danger index (mean across all active regions)
from collections import defaultdict
daily_danger = defaultdict(list)
for rec in all_records:
    daily_danger[rec['date']].append(rec['danger_level'])

daily_mean = {}
for date, levels in daily_danger.items():
    daily_mean[date] = np.mean(levels)

dates_sorted = sorted(daily_mean.keys())
print(f"Date range: {dates_sorted[0]} to {dates_sorted[-1]}")
print(f"Total unique dates: {len(dates_sorted)}")

# For each SSW event, compute danger level in ±15-day window vs control
print("\n=== Event-level SSW analysis ===")
ssw_results = []

for ssw_str in SSW_EVENTS:
    ssw_dt = datetime.strptime(ssw_str, "%Y-%m-%d")
    
    # SSW window: d-15 to d+15
    ssw_levels = []
    for d in range(-15, 16):
        dt = ssw_dt + timedelta(days=d)
        key = dt.strftime("%Y-%m-%d")
        if key in daily_mean:
            ssw_levels.append(daily_mean[key])
    
    if len(ssw_levels) < 10:
        print(f"  {ssw_str}: insufficient data ({len(ssw_levels)} days)")
        continue
    
    # Control: same calendar window in control seasons
    ctrl_levels = []
    for ctrl_yr_offset in [-2, -1, 1, 2]:
        for d in range(-15, 16):
            dt = ssw_dt + timedelta(days=d)
            ctrl_dt = dt.replace(year=dt.year + ctrl_yr_offset)
            key = ctrl_dt.strftime("%Y-%m-%d")
            if key in daily_mean:
                ctrl_levels.append(daily_mean[key])
    
    if len(ctrl_levels) < 10:
        # Use seasonal mean as control
        season_yr = ssw_dt.year if ssw_dt.month <= 6 else ssw_dt.year + 1
        season_levels = [v for k, v in daily_mean.items() 
                        if k[:4] in [str(season_yr-1), str(season_yr)]]
        ctrl_mean = np.mean(season_levels) if season_levels else 2.0
    else:
        ctrl_mean = np.mean(ctrl_levels)
    
    ssw_mean = np.mean(ssw_levels)
    diff = ssw_mean - ctrl_mean
    
    ssw_results.append({
        'ssw_date': ssw_str,
        'ssw_mean_danger': ssw_mean,
        'ctrl_mean_danger': ctrl_mean,
        'difference': diff,
        'n_ssw_days': len(ssw_levels),
        'n_ctrl_days': len(ctrl_levels),
        'decrease': diff < 0,
    })
    
    direction = "↓ DECREASE" if diff < 0 else "↑ increase"
    print(f"  {ssw_str}: SSW={ssw_mean:.2f} vs Ctrl={ctrl_mean:.2f} ({diff:+.3f}) {direction}")

# Summary
if ssw_results:
    n_decrease = sum(1 for r in ssw_results if r['decrease'])
    n_total = len(ssw_results)
    
    print(f"\n=== SUMMARY ===")
    print(f"SSW events with decreased danger: {n_decrease}/{n_total}")
    
    diffs = [r['difference'] for r in ssw_results]
    print(f"Mean danger difference: {np.mean(diffs):.4f}")
    
    # Sign test
    from scipy import stats
    if n_total >= 3:
        try:
            p_sign = stats.binomtest(n_decrease, n_total, 0.5).pvalue
        except:
            p_sign = float('nan')
        print(f"Sign test: P={p_sign:.4f}")
    
    # Effect by phase
    print("\n=== Phase-resolved analysis ===")
    for phase_name, d_start, d_end in [
        ("Pre-SSW (d-15 to d-1)", -15, -1),
        ("Post-SSW (d0 to d+15)", 0, 15),
        ("Late post (d+16 to d+30)", 16, 30),
    ]:
        phase_ssw = []
        phase_ctrl = []
        for ssw_str in SSW_EVENTS:
            ssw_dt = datetime.strptime(ssw_str, "%Y-%m-%d")
            for d in range(d_start, d_end+1):
                dt = ssw_dt + timedelta(days=d)
                key = dt.strftime("%Y-%m-%d")
                if key in daily_mean:
                    phase_ssw.append(daily_mean[key])
                # Control
                for offset in [-2, -1, 1, 2]:
                    ctrl_dt = dt.replace(year=dt.year + offset)
                    ctrl_key = ctrl_dt.strftime("%Y-%m-%d")
                    if ctrl_key in daily_mean:
                        phase_ctrl.append(daily_mean[ctrl_key])
        
        if phase_ssw and phase_ctrl:
            ssw_m = np.mean(phase_ssw)
            ctrl_m = np.mean(phase_ctrl)
            try:
                u_stat, u_p = stats.mannwhitneyu(phase_ssw, phase_ctrl, alternative='two-sided')
            except:
                u_p = float('nan')
            print(f"  {phase_name}: SSW={ssw_m:.3f} vs Ctrl={ctrl_m:.3f} (diff={ssw_m-ctrl_m:+.3f}, MW P={u_p:.4f})")

    # Region-level analysis
    print("\n=== Region-level analysis ===")
    for rid, rname in REGIONS.items():
        region_recs = [r for r in all_records if r['region_id'] == rid]
        if not region_recs:
            continue
        
        reg_daily = defaultdict(float)
        for r in region_recs:
            reg_daily[r['date']] = r['danger_level']
        
        reg_decrease = 0
        reg_total = 0
        for ssw_str in SSW_EVENTS:
            ssw_dt = datetime.strptime(ssw_str, "%Y-%m-%d")
            ssw_vals = []
            ctrl_vals = []
            for d in range(-15, 16):
                dt = ssw_dt + timedelta(days=d)
                key = dt.strftime("%Y-%m-%d")
                if key in reg_daily:
                    ssw_vals.append(reg_daily[key])
                for offset in [-2, -1, 1, 2]:
                    ctrl_dt = dt.replace(year=dt.year + offset)
                    ctrl_key = ctrl_dt.strftime("%Y-%m-%d")
                    if ctrl_key in reg_daily:
                        ctrl_vals.append(reg_daily[ctrl_key])
            
            if ssw_vals and ctrl_vals:
                if np.mean(ssw_vals) < np.mean(ctrl_vals):
                    reg_decrease += 1
                reg_total += 1
        
        if reg_total > 0:
            print(f"  {rname}: {reg_decrease}/{reg_total} SSW events show decreased danger")

print("\nDone!")
