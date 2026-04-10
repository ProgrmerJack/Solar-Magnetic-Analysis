"""
Download ALBINA (Austria/South Tyrol/Trentino) danger ratings for SSW analysis.
Public data from avalanche.report/albina_files/ (available 2018-12-04 onwards)
"""
import urllib.request, json, os, time, sys
from datetime import datetime, timedelta

os.makedirs('data/cryosphere/albina', exist_ok=True)

danger_map = {
    'low': 1, 'moderate': 2, 'considerable': 3, 'high': 4, 'very_high': 5,
    'no_rating': 0, 'no_snow': 0, 'missing': 0
}

# SSW events in ALBINA range (2018-12-04+)
# 2019-01-02, 2021-01-05, 2023-02-16
ssw_events = ['2019-01-02', '2021-01-05', '2023-02-16']

# Generate date ranges: 45 days before to 45 days after each SSW
all_dates = set()
for ssw in ssw_events:
    ssw_dt = datetime.strptime(ssw, '%Y-%m-%d')
    for delta in range(-45, 46):
        d = ssw_dt + timedelta(days=delta)
        if d >= datetime(2018, 12, 4):
            all_dates.add(d.strftime('%Y-%m-%d'))

# Also add control periods (non-SSW winters)
control_periods = [
    ('2019-12-01', '2020-03-31'),  # 2019-2020 winter (no SSW)
    ('2021-12-01', '2022-03-31'),  # 2021-2022 winter (no SSW)
]
for start_str, end_str in control_periods:
    start = datetime.strptime(start_str, '%Y-%m-%d')
    end = datetime.strptime(end_str, '%Y-%m-%d')
    d = start
    while d <= end:
        all_dates.add(d.strftime('%Y-%m-%d'))
        d += timedelta(days=1)

all_dates = sorted(all_dates)
print(f"Total dates to download: {len(all_dates)}")

all_records = []
success = 0
fail = 0

for date_str in all_dates:
    url = f'https://avalanche.report/albina_files/{date_str}/avalanche_report.json'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 Research'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.load(resp)
        
        for item in data:
            regions = item.get('regions', [])
            forenoon = item.get('forenoon', {})
            
            dr_above = forenoon.get('dangerRatingAbove', 'missing')
            dr_below = forenoon.get('dangerRatingBelow', 'missing')
            
            dl_above = danger_map.get(dr_above, 0)
            dl_below = danger_map.get(dr_below, 0)
            dl_max = max(dl_above, dl_below)
            
            # Identify country from region codes
            for region in regions:
                country = 'AT' if region.startswith('AT') else ('IT' if region.startswith('IT') else 'unknown')
                all_records.append({
                    'date': date_str,
                    'region': region,
                    'country': country,
                    'danger_above': dl_above,
                    'danger_below': dl_below,
                    'danger_max': dl_max,
                })
        
        success += 1
        if success % 50 == 0:
            sys.stdout.write(f"  Downloaded {success} days ({date_str})...\n")
            sys.stdout.flush()
    except Exception as e:
        fail += 1
    
    time.sleep(0.1)  # Rate limit

print(f"\nDone: {success} success, {fail} fail, {len(all_records)} total records")

# Save
with open('data/cryosphere/albina/albina_danger.json', 'w') as f:
    json.dump(all_records, f)

# Analyze
import pandas as pd
import numpy as np
from scipy import stats

df = pd.DataFrame(all_records)
df['date'] = pd.to_datetime(df['date'])
df = df[df['danger_max'] > 0]

print(f"\nValid records: {len(df)}")
print(f"Countries: {df['country'].value_counts().to_dict()}")
print(f"Date range: {df['date'].min().date()} to {df['date'].max().date()}")

# SSW analysis
for ssw_date_str in ssw_events:
    ssw_dt = pd.Timestamp(ssw_date_str)
    
    # SSW window: ±15 days
    ssw_win = df[(df['date'] >= ssw_dt - pd.Timedelta(days=15)) & 
                 (df['date'] <= ssw_dt + pd.Timedelta(days=15))]
    
    # DOY-matched control from non-SSW winters
    ssw_doys = ssw_win['date'].dt.dayofyear.unique()
    ctrl = df[df['date'].dt.dayofyear.isin(ssw_doys)]
    # Remove SSW winter from control
    ssw_year = ssw_dt.year
    ssw_winter_start = pd.Timestamp(f'{ssw_year-1}-11-01') if ssw_dt.month <= 6 else pd.Timestamp(f'{ssw_year}-11-01')
    ssw_winter_end = pd.Timestamp(f'{ssw_year}-06-30') if ssw_dt.month <= 6 else pd.Timestamp(f'{ssw_year+1}-06-30')
    ctrl = ctrl[(ctrl['date'] < ssw_winter_start) | (ctrl['date'] > ssw_winter_end)]
    
    if len(ssw_win) > 0 and len(ctrl) > 0:
        ssw_mean = ssw_win['danger_max'].mean()
        ctrl_mean = ctrl['danger_max'].mean()
        u, p = stats.mannwhitneyu(ssw_win['danger_max'], ctrl['danger_max'], alternative='less')
        print(f"\nSSW {ssw_date_str}: SSW={ssw_mean:.3f} (n={len(ssw_win)}), Ctrl={ctrl_mean:.3f} (n={len(ctrl)}), P={p:.4f}")
        
        # By country
        for country in ['AT', 'IT']:
            c_ssw = ssw_win[ssw_win['country'] == country]
            c_ctrl = ctrl[ctrl['country'] == country]
            if len(c_ssw) > 0 and len(c_ctrl) > 0:
                c_mean_s = c_ssw['danger_max'].mean()
                c_mean_c = c_ctrl['danger_max'].mean()
                _, c_p = stats.mannwhitneyu(c_ssw['danger_max'], c_ctrl['danger_max'], alternative='less')
                print(f"  {country}: SSW={c_mean_s:.3f}, Ctrl={c_mean_c:.3f}, P={c_p:.4f}")

# Overall SSW vs control
print("\n--- OVERALL ---")
ssw_all = pd.DataFrame()
ctrl_all = pd.DataFrame()
for ssw_date_str in ssw_events:
    ssw_dt = pd.Timestamp(ssw_date_str)
    win = df[(df['date'] >= ssw_dt - pd.Timedelta(days=15)) & 
             (df['date'] <= ssw_dt + pd.Timedelta(days=15))]
    ssw_all = pd.concat([ssw_all, win])

# Control = non-SSW winter data
for start_str, end_str in control_periods:
    ctrl_period = df[(df['date'] >= start_str) & (df['date'] <= end_str)]
    ctrl_all = pd.concat([ctrl_all, ctrl_period])

if len(ssw_all) > 0 and len(ctrl_all) > 0:
    print(f"SSW periods: {ssw_all['danger_max'].mean():.3f} (n={len(ssw_all)})")
    print(f"Control: {ctrl_all['danger_max'].mean():.3f} (n={len(ctrl_all)})")
    u, p = stats.mannwhitneyu(ssw_all['danger_max'], ctrl_all['danger_max'], alternative='less')
    print(f"Mann-Whitney P: {p:.6f}")
    
    pooled_std = np.sqrt((ssw_all['danger_max'].var() + ctrl_all['danger_max'].var()) / 2)
    d = (ssw_all['danger_max'].mean() - ctrl_all['danger_max'].mean()) / pooled_std
    print(f"Cohen's d: {d:.3f}")
