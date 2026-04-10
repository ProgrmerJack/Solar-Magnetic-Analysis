"""
Download US avalanche danger ratings from avalanche.org API
Covers ALL US avalanche centers: CAIC, NWAC, SAC, UAC, etc.
"""
import urllib.request, json, time, os
import pandas as pd

OUT_DIR = 'data/cryosphere/us_danger_ratings'
os.makedirs(OUT_DIR, exist_ok=True)

all_forecasts = []
centers_seen = set()

for year in range(2011, 2026):
    for month in [11, 12, 1, 2, 3, 4]:
        y = year if month >= 11 else year + 1
        if y > 2025: continue
        
        start = f"{y}-{month:02d}-01"
        if month in [1, 3]: end = f"{y}-{month:02d}-31"
        elif month == 12: end = f"{y}-12-31"
        elif month in [4, 11]: end = f"{y}-{month:02d}-30"
        elif month == 2: end = f"{y}-02-28"
        
        url = f'https://api.avalanche.org/v2/public/products?type=forecast&date_start={start}&date_end={end}'
        
        try:
            req = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (Academic Research)',
                'Accept': 'application/json'
            })
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read())
                
            if isinstance(data, list):
                for item in data:
                    center = item.get('avalanche_center', {}).get('name', 'Unknown')
                    centers_seen.add(center)
                    all_forecasts.append({
                        'date': item.get('start_date', '')[:10],
                        'center': center,
                        'danger_rating': item.get('danger_rating'),
                    })
                print(f"  {start}: {len(data)} forecasts")
        except Exception as e:
            print(f"  {start}: Error - {e}")
        
        time.sleep(0.5)

print(f"\nTotal: {len(all_forecasts)}")
print(f"Centers: {sorted(centers_seen)}")

df = pd.DataFrame(all_forecasts)
df['date'] = pd.to_datetime(df['date'])
df.to_csv(f'{OUT_DIR}/us_danger_ratings_all.csv', index=False)
print(f"Saved {len(df)} records")
print(df.groupby('center').size().sort_values(ascending=False).to_string())
