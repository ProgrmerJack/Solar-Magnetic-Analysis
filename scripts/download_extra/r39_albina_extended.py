"""Download ALBINA CAAMLv6 data (2022-2024) and merge with existing records."""
import json, urllib.request, time, datetime, os

OUT = "data/cryosphere/albina/albina_danger.json"
with open(OUT) as f:
    existing = json.load(f)
existing_dates = set(r['date'] for r in existing)
print(f"Existing: {len(existing)} records, {len(existing_dates)} dates")

danger_map = {'low':1,'moderate':2,'considerable':3,'high':4,'very_high':5}
new_records = []

# Download CAAMLv6 files for winter 2022-23 and 2023-24
start = datetime.date(2022, 11, 1)
end = datetime.date(2024, 4, 30)
current = start
success = fail = skip = 0

while current <= end:
    ds = current.isoformat()
    if ds in existing_dates:
        skip += 1
        current += datetime.timedelta(days=1)
        continue
    
    # Try CAAMLv6 JSON for AT-07 (English)
    url = f"https://avalanche.report/albina_files/{ds}/{ds}_AT-07_en_CAAMLv6.json"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
        
        bulletins = data.get('bulletins', data) if isinstance(data, dict) else data
        if isinstance(bulletins, dict):
            bulletins = bulletins.get('bulletins', [bulletins])
        
        for b in bulletins:
            # Get danger ratings
            ratings = b.get('dangerRatings', [])
            max_danger = 0
            above = below = 0
            for r in ratings:
                val = danger_map.get(r.get('mainValue', ''), 0)
                if val > max_danger:
                    max_danger = val
                elev = r.get('elevation', {})
                if 'lowerBound' in elev:
                    above = max(above, val)
                elif 'upperBound' in elev:
                    below = max(below, val)
                else:
                    above = below = max(above, val)
            
            if max_danger == 0:
                continue
            
            regions = b.get('regions', [])
            for reg in regions:
                rid = reg.get('regionID', '') if isinstance(reg, dict) else str(reg)
                country = 'AT' if rid.startswith('AT') else ('IT' if rid.startswith('IT') else 'other')
                if country in ('AT', 'IT'):
                    new_records.append({
                        'date': ds, 'region': rid, 'country': country,
                        'danger_above': above if above else max_danger,
                        'danger_below': below if below else max_danger,
                        'danger_max': max_danger
                    })
        success += 1
    except Exception as e:
        fail += 1
    
    time.sleep(0.12)
    current += datetime.timedelta(days=1)
    elapsed = (current - start).days
    if elapsed % 30 == 0:
        print(f"  Progress: {current} | new={len(new_records)} | ok={success} fail={fail} skip={skip}")

print(f"\nDone: {len(new_records)} new records | ok={success} fail={fail} skip={skip}")

# Merge
all_records = existing + new_records
with open(OUT, 'w') as f:
    json.dump(all_records, f)
print(f"Total saved: {len(all_records)} records")

# Summary
all_dates = sorted(set(r['date'] for r in all_records))
print(f"Date range: {all_dates[0]} to {all_dates[-1]}")
print(f"Unique dates: {len(all_dates)}")
