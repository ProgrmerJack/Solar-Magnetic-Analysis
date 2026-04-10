"""
R34 Multi-Country Replication Analysis
Parse EAWS 142K-row dataset properly and expand analysis:
1. Parse semicolon-delimited data
2. SSW events in the 2011-2015 window
3. Per-country danger level analysis
4. Geographic gradient with proper statistics
5. Try to fetch French BRA data from API
"""
import pandas as pd
import numpy as np
from scipy import stats
from collections import Counter
import json, os

# ============================================================
# 1. PARSE EAWS 5-COUNTRY DATA
# ============================================================
print("="*60)
print("1. PARSING EAWS MULTI-COUNTRY DATA")
print("="*60)

df = pd.read_csv('data/cryosphere/european_alps/data_dmax.csv', sep=';')
print(f"Columns: {list(df.columns)}")
print(f"Shape: {df.shape}")
print(f"Date range: {df['date'].min()} to {df['date'].max()}")
print(f"\nCountries: {df['country'].unique()}")
print(f"Forecast centers: {df['forecastCenter'].unique()}")
df = df.dropna(subset=['country'])
print(f"\nRegions per country:")
for country in sorted(df['country'].unique()):
    n_reg = df[df['country']==country]['warningRegion'].nunique()
    n_rows = len(df[df['country']==country])
    print(f"  {country}: {n_reg} regions, {n_rows} rows")

df['date'] = pd.to_datetime(df['date'])
df['danger'] = pd.to_numeric(df['dangerLevelMax'], errors='coerce')

# ============================================================
# 2. SSW EVENTS IN COVERAGE (2011-12 to 2015-04)
# ============================================================
print("\n" + "="*60)
print("2. SSW EVENTS IN EAWS COVERAGE")
print("="*60)

# SSW events from catalog
ssw_in_range = [
    '2012-01-06',   # January 2012 SSW
    '2013-01-06',   # January 2013 SSW
]

# Also check if there are SSWs in early 2012 or late 2014
# Butler catalog SSW dates near this period
additional_ssw = [
    '2014-01-06',   # Possible - check if in catalog
]

print(f"SSW events in range: {ssw_in_range}")

results = {}

for ssw_str in ssw_in_range:
    ssw_date = pd.Timestamp(ssw_str)
    print(f"\n--- SSW Event: {ssw_str} ---")
    
    # SSW window: ±15 days
    ssw_mask = (df['date'] >= ssw_date - pd.Timedelta(days=15)) & \
               (df['date'] <= ssw_date + pd.Timedelta(days=15))
    
    # DOY-matched control
    doy_center = ssw_date.dayofyear
    ctrl_mask = (df['date'].dt.dayofyear >= doy_center - 15) & \
                (df['date'].dt.dayofyear <= doy_center + 15) & ~ssw_mask
    
    ssw_df = df[ssw_mask].dropna(subset=['danger'])
    ctrl_df = df[ctrl_mask].dropna(subset=['danger'])
    
    # Overall
    ssw_mean = ssw_df['danger'].mean()
    ctrl_mean = ctrl_df['danger'].mean()
    stat, p = stats.mannwhitneyu(ssw_df['danger'], ctrl_df['danger'], alternative='two-sided')
    d = (ssw_mean - ctrl_mean) / np.sqrt((ssw_df['danger'].std()**2 + ctrl_df['danger'].std()**2)/2)
    
    print(f"  Overall: SSW={ssw_mean:.3f} (n={len(ssw_df)}), Ctrl={ctrl_mean:.3f} (n={len(ctrl_df)})")
    print(f"  Diff={ssw_mean-ctrl_mean:+.3f}, MW P={p:.4f}, d={d:.3f}")
    
    # Per country
    event_results = {'overall': {'ssw_mean': round(ssw_mean, 3), 'ctrl_mean': round(ctrl_mean, 3),
                                  'diff': round(ssw_mean-ctrl_mean, 3), 'p': round(p, 4)}}
    
    print(f"\n  Per country:")
    for country in sorted(df['country'].unique()):
        ssw_c = ssw_df[ssw_df['country']==country]['danger']
        ctrl_c = ctrl_df[ctrl_df['country']==country]['danger']
        if len(ssw_c) > 10 and len(ctrl_c) > 10:
            stat_c, p_c = stats.mannwhitneyu(ssw_c, ctrl_c, alternative='two-sided')
            d_c = (ssw_c.mean() - ctrl_c.mean()) / np.sqrt((ssw_c.std()**2 + ctrl_c.std()**2)/2)
            direction = 'DECREASE' if ssw_c.mean() < ctrl_c.mean() else 'INCREASE'
            print(f"    {country}: SSW={ssw_c.mean():.2f}, Ctrl={ctrl_c.mean():.2f}, diff={ssw_c.mean()-ctrl_c.mean():+.2f}, P={p_c:.4f} [{direction}]")
            event_results[country] = {
                'ssw_mean': round(ssw_c.mean(), 3),
                'ctrl_mean': round(ctrl_c.mean(), 3),
                'diff': round(ssw_c.mean()-ctrl_c.mean(), 3),
                'p': round(p_c, 4),
                'direction': direction,
                'n_ssw': len(ssw_c),
                'n_ctrl': len(ctrl_c)
            }
    
    results[ssw_str] = event_results

# ============================================================
# 3. GEOGRAPHIC GRADIENT BY FORECAST CENTER
# ============================================================
print("\n" + "="*60)
print("3. GEOGRAPHIC GRADIENT BY FORECAST CENTER")
print("="*60)

# Combine both SSW events
all_ssw_mask = pd.Series(False, index=df.index)
all_ctrl_mask = pd.Series(False, index=df.index)
for ssw_str in ssw_in_range:
    ssw_date = pd.Timestamp(ssw_str)
    ssw_m = (df['date'] >= ssw_date - pd.Timedelta(days=15)) & \
            (df['date'] <= ssw_date + pd.Timedelta(days=15))
    doy_center = ssw_date.dayofyear
    ctrl_m = (df['date'].dt.dayofyear >= doy_center - 15) & \
             (df['date'].dt.dayofyear <= doy_center + 15) & ~ssw_m
    all_ssw_mask |= ssw_m
    all_ctrl_mask |= ctrl_m

# By forecast center
gradient_results = {}
print(f"\n{'Center':>10} | {'Country':>3} | {'SSW':>6} | {'Ctrl':>6} | {'Diff':>6} | {'P':>8} | {'Dir':>10}")
print("-" * 70)

for center in sorted(df['forecastCenter'].unique()):
    ssw_c = df[all_ssw_mask & (df['forecastCenter']==center)]['danger'].dropna()
    ctrl_c = df[all_ctrl_mask & (df['forecastCenter']==center)]['danger'].dropna()
    country = df[df['forecastCenter']==center]['country'].iloc[0]
    
    if len(ssw_c) > 20 and len(ctrl_c) > 20:
        stat, p = stats.mannwhitneyu(ssw_c, ctrl_c, alternative='two-sided')
        diff = ssw_c.mean() - ctrl_c.mean()
        direction = 'DECREASE' if diff < 0 else 'INCREASE'
        print(f"{center:>10} | {country:>3} | {ssw_c.mean():6.2f} | {ctrl_c.mean():6.2f} | {diff:+6.2f} | {p:8.4f} | {direction:>10}")
        gradient_results[center] = {
            'country': country,
            'ssw_mean': round(ssw_c.mean(), 3),
            'ctrl_mean': round(ctrl_c.mean(), 3),
            'diff': round(diff, 3),
            'p': round(p, 4),
            'direction': direction
        }

# Count directions
decreases = sum(1 for v in gradient_results.values() if v['direction'] == 'DECREASE')
increases = sum(1 for v in gradient_results.values() if v['direction'] == 'INCREASE')
total = len(gradient_results)
print(f"\nGradient summary: {decreases}/{total} decrease, {increases}/{total} increase")

# Geographic pattern: Western Alps vs Northern Alps
western = ['VDA', 'PIE', 'LOM']  # Italy west
central = ['TAA', 'TRE', 'BZL', 'VEN']  # Italy central/east
swiss = ['VS', 'GR', 'BE', 'SG']  # Swiss cantons (if present)
austrian = ['TIR', 'SBG', 'VBG', 'KTN', 'STM']  # Austria

for group_name, centers in [('Western IT', western), ('Central/East IT', central), ('Austrian', austrian)]:
    ssw_vals = []
    ctrl_vals = []
    for c in centers:
        if c in gradient_results:
            ssw_v = df[all_ssw_mask & (df['forecastCenter']==c)]['danger'].dropna().values
            ctrl_v = df[all_ctrl_mask & (df['forecastCenter']==c)]['danger'].dropna().values
            ssw_vals.extend(ssw_v)
            ctrl_vals.extend(ctrl_v)
    
    if ssw_vals and ctrl_vals:
        ssw_arr = np.array(ssw_vals)
        ctrl_arr = np.array(ctrl_vals)
        stat, p = stats.mannwhitneyu(ssw_arr, ctrl_arr, alternative='two-sided')
        diff = ssw_arr.mean() - ctrl_arr.mean()
        print(f"\n{group_name}: SSW={ssw_arr.mean():.3f}, Ctrl={ctrl_arr.mean():.3f}, diff={diff:+.3f}, P={p:.4f}")

results['gradient'] = gradient_results

# ============================================================
# 4. COMBINED CONCORDANCE ACROSS ALL COUNTRIES
# ============================================================
print("\n" + "="*60)
print("4. COMBINED CONCORDANCE")
print("="*60)

# Count event-country pairs with decrease
concordance_pairs = []
for ssw_str in ssw_in_range:
    ssw_date = pd.Timestamp(ssw_str)
    ssw_mask = (df['date'] >= ssw_date - pd.Timedelta(days=15)) & \
               (df['date'] <= ssw_date + pd.Timedelta(days=15))
    doy_center = ssw_date.dayofyear
    ctrl_mask = (df['date'].dt.dayofyear >= doy_center - 15) & \
                (df['date'].dt.dayofyear <= doy_center + 15) & ~ssw_mask
    
    for country in sorted(df['country'].unique()):
        ssw_c = df[ssw_mask & (df['country']==country)]['danger'].dropna()
        ctrl_c = df[ctrl_mask & (df['country']==country)]['danger'].dropna()
        if len(ssw_c) > 10:
            decrease = ssw_c.mean() < ctrl_c.mean()
            concordance_pairs.append({
                'ssw': ssw_str,
                'country': country,
                'ssw_mean': ssw_c.mean(),
                'ctrl_mean': ctrl_c.mean(),
                'decrease': decrease
            })

n_decrease = sum(p['decrease'] for p in concordance_pairs)
n_total = len(concordance_pairs)
print(f"Event-country pairs: {n_total}")
print(f"Showing decrease: {n_decrease}/{n_total} ({n_decrease/n_total:.1%})")
binom_p = stats.binomtest(n_decrease, n_total, 0.5, alternative='greater').pvalue
print(f"Binomial test P = {binom_p:.4f}")

for p in concordance_pairs:
    dir_str = 'DOWN' if p['decrease'] else 'UP'
    print(f"  {p['ssw']} {p['country']}: {p['ssw_mean']:.2f} vs {p['ctrl_mean']:.2f} [{dir_str}]")

results['concordance'] = {
    'n_decrease': n_decrease,
    'n_total': n_total,
    'pct_decrease': round(n_decrease/n_total*100, 1),
    'binom_p': round(binom_p, 4),
    'pairs': concordance_pairs
}

# ============================================================
# 5. TRY FRENCH BRA DATA VIA METEOFRANCE API
# ============================================================
print("\n" + "="*60)
print("5. ATTEMPTING FRENCH BRA DATA DOWNLOAD")
print("="*60)

import requests

# Try the Météo-France public BRA API
urls_to_try = [
    'https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA/',
    'https://public-api.meteofrance.fr/public/DPBRA/v1/massif/liste',
    'https://rpcache-aa.meteofrance.com/internet2018client/2.0/report?domain=BRA&report_type=BULLETINBRA&report_subtype=BRA&massif=MONT-BLANC',
]

for url in urls_to_try:
    try:
        r = requests.get(url, timeout=15, headers={'User-Agent': 'AcademicResearch/1.0'})
        print(f"  {url[:60]}... : status={r.status_code}")
        if r.status_code == 200 and len(r.content) > 100:
            print(f"    Content type: {r.headers.get('content-type', 'unknown')}")
            if 'json' in r.headers.get('content-type', ''):
                print(f"    JSON preview: {r.text[:200]}")
            elif 'xml' in r.headers.get('content-type', ''):
                print(f"    XML preview: {r.text[:200]}")
            else:
                print(f"    Content preview: {r.text[:200]}")
    except Exception as e:
        print(f"  {url[:60]}... : ERROR: {e}")

# Try data.gouv.fr for French avalanche data
try:
    r = requests.get('https://www.data.gouv.fr/api/1/datasets/?q=avalanche&page_size=5', timeout=15)
    if r.status_code == 200:
        datasets = r.json().get('data', [])
        print(f"\ndata.gouv.fr avalanche datasets: {len(datasets)}")
        for ds in datasets:
            title = ds.get('title', 'Unknown')
            n_resources = len(ds.get('resources', []))
            print(f"  - {title} ({n_resources} resources)")
            for res in ds.get('resources', [])[:3]:
                fmt = res.get('format', 'unknown')
                url = res.get('url', '')
                print(f"    [{fmt}] {url[:80]}")
except Exception as e:
    print(f"data.gouv.fr error: {e}")

# Save results
os.makedirs('data/results', exist_ok=True)
with open('data/results/r34_multicountry.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to data/results/r34_multicountry.json")
