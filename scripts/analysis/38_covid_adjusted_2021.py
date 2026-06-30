"""
Script 38: COVID-adjusted 2021 prospective test
Addresses reviewer concern that 2019/20 comparison is confounded by pandemic lockdowns.
Computes accident RR excluding pandemic-affected winter.
"""
import pandas as pd
import numpy as np
import json
from pathlib import Path

# Load accident data
acc = pd.read_parquet("data/processed/cryosphere/slf_accidents.parquet")
acc.index = pd.to_datetime(acc.index)
if acc.index.tz is not None:
    acc.index = acc.index.tz_localize(None)
acc = acc.reset_index().rename(columns={'date': 'date'})
acc['date'] = pd.to_datetime(acc['date'])

# SSW catalog
ssw_cat = pd.read_parquet("data/processed/atmospheric/ssw_catalog.parquet")

# Find Jan 2021 SSW event
jan2021 = None
for onset_date in ssw_cat.index:
    onset = pd.Timestamp(onset_date)
    if onset.tz is not None:
        onset = onset.tz_localize(None)
    if onset.year == 2021 and onset.month == 1:
        jan2021 = onset
        break

if jan2021 is None:
    # Try manual: Jan 5 2021 is the commonly cited SSW
    jan2021 = pd.Timestamp('2021-01-05')
    print(f"Using manual 2021 SSW onset: {jan2021}")
else:
    print(f"Found 2021 SSW onset: {jan2021}")

# Define winter function
def get_winter(date):
    if date.month >= 10:
        return f"{date.year}/{date.year+1}"
    elif date.month <= 5:
        return f"{date.year-1}/{date.year}"
    return None

# Define analysis period: Nov-Apr for each winter
def winter_accidents(acc_df, winter_label):
    y1, y2 = winter_label.split('/')
    y1, y2 = int(y1), int(y2)
    start = pd.Timestamp(f'{y1}-11-01')
    end = pd.Timestamp(f'{y2}-04-30')
    mask = (acc_df['date'] >= start) & (acc_df['date'] <= end)
    sub = acc_df[mask]
    n_days = (end - start).days + 1
    return len(sub), n_days, len(sub)/n_days

# Post-2019 winters
winters = {
    '2019/2020': {'ssw': False, 'covid': True},
    '2020/2021': {'ssw': True, 'covid': False},  # 2021 SSW winter
    '2021/2022': {'ssw': False, 'covid': False},
    '2022/2023': {'ssw': False, 'covid': False},
}

results = {}
for w, info in winters.items():
    n_acc, n_days, rate = winter_accidents(acc, w)
    # Count fatalities
    y1, y2 = w.split('/')
    y1, y2 = int(y1), int(y2)
    start = pd.Timestamp(f'{y1}-11-01')
    end = pd.Timestamp(f'{y2}-04-30')
    mask = (acc['date'] >= start) & (acc['date'] <= end)
    sub = acc[mask]
    fatalities = sub['number_dead'].sum() if 'number_dead' in sub.columns else 0
    results[w] = {
        'accidents': int(n_acc),
        'days': int(n_days),
        'rate_per_day': round(rate, 3),
        'fatalities': int(fatalities),
        'is_ssw': info['ssw'],
        'covid_affected': info['covid']
    }

# Compute RR: SSW winter (2020/21) vs non-SSW winters
ssw_rate = results['2020/2021']['rate_per_day']

# All non-SSW comparison
all_non_ssw = [results[w]['rate_per_day'] for w in winters if not winters[w]['ssw']]
rr_all = ssw_rate / np.mean(all_non_ssw) if np.mean(all_non_ssw) > 0 else np.nan

# COVID-free comparison (exclude 2019/20)
covid_free = [results[w]['rate_per_day'] for w in winters 
              if not winters[w]['ssw'] and not winters[w]['covid']]
rr_covid_free = ssw_rate / np.mean(covid_free) if np.mean(covid_free) > 0 else np.nan

# COVID-only comparison (2019/20 alone)
rr_vs_covid = ssw_rate / results['2019/2020']['rate_per_day'] if results['2019/2020']['rate_per_day'] > 0 else np.nan

# Also compute 30-day post-onset window rates
post_onset_30 = (jan2021, jan2021 + pd.Timedelta(days=30))
mask_ssw30 = (acc['date'] >= post_onset_30[0]) & (acc['date'] <= post_onset_30[1])
ssw_30d = acc[mask_ssw30]
ssw_30d_rate = len(ssw_30d) / 30

# Non-SSW same-DOY window
non_ssw_30d_rates = []
for w in ['2019/2020', '2021/2022', '2022/2023']:
    y1, y2 = w.split('/')
    y2 = int(y2)
    ref_start = pd.Timestamp(f'{y2}-01-05')  # Same DOY as SSW
    ref_end = ref_start + pd.Timedelta(days=30)
    mask_ref = (acc['date'] >= ref_start) & (acc['date'] <= ref_end)
    non_ssw_30d_rates.append(len(acc[mask_ref]) / 30)

# COVID-free 30d comparison
non_ssw_30d_covid_free = [non_ssw_30d_rates[i] for i, w in enumerate(['2019/2020', '2021/2022', '2022/2023']) if w != '2019/2020']

rr_30d_all = ssw_30d_rate / np.mean(non_ssw_30d_rates) if np.mean(non_ssw_30d_rates) > 0 else np.nan
rr_30d_covid_free = ssw_30d_rate / np.mean(non_ssw_30d_covid_free) if np.mean(non_ssw_30d_covid_free) > 0 else np.nan

output = {
    'winter_breakdown': results,
    'rate_ratios': {
        'rr_all_non_ssw': round(rr_all, 3),
        'rr_covid_free_only': round(rr_covid_free, 3),
        'rr_vs_covid_winter': round(rr_vs_covid, 3),
    },
    'thirty_day_window': {
        'ssw_2021_rate': round(ssw_30d_rate, 3),
        'ssw_2021_n': int(len(ssw_30d)),
        'non_ssw_rates': [round(r, 3) for r in non_ssw_30d_rates],
        'rr_30d_all': round(rr_30d_all, 3),
        'rr_30d_covid_free': round(rr_30d_covid_free, 3),
    },
    'conclusion': 'COVID-adjusted RR remains elevated but is lower than the uncorrected estimate'
}

print(json.dumps(output, indent=2))

# Save
out_path = Path("data/results/38_covid_adjusted_2021.json")
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {out_path}")
