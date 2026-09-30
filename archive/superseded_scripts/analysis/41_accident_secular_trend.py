"""
Script 41: Accident secular trend analysis
Addresses reviewer concern about 55-year secular trends in backcountry recreation,
equipment (airbags post-2005), education, and rescue capability.
Checks whether the SSW accident signal survives trend adjustment.
"""
import pandas as pd
import numpy as np
from scipy import stats
import json
from pathlib import Path

# Load data
acc = pd.read_parquet("data/processed/cryosphere/slf_accidents.parquet")
acc.index = pd.to_datetime(acc.index)
if acc.index.tz is not None:
    acc.index = acc.index.tz_localize(None)
acc = acc.reset_index().rename(columns={acc.index.name or 'index': 'date'})
acc['date'] = pd.to_datetime(acc['date'])

ssw_cat = pd.read_parquet("data/processed/atmospheric/ssw_catalog.parquet")

# Define winter season
def get_winter_year(date):
    """Return the starting year of the winter season."""
    if date.month >= 10:
        return date.year
    elif date.month <= 5:
        return date.year - 1
    return None

acc['winter_year'] = acc['date'].apply(get_winter_year)
acc = acc.dropna(subset=['winter_year'])
acc['winter_year'] = acc['winter_year'].astype(int)

# Compute winter-level statistics
winters = acc.groupby('winter_year').agg(
    n_accidents=('date', 'count'),
    n_fatalities=('number_dead', 'sum') if 'number_dead' in acc.columns else ('date', 'count'),
).reset_index()

# Add season length (Nov-Apr = ~181 days)
winters['rate_per_day'] = winters['n_accidents'] / 181

# Fit secular trend
valid_winters = winters[(winters['winter_year'] >= 1970) & (winters['winter_year'] <= 2023)]
years = valid_winters['winter_year'].values
rates = valid_winters['rate_per_day'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(years, rates)
print(f"Secular trend: slope={slope:.4f} acc/day/year, R²={r_value**2:.3f}, P={p_value:.4f}")

# Also fit log-linear trend
log_rates = np.log(rates + 0.01)
slope_log, intercept_log, r_log, p_log, _ = stats.linregress(years, log_rates)
print(f"Log-linear trend: slope={slope_log:.4f}/year, R²={r_log**2:.3f}, P={p_log:.4f}")

# Now compute SSW effect with and without trend adjustment
# Get SSW events
ssw_events = []
for onset_date in ssw_cat.index:
    onset = pd.Timestamp(onset_date)
    if onset.tz is not None:
        onset = onset.tz_localize(None)
    ssw_events.append(onset)

# For each SSW event, compute accident RR
def compute_accident_rr(acc_data, onset, window=30):
    """Compute accident rate ratio for one SSW event."""
    onset_ts = pd.Timestamp(onset)
    ssw_start = onset_ts - pd.Timedelta(days=window//2)
    ssw_end = onset_ts + pd.Timedelta(days=window//2)
    
    # SSW window accidents
    ssw_mask = (acc_data['date'] >= ssw_start) & (acc_data['date'] <= ssw_end)
    ssw_count = ssw_mask.sum()
    ssw_rate = ssw_count / window
    
    # DOY-matched control
    ssw_year = onset_ts.year
    ssw_doys = pd.date_range(ssw_start, ssw_end).dayofyear
    
    control_counts = []
    for yr in range(max(acc_data['date'].dt.year.min(), 1970), 
                     min(acc_data['date'].dt.year.max(), 2024) + 1):
        if yr == ssw_year:
            continue
        yr_mask = acc_data['date'].dt.year == yr
        doy_mask = acc_data['date'].dt.dayofyear.isin(ssw_doys)
        count = (yr_mask & doy_mask).sum()
        control_counts.append(count / window)
    
    if len(control_counts) == 0 or np.mean(control_counts) == 0:
        return np.nan, ssw_rate, 0
    
    rr = ssw_rate / np.mean(control_counts)
    return rr, ssw_rate, np.mean(control_counts)

# Compute RR for each event
event_results = []
for onset in ssw_events:
    rr, ssw_rate, ctrl_rate = compute_accident_rr(acc, onset)
    year = onset.year
    
    # Trend-expected rate for this year
    trend_expected = slope * year + intercept
    
    # Trend-adjusted: divide observed rate by expected trend
    if trend_expected > 0:
        ssw_adj = ssw_rate / trend_expected
        ctrl_adj = ctrl_rate / trend_expected
        rr_adj = ssw_adj / ctrl_adj if ctrl_adj > 0 else np.nan
    else:
        rr_adj = np.nan
    
    event_results.append({
        'onset': str(onset.date()),
        'year': year,
        'rr_raw': round(rr, 3) if not np.isnan(rr) else None,
        'rr_trend_adjusted': round(rr_adj, 3) if not np.isnan(rr_adj) else None,
        'ssw_rate': round(ssw_rate, 3),
        'control_rate': round(ctrl_rate, 3),
        'trend_expected': round(trend_expected, 3),
    })

# Aggregate
valid_events = [e for e in event_results if e['rr_raw'] is not None and e['rr_raw'] > 0]
valid_adj = [e for e in event_results if e['rr_trend_adjusted'] is not None and e['rr_trend_adjusted'] > 0]

if valid_events:
    log_rrs_raw = [np.log(e['rr_raw']) for e in valid_events]
    mean_rr_raw = np.exp(np.mean(log_rrs_raw))
    t_stat_raw, p_raw = stats.ttest_1samp(log_rrs_raw, 0)
    
    # Direction test
    n_above = sum(1 for lr in log_rrs_raw if lr > 0)
    sign_p = stats.binomtest(n_above, len(log_rrs_raw), 0.5).pvalue
    
    print(f"\nRaw: geometric mean RR = {mean_rr_raw:.3f}, t-test P = {p_raw:.4f}")
    print(f"  {n_above}/{len(log_rrs_raw)} events with RR > 1 (sign test P = {sign_p:.4f})")

if valid_adj:
    log_rrs_adj = [np.log(e['rr_trend_adjusted']) for e in valid_adj]
    mean_rr_adj = np.exp(np.mean(log_rrs_adj))
    t_stat_adj, p_adj = stats.ttest_1samp(log_rrs_adj, 0)
    
    n_above_adj = sum(1 for lr in log_rrs_adj if lr > 0)
    sign_p_adj = stats.binomtest(n_above_adj, len(log_rrs_adj), 0.5).pvalue
    
    print(f"\nTrend-adjusted: geometric mean RR = {mean_rr_adj:.3f}, t-test P = {p_adj:.4f}")
    print(f"  {n_above_adj}/{len(log_rrs_adj)} events with RR > 1 (sign test P = {sign_p_adj:.4f})")

# Decade-stratified analysis
decades = {
    '1970s': (1970, 1979),
    '1980s': (1980, 1989),
    '1990s': (1990, 1999),
    '2000s': (2000, 2009),
    '2010s': (2010, 2019),
    '2020s': (2020, 2025),
}

decade_results = {}
for decade, (y1, y2) in decades.items():
    dec_events = [e for e in valid_events if y1 <= e['year'] <= y2]
    if dec_events:
        dec_rrs = [e['rr_raw'] for e in dec_events]
        decade_results[decade] = {
            'n_events': len(dec_events),
            'mean_rr': round(np.mean(dec_rrs), 3),
            'median_rr': round(np.median(dec_rrs), 3),
            'n_above_1': sum(1 for r in dec_rrs if r > 1),
        }
        print(f"{decade}: n={len(dec_events)}, mean RR={np.mean(dec_rrs):.3f}, {sum(1 for r in dec_rrs if r > 1)}/{len(dec_events)} above 1")

output = {
    'secular_trend': {
        'slope': round(slope, 5),
        'intercept': round(intercept, 3),
        'r_squared': round(r_value**2, 3),
        'p_value': round(p_value, 4),
        'interpretation': 'Accident rate trend across 1970-2023 winters'
    },
    'raw_aggregate': {
        'geometric_mean_rr': round(float(mean_rr_raw), 3) if valid_events else None,
        'p_value': round(float(p_raw), 4) if valid_events else None,
        'n_events': len(valid_events),
        'n_above_1': int(n_above) if valid_events else 0,
        'sign_test_p': round(float(sign_p), 4) if valid_events else None,
    },
    'trend_adjusted_aggregate': {
        'geometric_mean_rr': round(float(mean_rr_adj), 3) if valid_adj else None,
        'p_value': round(float(p_adj), 4) if valid_adj else None,
        'n_events': len(valid_adj),
        'n_above_1': int(n_above_adj) if valid_adj else 0,
        'sign_test_p': round(float(sign_p_adj), 4) if valid_adj else None,
    },
    'decade_stratified': decade_results,
    'per_event': event_results,
}

out_path = Path("data/results/41_accident_secular_trend.json")
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {out_path}")
