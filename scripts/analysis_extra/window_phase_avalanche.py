"""
Test whether the window-width dependence of the SSW-avalanche association
can be explained by the ERA5 phase decomposition.

If the 20-day window weakens because days 16-20 include the cold outbreak,
we should see:
1. Days 0-15: avalanche decrease (confirmed)
2. Days 16-20: avalanche increase or neutral (cold reversal)
3. This dilutes the 20-day signal
"""
import pandas as pd
import numpy as np
from scipy import stats

# Load data
panel = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\analysis_panel_v2.parquet')
ssw = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')

# Ensure compatible datetime
panel_dates = panel.index  # DatetimeIndex named 'time'
ssw_dates = ssw.index.tz_localize(None)

# Only use SSW events in study period
mask = (ssw_dates >= panel_dates.min()) & (ssw_dates <= panel_dates.max())
ssw_study = ssw_dates[mask]
print(f"SSW events in study period: {len(ssw_study)}")

def get_ssw_lag(date, ssw_dates):
    diffs = (date - ssw_dates).days
    past = diffs[diffs >= 0]
    future = diffs[diffs < 0]
    if len(past) > 0 and past.min() <= 30:
        return past.min()
    if len(future) > 0 and abs(future.max()) <= 15:
        return future.max()
    return None

y_col = 'dry_natural_size_1234'
print(f"Using column: {y_col}")

lags = []
for date in panel_dates:
    lag = get_ssw_lag(date, ssw_study)
    lags.append(lag)
panel['ssw_lag'] = lags

phases = {
    'Pre-SSW (d-15 to d-1)': (-15, -1),
    'Post-SSW (d0 to d+15)': (0, 15),
    'Late post (d+16 to d+20)': (16, 20),
    'Late post (d+16 to d+30)': (16, 30),
}

ctrl = panel[panel['ssw_lag'].isna()][y_col].dropna()
ctrl_mean = ctrl.mean()
print(f"\nControl mean ({y_col}): {ctrl_mean:.3f} (n={len(ctrl)})")

print(f"\n{'Phase':<30} {'Mean':>8} {'Diff':>8} {'n':>6} {'P(MW)':>10}")
print("-" * 72)
for name, (lo, hi) in phases.items():
    mask = (panel['ssw_lag'] >= lo) & (panel['ssw_lag'] <= hi)
    phase_data = panel[mask][y_col].dropna()
    if len(phase_data) > 5:
        diff = phase_data.mean() - ctrl_mean
        _, p = stats.mannwhitneyu(phase_data, ctrl, alternative='two-sided')
        print(f"{name:<30} {phase_data.mean():>8.3f} {diff:>+8.3f} {len(phase_data):>6} {p:>10.4f}")
    else:
        print(f"{name:<30} {'n<5':>8}")
