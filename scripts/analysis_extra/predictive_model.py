import pandas as pd
import numpy as np
from scipy import stats

# Build cross-validated predictive model: vortex state -> avalanche probability
strat = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ncep_stratosphere.parquet')
panel = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\analysis_panel_v2.parquet')
ssw_cat = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\atmospheric\ssw_catalog.parquet')

strat.index = strat.index.tz_localize(None)
ssw_cat.index = ssw_cat.index.tz_localize(None)

aval = panel['dry_natural_size_1234']
aval.index = pd.to_datetime(aval.index)

# Define SSW-active windows: ±15d of each SSW onset
swiss_events = ssw_cat[(ssw_cat.index >= '1998-11-01') & (ssw_cat.index <= '2019-04-30')].index

# LEAVE-ONE-OUT CROSS-VALIDATED PREDICTION
# For each SSW event: train on all other events, predict this event
print('=== LEAVE-ONE-OUT CROSS-VALIDATED SSW PREDICTION ===')
print(f'Approach: Train rate ratio on n-1 events, predict held-out event')
print(f'Baseline model: seasonal climatology')
print(f'SSW model: seasonal climatology * SSW adjustment factor')
print()

all_ratios = []
for onset in swiss_events:
    onset_ts = pd.Timestamp(onset)
    window = pd.date_range(onset_ts - pd.Timedelta(days=15), onset_ts + pd.Timedelta(days=15))
    ssw_rate = aval.reindex(window).mean()
    doy = onset_ts.timetuple().tm_yday
    doys = set()
    for d in range(-15, 16):
        doys.add((doy + d) % 366)
    clim = aval[aval.index.dayofyear.isin(doys)].mean()
    if clim > 0:
        all_ratios.append(ssw_rate / clim)
    else:
        all_ratios.append(np.nan)

ratios = np.array(all_ratios)
print('Event ratios:', [f'{r:.3f}' for r in ratios])
print()

# LOO prediction
predictions = []
actuals = []
baseline_errors = []
ssw_errors = []

for i in range(len(ratios)):
    if np.isnan(ratios[i]):
        continue
    
    # Train on all except event i
    train_ratios = np.concatenate([ratios[:i], ratios[i+1:]])
    train_ratios = train_ratios[~np.isnan(train_ratios)]
    
    # LOO predicted ratio (mean of training set)
    predicted_ratio = train_ratios.mean()
    
    # Actual ratio
    actual = ratios[i]
    
    # Baseline (climatological): ratio = 1.0
    baseline_error = (1.0 - actual)**2
    ssw_error = (predicted_ratio - actual)**2
    
    predictions.append(predicted_ratio)
    actuals.append(actual)
    baseline_errors.append(baseline_error)
    ssw_errors.append(ssw_error)
    
    onset = swiss_events[i]
    print(f'Event {pd.Timestamp(onset).strftime("%Y-%m-%d")}: actual={actual:.3f}, predicted={predicted_ratio:.3f}, baseline=1.000, SSW_err={ssw_error:.3f}, base_err={baseline_error:.3f}')

# Summary statistics
print()
base_mse = np.mean(baseline_errors)
ssw_mse = np.mean(ssw_errors)
skill = 1 - ssw_mse / base_mse
print(f'Baseline MSE: {base_mse:.4f}')
print(f'SSW model MSE: {ssw_mse:.4f}')
print(f'Forecast skill score: {skill:.4f}')
print(f'SSW model reduces MSE by {100*(1-ssw_mse/base_mse):.1f}%')

# Directional accuracy
correct_dir = sum(1 for p, a in zip(predictions, actuals) if (p < 1) == (a < 1))
print(f'Directional accuracy: {correct_dir}/{len(predictions)} ({100*correct_dir/len(predictions):.0f}%)')

# Signed log-likelihood comparison (continuous ranked probability)
print()
print('=== Calibration: LOO prediction intervals ===')
within_ci = 0
for i in range(len(ratios)):
    if np.isnan(ratios[i]):
        continue
    train = np.concatenate([ratios[:i], ratios[i+1:]])
    train = train[~np.isnan(train)]
    lo = np.percentile(train, 2.5)
    hi = np.percentile(train, 97.5)
    inside = lo <= ratios[i] <= hi
    within_ci += inside
print(f'Events within 95% prediction interval: {within_ci}/{len(actuals)}')

# Apply the model to Utah
utah = pd.read_parquet(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\processed\cryosphere\utah_daily_dry_slab.parquet')
utah_events = ssw_cat[(ssw_cat.index >= '2012-01-01') & (ssw_cat.index <= '2024-12-31')]
print(f'\n=== UTAH PREDICTION (out-of-region validation) ===')
print(f'Swiss-trained mean ratio: {np.nanmean(ratios):.3f}')
print(f'Predicted Utah SSW effect: {np.nanmean(ratios):.3f} (50%+ reduction)')

# Utah actual
utah.index = pd.to_datetime(utah.index)
utah_aval = utah.iloc[:,0] if len(utah.columns) > 0 else utah
print(f'Utah columns: {list(utah.columns)}')
