"""
Process extended ERA5 netCDF files (1998-2003, 2018-2019) and merge with existing
daily parquet (2004-2013) to create unified ERA5 dataset covering 1998-2019.
Then rerun sintering model with all SSW events in the extended period.
"""
import xarray as xr
import pandas as pd
import numpy as np
import os, json, glob

BASE = os.path.join(os.path.dirname(__file__), '..')

# === Step 1: Process new netCDF files ===
print("=" * 60)
print("Step 1: Processing extended ERA5 netCDF files")
print("=" * 60)

nc_files = sorted(glob.glob(os.path.join(BASE, 'data', 'raw', 'era5_extended', '*.nc')))
print(f"Found {len(nc_files)} netCDF files")

new_frames = []
for f in nc_files:
    year = os.path.basename(f).replace('era5_swiss_', '').replace('.nc', '')
    
    import zipfile, tempfile
    if zipfile.is_zipfile(f):
        with zipfile.ZipFile(f) as zf:
            nc_names = [n for n in zf.namelist() if n.endswith('.nc')]
            print(f"  {year}: ZIP with {nc_names}")
            
            dfs_year = []
            for nc_name in nc_names:
                with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as tmp:
                    tmp.write(zf.read(nc_name))
                    tmp_path = tmp.name
                try:
                    ds = xr.open_dataset(tmp_path)
                    spatial_dims = [d for d in ds.dims if d not in ['time', 'valid_time', 'date', 'number']]
                    time_dim = 'valid_time' if 'valid_time' in ds.dims else 'time'
                    var_names = list(ds.data_vars)
                    print(f"    {nc_name}: vars={var_names[:6]}, dims={list(ds.dims)}")
                    
                    if spatial_dims:
                        ds_mean = ds.mean(dim=spatial_dims)
                    else:
                        ds_mean = ds
                    
                    df = ds_mean.to_dataframe().reset_index()
                    if time_dim in df.columns:
                        df['date'] = pd.to_datetime(df[time_dim]).dt.date
                    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
                    daily = df.groupby('date')[numeric_cols].mean().reset_index()
                    daily['date'] = pd.to_datetime(daily['date'])
                    daily = daily.set_index('date')
                    dfs_year.append(daily)
                    ds.close()
                finally:
                    os.unlink(tmp_path)
            
            # Merge instant and accumulated
            if dfs_year:
                merged_year = pd.concat(dfs_year, axis=1)
                merged_year = merged_year.loc[:, ~merged_year.columns.duplicated()]
    else:
        ds = xr.open_dataset(f)
        spatial_dims = [d for d in ds.dims if d not in ['time', 'valid_time']]
        time_dim = 'valid_time' if 'valid_time' in ds.dims else 'time'
        if spatial_dims:
            ds_mean = ds.mean(dim=spatial_dims)
        else:
            ds_mean = ds
        df = ds_mean.to_dataframe().reset_index()
        df['date'] = pd.to_datetime(df[time_dim]).dt.date
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        merged_year = df.groupby('date')[numeric_cols].mean().reset_index()
        merged_year['date'] = pd.to_datetime(merged_year['date'])
        merged_year = merged_year.set_index('date')
        ds.close()
    
    # Rename variables to match existing schema
    rename_map = {}
    for col in merged_year.columns:
        cl = col.lower()
        if 'tp' in cl or 'total_precipitation' in cl:
            rename_map[col] = 'tp_mm'
        elif 'sf' in cl or 'snowfall' in cl:
            rename_map[col] = 'sf_mm'
        elif 't2m' in cl or '2m_temperature' in cl:
            rename_map[col] = 't2m_K'
        elif 'sd' in cl or 'snow_depth' in cl:
            rename_map[col] = 'sd_m'
        elif 'u10' in cl or '10m_u' in cl:
            rename_map[col] = 'u10'
        elif 'v10' in cl or '10m_v' in cl:
            rename_map[col] = 'v10'
    
    merged_year = merged_year.rename(columns=rename_map)
    
    # Convert precipitation from m to mm if needed
    if 'tp_mm' in merged_year.columns and merged_year['tp_mm'].max() < 1:
        merged_year['tp_mm'] = merged_year['tp_mm'] * 1000
    if 'sf_mm' in merged_year.columns and merged_year['sf_mm'].max() < 1:
        merged_year['sf_mm'] = merged_year['sf_mm'] * 1000
    
    # Compute wind speed
    if 'u10' in merged_year.columns and 'v10' in merged_year.columns:
        merged_year['wind_speed'] = np.sqrt(merged_year['u10']**2 + merged_year['v10']**2)
    
    # Add day of year
    merged_year['doy'] = merged_year.index.dayofyear
    
    new_frames.append(merged_year)
    print(f"    -> {len(merged_year)} days, cols: {list(merged_year.columns)[:8]}")

new_era5 = pd.concat(new_frames)
print(f"\nNew ERA5: {len(new_era5)} days, {new_era5.index.min()} to {new_era5.index.max()}")

# === Step 2: Load and merge with existing ===
print("\n" + "=" * 60)
print("Step 2: Merging with existing ERA5 (2004-2013)")
print("=" * 60)

existing = pd.read_parquet(os.path.join(BASE, 'data', 'processed', 'era5_swiss_alps_daily.parquet'))
existing.index = pd.to_datetime(existing.index)
print(f"Existing: {len(existing)} days, {existing.index.min()} to {existing.index.max()}")

# Keep only common columns
common_cols = sorted(set(new_era5.columns) & set(existing.columns))
print(f"Common columns: {common_cols}")

merged = pd.concat([new_era5[common_cols], existing[common_cols]])
merged = merged.sort_index()
merged = merged[~merged.index.duplicated(keep='first')]

# Recompute anomalies based on DOY climatology
merged['doy'] = merged.index.dayofyear
for var in ['tp_mm', 'sf_mm']:
    if var in merged.columns:
        clim = merged.groupby('doy')[var].transform('mean')
        merged[f'{var}_anom'] = merged[var] - clim

print(f"Merged: {len(merged)} days, {merged.index.min()} to {merged.index.max()}")

# Save
out_path = os.path.join(BASE, 'data', 'processed', 'era5_swiss_alps_extended.parquet')
merged.to_parquet(out_path)
print(f"Saved to {out_path}")

# === Step 3: Identify SSW events in extended period ===
print("\n" + "=" * 60)
print("Step 3: SSW events in extended ERA5 period")
print("=" * 60)

ssw_cat = pd.read_parquet(os.path.join(BASE, 'data', 'processed', 'atmospheric', 'ssw_catalog.parquet'))
ssw_dates = pd.to_datetime(ssw_cat.index).tz_localize(None)

era5_start = merged.index.min()
era5_end = merged.index.max()
ssw_in_era5 = [d for d in ssw_dates if era5_start <= d <= era5_end]
print(f"SSW events in ERA5 coverage ({era5_start.date()} to {era5_end.date()}):")
for d in sorted(ssw_in_era5):
    print(f"  {d.date()}")
print(f"Total: {len(ssw_in_era5)} events (was 8 with 2004-2013 only)")

# === Step 4: Rerun sintering model with extended data ===
print("\n" + "=" * 60)
print("Step 4: Sintering model with extended ERA5")
print("=" * 60)

# Arrhenius sintering model parameters
Ea = 0.5  # eV activation energy
kB = 8.617e-5  # eV/K
T_ref = 263.15  # -10°C reference
tau_sinter = 7  # days integration window

results = []
for ssw_d in sorted(ssw_in_era5):
    # Get temperature for ±25 day window
    window_start = ssw_d - pd.Timedelta(days=25)
    window_end = ssw_d + pd.Timedelta(days=25)
    
    era5_window = merged.loc[window_start:window_end]
    if 't2m_K' not in era5_window.columns or len(era5_window) < 30:
        print(f"  {ssw_d.date()}: insufficient data ({len(era5_window)} days)")
        continue
    
    T_series = era5_window['t2m_K'].dropna()
    if len(T_series) < 30:
        print(f"  {ssw_d.date()}: insufficient temperature data")
        continue
    
    # Compute sintering rate relative to reference
    rate = np.exp(-Ea / kB * (1/T_series.values - 1/T_ref))
    
    # SSW window: days -5 to +15
    ssw_mask = (era5_window.index >= ssw_d - pd.Timedelta(days=5)) & \
               (era5_window.index <= ssw_d + pd.Timedelta(days=15))
    ctrl_mask = ~ssw_mask
    
    ssw_idx = era5_window.index[ssw_mask]
    ctrl_idx = era5_window.index[ctrl_mask]
    
    ssw_temp = T_series.reindex(ssw_idx).dropna()
    ctrl_temp = T_series.reindex(ctrl_idx).dropna()
    
    if len(ssw_temp) < 5 or len(ctrl_temp) < 5:
        continue
    
    # Mean temperature anomaly
    delta_T = ssw_temp.mean() - ctrl_temp.mean()
    
    # Integrated sintering enhancement
    ssw_rate = np.exp(-Ea / kB * (1/ssw_temp.values - 1/T_ref))
    ctrl_rate = np.exp(-Ea / kB * (1/ctrl_temp.values - 1/T_ref))
    
    enhancement = (ssw_rate.mean() / ctrl_rate.mean() - 1) * 100  # percent
    
    results.append({
        'ssw_date': str(ssw_d.date()),
        'delta_T_K': round(delta_T, 3),
        'sintering_enhancement_pct': round(enhancement, 2),
        'ssw_mean_T': round(ssw_temp.mean(), 2),
        'ctrl_mean_T': round(ctrl_temp.mean(), 2),
        'n_ssw_days': len(ssw_temp),
        'n_ctrl_days': len(ctrl_temp)
    })
    
    sign = '+' if enhancement > 0 else ''
    print(f"  {ssw_d.date()}: ΔT={delta_T:+.2f}K, sintering {sign}{enhancement:.1f}%")

# Summary statistics
enhancements = [r['sintering_enhancement_pct'] for r in results]
positive = sum(1 for e in enhancements if e > 0)
mean_enh = np.mean(enhancements)
median_enh = np.median(enhancements)

# Statistical test
from scipy import stats
t_stat, p_val = stats.ttest_1samp(enhancements, 0)
sign_p = stats.binom_test(positive, len(enhancements), 0.5) if hasattr(stats, 'binom_test') else \
         stats.binomtest(positive, len(enhancements), 0.5).pvalue

print(f"\n{'='*60}")
print(f"SINTERING MODEL RESULTS (n={len(results)} SSW events)")
print(f"{'='*60}")
print(f"Mean enhancement: {mean_enh:.1f}%")
print(f"Median enhancement: {median_enh:.1f}%")
print(f"Positive: {positive}/{len(results)} ({100*positive/len(results):.0f}%)")
print(f"t-test vs 0: t={t_stat:.3f}, P={p_val:.6f}")
print(f"Sign test: P={sign_p:.6f}")
print(f"Range: [{min(enhancements):.1f}%, {max(enhancements):.1f}%]")

# Save results
output = {
    'n_events': len(results),
    'n_era5_years': f"{era5_start.year}-{era5_end.year}",
    'mean_enhancement_pct': round(mean_enh, 2),
    'median_enhancement_pct': round(median_enh, 2),
    'positive_fraction': f"{positive}/{len(results)}",
    't_test_p': round(p_val, 6),
    'sign_test_p': round(sign_p, 6),
    'per_event': results
}
out_json = os.path.join(BASE, 'data', 'results', 'sintering_extended.json')
with open(out_json, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to {out_json}")
