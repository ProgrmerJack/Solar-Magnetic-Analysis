import pandas as pd
ssw = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw.index = ssw.index.tz_localize(None)
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
era5 = pd.read_parquet('data/processed/era5_swiss_alps_daily.parquet')

print('Panel range:', panel.index.min(), 'to', panel.index.max())
print('ERA5 range:', era5.index.min(), 'to', era5.index.max())
print()
print('SSW events in panel range:')
in_panel = ssw[(ssw.index >= panel.index.min()) & (ssw.index <= panel.index.max())]
for d in in_panel.index:
    in_era5 = era5.index.min() <= d <= era5.index.max()
    status = "YES" if in_era5 else "NO"
    print(f'  {d.date()} - ERA5: {status}')
n_with = sum(1 for d in in_panel.index if era5.index.min() <= d <= era5.index.max())
print(f'Total: {len(in_panel)} SSW events, {n_with} with ERA5')

missing = [d for d in in_panel.index if not (era5.index.min() <= d <= era5.index.max())]
if missing:
    print(f'Missing ERA5 for: {[str(d.date()) for d in missing]}')
    years = sorted(set(d.year for d in missing))
    print(f'Need ERA5 for years: {years}')
