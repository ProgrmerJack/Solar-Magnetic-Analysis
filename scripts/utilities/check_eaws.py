import pandas as pd
eaws = pd.read_csv('data/cryosphere/european_alps/data_dmax.csv', sep=';')
print('Columns:', list(eaws.columns))
print('Shape:', eaws.shape)
eaws['date'] = pd.to_datetime(eaws['date'])
print('Date range:', eaws['date'].min(), 'to', eaws['date'].max())
print('Countries:', sorted(eaws['country'].unique()))
print('Centres:', sorted(eaws['forecastCenter'].unique()))
print('Years:', sorted(eaws['date'].dt.year.unique()))

ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_dates = ssw_cat.index.tz_localize(None)
print('\nSSW events in EAWS range:')
for d in ssw_dates:
    if d >= eaws['date'].min() and d <= eaws['date'].max():
        print(f'  {d.date()}')
