import pandas as pd
import json

# Check European Alps dmax data
df = pd.read_csv('C:/Users/Jack0/Solar-Magnetic-Analysis/data/cryosphere/european_alps/data_dmax.csv', nrows=10)
print('European Alps dmax columns:', list(df.columns))
print(df.head(3))
print(f'Shape: {df.shape}')

# Full file
df_full = pd.read_csv('C:/Users/Jack0/Solar-Magnetic-Analysis/data/cryosphere/european_alps/data_dmax.csv')
print(f'\nFull shape: {df_full.shape}')
if 'country' in [c.lower() for c in df_full.columns]:
    country_col = [c for c in df_full.columns if c.lower() == 'country'][0]
    print(f'Countries: {df_full[country_col].unique()}')
if 'region' in [c.lower() for c in df_full.columns]:
    region_col = [c for c in df_full.columns if c.lower() == 'region'][0]
    print(f'Regions: {df_full[region_col].nunique()} unique')
    print(f'Sample: {df_full[region_col].unique()[:15]}')

# Check date range
date_cols = [c for c in df_full.columns if 'date' in c.lower() or 'time' in c.lower() or 'day' in c.lower()]
print(f'Date columns: {date_cols}')
if date_cols:
    print(f'Date range: {df_full[date_cols[0]].min()} to {df_full[date_cols[0]].max()}')

# Check EAWS regions
with open('C:/Users/Jack0/Solar-Magnetic-Analysis/data/cryosphere/eaws/eaws_regions.geojson') as f:
    geo = json.load(f)
n_features = len(geo['features'])
print(f'\nEAWS regions: {n_features} features')
if geo['features']:
    props = geo['features'][0]['properties']
    print(f'Properties keys: {list(props.keys())}')
    # Show some region IDs
    ids = [f['properties'].get('id', f['properties'].get('ID', 'N/A')) for f in geo['features'][:20]]
    print(f'First 20 IDs: {ids}')
