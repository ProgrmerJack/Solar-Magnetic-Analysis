import pandas as pd

# Parse correctly with semicolon delimiter
df = pd.read_csv('data/cryosphere/european_alps/data_dmax.csv', sep=';')
print('Shape:', df.shape)
print('Columns:', list(df.columns))
print('\nFirst 5 rows:')
print(df.head())
print('\nCountry counts:')
print(df['country'].value_counts())
print('\nForecast center counts:')
print(df['forecastCenter'].value_counts())
print('\nDate range:', df['date'].min(), 'to', df['date'].max())
print('\nDanger level distribution:')
print(df['dangerLevelMax'].value_counts().sort_index())
print('\nUnique warning regions per country:')
for c in df['country'].unique():
    n = df[df['country']==c]['warningRegion'].nunique()
    regions = df[df['country']==c]['warningRegion'].unique()[:5]
    dates = df[df['country']==c]['date']
    print('  %s: %d regions, dates %s to %s' % (c, n, dates.min(), dates.max()))
    print('    Sample regions:', list(regions))

# Swiss snowpack data
print('\n\n=== Swiss Snowpack Data ===')
sf = pd.read_csv('data/cryosphere/swiss_snowpack/data_rf2_tidy.csv', nrows=5)
print('Columns:', list(sf.columns))
print(sf.head())
