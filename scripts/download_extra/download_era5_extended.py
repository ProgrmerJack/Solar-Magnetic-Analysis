"""
Download ERA5 daily data for Swiss Alps region for missing years.
Needed: 1998-2003 and 2018-2019 to complement existing 2004-2013 data.
Variables: 2m temperature, total precipitation, snowfall, snow depth, 10m winds.
Region: Swiss Alps box [45.5-47.5N, 6-11E]
"""
import cdsapi
import os
import sys

c = cdsapi.Client()

# Swiss Alps bounding box: [N, W, S, E]
area = [47.5, 6, 45.5, 11]

# Variables matching existing ERA5 data
variables = [
    '2m_temperature',
    'total_precipitation',
    'snowfall',
    'snow_depth',
    '10m_u_component_of_wind',
    '10m_v_component_of_wind',
]

outdir = 'data/raw/era5_extended'
os.makedirs(outdir, exist_ok=True)

# Years we need
years_needed = [1998, 1999, 2001, 2002, 2003, 2018, 2019]
# Also 2000 for control years
years_needed.append(2000)
years_needed.sort()

# Winter months only (Oct-Apr) to save download time
months = ['01', '02', '03', '04', '10', '11', '12']

for year in years_needed:
    outfile = os.path.join(outdir, f'era5_swiss_{year}.nc')
    if os.path.exists(outfile):
        print(f'{year}: already downloaded, skipping')
        continue
    
    print(f'Downloading ERA5 for {year}...')
    try:
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': variables,
                'year': str(year),
                'month': months,
                'day': [f'{d:02d}' for d in range(1, 32)],
                'time': '12:00',  # Daily snapshot at noon
                'area': area,
                'format': 'netcdf',
            },
            outfile
        )
        size_mb = os.path.getsize(outfile) / 1e6
        print(f'  {year}: downloaded ({size_mb:.1f} MB)')
    except Exception as e:
        print(f'  {year}: FAILED - {e}')

print('\nDone!')
