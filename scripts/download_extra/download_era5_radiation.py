"""Download ERA5 radiation and additional surface variables for Swiss Alps.
Variables: surface net solar radiation, surface net thermal radiation,
total cloud cover, skin temperature, boundary layer height.
"""
import cdsapi
import os

c = cdsapi.Client()

# Define Swiss Alps bounding box
AREA = [48, 6, 46, 11]  # N, W, S, E
YEARS = [str(y) for y in range(1998, 2020)]
MONTHS = ['01', '02', '03', '04', '11', '12']  # Winter months

output_dir = 'data/raw/era5_radiation'
os.makedirs(output_dir, exist_ok=True)

# Download in chunks by year to avoid timeout
for year in YEARS:
    outfile = os.path.join(output_dir, f'era5_radiation_{year}.nc')
    if os.path.exists(outfile):
        print(f"Skipping {year} (already exists)")
        continue
    
    print(f"Downloading ERA5 radiation variables for {year}...")
    try:
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': [
                    'surface_net_solar_radiation',
                    'surface_net_thermal_radiation', 
                    'total_cloud_cover',
                    'skin_temperature',
                    'boundary_layer_height',
                    'surface_sensible_heat_flux',
                    'surface_latent_heat_flux',
                    '2m_dewpoint_temperature',
                ],
                'year': year,
                'month': MONTHS,
                'day': [str(d).zfill(2) for d in range(1, 32)],
                'time': ['00:00', '06:00', '12:00', '18:00'],
                'area': AREA,
                'format': 'netcdf',
            },
            outfile
        )
        print(f"  Saved: {outfile}")
    except Exception as e:
        print(f"  Error for {year}: {e}")
        # Try with fewer variables
        try:
            c.retrieve(
                'reanalysis-era5-single-levels',
                {
                    'product_type': 'reanalysis',
                    'variable': [
                        'surface_net_solar_radiation',
                        'surface_net_thermal_radiation',
                        'total_cloud_cover',
                        'skin_temperature',
                    ],
                    'year': year,
                    'month': MONTHS,
                    'day': [str(d).zfill(2) for d in range(1, 32)],
                    'time': ['00:00', '06:00', '12:00', '18:00'],
                    'area': AREA,
                    'format': 'netcdf',
                },
                outfile
            )
            print(f"  Saved (reduced vars): {outfile}")
        except Exception as e2:
            print(f"  Failed completely: {e2}")

print("Done downloading ERA5 radiation data.")
