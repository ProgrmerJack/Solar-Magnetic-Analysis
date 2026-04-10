"""Save Butler SSW Compendium data as structured CSV.
Source: https://csl.noaa.gov/groups/csl8/sswcompendium/majorevents.html
Butler et al. (2017), updated Sep 2023.
"""
import pandas as pd

# ERA5-based SSW dates from the compendium
ssw_events = [
    # (event_name, era5_date, enso, qbo)
    ("FEB 1979", "1979-02-22", "N", "W"),
    ("FEB 1980", "1980-02-29", "E", "E"),
    ("MAR 1981", "1981-03-04", "N", "W"),
    ("DEC 1981", "1981-12-04", "N", "E"),
    ("FEB 1984", "1984-02-24", "L", "W"),
    ("JAN 1985", "1985-01-01", "L", "E"),
    ("JAN 1987", "1987-01-23", "E", "W"),
    ("DEC 1987", "1987-12-08", "E", "W"),
    ("MAR 1988", "1988-03-14", "E", "W"),
    ("FEB 1989", "1989-02-21", "L", "W"),
    ("DEC 1998", "1998-12-15", "L", "E"),
    ("FEB 1999", "1999-02-26", "L", "E"),
    ("MAR 2000", "2000-03-20", "L", "W"),
    ("FEB 2001", "2001-02-11", "L", "W"),
    ("DEC 2001", "2001-12-30", "N", "E"),
    ("FEB 2002", "2002-02-17", "N", "E"),
    ("JAN 2003", "2003-01-18", "E", "W"),
    ("JAN 2004", "2004-01-05", "N", "E"),
    ("JAN 2006", "2006-01-21", "L", "E"),
    ("FEB 2007", "2007-02-24", "E", "W"),
    ("FEB 2008", "2008-02-22", "L", "E"),
    ("JAN 2009", "2009-01-24", "L", "W"),
    ("FEB 2010", "2010-02-09", "E", "W"),
    ("MAR 2010", "2010-03-24", "E", "W"),
    ("JAN 2013", "2013-01-06", "N", "E"),
    ("FEB 2018", "2018-02-12", "L", "W"),
    ("JAN 2019", "2019-01-01", "E", "E"),
    ("JAN 2021", "2021-01-05", "L", "W"),
    ("FEB 2023", "2023-02-16", "L", "W"),
]

df = pd.DataFrame(ssw_events, columns=['event_name', 'date', 'enso_phase', 'qbo_phase'])
df['date'] = pd.to_datetime(df['date'])

# Save
outfile = 'data/processed/atmospheric/butler_ssw_compendium_era5.csv'
df.to_csv(outfile, index=False)
print("Saved %d ERA5 SSW events to %s" % (len(df), outfile))
print(df.to_string())
