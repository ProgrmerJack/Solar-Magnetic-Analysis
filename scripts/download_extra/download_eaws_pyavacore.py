"""
Download Austrian + French + Italian avalanche danger data using pyAvaCore.
Focus on SSW event windows and matched controls.
"""
import subprocess, json, os, sys
from datetime import datetime, timedelta
import pandas as pd

# SSW events from our catalog  
ssw_dates = [
    '2000-03-20', '2001-02-11', '2001-12-30', '2003-01-18',
    '2004-01-05', '2006-01-21', '2007-02-24', '2008-02-22',
    '2009-01-24', '2010-02-09', '2010-03-24', '2013-01-07',
    '2018-02-12', '2019-01-01', '2021-01-05',
    '2012-01-11', '2023-02-16'
]

# Regions to download (focus on Alpine countries)
regions = [
    'AT-02',  # Tirol
    'AT-03',  # Salzburg
    'AT-04',  # Oberösterreich
    'AT-05',  # Steiermark
    'AT-06',  # Kärnten
    'AT-07',  # Vorarlberg
    'AT-08',  # Niederösterreich
    'FR-01',  # France
    'IT-23',  # Aosta
    'IT-25',  # Lombardia
    'IT-32-BZ',  # Südtirol
    'IT-32-TN',  # Trentino
    'IT-34',  # Friuli
    'IT-36',  # Veneto
    'IT-57',  # Piemonte
    'CH',     # Switzerland (for comparison)
    'DE-BY',  # Bavaria
    'SI',     # Slovenia
    'ES-CT-L',  # Catalonia/Aran
]

output_dir = 'data/cryosphere/eaws_bulletins'
os.makedirs(output_dir, exist_ok=True)

# For each SSW event, download ±20 day window
all_records = []
errors = []

for ssw_date_str in ssw_dates:
    ssw_date = datetime.strptime(ssw_date_str, '%Y-%m-%d')
    
    # Download ±20 day window
    start = ssw_date - timedelta(days=20)
    end = ssw_date + timedelta(days=20)
    
    date_range = '{}/{}'.format(start.strftime('%Y-%m-%d'), end.strftime('%Y-%m-%d'))
    
    print("Downloading SSW {} ({} to {})...".format(ssw_date_str, start.strftime('%Y-%m-%d'), end.strftime('%Y-%m-%d')))
    
    # Use pyAvaCore CLI
    cmd = [
        sys.executable, '-m', 'avacore',
        '--date', date_range,
        '--output', os.path.join(output_dir, 'ssw_{}'.format(ssw_date_str)),
        '--cli', 'n',
        '--lang', 'en'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print("  OK")
        else:
            print("  Error: {}".format(result.stderr[:200]))
            errors.append({'ssw': ssw_date_str, 'error': result.stderr[:200]})
    except subprocess.TimeoutExpired:
        print("  Timeout")
        errors.append({'ssw': ssw_date_str, 'error': 'timeout'})
    except Exception as e:
        print("  Exception: {}".format(str(e)[:100]))
        errors.append({'ssw': ssw_date_str, 'error': str(e)[:100]})

# Also download control windows (adjacent non-SSW winters)
control_dates = [
    '2005-01-15', '2011-01-15', '2014-01-15', '2015-01-15', 
    '2016-01-15', '2017-01-15', '2020-01-15', '2022-01-15',
    '2024-01-15',
]

for ctrl_date_str in control_dates:
    ctrl_date = datetime.strptime(ctrl_date_str, '%Y-%m-%d')
    start = ctrl_date - timedelta(days=20)
    end = ctrl_date + timedelta(days=20)
    date_range = '{}/{}'.format(start.strftime('%Y-%m-%d'), end.strftime('%Y-%m-%d'))
    
    print("Downloading control {} ({} to {})...".format(ctrl_date_str, start.strftime('%Y-%m-%d'), end.strftime('%Y-%m-%d')))
    
    cmd = [
        sys.executable, '-m', 'avacore',
        '--date', date_range,
        '--output', os.path.join(output_dir, 'ctrl_{}'.format(ctrl_date_str)),
        '--cli', 'n',
        '--lang', 'en'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print("  OK")
        else:
            print("  Error: {}".format(result.stderr[:200]))
    except Exception as e:
        print("  Exception: {}".format(str(e)[:100]))

print("\nErrors: {}".format(len(errors)))
for e in errors:
    print("  {} -> {}".format(e['ssw'], e['error']))

# List what was downloaded
print("\nDownloaded files:")
for root, dirs, files in os.walk(output_dir):
    for f in files[:5]:
        path = os.path.join(root, f)
        size = os.path.getsize(path)
        print("  {} ({} bytes)".format(path, size))
    if len(files) > 5:
        print("  ... and {} more files".format(len(files) - 5))
