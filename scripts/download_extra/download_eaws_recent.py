"""
Download EAWS bulletins for RECENT SSW events only (2018+).
The EAWS/ALBINA system was launched around 2017-2018.
Download single dates to avoid timeouts.
"""
import subprocess, json, os, sys, glob
from datetime import datetime, timedelta

output_dir = 'data/cryosphere/eaws_bulletins'
os.makedirs(output_dir, exist_ok=True)

# Recent SSW events in EAWS era
ssw_events = [
    '2018-02-12',
    '2019-01-01', 
    '2021-01-05',
    '2023-02-16',
]

# Also some controls
control_periods = [
    '2020-01-15',
    '2022-01-15',  
    '2024-01-15',
]

all_dates = []
for ssw in ssw_events:
    base = datetime.strptime(ssw, '%Y-%m-%d')
    for offset in range(-20, 21):
        d = base + timedelta(days=offset)
        all_dates.append(('ssw_' + ssw, d))

for ctrl in control_periods:
    base = datetime.strptime(ctrl, '%Y-%m-%d')
    for offset in range(-20, 21):
        d = base + timedelta(days=offset)
        all_dates.append(('ctrl_' + ctrl, d))

print("Total date-downloads needed: {}".format(len(all_dates)))
print("Starting download...")

# Download in batches using date ranges (one week at a time)
def download_range(label, start, end, outdir):
    date_range = '{}/{}'.format(start.strftime('%Y-%m-%d'), end.strftime('%Y-%m-%d'))
    cmd = [
        sys.executable, '-m', 'avacore',
        '--date', date_range,
        '--output', outdir,
        '--cli', 'n',
        '--lang', 'en'
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
        return result.returncode == 0, result.stderr[:100] if result.returncode != 0 else ''
    except subprocess.TimeoutExpired:
        return False, 'timeout'
    except Exception as e:
        return False, str(e)[:80]

# Group by event and download in weekly chunks
events = ssw_events + control_periods
labels = ['ssw_' + s for s in ssw_events] + ['ctrl_' + c for c in control_periods]

for label, center_str in zip(labels, events):
    center = datetime.strptime(center_str, '%Y-%m-%d')
    outdir = os.path.join(output_dir, label)
    os.makedirs(outdir, exist_ok=True)
    
    # Download in 7-day chunks
    for week_start_offset in range(-21, 21, 7):
        week_start = center + timedelta(days=week_start_offset)
        week_end = week_start + timedelta(days=6)
        
        date_range = '{}/{}'.format(week_start.strftime('%Y-%m-%d'), week_end.strftime('%Y-%m-%d'))
        sys.stdout.write("  {} [{} to {}]... ".format(label, week_start.strftime('%m-%d'), week_end.strftime('%m-%d')))
        sys.stdout.flush()
        
        ok, err = download_range(label, week_start, week_end, outdir)
        if ok:
            # Count downloaded files
            files = glob.glob(os.path.join(outdir, '*.json'))
            sys.stdout.write("OK ({} files total)\n".format(len(files)))
        else:
            sys.stdout.write("ERROR: {}\n".format(err))
        sys.stdout.flush()

# Summary
print("\n=== Summary ===")
total_files = 0
for label in labels:
    outdir = os.path.join(output_dir, label)
    files = glob.glob(os.path.join(outdir, '*.json'))
    total_files += len(files)
    if files:
        # Read a sample to check structure
        with open(files[0], 'r', encoding='utf-8') as f:
            try:
                data = json.load(f)
                if isinstance(data, dict):
                    print("  {} -> {} files, keys: {}".format(label, len(files), list(data.keys())[:5]))
                elif isinstance(data, list):
                    print("  {} -> {} files, {} items per file".format(label, len(files), len(data)))
            except:
                print("  {} -> {} files (not valid JSON)".format(label, len(files)))
    else:
        print("  {} -> 0 files".format(label))

print("\nTotal downloaded: {} JSON files".format(total_files))
