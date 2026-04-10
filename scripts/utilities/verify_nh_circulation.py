"""Re-run Phase 7 NH circulation analysis to verify paper numbers."""
import pandas as pd, numpy as np
from scipy import stats

# Load NCEP stratosphere data
ncep = pd.read_parquet('data/processed/ncep_stratosphere.parquet')
print("NCEP shape:", ncep.shape)
print("NCEP columns:", list(ncep.columns))
print("NCEP date range:", ncep.index.min(), "to", ncep.index.max())

# Load SSW catalog
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_cat.index = ssw_cat.index.tz_localize(None)

# Displacement events only
disp_events = ssw_cat[ssw_cat['type'] == 'displacement'].index
print("\nDisplacement events:", len(disp_events))
for d in disp_events:
    print(" ", d.date())

# Check for NH circulation variables
nh_vars = [c for c in ncep.columns if 'hgt' in c.lower() or 'slp' in c.lower() or '500' in c or '850' in c]
print("\nPossible NH vars:", nh_vars[:20])

# Also check for Z500, SLP, U850 in other datasets
import glob
all_parquets = glob.glob('data/processed/**/*.parquet', recursive=True)
print("\nAll parquet files:")
for f in all_parquets[:30]:
    print(" ", f)

# Check if panel has these
panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
nh_panel = [c for c in panel.columns if 'z500' in c.lower() or 'slp' in c.lower() or 'u850' in c.lower() or 'hgt' in c.lower()]
print("\nPanel NH vars:", nh_panel)

# Check NCEP for 500hPa height and 850hPa wind
ncep_hgt = [c for c in ncep.columns if 'hgt' in c.lower()]
ncep_uwnd = [c for c in ncep.columns if 'uwnd' in c.lower()]
print("\nNCEP hgt columns:", ncep_hgt)
print("NCEP uwnd columns:", ncep_uwnd)

# Compute composites using NCEP data
# 500 hPa geopotential height and 850 hPa zonal wind
if 'hgt_m_500hPa' in ncep.columns and 'uwnd_ms_850hPa' in ncep.columns:
    # Add DOY climatology
    ncep['doy'] = ncep.index.dayofyear
    for var in ['hgt_m_500hPa', 'uwnd_ms_850hPa']:
        clim = ncep.groupby('doy')[var].mean()
        ncep[var + '_anom'] = ncep[var] - ncep['doy'].map(clim)
    
    # Composite for displacement events: days 0 to +30
    all_ssw_days = []
    all_ctrl_days = []
    for onset in disp_events:
        # SSW window: days 0 to +30
        ssw_start = onset
        ssw_end = onset + pd.Timedelta(days=30)
        mask = (ncep.index >= ssw_start) & (ncep.index <= ssw_end)
        ssw_days = ncep.loc[mask]
        all_ssw_days.append(ssw_days)
        
        # Control: same DOY in other years
        for yr_offset in [-2, -1, 1, 2]:
            ctrl_start = onset + pd.DateOffset(years=yr_offset)
            ctrl_end = ctrl_start + pd.Timedelta(days=30)
            ctrl_mask = (ncep.index >= ctrl_start) & (ncep.index <= ctrl_end)
            ctrl_days = ncep.loc[ctrl_mask]
            all_ctrl_days.append(ctrl_days)
    
    ssw_all = pd.concat(all_ssw_days)
    ctrl_all = pd.concat(all_ctrl_days)
    
    print("\n=== NH CIRCULATION: DISPLACEMENT SSW (days 0-30) ===")
    for var in ['hgt_m_500hPa_anom', 'uwnd_ms_850hPa_anom']:
        ssw_mean = ssw_all[var].mean()
        ctrl_mean = ctrl_all[var].mean()
        diff = ssw_mean - ctrl_mean
        t, p = stats.ttest_ind(ssw_all[var].dropna(), ctrl_all[var].dropna())
        label = var.replace('_anom','').replace('_',' ')
        print("  %s: SSW=%.3f, Ctrl=%.3f, diff=%.3f, t=%.3f, P=%.6f" % (
            label, ssw_mean, ctrl_mean, diff, t, p))
    
    # Also compute simple SSW window anomaly (deseasonalised)
    print("\n=== DESEASONALISED SSW ANOMALY (days 0-30, displacement) ===")
    for var in ['hgt_m_500hPa_anom', 'uwnd_ms_850hPa_anom']:
        ssw_mean = ssw_all[var].mean()
        # Test against zero (the climatology is zero by construction)
        t, p = stats.ttest_1samp(ssw_all[var].dropna(), 0)
        print("  %s: mean anom=%.3f, t=%.3f, P=%.6f" % (var, ssw_mean, t, p))

else:
    print("\n500hPa/850hPa not found. Available:", list(ncep.columns))
    # Try computing from what we have
    for c in ncep.columns:
        print("  ", c)
