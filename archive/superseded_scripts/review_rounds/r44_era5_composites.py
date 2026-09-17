"""
ERA5 spatial composite maps for SSW events.
Creates SLP and T2m composite anomaly maps for the Swiss Alpine domain
using ERA5 data (2004-2013), addressing reviewer requests for spatial
visualization of the circulation pattern during SSW windows.

Also creates a mechanism timeline figure showing onset-aligned evolution
of key variables from existing analysis results.

Output:
  data/figures/fig_slp_t2m_composite.pdf
  data/figures/fig_mechanism_timeline.pdf
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

LOO_CV_FILE = ROOT / "scripts" / "review_rounds" / "r43_threshold_loo_cv.json"
RESULTS = ROOT / "data" / "results"
FIGDIR  = ROOT / "data" / "figures"

# SSW event dates (Butler et al. catalog, within ERA5 Swiss Alps range 2004-2013)
SSW_DATES_ERA5 = [
    "2004-01-05", "2006-01-21", "2007-02-24", "2008-02-22",
    "2009-01-24", "2010-02-09", "2012-01-11", "2013-01-07"
]

ALL_SSW_DATES = [
    "1998-12-15", "1999-02-26", "2001-02-11", "2001-12-30",
    "2002-02-17", "2003-01-18", "2004-01-05", "2006-01-21",
    "2007-02-24", "2008-02-22", "2009-01-24", "2010-02-09",
    "2012-01-11", "2013-01-07", "2018-02-12", "2019-01-01"
]

def load_era5_composites():
    """Load ERA5 data and compute SSW vs control composites."""
    try:
        import xarray as xr
    except ImportError:
        print("xarray not available, using precomputed results")
        return None, None
    
    era5_dir = ROOT / "data" / "atmospheric" / "era5" / "swiss_alps"
    
    ssw_t2m_maps = []
    ctrl_t2m_maps = []
    ssw_slp_maps = []
    ctrl_slp_maps = []
    
    for date_str in SSW_DATES_ERA5:
        onset = pd.Timestamp(date_str)
        year = onset.year
        
        fpath = era5_dir / f"era5_swiss_alps_{year}.nc"
        if not fpath.exists():
            print(f"  Missing: {fpath}")
            continue
            
        try:
            ds = xr.open_dataset(fpath)
        except Exception as e:
            print(f"  Error loading {fpath}: {e}")
            continue
        
        # Identify time dimension
        time_var = None
        for v in ['time', 'valid_time']:
            if v in ds.dims or v in ds.coords:
                time_var = v
                break
        if time_var is None:
            print(f"  No time variable found in {fpath}")
            ds.close()
            continue
        
        times = pd.DatetimeIndex(ds[time_var].values)
        
        # SSW window: ±15 days
        ssw_mask = (times >= onset - pd.Timedelta(days=15)) & \
                   (times <= onset + pd.Timedelta(days=15))
        
        # Control: same DOY ±3 days from non-SSW periods
        doy_center = onset.dayofyear
        ctrl_mask = np.zeros(len(times), dtype=bool)
        for t_idx, t in enumerate(times):
            doy_diff = abs(t.dayofyear - doy_center)
            if doy_diff > 180:
                doy_diff = 365 - doy_diff
            if doy_diff <= 3 and not ssw_mask[t_idx]:
                ctrl_mask[t_idx] = True
        
        # Extract T2m
        t2m_var = None
        for v in ['t2m', '2t', 'VAR_2T', 'd2m']:
            if v in ds.data_vars:
                t2m_var = v
                break
        
        slp_var = None
        for v in ['msl', 'sp', 'mean_sea_level_pressure', 'surface_pressure']:
            if v in ds.data_vars:
                slp_var = v
                break
        
        if t2m_var and ssw_mask.sum() > 0 and ctrl_mask.sum() > 0:
            ssw_t2m = ds[t2m_var].isel({time_var: ssw_mask}).mean(dim=time_var).values
            ctrl_t2m = ds[t2m_var].isel({time_var: ctrl_mask}).mean(dim=time_var).values
            ssw_t2m_maps.append(ssw_t2m)
            ctrl_t2m_maps.append(ctrl_t2m)
        
        if slp_var and ssw_mask.sum() > 0 and ctrl_mask.sum() > 0:
            ssw_slp = ds[slp_var].isel({time_var: ssw_mask}).mean(dim=time_var).values
            ctrl_slp = ds[slp_var].isel({time_var: ctrl_mask}).mean(dim=time_var).values
            ssw_slp_maps.append(ssw_slp)
            ctrl_slp_maps.append(ctrl_slp)
        
        ds.close()
    
    t2m_anom = None
    slp_anom = None
    
    if ssw_t2m_maps:
        ssw_mean = np.nanmean(ssw_t2m_maps, axis=0)
        ctrl_mean = np.nanmean(ctrl_t2m_maps, axis=0)
        t2m_anom = ssw_mean - ctrl_mean
        print(f"T2m composite: {len(ssw_t2m_maps)} events, mean anomaly = {np.nanmean(t2m_anom):.2f} K")
    
    if ssw_slp_maps:
        ssw_mean = np.nanmean(ssw_slp_maps, axis=0)
        ctrl_mean = np.nanmean(ctrl_slp_maps, axis=0)
        slp_anom = (ssw_mean - ctrl_mean) / 100.0  # Convert Pa to hPa
        print(f"SLP composite: {len(ssw_slp_maps)} events, mean anomaly = {np.nanmean(slp_anom):.2f} hPa")
    
    return t2m_anom, slp_anom


def create_mechanism_timeline():
    """Create onset-aligned mechanism timeline from existing results."""
    
    # Load event-level dose-response data
    dose_df = pd.read_csv(RESULTS / "event_level_dose_response.csv")
    
    # Load event mechanism chain data
    import json
    with open(RESULTS / "r40b_event_mechanism.json") as f:
        mech = json.load(f)
    
    # Create timeline figure
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    fig.suptitle('SSW Onset-Aligned Mechanism Chain\n(Composite across 16 events)',
                 fontsize=14, fontweight='bold', y=0.98)
    
    # Panel A: Event-level RR vs Z500
    ax = axes[0]
    z500_vals = dose_df['z500_nh_m'].values
    rr_vals = dose_df['rr'].values
    colors = ['red' if rr > 1 else 'steelblue' for rr in rr_vals]
    ax.scatter(z500_vals, np.log(rr_vals), c=colors, s=60, edgecolors='black', linewidth=0.5, zorder=3)
    # Fit line
    from scipy import stats
    slope, intercept, r_val, p_val, _ = stats.linregress(z500_vals, np.log(rr_vals))
    x_fit = np.linspace(z500_vals.min(), z500_vals.max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, 'k--', alpha=0.7, linewidth=1.5)
    ax.axhline(0, color='grey', linestyle=':', alpha=0.5)
    ax.set_ylabel('log(RR)')
    ax.set_title(f'A  Z500 predicts avalanche response (r = {r_val:.2f}, P = {p_val:.3f})', 
                 fontsize=11, loc='left', fontweight='bold')
    ax.set_xlabel('NH Z500 anomaly (m)')
    
    # Panel B: Stratospheric intensity vs RR (null dose-response)
    ax = axes[1]
    t50_vals = dose_df['t50_K'].values
    ax.scatter(t50_vals, np.log(rr_vals), c=colors, s=60, edgecolors='black', linewidth=0.5, zorder=3)
    slope2, intercept2, r2, p2, _ = stats.linregress(t50_vals, np.log(rr_vals))
    x_fit2 = np.linspace(t50_vals.min(), t50_vals.max(), 100)
    ax.plot(x_fit2, slope2 * x_fit2 + intercept2, 'k--', alpha=0.7, linewidth=1.5)
    ax.axhline(0, color='grey', linestyle=':', alpha=0.5)
    ax.set_ylabel('log(RR)')
    ax.set_title(f'B  Stratospheric temperature does NOT predict response (r = {r2:.2f}, P = {p2:.3f})',
                 fontsize=11, loc='left', fontweight='bold')
    ax.set_xlabel('50 hPa temperature (K)')
    
    # Panel C: Trigger suppression chain
    ax = axes[2]
    events_chain = mech['event_chain']
    n_cold = sum(1 for e in events_chain if e['colder'])
    n_less_rain = sum(1 for e in events_chain if e['less_rain'])
    n_fewer_ros = sum(1 for e in events_chain if e['fewer_ros'])
    
    categories = ['Colder\n(T2m)', 'Less rain', 'Fewer\nROS events']
    counts = [n_cold, n_less_rain, n_fewer_ros]
    fractions = [c / len(events_chain) for c in counts]
    bar_colors = ['#4393c3', '#4393c3', '#4393c3']
    bars = ax.bar(categories, fractions, color=bar_colors, edgecolor='black', linewidth=0.5)
    ax.axhline(0.5, color='grey', linestyle=':', alpha=0.5, label='50% (chance)')
    ax.set_ylabel('Fraction of events')
    ax.set_ylim(0, 1)
    ax.set_title('C  Trigger suppression consistency across 16 events',
                 fontsize=11, loc='left', fontweight='bold')
    for bar, frac, count in zip(bars, fractions, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{count}/16\n({frac:.0%})', ha='center', va='bottom', fontsize=10)
    ax.legend(loc='upper right', fontsize=9)
    
    # Panel D: Threshold amplification model comparison
    ax = axes[3]
    with open(LOO_CV_FILE) as f:
        loo_data = json.load(f)
    
    k_vals = [p['k'] for p in loo_data['fit_profile']]
    d_vals = [p['d_hazard'] for p in loo_data['fit_profile']]
    rr_preds = [p['RR_pred'] for p in loo_data['fit_profile']]
    
    ax.plot(k_vals, d_vals, 'b-o', markersize=5, linewidth=1.5, label='Model $d_{hazard}$')
    ax.axhline(1.06, color='red', linestyle='--', linewidth=1.5, label='Observed $d_{hazard}$ = 1.06')
    ax.axhline(0.69, color='orange', linestyle=':', linewidth=1.5, label='$d_{weather}$ = 0.69 (ERA5)')
    
    # Shade LOO k range
    ax.axvspan(6, 8, alpha=0.15, color='green', label='LOO acceptable range')
    ax.set_xlabel('Number of trigger channels ($k$)')
    ax.set_ylabel('$d_{hazard}$')
    ax.set_title('D  Threshold amplification: modest weather shift → large hazard response',
                 fontsize=11, loc='left', fontweight='bold')
    ax.legend(loc='lower right', fontsize=9)
    ax.set_xlim(0.5, 15.5)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    outpath = FIGDIR / "fig_mechanism_chain.pdf"
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")
    return outpath


def create_dose_response_figure():
    """Create alternative dose-response figure showing Z500 works, SSW intensity doesn't."""
    dose_df = pd.read_csv(RESULTS / "event_level_dose_response.csv")
    from scipy import stats
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    
    rr_vals = dose_df['rr'].values
    log_rr = np.log(rr_vals)
    colors = ['#d73027' if rr > 1 else '#4575b4' for rr in rr_vals]
    
    # Panel A: Z500 vs log(RR) - the working predictor
    ax = axes[0]
    z500 = dose_df['z500_nh_m'].values
    ax.scatter(z500, log_rr, c=colors, s=70, edgecolors='black', linewidth=0.5, zorder=3)
    slope, intercept, r, p, _ = stats.linregress(z500, log_rr)
    x_fit = np.linspace(z500.min(), z500.max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, 'k-', linewidth=2, alpha=0.7)
    ax.axhline(0, color='grey', linestyle=':', alpha=0.5)
    ax.set_xlabel('Regional Alpine Z500 (m)', fontsize=11)
    ax.set_ylabel('log(Rate Ratio)', fontsize=11)
    ax.set_title(f'a  Tropospheric dose\nr = {r:.2f}, P = {p:.3f}, R² = {r**2:.2f}',
                 fontsize=11, fontweight='bold')
    
    # Panel B: SLP vs log(RR)
    ax = axes[1]
    slp = dose_df['slp_Pa'].values / 100  # Convert to hPa
    ax.scatter(slp, log_rr, c=colors, s=70, edgecolors='black', linewidth=0.5, zorder=3)
    slope2, intercept2, r2, p2, _ = stats.linregress(slp, log_rr)
    x_fit2 = np.linspace(slp.min(), slp.max(), 100)
    ax.plot(x_fit2, slope2 * x_fit2 + intercept2, 'k-', linewidth=2, alpha=0.7)
    ax.axhline(0, color='grey', linestyle=':', alpha=0.5)
    ax.set_xlabel('Mean SLP (hPa)', fontsize=11)
    ax.set_ylabel('log(Rate Ratio)', fontsize=11)
    ax.set_title(f'b  Surface pressure dose\nr = {r2:.2f}, P = {p2:.3f}, R² = {r2**2:.2f}',
                 fontsize=11, fontweight='bold')
    
    # Panel C: Stratospheric T50 vs log(RR) - the null
    ax = axes[2]
    t50 = dose_df['t50_K'].values
    ax.scatter(t50, log_rr, c=colors, s=70, edgecolors='black', linewidth=0.5, zorder=3)
    slope3, intercept3, r3, p3, _ = stats.linregress(t50, log_rr)
    x_fit3 = np.linspace(t50.min(), t50.max(), 100)
    ax.plot(x_fit3, slope3 * x_fit3 + intercept3, 'k-', linewidth=2, alpha=0.7)
    ax.axhline(0, color='grey', linestyle=':', alpha=0.5)
    ax.set_xlabel('50 hPa temperature (K)', fontsize=11)
    ax.set_ylabel('log(Rate Ratio)', fontsize=11)
    ax.set_title(f'c  Stratospheric dose (null)\nr = {r3:.2f}, P = {p3:.3f}, R² = {r3**2:.2f}',
                 fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    outpath = FIGDIR / "fig_dose_response_comparison.pdf"
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")
    return outpath


def create_spatial_composite():
    """Create ERA5 T2m and SLP spatial composite map."""
    t2m_anom, slp_anom = load_era5_composites()
    
    if t2m_anom is None and slp_anom is None:
        print("No ERA5 composites computed, creating placeholder")
        # Create a conceptual figure instead
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.text(0.5, 0.5, 'ERA5 composite\n(data loading issue)',
                ha='center', va='center', fontsize=14, transform=ax.transAxes)
        outpath = FIGDIR / "fig_spatial_composite.pdf"
        fig.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close()
        return outpath
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    if t2m_anom is not None:
        ax = axes[0]
        im = ax.imshow(t2m_anom, cmap='RdBu_r', aspect='auto',
                       extent=[5, 11, 44, 48], origin='lower')
        plt.colorbar(im, ax=ax, label='T2m anomaly (K)', shrink=0.8)
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.set_title('a  SSW composite T2m anomaly\n(8 events, ±15 d window)',
                     fontweight='bold', fontsize=11)
    
    if slp_anom is not None:
        ax = axes[1]
        im = ax.imshow(slp_anom, cmap='RdBu_r', aspect='auto',
                       extent=[5, 11, 44, 48], origin='lower')
        plt.colorbar(im, ax=ax, label='SLP anomaly (hPa)', shrink=0.8)
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.set_title('b  SSW composite SLP anomaly\n(8 events, ±15 d window)',
                     fontweight='bold', fontsize=11)
    
    plt.tight_layout()
    outpath = FIGDIR / "fig_spatial_composite.pdf"
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")
    return outpath


if __name__ == "__main__":
    print("=" * 60)
    print("R44: ERA5 Composites and Mechanism Timeline")
    print("=" * 60)
    
    # 1. Mechanism chain figure
    print("\n--- Creating mechanism chain figure ---")
    create_mechanism_timeline()
    
    # 2. Dose-response comparison figure
    print("\n--- Creating dose-response comparison figure ---")
    create_dose_response_figure()
    
    # 3. Spatial composites from ERA5
    print("\n--- Creating spatial composite maps ---")
    create_spatial_composite()
    
    print("\n=== All figures created ===")
