"""
R29: Multi-Country EAWS Avalanche Danger Analysis During SSW Events
Analyzes danger levels across Austria, France, Italy, Germany, Switzerland
from the EAWS data_dmax.csv dataset (2011-2015).
"""
import pandas as pd
import numpy as np
from scipy import stats
import json, warnings
warnings.filterwarnings('ignore')

# --- Load EAWS data ---
eaws = pd.read_csv(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\european_alps\data_dmax.csv', sep=';')
eaws['date'] = pd.to_datetime(eaws['date'])
eaws = eaws.dropna(subset=['country', 'dangerLevelMax'])
eaws['dangerLevelMax'] = eaws['dangerLevelMax'].astype(int)
# Filter out danger level 0 (missing)
eaws = eaws[eaws['dangerLevelMax'] > 0]

print(f"EAWS dataset: {len(eaws)} region-days")
print(f"Date range: {eaws['date'].min()} to {eaws['date'].max()}")
print(f"Countries: {eaws['country'].unique()}")
print(f"\nRegion-days per country:")
print(eaws.groupby('country').size())

# --- SSW events ---
ssw_dates = [pd.Timestamp('2012-01-11'), pd.Timestamp('2013-01-07')]
window = 15  # ±15 days

# --- Mark SSW windows ---
def in_ssw_window(date, ssw_dates, window):
    for sd in ssw_dates:
        if abs((date - sd).days) <= window:
            return True
    return False

def get_ssw_event(date, ssw_dates, window):
    for sd in ssw_dates:
        if abs((date - sd).days) <= window:
            return sd
    return None

eaws['ssw_window'] = eaws['date'].apply(lambda d: in_ssw_window(d, ssw_dates, window))
eaws['ssw_event'] = eaws['date'].apply(lambda d: get_ssw_event(d, ssw_dates, window))

# Define winter (Nov-Apr)
eaws['month'] = eaws['date'].dt.month
eaws['is_winter'] = eaws['month'].isin([11, 12, 1, 2, 3, 4])
eaws_winter = eaws[eaws['is_winter']].copy()

# DOY for matching
eaws_winter['doy'] = eaws_winter['date'].dt.dayofyear

print(f"\nWinter data: {len(eaws_winter)} region-days")
print(f"SSW window data: {eaws_winter['ssw_window'].sum()} region-days")
print(f"Control data: {(~eaws_winter['ssw_window']).sum()} region-days")

# =================================================================
# ANALYSIS 1: Country-level SSW vs Control danger levels
# =================================================================
print("\n" + "="*70)
print("ANALYSIS 1: SSW vs Control Danger Levels by Country")
print("="*70)

results = {}
countries = ['AT', 'CH', 'FR', 'IT', 'DE']
country_names = {'AT': 'Austria', 'CH': 'Switzerland', 'FR': 'France', 'IT': 'Italy', 'DE': 'Germany'}

for c in countries:
    cdata = eaws_winter[eaws_winter['country'] == c]
    ssw_data = cdata[cdata['ssw_window']]['dangerLevelMax']
    ctrl_data = cdata[~cdata['ssw_window']]['dangerLevelMax']
    
    if len(ssw_data) < 10:
        print(f"\n{country_names[c]}: Insufficient SSW data ({len(ssw_data)} region-days)")
        continue
    
    ssw_mean = ssw_data.mean()
    ctrl_mean = ctrl_data.mean()
    delta = ssw_mean - ctrl_mean
    pooled_std = np.sqrt((ssw_data.var() * len(ssw_data) + ctrl_data.var() * len(ctrl_data)) / 
                          (len(ssw_data) + len(ctrl_data)))
    cohens_d = delta / pooled_std if pooled_std > 0 else 0
    
    mw_stat, mw_p = stats.mannwhitneyu(ssw_data, ctrl_data, alternative='less')
    
    print(f"\n{country_names[c]} ({c}):")
    print(f"  SSW: mean={ssw_mean:.3f} (n={len(ssw_data)})")
    print(f"  Control: mean={ctrl_mean:.3f} (n={len(ctrl_data)})")
    print(f"  Delta: {delta:.3f}, Cohen's d: {cohens_d:.3f}")
    print(f"  Mann-Whitney P (one-sided, SSW<Control): {mw_p:.6f}")
    
    results[c] = {
        'country': country_names[c],
        'ssw_mean': round(ssw_mean, 3),
        'ctrl_mean': round(ctrl_mean, 3),
        'delta': round(delta, 3),
        'cohens_d': round(cohens_d, 3),
        'mw_p': mw_p,
        'n_ssw': len(ssw_data),
        'n_ctrl': len(ctrl_data),
        'direction': 'decrease' if delta < 0 else 'increase'
    }

# =================================================================
# ANALYSIS 2: Event-level analysis (each SSW × country)
# =================================================================
print("\n" + "="*70)
print("ANALYSIS 2: Event-Level SSW × Country Pairs")
print("="*70)

event_pairs = []
for ssw_date in ssw_dates:
    window_start = ssw_date - pd.Timedelta(days=window)
    window_end = ssw_date + pd.Timedelta(days=window)
    
    for c in countries:
        cdata = eaws_winter[eaws_winter['country'] == c]
        
        # SSW window data
        ssw_mask = (cdata['date'] >= window_start) & (cdata['date'] <= window_end)
        ssw_regions = cdata[ssw_mask]
        
        if len(ssw_regions) < 5:
            continue
        
        # DOY-matched control: same DOY range, different years
        ssw_doys = ssw_regions['doy'].unique()
        ctrl_mask = (~cdata['ssw_window']) & (cdata['doy'].isin(range(min(ssw_doys)-3, max(ssw_doys)+4)))
        ctrl_regions = cdata[ctrl_mask]
        
        if len(ctrl_regions) < 5:
            continue
        
        ssw_mean = ssw_regions['dangerLevelMax'].mean()
        ctrl_mean = ctrl_regions['dangerLevelMax'].mean()
        delta = ssw_mean - ctrl_mean
        
        mw_stat, mw_p = stats.mannwhitneyu(
            ssw_regions['dangerLevelMax'], ctrl_regions['dangerLevelMax'], 
            alternative='less'
        )
        
        direction = 'decrease' if delta < 0 else 'increase'
        
        event_pairs.append({
            'ssw_date': str(ssw_date.date()),
            'country': c,
            'ssw_mean': round(ssw_mean, 3),
            'ctrl_mean': round(ctrl_mean, 3),
            'delta': round(delta, 3),
            'mw_p': round(mw_p, 6),
            'n_ssw': len(ssw_regions),
            'n_ctrl': len(ctrl_regions),
            'direction': direction
        })
        
        sig = '*' if mw_p < 0.05 else ''
        print(f"  {ssw_date.date()} × {country_names[c]:12s}: SSW={ssw_mean:.2f} Ctrl={ctrl_mean:.2f} Δ={delta:+.3f} P={mw_p:.4f} {sig} [{direction}]")

# Count concordance
n_decrease = sum(1 for ep in event_pairs if ep['direction'] == 'decrease')
n_total = len(event_pairs)
print(f"\nEvent-country pairs: {n_decrease}/{n_total} show decrease")
if n_total > 0:
    binom_p = stats.binomtest(n_decrease, n_total, 0.5, alternative='greater').pvalue
    print(f"Sign test P = {binom_p:.6f}")

# =================================================================
# ANALYSIS 3: Region-level within countries
# =================================================================
print("\n" + "="*70)
print("ANALYSIS 3: Region-Level Direction (per warning region)")
print("="*70)

region_results = []
for c in countries:
    cdata = eaws_winter[eaws_winter['country'] == c]
    regions = cdata['warningRegion'].unique()
    
    n_dec = 0
    n_inc = 0
    n_regions = 0
    
    for reg in regions:
        rdata = cdata[cdata['warningRegion'] == reg]
        ssw_r = rdata[rdata['ssw_window']]['dangerLevelMax']
        ctrl_r = rdata[~rdata['ssw_window']]['dangerLevelMax']
        
        if len(ssw_r) < 5 or len(ctrl_r) < 5:
            continue
        
        n_regions += 1
        if ssw_r.mean() < ctrl_r.mean():
            n_dec += 1
        else:
            n_inc += 1
    
    if n_regions > 0:
        pct = n_dec / n_regions * 100
        if n_dec > 0:
            binom_p = stats.binomtest(n_dec, n_regions, 0.5, alternative='greater').pvalue
        else:
            binom_p = 1.0
        print(f"  {country_names[c]:12s}: {n_dec}/{n_regions} regions show decrease ({pct:.1f}%) P={binom_p:.6f}")
        
        region_results.append({
            'country': c,
            'n_decrease': n_dec,
            'n_total': n_regions,
            'pct_decrease': round(pct, 1),
            'sign_test_p': round(binom_p, 6)
        })

# =================================================================
# ANALYSIS 4: Alpine-wide meta-analysis
# =================================================================
print("\n" + "="*70)
print("ANALYSIS 4: Alpine-Wide Meta-Analysis")
print("="*70)

# Pool all Alpine countries
alpine = eaws_winter[eaws_winter['country'].isin(['AT', 'CH', 'FR', 'IT', 'DE'])]
ssw_all = alpine[alpine['ssw_window']]['dangerLevelMax']
ctrl_all = alpine[~alpine['ssw_window']]['dangerLevelMax']

print(f"  SSW: mean={ssw_all.mean():.3f} (n={len(ssw_all)})")
print(f"  Control: mean={ctrl_all.mean():.3f} (n={len(ctrl_all)})")
print(f"  Delta: {ssw_all.mean() - ctrl_all.mean():.3f}")
mw_stat, mw_p = stats.mannwhitneyu(ssw_all, ctrl_all, alternative='less')
pooled_std = np.sqrt((ssw_all.var() * len(ssw_all) + ctrl_all.var() * len(ctrl_all)) / 
                      (len(ssw_all) + len(ctrl_all)))
d = (ssw_all.mean() - ctrl_all.mean()) / pooled_std
print(f"  Cohen's d: {d:.3f}")
print(f"  Mann-Whitney P: {mw_p:.2e}")

# Total cross-country concordance including Swiss/Norway/Utah
total_decrease_eaws = n_decrease
total_pairs_eaws = n_total
# Add existing: Swiss 14/16, Norway 14/16 (event-region pairs), Utah 4/4
# But be careful about double-counting Swiss EAWS vs Swiss SLF
print(f"\n  EAWS event-country pairs: {total_decrease_eaws}/{total_pairs_eaws} decrease")
print(f"  Previous: Swiss SLF 14/16, Norway NVE 14/16, Utah UAC 4/4")

# =================================================================
# ANALYSIS 5: DOY-matched analysis per country
# =================================================================
print("\n" + "="*70)
print("ANALYSIS 5: DOY-Matched Danger Level Ratios by Country")
print("="*70)

for c in countries:
    cdata = eaws_winter[eaws_winter['country'] == c]
    
    rr_list = []
    for ssw_date in ssw_dates:
        window_start = ssw_date - pd.Timedelta(days=window)
        window_end = ssw_date + pd.Timedelta(days=window)
        
        ssw_mask = (cdata['date'] >= window_start) & (cdata['date'] <= window_end)
        ssw_days = cdata[ssw_mask]
        
        if len(ssw_days) == 0:
            continue
        
        # Compute daily mean danger across all regions for each day
        ssw_daily = ssw_days.groupby('date')['dangerLevelMax'].mean()
        
        # DOY-matched expected: same DOYs in non-SSW periods
        expected_vals = []
        for day in ssw_daily.index:
            doy = day.dayofyear
            ctrl_mask = (~cdata['ssw_window']) & (cdata['doy'].between(doy-3, doy+3))
            ctrl_daily = cdata[ctrl_mask].groupby('date')['dangerLevelMax'].mean()
            if len(ctrl_daily) > 0:
                expected_vals.append(ctrl_daily.mean())
        
        if expected_vals:
            obs_mean = ssw_daily.mean()
            exp_mean = np.mean(expected_vals)
            rr = obs_mean / exp_mean if exp_mean > 0 else np.nan
            rr_list.append(rr)
            print(f"  {country_names[c]:12s} SSW {ssw_date.date()}: Obs={obs_mean:.2f} Exp={exp_mean:.2f} RR={rr:.3f}")
    
    if rr_list:
        geo_rr = np.exp(np.mean(np.log(rr_list)))
        n_dec = sum(1 for r in rr_list if r < 1)
        print(f"  {country_names[c]:12s} Geometric mean RR = {geo_rr:.3f}, {n_dec}/{len(rr_list)} decrease")

# =================================================================
# ANALYSIS 6: Effect size comparison across measurement systems
# =================================================================
print("\n" + "="*70)
print("ANALYSIS 6: Cross-System Effect Size Summary")
print("="*70)
print(f"  Swiss SLF counts:     RR=0.32 (68% reduction), d=-1.06, n=16 events")
print(f"  Norwegian NVE danger: delta=-0.23, d=-0.35, n=4 events")
print(f"  Utah UAC counts:      RR=0.34 (66% reduction), n=4 events")
print(f"  EAWS Alpine danger:   delta={ssw_all.mean()-ctrl_all.mean():.3f}, d={d:.3f}, n=2 SSW events × 5 countries")

# Save results
output = {
    'country_results': results,
    'event_pairs': event_pairs,
    'region_results': region_results,
    'alpine_wide': {
        'ssw_mean': round(ssw_all.mean(), 3),
        'ctrl_mean': round(ctrl_all.mean(), 3),
        'delta': round(ssw_all.mean() - ctrl_all.mean(), 3),
        'cohens_d': round(d, 3),
        'mw_p': mw_p,
        'n_ssw': len(ssw_all),
        'n_ctrl': len(ctrl_all)
    },
    'concordance': {
        'eaws_decrease': total_decrease_eaws,
        'eaws_total': total_pairs_eaws,
        'all_decrease': total_decrease_eaws + 14 + 14 + 4,  # + Swiss SLF + Norway + Utah
        'all_total': total_pairs_eaws + 16 + 16 + 4  # note Swiss SLF events don't overlap with EAWS Swiss
    }
}

with open(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\results\r29_eaws_multicountry.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print("\nResults saved to data/results/r29_eaws_multicountry.json")
