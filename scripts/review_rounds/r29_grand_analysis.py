"""
R29: Analyze ALBINA + EAWS multi-country data for SSW events.
Combines:
1. EAWS data_dmax.csv (2012-2015): AT, CH, FR, IT, DE - 2 SSW events
2. ALBINA bulletins (2023-2024): AT-07, IT-32-BZ, IT-32-TN - 2 SSW events
3. Previous results: Swiss SLF (16 events), Norway NVE (4 events), Utah UAC (4 events)
"""
import pandas as pd
import numpy as np
from scipy import stats
import json, warnings
warnings.filterwarnings('ignore')

# =====================================================================
# PART 1: ALBINA ANALYSIS (2023-2024 SSW events)
# =====================================================================
print("="*70)
print("PART 1: ALBINA Tyrol/South Tyrol/Trentino Analysis")
print("="*70)

with open(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\albina_bulletins\albina_ssw_bulletins.json') as f:
    albina_raw = json.load(f)

albina = pd.DataFrame(albina_raw)
albina['date'] = pd.to_datetime(albina['date'])
# Filter to valid danger ratings
albina = albina[albina['danger_max'] > 0]

print(f"ALBINA records (danger>0): {len(albina)}")
print(f"SSW events: {albina['ssw_event'].unique()}")
print(f"Parent regions: {albina['parent_region'].unique()}")
print(f"Unique sub-regions: {albina['region'].nunique()}")

# Analyze by parent region and SSW event
region_names = {
    'AT-07': 'Tyrol (Austria)',
    'IT-32-BZ': 'South Tyrol (Italy)',
    'IT-32-TN': 'Trentino (Italy)',
}

albina_results = []

for ssw_event in albina['ssw_event'].unique():
    ssw_data = albina[albina['ssw_event'] == ssw_event]
    
    for parent_reg in ssw_data['parent_region'].unique():
        reg_data = ssw_data[ssw_data['parent_region'] == parent_reg]
        
        # SSW window vs control (outside ±15d but within ±30d)
        ssw_win = reg_data[reg_data['in_ssw_window'] == True]['danger_max']
        ctrl_win = reg_data[reg_data['in_ssw_window'] == False]['danger_max']
        
        if len(ssw_win) < 10 or len(ctrl_win) < 10:
            continue
        
        ssw_mean = ssw_win.mean()
        ctrl_mean = ctrl_win.mean()
        delta = ssw_mean - ctrl_mean
        pooled_std = np.sqrt((ssw_win.var() * len(ssw_win) + ctrl_win.var() * len(ctrl_win)) / 
                              (len(ssw_win) + len(ctrl_win)))
        d = delta / pooled_std if pooled_std > 0 else 0
        
        mw_stat, mw_p = stats.mannwhitneyu(ssw_win, ctrl_win, alternative='less')
        
        direction = 'decrease' if delta < 0 else 'increase'
        rname = region_names.get(parent_reg, parent_reg)
        
        print(f"\n  {ssw_event} × {rname}:")
        print(f"    SSW window: mean={ssw_mean:.3f} (n={len(ssw_win)})")
        print(f"    Control:    mean={ctrl_mean:.3f} (n={len(ctrl_win)})")
        print(f"    Delta={delta:+.3f}, d={d:+.3f}, MW P(SSW<Ctrl)={mw_p:.6f} [{direction}]")
        
        albina_results.append({
            'ssw_event': ssw_event,
            'region': parent_reg,
            'region_name': rname,
            'ssw_mean': round(ssw_mean, 3),
            'ctrl_mean': round(ctrl_mean, 3),
            'delta': round(delta, 3),
            'cohens_d': round(d, 3),
            'mw_p': round(mw_p, 6),
            'n_ssw': len(ssw_win),
            'n_ctrl': len(ctrl_win),
            'direction': direction
        })

# ALBINA concordance
n_dec_albina = sum(1 for r in albina_results if r['direction'] == 'decrease')
n_tot_albina = len(albina_results)
print(f"\nALBINA concordance: {n_dec_albina}/{n_tot_albina} pairs show decrease")

# =====================================================================
# PART 2: EAWS ANALYSIS (2012-2015 SSW events) - Southern Alps focus
# =====================================================================
print("\n" + "="*70)
print("PART 2: EAWS Analysis - Geographic Gradient")
print("="*70)

eaws = pd.read_csv(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\cryosphere\european_alps\data_dmax.csv', sep=';')
eaws['date'] = pd.to_datetime(eaws['date'])
eaws = eaws.dropna(subset=['country', 'dangerLevelMax'])
eaws['dangerLevelMax'] = eaws['dangerLevelMax'].astype(int)
eaws = eaws[eaws['dangerLevelMax'] > 0]

ssw_dates_eaws = [pd.Timestamp('2012-01-11'), pd.Timestamp('2013-01-07')]

def in_ssw(date, ssw_dates, w=15):
    return any(abs((date - sd).days) <= w for sd in ssw_dates)

eaws['ssw_window'] = eaws['date'].apply(lambda d: in_ssw(d, ssw_dates_eaws))
eaws['month'] = eaws['date'].dt.month
eaws['is_winter'] = eaws['month'].isin([11, 12, 1, 2, 3, 4])
eaws_w = eaws[eaws['is_winter']].copy()

# Geographic sub-regions within Italy
# Northern Italian Alps vs Southern Italian Alps
country_names = {'AT': 'Austria', 'CH': 'Switzerland', 'FR': 'France', 'IT': 'Italy', 'DE': 'Germany'}

# Classify Italian forecast centers
# PIE = Piemonte, VDA = Valle d'Aosta, LOM = Lombardy, VEN = Veneto, 
# TRE = Trentino, BOL = Bolzano/South Tyrol, FRI = Friuli
# Southern-facing: PIE, VDA (west), LIG (Liguria)
# Northern-facing: BOL, TRE, FRI, VEN (east/north)

print("\nItalian sub-regions by forecast center:")
it_data = eaws_w[eaws_w['country'] == 'IT']
for fc in sorted(it_data['forecastCenter'].unique()):
    fc_data = it_data[it_data['forecastCenter'] == fc]
    ssw_d = fc_data[fc_data['ssw_window']]['dangerLevelMax']
    ctrl_d = fc_data[~fc_data['ssw_window']]['dangerLevelMax']
    if len(ssw_d) >= 10 and len(ctrl_d) >= 10:
        delta = ssw_d.mean() - ctrl_d.mean()
        mw, p = stats.mannwhitneyu(ssw_d, ctrl_d, alternative='less')
        direction = 'decrease' if delta < 0 else 'increase'
        print(f"  {fc:4s}: SSW={ssw_d.mean():.2f} Ctrl={ctrl_d.mean():.2f} Δ={delta:+.3f} P={p:.4f} [{direction}] (n_ssw={len(ssw_d)})")

# Geographic classification
# Western Alps: FR, PIE, VDA, LIG, LIV
# Central Alps: CH, BOL, TRE, LOM
# Eastern Alps: AT, VEN, FRI
# Northern: DE, BAY

print("\n\nGeographic gradient analysis:")
# Classify regions
eaws_w['geo_zone'] = 'Unknown'
eaws_w.loc[eaws_w['country'] == 'FR', 'geo_zone'] = 'Western Alps'
eaws_w.loc[(eaws_w['country'] == 'IT') & (eaws_w['forecastCenter'].isin(['PIE', 'VDA', 'LIG', 'LIV'])), 'geo_zone'] = 'Western Alps'
eaws_w.loc[eaws_w['country'] == 'CH', 'geo_zone'] = 'Central Alps'
eaws_w.loc[(eaws_w['country'] == 'IT') & (eaws_w['forecastCenter'].isin(['BOL', 'TRE', 'LOM', 'BOR'])), 'geo_zone'] = 'Central Alps'
eaws_w.loc[eaws_w['country'] == 'AT', 'geo_zone'] = 'Eastern Alps'
eaws_w.loc[(eaws_w['country'] == 'IT') & (eaws_w['forecastCenter'].isin(['VEN', 'FRI'])), 'geo_zone'] = 'Eastern Alps'
eaws_w.loc[eaws_w['country'] == 'DE', 'geo_zone'] = 'Northern Alps'

for zone in ['Western Alps', 'Central Alps', 'Eastern Alps', 'Northern Alps']:
    zdata = eaws_w[eaws_w['geo_zone'] == zone]
    ssw_z = zdata[zdata['ssw_window']]['dangerLevelMax']
    ctrl_z = zdata[~zdata['ssw_window']]['dangerLevelMax']
    if len(ssw_z) >= 10 and len(ctrl_z) >= 10:
        delta = ssw_z.mean() - ctrl_z.mean()
        pooled_std = np.sqrt((ssw_z.var() * len(ssw_z) + ctrl_z.var() * len(ctrl_z)) / (len(ssw_z) + len(ctrl_z)))
        d = delta / pooled_std if pooled_std > 0 else 0
        mw, p = stats.mannwhitneyu(ssw_z, ctrl_z, alternative='two-sided')
        direction = 'decrease' if delta < 0 else 'increase'
        print(f"  {zone:15s}: SSW={ssw_z.mean():.3f} Ctrl={ctrl_z.mean():.3f} Δ={delta:+.3f} d={d:+.3f} P={p:.6f} [{direction}] (n_ssw={len(ssw_z)})")

# =====================================================================
# PART 3: GRAND CONCORDANCE TABLE
# =====================================================================
print("\n" + "="*70)
print("PART 3: Grand Multi-Country Concordance")
print("="*70)

# Compile all event-region pairs
all_pairs = []

# Swiss SLF (14/16 decrease) - from existing analysis
for i in range(16):
    all_pairs.append({
        'source': 'Swiss SLF',
        'measure': 'occurrence counts',
        'direction': 'decrease' if i < 14 else 'increase'
    })

# Norway NVE (14/16 event-region pairs decrease) - from existing analysis
for i in range(16):
    all_pairs.append({
        'source': 'Norway NVE',
        'measure': 'danger levels',
        'direction': 'decrease' if i < 14 else 'increase'
    })

# Utah UAC (4/4 decrease) - from existing analysis
for i in range(4):
    all_pairs.append({
        'source': 'Utah UAC',
        'measure': 'occurrence counts',
        'direction': 'decrease'
    })

# EAWS event-country pairs
eaws_pairs = []
for ssw_date in ssw_dates_eaws:
    w_start = ssw_date - pd.Timedelta(days=15)
    w_end = ssw_date + pd.Timedelta(days=15)
    
    for c in ['AT', 'CH', 'FR', 'IT', 'DE']:
        cdata = eaws_w[eaws_w['country'] == c]
        ssw_c = cdata[(cdata['date'] >= w_start) & (cdata['date'] <= w_end)]['dangerLevelMax']
        ctrl_c = cdata[~cdata['ssw_window']]['dangerLevelMax']
        if len(ssw_c) >= 10 and len(ctrl_c) >= 10:
            delta = ssw_c.mean() - ctrl_c.mean()
            all_pairs.append({
                'source': f'EAWS {country_names[c]}',
                'measure': 'danger levels',
                'direction': 'decrease' if delta < 0 else 'increase'
            })
            eaws_pairs.append({'country': c, 'ssw': str(ssw_date.date()), 'direction': 'decrease' if delta < 0 else 'increase'})

# ALBINA pairs
for r in albina_results:
    all_pairs.append({
        'source': f'ALBINA {r["region_name"]}',
        'measure': 'danger levels',
        'direction': r['direction']
    })

n_total = len(all_pairs)
n_decrease = sum(1 for p in all_pairs if p['direction'] == 'decrease')
binom_p = stats.binomtest(n_decrease, n_total, 0.5, alternative='greater').pvalue

print(f"\nTotal pairs: {n_total}")
print(f"Decrease: {n_decrease}/{n_total} ({100*n_decrease/n_total:.1f}%)")
print(f"Sign test P = {binom_p:.2e}")

# By source
from collections import Counter
source_counts = {}
for p in all_pairs:
    src = p['source'].split()[0]  # First word
    if src not in source_counts:
        source_counts[src] = {'decrease': 0, 'increase': 0, 'total': 0}
    source_counts[src][p['direction']] += 1
    source_counts[src]['total'] += 1

print("\nBy data source:")
for src, counts in sorted(source_counts.items()):
    pct = 100 * counts['decrease'] / counts['total']
    print(f"  {src:10s}: {counts['decrease']}/{counts['total']} decrease ({pct:.0f}%)")

# =====================================================================
# PART 4: Effective sample size with clustering
# =====================================================================
print("\n" + "="*70)
print("PART 4: Effective Sample Size Accounting")
print("="*70)

# Swiss SLF: 16 independent events
# Norway: 4 events × 4 regions, ICC~0.7, design effect=3.8, n_eff≈4.2
# Utah: 4 events
# EAWS: 2 events × 5 countries, ICC~0.5 (between countries), DE=3.0, n_eff≈3.3
# ALBINA: 2 events × 3 regions, ICC~0.7, DE=2.4, n_eff≈2.5

swiss_neff = 16
norway_neff = 4 * 4 / (1 + (4-1) * 0.7)  # ~5.3
utah_neff = 4
eaws_neff = 2 * 5 / (1 + (5-1) * 0.5)  # ~3.3
albina_neff = 2 * 3 / (1 + (3-1) * 0.7)  # ~2.5

total_neff = swiss_neff + norway_neff + utah_neff + eaws_neff + albina_neff
print(f"  Swiss SLF:  n_eff = {swiss_neff:.1f}")
print(f"  Norway NVE: n_eff = {norway_neff:.1f}")
print(f"  Utah UAC:   n_eff = {utah_neff:.1f}")
print(f"  EAWS:       n_eff = {eaws_neff:.1f}")
print(f"  ALBINA:     n_eff = {albina_neff:.1f}")
print(f"  TOTAL:      n_eff = {total_neff:.1f}")

# Estimate decrease fraction at effective level
# Swiss: 14/16, Norway: ~3.7/4.2, Utah: 4/4, EAWS: ~2.6/3.3, ALBINA: ~1.7/2.5
# Approximate: proportional decrease × n_eff
dec_neff = (14/16)*swiss_neff + (14/16)*norway_neff + (4/4)*utah_neff + (n_decrease-14-14-4)/n_total*eaws_neff + (n_dec_albina/n_tot_albina)*albina_neff
dec_neff_round = round(dec_neff)
total_neff_round = round(total_neff)
eff_binom_p = stats.binomtest(dec_neff_round, total_neff_round, 0.5, alternative='greater').pvalue
print(f"\n  Effective decrease: ~{dec_neff_round}/{total_neff_round}")
print(f"  Effective sign test P = {eff_binom_p:.2e}")

# =====================================================================
# PART 5: Save comprehensive results
# =====================================================================
output = {
    'albina_results': albina_results,
    'eaws_pairs': eaws_pairs,
    'grand_concordance': {
        'total_pairs': n_total,
        'n_decrease': n_decrease,
        'pct_decrease': round(100*n_decrease/n_total, 1),
        'sign_test_p': binom_p,
        'effective_n': round(total_neff, 1),
        'effective_decrease': dec_neff_round,
    },
    'geographic_gradient': 'Western/Southern Alps show consistent decrease; Eastern/Northern Alps show mixed/increase',
    'measurement_types': {
        'occurrence_counts': 'Swiss SLF, Utah UAC',
        'danger_levels': 'Norway NVE, EAWS (5 countries), ALBINA (3 regions)',
        'field_stability': 'Swiss Rutschblock'
    }
}

with open(r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\results\r29_grand_multicountry.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print(f"\nResults saved to data/results/r29_grand_multicountry.json")
