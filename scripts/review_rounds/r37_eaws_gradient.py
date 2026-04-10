"""
EAWS multi-country spatial gradient analysis for SSW events.
Tests whether avalanche danger reduction during SSW varies by geography.
Uses 5-country EAWS data (2011-2015), 2 SSW events (Jan 2012, Jan 2013).
"""
import pandas as pd, numpy as np, json
from scipy import stats

eaws = pd.read_csv('data/cryosphere/european_alps/data_dmax.csv', sep=';', low_memory=False)
eaws['date'] = pd.to_datetime(eaws['date'])
eaws = eaws.dropna(subset=['dangerLevelMax', 'forecastCenter'])

# Approximate centre coordinates (lat, lon) for spatial regression
centre_coords = {
    'PIE': (45.0, 7.5, 'IT'),   # Piemonte
    'BOL': (46.7, 11.3, 'IT'),  # Bolzano/South Tyrol
    'SWI': (46.8, 8.2, 'CH'),   # Switzerland
    'LOM': (46.0, 10.0, 'IT'),  # Lombardia
    'VOR': (47.2, 10.0, 'AT'),  # Vorarlberg
    'KAE': (46.7, 13.8, 'AT'),  # Kärnten
    'VEN': (46.2, 11.8, 'IT'),  # Veneto
    'FRI': (46.4, 13.0, 'IT'),  # Friuli
    'VDA': (45.7, 7.3, 'IT'),   # Valle d'Aosta
    'CHX': (45.9, 6.9, 'FR'),   # Chamonix
    'BSM': (45.5, 6.5, 'FR'),   # Bourg-Saint-Maurice
    'BAY': (47.5, 11.0, 'DE'),  # Bayern
    'GRE': (45.2, 5.7, 'FR'),   # Grenoble/Isère
    'BRI': (44.9, 6.6, 'FR'),   # Briançon
    'LIV': (46.5, 10.1, 'IT'),  # Livigno
    'LIG': (44.3, 8.5, 'IT'),   # Liguria
    'BOR': (45.0, 6.0, 'FR'),   # Bourg d'Oisans
    'NIE': (47.5, 15.8, 'AT'),  # Niederösterreich
    'OBE': (47.8, 14.0, 'AT'),  # Oberösterreich
    'TIR': (47.0, 11.4, 'AT'),  # Tirol
    'SAL': (47.2, 13.0, 'AT'),  # Salzburg
    'STE': (47.3, 14.8, 'AT'),  # Steiermark
    'TRE': (46.1, 11.1, 'IT'),  # Trentino
}

ssw_events = [pd.Timestamp('2012-01-11'), pd.Timestamp('2013-01-07')]

results_by_centre = []

for centre in eaws['forecastCenter'].unique():
    if centre not in centre_coords:
        continue
    
    lat, lon, country = centre_coords[centre]
    cdata = eaws[eaws['forecastCenter'] == centre]
    
    centre_results = []
    for ssw_date in ssw_events:
        # SSW window
        w = (cdata['date'] >= ssw_date - pd.Timedelta(days=15)) & \
            (cdata['date'] <= ssw_date + pd.Timedelta(days=15))
        ssw_danger = cdata.loc[w, 'dangerLevelMax'].mean()
        
        # DOY-matched control (same DOY range, other years)
        doy_start = (ssw_date - pd.Timedelta(days=15)).dayofyear
        doy_end = (ssw_date + pd.Timedelta(days=15)).dayofyear
        if doy_start <= doy_end:
            doy_match = (cdata['date'].dt.dayofyear >= doy_start) & \
                        (cdata['date'].dt.dayofyear <= doy_end)
        else:  # year boundary crossing
            doy_match = (cdata['date'].dt.dayofyear >= doy_start) | \
                        (cdata['date'].dt.dayofyear <= doy_end)
        ctrl_mask = doy_match & ~w
        ctrl_danger = cdata.loc[ctrl_mask, 'dangerLevelMax'].mean()
        
        if not np.isnan(ssw_danger) and not np.isnan(ctrl_danger) and ctrl_danger > 0:
            rr = ssw_danger / ctrl_danger
            anomaly = ssw_danger - ctrl_danger
            centre_results.append({
                'ssw_date': ssw_date,
                'ssw_danger': ssw_danger,
                'ctrl_danger': ctrl_danger,
                'rr': rr,
                'anomaly': anomaly,
            })
    
    if len(centre_results) > 0:
        mean_rr = np.mean([r['rr'] for r in centre_results])
        mean_anom = np.mean([r['anomaly'] for r in centre_results])
        n_decrease = sum(1 for r in centre_results if r['rr'] < 1)
        
        results_by_centre.append({
            'centre': centre,
            'country': country,
            'lat': lat,
            'lon': lon,
            'mean_rr': mean_rr,
            'mean_anomaly': mean_anom,
            'n_events': len(centre_results),
            'n_decrease': n_decrease,
        })

rdf = pd.DataFrame(results_by_centre)

print(f"{'='*70}")
print(f"EAWS SPATIAL GRADIENT ANALYSIS")
print(f"{'='*70}")
print(f"Centres analyzed: {len(rdf)}")
print(f"Countries: {rdf['country'].nunique()}")
print(f"\nCentre-level results (sorted by latitude):")
print(f"{'Centre':>6} {'Country':>3} {'Lat':>6} {'Lon':>6} {'RR':>6} {'Anom':>6} {'Dir':>4}")
print("-" * 50)

for _, r in rdf.sort_values('lat').iterrows():
    d = "↓" if r['mean_rr'] < 1 else "↑"
    print(f"{r['centre']:>6} {r['country']:>3} {r['lat']:>6.1f} {r['lon']:>6.1f} "
          f"{r['mean_rr']:>6.3f} {r['mean_anomaly']:>+6.3f} {d:>4}")

n_dec = (rdf['mean_rr'] < 1).sum()
n_total = len(rdf)
print(f"\nOverall: {n_dec}/{n_total} centres show decrease ({n_dec/n_total*100:.1f}%)")
print(f"Sign test vs 50%: P = {stats.binomtest(n_dec, n_total, 0.5).pvalue:.4f}")

# Spatial regression: RR vs latitude
r_lat, p_lat = stats.pearsonr(rdf['lat'], rdf['mean_rr'])
print(f"\n--- Latitude Gradient ---")
print(f"  Pearson r(lat, RR) = {r_lat:.3f}, P = {p_lat:.4f}")
print(f"  Interpretation: {'Higher lat = MORE suppression' if r_lat < 0 else 'Higher lat = LESS suppression'}")

# RR vs longitude (west-east)
r_lon, p_lon = stats.pearsonr(rdf['lon'], rdf['mean_rr'])
print(f"\n--- Longitude Gradient ---")
print(f"  Pearson r(lon, RR) = {r_lon:.3f}, P = {p_lon:.4f}")
print(f"  Interpretation: {'Eastern = MORE suppression' if r_lon < 0 else 'Eastern = LESS suppression'}")

# Mediterranean proximity (distance from ~43°N, 8°E)
med_lat, med_lon = 43.0, 8.0
rdf['med_dist'] = np.sqrt((rdf['lat'] - med_lat)**2 + (rdf['lon'] - med_lon)**2)
r_med, p_med = stats.pearsonr(rdf['med_dist'], rdf['mean_rr'])
print(f"\n--- Mediterranean Distance Gradient ---")
print(f"  Pearson r(dist, RR) = {r_med:.3f}, P = {p_med:.4f}")
print(f"  Interpretation: {'Closer to Med = MORE suppression' if r_med > 0 else 'Closer to Med = LESS suppression'}")

# Country-level summary
print(f"\n--- Country-Level Summary ---")
for country in ['FR', 'CH', 'IT', 'AT', 'DE']:
    cc = rdf[rdf['country'] == country]
    if len(cc) > 0:
        n_d = (cc['mean_rr'] < 1).sum()
        print(f"  {country}: {n_d}/{len(cc)} decrease, mean RR = {cc['mean_rr'].mean():.3f}")

# SSW-type stratification
print(f"\n--- SSW-Type Stratification ---")
for ssw_date, ssw_type in [(pd.Timestamp('2012-01-11'), 'displacement'),
                            (pd.Timestamp('2013-01-07'), 'displacement')]:
    event_results = []
    for centre in eaws['forecastCenter'].unique():
        if centre not in centre_coords:
            continue
        cdata = eaws[eaws['forecastCenter'] == centre]
        w = (cdata['date'] >= ssw_date - pd.Timedelta(days=15)) & \
            (cdata['date'] <= ssw_date + pd.Timedelta(days=15))
        ssw_danger = cdata.loc[w, 'dangerLevelMax'].mean()
        doy_start = (ssw_date - pd.Timedelta(days=15)).dayofyear
        doy_end = (ssw_date + pd.Timedelta(days=15)).dayofyear
        if doy_start <= doy_end:
            doy_match = (cdata['date'].dt.dayofyear >= doy_start) & \
                        (cdata['date'].dt.dayofyear <= doy_end)
        else:
            doy_match = (cdata['date'].dt.dayofyear >= doy_start) | \
                        (cdata['date'].dt.dayofyear <= doy_end)
        ctrl_mask = doy_match & ~w
        ctrl_danger = cdata.loc[ctrl_mask, 'dangerLevelMax'].mean()
        if not np.isnan(ssw_danger) and not np.isnan(ctrl_danger) and ctrl_danger > 0:
            event_results.append({
                'centre': centre,
                'rr': ssw_danger / ctrl_danger,
            })
    n_d = sum(1 for r in event_results if r['rr'] < 1)
    print(f"  {ssw_date.date()} ({ssw_type}): {n_d}/{len(event_results)} decrease")

# Multivariate: RR ~ lat + lon
from numpy.linalg import lstsq
X = np.column_stack([np.ones(len(rdf)), rdf['lat'].values, rdf['lon'].values])
y = rdf['mean_rr'].values
beta, _, _, _ = lstsq(X, y, rcond=None)
yhat = X @ beta
r2 = np.corrcoef(y, yhat)[0,1]**2
print(f"\n--- Multivariate: RR ~ lat + lon ---")
print(f"  β_lat = {beta[1]:.4f}, β_lon = {beta[2]:.4f}, R² = {r2:.3f}")

# Save
output = {
    'n_centres': int(len(rdf)),
    'n_decrease': int(n_dec),
    'r_latitude': float(r_lat),
    'p_latitude': float(p_lat),
    'r_longitude': float(r_lon),
    'p_longitude': float(p_lon),
    'r_med_distance': float(r_med),
    'p_med_distance': float(p_med),
    'r2_multivariate': float(r2),
    'centres': rdf.to_dict(orient='records'),
}
with open('data/results/r37_eaws_gradient.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print("\nSaved to data/results/r37_eaws_gradient.json")
