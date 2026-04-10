"""Extract peak propagation values from R14 results."""
import json, numpy as np

with open('data/results/r14_mechanism_breakthrough.json') as f:
    d = json.load(f)

p4 = d['phase4_propagation']
lags = p4['lags']  # -20 to +30
levels_data = p4['levels']

print("=== DISPLACEMENT SSW PROPAGATION (n=%d events) ===" % p4['n_events'])
print()

for level_name in ['ncep_t_10hpa_anom', 'ncep_t_20hpa_anom', 'ncep_t_30hpa_anom', 
                    'ncep_t_50hpa_anom', 'ncep_t_70hpa_anom', 'ncep_t_100hpa_anom', 't2m_anom']:
    if level_name in levels_data:
        means = levels_data[level_name]['means']
        pvals = levels_data[level_name]['p_values']
        
        # Find peak (maximum mean) for positive lags only (lag >= 0)
        pos_mask = [i for i, lag in enumerate(lags) if lag >= 0]
        pos_means = [means[i] for i in pos_mask]
        pos_lags = [lags[i] for i in pos_mask]
        pos_pvals = [pvals[i] for i in pos_mask]
        
        peak_idx = np.argmax(pos_means)
        peak_lag = pos_lags[peak_idx]
        peak_mean = pos_means[peak_idx]
        peak_p = pos_pvals[peak_idx]
        
        # Find minimum P value for positive lags
        min_p_idx = np.argmin(pos_pvals)
        min_p_lag = pos_lags[min_p_idx]
        min_p = pos_pvals[min_p_idx]
        
        label = level_name.replace('ncep_t_', '').replace('_anom', '').replace('hpa', ' hPa')
        print("%s:" % label)
        print("  Peak anomaly: +%.1f K at lag +%d d (P=%.6f)" % (peak_mean, peak_lag, peak_p))
        print("  Min P value:  P=%.2e at lag +%d d (anom=+%.1f K)" % (min_p, min_p_lag, pos_means[min_p_idx]))
        print()

# Also check the NH circulation results (Phase 7 - if it's in the console output)
print("=== PHASE 1: ERA5 COMPOSITES ===")
p1 = d['phase1_composites']
for ssw_type, data in p1.items():
    if isinstance(data, dict):
        print("%s: n=%s, T2m post=%.3f K, P=%.4f" % (
            ssw_type, data.get('n_events','?'), 
            data.get('t2m_post_mean_K', 0), data.get('t2m_post_p', 1)))
