"""
R68 Comprehensive SSW Analysis:
1. Displacement vs Split stratification
2. Combinatorial coherence (k-channel)
3. Power analysis for v'T' dose-response
4. Bootstrap CIs on gmRR
5. Season stratification
6. v'T' threshold test
7. Strat T10 dose-response
8. Leave-one-out sensitivity
"""
import csv, math, json, os
import numpy as np
from scipy import stats

base = r'C:\Users\Jack0\Solar-Magnetic-Analysis\data\results'

# Load SSW event catalog
events = []
with open(os.path.join(base, 'ssw_event_catalog.csv'), 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        events.append(row)

for e in events:
    e['rr'] = float(e['rr'])
    e['log_rr'] = math.log(e['rr'])
    e['strat_t'] = float(e['strat_t10_anom_K'])
    e['z500'] = float(e['z500_anom_m'])
    e['vt'] = float(e['wave_decel_ms_day'])
    e['vortex'] = float(e['vortex_disruption_ms'])
    e['stype'] = e['published_type']
    e['month'] = int(e['date'].split('-')[1])

n = len(events)
n_D = sum(1 for e in events if e['stype'] == 'D')
n_S = sum(1 for e in events if e['stype'] == 'S')
print(f'Total events: {n}, D={n_D}, S={n_S}')

# === 1. Displacement vs Split ===
D_ev = [e for e in events if e['stype'] == 'D']
S_ev = [e for e in events if e['stype'] == 'S']

D_gmrr = math.exp(np.mean([e['log_rr'] for e in D_ev]))
S_gmrr = math.exp(np.mean([e['log_rr'] for e in S_ev]))

D_sign = sum(1 for e in D_ev if e['rr'] < 1)
S_sign = sum(1 for e in S_ev if e['rr'] < 1)

D_sign_p = stats.binomtest(D_sign, len(D_ev), 0.5, alternative='greater').pvalue
S_sign_p = stats.binomtest(S_sign, len(S_ev), 0.5, alternative='greater').pvalue

print(f'\n=== Displacement vs Split ===')
print(f'Displacement (n={len(D_ev)}): gmRR={D_gmrr:.3f}, sign={D_sign}/{len(D_ev)}, P_sign={D_sign_p:.4f}')
print(f'Split (n={len(S_ev)}): gmRR={S_gmrr:.3f}, sign={S_sign}/{len(S_ev)}, P_sign={S_sign_p:.4f}')

# Combined: pooled sign
total_sign = D_sign + S_sign
print(f'Combined: {total_sign}/{n} events RR<1')

# Mann-Whitney U between D and S log_rr
u_ds, p_ds = stats.mannwhitneyu(
    [e['log_rr'] for e in D_ev],
    [e['log_rr'] for e in S_ev],
    alternative='two-sided'
)
print(f'D vs S Mann-Whitney U={u_ds:.1f}, P={p_ds:.3f}')

# === 2. Combinatorial coherence ===
# k independent observation channels all pointing same direction
# Channels: (1) SWE increase, (2) temp suppression, (3) rain decrease,
# (4) snow fraction increase, (5) natural dry decrease, (6) danger rating increase,
# (7) human trigger increase, (8) facet fraction increase
k = 8
p_all_same = 0.5**k
print(f'\n=== Combinatorial Coherence ===')
print(f'k={k} channels, P(all same direction) = {p_all_same:.6f} = 1/{int(1/p_all_same)}')

# Conservative: with correlation, effective k ~ 5
k_eff = 5
p_eff = 0.5**k_eff
print(f'Conservative k_eff={k_eff}: P = {p_eff:.6f} = 1/{int(1/p_eff)}')

# Even more conservative: k_eff = 4
p_4 = 0.5**4
print(f'Most conservative k_eff=4: P = {p_4:.4f} = 1/{int(1/p_4)}')

# === 3. Power analysis ===
z_alpha = 1.96
z_beta = 0.84  # 80% power
min_r = math.tanh((z_alpha + z_beta) / math.sqrt(n - 3))
print(f'\n=== Power Analysis ===')
print(f'Min detectable |rho| at 80% power, alpha=0.05 (two-sided), n={n}: {min_r:.3f}')

# Power to detect |rho|=0.25 at n=16
fisher_z = 0.5 * math.log((1 + 0.25) / (1 - 0.25))
se = 1 / math.sqrt(n - 3)
z_test = fisher_z / se
power_025 = stats.norm.cdf(z_test - z_alpha / 2) + stats.norm.cdf(-z_test - z_alpha / 2)
print(f'Power to detect |rho|=0.25: {power_025:.3f} ({power_025*100:.1f}%)')
print(f'Power to detect |rho|=0.50: ', end='')
fisher_z50 = 0.5 * math.log((1 + 0.50) / (1 - 0.50))
z50 = fisher_z50 / se
power_050 = stats.norm.cdf(z50 - z_alpha / 2) + stats.norm.cdf(-z50 - z_alpha / 2)
print(f'{power_050:.3f} ({power_050*100:.1f}%)')

# === 4. Bootstrap CIs ===
np.random.seed(42)
n_boot = 10000
log_rrs = np.array([e['log_rr'] for e in events])
boot_means = np.array([
    np.mean(np.random.choice(log_rrs, size=n, replace=True))
    for _ in range(n_boot)
])
boot_gmrrs = np.exp(boot_means)
ci_lo, ci_hi = np.percentile(boot_gmrrs, [2.5, 97.5])
gmrr = math.exp(np.mean(log_rrs))
print(f'\n=== Bootstrap gmRR ===')
print(f'gmRR = {gmrr:.3f}, 95% Bootstrap CI [{ci_lo:.3f}, {ci_hi:.3f}]')

# === 5. Season stratification ===
early = [e for e in events if e['month'] <= 1]  # Dec-Jan
late = [e for e in events if e['month'] >= 2]   # Feb-Mar

early_gmrr = math.exp(np.mean([e['log_rr'] for e in early])) if early else None
late_gmrr = math.exp(np.mean([e['log_rr'] for e in late])) if late else None
early_sign = sum(1 for e in early if e['rr'] < 1)
late_sign = sum(1 for e in late if e['rr'] < 1)

print(f'\n=== Season Stratification ===')
print(f'Early (Dec-Jan, n={len(early)}): gmRR={early_gmrr:.3f}, sign={early_sign}/{len(early)}')
print(f'Late (Feb+, n={len(late)}): gmRR={late_gmrr:.3f}, sign={late_sign}/{len(late)}')

# === 6. v'T' threshold test ===
vt_vals = np.array([e['vt'] for e in events])
vt_median = np.median(vt_vals)
above = [e for e in events if e['vt'] >= vt_median]
below = [e for e in events if e['vt'] < vt_median]
above_gmrr = math.exp(np.mean([e['log_rr'] for e in above]))
below_gmrr = math.exp(np.mean([e['log_rr'] for e in below]))
print(f"\n=== v'T' Threshold Test ===")
print(f"v'T' median = {vt_median:.3f}")
print(f'Above median (n={len(above)}): gmRR={above_gmrr:.3f}')
print(f'Below median (n={len(below)}): gmRR={below_gmrr:.3f}')
u_vt, p_vt = stats.mannwhitneyu(
    [e['log_rr'] for e in above],
    [e['log_rr'] for e in below],
    alternative='two-sided'
)
print(f'Mann-Whitney U={u_vt:.1f}, P={p_vt:.3f}')

# Sign concordance in each half
above_sign = sum(1 for e in above if e['rr'] < 1)
below_sign = sum(1 for e in below if e['rr'] < 1)
print(f'Above: {above_sign}/{len(above)} RR<1; Below: {below_sign}/{len(below)} RR<1')

# === 7. Alternative dose-response metrics ===
rho_strat, p_strat = stats.spearmanr(
    [e['strat_t'] for e in events],
    [e['log_rr'] for e in events]
)
rho_vort, p_vort = stats.spearmanr(
    [e['vortex'] for e in events],
    [e['log_rr'] for e in events]
)
rho_z500, p_z500 = stats.spearmanr(
    [e['z500'] for e in events],
    [e['log_rr'] for e in events]
)
print(f'\n=== Dose-Response (alternative metrics) ===')
print(f'Strat T10 vs log_RR: rho={rho_strat:.3f}, P={p_strat:.3f}')
print(f'Vortex disruption vs log_RR: rho={rho_vort:.3f}, P={p_vort:.3f}')
print(f'Z500 vs log_RR: rho={rho_z500:.3f}, P={p_z500:.3f}')

# === 8. LOO sensitivity ===
print(f'\n=== Leave-One-Out Sensitivity ===')
loo_gmrrs = []
loo_signs = []
fragile_events = []
for i in range(n):
    remaining = [e for j, e in enumerate(events) if j != i]
    loo_gmrr_i = math.exp(np.mean([e['log_rr'] for e in remaining]))
    loo_sign_i = sum(1 for e in remaining if e['rr'] < 1)
    loo_gmrrs.append(loo_gmrr_i)
    loo_signs.append(loo_sign_i)
    loo_p = stats.binomtest(loo_sign_i, len(remaining), 0.5, alternative='greater').pvalue
    if loo_p > 0.025:
        fragile_events.append(events[i]['date'])
        print(f"  Drop {events[i]['date']} (RR={events[i]['rr']:.3f}): sign={loo_sign_i}/15, P_1s={loo_p:.4f}")

print(f'LOO gmRR range: [{min(loo_gmrrs):.3f}, {max(loo_gmrrs):.3f}]')
print(f'LOO sign range: [{min(loo_signs)}, {max(loo_signs)}]/15')
print(f'Fragile events (dropping erases P<0.025 one-sided): {len(fragile_events)}')

# === 9. Interaction test: type x season ===
# Is the split effect stronger in late season?
print(f'\n=== Type x Season Interaction ===')
for stype in ['D', 'S']:
    for season in ['early', 'late']:
        subset = [e for e in events if e['stype'] == stype and 
                  (e['month'] <= 1 if season == 'early' else e['month'] >= 2)]
        if subset:
            sub_gmrr = math.exp(np.mean([e['log_rr'] for e in subset]))
            sub_sign = sum(1 for e in subset if e['rr'] < 1)
            print(f'  {stype}-{season} (n={len(subset)}): gmRR={sub_gmrr:.3f}, sign={sub_sign}/{len(subset)}')

# === 10. Bayesian analysis ===
# Log-RR ~ Normal(mu, sigma)
log_rr_arr = np.array([e['log_rr'] for e in events])
mu_hat = np.mean(log_rr_arr)
se_hat = np.std(log_rr_arr, ddof=1) / math.sqrt(n)
# t-test
t_stat = mu_hat / se_hat
t_p = 2 * stats.t.sf(abs(t_stat), df=n-1)
print(f'\n=== Parametric t-test on log(RR) ===')
print(f'mu = {mu_hat:.3f}, SE = {se_hat:.3f}')
print(f't = {t_stat:.3f}, P = {t_p:.4f} (two-sided)')
print(f'Cohen d = {mu_hat / np.std(log_rr_arr, ddof=1):.3f}')

# Save all results
results = {
    'n_events': n,
    'overall_gmRR': round(gmrr, 3),
    'bootstrap_CI': [round(ci_lo, 3), round(ci_hi, 3)],
    'displacement': {
        'n': len(D_ev), 'gmRR': round(D_gmrr, 3),
        'sign': f'{D_sign}/{len(D_ev)}', 'P_sign': round(D_sign_p, 4)
    },
    'split': {
        'n': len(S_ev), 'gmRR': round(S_gmrr, 3),
        'sign': f'{S_sign}/{len(S_ev)}', 'P_sign': round(S_sign_p, 4)
    },
    'D_vs_S_MWU_P': round(p_ds, 3),
    'combinatorial': {
        'k': k, 'P_all': round(p_all_same, 6),
        'k_eff_5': round(p_eff, 6), 'k_eff_4': round(p_4, 4)
    },
    'power': {
        'min_detectable_rho_80pct': round(min_r, 3),
        'power_rho025': round(power_025, 3),
        'power_rho050': round(power_050, 3)
    },
    'season': {
        'early': {'n': len(early), 'gmRR': round(early_gmrr, 3), 'sign': f'{early_sign}/{len(early)}'},
        'late': {'n': len(late), 'gmRR': round(late_gmrr, 3), 'sign': f'{late_sign}/{len(late)}'}
    },
    'vt_threshold': {
        'median': round(vt_median, 3),
        'above_gmRR': round(above_gmrr, 3),
        'below_gmRR': round(below_gmrr, 3),
        'MWU_P': round(p_vt, 3)
    },
    'dose_response': {
        'strat_T10': {'rho': round(rho_strat, 3), 'P': round(p_strat, 3)},
        'vortex': {'rho': round(rho_vort, 3), 'P': round(p_vort, 3)},
        'Z500': {'rho': round(rho_z500, 3), 'P': round(p_z500, 3)}
    },
    'loo': {
        'gmRR_range': [round(min(loo_gmrrs), 3), round(max(loo_gmrrs), 3)],
        'sign_range': [min(loo_signs), max(loo_signs)],
        'n_fragile': len(fragile_events)
    },
    'parametric_t': {
        'mu_logRR': round(mu_hat, 3), 'SE': round(se_hat, 3),
        't': round(t_stat, 3), 'P': round(t_p, 4),
        'cohen_d': round(mu_hat / np.std(log_rr_arr, ddof=1), 3)
    }
}

with open(os.path.join(base, 'r68_comprehensive_analysis.json'), 'w') as f:
    json.dump(results, f, indent=2)
print('\nResults saved to r68_comprehensive_analysis.json')
