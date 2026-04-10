"""
23_resolve_limitations.py — Resolve all paper limitations with new analyses.
Run parts sequentially: python 23_resolve_limitations.py --part N
"""
import sys, json, warnings, pathlib
import numpy as np
import pandas as pd
from itertools import product
warnings.filterwarnings('ignore')

PROCESSED = pathlib.Path('data/processed')
RESULTS   = pathlib.Path('data/results')
FIGURES   = pathlib.Path('data/figures')
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

panel = pd.read_parquet(PROCESSED / 'analysis_panel_v2.parquet')
ssw   = pd.read_parquet(PROCESSED / 'atmospheric' / 'ssw_catalog.parquet')
ssw.index = pd.to_datetime(ssw.index).tz_localize(None)
winter = panel[panel['is_winter'] == 1].copy()

# Get SSW dates in study period
ssw_dates = ssw.index[(ssw.index >= winter.index.min()) & (ssw.index <= winter.index.max())]
print(f"Panel: {len(panel)} rows, Winter: {len(winter)} rows, SSW events: {len(ssw_dates)}")

def get_part():
    for a in sys.argv[1:]:
        if a.startswith('--part'):
            return int(sys.argv[sys.argv.index(a)+1])
    return int(sys.argv[-1]) if len(sys.argv) > 1 else 1

part = get_part()

# ============================================================
# PART 1: Formal pre-vs-post SSW test (exact event-level)
# ============================================================
if part == 1:
    print("\n" + "="*60)
    print("PART 1: FORMAL PRE-VS-POST SSW TEST")
    print("="*60)
    
    results = {}
    
    for outcome_col, label in [('dry_natural_size_1234', 'dry_natural'),
                                ('norway_aval_count', 'norway_total')]:
        print(f"\n--- {label} ---")
        
        pre_deltas = []
        post_deltas = []
        
        for ssw_date in ssw_dates:
            # Post window: 0–15d after SSW
            post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
            post_vals = winter.loc[post_mask, outcome_col].dropna()
            
            # Pre window: 15–0d before SSW
            pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
            pre_vals = winter.loc[pre_mask, outcome_col].dropna()
            
            # Matched control: same calendar days, nearest non-SSW winters
            ssw_year = ssw_date.year if ssw_date.month < 7 else ssw_date.year + 1
            control_means_post = []
            control_means_pre = []
            
            for offset in [-2, -1, 1, 2]:
                ctrl_year = ssw_year + offset
                # Check this isn't also an SSW winter
                is_ssw_winter = False
                for sd in ssw_dates:
                    sd_year = sd.year if sd.month < 7 else sd.year + 1
                    if sd_year == ctrl_year:
                        is_ssw_winter = True
                        break
                if is_ssw_winter:
                    continue
                
                # Same calendar days shifted by offset years
                post_ctrl_start = ssw_date.replace(year=ssw_date.year + offset)
                post_ctrl_end = post_ctrl_start + pd.Timedelta(days=15)
                pre_ctrl_start = post_ctrl_start - pd.Timedelta(days=15)
                
                try:
                    post_ctrl = winter.loc[(winter.index >= post_ctrl_start) & 
                                          (winter.index < post_ctrl_end), outcome_col].dropna()
                    pre_ctrl = winter.loc[(winter.index >= pre_ctrl_start) & 
                                         (winter.index < post_ctrl_start), outcome_col].dropna()
                    if len(post_ctrl) > 5:
                        control_means_post.append(post_ctrl.mean())
                    if len(pre_ctrl) > 5:
                        control_means_pre.append(pre_ctrl.mean())
                except:
                    continue
            
            if len(post_vals) > 5 and len(pre_vals) > 5 and control_means_post and control_means_pre:
                post_delta = post_vals.mean() - np.mean(control_means_post)
                pre_delta = pre_vals.mean() - np.mean(control_means_pre)
                post_deltas.append(post_delta)
                pre_deltas.append(pre_delta)
        
        post_deltas = np.array(post_deltas)
        pre_deltas = np.array(pre_deltas)
        n = len(post_deltas)
        diff = post_deltas - pre_deltas  # post - pre (negative means post is MORE suppressed)
        
        print(f"  N events with both windows: {n}")
        print(f"  Mean post Δ: {post_deltas.mean():.3f}")
        print(f"  Mean pre Δ:  {pre_deltas.mean():.3f}")
        print(f"  Mean diff (post-pre): {diff.mean():.3f}")
        print(f"  Median diff: {np.median(diff):.3f}")
        
        # Exact sign-flip permutation on differences (2^n permutations)
        obs_mean_diff = diff.mean()
        n_perm = 2**n
        count_extreme = 0
        for i in range(n_perm):
            signs = np.array([1 if (i >> j) & 1 else -1 for j in range(n)])
            perm_mean = (signs * np.abs(diff)).mean()
            if perm_mean <= obs_mean_diff:  # one-sided: post more negative than pre
                count_extreme += 1
        perm_p_onesided = count_extreme / n_perm
        perm_p_twosided = min(1.0, 2 * perm_p_onesided)
        
        print(f"  Exact permutation P (one-sided, post < pre): {perm_p_onesided:.4f}")
        print(f"  Exact permutation P (two-sided): {perm_p_twosided:.4f}")
        
        # Bootstrap CI on the difference
        np.random.seed(42)
        boot_diffs = []
        for _ in range(10000):
            idx = np.random.choice(n, n, replace=True)
            boot_diffs.append(diff[idx].mean())
        boot_diffs = np.array(boot_diffs)
        ci_lo, ci_hi = np.percentile(boot_diffs, [2.5, 97.5])
        print(f"  Bootstrap 95% CI on diff: [{ci_lo:.3f}, {ci_hi:.3f}]")
        
        # Leave-one-out influence
        loo_means = []
        for i in range(n):
            loo = np.delete(diff, i)
            loo_means.append(loo.mean())
        loo_means = np.array(loo_means)
        most_influential = np.argmax(np.abs(loo_means - diff.mean()))
        print(f"  LOO range: [{loo_means.min():.3f}, {loo_means.max():.3f}]")
        print(f"  Most influential event #{most_influential}: removing it gives mean={loo_means[most_influential]:.3f}")
        
        # Paired t-test
        from scipy import stats
        t_stat, t_p = stats.ttest_rel(post_deltas, pre_deltas)
        print(f"  Paired t-test: t={t_stat:.3f}, P={t_p:.4f}")
        
        # Wilcoxon signed-rank on differences
        w_stat, w_p = stats.wilcoxon(diff, alternative='two-sided')
        print(f"  Wilcoxon on diff: W={w_stat:.1f}, P={w_p:.4f}")
        
        results[label] = {
            'n_events': int(n),
            'mean_post_delta': float(post_deltas.mean()),
            'mean_pre_delta': float(pre_deltas.mean()),
            'mean_diff': float(diff.mean()),
            'median_diff': float(np.median(diff)),
            'exact_perm_p_onesided': float(perm_p_onesided),
            'exact_perm_p_twosided': float(perm_p_twosided),
            'bootstrap_ci': [float(ci_lo), float(ci_hi)],
            'loo_range': [float(loo_means.min()), float(loo_means.max())],
            'paired_t_p': float(t_p),
            'wilcoxon_p': float(w_p),
            'post_deltas': post_deltas.tolist(),
            'pre_deltas': pre_deltas.tolist()
        }
    
    with open(RESULTS / 'lim_part1_prepost.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\n✓ Part 1 saved to", RESULTS / 'lim_part1_prepost.json')

# ============================================================
# PART 2: Distributed-lag event-study model
# ============================================================
elif part == 2:
    print("\n" + "="*60)
    print("PART 2: DISTRIBUTED-LAG EVENT-STUDY MODEL")
    print("="*60)
    
    results = {}
    
    for outcome_col, label in [('dry_natural_size_1234', 'dry_natural'),
                                ('norway_aval_count', 'norway_total')]:
        print(f"\n--- {label} ---")
        
        # Build seasonal expectation (day-of-season mean across all winters)
        winter_copy = winter.copy()
        winter_copy['dos'] = winter_copy['day_of_season']
        seasonal = winter_copy.groupby('dos')[outcome_col].mean()
        winter_copy['seasonal_expected'] = winter_copy['dos'].map(seasonal)
        winter_copy['anomaly'] = winter_copy[outcome_col] - winter_copy['seasonal_expected']
        
        # Event-study: stack anomalies around each SSW event
        lag_range = np.arange(-30, 61)  # -30 to +60 days
        event_anomalies = {lag: [] for lag in lag_range}
        
        for ssw_date in ssw_dates:
            for lag in lag_range:
                target_date = ssw_date + pd.Timedelta(days=int(lag))
                if target_date in winter_copy.index:
                    val = winter_copy.loc[target_date, 'anomaly']
                    if not np.isnan(val):
                        event_anomalies[lag].append(float(val))
        
        # Compute mean anomaly and CI at each lag
        lag_means = []
        lag_cis_lo = []
        lag_cis_hi = []
        lag_ns = []
        lag_ps = []
        
        for lag in lag_range:
            vals = np.array(event_anomalies[lag])
            n = len(vals)
            if n >= 5:
                m = vals.mean()
                se = vals.std(ddof=1) / np.sqrt(n)
                from scipy import stats
                ci_lo = m - 1.96 * se
                ci_hi = m + 1.96 * se
                _, p = stats.ttest_1samp(vals, 0)
            else:
                m, ci_lo, ci_hi, p = np.nan, np.nan, np.nan, np.nan
            lag_means.append(float(m))
            lag_cis_lo.append(float(ci_lo))
            lag_cis_hi.append(float(ci_hi))
            lag_ns.append(int(n))
            lag_ps.append(float(p) if not np.isnan(p) else None)
        
        # Compute 7-day rolling average for smoothing
        lag_means_arr = np.array(lag_means)
        smoothed = pd.Series(lag_means_arr).rolling(7, center=True, min_periods=3).mean().values
        
        # Summary: mean anomaly in pre/post windows
        pre_idx = [i for i, l in enumerate(lag_range) if -15 <= l < 0]
        post_idx = [i for i, l in enumerate(lag_range) if 0 <= l < 15]
        late_idx = [i for i, l in enumerate(lag_range) if 15 <= l < 30]
        
        pre_mean = np.nanmean([lag_means[i] for i in pre_idx])
        post_mean = np.nanmean([lag_means[i] for i in post_idx])
        late_mean = np.nanmean([lag_means[i] for i in late_idx])
        
        print(f"  Pre-SSW (-15 to 0d) mean anomaly: {pre_mean:.3f}")
        print(f"  Post-SSW (0 to 15d) mean anomaly: {post_mean:.3f}")
        print(f"  Late post (15 to 30d) mean anomaly: {late_mean:.3f}")
        
        # Count significant lags
        sig_pre = sum(1 for i in pre_idx if lag_ps[i] is not None and lag_ps[i] < 0.05)
        sig_post = sum(1 for i in post_idx if lag_ps[i] is not None and lag_ps[i] < 0.05)
        print(f"  Significant lags (P<0.05): pre={sig_pre}/{len(pre_idx)}, post={sig_post}/{len(post_idx)}")
        
        results[label] = {
            'lags': lag_range.tolist(),
            'mean_anomaly': lag_means,
            'ci_lo': lag_cis_lo,
            'ci_hi': lag_cis_hi,
            'smoothed_7d': smoothed.tolist(),
            'n_events_per_lag': lag_ns,
            'p_values': lag_ps,
            'summary': {
                'pre_15_0_mean': float(pre_mean),
                'post_0_15_mean': float(post_mean),
                'post_15_30_mean': float(late_mean),
                'sig_lags_pre': sig_pre,
                'sig_lags_post': sig_post
            }
        }
    
    # Generate figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    for idx, (label, title) in enumerate([('dry_natural', 'Swiss Dry Slab Avalanches'),
                                           ('norway_total', 'Norwegian Total Avalanches')]):
        ax = axes[idx]
        r = results[label]
        lags = r['lags']
        means = np.array(r['mean_anomaly'])
        ci_lo = np.array(r['ci_lo'])
        ci_hi = np.array(r['ci_hi'])
        smoothed = np.array(r['smoothed_7d'])
        
        ax.fill_between(lags, ci_lo, ci_hi, alpha=0.2, color='steelblue')
        ax.plot(lags, means, 'o', ms=2, alpha=0.4, color='steelblue')
        ax.plot(lags, smoothed, '-', lw=2, color='navy', label='7-day rolling mean')
        ax.axhline(0, color='gray', ls='--', lw=0.8)
        ax.axvline(0, color='red', ls='-', lw=1.5, alpha=0.7, label='SSW onset')
        ax.axvspan(-15, 0, alpha=0.05, color='orange', label='Pre-SSW window')
        ax.axvspan(0, 15, alpha=0.05, color='green', label='Post-SSW window')
        ax.set_ylabel('Anomaly (events/day)')
        ax.set_title(title)
        ax.legend(fontsize=8)
    
    axes[1].set_xlabel('Days relative to SSW onset')
    plt.suptitle('Event-Study: Daily Avalanche Anomaly Around SSW Events\n(Seasonal expectation subtracted)', 
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES / 'fig_event_study_lag.png', dpi=200, bbox_inches='tight')
    plt.savefig(FIGURES / 'fig_event_study_lag.pdf', bbox_inches='tight')
    plt.close()
    print("\n✓ Event-study figure saved")
    
    with open(RESULTS / 'lim_part2_eventlag.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("✓ Part 2 saved to", RESULTS / 'lim_part2_eventlag.json')

# ============================================================
# PART 3: Pre-specified meteorological chain
# ============================================================
elif part == 3:
    print("\n" + "="*60)
    print("PART 3: PRE-SPECIFIED METEOROLOGICAL CHAIN")
    print("="*60)
    print("Chain: 10hPa Temp → AO (u850) → Z500/SLP → Avalanche magnitude")
    
    from scipy import stats
    results = {}
    
    # Step 1: NCEP 10hPa temperature response around SSW
    print("\n--- Step 1: Stratospheric temperature (10hPa) ---")
    pre_temps = []
    post_temps = []
    for ssw_date in ssw_dates:
        pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
        post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
        pre_t = winter.loc[pre_mask, 'ncep_t_10hpa'].dropna()
        post_t = winter.loc[post_mask, 'ncep_t_10hpa'].dropna()
        if len(pre_t) > 5 and len(post_t) > 5:
            pre_temps.append(pre_t.mean())
            post_temps.append(post_t.mean())
    
    pre_temps = np.array(pre_temps)
    post_temps = np.array(post_temps)
    temp_diff = post_temps - pre_temps
    t_stat, t_p = stats.ttest_rel(post_temps, pre_temps)
    print(f"  N events: {len(temp_diff)}")
    print(f"  Mean pre: {pre_temps.mean():.1f} K, Mean post: {post_temps.mean():.1f} K")
    print(f"  Mean warming: {temp_diff.mean():.1f} K")
    print(f"  Paired t: t={t_stat:.2f}, P={t_p:.6f}")
    
    results['step1_strat_temp'] = {
        'variable': 'ncep_t_10hpa',
        'n': len(temp_diff),
        'mean_pre': float(pre_temps.mean()),
        'mean_post': float(post_temps.mean()),
        'mean_change': float(temp_diff.mean()),
        'paired_t_stat': float(t_stat),
        'paired_t_p': float(t_p),
        'n_positive': int(np.sum(temp_diff > 0)),
        'event_changes': temp_diff.tolist()
    }
    
    # Step 2: AO proxy — use ncep_u_10hpa (zonal wind, captures vortex weakening)
    # and ncep_u850_nh (low-level jet, captures tropospheric response)
    print("\n--- Step 2: Zonal wind (vortex) and tropospheric response ---")
    for var, var_label in [('ncep_u_10hpa', 'Strat zonal wind 10hPa'),
                           ('ncep_u850_nh', 'Trop U850'),
                           ('ncep_z500_nh', 'Z500'),
                           ('ncep_slp_nh', 'SLP')]:
        pre_vals = []
        post_vals = []
        for ssw_date in ssw_dates:
            pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
            post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
            pre_v = winter.loc[pre_mask, var].dropna()
            post_v = winter.loc[post_mask, var].dropna()
            if len(pre_v) > 5 and len(post_v) > 5:
                pre_vals.append(pre_v.mean())
                post_vals.append(post_v.mean())
        
        pre_vals = np.array(pre_vals)
        post_vals = np.array(post_vals)
        changes = post_vals - pre_vals
        t_stat, t_p = stats.ttest_rel(post_vals, pre_vals)
        w_stat, w_p = stats.wilcoxon(changes)
        
        print(f"  {var_label}: mean change = {changes.mean():.3f}, t-P={t_p:.4f}, W-P={w_p:.4f}")
        
        results[f'step2_{var}'] = {
            'variable': var,
            'label': var_label,
            'n': len(changes),
            'mean_pre': float(pre_vals.mean()),
            'mean_post': float(post_vals.mean()),
            'mean_change': float(changes.mean()),
            'paired_t_p': float(t_p),
            'wilcoxon_p': float(w_p),
            'n_positive': int(np.sum(changes > 0)),
            'event_changes': changes.tolist()
        }
    
    # Step 3: Event-by-event correlation — strat warming magnitude vs avalanche decrease
    print("\n--- Step 3: Event-by-event correlation ---")
    strat_warmings = []
    aval_changes = []
    
    for ssw_date in ssw_dates:
        pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
        post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
        
        pre_t = winter.loc[pre_mask, 'ncep_t_10hpa'].dropna()
        post_t = winter.loc[post_mask, 'ncep_t_10hpa'].dropna()
        
        pre_a = winter.loc[pre_mask, 'dry_natural_size_1234'].dropna()
        post_a = winter.loc[post_mask, 'dry_natural_size_1234'].dropna()
        
        if len(pre_t) > 5 and len(post_t) > 5 and len(pre_a) > 5 and len(post_a) > 5:
            strat_warmings.append(post_t.mean() - pre_t.mean())
            aval_changes.append(post_a.mean() - pre_a.mean())
    
    strat_warmings = np.array(strat_warmings)
    aval_changes = np.array(aval_changes)
    
    r_pearson, p_pearson = stats.pearsonr(strat_warmings, aval_changes)
    r_spearman, p_spearman = stats.spearmanr(strat_warmings, aval_changes)
    
    print(f"  N events: {len(strat_warmings)}")
    print(f"  Pearson r={r_pearson:.3f}, P={p_pearson:.4f}")
    print(f"  Spearman ρ={r_spearman:.3f}, P={p_spearman:.4f}")
    
    results['step3_event_correlation'] = {
        'n': len(strat_warmings),
        'pearson_r': float(r_pearson),
        'pearson_p': float(p_pearson),
        'spearman_rho': float(r_spearman),
        'spearman_p': float(p_spearman),
        'strat_warmings': strat_warmings.tolist(),
        'aval_changes': aval_changes.tolist()
    }
    
    # Step 4: AO index (nao_daily as proxy, also check if we have AO)
    # Test direct avalanche prediction: does AO/NAO predict dry aval rate?
    print("\n--- Step 4: NAO/AO → avalanche direct relationship ---")
    valid = winter[['nao_daily', 'dry_natural_size_1234']].dropna()
    r_nao, p_nao = stats.pearsonr(valid['nao_daily'], valid['dry_natural_size_1234'])
    print(f"  NAO-daily vs dry aval: r={r_nao:.3f}, P={p_nao:.6f}, N={len(valid)}")
    
    # Bin NAO into terciles
    nao_terciles = pd.qcut(valid['nao_daily'], 3, labels=['NAO-', 'NAO0', 'NAO+'])
    for t in ['NAO-', 'NAO0', 'NAO+']:
        subset = valid.loc[nao_terciles == t, 'dry_natural_size_1234']
        print(f"    {t}: mean dry aval = {subset.mean():.3f} (n={len(subset)})")
    
    results['step4_nao_aval'] = {
        'n': len(valid),
        'pearson_r': float(r_nao),
        'pearson_p': float(p_nao),
        'tercile_means': {
            str(t): float(valid.loc[nao_terciles == t, 'dry_natural_size_1234'].mean())
            for t in ['NAO-', 'NAO0', 'NAO+']
        }
    }
    
    # Bonferroni correction for 4 main tests
    main_pvals = [
        results['step1_strat_temp']['paired_t_p'],
        results['step2_ncep_u_10hpa']['paired_t_p'],
        results['step2_ncep_z500_nh']['paired_t_p'],
        results['step3_event_correlation']['pearson_p']
    ]
    bonferroni = [min(1.0, p * 4) for p in main_pvals]
    results['bonferroni_correction'] = {
        'raw_pvals': main_pvals,
        'corrected_pvals': bonferroni,
        'labels': ['strat_temp', 'u10hpa', 'z500', 'event_corr']
    }
    print(f"\n  Bonferroni-corrected P-values: {dict(zip(['strat_temp','u10hpa','z500','event_corr'], bonferroni))}")
    
    # Generate correlation scatter figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    ax.scatter(strat_warmings, aval_changes, s=60, c='steelblue', edgecolors='navy', zorder=3)
    # Regression line
    slope, intercept = np.polyfit(strat_warmings, aval_changes, 1)
    x_line = np.linspace(strat_warmings.min(), strat_warmings.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, '--', color='red', lw=1.5)
    ax.axhline(0, color='gray', ls=':', lw=0.8)
    ax.axvline(0, color='gray', ls=':', lw=0.8)
    ax.set_xlabel('Stratospheric warming (ΔT 10hPa, K)')
    ax.set_ylabel('Dry slab avalanche change (Δ events/day)')
    ax.set_title(f'Event-by-Event: Stratospheric Warming vs Avalanche Response\n'
                 f'Pearson r={r_pearson:.2f} (P={p_pearson:.3f}), Spearman ρ={r_spearman:.2f} (P={p_spearman:.3f})')
    plt.tight_layout()
    plt.savefig(FIGURES / 'fig_event_correlation.png', dpi=200, bbox_inches='tight')
    plt.savefig(FIGURES / 'fig_event_correlation.pdf', bbox_inches='tight')
    plt.close()
    print("✓ Correlation figure saved")
    
    with open(RESULTS / 'lim_part3_metchain.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("✓ Part 3 saved")

# ============================================================
# PART 4: Norwegian cold-period sensitivity
# ============================================================
elif part == 4:
    print("\n" + "="*60)
    print("PART 4: NORWEGIAN COLD-PERIOD SENSITIVITY")
    print("="*60)
    print("(Supportive sensitivity analysis — NOT claiming dry/wet recovery)")
    
    from scipy import stats
    results = {}
    
    # Use stratospheric temperature as proxy for cold/warm regime
    # and ncep_t_100hpa (lower stratosphere, closer to surface)
    # Split Norwegian winter days into cold vs warm halves
    norway_valid = winter[['norway_aval_count', 'ncep_t_100hpa', 'ssw_within_15d']].dropna()
    
    median_t = norway_valid['ncep_t_100hpa'].median()
    norway_valid['cold_regime'] = (norway_valid['ncep_t_100hpa'] < median_t).astype(int)
    
    print(f"  Total days: {len(norway_valid)}")
    print(f"  Median 100hPa temp: {median_t:.1f} K")
    print(f"  Cold regime days: {norway_valid['cold_regime'].sum()}")
    print(f"  Warm regime days: {(~norway_valid['cold_regime'].astype(bool)).sum()}")
    
    # SSW effect in cold vs warm regimes
    for regime, label in [(1, 'Cold (below-median 100hPa T)'), (0, 'Warm (above-median 100hPa T)')]:
        subset = norway_valid[norway_valid['cold_regime'] == regime]
        ssw_days = subset[subset['ssw_within_15d'] == 1]['norway_aval_count']
        non_ssw_days = subset[subset['ssw_within_15d'] == 0]['norway_aval_count']
        
        if len(ssw_days) > 10 and len(non_ssw_days) > 10:
            t_stat, t_p = stats.ttest_ind(ssw_days, non_ssw_days)
            mw_stat, mw_p = stats.mannwhitneyu(ssw_days, non_ssw_days, alternative='two-sided')
            
            print(f"\n  {label}:")
            print(f"    SSW days: mean={ssw_days.mean():.2f}, n={len(ssw_days)}")
            print(f"    Non-SSW: mean={non_ssw_days.mean():.2f}, n={len(non_ssw_days)}")
            print(f"    Diff: {ssw_days.mean() - non_ssw_days.mean():.2f}")
            print(f"    t-test P={t_p:.4f}, Mann-Whitney P={mw_p:.4f}")
            
            results[f'regime_{["warm","cold"][regime]}'] = {
                'label': label,
                'ssw_mean': float(ssw_days.mean()),
                'ssw_n': int(len(ssw_days)),
                'non_ssw_mean': float(non_ssw_days.mean()),
                'non_ssw_n': int(len(non_ssw_days)),
                'diff': float(ssw_days.mean() - non_ssw_days.mean()),
                't_p': float(t_p),
                'mw_p': float(mw_p)
            }
    
    # Also stratify by season: early winter (Nov-Jan, more dry) vs late (Feb-Apr, more wet)
    print("\n--- Season stratification ---")
    norway_valid_s = winter[['norway_aval_count', 'ssw_within_15d', 'month']].dropna()
    norway_valid_s['early_winter'] = norway_valid_s['month'].isin([11, 12, 1]).astype(int)
    
    for season, label in [(1, 'Early winter (Nov-Jan, dry-dominant)'), (0, 'Late winter (Feb-Apr, wet-dominant)')]:
        subset = norway_valid_s[norway_valid_s['early_winter'] == season]
        ssw_days = subset[subset['ssw_within_15d'] == 1]['norway_aval_count']
        non_ssw_days = subset[subset['ssw_within_15d'] == 0]['norway_aval_count']
        
        if len(ssw_days) > 5 and len(non_ssw_days) > 10:
            t_stat, t_p = stats.ttest_ind(ssw_days, non_ssw_days)
            print(f"\n  {label}:")
            print(f"    SSW days: mean={ssw_days.mean():.2f}, n={len(ssw_days)}")
            print(f"    Non-SSW: mean={non_ssw_days.mean():.2f}, n={len(non_ssw_days)}")
            print(f"    Diff: {ssw_days.mean() - non_ssw_days.mean():.2f}")
            print(f"    t-test P={t_p:.4f}")
            
            results[f'season_{["late","early"][season]}'] = {
                'label': label,
                'ssw_mean': float(ssw_days.mean()),
                'ssw_n': int(len(ssw_days)),
                'non_ssw_mean': float(non_ssw_days.mean()),
                'non_ssw_n': int(len(non_ssw_days)),
                'diff': float(ssw_days.mean() - non_ssw_days.mean()),
                't_p': float(t_p)
            }
    
    with open(RESULTS / 'lim_part4_norway.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\n✓ Part 4 saved")

# ============================================================
# PART 5: SNOTEL western US extension
# ============================================================
elif part == 5:
    print("\n" + "="*60)
    print("PART 5: SNOTEL WESTERN US EXTENSION")
    print("="*60)
    
    from scipy import stats
    results = {}
    
    for var, label in [('snotel_swe_mean', 'SWE (snow water equivalent)'),
                       ('snotel_prec_mean', 'Precipitation'),
                       ('snotel_temp_mean', 'Temperature')]:
        print(f"\n--- {label} ({var}) ---")
        
        pre_vals = []
        post_vals = []
        
        for ssw_date in ssw_dates:
            pre_mask = (winter.index >= ssw_date - pd.Timedelta(days=15)) & (winter.index < ssw_date)
            post_mask = (winter.index >= ssw_date) & (winter.index < ssw_date + pd.Timedelta(days=15))
            
            pre_v = winter.loc[pre_mask, var].dropna()
            post_v = winter.loc[post_mask, var].dropna()
            
            if len(pre_v) > 5 and len(post_v) > 5:
                pre_vals.append(pre_v.mean())
                post_vals.append(post_v.mean())
        
        if len(pre_vals) >= 5:
            pre_arr = np.array(pre_vals)
            post_arr = np.array(post_vals)
            changes = post_arr - pre_arr
            
            t_stat, t_p = stats.ttest_rel(post_arr, pre_arr)
            n_pos = int(np.sum(changes > 0))
            
            print(f"  N events: {len(changes)}")
            print(f"  Mean pre: {pre_arr.mean():.3f}")
            print(f"  Mean post: {post_arr.mean():.3f}")
            print(f"  Mean change: {changes.mean():.3f}")
            print(f"  {n_pos}/{len(changes)} positive")
            print(f"  Paired t: P={t_p:.4f}")
            
            # Also compare to seasonal expectation
            winter_copy = winter.copy()
            winter_copy['dos'] = winter_copy['day_of_season']
            seasonal = winter_copy.groupby('dos')[var].mean()
            
            results[var] = {
                'label': label,
                'n_events': len(changes),
                'mean_pre': float(pre_arr.mean()),
                'mean_post': float(post_arr.mean()),
                'mean_change': float(changes.mean()),
                'n_positive': n_pos,
                'paired_t_p': float(t_p),
                'event_changes': changes.tolist()
            }
        else:
            print(f"  Insufficient data (only {len(pre_vals)} events)")
            results[var] = {'label': label, 'insufficient_data': True}
    
    # SWE change rate (difference in daily SWE change pre vs post SSW)
    print("\n--- SWE change rate (dSWE/dt) ---")
    winter_copy = winter.copy()
    winter_copy['dswe'] = winter_copy['snotel_swe_mean'].diff()
    
    pre_rates = []
    post_rates = []
    for ssw_date in ssw_dates:
        pre_mask = (winter_copy.index >= ssw_date - pd.Timedelta(days=15)) & (winter_copy.index < ssw_date)
        post_mask = (winter_copy.index >= ssw_date) & (winter_copy.index < ssw_date + pd.Timedelta(days=15))
        
        pre_r = winter_copy.loc[pre_mask, 'dswe'].dropna()
        post_r = winter_copy.loc[post_mask, 'dswe'].dropna()
        
        if len(pre_r) > 5 and len(post_r) > 5:
            pre_rates.append(pre_r.mean())
            post_rates.append(post_r.mean())
    
    if len(pre_rates) >= 5:
        pre_arr = np.array(pre_rates)
        post_arr = np.array(post_rates)
        changes = post_arr - pre_arr
        t_stat, t_p = stats.ttest_rel(post_arr, pre_arr)
        print(f"  N events: {len(changes)}")
        print(f"  Mean pre dSWE/dt: {pre_arr.mean():.4f}")
        print(f"  Mean post dSWE/dt: {post_arr.mean():.4f}")
        print(f"  Change in accumulation rate: {changes.mean():.4f}")
        print(f"  Paired t: P={t_p:.4f}")
        
        results['swe_rate'] = {
            'label': 'SWE accumulation rate (dSWE/dt)',
            'n_events': len(changes),
            'mean_pre_rate': float(pre_arr.mean()),
            'mean_post_rate': float(post_arr.mean()),
            'mean_change': float(changes.mean()),
            'paired_t_p': float(t_p)
        }
    
    with open(RESULTS / 'lim_part5_snotel.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\n✓ Part 5 saved")

# ============================================================
# PART 6: Improved prediction model
# ============================================================
elif part == 6:
    print("\n" + "="*60)
    print("PART 6: IMPROVED PREDICTION MODEL")
    print("="*60)
    print("Leave-one-winter-out blocked CV with SSW + met covariates")
    
    from scipy import stats
    import statsmodels.api as sm
    
    results = {}
    
    # Prepare data
    model_data = winter[['dry_natural_size_1234', 'day_of_season', 'day_of_season_sq',
                         'ssw_within_15d', 'nao_daily', 'ncep_t_10hpa',
                         'post_event_1_3d', 'winter_id']].dropna()
    
    print(f"  Model data: {len(model_data)} rows, {model_data['winter_id'].nunique()} winters")
    
    winters = sorted(model_data['winter_id'].unique())
    
    # Model comparison: baseline (seasonal) vs enhanced (+ SSW + strat temp + geomag)
    baseline_rmse = []
    enhanced_rmse = []
    baseline_mae = []
    enhanced_mae = []
    baseline_ll = []
    enhanced_ll = []
    
    for held_out in winters:
        train = model_data[model_data['winter_id'] != held_out]
        test = model_data[model_data['winter_id'] == held_out]
        
        if len(test) < 10:
            continue
        
        y_train = train['dry_natural_size_1234']
        y_test = test['dry_natural_size_1234']
        
        # Baseline: seasonal only
        X_base_train = sm.add_constant(train[['day_of_season', 'day_of_season_sq']])
        X_base_test = sm.add_constant(test[['day_of_season', 'day_of_season_sq']])
        
        try:
            m_base = sm.GLM(y_train, X_base_train, family=sm.families.Poisson()).fit()
            pred_base = m_base.predict(X_base_test)
            
            resid_base = y_test - pred_base
            rmse_base = np.sqrt((resid_base**2).mean())
            mae_base = np.abs(resid_base).mean()
            
            # Poisson log-likelihood on test set
            ll_base = float(stats.poisson.logpmf(y_test.astype(int), pred_base.clip(0.01)).sum())
        except:
            continue
        
        # Enhanced: seasonal + SSW + strat temp + geomag
        enhanced_cols = ['day_of_season', 'day_of_season_sq', 'ssw_within_15d', 
                        'ncep_t_10hpa', 'post_event_1_3d']
        X_enh_train = sm.add_constant(train[enhanced_cols])
        X_enh_test = sm.add_constant(test[enhanced_cols])
        
        try:
            m_enh = sm.GLM(y_train, X_enh_train, family=sm.families.Poisson()).fit()
            pred_enh = m_enh.predict(X_enh_test)
            
            resid_enh = y_test - pred_enh
            rmse_enh = np.sqrt((resid_enh**2).mean())
            mae_enh = np.abs(resid_enh).mean()
            ll_enh = float(stats.poisson.logpmf(y_test.astype(int), pred_enh.clip(0.01)).sum())
        except:
            continue
        
        baseline_rmse.append(rmse_base)
        enhanced_rmse.append(rmse_enh)
        baseline_mae.append(mae_base)
        enhanced_mae.append(mae_enh)
        baseline_ll.append(ll_base)
        enhanced_ll.append(ll_enh)
    
    baseline_rmse = np.array(baseline_rmse)
    enhanced_rmse = np.array(enhanced_rmse)
    baseline_mae = np.array(baseline_mae)
    enhanced_mae = np.array(enhanced_mae)
    baseline_ll = np.array(baseline_ll)
    enhanced_ll = np.array(enhanced_ll)
    
    # Improvement
    rmse_diff = enhanced_rmse - baseline_rmse  # negative = improvement
    mae_diff = enhanced_mae - baseline_mae
    ll_diff = enhanced_ll - baseline_ll  # positive = improvement
    
    t_rmse, p_rmse = stats.ttest_rel(enhanced_rmse, baseline_rmse)
    t_mae, p_mae = stats.ttest_rel(enhanced_mae, baseline_mae)
    t_ll, p_ll = stats.ttest_rel(enhanced_ll, baseline_ll)
    
    n_rmse_improved = int(np.sum(rmse_diff < 0))
    n_ll_improved = int(np.sum(ll_diff > 0))
    
    print(f"\n  Leave-one-winter-out results ({len(baseline_rmse)} folds):")
    print(f"  RMSE: baseline={baseline_rmse.mean():.3f}, enhanced={enhanced_rmse.mean():.3f}")
    print(f"    {n_rmse_improved}/{len(rmse_diff)} folds improved, paired t P={p_rmse:.4f}")
    print(f"  MAE: baseline={baseline_mae.mean():.3f}, enhanced={enhanced_mae.mean():.3f}")
    print(f"    paired t P={p_mae:.4f}")
    print(f"  Log-lik: baseline={baseline_ll.mean():.1f}, enhanced={enhanced_ll.mean():.1f}")
    print(f"    {n_ll_improved}/{len(ll_diff)} folds improved, paired t P={p_ll:.4f}")
    
    # Full model coefficients
    print("\n  Full-sample enhanced model coefficients:")
    full_data = model_data.dropna()
    y_full = full_data['dry_natural_size_1234']
    enhanced_cols = ['day_of_season', 'day_of_season_sq', 'ssw_within_15d', 
                    'ncep_t_10hpa', 'post_event_1_3d']
    X_full = sm.add_constant(full_data[enhanced_cols])
    m_full = sm.GLM(y_full, X_full, family=sm.families.Poisson()).fit()
    
    for name, coef, pval in zip(m_full.params.index, m_full.params, m_full.pvalues):
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
        print(f"    {name}: {coef:.4f} (P={pval:.4f}) {sig}")
    
    results = {
        'n_folds': len(baseline_rmse),
        'rmse_baseline': float(baseline_rmse.mean()),
        'rmse_enhanced': float(enhanced_rmse.mean()),
        'rmse_improvement_pct': float((1 - enhanced_rmse.mean()/baseline_rmse.mean()) * 100),
        'rmse_paired_t_p': float(p_rmse),
        'n_rmse_improved': n_rmse_improved,
        'mae_baseline': float(baseline_mae.mean()),
        'mae_enhanced': float(enhanced_mae.mean()),
        'mae_paired_t_p': float(p_mae),
        'loglik_baseline': float(baseline_ll.mean()),
        'loglik_enhanced': float(enhanced_ll.mean()),
        'loglik_paired_t_p': float(p_ll),
        'n_loglik_improved': n_ll_improved,
        'full_model_coefficients': {
            name: {'coef': float(coef), 'pval': float(pval)}
            for name, coef, pval in zip(m_full.params.index, m_full.params, m_full.pvalues)
        }
    }
    
    with open(RESULTS / 'lim_part6_prediction.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\n✓ Part 6 saved")

print("\n" + "="*60)
print(f"PART {part} COMPLETE")
print("="*60)
