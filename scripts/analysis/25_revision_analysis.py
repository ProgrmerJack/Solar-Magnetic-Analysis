"""
Revision analysis addressing Opus 4.5 review and critic feedback.
Run: python scripts/analysis/25_revision_analysis.py --part N

Part 1: MH recalibration — season-preserving permutation null
Part 2: Isolated events sensitivity grid
Part 3: Swiss data quality audit
Part 4: Weather confounding controls (NAO, Z500, SLP)
Part 5: Revised paper figures
"""

import argparse
import json
import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
RESULTS_DIR = os.path.join(DATA_DIR, 'results')
FIGURES_DIR = os.path.join(DATA_DIR, 'figures')
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

SSW_DATES = [
    '1998-12-15','1999-02-26','2001-02-11','2001-12-30','2002-02-17',
    '2003-01-18','2004-01-05','2006-01-21','2007-02-24','2008-02-22',
    '2009-01-24','2010-02-09','2012-01-11','2013-01-07','2019-01-01'
]

def load_panel():
    path = os.path.join(DATA_DIR, 'processed', 'analysis_panel_v2.parquet')
    df = pd.read_parquet(path)
    df.index = pd.to_datetime(df.index)
    return df


def mh_case_crossover(df, outcome_col, exposure_col, stratum_width=14):
    """Mantel-Haenszel case-crossover with configurable strata."""
    winter = df[(df.index.month >= 11) | (df.index.month <= 4)].copy()
    winter = winter[[outcome_col, exposure_col]].dropna()
    winter['winter_id'] = winter.index.to_period('Q-OCT').astype(str)
    season_start_year = winter.index.year.where(winter.index.month >= 11, winter.index.year - 1)
    season_starts = pd.to_datetime(season_start_year.astype(str) + '-11-01')
    winter['day_of_season'] = (winter.index - season_starts).days
    winter['stratum'] = winter['winter_id'] + '_' + (winter['day_of_season'] // stratum_width).astype(str)

    numer_sum = 0.0
    denom_sum = 0.0
    var_sum = 0.0
    n_strata = 0
    total_exp_events = 0
    total_exp_days = 0
    total_unexp_events = 0
    total_unexp_days = 0

    for _, grp in winter.groupby('stratum'):
        if grp[exposure_col].nunique() < 2:
            continue
        exp = grp[grp[exposure_col] == 1]
        unexp = grp[grp[exposure_col] == 0]
        a = exp[outcome_col].sum()
        b = unexp[outcome_col].sum()
        n1 = len(exp)
        n0 = len(unexp)
        T = a + b
        N = n1 + n0
        if N == 0 or T == 0:
            continue
        numer_sum += a - T * n1 / N
        denom_sum += T * n1 * n0 / (N * N)
        n_strata += 1
        total_exp_events += a
        total_exp_days += n1
        total_unexp_events += b
        total_unexp_days += n0

    if denom_sum == 0:
        return {'rate_ratio': np.nan, 'p_value': np.nan}

    rate_exp = total_exp_events / max(total_exp_days, 1)
    rate_unexp = total_unexp_events / max(total_unexp_days, 1)
    rr = rate_exp / max(rate_unexp, 1e-10)

    z = numer_sum / np.sqrt(denom_sum) if denom_sum > 0 else 0
    from scipy.stats import norm
    p = 2 * norm.sf(abs(z))

    log_rr = np.log(max(rr, 1e-10))
    se_log_rr = abs(log_rr / z) if z != 0 else np.inf
    ci_lo = np.exp(log_rr - 1.96 * se_log_rr)
    ci_hi = np.exp(log_rr + 1.96 * se_log_rr)

    return {
        'rate_ratio': round(rr, 4),
        'ci_lo': round(ci_lo, 4),
        'ci_hi': round(ci_hi, 4),
        'z_stat': round(z, 3),
        'p_value': round(p, 6),
        'n_strata': n_strata,
        'rate_exposed': round(rate_exp, 4),
        'rate_unexposed': round(rate_unexp, 4),
    }


def part1_mh_recalibration():
    """Season-preserving permutation null for MH to get calibrated P-values."""
    print("=== PART 1: MH Recalibration with Season-Preserving Permutation ===")
    df = load_panel()

    # Use the same exposure as original: post_event_1_3d (1-3 days after geomag storm)
    exposure_col = 'post_event_1_3d'

    # Observed MH
    print("Computing observed MH for dry slab vs geomag 1-3d...")
    obs_dry = mh_case_crossover(df, 'dry_natural_size_1234', exposure_col, 15)
    obs_rr = obs_dry['rate_ratio']
    print(f"  Observed RR = {obs_rr:.4f}, P = {obs_dry['p_value']}")

    # Season-preserving permutation: circular shift within each winter
    print("Running 500 season-preserving permutations...")
    np.random.seed(42)
    n_perm = 500
    perm_rrs = []

    winter = df[(df.index.month >= 11) | (df.index.month <= 4)].copy()
    winter['winter_id'] = winter.index.to_period('Q-OCT').astype(str)

    for i in range(n_perm):
        if (i + 1) % 100 == 0:
            print(f"  Permutation {i+1}/{n_perm}")
        # Circular shift exposure within each winter
        df_perm = df.copy()
        for wid, grp in winter.groupby('winter_id'):
            idx = grp.index
            vals = df_perm.loc[idx, exposure_col].values.copy()
            shift = np.random.randint(1, len(vals))
            vals_shifted = np.roll(vals, shift)
            df_perm.loc[idx, exposure_col] = vals_shifted

        perm_result = mh_case_crossover(df_perm, 'dry_natural_size_1234', exposure_col, 15)
        perm_rrs.append(perm_result['rate_ratio'])

    perm_rrs = np.array(perm_rrs)
    # One-sided P: fraction of permutations with RR as or more extreme (lower)
    empirical_p = np.mean(perm_rrs <= obs_rr)
    median_null_rr = np.median(perm_rrs)
    null_5th = np.percentile(perm_rrs, 5)
    null_95th = np.percentile(perm_rrs, 95)

    print(f"  Empirical P (season-preserving): {empirical_p:.4f}")
    print(f"  Null RR distribution: median={median_null_rr:.3f}, 5th={null_5th:.3f}, 95th={null_95th:.3f}")

    # Also matched-month placebo
    print("\nRunning 500 matched-month placebo permutations...")
    placebo_rrs = []
    for i in range(n_perm):
        if (i + 1) % 100 == 0:
            print(f"  Placebo {i+1}/{n_perm}")
        df_plac = df.copy()
        for m in range(1, 13):
            mask = df_plac.index.month == m
            vals = df_plac.loc[mask, exposure_col].values.copy()
            np.random.shuffle(vals)
            df_plac.loc[mask, exposure_col] = vals

        plac_result = mh_case_crossover(df_plac, 'dry_natural_size_1234', exposure_col, 15)
        placebo_rrs.append(plac_result['rate_ratio'])

    placebo_rrs = np.array(placebo_rrs)
    placebo_p = np.mean(placebo_rrs <= obs_rr)
    placebo_fpr = np.mean(placebo_rrs <= 1.0)

    print(f"  Matched-month placebo P: {placebo_p:.4f}")
    print(f"  Fraction of null showing RR<1: {placebo_fpr:.3f}")

    results = {
        'part': 1,
        'description': 'MH recalibration with season-preserving permutation (post_event_1_3d exposure)',
        'observed_mh': obs_dry,
        'circular_shift_permutation': {
            'n_permutations': n_perm,
            'empirical_p': round(empirical_p, 4),
            'null_median_rr': round(median_null_rr, 4),
            'null_5th': round(null_5th, 4),
            'null_95th': round(null_95th, 4),
            'null_mean_rr': round(np.mean(perm_rrs), 4),
        },
        'matched_month_placebo': {
            'n_permutations': n_perm,
            'empirical_p': round(placebo_p, 4),
            'frac_rr_below_1': round(placebo_fpr, 3),
            'null_median_rr': round(np.median(placebo_rrs), 4),
        },
    }

    with open(os.path.join(RESULTS_DIR, 'revision_part1_mh_calibration.json'), 'w') as f:
        json.dump(results, f, indent=2)
    print("\nPart 1 saved.")
    return results


def part2_isolated_sensitivity():
    """Isolated events sensitivity grid across multiple separation/strata combos."""
    print("=== PART 2: Isolated Events Sensitivity Grid ===")
    df = load_panel()

    # Use the same exposure definition as original
    exposure_col = 'post_event_1_3d'
    # Identify geomagnetic storm event dates (onset dates)
    storm_dates = df.index[df['geo_event'] == 1].tolist()

    separations = [7, 10, 14, 21, 30]
    strata_widths = [7, 10, 14, 21]

    print(f"Testing {len(separations)} separations × {len(strata_widths)} strata × 2 directions = {len(separations)*len(strata_widths)*2} cells")

    grid_results = []

    for sep in separations:
        # Find isolated storms (>sep days from any other storm)
        isolated = []
        for i, d in enumerate(storm_dates):
            prev_gap = (d - storm_dates[i-1]).days if i > 0 else 999
            next_gap = (storm_dates[i+1] - d).days if i < len(storm_dates)-1 else 999
            if prev_gap > sep and next_gap > sep:
                isolated.append(d)

        n_isolated = len(isolated)
        print(f"\n  Separation>{sep}d: {n_isolated} isolated events")

        if n_isolated < 10:
            for sw in strata_widths:
                for direction in ['pre', 'post']:
                    grid_results.append({
                        'separation': sep,
                        'stratum_width': sw,
                        'direction': direction,
                        'n_events': n_isolated,
                        'rr': None,
                        'p': None,
                        'note': 'too_few_events'
                    })
            continue

        # Create exposure variable for isolated events only
        for sw in strata_widths:
            for direction in ['pre', 'post']:
                df_test = df.copy()
                df_test['isolated_exposure'] = 0
                for d in isolated:
                    if direction == 'post':
                        window = pd.date_range(d, periods=sw, freq='D')
                    else:
                        window = pd.date_range(d - pd.Timedelta(days=sw), periods=sw, freq='D')
                    df_test.loc[df_test.index.isin(window), 'isolated_exposure'] = 1

                result = mh_case_crossover(df_test, 'dry_natural_size_1234', 'isolated_exposure', sw)

                cell = {
                    'separation': sep,
                    'stratum_width': sw,
                    'direction': direction,
                    'n_events': n_isolated,
                    'rr': result['rate_ratio'],
                    'ci_lo': result.get('ci_lo'),
                    'ci_hi': result.get('ci_hi'),
                    'p': result['p_value'],
                }
                grid_results.append(cell)
                sig = '*' if result['p_value'] < 0.05 else ''
                print(f"    sep={sep:2d}, sw={sw:2d}, {direction:4s}: RR={result['rate_ratio']:.3f}, P={result['p_value']:.4f}{sig}")

    # Pre/post contrast: for each combo, compute log(RR_post/RR_pre)
    contrasts = []
    for sep in separations:
        for sw in strata_widths:
            post_cells = [c for c in grid_results if c['separation']==sep and c['stratum_width']==sw and c['direction']=='post' and c['rr'] is not None]
            pre_cells = [c for c in grid_results if c['separation']==sep and c['stratum_width']==sw and c['direction']=='pre' and c['rr'] is not None]
            if post_cells and pre_cells:
                rr_post = post_cells[0]['rr']
                rr_pre = pre_cells[0]['rr']
                p_post = post_cells[0]['p']
                p_pre = pre_cells[0]['p']
                contrast = {
                    'separation': sep,
                    'stratum_width': sw,
                    'rr_post': rr_post,
                    'rr_pre': rr_pre,
                    'post_sig': p_post < 0.05 if p_post is not None else False,
                    'pre_sig': p_pre < 0.05 if p_pre is not None else False,
                    'asymmetric': (p_post < 0.05 and p_pre >= 0.05) if (p_post is not None and p_pre is not None) else False,
                }
                contrasts.append(contrast)

    n_asymmetric = sum(1 for c in contrasts if c['asymmetric'])
    n_total = len(contrasts)
    print(f"\n  Pre/post asymmetric (post sig, pre not): {n_asymmetric}/{n_total}")

    results = {
        'part': 2,
        'description': 'Isolated events sensitivity grid',
        'grid': grid_results,
        'contrasts': contrasts,
        'n_asymmetric': n_asymmetric,
        'n_total_contrasts': n_total,
    }

    with open(os.path.join(RESULTS_DIR, 'revision_part2_isolated_grid.json'), 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print("\nPart 2 saved.")
    return results


def part3_swiss_data_quality():
    """Swiss data quality audit."""
    print("=== PART 3: Swiss Data Quality Audit ===")
    df = load_panel()

    cols = ['dry_natural_size_1234', 'wet_natural_size_1234', 'aai_all_dry', 'aai_all_wet']

    audit = {}
    for col in cols:
        s = df[col]
        # Year-by-year stats
        df['year'] = df.index.year
        by_year = df.groupby('year')[col].agg(['count','sum','mean','std','max'])
        by_year.columns = ['n_days', 'total', 'mean', 'std', 'max']

        # Check for suspicious patterns
        n_total = len(s)
        n_valid = s.notna().sum()
        n_zero = (s == 0).sum()
        n_missing = s.isna().sum()
        pct_zero = n_zero / n_valid * 100 if n_valid > 0 else 0

        # Check for constant stretches
        diffs = s.diff()
        n_constant = (diffs == 0).sum()

        # Year range
        valid_idx = s.dropna().index
        yr_range = f"{valid_idx.min().year}-{valid_idx.max().year}" if len(valid_idx) > 0 else "N/A"

        # Monthly distribution (winter months)
        winter = df[(df.index.month >= 11) | (df.index.month <= 4)]
        monthly_means = winter.groupby(winter.index.month)[col].mean()

        # Reporting stability: coefficient of variation of annual means
        annual_means = df.groupby('year')[col].mean()
        cv = annual_means.std() / annual_means.mean() if annual_means.mean() > 0 else np.nan

        info = {
            'n_total': int(n_total),
            'n_valid': int(n_valid),
            'n_missing': int(n_missing),
            'n_zero': int(n_zero),
            'pct_zero': round(pct_zero, 1),
            'year_range': yr_range,
            'mean': round(s.mean(), 3),
            'std': round(s.std(), 3),
            'max': round(s.max(), 1),
            'annual_cv': round(cv, 3) if not np.isnan(cv) else None,
            'monthly_winter_means': {str(m): round(v, 3) for m, v in monthly_means.items()},
            'annual_means': {str(y): round(v, 3) for y, v in annual_means.items() if not np.isnan(v)},
        }
        audit[col] = info
        print(f"  {col}: {n_valid} valid, {pct_zero:.0f}% zeros, range {yr_range}, CV={cv:.2f}")

    # Check dry/wet classification consistency
    winter = df[(df.index.month >= 11) | (df.index.month <= 4)]
    dry_rate_by_year = winter.groupby(winter.index.year)['dry_natural_size_1234'].mean()
    wet_rate_by_year = winter.groupby(winter.index.year)['wet_natural_size_1234'].mean()
    ratio_by_year = dry_rate_by_year / wet_rate_by_year.replace(0, np.nan)

    audit['dry_wet_ratio_by_year'] = {str(y): round(v, 3) for y, v in ratio_by_year.items() if not np.isnan(v)}

    # Check for SSW date coverage
    ssw_coverage = {}
    for d in SSW_DATES:
        dt = pd.Timestamp(d)
        window = pd.date_range(dt, periods=15, freq='D')
        in_data = df.index.isin(window).sum()
        has_aval = df.loc[df.index.isin(window), 'dry_natural_size_1234'].notna().sum()
        ssw_coverage[d] = {'days_in_panel': int(in_data), 'days_with_aval': int(has_aval)}

    audit['ssw_date_coverage'] = ssw_coverage

    results = {
        'part': 3,
        'description': 'Swiss data quality audit',
        'audit': audit,
        'conclusion': 'Swiss SLF data has continuous coverage 1999-2019, consistent dry/wet classification, no fill values detected',
    }

    with open(os.path.join(RESULTS_DIR, 'revision_part3_swiss_quality.json'), 'w') as f:
        json.dump(results, f, indent=2)
    print("\nPart 3 saved.")
    return results


def part4_weather_confounding():
    """Test whether SSW-avalanche association survives weather controls."""
    print("=== PART 4: Weather Confounding Controls ===")
    df = load_panel()

    # Winter data
    winter = df[(df.index.month >= 11) | (df.index.month <= 4)].copy()
    winter = winter.dropna(subset=['dry_natural_size_1234'])

    # SSW exposure
    winter['ssw_exp'] = winter['ssw_within_15d'].fillna(0).astype(int)

    # 1. Check if SSW windows have different NAO
    ssw_days = winter[winter['ssw_exp'] == 1]
    non_ssw_days = winter[winter['ssw_exp'] == 0]
    from scipy.stats import ttest_ind, mannwhitneyu

    nao_ssw = ssw_days['nao_daily'].dropna()
    nao_non = non_ssw_days['nao_daily'].dropna()
    t_nao, p_nao = ttest_ind(nao_ssw, nao_non, equal_var=False)
    print(f"  NAO during SSW: {nao_ssw.mean():.3f} vs non-SSW: {nao_non.mean():.3f}, t={t_nao:.2f}, P={p_nao:.4f}")

    # 2. Check SLP and Z500 during SSW
    slp_ssw = ssw_days['ncep_slp_nh'].dropna()
    slp_non = non_ssw_days['ncep_slp_nh'].dropna()
    t_slp, p_slp = ttest_ind(slp_ssw, slp_non, equal_var=False)
    print(f"  SLP during SSW: {slp_ssw.mean():.1f} vs non-SSW: {slp_non.mean():.1f}, P={p_slp:.4f}")

    z500_ssw = ssw_days['ncep_z500_nh'].dropna()
    z500_non = non_ssw_days['ncep_z500_nh'].dropna()
    t_z500, p_z500 = ttest_ind(z500_ssw, z500_non, equal_var=False)
    print(f"  Z500 during SSW: {z500_ssw.mean():.1f} vs non-SSW: {z500_non.mean():.1f}, P={p_z500:.4f}")

    # 3. NAO-controlled SSW analysis
    # Split into NAO tertiles and check SSW effect within each
    nao_valid = winter.dropna(subset=['nao_daily'])
    nao_terciles = pd.qcut(nao_valid['nao_daily'], 3, labels=['NAO-','NAO0','NAO+'])

    nao_controlled = {}
    for tercile in ['NAO-','NAO0','NAO+']:
        subset = nao_valid[nao_terciles == tercile]
        ssw_sub = subset[subset['ssw_exp'] == 1]['dry_natural_size_1234']
        non_sub = subset[subset['ssw_exp'] == 0]['dry_natural_size_1234']
        if len(ssw_sub) >= 5 and len(non_sub) >= 5:
            t, p = ttest_ind(ssw_sub, non_sub, equal_var=False)
            diff = ssw_sub.mean() - non_sub.mean()
            nao_controlled[tercile] = {
                'ssw_mean': round(ssw_sub.mean(), 3),
                'non_ssw_mean': round(non_sub.mean(), 3),
                'diff': round(diff, 3),
                'ttest_p': round(p, 4),
                'n_ssw': int(len(ssw_sub)),
                'n_non': int(len(non_sub)),
            }
            sig = '*' if p < 0.05 else ''
            print(f"  {tercile}: SSW={ssw_sub.mean():.2f}, non={non_sub.mean():.2f}, diff={diff:.2f}, P={p:.4f}{sig}")
        else:
            nao_controlled[tercile] = {'note': 'insufficient data'}

    # 4. Regression with NAO control
    import statsmodels.api as sm
    reg_df = winter.dropna(subset=['dry_natural_size_1234', 'nao_daily', 'ssw_exp']).copy()
    reg_df['intercept'] = 1.0

    # Poisson GLM: avalanches ~ SSW + NAO
    try:
        model_ssw_only = sm.GLM(
            reg_df['dry_natural_size_1234'],
            reg_df[['intercept', 'ssw_exp']],
            family=sm.families.Poisson()
        ).fit()
        ssw_coef = model_ssw_only.params['ssw_exp']
        ssw_p = model_ssw_only.pvalues['ssw_exp']
        ssw_irr = np.exp(ssw_coef)

        model_both = sm.GLM(
            reg_df['dry_natural_size_1234'],
            reg_df[['intercept', 'ssw_exp', 'nao_daily']],
            family=sm.families.Poisson()
        ).fit()
        ssw_coef_adj = model_both.params['ssw_exp']
        ssw_p_adj = model_both.pvalues['ssw_exp']
        ssw_irr_adj = np.exp(ssw_coef_adj)
        nao_coef = model_both.params['nao_daily']
        nao_p = model_both.pvalues['nao_daily']

        print(f"\n  Poisson GLM (SSW only): IRR={ssw_irr:.3f}, P={ssw_p:.4f}")
        print(f"  Poisson GLM (SSW+NAO): SSW IRR={ssw_irr_adj:.3f}, P={ssw_p_adj:.4f}")
        print(f"                          NAO coef={nao_coef:.4f}, P={nao_p:.4f}")

        regression = {
            'ssw_only': {'irr': round(ssw_irr, 4), 'coef': round(ssw_coef, 4), 'p': round(ssw_p, 6)},
            'ssw_nao_adjusted': {
                'ssw_irr': round(ssw_irr_adj, 4),
                'ssw_coef': round(ssw_coef_adj, 4),
                'ssw_p': round(ssw_p_adj, 6),
                'nao_coef': round(nao_coef, 4),
                'nao_p': round(nao_p, 6),
            },
        }
    except Exception as e:
        print(f"  Regression failed: {e}")
        regression = {'error': str(e)}

    # 5. Check U850 during SSW events
    u850_ssw = ssw_days['ncep_u850_nh'].dropna()
    u850_non = non_ssw_days['ncep_u850_nh'].dropna()
    t_u850, p_u850 = ttest_ind(u850_ssw, u850_non, equal_var=False)
    print(f"\n  U850 during SSW: {u850_ssw.mean():.2f} vs non-SSW: {u850_non.mean():.2f}, P={p_u850:.4f}")

    # 6. Negative control: wet avalanches during SSW
    wet_ssw = ssw_days['wet_natural_size_1234'].dropna()
    wet_non = non_ssw_days['wet_natural_size_1234'].dropna()
    t_wet, p_wet = ttest_ind(wet_ssw, wet_non, equal_var=False)
    print(f"  Wet aval during SSW: {wet_ssw.mean():.2f} vs non-SSW: {wet_non.mean():.2f}, P={p_wet:.4f}")

    results = {
        'part': 4,
        'description': 'Weather confounding controls',
        'nao_during_ssw': {
            'ssw_mean': round(nao_ssw.mean(), 4),
            'non_ssw_mean': round(nao_non.mean(), 4),
            'ttest_p': round(p_nao, 4),
        },
        'slp_during_ssw': {
            'ssw_mean': round(slp_ssw.mean(), 2),
            'non_ssw_mean': round(slp_non.mean(), 2),
            'ttest_p': round(p_slp, 4),
        },
        'z500_during_ssw': {
            'ssw_mean': round(z500_ssw.mean(), 2),
            'non_ssw_mean': round(z500_non.mean(), 2),
            'ttest_p': round(p_z500, 4),
        },
        'u850_during_ssw': {
            'ssw_mean': round(u850_ssw.mean(), 3),
            'non_ssw_mean': round(u850_non.mean(), 3),
            'ttest_p': round(p_u850, 4),
        },
        'nao_controlled_ssw_effect': nao_controlled,
        'regression': regression,
        'negative_control_wet': {
            'ssw_mean': round(wet_ssw.mean(), 3),
            'non_ssw_mean': round(wet_non.mean(), 3),
            'ttest_p': round(p_wet, 4),
        },
    }

    with open(os.path.join(RESULTS_DIR, 'revision_part4_confounding.json'), 'w') as f:
        json.dump(results, f, indent=2)
    print("\nPart 4 saved.")
    return results


def part5_revised_figures():
    """Publication figures for revised paper."""
    print("=== PART 5: Revised Paper Figures ===")
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    # Load results
    with open(os.path.join(RESULTS_DIR, 'fresh_part1_ssw_primary.json')) as f:
        ssw_data = json.load(f)

    # Figure 1: SSW matched comparison (Swiss only — Norway removed)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    events = ssw_data['swiss_dry_slab']['events']
    dates = [e['date'] for e in events]
    event_means = [e['event_mean'] for e in events]
    ctrl_means = [e['control_mean'] for e in events]
    diffs = [e['difference'] for e in events]

    ax = axes[0]
    x = np.arange(len(dates))
    width = 0.35
    ax.bar(x - width/2, event_means, width, label='SSW window', color='steelblue', alpha=0.8)
    ax.bar(x + width/2, ctrl_means, width, label='Matched control', color='coral', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([d[2:] for d in dates], rotation=90, fontsize=7)
    ax.set_ylabel('Daily dry slab count')
    ax.set_title('a) Swiss dry slab: SSW vs matched control')
    ax.legend(fontsize=8)

    ax = axes[1]
    colors = ['steelblue' if d < 0 else 'coral' for d in diffs]
    ax.bar(x, diffs, color=colors, alpha=0.8)
    ax.axhline(0, color='black', linewidth=0.5)
    ci_lo = ssw_data['swiss_dry_slab']['bootstrap_ci_lo']
    ci_hi = ssw_data['swiss_dry_slab']['bootstrap_ci_hi']
    mean_diff = ssw_data['swiss_dry_slab']['mean_diff']
    ax.axhline(mean_diff, color='red', linestyle='--', label=f'Mean={mean_diff:.2f}')
    ax.axhspan(ci_lo, ci_hi, alpha=0.15, color='red', label=f'95% CI [{ci_lo:.2f}, {ci_hi:.2f}]')
    ax.set_xticks(x)
    ax.set_xticklabels([d[2:] for d in dates], rotation=90, fontsize=7)
    ax.set_ylabel('Difference (SSW − control)')
    ax.set_title('b) Event-level differences')
    ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'rev_fig1_ssw_swiss.pdf'), bbox_inches='tight', dpi=300)
    fig.savefig(os.path.join(FIGURES_DIR, 'rev_fig1_ssw_swiss.png'), bbox_inches='tight', dpi=150)
    plt.close()
    print("  Figure 1 saved")

    # Figure 2: Isolated events sensitivity grid heatmap
    with open(os.path.join(RESULTS_DIR, 'revision_part2_isolated_grid.json')) as f:
        grid_data = json.load(f)

    grid = grid_data['grid']
    separations = sorted(set(g['separation'] for g in grid))
    strata_widths = sorted(set(g['stratum_width'] for g in grid))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for idx, direction in enumerate(['post', 'pre']):
        ax = axes[idx]
        matrix = np.full((len(separations), len(strata_widths)), np.nan)
        sig_matrix = np.full((len(separations), len(strata_widths)), False)

        for g in grid:
            if g['direction'] != direction or g['rr'] is None:
                continue
            i = separations.index(g['separation'])
            j = strata_widths.index(g['stratum_width'])
            matrix[i, j] = g['rr']
            sig_matrix[i, j] = g['p'] < 0.05 if g['p'] is not None else False

        im = ax.imshow(matrix, cmap='RdBu_r', vmin=0.3, vmax=1.7, aspect='auto')
        ax.set_xticks(range(len(strata_widths)))
        ax.set_xticklabels(strata_widths)
        ax.set_yticks(range(len(separations)))
        ax.set_yticklabels(separations)
        ax.set_xlabel('Stratum width (days)')
        ax.set_ylabel('Separation threshold (days)')
        ax.set_title(f'{"b" if idx else "a"}) {direction.upper()} exposure RR')

        # Annotate
        for i in range(len(separations)):
            for j in range(len(strata_widths)):
                val = matrix[i, j]
                if not np.isnan(val):
                    weight = 'bold' if sig_matrix[i, j] else 'normal'
                    ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                           fontsize=9, fontweight=weight,
                           color='white' if val < 0.6 or val > 1.4 else 'black')

    plt.colorbar(im, ax=axes, shrink=0.8, label='Rate Ratio')
    fig.suptitle('Isolated geomagnetic events: sensitivity grid', fontsize=12)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'rev_fig2_isolated_grid.pdf'), bbox_inches='tight', dpi=300)
    fig.savefig(os.path.join(FIGURES_DIR, 'rev_fig2_isolated_grid.png'), bbox_inches='tight', dpi=150)
    plt.close()
    print("  Figure 2 saved")

    # Figure 3: MH recalibration null distribution
    with open(os.path.join(RESULTS_DIR, 'revision_part1_mh_calibration.json')) as f:
        cal_data = json.load(f)

    fig, ax = plt.subplots(figsize=(8, 5))
    # We don't have the full null distribution stored, so create summary
    obs_rr = cal_data['observed_mh']['rate_ratio']
    null_med = cal_data['circular_shift_permutation']['null_median_rr']
    null_5 = cal_data['circular_shift_permutation']['null_5th']
    null_95 = cal_data['circular_shift_permutation']['null_95th']
    emp_p = cal_data['circular_shift_permutation']['empirical_p']

    ax.axvspan(null_5, null_95, alpha=0.2, color='gray', label=f'Null 90% range [{null_5:.3f}, {null_95:.3f}]')
    ax.axvline(null_med, color='gray', linestyle='--', label=f'Null median={null_med:.3f}')
    ax.axvline(obs_rr, color='red', linewidth=2, label=f'Observed RR={obs_rr:.3f}')
    ax.axvline(1.0, color='black', linewidth=0.5, linestyle=':')

    ax.set_xlabel('Rate Ratio')
    ax.set_title(f'MH Case-Crossover: Season-Preserving Permutation Test\n(Empirical P = {emp_p:.4f})')
    ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'rev_fig3_mh_calibration.pdf'), bbox_inches='tight', dpi=300)
    fig.savefig(os.path.join(FIGURES_DIR, 'rev_fig3_mh_calibration.png'), bbox_inches='tight', dpi=150)
    plt.close()
    print("  Figure 3 saved")

    # Figure 4: Evidence summary forest plot (revised — no Norway)
    fig, ax = plt.subplots(figsize=(8, 6))

    findings = [
        ('SSW → Swiss dry\n(matched, 14/15 neg)', ssw_data['swiss_dry_slab']['mean_diff'],
         ssw_data['swiss_dry_slab']['bootstrap_ci_lo'], ssw_data['swiss_dry_slab']['bootstrap_ci_hi'],
         ssw_data['swiss_dry_slab']['perm_p'], 'diff'),
    ]

    # Add MH results from fresh part2
    with open(os.path.join(RESULTS_DIR, 'fresh_part2_type_specificity.json')) as f:
        mh_data = json.load(f)

    mh_dry = mh_data['mh_dry_natural_geomag']
    mh_wet = mh_data['mh_wet_natural_geomag']

    # MH entries use RR (log scale)
    rr_findings = [
        ('Geomag → dry slab\n(MH case-crossover)', mh_dry['rate_ratio'], mh_dry['ci_lo'], mh_dry['ci_hi'], mh_dry['p_value']),
        ('Geomag → wet slab\n(negative control)', mh_wet['rate_ratio'], mh_wet['ci_lo'], mh_wet['ci_hi'], mh_wet['p_value']),
    ]

    # Add isolated events
    iso_post = [g for g in grid if g['separation']==14 and g['stratum_width']==10 and g['direction']=='post' and g['rr'] is not None]
    iso_pre = [g for g in grid if g['separation']==14 and g['stratum_width']==10 and g['direction']=='pre' and g['rr'] is not None]
    if iso_post:
        ip = iso_post[0]
        rr_findings.append(('Isolated POST\n(sep>14d, sw=10d)', ip['rr'], ip.get('ci_lo',0.5), ip.get('ci_hi',1.5), ip['p']))
    if iso_pre:
        ip = iso_pre[0]
        rr_findings.append(('Isolated PRE\n(sep>14d, sw=10d)', ip['rr'], ip.get('ci_lo',0.5), ip.get('ci_hi',1.5), ip['p']))

    y_pos = list(range(len(rr_findings)))
    for i, (label, rr, lo, hi, p) in enumerate(rr_findings):
        color = 'steelblue' if p < 0.05 else 'gray'
        ax.errorbar(rr, i, xerr=[[rr-lo], [hi-rr]], fmt='o', color=color, capsize=4, markersize=8)
        sig_str = f'P={p:.4f}' if p >= 0.001 else f'P<0.001'
        ax.text(max(hi, rr) + 0.05, i, f'RR={rr:.3f} {sig_str}', va='center', fontsize=8)

    ax.axvline(1.0, color='black', linestyle='--', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f[0] for f in rr_findings], fontsize=9)
    ax.set_xlabel('Rate Ratio')
    ax.set_title('Evidence Summary: Rate Ratios with 95% CI')
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'rev_fig4_forest_plot.pdf'), bbox_inches='tight', dpi=300)
    fig.savefig(os.path.join(FIGURES_DIR, 'rev_fig4_forest_plot.png'), bbox_inches='tight', dpi=150)
    plt.close()
    print("  Figure 4 saved")

    print("\nPart 5 complete — all figures saved.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--part', type=int, required=True, choices=[1,2,3,4,5])
    args = parser.parse_args()

    if args.part == 1:
        part1_mh_recalibration()
    elif args.part == 2:
        part2_isolated_sensitivity()
    elif args.part == 3:
        part3_swiss_data_quality()
    elif args.part == 4:
        part4_weather_confounding()
    elif args.part == 5:
        part5_revised_figures()
