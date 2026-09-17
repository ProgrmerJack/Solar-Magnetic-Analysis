"""
Winter-blocked bootstrap for dependence-aware confidence intervals.

Reviewer concern: standard bootstrap ignores within-winter temporal
correlation. This script resamples ENTIRE winters (not individual events),
producing conservative CIs that respect the dependence structure.

Output: data/results/r44_winter_blocked_bootstrap.json
"""
import sys
from pathlib import Path
import json
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

RESULTS = ROOT / "data" / "results"
RESULTS.mkdir(exist_ok=True)

# All 16 SSW events with their rate ratios and winter assignments
EVENTS = [
    {"date": "1998-12-15", "rr": 0.36, "winter": "1998-99"},
    {"date": "1999-02-26", "rr": 0.21, "winter": "1998-99"},
    {"date": "2001-02-11", "rr": 0.42, "winter": "2000-01"},
    {"date": "2001-12-30", "rr": 0.39, "winter": "2001-02"},
    {"date": "2002-02-17", "rr": 0.28, "winter": "2001-02"},
    {"date": "2003-01-18", "rr": 0.18, "winter": "2002-03"},
    {"date": "2004-01-05", "rr": 0.32, "winter": "2003-04"},
    {"date": "2006-01-21", "rr": 0.15, "winter": "2005-06"},
    {"date": "2007-02-24", "rr": 0.22, "winter": "2006-07"},
    {"date": "2008-02-22", "rr": 0.19, "winter": "2007-08"},
    {"date": "2009-01-24", "rr": 0.38, "winter": "2008-09"},
    {"date": "2010-02-09", "rr": 0.44, "winter": "2009-10"},
    {"date": "2012-01-11", "rr": 0.28, "winter": "2011-12"},
    {"date": "2013-01-07", "rr": 0.31, "winter": "2012-13"},
    {"date": "2018-02-12", "rr": 0.85, "winter": "2017-18"},
    {"date": "2019-01-01", "rr": 3.43, "winter": "2018-19"},
]


def standard_bootstrap(log_rr_vals, n_boot=10000, seed=42):
    """Standard i.i.d. bootstrap (assumes independence)."""
    rng = np.random.RandomState(seed)
    n = len(log_rr_vals)
    boot_means = []
    for _ in range(n_boot):
        sample = rng.choice(log_rr_vals, size=n, replace=True)
        boot_means.append(np.mean(sample))
    return np.array(boot_means)


def winter_blocked_bootstrap(events, n_boot=10000, seed=42):
    """Block bootstrap resampling by winter season."""
    rng = np.random.RandomState(seed)
    
    # Group events by winter
    winters = {}
    for e in events:
        w = e['winter']
        if w not in winters:
            winters[w] = []
        winters[w].append(np.log(e['rr']))
    
    winter_keys = list(winters.keys())
    n_winters = len(winter_keys)
    
    boot_means = []
    for _ in range(n_boot):
        # Resample entire winters
        sampled_winters = rng.choice(winter_keys, size=n_winters, replace=True)
        all_log_rr = []
        for w in sampled_winters:
            all_log_rr.extend(winters[w])
        boot_means.append(np.mean(all_log_rr))
    
    return np.array(boot_means)


def compute_all():
    """Compute both standard and blocked bootstrap CIs."""
    log_rr = np.array([np.log(e['rr']) for e in EVENTS])
    
    print(f"N events: {len(EVENTS)}")
    print(f"Unique winters: {len(set(e['winter'] for e in EVENTS))}")
    print(f"Mean log(RR): {np.mean(log_rr):.4f}")
    print(f"Geometric mean RR: {np.exp(np.mean(log_rr)):.4f}")
    
    # Standard bootstrap
    std_boots = standard_bootstrap(log_rr)
    std_ci = np.percentile(std_boots, [2.5, 97.5])
    std_geo_ci = np.exp(std_ci)
    
    # Winter-blocked bootstrap
    block_boots = winter_blocked_bootstrap(EVENTS)
    block_ci = np.percentile(block_boots, [2.5, 97.5])
    block_geo_ci = np.exp(block_ci)
    
    # Parametric t-test (for comparison)
    t_stat, p_val = stats.ttest_1samp(log_rr, 0)
    t_ci = stats.t.interval(0.95, df=len(log_rr)-1, 
                            loc=np.mean(log_rr), scale=stats.sem(log_rr))
    t_geo_ci = np.exp(np.array(t_ci))
    
    # Sign test (non-parametric)
    n_decrease = sum(1 for r in log_rr if r < 0)
    sign_p = stats.binomtest(n_decrease, len(log_rr), 0.5).pvalue
    
    results = {
        "n_events": len(EVENTS),
        "n_winters": len(set(e['winter'] for e in EVENTS)),
        "n_multi_event_winters": sum(1 for w in set(e['winter'] for e in EVENTS) 
                                     if sum(1 for e2 in EVENTS if e2['winter'] == w) > 1),
        "mean_log_rr": float(np.mean(log_rr)),
        "geo_mean_rr": float(np.exp(np.mean(log_rr))),
        "standard_bootstrap": {
            "ci_95_log": [float(std_ci[0]), float(std_ci[1])],
            "ci_95_rr": [float(std_geo_ci[0]), float(std_geo_ci[1])],
            "se": float(np.std(std_boots)),
        },
        "winter_blocked_bootstrap": {
            "ci_95_log": [float(block_ci[0]), float(block_ci[1])],
            "ci_95_rr": [float(block_geo_ci[0]), float(block_geo_ci[1])],
            "se": float(np.std(block_boots)),
        },
        "parametric_t": {
            "t_stat": float(t_stat),
            "p_value": float(p_val),
            "ci_95_log": [float(t_ci[0]), float(t_ci[1])],
            "ci_95_rr": [float(t_geo_ci[0]), float(t_geo_ci[1])],
        },
        "sign_test": {
            "n_decrease": int(n_decrease),
            "n_total": len(log_rr),
            "fraction_decrease": float(n_decrease / len(log_rr)),
            "p_value": float(sign_p),
        },
        "ci_width_comparison": {
            "standard_width_rr": float(std_geo_ci[1] - std_geo_ci[0]),
            "blocked_width_rr": float(block_geo_ci[1] - block_geo_ci[0]),
            "ratio_blocked_to_standard": float((block_geo_ci[1] - block_geo_ci[0]) / 
                                                (std_geo_ci[1] - std_geo_ci[0])),
            "interpretation": "Blocked CI width / standard CI width. >1 means blocked is wider (more conservative)."
        }
    }
    
    print("\n=== Results ===")
    print(f"Standard bootstrap 95% CI (RR): [{std_geo_ci[0]:.3f}, {std_geo_ci[1]:.3f}]")
    print(f"Blocked bootstrap 95% CI (RR):  [{block_geo_ci[0]:.3f}, {block_geo_ci[1]:.3f}]")
    print(f"Parametric t CI (RR):           [{t_geo_ci[0]:.3f}, {t_geo_ci[1]:.3f}]")
    print(f"Blocked/Standard CI width ratio: {results['ci_width_comparison']['ratio_blocked_to_standard']:.2f}")
    print(f"Sign test: {n_decrease}/{len(log_rr)} decrease, P = {sign_p:.4f}")
    
    both_below_1 = block_geo_ci[1] < 1.0
    print(f"\nBlocked CI entirely below 1.0: {both_below_1}")
    results["blocked_ci_below_unity"] = bool(both_below_1)
    
    outpath = RESULTS / "r44_winter_blocked_bootstrap.json"
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {outpath}")
    
    return results


if __name__ == "__main__":
    compute_all()
