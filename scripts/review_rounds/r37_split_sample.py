"""
Split-sample validation: divide 16 Swiss SSW events chronologically into
training (first 8: 1998-2006) and validation (last 8: 2007-2019).
If both halves show the effect, this provides a form of temporal replication
that addresses concerns about analytical flexibility.
"""
import pandas as pd, numpy as np, json
from scipy import stats
from scipy.stats import binomtest

panel = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
ssw_cat = pd.read_parquet('data/processed/atmospheric/ssw_catalog.parquet')
ssw_dates = ssw_cat.index.tz_localize(None)

# Get SSW events in panel range
ssw_in = ssw_dates[(ssw_dates >= panel.index.min()) & (ssw_dates <= panel.index.max())]
print(f"Total SSW events in panel: {len(ssw_in)}")

# Compute RR for each event
events = []
for d in ssw_in:
    w = (panel.index >= d - pd.Timedelta(days=15)) & (panel.index <= d + pd.Timedelta(days=15))
    obs = panel.loc[w, 'dry_natural_size_1234'].mean()
    
    doy = d.dayofyear
    ctrl_mask = (panel.index.dayofyear >= doy - 3) & (panel.index.dayofyear <= doy + 3) & ~w
    exp = panel.loc[ctrl_mask, 'dry_natural_size_1234'].mean()
    
    rr = obs / max(exp, 0.01)
    events.append({
        'onset': d,
        'obs': obs,
        'exp': exp,
        'rr': rr,
        'log_rr': np.log(rr),
        'decrease': rr < 1,
    })

edf = pd.DataFrame(events).sort_values('onset')
print(f"\nAll {len(edf)} events:")
for _, r in edf.iterrows():
    sign = "v" if r['decrease'] else "^"
    print(f"  {r['onset'].date()} RR={r['rr']:.3f} log(RR)={r['log_rr']:+.2f} {sign}")

# Split chronologically
n = len(edf)
half = n // 2
train = edf.iloc[:half]
test = edf.iloc[half:]

print(f"\n{'='*60}")
print(f"SPLIT-SAMPLE VALIDATION")
print(f"{'='*60}")

for name, subset in [('Training (first 8: 1998-2006)', train), 
                      ('Validation (last 8: 2007-2019)', test),
                      ('Full sample (n=16)', edf)]:
    n_events = len(subset)
    n_decrease = subset['decrease'].sum()
    gm_rr = np.exp(subset['log_rr'].mean())
    d = subset['log_rr'].mean() / subset['log_rr'].std()
    
    # Sign test
    p_sign = binomtest(n_decrease, n_events, 0.5, alternative='greater').pvalue
    
    # t-test
    t_stat, p_t = stats.ttest_1samp(subset['log_rr'], 0)
    
    # Wilcoxon
    try:
        w_stat, p_w = stats.wilcoxon(subset['log_rr'], alternative='less')
    except:
        p_w = np.nan
    
    # Bootstrap CI
    np.random.seed(42)
    boot_means = [np.exp(np.random.choice(subset['log_rr'].values, n_events, replace=True).mean()) for _ in range(10000)]
    ci_lo, ci_hi = np.percentile(boot_means, [2.5, 97.5])
    
    print(f"\n--- {name} ---")
    print(f"  Events: {n_events}")
    print(f"  Decrease: {n_decrease}/{n_events} ({n_decrease/n_events*100:.1f}%)")
    print(f"  Geometric mean RR: {gm_rr:.3f} [95% CI: {ci_lo:.3f}, {ci_hi:.3f}]")
    print(f"  Cohen's d: {d:.3f}")
    print(f"  Sign test P: {p_sign:.4f}")
    print(f"  t-test P: {p_t:.4f}")
    print(f"  Wilcoxon P: {p_w:.4f}")

# Cross-validation: train on first half, predict second half
print(f"\n{'='*60}")
print(f"CROSS-TEMPORAL PREDICTION")
print(f"{'='*60}")

# Train: learn that SSW → decrease
train_direction = (train['decrease'].sum() / len(train))
print(f"Training base rate: {train_direction:.3f} decrease")
print(f"Validation observed: {test['decrease'].sum()}/{len(test)} decrease")

# Predict "decrease" for all validation events
correct = test['decrease'].sum()  # All decrease predictions
accuracy = correct / len(test)
print(f"Prediction accuracy: {correct}/{len(test)} = {accuracy:.3f}")

# Is the validation effect consistent with training?
# Two-sample test of log(RR) between halves
t_2s, p_2s = stats.ttest_ind(train['log_rr'], test['log_rr'])
mw_2s, p_mw = stats.mannwhitneyu(train['log_rr'], test['log_rr'])
print(f"\nTraining vs validation difference:")
print(f"  Training mean log(RR): {train['log_rr'].mean():.3f}")
print(f"  Validation mean log(RR): {test['log_rr'].mean():.3f}")
print(f"  Two-sample t-test P: {p_2s:.3f}")
print(f"  Mann-Whitney P: {p_mw:.3f}")
print(f"  Interpretation: {'NO significant difference' if p_2s > 0.05 else 'Significant difference'}")

# Save results
results = {
    'training': {
        'n': int(len(train)),
        'n_decrease': int(train['decrease'].sum()),
        'gm_rr': float(np.exp(train['log_rr'].mean())),
        'd': float(train['log_rr'].mean() / train['log_rr'].std()),
        'p_sign': float(binomtest(int(train['decrease'].sum()), len(train), 0.5, alternative='greater').pvalue),
    },
    'validation': {
        'n': int(len(test)),
        'n_decrease': int(test['decrease'].sum()),
        'gm_rr': float(np.exp(test['log_rr'].mean())),
        'd': float(test['log_rr'].mean() / test['log_rr'].std()),
        'p_sign': float(binomtest(int(test['decrease'].sum()), len(test), 0.5, alternative='greater').pvalue),
    },
    'halves_consistent': float(p_2s),
}

with open('data/results/r37_split_sample.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to data/results/r37_split_sample.json")
