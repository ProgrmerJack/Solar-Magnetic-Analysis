"""
R60 Mixed-Effects NB2 Analysis - Swiss Natural Dry Slab Counts
Fits GEE NB2 with winter-level exchangeable correlation and covariates.
"""
import pandas as pd
import numpy as np
import json
import warnings
warnings.filterwarnings('ignore')
import statsmodels.api as sm
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import NegativeBinomial, Poisson
from statsmodels.genmod.cov_struct import Exchangeable

df = pd.read_parquet('data/processed/analysis_panel_v2.parquet')
w = df[df['winter_id'].notna()].copy()
d = w[w['dry_natural_size_1234'].notna()].copy()

d['z500_std'] = (d['ncep_z500_nh'] - d['ncep_z500_nh'].mean()) / d['ncep_z500_nh'].std()
d['precip_std'] = (d['snotel_prec_mean'] - d['snotel_prec_mean'].mean()) / d['snotel_prec_mean'].std()

winter_map = {wid: i for i, wid in enumerate(sorted(d['winter_id'].unique()))}
d['winter_int'] = d['winter_id'].map(winter_map)

endog = d['dry_natural_size_1234'].values
groups = d['winter_int'].values

nw = d['winter_id'].nunique()
nssw = int(d['ssw_within_15d'].sum())
print(f"DATA: N={len(d)}, winters={nw}, SSW_days={nssw}")
print(f"Mean={endog.mean():.3f}, Var={endog.var():.3f}, Zeros={(endog==0).sum()}")

results = {}

# Model 1: GEE NB2 + covariates (z500 + precip)
print("\n--- GEE NB2 + covariates ---")
exog1 = sm.add_constant(d[['ssw_within_15d', 'z500_std', 'precip_std']].values)
gee1 = GEE(endog, exog1, groups=groups, family=NegativeBinomial(alpha=1.0), cov_struct=Exchangeable())
r1 = gee1.fit()
irr1 = np.exp(r1.params[1])
ci1 = [np.exp(r1.params[1] - 1.96*r1.bse[1]), np.exp(r1.params[1] + 1.96*r1.bse[1])]
p1 = r1.pvalues[1]
print(f"SSW IRR={irr1:.3f} [{ci1[0]:.3f}, {ci1[1]:.3f}] P={p1:.6f}")
print(f"Z500 P={r1.pvalues[2]:.6f}, Precip P={r1.pvalues[3]:.6f}")
results['gee_nb2_covariates'] = {
    'irr': round(float(irr1), 3),
    'ci': [round(float(ci1[0]), 3), round(float(ci1[1]), 3)],
    'p': float(p1),
    'z500_p': float(r1.pvalues[2]),
    'precip_p': float(r1.pvalues[3])
}

# Model 2: GEE NB2 SSW only (matches manuscript's GEE)
print("\n--- GEE NB2 SSW only ---")
exog2 = sm.add_constant(d[['ssw_within_15d']].values)
gee2 = GEE(endog, exog2, groups=groups, family=NegativeBinomial(alpha=1.0), cov_struct=Exchangeable())
r2 = gee2.fit()
irr2 = np.exp(r2.params[1])
ci2 = [np.exp(r2.params[1] - 1.96*r2.bse[1]), np.exp(r2.params[1] + 1.96*r2.bse[1])]
p2 = r2.pvalues[1]
print(f"SSW IRR={irr2:.3f} [{ci2[0]:.3f}, {ci2[1]:.3f}] P={p2:.6f}")
results['gee_nb2_ssw_only'] = {
    'irr': round(float(irr2), 3),
    'ci': [round(float(ci2[0]), 3), round(float(ci2[1]), 3)],
    'p': float(p2)
}

# Model 3: GLM NB2 + covariates (no clustering)
print("\n--- GLM NB2 + covariates (no clustering) ---")
glm3 = sm.GLM(endog, exog1, family=sm.families.NegativeBinomial(alpha=1.0))
r3 = glm3.fit()
irr3 = np.exp(r3.params[1])
ci3 = [np.exp(r3.params[1] - 1.96*r3.bse[1]), np.exp(r3.params[1] + 1.96*r3.bse[1])]
p3 = r3.pvalues[1]
print(f"SSW IRR={irr3:.3f} [{ci3[0]:.3f}, {ci3[1]:.3f}] P={p3:.6f}")
results['glm_nb2_covariates'] = {
    'irr': round(float(irr3), 3),
    'ci': [round(float(ci3[0]), 3), round(float(ci3[1]), 3)],
    'p': float(p3)
}

# Model 4: GEE Poisson + covariates 
print("\n--- GEE Poisson + covariates ---")
gee4 = GEE(endog, exog1, groups=groups, family=Poisson(), cov_struct=Exchangeable())
r4 = gee4.fit()
irr4 = np.exp(r4.params[1])
ci4 = [np.exp(r4.params[1] - 1.96*r4.bse[1]), np.exp(r4.params[1] + 1.96*r4.bse[1])]
p4 = r4.pvalues[1]
print(f"SSW IRR={irr4:.3f} [{ci4[0]:.3f}, {ci4[1]:.3f}] P={p4:.6f}")
results['gee_poisson_covariates'] = {
    'irr': round(float(irr4), 3),
    'ci': [round(float(ci4[0]), 3), round(float(ci4[1]), 3)],
    'p': float(p4)
}

# Model 5: GEE NB2 + Z500 only (no precip)
print("\n--- GEE NB2 + Z500 only ---")
exog5 = sm.add_constant(d[['ssw_within_15d', 'z500_std']].values)
gee5 = GEE(endog, exog5, groups=groups, family=NegativeBinomial(alpha=1.0), cov_struct=Exchangeable())
r5 = gee5.fit()
irr5 = np.exp(r5.params[1])
ci5 = [np.exp(r5.params[1] - 1.96*r5.bse[1]), np.exp(r5.params[1] + 1.96*r5.bse[1])]
p5 = r5.pvalues[1]
print(f"SSW IRR={irr5:.3f} [{ci5[0]:.3f}, {ci5[1]:.3f}] P={p5:.6f}")
results['gee_nb2_z500_only'] = {
    'irr': round(float(irr5), 3),
    'ci': [round(float(ci5[0]), 3), round(float(ci5[1]), 3)],
    'p': float(p5)
}

results['n_obs'] = len(d)
results['n_winters'] = int(nw)
results['n_ssw_days'] = nssw
results['outcome'] = 'dry_natural_size_1234 (Swiss natural dry slab counts)'

with open('data/results/r60_mixed_effects_swiss.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to data/results/r60_mixed_effects_swiss.json")
