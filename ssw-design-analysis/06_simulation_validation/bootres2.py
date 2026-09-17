"""Decisive test: does coverage recover with more bootstrap replicates?"""
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path('../02_event_catalogues').resolve()))
import calibrate_seasonality as C, gate3_clean_null as G
from build_catalogue import load_catalogue
from scipy import stats
import os
SPEC=os.environ.get('SPEC','harm6')
ao = C.read_cpc(C.ROOT/'data/processed/atmospheric/ao_daily_cpc.txt')
d = pd.DataFrame({'y': ao.values}, index=ao.index)
d = d[np.isin(d.index.month, C.SEASON)]
d['winter'] = C.winter_of(d.index); d['doy'] = d.index.dayofyear
on = load_catalogue('primary'); on = on[(on >= d.index.min()) & (on <= d.index.max())]
clean = d[~G.real_influence_mask(d.index, on)].copy()
doys = pd.DatetimeIndex(on).dayofyear.values
def cp(k,n,a=0.05):
    lo=stats.beta.ppf(a/2,k,n-k+1) if k else 0.0
    hi=stats.beta.ppf(1-a/2,k+1,n-k) if k<n else 1.0
    return round(float(lo),3),round(float(hi),3)
print(f"{'n_boot':>7} {'draws':>6} {'coverage':>9} {'CI95':>16} {'FPR':>6}")
for nb, nd in ((1000, 300),):
    pass
for nb, nd in ((1000, 300),):
    rng=np.random.default_rng(23); cov=rej=0
    for _ in range(nd):
        fk=G.draw_clean(clean.index,doys,rng)
        ci=C.bootstrap_ci(clean,fk,SPEC,n_boot=nb,seed=int(rng.integers(1e6)))
        if ci:
            if ci[0]<=0<=ci[1]: cov+=1
            else: rej+=1
    n=cov+rej
    print(f"{nb:>7} {nd:>6} {cov/n:>9.3f} {str(cp(cov,n)):>16} {rej/n:>6.3f}", flush=True)
