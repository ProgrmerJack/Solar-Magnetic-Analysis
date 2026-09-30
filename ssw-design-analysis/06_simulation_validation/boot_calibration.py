#!/usr/bin/env python3
"""
boot_calibration.py
===================
CORRECTS ao_null_perbin.py, which was confounded and must not be used.

That script compared the winter-block bootstrap SE (computed on all 13,864 days,
real onsets) against the pseudo-onset null SD (computed on 8,862 CLEAN days,
after the real-event zone was deleted). Different sample sizes:
sqrt(13864/8862) = 1.251, so a ~25% larger null SD is expected from volume
alone. The observed median ratio was 0.80 -- i.e. exactly 25%. It measured the
sample-size difference, not calibration, and its "bootstrap understates
uncertainty by 20%" reading was an artefact.

The valid test keeps everything on the SAME data. For each pseudo-onset draw on
clean data, fit the estimate AND its winter-block bootstrap interval. Then:

  null SD        = SD of the point estimates across draws  (truth: sampling
                   variability of the estimator under no effect)
  mean boot SE   = average of (CI width)/(2*1.96) across those same draws
  ratio          = mean boot SE / null SD.  ~1.0 is calibrated.
  coverage       = fraction of intervals containing 0, the true null value.

Because both quantities now come from identical data under a known truth, the
ratio is interpretable. This is the per-bin version of the Gate 3 coverage
check (which returned 0.947 for one window under one spec).

Output: boot_calibration.json
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "3_calibration"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue
import gate3_clean_null as G
import canonical_event_study as K
import multi_index_event_study as M

N_DRAW = 200
N_BOOT = 1000        # Gate 3 standing floor


def boot_se(dd, labels, fitfn, buildfn, rng):
    """Winter-block bootstrap SE for every bin, on the data as given."""
    winters = np.array(sorted(dd["winter"].unique()))
    idx = {w: np.flatnonzero((dd["winter"] == w).values) for w in winters}
    out = []
    for _ in range(N_BOOT):
        pick = rng.choice(winters, len(winters), replace=True)
        sub = dd.iloc[np.concatenate([idx[w] for w in pick])].copy()
        sub["winter"] = np.concatenate([np.full(len(idx[w]), j) for j, w in enumerate(pick)])
        b = fitfn(sub)
        if b is not None and np.all(np.isfinite(b)):
            out.append(b)
    a = np.array(out)
    lo = np.percentile(a, 2.5, axis=0)
    hi = np.percentile(a, 97.5, axis=0)
    return (hi - lo) / (2 * 1.96), lo, hi


def run(tag, d, labels, fitfn, buildfn):
    real = load_catalogue("primary")
    real = real[(real >= d.index.min()) & (real <= d.index.max())]
    clean = d[~G.real_influence_mask(d.index, real)]
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    rng = np.random.default_rng(20260730)
    ests, ses, covs = [], [], []
    for i in range(N_DRAW):
        fake = G.draw_clean(clean.index, doys, rng)
        if len(fake) < 10:
            continue
        dd = buildfn(clean, fake)
        e = fitfn(dd)
        if e is None or not np.all(np.isfinite(e)):
            continue
        se, lo, hi = boot_se(dd, labels, fitfn, buildfn, rng)
        ests.append(e); ses.append(se); covs.append((lo <= 0) & (hi >= 0))
        if (i + 1) % 50 == 0:
            print(f"  {tag} {i+1}/{N_DRAW}", flush=True)
    E = np.array(ests); S = np.array(ses); C_ = np.array(covs)
    res = {}
    print(f"\n=== {tag}: {len(E)} draws x {N_BOOT} bootstrap reps (same clean data) ===")
    print(f"{'bin':10s} {'nullSD':>8s} {'meanBootSE':>11s} {'ratio':>7s} {'coverage':>9s}")
    print("-" * 50)
    for j, l in enumerate(labels):
        sd = float(E[:, j].std()); mse = float(S[:, j].mean())
        cov = float(C_[:, j].mean())
        res[l] = {"null_SD": round(sd, 4), "mean_bootstrap_SE": round(mse, 4),
                  "SE_ratio": round(mse / sd, 3), "coverage_of_zero": round(cov, 4)}
        flag = "  <-- anti-conservative" if cov < 0.90 else ""
        print(f"{l:10s} {sd:8.4f} {mse:11.4f} {mse/sd:7.2f} {cov:9.3f}{flag}")
    r = [res[l]["SE_ratio"] for l in labels]
    c = [res[l]["coverage_of_zero"] for l in labels]
    print(f"  SE ratio median {np.median(r):.2f}; coverage median {np.median(c):.3f} "
          f"(nominal 0.95)")
    return res


def main():
    out = {"n_draw": N_DRAW, "n_boot": N_BOOT, "series": {}}
    out["series"]["AO"] = run("AO", K.load_series("AO"), K.LABELS, K.fit, K.build)
    out["series"]["AAO"] = run("AAO", M.load("aao"), M.LABELS, M.fit, M.build)
    (RESULTS / "boot_calibration.json").write_text(json.dumps(out, indent=2), encoding="utf8")
    print("\nSaved -> boot_calibration.json")


if __name__ == "__main__":
    main()
