#!/usr/bin/env python3
"""
ao_null_perbin.py  --  VOID. DO NOT CITE. Superseded by boot_calibration.py.
=============================================================================
CONFOUNDED BY SAMPLE SIZE. This script compares a bootstrap SE computed on all
13,864 winter days (real onsets, full data) against a pseudo-onset null SD
computed on 8,862 CLEAN days (real-event zone deleted). Those are different
sample sizes: sqrt(13864/8862) = 1.251, so the null SD is expected to be ~25%
larger from data volume ALONE. The observed median SE ratio was 0.80 -- i.e.
exactly the 25% predicted. The apparent conclusion, "the winter-block bootstrap
understates uncertainty by 20%", does not follow from this comparison at all.

The randomisation p-values it printed are similarly not comparable to the
bootstrap p-values beside them, for the same reason.

Kept only so the mistake is on the record. The valid like-for-like test --
bootstrap SE and null SD from identical clean data under a known truth -- is
boot_calibration.py.

ORIGINAL HEADER FOLLOWS
-----------------------
aao_null.py found the winter-block bootstrap SE to be 27% SMALLER than the clean
pseudo-onset null SD for the AAO at +45..+60. That made a negative control look
like a failure (bootstrap p=0.015) when the randomisation null says it passes
(p=0.081, observed inside the null 95% range).

If the same understatement applies to the AO, every bootstrap p-value in this
project is anti-conservative and the headline intervals are too narrow.

Gate 3 measured AO coverage at 0.947, but for one window under one spec. This
compares bootstrap SE against clean-null SD bin by bin, which is the quantity
that actually failed for the AAO.

  ratio = bootstrap_SE / null_SD.  1.0 is calibrated; <1 means the bootstrap
  understates uncertainty and the reported p-values are too small.

Same clean-null construction as gate3_clean_null.py: real-event influence
(-60..+75 d) deleted from the data first, then pseudo-onsets drawn from the
surviving days at the observed day-of-year set.

Output: ao_null_perbin.json
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
import calibrate_seasonality as C
from build_catalogue import load_catalogue
import gate3_clean_null as G
import canonical_event_study as K

N_DRAW = 2000
K.N_BOOT = 6000


def main():
    d = K.load_series("AO")
    real = load_catalogue("primary")
    real = real[(real >= d.index.min()) & (real <= d.index.max())]
    prof, _ = K.run(d, real, seed=20260730)

    mask = G.real_influence_mask(d.index, real)
    clean = d[~mask]
    print(f"AO: {len(d)} winter days, {len(real)} onsets; "
          f"clean days {len(clean)} ({100*len(clean)/len(d):.0f}%)")

    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])
    rng = np.random.default_rng(20260730)
    draws = {l: [] for l in K.LABELS}
    for i in range(N_DRAW):
        fake = G.draw_clean(clean.index, doys, rng)
        if len(fake) < 10:
            continue
        b = K.fit(K.build(clean, fake))
        if b is None or not np.all(np.isfinite(b)):
            continue
        for j, l in enumerate(K.LABELS):
            draws[l].append(b[j])
        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{N_DRAW} draws", flush=True)

    res = {"n_draw": len(draws[K.LABELS[0]]), "n_boot": K.N_BOOT, "bins": {}}
    print(f"\n{'bin':10s} {'effect':>8s} {'bootSE':>8s} {'nullSD':>8s} "
          f"{'ratio':>7s} {'nullmean':>9s} {'boot p':>8s} {'rand p':>8s}")
    print("-" * 72)
    for l in K.LABELS:
        v = np.array(draws[l])
        e = prof[l]
        lo, hi = e["CI95"]
        se = (hi - lo) / (2 * 1.96)
        ratio = se / v.std()
        rp = float(2 * min((v <= e["effect"]).mean(), (v >= e["effect"]).mean()))
        res["bins"][l] = {
            "effect": e["effect"], "bootstrap_SE": round(float(se), 4),
            "null_SD": round(float(v.std()), 4), "null_mean": round(float(v.mean()), 4),
            "SE_ratio": round(float(ratio), 3),
            "bootstrap_p": e["p_two_sided"], "randomisation_p": rp}
        flag = "  <-- understated" if ratio < 0.85 else ""
        print(f"{l:10s} {e['effect']:+8.3f} {se:8.3f} {v.std():8.3f} {ratio:7.2f} "
              f"{v.mean():+9.3f} {e['p_two_sided']:8.4f} {rp:8.4f}{flag}")

    r = [res["bins"][l]["SE_ratio"] for l in K.LABELS]
    print(f"\n  SE ratio across bins: min {min(r):.2f}, median {np.median(r):.2f}, "
          f"max {max(r):.2f}   (AAO at +45..+60 was 0.73)")
    print("  ratio near 1 -> AO bootstrap is calibrated and the headline p-values hold")
    (RESULTS / "ao_null_perbin.json").write_text(json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> ao_null_perbin.json")


if __name__ == "__main__":
    main()
