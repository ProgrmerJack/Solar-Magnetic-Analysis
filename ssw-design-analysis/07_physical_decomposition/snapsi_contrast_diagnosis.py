#!/usr/bin/env python3
"""
snapsi_contrast_diagnosis.py
============================
WHY IS THE NUDGED DW-NDW CONTRAST SMALLER THAN ONE SHIFTED GAUSSIAN PREDICTS, AND
IS THE SOUTHERN HEMISPHERE CASE A TRANSLATION TOO?

Plan approved 2026-09-29. In result L the nudged contrast falls 0.19 sigma
[0.11, 0.30] short of the contrast a cut at zero gives on one Gaussian with the
ensemble's own mean and s.d. That expectation (a) is Gaussian and (b) models
only condition 1 of the criterion. Written before these diagnostics were run.

TESTS (NH pairs of result L: same centres, initialisations, window days +8..+25,
control standardisation; centre-cluster bootstrap, 10,000 resamples)
  D1 CONDITION 1 ONLY: the residual from the Gaussian expectation when DW is
     defined by the window mean alone. If the shortfall comes from condition 2,
     D1 is ~0.
  D2 EMPIRICAL SHIFTED NULL: the control members of the same pair, every 6-hourly
     value shifted by the pair's imposed effect (nudged mean minus control mean of
     the window-mean NAM proxy), classified with the full criterion (conditions
     1-2); residual = nudged contrast minus this contrast. It keeps the shape of
     the SSW-free distribution and condition 2. If one shifted population
     suffices, D2 is ~0.
  D3 SKEWNESS of member window means, nudged minus control.
  D4 SOUTHERN HEMISPHERE variance ratio (nudged/control, members re-centred on
     their ensemble mean and pooled, s20190829), 95% interval from 10,000
     resamples of members within ensembles, as for the NH ratio of result M.

Output: results/current/8_experiment/snapsi_contrast_diagnosis.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
from scipy.stats import skew

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
sys.path.insert(0, str(HERE))
import snapsi_selection_test as L                   # noqa: E402

NAME = "snapsi_contrast_diagnosis"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT = 10000
MIN_CLASS = 3


def member_nam(centre, init):
    """{arm: {member: 6-hourly NAM proxy over the window}} with L's standardisation."""
    con = L.load(centre, "control", init)
    if con is None or not L.spans_window(con, init):
        return None
    cm = L.window_means(con, init)
    base, sd = float(cm.mean()), float(cm.std(ddof=1))
    off = L.post_onset_offset(init)
    lo, hi = L.WINDOW[0] + off, L.WINDOW[1] + off
    out = {}
    for arm in ("nudged", "control"):
        d = L.load(centre, arm, init)
        if d is None or not L.spans_window(d, init):
            return None
        w = d[(d["lead_days"] >= lo) & (d["lead_days"] <= hi)]
        out[arm] = {m: (-(g["psl_cap"].values - base) / sd) for m, g in w.groupby("member")}
    return out


def contrast(series, cond2=True):
    means = np.array([s.mean() for s in series])
    dw = means < 0
    if cond2:
        dw &= np.array([(s < 0).mean() > 0.5 for s in series])
    if dw.sum() < MIN_CLASS or (~dw).sum() < MIN_CLASS:
        return None
    return float(means[dw].mean() - means[~dw].mean())


def cboot(vals, cens, rng):
    vals, cens = np.asarray(vals), np.asarray(cens)
    uc = np.unique(cens)
    bs = [np.mean(np.concatenate([vals[cens == c] for c in rng.choice(uc, len(uc))])) for _ in range(N_BOOT)]
    return {"mean": round(float(vals.mean()), 4),
            "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
            "n_pairs": int(len(vals)), "n_centres": int(len(uc))}


def main():
    rng = np.random.default_rng(SEED)
    centres = sorted({p.name.split("_")[0] for p in L.RED.glob("*.parquet")})
    usable = [c for c in centres if not L.corruption_guard(c)[0]]
    rows = []
    for c in usable:
        for init in L.ONSET:
            nm = member_nam(c, init)
            if nm is None:
                continue
            nud, con = list(nm["nudged"].values()), list(nm["control"].values())
            mn = np.array([s.mean() for s in nud]); mc = np.array([s.mean() for s in con])
            cn2, cc2 = contrast(nud), contrast(con)
            if cn2 is None or cc2 is None:
                continue                                   # result L's eligible pairs
            shift = float(mn.mean() - mc.mean())
            cn1 = contrast(nud, cond2=False)
            shifted = contrast([s + shift for s in con])
            rows.append({"centre": c, "init": init,
                         "D1_resid_cond1": (cn1 - L.threshold_contrast(mn.mean(), mn.std(ddof=1)))
                         if cn1 is not None else np.nan,
                         "D2_resid_empirical": (cn2 - shifted) if shifted is not None else np.nan,
                         "resid_gaussian_full": cn2 - L.threshold_contrast(mn.mean(), mn.std(ddof=1)),
                         "D3_skew_diff": float(skew(mn) - skew(mc)),
                         "sd_ratio": float(mn.std(ddof=1) / mc.std(ddof=1))})
    import pandas as pd
    df = pd.DataFrame(rows)
    res = {"plan_approved": "2026-09-29", "seed": SEED, "n_boot": N_BOOT, "window": list(L.WINDOW),
           "n_pairs": int(len(df)), "centres": sorted(df.centre.unique())}
    for k in ("resid_gaussian_full", "D1_resid_cond1", "D2_resid_empirical", "D3_skew_diff", "sd_ratio"):
        ok = df[k].notna()
        res[k] = cboot(df.loc[ok, k].values, df.loc[ok, "centre"].values, rng)
    # D4: Southern Hemisphere variance ratio
    nud_c, con_c = [], []
    ens = []
    for c in usable:
        for init in ("s20190829",):
            n, k = L.load(c, "nudged", init), L.load(c, "control", init)
            if n is None or k is None or not (L.spans_window(n, init) and L.spans_window(k, init)):
                continue
            wn, wk = L.window_means(n, init).values, L.window_means(k, init).values
            sd = wk.std(ddof=1)
            ens.append(((wn - wn.mean()) / sd, (wk - wk.mean()) / sd))
    vr = float(np.var(np.concatenate([e[0] for e in ens]), ddof=1) /
               np.var(np.concatenate([e[1] for e in ens]), ddof=1))
    bs = []
    for _ in range(N_BOOT):
        a = np.concatenate([rng.choice(e[0], len(e[0])) for e in ens])
        b = np.concatenate([rng.choice(e[1], len(e[1])) for e in ens])
        bs.append(np.var(a, ddof=1) / np.var(b, ddof=1))
    res["D4_SH_variance_ratio"] = {"ratio": round(vr, 4),
                                  "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
                                  "n_ensembles": len(ens),
                                  "n_nudged": int(sum(len(e[0]) for e in ens)),
                                  "n_control": int(sum(len(e[1]) for e in ens))}
    res["per_pair"] = df.round(4).to_dict(orient="records")
    for k in ("resid_gaussian_full", "D1_resid_cond1", "D2_resid_empirical", "D3_skew_diff", "sd_ratio"):
        print(f"{k:22s} {res[k]['mean']:+.3f} {res[k]['ci95']} (n={res[k]['n_pairs']})")
    print("D4 SH variance ratio", res["D4_SH_variance_ratio"])
    (RESULTS / "snapsi_contrast_diagnosis.json").write_text(json.dumps(res, indent=2),
                                                            encoding="utf8", newline="\n")
    print("Saved -> snapsi_contrast_diagnosis.json")


if __name__ == "__main__":
    main()
