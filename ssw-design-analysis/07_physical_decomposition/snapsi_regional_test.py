#!/usr/bin/env python3
"""
snapsi_regional_test.py
=======================
REGIONAL COLD RISK UNDER AN IMPOSED SSW: IS IT THE CIRCULATION SHIFT, OR MORE?

Plan approved 2026-09-26 ("I approve the new plan"). Everything below -- regions,
window, metrics, tests, falsifiers -- was written before any regional
temperature was analysed; the tas download (acquire_snapsi_tas.py) stores a
generic 2.5-degree daily field, so it did not steer these choices.

QUESTION
  After an imposed SSW, (a) how much does regional near-surface temperature and
  the risk of a cold fortnight change, and (b) is that change what the SSW's
  shift of the polar-cap circulation implies through the ordinary, SSW-free
  relation between that circulation and regional temperature -- or does the SSW
  (or the "downward" label) carry regional information beyond it? No study
  located in the literature review of 2026-09-26 tests (b) (research_notes/.../
  regional_cold_impacts.md, Q2).

DATA
  SNAPSI nudged and control, the four NH initialisations, every centre passing
  snapsi_selection_test.corruption_guard (NRL excluded). Member tas: daily means
  on a 2.5-degree grid (_snapsi_tas/). Member circulation: the polar-cap
  sea-level-pressure NAM proxy of result L, -(psl_cap - control mean)/control s.d.
  over post-onset days +8..+25.

REGIONS (cos-latitude-weighted means over all grid points, land and sea)
  NEURASIA   50-65N, 10-130E   Kretschmer et al. 2018 (npj Clim Atmos Sci)  PRIMARY
  HI_EUROPE  55-70N,  0-60E    Huang et al. 2021 (Commun Earth Environ)
  MID_EASIA  35-55N, 90-150E   Huang et al. 2021
  MID_NAMER  35-55N, 120-60W   Huang et al. 2021
  WINDOW: daily means of post-onset days +8..+24 (00 UTC day 8 to 00 UTC day 25,
  the span of result L's window).
  T = window-mean regional tas, standardised by the control ensemble of the same
  centre and initialisation (mean 0, s.d. 1): sigma units as for the NAM proxy.

TESTS (unit: centre x initialisation pair; intervals from 10,000 bootstrap
resamples of CENTRES, as in result L)
  T1 regional shift     S_T = mean(T nudged) - mean(T control); and the cold
                        risk P(T nudged < control 10th percentile) against 0.10.
  T2 MEDIATION (primary) fit T = a + b N on the CONTROL members of each pair;
                        predict each nudged member from its own N; residual
                        R = mean over nudged members of (T - a - b N).
                        H0: R = 0 (the regional change is the circulation shift
                        acting through the SSW-free relation).
                        FALSIFIER: the 95% interval of R excludes 0 AND |R| > 0.2
                        sigma in NEURASIA -> the SSW has a regional effect beyond
                        the polar-cap circulation, and the paper must say so.
                        Reported also: the share of S_T the shift explains,
                        1 - R / S_T.
  T3 LABEL              within each nudged ensemble, T regressed on N and the
                        downward indicator (Karpechko conditions 1-2 as in L),
                        both demeaned within ensemble; pooled coefficient of the
                        indicator. H0: 0 (the label adds nothing beyond N).
                        Same regression on control members for comparison.
  T4 COLD RISK          the probability of T < control 10th percentile that the
                        shift predicts: Gaussian with mean a + b N_i and the
                        control residual s.d., averaged over nudged members i;
                        against the observed nudged frequency.
  Secondary: all four regions; T2 per initialisation.
  Sensitivity added 2026-09-26 AFTER the primary result: every test without ECCC,
  the Tukey outlier of the causal shift (result K), as K reports it.

VERIFICATION
  --selftest: synthetic members with T = b N + noise in both arms (R must be ~0,
  interval covering 0 in >= 90% of 200 draws) and with an added +0.5 sigma direct
  effect in nudged (must be detected in >= 90%).

Output: results/current/8_experiment/snapsi_regional_test.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
TAS = HERE.parents[0] / "03_data_ingestion" / "_snapsi_tas"
sys.path.insert(0, str(HERE))
import snapsi_selection_test as L                   # noqa: E402

NAME = "snapsi_regional_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT = 10000
REGIONS = {"NEURASIA": (50, 65, 10, 130), "HI_EUROPE": (55, 70, 0, 60),
           "MID_EASIA": (35, 55, 90, 150), "MID_NAMER": (35, 55, 240, 300)}
PRIMARY = "NEURASIA"
TAS_DAYS = (8, 24)
COLD_Q = 0.10


def region_series(centre, exp, init):
    """member -> window-mean regional tas (K) for every region."""
    off = L.post_onset_offset(init)
    lo, hi = TAS_DAYS[0] + off, TAS_DAYS[1] + off
    out = {}
    for f in sorted(TAS.glob(f"{centre}_{exp}_{init}_*.npz")):
        member = f.stem.split("_", 3)[3]
        with np.load(f) as z:
            day, tas, n = z["lead_day"], z["tas"], z["n_steps"]
            lat, lon = z["lat"], z["lon"]
        sel = (day >= lo) & (day <= hi)
        if sel.sum() != hi - lo + 1 or (n[sel] != 4).any():
            continue                                     # window not fully covered
        field = tas[sel].mean(axis=0)
        row = {}
        for r, (la0, la1, lo0, lo1) in REGIONS.items():
            mi = (lat >= la0) & (lat <= la1)
            mj = (lon >= lo0) & (lon <= lo1)
            w = np.cos(np.deg2rad(lat[mi]))[:, None] * np.ones(mj.sum())[None, :]
            row[r] = float((field[np.ix_(mi, mj)] * w).sum() / w.sum())
        out[member] = row
    return pd.DataFrame.from_dict(out, orient="index")


def nam_and_label(centre, init):
    """Per arm: member -> (N, DW) with the control standardisation of result L."""
    con = L.load(centre, "control", init)
    if con is None or not L.spans_window(con, init):
        return None
    cm = L.window_means(con, init)
    base, sd = float(cm.mean()), float(cm.std(ddof=1))
    off = L.post_onset_offset(init)
    lo, hi = L.WINDOW[0] + off, L.WINDOW[1] + off
    res = {}
    for arm in ("nudged", "control"):
        d = L.load(centre, arm, init)
        if d is None or not L.spans_window(d, init):
            return None
        w = d[(d["lead_days"] >= lo) & (d["lead_days"] <= hi)].copy()
        w["nam"] = -(w["psl_cap"] - base) / sd
        g = w.groupby("member")["nam"]
        mean, frac = g.mean(), g.apply(lambda s: float((s < 0).mean()))
        res[arm] = pd.DataFrame({"N": mean, "DW": ((mean < 0) & (frac > 0.5)).astype(float)})
    return res


def pair_stats(Tn, Tc, Nn, Nc, DWn, DWc):
    """All per-pair quantities for one region (T already standardised)."""
    b, a = np.polyfit(Nc, Tc, 1)
    resid_c = Tc - (a + b * Nc)
    s_res = float(np.std(resid_c, ddof=2))
    q10 = float(np.quantile(Tc, COLD_Q))
    from scipy.stats import norm
    pred_cold = float(np.mean(norm.cdf((q10 - (a + b * Nn)) / s_res)))
    out = {"S_T": float(Tn.mean() - Tc.mean()), "R": float(np.mean(Tn - (a + b * Nn))),
           "slope_control": float(b), "cold_obs": float(np.mean(Tn < q10)),
           "cold_pred": pred_cold, "cold_control": float(np.mean(Tc < q10))}
    for arm, T, N, D in (("nudged", Tn, Nn, DWn), ("control", Tc, Nc, DWc)):
        if 3 <= D.sum() <= len(D) - 3:
            X = np.column_stack([np.ones(len(T)), N - N.mean(), D - D.mean()])
            coef = np.linalg.lstsq(X, T - T.mean(), rcond=None)[0]
            out[f"label_coef_{arm}"] = float(coef[2])
    return out


def centre_boot(df, col, rng):
    """Mean over pairs and 95% interval from resampling centres."""
    df = df.dropna(subset=[col])
    cs = df["centre"].unique()
    by = {c: df.loc[df.centre == c, col].values for c in cs}
    bs = np.empty(N_BOOT)
    for k in range(N_BOOT):
        pick = rng.choice(cs, len(cs))
        bs[k] = np.mean(np.concatenate([by[c] for c in pick]))
    return {"mean": round(float(df[col].mean()), 4),
            "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
            "n_pairs": int(len(df)), "n_centres": int(len(cs))}


def build_pairs():
    centres = sorted({p.name.split("_")[0] for p in TAS.glob("*.npz")})
    rows = []
    for c in centres:
        bad, _ = L.corruption_guard(c)
        if bad:
            continue
        for init in L.ONSET:
            nl = nam_and_label(c, init)
            if nl is None:
                continue
            tn, tc = region_series(c, "nudged", init), region_series(c, "control", init)
            n_, c_ = nl["nudged"].join(tn, how="inner"), nl["control"].join(tc, how="inner")
            if len(n_) < 10 or len(c_) < 10:
                continue
            for r in REGIONS:
                mu, sd = c_[r].mean(), c_[r].std(ddof=1)
                st = pair_stats(((n_[r] - mu) / sd).values, ((c_[r] - mu) / sd).values,
                                n_.N.values, c_.N.values, n_.DW.values, c_.DW.values)
                rows.append({"centre": c, "init": init, "region": r,
                             "n_nudged": len(n_), "n_control": len(c_), **st})
    return pd.DataFrame(rows)


def summarise(pairs, rng):
    out = {}
    for r in REGIONS:
        d = pairs[pairs.region == r]
        s = {k: centre_boot(d, k, rng) for k in ("S_T", "R", "slope_control", "cold_obs",
                                                 "cold_pred", "cold_control")}
        for arm in ("nudged", "control"):
            col = f"label_coef_{arm}"
            if col in d:
                s[col] = centre_boot(d, col, rng)
        s["share_of_shift_explained"] = round(1 - s["R"]["mean"] / s["S_T"]["mean"], 3) \
            if abs(s["S_T"]["mean"]) > 1e-9 else None
        s["R_by_init"] = {i: round(float(g.R.mean()), 3) for i, g in d.groupby("init")}
        out[r] = s
    return out


def selftest():
    rng = np.random.default_rng(1)
    hits0 = hits1 = 0
    for k in range(200):
        rows = []
        for c in range(8):
            for i in range(4):
                Nc = rng.normal(0, 1, 50); Nn = rng.normal(-1.1, 1, 50)
                b = rng.uniform(0.3, 0.8)
                Tc = b * Nc + rng.normal(0, 0.7, 50)
                Tn0 = b * Nn + rng.normal(0, 0.7, 50)
                for eff, lab in ((0.0, "null"), (-0.5, "eff")):
                    Tn = Tn0 + eff
                    mu, sd = Tc.mean(), Tc.std(ddof=1)
                    st = pair_stats((Tn - mu) / sd, (Tc - mu) / sd, Nn, Nc,
                                    (Nn < 0).astype(float), (Nc < 0).astype(float))
                    rows.append({"centre": f"c{c}", "case": lab, **st})
        df = pd.DataFrame(rows)
        sub = np.random.default_rng(k)
        global N_BOOT
        nb, N_BOOT = N_BOOT, 500
        z = centre_boot(df[df.case == "null"], "R", sub)
        e = centre_boot(df[df.case == "eff"], "R", sub)
        N_BOOT = nb
        hits0 += z["ci95"][0] <= 0 <= z["ci95"][1]
        hits1 += e["ci95"][1] < 0
    print(f"selftest: null interval covers 0 in {hits0}/200; -0.5 sigma effect detected in {hits1}/200")
    return hits0 >= 180 and hits1 >= 180


def main():
    if "--selftest" in sys.argv:
        return 0 if selftest() else 1
    rng = np.random.default_rng(SEED)
    pairs = build_pairs()
    if pairs.empty:
        sys.exit("no pairs: is _snapsi_tas populated?")
    res = {"plan_approved": "2026-09-26", "seed": SEED, "n_boot": N_BOOT,
           "regions": REGIONS, "primary": PRIMARY, "tas_days": list(TAS_DAYS),
           "nam_window_days": list(L.WINDOW), "cold_quantile": COLD_Q,
           "inputs": ["ssw-design-analysis/03_data_ingestion/_snapsi_tas",
                      "ssw-design-analysis/03_data_ingestion/_snapsi_reduced"],
           "n_pairs": int(pairs[pairs.region == PRIMARY].shape[0]),
           "centres": sorted(pairs.centre.unique()),
           "results": summarise(pairs, rng),
           "sensitivity_excl_ECCC": summarise(pairs[pairs.centre != "ECCC"], rng),
           "pairs": pairs.round(4).to_dict(orient="records")}
    p = res["results"][PRIMARY]
    res["falsifier_primary"] = {
        "rule": "R interval excludes 0 AND |R| > 0.2 sigma in NEURASIA",
        "met": bool((p["R"]["ci95"][0] > 0 or p["R"]["ci95"][1] < 0) and abs(p["R"]["mean"]) > 0.2)}
    for r, s in res["results"].items():
        print(f"{r:10s} S_T {s['S_T']['mean']:+.3f} {s['S_T']['ci95']}  R {s['R']['mean']:+.3f} "
              f"{s['R']['ci95']}  explained {s['share_of_shift_explained']}  cold obs "
              f"{s['cold_obs']['mean']:.3f} pred {s['cold_pred']['mean']:.3f} "
              f"(control {s['cold_control']['mean']:.3f})  label(nudged) "
              f"{s.get('label_coef_nudged', {}).get('mean')} {s.get('label_coef_nudged', {}).get('ci95')}",
              flush=True)
    print("falsifier met:", res["falsifier_primary"]["met"])
    (RESULTS / "snapsi_regional_test.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
