#!/usr/bin/env python3
"""
damip_class_law.py
==================
IS THE THRESHOLD-CLASS LAW A PROPERTY OF OUTCOME-DEFINED CLASSES IN GEOSCIENCE,
NOT OF THE SSW LABEL? A SECOND DOMAIN: ANTHROPOGENIC FORCING AND EXTREME SEASONS.

Registered 2026-10-03, committed BEFORE any member of these ensembles was fetched
(03_data_ingestion/acquire_damip_hadgem.py writes the reduced data).

INTERVENTION  CMIP6 DAMIP, HadGEM3-GC31-LL: `historical` (all forcings, 55 members,
  r*i1p1f3) against `hist-nat` (natural forcings only, 60 members, r*i1p1f3): the
  same coupled model; the arms differ only in human forcing (research note
  cross_domain_interventions.md). Years 1985-2014 (30), 360-day calendar.
OUTCOMES  area-weighted (cos-latitude) box means of daily 2 m temperature, seasonal
  means per member and year:
  PRIMARY    central Europe, 45-55N, 5-15E, June-August (heat; the box of the
             heatwave-class literature, e.g. Stefanon et al. 2012)
  SECONDARY  western Russia 50-60N, 30-60E, JJA; central North America 35-50N,
             100-80W, JJA; central Eurasia 45-60N, 60-120E, December-February
             (cold; approximately the region of Mori et al. 2014)
CLASS  "extreme season" = beyond the hist-nat arm's tercile of that year (upper
  tercile for JJA heat, lower tercile for DJF cold): a fixed threshold, as the
  downward label's zero. Secondary: the 90th (10th) percentile.
UNIT  a year: its 55 historical members against its 60 hist-nat members, values
  standardised by the hist-nat members of that year (the SNAPSI "pair").
TESTS (intervals from 10,000 resamples of years)
  D1 rate: share in the extreme class, historical vs hist-nat; shift; s.d. ratio.
  D2 PRIMARY, the law: residual = C_hist - C_LS, where C_LS is the class contrast
     of the hist-nat members mapped x -> m_h + k (x - m_n) (k = ratio of s.d.),
     classified with the same threshold. READING, fixed now: interval inside
     +-0.25 sigma -> "the law holds in a second domain"; interval excluding 0
     with |mean| > 0.25 -> "it fails"; otherwise inconclusive.
  D3 shift-only residual (no rescaling), and the share of years with |residual|
     < 0.25 for the law and for shift only.
  D4 shape: within each year both arms standardised by their own mean and s.d.;
     pooled residual quantiles (q 0.05..0.95) and the Anderson-Darling k-sample
     statistic with 2,000 within-year label permutations; |dq| <= 0.25 and p >= 0.05
     -> "shape unchanged by the forcing".
  D5 elasticity: slope of (C_hist - C_nat) on the law's prediction C_nat (k - 1)
     over years, as F3 of snapsi_free_natural_experiment (mixture < 0.60).
CAVEATS stated now: one model; seasonal means, not daily events; the forcing
  difference grows over 1985-2014, so years differ in shift (that is the dose).

Output: results/current/8_experiment/damip_class_law.json
"""
import json
import sys
import warnings
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", message=".*(variant|p-value (capped|floored)).*")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
DATA = HERE.parents[0] / "03_data_ingestion" / "damip_hadgem_boxes.parquet"
NAME = "damip_class_law"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT, N_PERM, TOL, MIN_CLASS = 10000, 2000, 0.25, 3
YEARS = range(1985, 2015)
TARGETS = {"CEUR_JJA": ("CEUR", "JJA", "hot"), "WRUS_JJA": ("WRUS", "JJA", "hot"),
           "CNAM_JJA": ("CNAM", "JJA", "hot"), "CEAS_DJF": ("CEAS", "DJF", "cold")}
PRIMARY = "CEUR_JJA"
QS = np.round(np.arange(0.05, 0.951, 0.05), 2)


def seasonal(d, box, season):
    """member x year seasonal means for one arm (360-day calendar: months 30 days)."""
    d = d[["arm", "member", "year", "month", box]].copy()
    if season == "JJA":
        d = d[d.month.isin([6, 7, 8])]
        d["syear"] = d.year
    else:                                   # DJF labelled by the January year
        d = d[d.month.isin([12, 1, 2])]
        d["syear"] = np.where(d.month == 12, d.year + 1, d.year)
    g = d.groupby(["arm", "member", "syear"])[box].agg(["mean", "count"]).reset_index()
    return g[g["count"] >= 85]               # 90 days in a 360-day season; allow none missing beyond 5


def contrast(x, thr, side):
    cls = x > thr if side == "hot" else x < thr
    if cls.sum() < MIN_CLASS or (~cls).sum() < MIN_CLASS:
        return np.nan
    return float(x[cls].mean() - x[~cls].mean())


def year_stats(xh, xn, q, side):
    mu, sd = xn.mean(), xn.std(ddof=1)
    h, n = (xh - mu) / sd, (xn - mu) / sd
    thr = float(np.quantile(n, q))
    m_h, m_n, k = h.mean(), n.mean(), h.std(ddof=1) / n.std(ddof=1)
    ch, cn = contrast(h, thr, side), contrast(n, thr, side)
    rate = (lambda x: float((x > thr).mean() if side == "hot" else (x < thr).mean()))
    return {"rate_hist": rate(h), "rate_nat": rate(n), "shift": float(m_h - m_n), "sd_ratio": float(k),
            "C_hist": ch, "C_nat": cn, "resid_law": ch - contrast(m_n + k * (n - m_n), thr, side),
            "resid_shift": ch - contrast(n + (m_h - m_n), thr, side), "diff": ch - cn,
            "law_pred_diff": cn * (k - 1), "_z": ((h - m_h) / h.std(ddof=1), (n - m_n) / n.std(ddof=1))}


def boot(v, rng):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    bs = [rng.choice(v, len(v)).mean() for _ in range(N_BOOT)]
    return {"mean": round(float(v.mean()), 4), "ci95": [round(float(x), 4) for x in np.quantile(bs, [0.025, 0.975])],
            "n_years": int(len(v))}


def reading(r):
    lo, hi = r["ci95"]
    if -TOL <= lo and hi <= TOL:
        return "the law holds in a second domain"
    if (lo > 0 or hi < 0) and abs(r["mean"]) > TOL:
        return "it fails"
    return "inconclusive"


def shape(zl, rng):
    from scipy.stats import anderson_ksamp
    zh = np.concatenate([a for a, _ in zl]); zn = np.concatenate([b for _, b in zl])
    dq = np.quantile(zh, QS) - np.quantile(zn, QS)
    ad0 = anderson_ksamp([zh, zn]).statistic; ex = 0
    for _ in range(N_PERM):
        a, b = [], []
        for x, y in zl:
            z = np.concatenate([x, y]); p = rng.permutation(len(z)); a.append(z[p[:len(x)]]); b.append(z[p[len(x):]])
        ex += anderson_ksamp([np.concatenate(a), np.concatenate(b)]).statistic >= ad0
    p = (ex + 1) / (N_PERM + 1)
    return {"max_abs_quantile_diff": round(float(np.max(np.abs(dq))), 4), "AD_perm_p": round(p, 4),
            "quantile_diff": dict(zip([str(x) for x in QS], np.round(dq, 4).tolist())),
            "reading": "shape unchanged by the forcing" if np.max(np.abs(dq)) <= TOL and p >= 0.05 else "shape changed"}


def analyse(d, target, q, rng):
    box, season, side = TARGETS[target]
    s = seasonal(d, box, season)
    rows, zl = [], []
    for y in YEARS:
        xh = s[(s.arm == "historical") & (s.syear == y)]["mean"].values
        xn = s[(s.arm == "hist-nat") & (s.syear == y)]["mean"].values
        if len(xh) < 20 or len(xn) < 20:
            continue
        st = year_stats(xh, xn, q, side); zl.append(st.pop("_z")); st["year"] = y; rows.append(st)
    df = pd.DataFrame(rows)
    out = {"n_years": int(len(df)), "members": [int(s[s.arm == "historical"].member.nunique()),
                                                 int(s[s.arm == "hist-nat"].member.nunique())]}
    for col in ("rate_hist", "rate_nat", "shift", "sd_ratio", "resid_law", "resid_shift", "diff"):
        out[col] = boot(df[col], rng)
    out["reading_law"] = reading(out["resid_law"])
    ok = df.dropna(subset=["resid_law", "resid_shift"])
    out["share_abs_lt_tol"] = {"law": round(float((ok.resid_law.abs() < TOL).mean()), 4),
                               "shift_only": round(float((ok.resid_shift.abs() < TOL).mean()), 4)}
    d5 = df.dropna(subset=["diff", "law_pred_diff"])
    if len(d5) > 5:
        b = np.polyfit(d5.law_pred_diff, d5["diff"], 1)[0]
        bs = []
        for _ in range(N_BOOT):
            i = rng.integers(len(d5), size=len(d5))
            if d5.law_pred_diff.values[i].std() > 0:
                bs.append(np.polyfit(d5.law_pred_diff.values[i], d5["diff"].values[i], 1)[0])
        out["elasticity_slope"] = {"slope": round(float(b), 4),
                                   "ci95": [round(float(x), 4) for x in np.quantile(bs, [0.025, 0.975])]}
    out["shape"] = shape(zl, rng)
    out["years"] = df.round(4).to_dict(orient="records")
    return out


def main():
    if not DATA.exists():
        sys.exit("no data: run 03_data_ingestion/acquire_damip_hadgem.py")
    d = pd.read_parquet(DATA)
    res = {"registered": "2026-10-03, committed before any member was fetched", "seed": SEED, "n_boot": N_BOOT,
           "tolerance_sigma": TOL, "model": "HadGEM3-GC31-LL", "years": [YEARS[0], YEARS[-1]], "targets": {}}
    for i, t in enumerate(TARGETS):
        for q_lab, q in (("tercile", None), ("decile", None)):
            side = TARGETS[t][2]
            qq = (2 / 3 if side == "hot" else 1 / 3) if q_lab == "tercile" else (0.9 if side == "hot" else 0.1)
            r = analyse(d, t, qq, np.random.default_rng(SEED + 10 * i + (q_lab == "decile")))
            res["targets"][f"{t}_{q_lab}"] = r
            print(t, q_lab, {k: r[k] for k in ("n_years", "rate_hist", "rate_nat", "resid_law", "resid_shift",
                                               "reading_law", "share_abs_lt_tol")},
                  "| slope", r.get("elasticity_slope"), "| shape", {k: r["shape"][k] for k in ("max_abs_quantile_diff", "AD_perm_p", "reading")},
                  flush=True)
    res["primary"] = res["targets"][f"{PRIMARY}_tercile"]["reading_law"]
    # POST HOC (added 2026-10-04 after the registered run): the estimator's own bias
    # when the law holds EXACTLY (one Gaussian population, 55 vs 60 members, shift 0.6),
    # at the two thresholds; separate random stream, registered numbers unchanged.
    rng_b = np.random.default_rng(SEED + 99)
    bias = {}
    for q_lab, qq in (("tercile", 2 / 3), ("decile", 0.9)):
        for k in (1.0, 1.1):
            r, sft = [], []
            for _ in range(3000):
                st = year_stats(rng_b.normal(0.6, k, 55), rng_b.normal(0, 1, 60), qq, "hot")
                r.append(st["resid_law"]); sft.append(st["resid_shift"])
            bias[f"{q_lab}_k{k}"] = {"law_resid_mean": round(float(np.nanmean(r)), 4),
                                     "law_resid_sd": round(float(np.nanstd(r)), 4),
                                     "shift_resid_mean": round(float(np.nanmean(sft)), 4)}
    res["post_hoc_estimator_bias"] = {"note": "post hoc; law holds exactly in these simulations", **bias}
    print("post hoc estimator bias:", bias)
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
