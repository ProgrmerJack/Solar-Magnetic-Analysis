#!/usr/bin/env python3
"""
damip_class_law_multimodel.py
=============================
THE THRESHOLD-CLASS LAW UNDER HUMAN FORCING IN SEVERAL MODELS, AND WHERE FORCING
CHANGES THE SPREAD (so that the elasticity is identified).

Registered 2026-10-04, committed BEFORE any data for it were fetched
(03_data_ingestion/acquire_cmip6_monthly_boxes.py). Extends damip_class_law.py
(one model, spread barely changed) on the two weaknesses a reviewer named: one
model, and weakly determined elasticities.

MODELS  CMIP6 historical vs DAMIP hist-nat, monthly 2 m temperature (Amon tas),
  historical members restricted to the physics and forcing variant of the hist-nat
  members (CanESM5 historical also has a p2 physics, excluded):
  PRIMARY    CanESM5 (Pangeo), MIROC6 (Pangeo), HadGEM3-GC31-LL (ESGF CEDA, Amon)
  SECONDARY  GISS-E2-1-G (Pangeo; 20 hist-nat members)
AMENDED 2026-10-04 before any analysis (one test member fetched; no statistic
  computed): the members on Pangeo mix variants. CanESM5 hist-nat has 25 p1 + 25 p2
  and historical 25 p1 + 40 p2 physics members, so the arms would mix two model
  physics in different proportions -- a genuine mixture, the very thing the test
  detects. CanESM5 is therefore analysed as two units, p1 (25 vs 25) and p2 (40 vs
  25), each within one physics, and counts as ONE model in M1 (mean of its two
  units). GISS-E2-1-G mixes forcing variants f2/f3/f4 with at most 11 members per
  arm in any one, and is not analysed.
YEARS 1985-2014; seasons from monthly means (all three months required; DJF by the
  January year).
TARGETS (cos-latitude box means; classes = the hist-nat arm's tercile of that year)
  CEUR JJA 45-55N 5-15E (hot)     WRUS JJA 50-60N 30-60E (hot)
  CNAM JJA 35-50N 260-280E (hot)  CEAS DJF 45-60N 60-120E (cold)
  NSIB DJF 60-75N 60-140E (cold)  NWNA DJF 55-70N 200-260E (cold)
  The two high-latitude winter boxes are added because human forcing is expected to
  reduce winter temperature variance there (research note novelty_positioning.md),
  which is what identifies the elasticity.
TESTS (unit = a year; per model-target as damip_class_law: law residual, shift-only
  residual, s.d. ratio, elasticity, shape)
  M1 PRIMARY (second domain beyond one model): law residual for CEUR JJA averaged
     over the three primary models with equal weight; interval from resampling
     years within each model. READING: inside +-0.25 -> "the law holds in three
     models"; excluding 0 with |mean| > 0.25 -> "fails"; else inconclusive.
  M2 PRIMARY (identification): elasticity of (C_hist - C_nat) on C_nat (k - 1),
     pooled over all primary model-target-years; interval from resampling years
     within model-target. Also over the two high-latitude winter targets alone.
     READING: interval entirely above 0.60 -> "forced mixing of two symmetric
     types >= 2 s.d. apart is excluded in the second domain"; whether 1 is inside
     is reported; the pooled s.d. ratio is reported to show identification.
  M3 per model-target readings, the law residual pooled over the six targets per
     model, and GISS-E2-1-G as a secondary.
Output: results/current/8_experiment/damip_class_law_multimodel.json
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
DATA = HERE.parents[0] / "03_data_ingestion" / "cmip6_monthly_boxes.parquet"
sys.path.insert(0, str(HERE))
import damip_class_law as D1                         # noqa: E402  year_stats, boot, shape

NAME = "damip_class_law_multimodel"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT, TOL = 10000, 0.25
YEARS = range(1985, 2015)
PRIMARY_MODELS = ["CanESM5", "MIROC6", "HadGEM3-GC31-LL"]
SECONDARY_MODELS = ["GISS-E2-1-G"]
TARGETS = {"CEUR_JJA": ("CEUR", "JJA", "hot"), "WRUS_JJA": ("WRUS", "JJA", "hot"), "CNAM_JJA": ("CNAM", "JJA", "hot"),
           "CEAS_DJF": ("CEAS", "DJF", "cold"), "NSIB_DJF": ("NSIB", "DJF", "cold"), "NWNA_DJF": ("NWNA", "DJF", "cold")}
HIGHLAT = ["NSIB_DJF", "NWNA_DJF"]


def seasonal(d, box, season):
    d = d[["arm", "member", "year", "month", box]].copy()
    if season == "JJA":
        d = d[d.month.isin([6, 7, 8])]; d["syear"] = d.year
    else:
        d = d[d.month.isin([12, 1, 2])]; d["syear"] = np.where(d.month == 12, d.year + 1, d.year)
    g = d.groupby(["arm", "member", "syear"])[box].agg(["mean", "count"]).reset_index()
    return g[g["count"] == 3]


def unit_rows(d, target):
    box, season, side = TARGETS[target]
    s = seasonal(d, box, season)
    q = 2 / 3 if side == "hot" else 1 / 3
    rows, zl = [], []
    for y in YEARS:
        xh = s[(s.arm == "historical") & (s.syear == y)]["mean"].values
        xn = s[(s.arm == "hist-nat") & (s.syear == y)]["mean"].values
        if len(xh) < 15 or len(xn) < 15:
            continue
        st = D1.year_stats(xh, xn, q, side); zl.append(st.pop("_z")); st["year"] = y; rows.append(st)
    return pd.DataFrame(rows), zl


def pooled_mean(frames, col, rng):
    """Equal weight per frame; resample years within each frame."""
    vals = [f[col].dropna().values for f in frames]
    pt = float(np.mean([v.mean() for v in vals]))
    bs = [np.mean([rng.choice(v, len(v)).mean() for v in vals]) for _ in range(N_BOOT)]
    return {"mean": round(pt, 4), "ci95": [round(float(x), 4) for x in np.quantile(bs, [0.025, 0.975])],
            "n_frames": len(vals)}


def pooled_by_model(by_model, col, rng):
    """Equal weight per MODEL; a model with two units (CanESM5 p1, p2) is the mean of
    its units; years resampled within each unit."""
    by_model = {m: [f[col].dropna().values for f in fs] for m, fs in by_model.items() if fs}
    pt = float(np.mean([np.mean([v.mean() for v in vs]) for vs in by_model.values()]))
    bs = [np.mean([np.mean([rng.choice(v, len(v)).mean() for v in vs]) for vs in by_model.values()])
          for _ in range(N_BOOT)]
    return {"mean": round(pt, 4), "ci95": [round(float(x), 4) for x in np.quantile(bs, [0.025, 0.975])],
            "n_models": len(by_model)}


def pooled_slope(frames, rng):
    fr = [f.dropna(subset=["diff", "law_pred_diff"]) for f in frames]
    allx = np.concatenate([f.law_pred_diff.values for f in fr]); ally = np.concatenate([f["diff"].values for f in fr])
    b = float(np.polyfit(allx, ally, 1)[0]); bs = []
    for _ in range(N_BOOT):
        xs, ys = [], []
        for f in fr:
            i = rng.integers(len(f), size=len(f)); xs.append(f.law_pred_diff.values[i]); ys.append(f["diff"].values[i])
        bs.append(np.polyfit(np.concatenate(xs), np.concatenate(ys), 1)[0])
    ci = [round(float(x), 4) for x in np.quantile(bs, [0.025, 0.975])]
    return {"slope": round(b, 4), "ci95": ci, "contains_1": bool(ci[0] <= 1 <= ci[1]),
            "reading": ("forced mixing of two symmetric types >= 2 s.d. apart is excluded in the second domain"
                        if ci[0] > 0.60 else "consistent with forced mixing" if ci[1] < 0.60 else "inconclusive"),
            "pooled_sd_ratio_mean": round(float(np.mean(np.concatenate([f.sd_ratio.values for f in fr]))), 4),
            "pooled_sd_ratio_range": [round(float(np.min(np.concatenate([f.sd_ratio.values for f in fr]))), 4),
                                      round(float(np.max(np.concatenate([f.sd_ratio.values for f in fr]))), 4)]}


def calibrate(frames, d, units, n_rep=200):
    """POST HOC (added 2026-10-04 after the registered run): what the same estimators
    give when the law holds EXACTLY. For every model-target-year, one Gaussian
    population with the observed member counts, mean shift and s.d. ratio (in units
    of the natural-forcing arm); n_rep replicates. Reports the expected law residual
    per model-target (and the observed minus expected), and the expected pooled
    elasticity. Own random stream; the registered numbers are unchanged."""
    rc = np.random.default_rng(SEED + 500)
    out = {"note": "post hoc calibration under the exact law (Gaussian, observed n, shift, s.d. ratio)", "per_model_target": {}}
    sims_x, sims_y = [[] for _ in range(n_rep)], [[] for _ in range(n_rep)]
    for (m, t), f in frames.items():
        if m not in units:
            continue
        box, season, side = TARGETS[t]
        s = seasonal(d[d.model == m], box, season)
        q = 2 / 3 if side == "hot" else 1 / 3
        reps = [[] for _ in range(n_rep)]   # per replicate, the years with a defined contrast (as observed)
        for _, row in f.iterrows():
            y = int(row.year)
            nh = int(((s.arm == "historical") & (s.syear == y)).sum()); nn = int(((s.arm == "hist-nat") & (s.syear == y)).sum())
            for r in range(n_rep):
                st = D1.year_stats(rc.normal(row["shift"], row.sd_ratio, nh), rc.normal(0, 1, nn), q, side)
                if np.isfinite(st["resid_law"]):
                    reps[r].append(st["resid_law"])
                if np.isfinite(st["diff"]) and np.isfinite(st["law_pred_diff"]):
                    sims_x[r].append(st["law_pred_diff"]); sims_y[r].append(st["diff"])
        # until 2026-10-05 an undefined year added 0 and the sum was divided by all years,
        # shrinking the range wherever the forced class is often near-empty (FAILURES.md)
        reps = np.array([np.mean(v) for v in reps if v])
        obs = float(f.resid_law.mean())
        out["per_model_target"][f"{m}|{t}"] = {
            "observed": round(obs, 4), "expected_under_law": round(float(reps.mean()), 4),
            "range95_under_law": [round(float(x), 4) for x in np.quantile(reps, [0.025, 0.975])],
            "observed_minus_expected": round(obs - float(reps.mean()), 4),
            "inside_range": bool(np.quantile(reps, 0.025) <= obs <= np.quantile(reps, 0.975))}
    slopes = [np.polyfit(x, y, 1)[0] for x, y in zip(sims_x, sims_y) if len(x) > 5]
    out["expected_elasticity_under_law"] = {"mean": round(float(np.mean(slopes)), 4),
                                            "range95": [round(float(x), 4) for x in np.quantile(slopes, [0.025, 0.975])]}
    inside = [v["inside_range"] for v in out["per_model_target"].values()]
    out["share_targets_inside_range"] = round(float(np.mean(inside)), 4)
    print("post hoc calibration:", out["expected_elasticity_under_law"], "inside", out["share_targets_inside_range"], flush=True)
    return out


def calibrate_empirical(frames, d, units, n_rep=200):
    """POST HOC (added 2026-10-04): as calibrate(), but keeping each year's actual
    outcome shape. Synthetic natural-forcing members are resampled from that year's
    natural-forcing members; synthetic historical members are resampled from them
    and mapped by the year's observed location and scale change, so the law holds
    exactly while the (possibly skewed) shape is the model's own."""
    rc = np.random.default_rng(SEED + 600)
    out = {"note": "post hoc calibration under the exact law with each year's empirical natural-forcing shape",
           "per_model_target": {}}
    for (m, t), f in frames.items():
        if m not in units:
            continue
        box, season, side = TARGETS[t]
        s = seasonal(d[d.model == m], box, season)
        q = 2 / 3 if side == "hot" else 1 / 3
        reps = [[] for _ in range(n_rep)]
        for _, row in f.iterrows():
            y = int(row.year)
            xn = s[(s.arm == "hist-nat") & (s.syear == y)]["mean"].values
            nh = int(((s.arm == "historical") & (s.syear == y)).sum())
            z = (xn - xn.mean()) / xn.std(ddof=1)
            for r in range(n_rep):
                syn_n = rc.choice(z, len(z)); syn_h = row["shift"] + row.sd_ratio * rc.choice(z, nh)
                st = D1.year_stats(syn_h, syn_n, q, side)
                if np.isfinite(st["resid_law"]):
                    reps[r].append(st["resid_law"])
        reps = np.array([np.mean(v) for v in reps if v])
        obs = float(f.resid_law.mean())
        out["per_model_target"][f"{m}|{t}"] = {
            "observed": round(obs, 4), "expected_under_law": round(float(reps.mean()), 4),
            "range95_under_law": [round(float(x), 4) for x in np.quantile(reps, [0.025, 0.975])],
            "inside_range": bool(np.quantile(reps, 0.025) <= obs <= np.quantile(reps, 0.975))}
    out["share_targets_inside_range"] = round(float(np.mean([v["inside_range"] for v in out["per_model_target"].values()])), 4)
    print("post hoc empirical-shape calibration: inside", out["share_targets_inside_range"], flush=True)
    return out


def reading(r):
    lo, hi = r["ci95"]
    if -TOL <= lo and hi <= TOL:
        return "the law holds"
    if (lo > 0 or hi < 0) and abs(r["mean"]) > TOL:
        return "fails"
    return "inconclusive"


def main():
    if not DATA.exists():
        sys.exit("no data: run 03_data_ingestion/acquire_cmip6_monthly_boxes.py")
    d = pd.read_parquet(DATA)
    rng = np.random.default_rng(SEED)
    res = {"registered": "2026-10-04, committed before any data were fetched", "seed": SEED, "n_boot": N_BOOT,
           "tolerance_sigma": TOL, "years": [YEARS[0], YEARS[-1]], "members": {}, "per_model_target": {}}
    frames = {}
    d = d.copy()
    phys = d.member.str.extract(r"i\d+(p\d+)f\d+")[0]
    d.loc[d.model == "CanESM5", "model"] = "CanESM5_" + phys[d.model == "CanESM5"]
    units = ["CanESM5_p1", "CanESM5_p2", "MIROC6", "HadGEM3-GC31-LL"]
    model_of = {"CanESM5_p1": "CanESM5", "CanESM5_p2": "CanESM5", "MIROC6": "MIROC6", "HadGEM3-GC31-LL": "HadGEM3-GC31-LL"}
    for m in units:
        dm = d[d.model == m]
        if dm.empty:
            continue
        res["members"][m] = {a: int(dm[dm.arm == a].member.nunique()) for a in ("historical", "hist-nat")}
        for t in TARGETS:
            f, zl = unit_rows(dm, t)
            if f.empty:
                continue
            frames[(m, t)] = f
            r = {"n_years": int(len(f)), "law": D1.boot(f.resid_law, rng), "shift_only": D1.boot(f.resid_shift, rng),
                 "sd_ratio": D1.boot(f.sd_ratio, rng), "rate_hist": round(float(f.rate_hist.mean()), 4),
                 "rate_nat": round(float(f.rate_nat.mean()), 4)}
            r["reading"] = reading(r["law"])
            r["shape"] = {k: v for k, v in D1.shape(zl, rng).items() if k != "quantile_diff"}
            res["per_model_target"][f"{m}|{t}"] = r
            print(m, t, r["n_years"], "rate", r["rate_nat"], "->", r["rate_hist"], "law", r["law"]["mean"], r["law"]["ci95"],
                  r["reading"], "| sd ratio", r["sd_ratio"]["mean"], "| shape p", r["shape"]["AD_perm_p"], flush=True)
    res["M1_CEUR_JJA_three_models"] = pooled_by_model(
        {mm: [frames[(u, "CEUR_JJA")] for u in units if model_of[u] == mm and (u, "CEUR_JJA") in frames]
         for mm in PRIMARY_MODELS}, "resid_law", rng)
    res["M1_CEUR_JJA_three_models"]["reading"] = reading(res["M1_CEUR_JJA_three_models"])
    allp = [f for (m, t), f in frames.items() if m in units]
    res["M2_elasticity_all"] = pooled_slope(allp, rng)
    res["M2_elasticity_highlat_winter"] = pooled_slope([f for (m, t), f in frames.items()
                                                        if m in units and t in HIGHLAT], rng)
    res["M3_law_all_primary"] = pooled_mean(allp, "resid_law", rng)
    res["M3_law_all_primary"]["reading"] = reading(res["M3_law_all_primary"])
    res["M3_per_model"] = {}
    for m in units:
        fm = [f for (mm, t), f in frames.items() if mm == m]
        if fm:
            res["M3_per_model"][m] = pooled_mean(fm, "resid_law", rng)
            res["M3_per_model"][m]["reading"] = reading(res["M3_per_model"][m])
    res["post_hoc_calibration"] = calibrate(frames, d, units)
    res["post_hoc_calibration_empirical_shape"] = calibrate_empirical(frames, d, units)
    for k in ("M1_CEUR_JJA_three_models", "M2_elasticity_all", "M2_elasticity_highlat_winter", "M3_law_all_primary", "M3_per_model"):
        print(k, res[k], flush=True)
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
