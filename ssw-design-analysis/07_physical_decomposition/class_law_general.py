#!/usr/bin/env python3
"""
class_law_general.py
====================
IS THE THRESHOLD-CLASS LAW GENERAL? PRECIPITATION AS WELL AS TEMPERATURE, FIVE
MODELS, AND FORCING STRONG ENOUGH TO DISTINGUISH ONE POPULATION FROM FORCED MIXING
OF RESPONSE TYPES.

Registered 2026-10-05, committed BEFORE any data for it were fetched
(03_data_ingestion/acquire_cmip6_general.py writes cmip6_general_boxes.parquet).
Already seen, so NOT new evidence for T1: monthly tas, historical vs hist-nat
1985-2014, in CanESM5 and MIROC6 (damip_class_law_multimodel.py). Everything else
(precipitation; EC-Earth3 and ACCESS-ESM1-5; ssp585; T3 everywhere) is unseen.

UNITS  model unit x box x season x variable. Model units: CanESM5_p1, CanESM5_p2
  (physics variants kept apart, 6346c9e), MIROC6, EC-Earth3, ACCESS-ESM1-5, and
  HadGEM3-GC31-LL (precipitation, T1 and T2 only). Boxes CEUR WRUS CNAM CEAS NSIB
  NWNA (acquire_cmip6_monthly_boxes.py); seasons JJA and DJF (January year), all
  three months required; a year needs >= 15 members in each arm.
FORCING LEVELS (30 years each; years paired by index within the period)
  L  CanESM5, MIROC6, HadGEM3: hist-nat 1985-2014; EC-Earth3, ACCESS: historical 1951-1980
  M  historical 1985-2014
  H  ssp585 2071-2100 (not HadGEM3)
CLASSES  beyond the tercile of the reference arm's members of that year:
  tas JJA hot (upper), tas DJF cold (lower), pr JJA dry (lower), pr DJF wet (upper).

T1 PRIMARY -- the location-scale law, L -> M. Per unit, per year: the reference
   members mapped x -> m_f + k (x - m_r) and classified at the same threshold;
   residual = forced contrast minus mapped contrast (damip_class_law.year_stats,
   corrected 2026-10-05); intervals from 10,000 resamples of years; calibration
   under the exact law with each year's own reference shape (resampled reference
   members, mapped by the year's location and scale change; 200 replicates; mean
   over the years with a defined contrast, as observed) -> 95% range per unit.
   READING per variable, fixed now: the pooled residual (equal weight per model;
   years resampled within unit) read against +-0.25 sigma as in damip_class_law;
   "the law holds for <variable> classes" if that reading is "holds" AND >= 85% of
   units lie inside their calibration range; "fails" if the pooled reading is
   "fails" OR < 70% inside; otherwise inconclusive.
T2 SECONDARY -- a location-scale-skew law, for units outside their T1 range. The
   standardised reference z is mapped z -> z + (dg / 6)(z^2 - 1) (first-order
   Cornish-Fisher change of skewness), restandardised, then located and scaled as
   in T1; dg = forced minus reference skewness, each the mean over years of the
   year's sample skewness (one dg per unit). Calibration as T1 under this law,
   re-estimating dg in every replicate. Reported: how many T1 departures lie inside
   their T2 range, and the share of map values that are non-monotone.
T3 PRIMARY -- one population or forced mixing of response types. For fixed types
   whose weights change with forcing (any number of types along a linear path in
   the weights; any two types), the ensemble variance is a quadratic in the
   ensemble mean with curvature exactly -1: V = E[X^2] - mu^2 and both E[X^2] and
   mu are linear in the weights. One population whose spread changes linearly with
   its mean gives curvature b^2 >= 0. Per unit: points (mu, V) of every arm-year of
   L, M, H in units of L (mean and pooled within-year s.d. of L); OLS fit
   V = c0 + c1 mu + c2 mu^2; kappa = c2; 95% interval from 2,000 resamples of year
   indices (all three arms of an index together).
   READING per unit: "mixing excluded" if the lower bound > -0.5; "consistent with
   mixing" if the upper bound < -0.5; otherwise "not identified". Per variable:
   "forced mixing of response types is excluded" if >= 10 units are identified and
   >= 90% of identified units read "mixing excluded"; "mixing present" if >= 50%
   read "consistent with mixing"; otherwise inconclusive. Also per model.
T4 SECONDARY -- the T1 law under strong forcing, M -> H, for units with a defined
   contrast in >= 15 years (shifts of several s.d. leave the forced class empty
   elsewhere; expected mainly for precipitation). Reading as T1, per unit.

Output: results/current/8_experiment/class_law_general.json
Run:    class_law_general.py [--selftest]
"""
import json
import sys
import warnings
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import skew

warnings.filterwarnings("ignore", category=RuntimeWarning)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
DATA = HERE.parents[0] / "03_data_ingestion" / "cmip6_general_boxes.parquet"
sys.path.insert(0, str(HERE))
import damip_class_law as D1                         # noqa: E402  year_stats (corrected), contrast, boot

NAME = "class_law_general"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT, N_BOOT_K, N_REP, TOL, MIN_MEMBERS = 10000, 2000, 200, 0.25, 15
BOXES = ["CEUR", "WRUS", "CNAM", "CEAS", "NSIB", "NWNA"]
SIDE = {("tas", "JJA"): "hot", ("tas", "DJF"): "cold", ("pr", "JJA"): "cold", ("pr", "DJF"): "hot"}
LEVEL_L = {"CanESM5_p1": ("hist-nat", 1985), "CanESM5_p2": ("hist-nat", 1985), "MIROC6": ("hist-nat", 1985),
           "HadGEM3-GC31-LL": ("hist-nat", 1985), "EC-Earth3": ("historical", 1951), "ACCESS-ESM1-5": ("historical", 1951)}
LEVEL_M, LEVEL_H = ("historical", 1985), ("ssp585", 2071)
MODEL_OF = {"CanESM5_p1": "CanESM5", "CanESM5_p2": "CanESM5"}


# ------------------------------------------------------------------ data
def seasonal(d, var, box, season):
    d = d[(d["var"] == var)][["experiment", "member", "year", "month", box]]
    if season == "JJA":
        d = d[d.month.isin([6, 7, 8])].assign(syear=lambda x: x.year)
    else:
        d = d[d.month.isin([12, 1, 2])].assign(syear=lambda x: np.where(x.month == 12, x.year + 1, x.year))
    g = d.groupby(["experiment", "member", "syear"])[box].agg(["mean", "count"]).reset_index()
    return g[g["count"] == 3]


def arm_years(s, level):
    """{year index: member values} for one forcing level (experiment, first year)."""
    exp, y0 = level
    a = s[(s.experiment == exp) & (s.syear >= y0) & (s.syear < y0 + 30)]
    return {int(y - y0): g["mean"].values for y, g in a.groupby("syear") if len(g) >= MIN_MEMBERS}


# ------------------------------------------------------------------ T1 / T4
def pair_frame(ref, frc, q, side):
    rows = []
    for i in sorted(set(ref) & set(frc)):
        st = D1.year_stats(frc[i], ref[i], q, side); st.pop("_z"); st["i"] = i; rows.append(st)
    return pd.DataFrame(rows)


def calib_ls(ref, frc, f, q, side, rng):
    """95% range of the mean law residual when the law holds exactly with each year's own
    reference shape; mean over the years with a defined contrast (as observed)."""
    reps = [[] for _ in range(N_REP)]
    for _, row in f.iterrows():
        xr = ref[int(row.i)]; z = (xr - xr.mean()) / xr.std(ddof=1); nf = len(frc[int(row.i)])
        for r in range(N_REP):
            st = D1.year_stats(row["shift"] + row.sd_ratio * rng.choice(z, nf), rng.choice(z, len(z)), q, side)
            if np.isfinite(st["resid_law"]):
                reps[r].append(st["resid_law"])
    reps = np.array([np.mean(v) for v in reps if v])
    if len(reps) < N_REP // 2:
        return None, None
    return [round(float(x), 4) for x in np.quantile(reps, [0.025, 0.975])], round(float(reps.mean()), 4)


# ------------------------------------------------------------------ T2
def cf_map(z, dg):
    t = z + dg / 6.0 * (z * z - 1.0)
    return t


def skew_resid(xf, xr, q, side, dg):
    mu, sd = xr.mean(), xr.std(ddof=1)
    f, r = (xf - mu) / sd, (xr - mu) / sd
    thr = float(np.quantile(r, q))
    z = (r - r.mean()) / r.std(ddof=1)
    t = cf_map(z, dg); t = (t - t.mean()) / t.std(ddof=1)
    return D1.contrast(f, thr, side) - D1.contrast(f.mean() + f.std(ddof=1) * t, thr, side)


def pooled_dg(ref, frc):
    ks = sorted(set(ref) & set(frc))
    return float(np.mean([skew(frc[i], bias=False) for i in ks]) - np.mean([skew(ref[i], bias=False) for i in ks]))


def t2_unit(ref, frc, q, side, rng):
    dg = pooled_dg(ref, frc)
    ks = sorted(set(ref) & set(frc))
    res = np.array([skew_resid(frc[i], ref[i], q, side, dg) for i in ks])
    obs = float(np.nanmean(res))
    nonmono = float(np.mean([np.mean(1 + dg / 3.0 * ((ref[i] - ref[i].mean()) / ref[i].std(ddof=1)) <= 0) for i in ks]))
    reps = []
    for _ in range(N_REP):
        syn_r, syn_f = {}, {}
        for i in ks:
            xr = ref[i]; z = (xr - xr.mean()) / xr.std(ddof=1)
            sh = (frc[i].mean() - xr.mean()) / xr.std(ddof=1); k = frc[i].std(ddof=1) / xr.std(ddof=1)
            tz = cf_map(z, dg); m_t, s_t = tz.mean(), tz.std(ddof=1)
            syn_r[i] = rng.choice(z, len(z))
            syn_f[i] = sh + k * (cf_map(rng.choice(z, len(frc[i])), dg) - m_t) / s_t
        dg_s = pooled_dg(syn_r, syn_f)
        v = [skew_resid(syn_f[i], syn_r[i], q, side, dg_s) for i in ks]
        v = [x for x in v if np.isfinite(x)]
        if v:
            reps.append(np.mean(v))
    lo, hi = np.quantile(reps, [0.025, 0.975])
    return {"dg": round(dg, 4), "observed": round(obs, 4), "range95_under_law": [round(float(lo), 4), round(float(hi), 4)],
            "inside_range": bool(lo <= obs <= hi), "share_nonmonotone": round(nonmono, 5)}


# ------------------------------------------------------------------ T3
def curvature_points(arms):
    """arms: list of {index: values}; returns per index the (mu, V) points of every arm,
    in units of the first arm (its mean, and its pooled within-year s.d.)."""
    L = arms[0]
    m0 = np.mean([v.mean() for v in L.values()]); s0 = np.sqrt(np.mean([v.var(ddof=1) for v in L.values()]))
    pts = {}
    for a in arms:
        for i, v in a.items():
            pts.setdefault(i, []).append(((v.mean() - m0) / s0, v.var(ddof=1) / s0 ** 2))
    return pts


def kappa(pts, idx):
    P = np.array([p for i in idx for p in pts[i]])
    return float(np.polyfit(P[:, 0], P[:, 1], 2)[0])


def t3_unit(arms, rng):
    pts = curvature_points(arms)
    idx = np.array(sorted(pts))
    k0 = kappa(pts, idx)
    bs = [kappa(pts, rng.choice(idx, len(idx))) for _ in range(N_BOOT_K)]
    lo, hi = np.quantile(bs, [0.025, 0.975])
    P = np.array([p for i in idx for p in pts[i]])
    read = "mixing excluded" if lo > -0.5 else "consistent with mixing" if hi < -0.5 else "not identified"
    return {"kappa": round(k0, 4), "ci95": [round(float(lo), 4), round(float(hi), 4)], "reading": read,
            "mean_range_sd": round(float(P[:, 0].max() - P[:, 0].min()), 3),
            "var_range": [round(float(P[:, 1].min()), 3), round(float(P[:, 1].max()), 3)], "n_index": int(len(idx))}


# ------------------------------------------------------------------ readings
def reading(r):
    lo, hi = r["ci95"]
    if -TOL <= lo and hi <= TOL:
        return "holds"
    if (lo > 0 or hi < 0) and abs(r["mean"]) > TOL:
        return "fails"
    return "inconclusive"


def pooled_by_model(by_model, rng):
    by_model = {m: vs for m, vs in by_model.items() if vs}
    pt = float(np.mean([np.mean([v.mean() for v in vs]) for vs in by_model.values()]))
    bs = [np.mean([np.mean([rng.choice(v, len(v)).mean() for v in vs]) for vs in by_model.values()]) for _ in range(N_BOOT)]
    return {"mean": round(pt, 4), "ci95": [round(float(x), 4) for x in np.quantile(bs, [0.025, 0.975])], "n_models": len(by_model)}


def t1_reading(pooled_read, share_inside):
    if pooled_read == "holds" and share_inside >= 0.85:
        return "the law holds"
    if pooled_read == "fails" or share_inside < 0.70:
        return "fails"
    return "inconclusive"


def t3_reading(reads):
    ident = [r for r in reads if r != "not identified"]
    if len(ident) >= 10 and np.mean([r == "mixing excluded" for r in ident]) >= 0.90:
        return "forced mixing of response types is excluded"
    if ident and np.mean([r == "consistent with mixing" for r in ident]) >= 0.50:
        return "mixing present"
    return "inconclusive"


# ------------------------------------------------------------------ main
def analyse(d):
    rng = np.random.default_rng(SEED)
    rng_c = np.random.default_rng(SEED + 1); rng_2 = np.random.default_rng(SEED + 2); rng_3 = np.random.default_rng(SEED + 3)
    d = d.copy()
    can = d.model == "CanESM5"
    d.loc[can, "model"] = "CanESM5_" + d.loc[can, "member"].str.extract(r"i\d+(p\d+)f\d+")[0]
    res = {"registered": "2026-10-05, committed before any data were fetched", "seed": SEED, "n_boot": N_BOOT,
           "n_boot_kappa": N_BOOT_K, "n_rep_calibration": N_REP, "tolerance_sigma": TOL, "min_members": MIN_MEMBERS,
           "members": {}, "units": {}}
    t1_by = {"tas": {}, "pr": {}}; t1_in = {"tas": [], "pr": []}; t3_reads = {"tas": [], "pr": []}; t3_by_model = {}
    for m in LEVEL_L:
        dm = d[d.model == m]
        if dm.empty:
            continue
        res["members"][m] = {f"{e}|{v}": int(g.member.nunique()) for (e, v), g in dm.groupby(["experiment", "var"])}
        for var in ("tas", "pr"):
            for season in ("JJA", "DJF"):
                side = SIDE[(var, season)]; q = 2 / 3 if side == "hot" else 1 / 3
                for box in BOXES:
                    s = seasonal(dm, var, box, season)
                    L, M, H = arm_years(s, LEVEL_L[m]), arm_years(s, LEVEL_M), arm_years(s, LEVEL_H)
                    if len(L) < 15 or len(M) < 15:
                        continue
                    key = f"{m}|{var}|{box}_{season}"; u = {"class": side, "n_years_LM": len(set(L) & set(M))}
                    f = pair_frame(L, M, q, side)
                    v = f.resid_law.dropna()
                    if len(v) >= 5:
                        u["T1"] = {"law": D1.boot(v, rng), "shift": D1.boot(f["shift"], rng), "sd_ratio": D1.boot(f.sd_ratio, rng),
                                   "rate_ref": round(float(f.rate_nat.mean()), 4), "rate_forced": round(float(f.rate_hist.mean()), 4)}
                        u["T1"]["reading"] = reading(u["T1"]["law"])
                        rg, ex = calib_ls(L, M, f.dropna(subset=["resid_law"]), q, side, rng_c)
                        u["T1"]["range95_under_law"] = rg; u["T1"]["expected_under_law"] = ex
                        u["T1"]["inside_range"] = bool(rg is not None and rg[0] <= u["T1"]["law"]["mean"] <= rg[1])
                        t1_by[var].setdefault(MODEL_OF.get(m, m), []).append(v.values); t1_in[var].append(u["T1"]["inside_range"])
                        if not u["T1"]["inside_range"]:
                            u["T2"] = t2_unit(L, M, q, side, rng_2)
                    if len(H) >= 15:
                        u["T3"] = t3_unit([L, M, H], rng_3)
                        t3_reads[var].append(u["T3"]["reading"])
                        t3_by_model.setdefault(f"{MODEL_OF.get(m, m)}|{var}", []).append(u["T3"]["reading"])
                        f4 = pair_frame(M, H, q, side); v4 = f4.resid_law.dropna()
                        if len(v4) >= 15:
                            u["T4"] = {"law": D1.boot(v4, rng), "sd_ratio": D1.boot(f4.sd_ratio, rng),
                                       "rate_ref": round(float(f4.rate_nat.mean()), 4), "rate_forced": round(float(f4.rate_hist.mean()), 4)}
                            u["T4"]["reading"] = reading(u["T4"]["law"])
                    res["units"][key] = u
                    print(key, {k: (u[k].get("reading"), u[k].get("law", {}).get("mean") if k != "T3" else u[k]["kappa"])
                                for k in ("T1", "T3", "T4") if k in u}, flush=True)
    res["T1_by_variable"] = {}
    for var in ("tas", "pr"):
        if not t1_in[var]:
            continue
        p = pooled_by_model(t1_by[var], rng); p["reading_pooled"] = reading(p)
        share = float(np.mean(t1_in[var]))
        t2 = [u["T2"]["inside_range"] for k, u in res["units"].items() if f"|{var}|" in k and "T2" in u]
        res["T1_by_variable"][var] = {"pooled": p, "n_units": len(t1_in[var]), "share_inside_calibration": round(share, 4),
                                      "reading": t1_reading(p["reading_pooled"], share),
                                      "T2_departures_inside_skew_law": f"{sum(t2)}/{len(t2)}"}
    res["T3_by_variable"] = {var: {"n_units": len(r), "n_identified": sum(x != "not identified" for x in r),
                                   "n_mixing_excluded": sum(x == "mixing excluded" for x in r),
                                   "n_consistent_with_mixing": sum(x == "consistent with mixing" for x in r),
                                   "reading": t3_reading(r)} for var, r in t3_reads.items() if r}
    res["T3_by_model"] = {k: {"n_identified": sum(x != "not identified" for x in r),
                              "n_mixing_excluded": sum(x == "mixing excluded" for x in r), "reading": t3_reading(r)}
                          for k, r in t3_by_model.items()}
    t4 = [u["T4"]["reading"] for u in res["units"].values() if "T4" in u]
    res["T4_summary"] = {"n_units": len(t4), **{r: t4.count(r) for r in ("holds", "inconclusive", "fails")}}
    for k in ("T1_by_variable", "T3_by_variable", "T3_by_model", "T4_summary"):
        print(k, json.dumps(res[k]), flush=True)
    return res


def selftest():
    """Synthetic checks of T1, T2 and T3 (no model data)."""
    rng = np.random.default_rng(1)
    out = {}
    def mix(n, w, D=3.0):
        return np.where(rng.random(n) < w, D + rng.normal(size=n), rng.normal(size=n))
    arms = [{i: mix(40, w) for i in range(30)} for w in (0.2, 0.5, 0.8)]
    out["T3_two_types_D3"] = t3_unit(arms, rng)
    arms = [{i: mu + s * rng.normal(size=40) for i in range(30)} for mu, s in ((0, 1), (1, 1), (4, 0.7))]
    out["T3_one_population"] = t3_unit(arms, rng)
    arms = [{i: (mu + s * rng.gamma(4, 1, 40)) for i in range(30)} for mu, s in ((0, 1), (0.5, 1.2), (2, 1.6))]
    out["T3_one_population_gamma"] = t3_unit(arms, rng)
    ref = {i: rng.normal(size=40) for i in range(30)}; frc = {i: 0.8 + 1.1 * rng.normal(size=40) for i in range(30)}
    f = pair_frame(ref, frc, 1 / 3, "cold")
    rg, ex = calib_ls(ref, frc, f.dropna(subset=["resid_law"]), 1 / 3, "cold", rng)
    out["T1_exact_law"] = {"observed": round(float(f.resid_law.mean()), 4), "range95": rg, "expected": ex}
    ref = {i: rng.gamma(4, 1, 40) for i in range(30)}; frc = {i: rng.gamma(16, 0.5, 40) - 3.5 for i in range(30)}   # skewness 1.0 -> 0.5
    f = pair_frame(ref, frc, 1 / 3, "cold")
    rg, ex = calib_ls(ref, frc, f.dropna(subset=["resid_law"]), 1 / 3, "cold", rng)
    out["T1_skew_change"] = {"observed": round(float(f.resid_law.mean()), 4), "range95": rg}
    out["T2_skew_change"] = t2_unit(ref, frc, 1 / 3, "cold", rng)
    print(json.dumps(out, indent=1))
    ok = (out["T3_two_types_D3"]["reading"] == "consistent with mixing" and out["T3_one_population"]["reading"] == "mixing excluded"
          and out["T3_one_population_gamma"]["reading"] == "mixing excluded"
          and out["T1_exact_law"]["range95"][0] <= out["T1_exact_law"]["observed"] <= out["T1_exact_law"]["range95"][1])
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    if "--selftest" in sys.argv:
        return selftest()
    if not DATA.exists():
        sys.exit("no data: run 03_data_ingestion/acquire_cmip6_general.py")
    res = analyse(pd.read_parquet(DATA))
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
