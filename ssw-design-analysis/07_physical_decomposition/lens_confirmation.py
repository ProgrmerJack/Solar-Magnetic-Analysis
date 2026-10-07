#!/usr/bin/env python3
"""
lens_confirmation.py
====================
CONFIRMATION ON UNTOUCHED ENSEMBLES: IS FORCED MIXING OF TWO FIXED RESPONSE TYPES
INCOMPATIBLE, AND DOES THE THRESHOLD-CLASS LAW HOLD, IN LARGE ENSEMBLES THIS PROJECT
HAS NEVER ANALYSED FOR SURFACE OUTCOMES?

Registered 2026-10-07, committed BEFORE any value of these data was fetched
(03_data_ingestion/acquire_lens_confirmation.py writes lens_confirmation_boxes.parquet).
EXPOSURE AUDIT: the project has used CESM2 only as CMIP6 zonal-mean stratospheric wind
(event catalogue) and as the CESM2(WACCM6) January-2021 forecasts; MPI-GE never. No
surface temperature or precipitation of either ensemble has been examined.

ENSEMBLES (units)  CESM2-cmip6 (CESM2-LE, standard biomass-burning forcing, 50 members),
  CESM2-smbb (smoothed biomass burning, 50 members), MPI-GE (100 members). The two CESM2
  variants differ in forcing, so they are separate units and count as ONE model when
  models are weighted. Scenarios: CESM2 SSP3-7.0, MPI-GE RCP8.5 (not equated; each model
  is compared only with itself).
LEVELS (30 years each; seasons as class_law_general: JJA, and DJF by January year)
  L historical 1951-1980;  M 2040-2069 (intermediate);  H 2070-2099
  A level-year needs >= 15 members.
CASES  3 units x {tas, pr} x 6 boxes (CEUR WRUS CNAM CEAS NSIB NWNA) x {JJA, DJF}; the
  family is the 36 cases of a variable. All are reported; no case is selected by power
  or by any observed shift.
C1 PRIMARY -- mixture_moment_test.tests exactly as frozen in 70df9fa (raw-moment
  collinearity T, two-way member x 5-year-block bootstrap, B = 1999; validation gate
  passed, 0e57742). Per case:
    "incompatible with fixed two-type mixing"  Holm-adjusted p_T < 0.05 (within variable)
    "compatible with fixed two-type mixing"    the 95% interval of rho lies inside
                                               [-0.25, +0.25] (rho = 0 under mixing)
    "unresolved"                               otherwise
  READING per variable, fixed now: "replicated: fixed two-type mixing rejected" if at
  least 12 of the 36 cases are incompatible and none is compatible; "fixed two-type
  mixing supported" if at least 12 are compatible and none incompatible; "mixed" if
  both kinds occur and either reaches 12; otherwise "unresolved".
C2 SECONDARY -- the CDF-collinearity test of the same function (p_cdf), same family,
  categories and reading.
C3 SECONDARY -- the threshold-class law on fresh data, M -> H (2040-69 -> 2070-99, where
  even temperature classes stay defined): class_law_general T1 machinery (location-scale
  null mapped to the forced mean; classes beyond the M arm's tercile of that year: tas
  JJA hot, DJF cold; pr JJA dry, DJF wet), cases with a defined residual in >= 15 years;
  calibration of 200 replicates with each year's own reference shape. READING per
  variable as class_law_general T1: "the law holds" if the pooled residual (models equal
  weight) lies inside +-0.25 sigma and >= 85% of cases lie inside their calibrated range;
  "fails" if the pooled reading fails or < 70% inside; otherwise inconclusive.
C4 SECONDARY -- the same law test L -> M, for cases with >= 15 defined years.
Reported descriptively: rho per case, its median per variable, forced change of the
level means in s.d. of L.

Output: results/current/8_experiment/lens_confirmation.json
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
DATA = HERE.parents[0] / "03_data_ingestion" / "lens_confirmation_boxes.parquet"
sys.path.insert(0, str(HERE))
import class_law_general as CL                      # noqa: E402  seasonal, pair_frame, calib_ls, pooled_by_model, reading
import damip_class_law as D1                        # noqa: E402  boot
import mixture_moment_test as MM                    # noqa: E402  tests, holm

NAME = "lens_confirmation"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
UNITS = {"CESM2-cmip6": "CESM2", "CESM2-smbb": "CESM2", "MPI-GE": "MPI-GE"}
LEVELS = (("historical", 1951), ("future", 2040), ("future", 2070))
TOL_RHO, MIN_FAMILY = 0.25, 12


def arm_matrix(s, level):
    exp, y0 = level
    a = s[(s.experiment == exp) & (s.syear >= y0) & (s.syear < y0 + 30)]
    piv = a.pivot_table(index="member", columns="syear", values="mean").reindex(columns=range(y0, y0 + 30))
    piv.loc[:, piv.notna().sum() < CL.MIN_MEMBERS] = np.nan
    return piv.index.values, piv.values


def law_block(ref, frc, q, side, rng, rng_c):
    f = CL.pair_frame(ref, frc, q, side)
    v = f.resid_law.dropna() if len(f) else pd.Series(dtype=float)
    if len(v) < 15:
        return None, None
    out = {"law": D1.boot(v, rng), "sd_ratio": D1.boot(f.sd_ratio, rng),
           "rate_ref": round(float(f.rate_nat.mean()), 4), "rate_forced": round(float(f.rate_hist.mean()), 4)}
    out["reading"] = CL.reading(out["law"])
    rg, ex = CL.calib_ls(ref, frc, f.dropna(subset=["resid_law"]), q, side, rng_c)
    out["range95_under_law"] = rg; out["expected_under_law"] = ex
    out["inside_range"] = bool(rg is not None and rg[0] <= out["law"]["mean"] <= rg[1])
    return out, v.values


def family_reading(cats):
    n_inc, n_com = cats.count("incompatible"), cats.count("compatible")
    if n_inc >= MIN_FAMILY and n_com == 0:
        return "replicated: fixed two-type mixing rejected"
    if n_com >= MIN_FAMILY and n_inc == 0:
        return "fixed two-type mixing supported"
    if (n_inc >= MIN_FAMILY or n_com >= MIN_FAMILY) and n_inc and n_com:
        return "mixed"
    return "unresolved"


def main():
    if not DATA.exists():
        sys.exit("no data: run 03_data_ingestion/acquire_lens_confirmation.py")
    d = pd.read_parquet(DATA)
    rng = np.random.default_rng(SEED); rng_b = np.random.default_rng(SEED + 1); rng_c = np.random.default_rng(SEED + 2)
    res = {"registered": "2026-10-07, committed before any value was fetched", "seed": SEED, "B": MM.B_REAL,
           "members": {m: {f"{e}|{v}": int(g.member.nunique()) for (e, v), g in d[d.model == m].groupby(["experiment", "var"])}
                       for m in UNITS if (d.model == m).any()},
           "cases": {}}
    law = {"C3": {"tas": {}, "pr": {}}, "C4": {"tas": {}, "pr": {}}}
    law_in = {"C3": {"tas": [], "pr": []}, "C4": {"tas": [], "pr": []}}
    for m in UNITS:
        dm = d[d.model == m]
        if dm.empty:
            continue
        for var in ("tas", "pr"):
            for season in ("JJA", "DJF"):
                side = CL.SIDE[(var, season)]; q = 2 / 3 if side == "hot" else 1 / 3
                for box in CL.BOXES:
                    s = CL.seasonal(dm, var, box, season)
                    key = f"{m}|{var}|{box}_{season}"
                    arms = [arm_matrix(s, lv) for lv in LEVELS]
                    c = MM.tests(arms, MM.B_REAL, rng)
                    sd0 = np.sqrt(np.nanmean(np.nanvar(arms[0][1], axis=0, ddof=1)))
                    mu = [np.nanmean(a[1]) for a in arms]
                    c["forced_change_sd"] = round(float(abs(mu[2] - mu[0]) / sd0), 3)
                    u = {"C1C2": c}
                    L, M, H = (CL.arm_years(s, lv) for lv in LEVELS)
                    for tag, (ref, frc) in (("C3", (M, H)), ("C4", (L, M))):
                        blk, vals = law_block(ref, frc, q, side, rng_b, rng_c)
                        if blk:
                            u[tag] = blk
                            law[tag][var].setdefault(UNITS[m], []).append(vals); law_in[tag][var].append(blk["inside_range"])
                    res["cases"][key] = u
                    print(key, {k: c[k] for k in ("p_T", "rho", "p_cdf", "forced_change_sd")},
                          {t: (u[t]["reading"], u[t]["law"]["mean"]) for t in ("C3", "C4") if t in u}, flush=True)
    summ = {}
    for var in ("tas", "pr"):
        ks = [k for k in res["cases"] if f"|{var}|" in k]
        for test in ("p_T", "p_cdf"):
            ps = [res["cases"][k]["C1C2"][test] for k in ks]; adj = MM.holm(ps); cats = []
            for k, a in zip(ks, adj):
                c = res["cases"][k]["C1C2"]; ci = c.get("rho_ci95")
                cat = ("incompatible" if a < 0.05 else
                       "compatible" if ci and -TOL_RHO <= ci[0] and ci[1] <= TOL_RHO else "unresolved")
                c[test + "_holm"] = round(float(a), 4); c[test.replace("p_", "") + "_category"] = cat; cats.append(cat)
            summ[f"{'C1' if test == 'p_T' else 'C2'}|{var}"] = {
                "n_cases": len(ks), "n_incompatible": cats.count("incompatible"), "n_compatible": cats.count("compatible"),
                "n_unresolved": cats.count("unresolved"), "n_raw_p_below_0.05": int(sum(p < 0.05 for p in ps)),
                "rho_median": round(float(np.nanmedian([res["cases"][k]["C1C2"]["rho"] for k in ks])), 3),
                "reading": family_reading(cats)}
        for tag in ("C3", "C4"):
            if law_in[tag][var]:
                p = CL.pooled_by_model(law[tag][var], rng_b); p["reading_pooled"] = CL.reading(p)
                share = float(np.mean(law_in[tag][var]))
                rds = [res["cases"][k][tag]["reading"] for k in ks if tag in res["cases"][k]]
                summ[f"{tag}|{var}"] = {"pooled": p, "n_cases": len(law_in[tag][var]), "share_inside_calibration": round(share, 4),
                                        "readings": {r: rds.count(r) for r in ("holds", "inconclusive", "fails")},
                                        "reading": CL.t1_reading(p["reading_pooled"], share)}
    res["summary"] = summ
    print(json.dumps(summ, indent=1))
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
