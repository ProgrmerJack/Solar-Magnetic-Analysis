#!/usr/bin/env python3
"""
vortex_threshold_continuity.py
==============================
IS AN SSW ITSELF A THRESHOLD ON A CONTINUUM? A MODEL-FREE TEST AT THE REVERSAL.

Plan approved 2026-10-01 (user: "I do approve any new kind of analyses"); design from
research_notes/SSW paper Nature Geo strengthening/C_observational_causal.md section 5.
Committed before any outcome was related to the running variable. Only the wind
(running variable) and the unconditional outcome s.d. had been computed (power).

QUESTION
  The SSW definition is a threshold: the zonal-mean wind at 10 hPa, 60N reversing
  below zero. If "SSW" marked a distinct kind of event, surface outcomes would jump
  at u_min = 0 among otherwise similar vortex disruptions; if SSWs are the extreme
  of a continuum of weakenings, the outcome changes smoothly with how far the wind
  falls. This is the observational counterpart of the label result (the downward
  label as a threshold on one shifted distribution), and needs no model.

UNITS (NCEP-R1 daily u(10 hPa, 60N), 1958-01-01 .. 2025-04-30; winter 2025-26
excluded to protect s2s_heldout2026_test.py)
  1. minimum day t in 1 Nov - 31 Mar where u(t) is the minimum over [t-20, t+20]
     (first day on ties);
  2. deceleration: max u over [t-30, t] minus u(t) >= 15 m/s;
  3. recovery: 10 consecutive days with u > max(u_min, 0) before 30 April (equal to
     the Charlton-Polvani final-warming rule for reversals, continuous at 0);
  running variable X = u_min (m/s), cutoff 0, D = 1(X < 0); key date t_min.

OUTCOMES (days t_min + 8 .. + 52 unless stated)
  Y1 (primary) CPC daily AO, mean;  Y2 ERA5 NAM1000 mean (units to 2022);
  Y3 ERA5 northern-Eurasian 2 m temperature anomaly (WB2 spliced with ARCO-ERA5 as
     in s2s_heldout_test.py; units to March 2024), K;
  Y4 the downward label (conditions 1-2 of the criterion on the AO: window mean
     < 0 and > 50% of days < 0), applied to every unit.
  Placebo outcome: AO mean over the 30 days before the pre-minimum maximum (before
  the deceleration began). Balance covariates: day of year (sin, cos), era (pre/post
  1979), deceleration size.

ESTIMATES
  E1 DOSE (primary for the size of the effect): OLS slope of Y on X over all units,
     with the placebo AO, day-of-year sin/cos and era as covariates; winter-cluster
     bootstrap (10,000) interval; per 10 m/s.
  E2 GLOBAL JUMP (primary for continuity): the coefficient of D in the same
     regression with separate slopes on each side (piecewise linear); interval as
     E1; reported with the equivalence bound |jump| < upper 90% limit.
  E3 LOCAL RD: rdrobust (Calonico-Cattaneo-Titiunik) local linear, triangular kernel,
     MSE-optimal bandwidth, robust bias-corrected CI, clustered by winter; and fixed
     bandwidths 5, 7.5, 10, 15, 20 m/s.
  E4 LOCAL RANDOMISATION: difference in means for |X| < w, w in {2.5, 5}; Fisher
     permutation p (10,000).
  Holm over Y1-Y4 within E2 (the continuity family). Placebo family separately.
FALSIFICATION (fixed now): placebo cutoffs at -10, -5, +5, +10 m/s (same-side units
only); donut excluding |X| < 1 and < 2 m/s; placebo outcome at 0; balance of
covariates at 0 (E3 on each covariate).
READING (fixed now): one population predicts E1 slope != 0 (outcome more negative /
colder as the wind falls further) and E2 jump consistent with 0. A jump with p <
0.05 (Holm) whose placebo and donut checks pass would mean the reversal marks a
distinct kind of event. A null jump is reported only as a bound, because the design
cannot detect jumps smaller than about the whole SSW effect in observations.

CMIP6 CALIBRATION (minimally model-dependent)
  The identical episode rule on each member's daily u(10 hPa, 60N) and annular-mode
  index (Y1 analogue, standardised by the member's daily s.d.) for the 20 CMIP6
  members: E1, E2 with thousands of units (the model's answer), and the size and
  power of the observational E2 from 1,000 random draws of 68 consecutive winters
  (power for a planted jump equal to the observed whole-SSW shift, -0.8 s.d.).

REVISION-3 ADDENDUM: ERA5 RUNNING VARIABLE, 1940-2025 (approved 2026-10-01;
committed BEFORE the ERA5 wind series and any ERA5 surface value before 1959 were
retrieved). Written to the same JSON under "era5_1940_2025"; the registered NCEP
numbers above are not changed and are reproduced first.
  Why: several NCEP units just above 0 (e.g. January 1977, February 1963, January
  1968, February/March 1981, February 2002) are reversals in other reanalyses
  (Butler et al. 2017 compendium), so measurement error in the running variable
  sits exactly at the cutoff; and ERA5 adds 1940-1957.
  Units: the identical episode rule on the ERA5 u(10 hPa, 60N) daily mean
  (acquire_era5_u10_cds.py), minima 1 November 1940 - 31 March 2025.
  Outcomes, days +8..+52: Y2E the polar-cap NAM proxy, -(ERA5 polar-cap mean
  sea-level pressure at 00 UTC minus its 31-day-smoothed day-of-year mean)
  divided by its November-March daily s.d. (WeatherBench 2 1959-April 2022,
  ARCO-ERA5 otherwise with the overlap offset of s2s_heldout_test.splice; the
  primary outcome, available for all years); Y1 CPC AO (from 1950); Y3
  northern-Eurasian 2 m temperature, anomaly from the smoothed day-of-year mean
  with a linear trend over 1940-2025 removed; Y4 the downward label on Y2E.
  Covariates as before plus a linear year term (trend). Estimates E1-E4, Holm over
  the outcomes, placebo cutoffs, donut, balance, exactly as registered; the
  falsification rule for the local estimates is unchanged.
  Also reported: the NCEP-ERA5 agreement of u_min for units in both, and how many
  units change side of 0. Sensitivity: minima from 1946 only. CMIP6 calibration is
  not repeated (it does not depend on the reanalysis).

Output: results/current/5_mechanism/vortex_threshold_continuity.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RESULTS = ROOT / "results" / "current" / "5_mechanism"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ING))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
import ensemble_precursor as EP                     # noqa: E402

NAME = "vortex_threshold_continuity"
SEED = zlib.crc32(NAME.encode()) % (2 ** 32)
N_BOOT, N_PERM, N_CAL = 10000, 10000, 1000
HALF, DROP, REC = 20, 15.0, 10
WIN = (8, 52)
BWS = (5.0, 7.5, 10.0, 15.0, 20.0)
NCEP = ING / "_ncep_reduced"
AO = ROOT / "data" / "processed" / "atmospheric" / "ao_daily_cpc.txt"
NAM = ING / "era5_nam_daily.parquet"
T2M = ING / "era5_t2m_regions_daily_spliced.parquet"
RAW = ING / "raw" / "cmip6"
END = pd.Timestamp("2025-04-30")


def episodes(u):
    """Weakening episodes from a daily series u (DatetimeIndex)."""
    vals, idx = u.values, u.index
    rows = []
    for i, t in enumerate(idx):
        if t.month not in (11, 12, 1, 2, 3):
            continue
        lo, hi = max(0, i - HALF), min(len(u), i + HALF + 1)
        w = vals[lo:hi]
        if not np.isfinite(vals[i]) or vals[i] != np.nanmin(w) or int(np.nanargmin(w)) != i - lo:
            continue
        j0 = max(0, i - 30)
        jmax = j0 + int(np.nanargmax(vals[j0:i + 1]))
        if vals[jmax] - vals[i] < DROP:
            continue
        end = pd.Timestamp(year=t.year if t.month <= 4 else t.year + 1, month=4, day=30)
        post = u[t + pd.Timedelta(days=1):end].values
        thr = max(vals[i], 0.0)
        run, rec = 0, False
        for v in post:
            run = run + 1 if v > thr else 0
            if run >= REC:
                rec = True
                break
        if not rec:
            continue
        rows.append({"date": t, "winter": t.year if t.month <= 4 else t.year + 1,
                     "X": float(vals[i]), "t_max": idx[jmax], "drop": float(vals[jmax] - vals[i])})
    return pd.DataFrame(rows)


def window_mean(s, t, a, b, need=0.9):
    w = s[(s.index >= t + pd.Timedelta(days=a)) & (s.index <= t + pd.Timedelta(days=b))]
    return float(w.mean()) if len(w) >= need * (b - a + 1) else np.nan


def label(s, t):
    w = s[(s.index >= t + pd.Timedelta(days=WIN[0])) & (s.index <= t + pd.Timedelta(days=WIN[1]))]
    if len(w) < 0.9 * (WIN[1] - WIN[0] + 1):
        return np.nan
    return float((w.mean() < 0) and ((w < 0).mean() > 0.5))


def design(df, y):
    d = df.dropna(subset=[y, "plac"]).copy()
    doy = d["date"].dt.dayofyear.values
    X = np.column_stack([np.ones(len(d)), d["X"].values, d["plac"].values,
                         np.sin(2 * np.pi * doy / 365.25), np.cos(2 * np.pi * doy / 365.25),
                         (d["winter"] >= 1979).astype(float)])
    return d, X


def ols(Xm, y):
    return np.linalg.lstsq(Xm, y, rcond=None)[0]


def dose_and_jump(df, y, rng, n_boot=N_BOOT):
    d, X = design(df, y)
    D = (d["X"].values < 0).astype(float)
    Xj = np.column_stack([X, D, D * d["X"].values])
    yv = d[y].values
    b_dose, b_jump = ols(X, yv), ols(Xj, yv)
    win = d["winter"].values
    uw = np.unique(win)
    bd, bj = [], []
    for _ in range(n_boot):
        pick = rng.choice(uw, len(uw))
        i = np.concatenate([np.flatnonzero(win == w) for w in pick])
        if len(np.unique(D[i])) < 2:
            continue
        bd.append(ols(X[i], yv[i])[1]); bj.append(ols(Xj[i], yv[i])[6])
    bd, bj = np.array(bd), np.array(bj)
    return {"n": int(len(d)), "n_below": int(D.sum()),
            "dose_slope_per_10ms": round(float(10 * b_dose[1]), 4),
            "dose_ci95": [round(float(10 * q), 4) for q in np.percentile(bd, [2.5, 97.5])],
            "dose_p_two_sided": round(float(min(1, 2 * min(np.mean(bd <= 0), np.mean(bd >= 0)))), 4),
            "jump": round(float(b_jump[6]), 4),
            "jump_ci95": [round(float(q), 4) for q in np.percentile(bj, [2.5, 97.5])],
            "jump_ci90": [round(float(q), 4) for q in np.percentile(bj, [5, 95])],
            "jump_p_two_sided": round(float(min(1, 2 * min(np.mean(bj <= 0), np.mean(bj >= 0)))), 4)}


def rd_local(df, y, h=None, cutoff=0.0, sub=None):
    from rdrobust import rdrobust
    d = df.dropna(subset=[y]) if sub is None else sub.dropna(subset=[y])
    try:
        kw = dict(y=d[y].values, x=d["X"].values, c=cutoff, kernel="triangular",
                  cluster=d["winter"].values, vce="hc1")
        if h is not None:
            kw["h"] = h
        r = rdrobust(**kw)
        coef = float(np.asarray(r.coef).ravel()[0])
        ci = np.asarray(r.ci)[2]                    # robust bias-corrected row
        pv = float(np.asarray(r.pv).ravel()[2])
        bw = float(np.asarray(r.bws).ravel()[0])
        # rdrobust reports the right-minus-left jump; D = 1 is the left (X < c) side
        return {"jump_left_minus_right": round(-coef, 4),
                "rbc_ci95_left_minus_right": [round(float(-ci[1]), 4), round(float(-ci[0]), 4)],
                "rbc_p": round(pv, 4), "h": round(bw, 2),
                "n_eff": [int(x) for x in np.asarray(r.N_h).ravel()]}
    except Exception as ex:
        return {"error": str(ex)[:160]}


def local_randomisation(df, y, w, rng):
    d = df.dropna(subset=[y])
    d = d[np.abs(d["X"]) < w]
    a, b = d[d["X"] < 0][y].values, d[d["X"] >= 0][y].values
    if len(a) < 3 or len(b) < 3:
        return {"w": w, "n": [len(a), len(b)], "note": "too few units"}
    obs = a.mean() - b.mean()
    allv = np.concatenate([a, b])
    perm = np.array([(lambda p: p[:len(a)].mean() - p[len(a):].mean())(rng.permutation(allv))
                     for _ in range(N_PERM)])
    return {"w": w, "n": [int(len(a)), int(len(b))], "diff_below_minus_above": round(float(obs), 4),
            "perm_p_two_sided": round(float(np.mean(np.abs(perm) >= abs(obs))), 4)}


def holm(ps):
    idx = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, run = [None] * len(ps), 0.0
    for r, i in enumerate(idx):
        run = max(run, min(1.0, (len(ps) - r) * ps[i]))
        adj[i] = run
    return adj


def observations(rng):
    u = pd.concat([pd.read_parquet(f)["uwnd_ms_60N_10hPa"] for f in sorted(NCEP.glob("uwnd.*.10.parquet"))])
    u.index = pd.to_datetime(u.index).tz_localize(None).normalize()
    u = u[~u.index.duplicated()].sort_index()
    u = u[u.index <= END]
    E = episodes(u)
    E = E[E["winter"] <= 2025].reset_index(drop=True)
    # fixed-width CPC file; one line (2003-04-30) fuses the day with the missing
    # code ("30-99.000"), so parse by regular expression and treat -99 as missing
    import re
    recs = [re.match(r"\s*(\d{4})\s+(\d+)\s+(\d+)\s*(-?\d+\.\d+)", ln) for ln in AO.read_text().splitlines()]
    recs = [m.groups() for m in recs if m]
    ao = pd.Series([float(v) for *_, v in recs],
                   index=pd.to_datetime([f"{y}-{m}-{d}" for y, m, d, _ in recs]))
    ao = ao.where(ao > -90)
    nam = pd.read_parquet(NAM)["nam_1000"]
    t2 = pd.read_parquet(T2M).set_index("date")["NEURASIA"]
    t2.index = pd.to_datetime(t2.index)
    t2a = t2 - t2.groupby(t2.index.dayofyear).transform("mean")
    E["Y1_AO"] = [window_mean(ao, t, *WIN) for t in E["date"]]
    E["Y2_NAM1000"] = [window_mean(nam, t, *WIN) for t in E["date"]]
    E["Y3_T_NEURASIA"] = [window_mean(t2a, t, *WIN) for t in E["date"]]
    E["Y4_downward"] = [label(ao, t) for t in E["date"]]
    E["plac"] = [window_mean(ao, tm, -30, -1) for tm in E["t_max"]]
    return E


def cmip6(rng):
    rows = []
    for k, f in enumerate(sorted(RAW.glob("*_zm.nc"))):
        try:
            m = EP.load_member(f)
        except Exception:
            continue
        u = m["u10"].dropna()
        if len(u) < 3650:
            continue
        am = m["am"]
        an = am - am.groupby(am.index.dayofyear).transform("mean")
        an = an / an.std()
        E = episodes(u)
        if E.empty:
            continue
        E["Y1"] = [window_mean(an, t, *WIN) for t in E["date"]]
        E["plac"] = [window_mean(an, tm, -30, -1) for tm in E["t_max"]]
        E["member"] = f.stem
        E["winter"] = E["winter"] + 100000 * k          # winters distinct across members
        rows.append(E)
    return pd.concat(rows, ignore_index=True)


def main():
    rng = np.random.default_rng(SEED)
    res = {"plan_approved": "2026-10-01", "seed": SEED, "n_boot": N_BOOT, "n_perm": N_PERM,
           "episode_rule": {"half_window_days": HALF, "deceleration_ms": DROP, "recovery_days": REC},
           "outcome_window": list(WIN), "inputs": ["03_data_ingestion/_ncep_reduced/uwnd.*.10.parquet",
                                                    str(AO.relative_to(ROOT)), str(NAM.relative_to(ROOT)),
                                                    str(T2M.relative_to(ROOT))]}
    E = observations(rng)
    res["units"] = {"n": int(len(E)), "below_0": int((E.X < 0).sum()), "above_0": int((E.X >= 0).sum()),
                    "years": [int(E.winter.min()), int(E.winter.max())]}
    print("units", res["units"], flush=True)
    # per-unit values, for the figure (09_figures reads, never computes)
    res["unit_table"] = [{"date": str(pd.Timestamp(r.date).date()), "winter": int(r.winter),
                          "X": round(float(r.X), 3),
                          **{k: (None if not np.isfinite(getattr(r, k)) else round(float(getattr(r, k)), 4))
                             for k in ("Y1_AO", "Y2_NAM1000", "Y3_T_NEURASIA", "Y4_downward")}}
                         for r in E.itertuples()]
    ys = ["Y1_AO", "Y2_NAM1000", "Y3_T_NEURASIA", "Y4_downward"]
    res["observations"] = {}
    for y in ys:
        r = {"E1_E2": dose_and_jump(E, y, rng), "E3_mse": rd_local(E, y),
             "E3_fixed": {str(h): rd_local(E, y, h=h) for h in BWS},
             "E4": [local_randomisation(E, y, w, rng) for w in (2.5, 5.0)],
             "donut": {str(dn): rd_local(E, y, sub=E[np.abs(E.X) >= dn]) for dn in (1.0, 2.0)},
             "placebo_cutoffs": {str(c): rd_local(E, y, cutoff=c, sub=E[(E.X < 0) if c < 0 else (E.X >= 0)])
                                 for c in (-10.0, -5.0, 5.0, 10.0)}}
        res["observations"][y] = r
        print(y, json.dumps(r["E1_E2"]), "| E3", r["E3_mse"], flush=True)
    pj = [res["observations"][y]["E1_E2"]["jump_p_two_sided"] for y in ys]
    res["E2_holm"] = dict(zip(ys, [round(p, 4) for p in holm(pj)]))
    # falsification: placebo outcome and covariate balance at 0
    E["doy_sin"] = np.sin(2 * np.pi * E["date"].dt.dayofyear / 365.25)
    E["era"] = (E["winter"] >= 1979).astype(float)
    res["falsification"] = {"placebo_outcome_AO_pre": rd_local(E, "plac"),
                            "balance_doy_sin": rd_local(E, "doy_sin"), "balance_era": rd_local(E, "era"),
                            "balance_drop": rd_local(E, "drop")}
    print("falsification", res["falsification"], flush=True)
    # CMIP6 calibration
    C = cmip6(rng)
    res["cmip6"] = {"n_units": int(len(C)), "below_0": int((C.X < 0).sum()),
                    "E1_E2": dose_and_jump(C.rename(columns={"Y1": "Y1"}), "Y1", rng, n_boot=2000)}
    print("CMIP6", res["cmip6"], flush=True)
    # observational design size and power: 68-winter blocks within members
    sizes, powers = [], []
    obs_sd = float(np.nanstd(E["Y1_AO"]))
    shift = -0.8
    members = C["member"].unique()
    for _ in range(N_CAL):
        mm = rng.choice(members)
        Cm = C[C.member == mm]
        ws = np.sort(Cm["winter"].unique())
        if len(ws) < 68:
            continue
        s0 = rng.integers(0, len(ws) - 67)
        blk = Cm[Cm["winter"].isin(ws[s0:s0 + 68])].copy()
        if (blk.X < 0).sum() < 5:
            continue
        r0 = dose_and_jump(blk, "Y1", rng, n_boot=200)
        sizes.append(r0["jump_p_two_sided"] < 0.05)
        blk["Y1"] = blk["Y1"] + shift * (blk["X"] < 0)
        r1 = dose_and_jump(blk, "Y1", rng, n_boot=200)
        powers.append(r1["jump_p_two_sided"] < 0.05)
    res["cmip6"]["obs_design_rejection_rate_model_jump"] = round(float(np.mean(sizes)), 3) if sizes else None
    res["cmip6"]["obs_design_power_planted_jump_minus0p8sd"] = round(float(np.mean(powers)), 3) if powers else None
    res["cmip6"]["n_calibration_blocks"] = len(sizes)
    print("CMIP6 calibration", {k: res["cmip6"][k] for k in res["cmip6"] if k.startswith("obs_") or k == "n_calibration_blocks"}, flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "vortex_threshold_continuity.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> vortex_threshold_continuity.json")


if __name__ == "__main__":
    main()
