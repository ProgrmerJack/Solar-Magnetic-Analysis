#!/usr/bin/env python3
"""
prospective_next_ssw.py
=======================
A REGISTERED PROSPECTIVE TEST OF THE SHIFTED NULL ON THE NEXT SSWs.

Registered 2026-10-02, committed (with the prediction table it writes) before the
onset of any SSW it scores. Every result in the paper that rests on observed
outcomes is retrospective; this is the first test whose outcomes do not yet exist.

WHAT IS PREDICTED, AND WHY IT IS FIXED NOW
  The shifted null forecasts an SSW's surface outcome from its onset date alone:
  event-free winter dates within +-10 days of the onset's day of year (1959-2022,
  ERA5 via WeatherBench 2; criterion_transport_cv's candidates and estimator, here
  trained on all 39 catalogued events), displaced by the mean SSW shift. Because
  nothing else enters, the forecast for every possible onset day (1 November - 31
  March) is computed and committed NOW, in the table this script writes
  (results/current/9_literature/prospective_next_ssw.json, "predictions"). The
  scorer recomputes the ensembles from the same frozen inputs and refuses to score
  if they differ from the committed table.

EVENTS
  Every major Northern Hemisphere SSW with onset from 1 November 2026 to 31 March
  2030: the project's Charlton-Polvani detector (ensemble_precursor.detect_ssw,
  validated on the NOAA compendium) on ERA5T daily-mean zonal-mean u(10 hPa, 60N)
  from the CDS (acquire_era5_realtime_cds.py), continued from
  era5_u10_60N_daily_cds.parquet. An onset is scored only once the detector's
  final-warming test can be applied (data to 30 April) and passes; candidates that
  turn out to be final warmings are listed and not scored. Southern Hemisphere
  events are not covered (the downward criterion is defined for the NH).

OUTCOMES (ERA5T, same reductions as the training data; acquire_era5_realtime_cds)
  NAM: polar-cap 1000 and 150 hPa height, 00 UTC, as anomalies from the training
  file's day-of-year mean and s.d. (exactly as acquire_era5_nam.py standardises).
  Temperature: northern-Eurasian 2 m temperature anomaly from the 1959-2022
  smoothed day-of-year mean (criterion_regional_contrast.t_anomalies).
  CONSISTENCY GATE, fixed now: on 2021-12-01..2022-03-31, covered by both sources,
  the CDS series must correlate >= 0.99 with the training series and differ by an
  RMS <= 0.10 of the day-of-year s.d., for nam_1000, nam_150 and NEURASIA;
  otherwise nothing is scored until the reduction is fixed, and that is reported.
  TREND, fixed now: temperature ensemble members are moved to the target winter by
  the linear trend of November-March event-free anomalies 1959-2022 (OLS on winter
  year), because the climate has warmed since most candidate winters. The NAM is
  not adjusted.

FORECASTS COMPARED (each event)
  shifted     candidates displaced by the SSW shift (NAM 1000/150, W = 52) and,
              for temperature, by the SSW shift of the days 8-24 mean
  climatology the same candidates, not displaced
  rate        the constant class rate of the 39 training events (label only)

SCORES, per event and cumulatively
  S1 PRIMARY  CRPS of the days 8-52 mean 1000 hPa NAM, climatology minus shifted
  S2          Brier score of the downward label (published surface conditions,
              days 8-52) for shifted, climatology and rate
  S3          CRPS of northern-Eurasian temperature, days 8-24, climatology minus
              shifted
  S4          PIT of the observed NAM in the shifted ensemble
READING, fixed now (cumulative, reported after every scored event; decisive only
  from five events): mean S1 gain with a 95% interval from resampling events that
  excludes zero -- positive: "prospectively supported"; negative: "prospectively
  refuted: one shifted population does not describe new SSW outcomes" --
  otherwise "not yet decided". The same reading is applied to S3 as a secondary.

USAGE
  prospective_next_ssw.py                     build and save the prediction table
  prospective_next_ssw.py --gate              check the CDS overlap (consistency gate)
  prospective_next_ssw.py --detect            list onsets in the real-time series
  prospective_next_ssw.py --score             score every confirmed onset with a
                                              complete outcome window

Output: results/current/9_literature/prospective_next_ssw.json
"""
import hashlib
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import criterion_sweep as CS                         # noqa: E402
import criterion_regional_contrast as RC             # noqa: E402

sys.path.insert(0, str(CS.PHYS))
from ensemble_precursor import detect_ssw            # noqa: E402

G, load_catalogue = CS.G, CS.load_catalogue
NAME = "prospective_next_ssw"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
RESULTS = CS.RESULTS
OUT = RESULTS / f"{NAME}.json"
ING = HERE.parents[0] / "03_data_ingestion"
RT_FILE = ING / "era5_realtime_cds.parquet"
CHECK_FILE = ING / "era5_realtime_cds_check.parquet"
U_HIST = ING / "era5_u10_60N_daily_cds.parquet"
N_SHIFT_SETS, HALFWIN, N_BOOT = 200, 10, 10000
FIRST, LAST = pd.Timestamp("2026-11-01"), pd.Timestamp("2030-03-31")
REF_DAYS = pd.date_range("2026-11-01", "2027-03-31")         # one prediction per onset day
GATE_R, GATE_RMS = 0.99, 0.10
T_DAYS = (8, 24)


def crps_ens(x, y):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.mean(np.abs(x - y)) - 0.5 * np.mean(np.abs(x[:, None] - x[None, :])))


def winter(t):
    return t.year if t.month <= 6 else t.year + 1


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]


# ------------------------------------------------------------------ the frozen model
def frozen_model():
    rng = np.random.default_rng(SEED)
    e = pd.read_parquet(CS.ERA5)
    tan = RC.t_anomalies()
    allev = load_catalogue("primary")
    real = allev[(allev >= e.index.min() + pd.Timedelta(days=70)) & (allev <= e.index.max() - pd.Timedelta(days=70))]
    clean = G.zone_free_index(e.index, allev)
    Mr = {c: CS.window_matrix(e[c], real) for c in ("nam_1000", "nam_150")}
    Mc = {c: CS.window_matrix(e[c], clean) for c in ("nam_1000", "nam_150")}
    Tr = RC.t_window(tan, real, *T_DAYS)[:, 0]
    Tc = RC.t_window(tan, clean, *T_DAYS)[:, 0]
    dr = np.array([t.dayofyear for t in real]); dc = np.array([t.dayofyear for t in clean])
    by = G.build_by_doy(clean)
    sets = [G.draw_clean(clean, dr, rng, by) for _ in range(N_SHIFT_SETS)]
    pos = {t: i for i, t in enumerate(clean)}
    idx = [np.array([pos[t] for t in s_]) for s_ in sets]
    shift = {c: float(np.nanmean(CS.wmean(Mr[c], 52)) - np.nanmean(np.concatenate([CS.wmean(Mc[c][j], 52) for j in idx])))
             for c in Mr}
    shift["T_NEURASIA_d8_24"] = float(np.nanmean(Tr) - np.nanmean(np.concatenate([Tc[j] for j in idx])))
    lab_r, ok_r = CS.classify(Mr["nam_1000"], Mr["nam_150"], 52, 0.0, 0.5, False)
    rate = float(lab_r[ok_r].mean())
    # linear trend of Nov-Mar event-free temperature anomalies on winter year
    nd = np.array([t.month in (11, 12, 1, 2, 3) for t in clean]) & np.isfinite(Tc)
    yrs = np.array([winter(t) for t in clean])
    beta = float(np.polyfit(yrs[nd], Tc[nd], 1)[0])
    return {"e": e, "tan": tan, "clean": clean, "Mc": Mc, "Tc": Tc, "dc": dc, "yrs": yrs,
            "shift": shift, "rate": rate, "beta_T_per_year": beta, "n_events": int(len(real))}


def ensemble(m, onset):
    """Shifted and climatological ensembles for an onset date."""
    d = pd.Timestamp(onset).dayofyear
    dd = np.abs(m["dc"] - d); dd = np.minimum(dd, 366 - dd)
    cand = np.flatnonzero(dd <= HALFWIN)
    sh = m["shift"]
    lab_s, ok_s = CS.classify(m["Mc"]["nam_1000"][cand], m["Mc"]["nam_150"][cand], 52, 0.0, 0.5, False, sh["nam_1000"])
    lab_c, ok_c = CS.classify(m["Mc"]["nam_1000"][cand], m["Mc"]["nam_150"][cand], 52, 0.0, 0.5, False)
    nam = CS.wmean(m["Mc"]["nam_1000"][cand], 52)
    T = m["Tc"][cand] + m["beta_T_per_year"] * (winter(pd.Timestamp(onset)) - m["yrs"][cand])
    return {"cand": cand, "p_dw_shifted": float(lab_s[ok_s].mean()), "p_dw_clim": float(lab_c[ok_c].mean()),
            "nam_clim": nam[np.isfinite(nam)], "nam_shifted": nam[np.isfinite(nam)] + sh["nam_1000"],
            "T_clim": T[np.isfinite(T)], "T_shifted": T[np.isfinite(T)] + sh["T_NEURASIA_d8_24"]}


def q(x):
    return [round(float(v), 3) for v in np.quantile(x, [0.1, 0.5, 0.9])]


def table(m):
    rows = {}
    for d in REF_DAYS:
        en = ensemble(m, d)
        rows[d.strftime("%m-%d")] = {"n": int(len(en["cand"])), "p_dw_shifted": round(en["p_dw_shifted"], 4),
                                     "p_dw_clim": round(en["p_dw_clim"], 4),
                                     "nam_shifted_q10_50_90": q(en["nam_shifted"]),
                                     "nam_clim_q10_50_90": q(en["nam_clim"]),
                                     "T_shifted_q10_50_90_K_winter2027": q(en["T_shifted"]),
                                     "T_clim_q10_50_90_K_winter2027": q(en["T_clim"])}
    return rows


# ------------------------------------------------------------------ real-time data
def standardise(rt, e):
    """NAM of real-time heights against the training file's day-of-year statistics."""
    out = pd.DataFrame(index=rt.index)
    doy_e = e.index.dayofyear
    for lev in (1000, 150):
        h = e[f"z{lev}_m"]
        mu = h.groupby(doy_e).mean(); sd = h.groupby(doy_e).std()
        d = rt.index.dayofyear
        out[f"nam_{lev}"] = -((rt[f"z{lev}_m"].values - mu.reindex(d).values) / sd.reindex(d).values)
    return out


def gate(m):
    if not CHECK_FILE.exists():
        return {"passed": False, "reason": "no overlap file: run acquire_era5_realtime_cds.py --check"}
    ck = pd.read_parquet(CHECK_FILE)
    nam = standardise(ck, m["e"])
    res, ok = {}, True
    for lev in ("nam_1000", "nam_150"):
        a = m["e"][lev].reindex(nam.index); b = nam[lev]
        k = a.notna() & b.notna()
        r = float(np.corrcoef(a[k], b[k])[0, 1]); rms = float(np.sqrt(np.mean((a[k] - b[k]) ** 2)))   # s.d. units
        res[lev] = {"r": round(r, 4), "rms_sd_units": round(rms, 4), "n_days": int(k.sum())}
        ok &= r >= GATE_R and rms <= GATE_RMS
    t_tr = pd.read_parquet(RC.T_FILE).set_index("date")["NEURASIA"]; t_tr.index = pd.to_datetime(t_tr.index)
    a = t_tr.reindex(ck.index); b = ck["NEURASIA"]; k = a.notna() & b.notna()
    tan = m["tan"]["NEURASIA"]
    sd_doy = tan.groupby(tan.index.dayofyear).std()
    rms = float(np.sqrt(np.mean(((a[k] - b[k]) / sd_doy.reindex(a[k].index.dayofyear).values) ** 2)))
    r = float(np.corrcoef(a[k], b[k])[0, 1])
    res["NEURASIA"] = {"r": round(r, 4), "rms_sd_units": round(rms, 4), "bias_K": round(float((b[k] - a[k]).mean()), 3),
                       "n_days": int(k.sum())}
    ok &= r >= GATE_R and rms <= GATE_RMS
    res["passed"] = bool(ok)
    return res


def onsets():
    u = pd.read_parquet(U_HIST).set_index("date")["u10_60N"]; u.index = pd.to_datetime(u.index)
    if RT_FILE.exists():
        rt = pd.read_parquet(RT_FILE)["u10_60N"].dropna()
        u = pd.concat([u[~u.index.isin(rt.index)], rt]).sort_index()
    found = [pd.Timestamp(t) for t in detect_ssw(u.values, u.index)]
    last = u.index.max()
    out = []
    for t in found:
        if FIRST <= t <= LAST:
            end = pd.Timestamp(year=winter(t), month=4, day=30)
            out.append({"onset": str(t.date()), "final_warming_test_applicable": bool(last >= end)})
    return out, str(last.date())


def score(m, committed):
    rt = pd.read_parquet(RT_FILE)
    nam = standardise(rt, m["e"])
    base = pd.read_parquet(RC.T_FILE).set_index("date")[["NEURASIA"]]; base.index = pd.to_datetime(base.index)
    b = base[(base.index.year >= 1959) & (base.index.year <= 2022)]
    doy = b.groupby(b.index.dayofyear).mean().reindex(range(1, 367))
    doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    doy.index = range(1, 367)
    tan_rt = rt[["NEURASIA"]] - doy.reindex(rt.index.dayofyear).values
    ev, _ = onsets()
    scored = []
    for o in ev:
        if not o["final_warming_test_applicable"]:
            continue
        t0 = pd.Timestamp(o["onset"])
        key = t0.strftime("%m-%d")
        en = ensemble(m, t0)
        if committed and committed[key]["p_dw_shifted"] != round(en["p_dw_shifted"], 4):
            sys.exit(f"recomputed prediction for {key} differs from the committed table -- refusing to score")
        win = pd.date_range(t0 + pd.Timedelta(days=CS.LO), t0 + pd.Timedelta(days=52))
        y = nam.reindex(win)
        tw = tan_rt["NEURASIA"].reindex(pd.date_range(t0 + pd.Timedelta(days=T_DAYS[0]), t0 + pd.Timedelta(days=T_DAYS[1])))
        if y.isna().any().any() or tw.isna().any():
            scored.append({"onset": o["onset"], "status": "outcome window incomplete"})
            continue
        yn = float(y["nam_1000"].mean())
        dw = bool(yn < 0 and (y["nam_1000"] < 0).mean() > 0.5)
        ty = float(tw.mean())
        scored.append({"onset": o["onset"], "status": "scored", "observed_nam_d8_52": round(yn, 4), "observed_dw": dw,
                       "observed_T_d8_24_K": round(ty, 3),
                       "S1_crps_gain": round(crps_ens(en["nam_clim"], yn) - crps_ens(en["nam_shifted"], yn), 4),
                       "S2_brier": {"shifted": round((en["p_dw_shifted"] - dw) ** 2, 4),
                                    "climatology": round((en["p_dw_clim"] - dw) ** 2, 4),
                                    "rate": round((m["rate"] - dw) ** 2, 4)},
                       "S3_crps_gain_T": round(crps_ens(en["T_clim"], ty) - crps_ens(en["T_shifted"], ty), 4),
                       "S4_pit": round(float((np.sum(en["nam_shifted"] < yn) + 0.5) / (len(en["nam_shifted"]) + 1)), 4)})
    done = [s for s in scored if s["status"] == "scored"]
    cum = {"n_scored": len(done)}
    rng = np.random.default_rng(SEED)
    for k in ("S1_crps_gain", "S3_crps_gain_T"):
        v = np.array([s[k] for s in done])
        if len(v):
            bs = [rng.choice(v, len(v)).mean() for _ in range(N_BOOT)]
            ci = [float(x) for x in np.quantile(bs, [0.025, 0.975])]
            verdict = ("not yet decided" if len(v) < 5 else
                       "prospectively supported" if ci[0] > 0 else
                       "prospectively refuted" if ci[1] < 0 else "not yet decided")
            cum[k] = {"mean": round(float(v.mean()), 4), "ci95": [round(x, 4) for x in ci], "reading": verdict}
    return {"events": scored, "cumulative": cum}


def main():
    argv = sys.argv[1:]
    m = frozen_model()
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    if "--gate" in argv:
        prev["consistency_gate"] = gate(m)
        print(json.dumps(prev["consistency_gate"], indent=1))
    elif "--detect" in argv:
        ev, last = onsets()
        prev["detected"] = {"through": last, "onsets": ev}
        print(prev["detected"])
    elif "--score" in argv:
        g = prev.get("consistency_gate") or gate(m)
        if not g.get("passed"):
            sys.exit(f"consistency gate not passed: {g}")
        prev["scores"] = score(m, prev.get("predictions"))
        print(json.dumps(prev["scores"], indent=1))
    else:
        prev = {"registered": "2026-10-02; prediction table committed before any scored onset",
                "seed": SEED, "n_shift_sets": N_SHIFT_SETS, "halfwin_days": HALFWIN, "n_boot": N_BOOT,
                "events_window": [str(FIRST.date()), str(LAST.date())],
                "inputs": {str(CS.ERA5.relative_to(CS.ROOT)): sha(CS.ERA5),
                           str(RC.T_FILE.relative_to(CS.ROOT)): sha(RC.T_FILE)},
                "n_training_events": m["n_events"],
                "shift": {k: round(v, 4) for k, v in m["shift"].items()},
                "constant_rate": round(m["rate"], 4), "beta_T_K_per_year": round(m["beta_T_per_year"], 5),
                "predictions": table(m)}
        p = prev["predictions"]
        for k in ("11-15", "12-15", "01-15", "02-15", "03-15"):
            print(k, p[k])
        print("shift", prev["shift"], "rate", prev["constant_rate"], "beta", prev["beta_T_K_per_year"])
    OUT.write_text(json.dumps(prev, indent=1), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
