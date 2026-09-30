#!/usr/bin/env python3
"""
s2s_heldout_test.py
===================
HELD-OUT EVALUATION: DO FORECASTS UNDER-PREDICT THE SHIFT AFTER THE 2023-2024 SSWs?

Plan approved 2026-09-30. This design was committed BEFORE any forecast of model
year 2025 (ECMWF, CMA) or any ERA5 value after 2023-01-10 was retrieved, and
before any statistic at these dates was computed from the CNRM files already on
disk. Whatever it shows is reported.

WHAT HAD BEEN SEEN
  - Published descriptions of the events: Lee, Butler & Manney (2025, Weather 80,
    45-53) describe January 2024 as not exerting a canonical surface influence and
    March 2024 as weakly coupled.
  - The CNRM reforecast files (s2s_cnrm_psl_cap.parquet, s2s_cnrm_t2m_regions.parquet;
    hindcast years 1999-2024) were retrieved for the confirmatory test; their
    2023/2024 starts entered only the leave-one-year-out forecast climatology of
    other years. No outcome, rank or ensemble statistic at the dates below was
    computed or looked at.

EVENTS (the frozen primary catalogue, load_catalogue("primary"); no new detection)
  2023-02-16, 2024-01-16, 2024-03-04: the catalogued onsets after the end of the
  ERA5 series used by the operational tests (polar cap to 2022-04-30, regional
  temperature to 2023-01-10). They were never part of any operational test.

FORECASTS (one model version per system, Dec-Mar starts, as in the P' design)
  cnrm       model version 1 June 2025, hindcast years 1999-2024 (on disk)
  ecmwf2025  ECMWF model year 2025, hindcast years 2005-2024, odd-day starts, 11 members
  cma2025    CMA model year 2025, hindcast years 2010-2024, Jan-Mar starts (none in
             December), 4 members
  Polar-cap msl leads 10-34 and daily-mean 2 m temperature windows 10-33, reduced
  exactly as in P' and the regional test (acquire_s2s_reforecasts.py).

OBSERVATIONS
  The existing ERA5 series (WeatherBench 2) wherever they exist; after their end,
  ARCO-ERA5 reduced to the same grid and regions (acquire_era5_arco_extension.py).
  Before splicing, the mean ARCO-minus-WB2 difference over the overlap (00 UTC
  polar cap: Nov-Apr 2020/21 and 2021/22; regional temperature: those winters plus
  2022-11-01..2023-01-10) is removed, one constant per series; the constant, the
  daily RMS difference and the correlation are reported. If the correlation of
  daily values over the overlap is below 0.99, the splice is reported as failed and
  the held-out test is not run.

STATISTICS (identical to P' and the regional test: s2s_forecast_test.start_stats,
leave-one-year-out anomalies, starts 2-9 days before onset, outcome days +8..+25
(polar cap) and +8..+24 (temperature), sign so that a LOW rank means the observed
outcome lies on the downward / cold side of the ensemble)
  HO1 (primary)   mean over the held-out events of the mean rank over the systems
                  that have starts 2-9 d before the event; null: 10,000 draws of
                  pseudo-onsets within +-21 calendar days of each event in any
                  hindcast year of those systems, > 135 d from every catalogued
                  onset, a pseudo-onset qualifying if at least half of the event's
                  systems have starts before it (the P' quorum rule);
                  p = P(null mean rank <= observed), one-sided.
  HO2 (secondary) the same for northern-Eurasian (NEURASIA) temperature.
  HO3 (secondary) pooled: the 17 P' events (nine confirmatory systems) together
                  with the held-out events (their systems), each event against its
                  own systems' null; p as above; for the polar cap and NEURASIA.
  Per-system and per-event ranks are reported, with the number of starts.
  Reading, fixed now: HO1 below the null mean is consistency with the registered
  result; p < 0.05 is an independent replication. With three events and about
  three systems the power is low, and a non-significant result is reported as
  uninformative, not as a failure to replicate, unless the mean rank lies above
  the null mean.

Output: results/current/6_predictability/s2s_heldout_test.json
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
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
import s2s_forecast_test as T                        # noqa: E402
import s2s_multimodel_test as MM                     # noqa: E402
import s2s_regional_test as SR                       # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

NAME = "s2s_heldout_test"
N_NULL = 10000
HELDOUT = [pd.Timestamp(d) for d in ("2023-02-16", "2024-01-16", "2024-03-04")]
SYSTEMS = ["cnrm", "ecmwf2025", "cma2025"]
WB2_PSL, ARCO_PSL = ING / "era5_psl_cap_6h.parquet", ING / "era5_psl_cap_00utc_arco.parquet"
WB2_T2M, ARCO_T2M = ING / "era5_t2m_regions_daily.parquet", ING / "era5_t2m_regions_daily_arco.parquet"
SPL_PSL, SPL_T2M = ING / "era5_psl_cap_00utc_spliced.parquet", ING / "era5_t2m_regions_daily_spliced.parquet"
OVERLAP = [("2020-11-01", "2021-04-30"), ("2021-11-01", "2022-04-30"), ("2022-11-01", "2023-01-10")]
MIN_CORR = 0.99


def in_overlap(ix):
    ix = pd.DatetimeIndex(ix)
    m = np.zeros(len(ix), bool)
    for a, b in OVERLAP:
        m |= (ix >= pd.Timestamp(a)) & (ix <= pd.Timestamp(b) + pd.Timedelta(hours=23))
    return m


def splice():
    """WB2 where it exists, ARCO minus the overlap-mean difference after it."""
    w = pd.read_parquet(WB2_PSL); w["time"] = pd.to_datetime(w["time"])
    w = w[w["time"].dt.hour == 0].set_index("time")["psl_cap_N"]
    a = pd.read_parquet(ARCO_PSL); a = a.set_index(pd.to_datetime(a["time"]))["psl_cap_N"]
    both = w.index.intersection(a.index); both = both[in_overlap(both)]
    d = a[both] - w[both]
    stats = {"psl_cap_N": {"n_overlap": int(len(both)), "offset": round(float(d.mean()), 3),
                           "rms_diff": round(float(np.sqrt((d ** 2).mean())), 3),
                           "corr": round(float(np.corrcoef(a[both], w[both])[0, 1]), 6)}}
    ext = (a[a.index > w.index.max()] - d.mean())
    psl = pd.concat([w, ext]).sort_index()
    pd.DataFrame({"time": psl.index, "psl_cap_N": psl.values}).to_parquet(SPL_PSL)
    wt = pd.read_parquet(WB2_T2M).set_index("date"); wt.index = pd.to_datetime(wt.index)
    at = pd.read_parquet(ARCO_T2M).set_index("date"); at.index = pd.to_datetime(at.index)
    both = wt.index.intersection(at.index); both = both[in_overlap(both)]
    out = wt.copy()
    ext = at[at.index > wt.index.max()].copy()
    for r in ("NEURASIA", "HI_EUROPE", "MID_EASIA", "MID_NAMER"):
        dd = at.loc[both, r] - wt.loc[both, r]
        stats[r] = {"n_overlap": int(len(both)), "offset": round(float(dd.mean()), 4),
                    "rms_diff": round(float(np.sqrt((dd ** 2).mean())), 4),
                    "corr": round(float(np.corrcoef(at.loc[both, r], wt.loc[both, r])[0, 1]), 6)}
        ext[r] = ext[r] - dd.mean()
    t2 = pd.concat([out, ext]).sort_index()
    t2.index.name = "date"
    t2.reset_index().to_parquet(SPL_T2M)
    ok = all(v["corr"] >= MIN_CORR for v in stats.values())
    return stats, ok


def files(outcome):
    tagf = {"cnrm": "cnrm", "ecmwf2025": "ecmwf2025", "cma2025": "cma2025"}
    if outcome == "psl":
        return {c: ING / f"s2s_{tagf[c]}_psl_cap.parquet" for c in SYSTEMS}
    return {c: ING / f"s2s_{tagf[c]}_t2m_regions.parquet" for c in SYSTEMS}


def setup(outcome):
    T.OBS = SPL_PSL
    if outcome == "psl":
        T.WIN, T.load = (8, 25), T.__dict__["_orig_load"]
    else:
        SR.OBS = SPL_T2M
        T.WIN, T.load = SR.WIN_T, SR.loader(outcome)


def centres(fmap):
    return {c: MM.Centre(c, f) for c, f in fmap.items() if f.exists()}


def test(cen, events, cat, rng):
    """Multi-system mean rank at each event (systems with starts), and the P'
    quorum calendar-window null. Returns a dict (None if nothing is covered)."""
    def zone_free(p):
        return bool(np.all(np.abs((cat - p).days) > T.ZONE_SEP))
    cover = {o: [c for c in cen if cen[c].at(o) is not None] for o in events}
    ev = [o for o in events if cover[o]]
    if not ev:
        return None
    years = sorted({y for C in cen.values() for y in C.hyears} | {y + 1 for C in cen.values() for y in C.hyears})
    cands = {o: [p for y in years for dd in range(-MM.NULL_WINDOW, MM.NULL_WINDOW + 1)
                 for p in [o + pd.DateOffset(years=y - o.year) + pd.Timedelta(days=dd)] if zone_free(p)]
             for o in ev}
    need = {o: int(np.ceil(len(cover[o]) / 2)) for o in ev}
    pools = {o: [p for p in cands[o] if sum(cen[c].at(p) is not None for c in cover[o]) >= need[o]] for o in ev}
    ev_pit = {o: float(np.mean([cen[c].at(o)["pit"] for c in cover[o]])) for o in ev}
    obs = float(np.mean(list(ev_pit.values())))
    null = []
    for _ in range(N_NULL):
        v = []
        for o in ev:
            p = pools[o][rng.integers(len(pools[o]))]
            v.append(np.mean([cen[c].at(p)["pit"] for c in cover[o] if cen[c].at(p) is not None]))
        null.append(np.mean(v))
    null = np.array(null)
    return {"n_events": len(ev),
            "events": {str(o.date()): {"systems": cover[o], "mean_rank": round(ev_pit[o], 4),
                                       "n_starts": {c: cen[c].at(o)["n_starts"] for c in cover[o]},
                                       "per_system_rank": {c: round(cen[c].at(o)["pit"], 4) for c in cover[o]},
                                       "A_ens": {c: round(cen[c].at(o)["A_ens"], 3) for c in cover[o]},
                                       "A_obs": round(cen[cover[o][0]].at(o)["A_obs"], 3),
                                       "n_pseudo_candidates": len(pools[o])} for o in ev},
            "mean_rank": round(obs, 4), "null_mean": round(float(null.mean()), 4),
            "null_q025_q975": [round(float(q), 4) for q in np.quantile(null, [0.025, 0.975])],
            "p": round(float(np.mean(null <= obs)), 4)}


def main():
    T.__dict__["_orig_load"] = T.load
    res = {"plan_approved": "2026-09-30", "registered_commit": "4c8b3f2", "n_null": N_NULL,
           "seed_rule": "crc32(NAME|outcome|set)", "systems": SYSTEMS,
           "heldout_events": [str(o.date()) for o in HELDOUT]}
    stats, ok = splice()
    res["splice"] = stats
    print(f"splice: {stats}", flush=True)
    if not ok:
        res["status"] = f"splice failed (overlap correlation below {MIN_CORR}); test not run"
        print(res["status"])
    else:
        cat = load_catalogue("primary")
        for outcome, lab in (("psl", "HO1_polar_cap"), ("NEURASIA", "HO2_NEURASIA")):
            setup(outcome)
            cen = centres(files(outcome))
            res[lab] = {"systems_found": list(cen),
                        **(test(cen, HELDOUT, cat, np.random.default_rng(zlib.crc32(f"{NAME}|{outcome}|heldout".encode())))
                           or {"note": "no held-out event covered"})}
            print(lab, res[lab], flush=True)
            # HO3 pooled: the P' events with the nine confirmatory systems, plus held-out
            conf = centres({c: (MM.FILES[c] if outcome == "psl" else SR.FILES[c]) for c in MM.CONFIRM})
            old_ev = [o for o in cat if o < pd.Timestamp("2022-06-01")]
            rng = np.random.default_rng(zlib.crc32(f"{NAME}|{outcome}|pooled".encode()))
            a = test(conf, old_ev, cat, rng)
            b = test(cen, HELDOUT, cat, rng)
            if a and b:
                n = a["n_events"] + b["n_events"]
                obs = (a["mean_rank"] * a["n_events"] + b["mean_rank"] * b["n_events"]) / n
                # pooled null: redraw both parts jointly
                na = test_null(conf, list(map(pd.Timestamp, a["events"])), cat, rng)
                nb = test_null(cen, list(map(pd.Timestamp, b["events"])), cat, rng)
                pooled = (na * a["n_events"] + nb * b["n_events"]) / n
                res[lab.replace("HO1", "HO3").replace("HO2", "HO3") + "_pooled"] = {
                    "n_events": n, "mean_rank": round(float(obs), 4),
                    "null_mean": round(float(pooled.mean()), 4),
                    "null_q025_q975": [round(float(q), 4) for q in np.quantile(pooled, [0.025, 0.975])],
                    "p": round(float(np.mean(pooled <= obs)), 4),
                    "note": "P' events re-evaluated with the spliced observation series, whose "
                            "leave-one-year-out observed climatology now includes 2022-2024",
                    "P_events_alone": {k: a[k] for k in ("n_events", "mean_rank", "null_mean", "p")}}
                print("HO3", outcome, res[lab.replace("HO1", "HO3").replace("HO2", "HO3") + "_pooled"], flush=True)
    (RESULTS / "s2s_heldout_test.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> s2s_heldout_test.json")


def test_null(cen, events, cat, rng):
    """Null draws of the multi-system mean rank for a fixed event list."""
    def zone_free(p):
        return bool(np.all(np.abs((cat - p).days) > T.ZONE_SEP))
    cover = {o: [c for c in cen if cen[c].at(o) is not None] for o in events}
    years = sorted({y for C in cen.values() for y in C.hyears} | {y + 1 for C in cen.values() for y in C.hyears})
    need = {o: int(np.ceil(len(cover[o]) / 2)) for o in events}
    pools = {}
    for o in events:
        c_ = [p for y in years for dd in range(-MM.NULL_WINDOW, MM.NULL_WINDOW + 1)
              for p in [o + pd.DateOffset(years=y - o.year) + pd.Timedelta(days=dd)] if zone_free(p)]
        pools[o] = [p for p in c_ if sum(cen[c].at(p) is not None for c in cover[o]) >= need[o]]
    out = np.empty(N_NULL)
    for k in range(N_NULL):
        out[k] = np.mean([np.mean([cen[c].at(p)["pit"] for c in cover[o] if cen[c].at(p) is not None])
                          for o in events for p in [pools[o][rng.integers(len(pools[o]))]]])
    return out


if __name__ == "__main__":
    main()
