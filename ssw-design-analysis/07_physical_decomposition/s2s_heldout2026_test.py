#!/usr/bin/env python3
"""
s2s_heldout2026_test.py
=======================
A FOURTH OUT-OF-SAMPLE SSW: 4 MARCH 2026, SCORED WITH REAL-TIME FORECASTS

Plan approved 2026-09-30. This design was committed BEFORE any forecast for 2026,
any 2026 reforecast of ECCC/HMCR or model year 2026 of ECMWF, and any surface
observation (sea-level pressure, 2 m temperature) after 30 April 2024 was
retrieved or examined.

WHAT HAD BEEN SEEN
  Only the stratospheric wind: the project's Charlton-Polvani detector
  (ensemble_precursor.detect_ssw, validated on the NCEP compendium) run on NCEP
  u(10 hPa, 60N) daily means to 17 March 2026 (where the NCEP series ends),
  continued with ERA5 (ARCO) 00 UTC values (correlation 0.976 with NCEP over the
  overlap), finds one onset in winter 2025-26: 4 March 2026 (easterly 4-10 March,
  westerly 11-23 March, final warming about 10 April). None in winter 2024-25. The
  onset is marginal (a late, short reversal), which is stated. The earlier held-out
  test (s2s_heldout_test.py) had failed to replicate; this design was written after
  that result and its exploratory diagnosis were known.

DISCLOSURE ADDED 2026-10-01 (after registration, before any test was run): while
  checking the date range of data/processed/atmospheric/ao_daily_cpc.txt, its last
  two lines were displayed: the CPC daily AO for 30 and 31 March 2026 (+2.53 and
  +1.96), i.e. days 26-27 after the 4 March onset, inside the outcome window. No
  other surface observation for 2026 has been examined. The design is unchanged.

FORECASTS (ECDS dataset s2s-forecasts, real time; starts 2-9 days before onset,
i.e. 23 February - 2 March 2026, each start used only if a reforecast of the same
system and model version exists for the same start month-day)
  PRIMARY systems, whose 2026 real-time model is the one of available reforecasts:
    ECMWF (reforecasts model year 2026, on the fly, 2006-2025), ECCC (model year
    2026), HMCR (model year 2026), KMA (model year 2026), NCEP (CFSv2, the only
    version, fixed hindcasts 1999-2010).
  SECONDARY (version match not verifiable): CMA (latest reforecasts model year
    2025), CNRM (2025), JMA (2022), CNR-ISAC (2023).
  Real-time ensembles are larger than reforecast ensembles; the rank of the observed
  outcome is (below + 0.5 equal + 0.5)/(n + 1) with the real-time members, whose
  expectation under calibration is 0.5 whatever n; this is stated as a difference
  from the null draws, which use reforecast members.

ANOMALIES AND OUTCOMES (as P' and the regional test)
  Forecast anomaly: real-time value minus the mean of the reforecasts of the same
  system, model version, start month-day and lead over all their hindcast years.
  Observed anomaly: ERA5 minus its mean over those hindcast years at the same
  calendar dates (WeatherBench 2 series spliced with ARCO-ERA5 exactly as in
  s2s_heldout_test.py, extended to winters 2024-25 and 2025-26 with the same
  reduction; overlap rule unchanged). Outcomes: polar-cap NAM proxy days +8..+25
  (primary) and northern-Eurasian 2 m temperature days +8..+24 (secondary); sign
  so that a low rank means downward / colder than forecast.

TESTS
  HO2026 (primary): multi-system mean rank of the primary systems at 4 March 2026
    against the calendar-window null of their reforecasts (pseudo-onsets within
    +-21 days in any hindcast year, > 135 days from every catalogued onset, quorum
    of half), 10,000 draws; one-sided p = P(null <= observed).
  Secondary: the same for temperature; with the secondary systems added; and the
    four held-out events (2023, Jan 2024, Mar 2024, 2026) together, each against its
    own systems' null.
  Reading, fixed now: a rank below the null mean is consistent with the 1998-2021
  under-forecast, above it inconsistent; with one event neither is decisive, and
  the result is reported as one more event, not as a replication or its failure.

IMPLEMENTATION NOTES (2026-10-01, written from the data inventory -- start dates,
  member counts, model years -- before any 2026 forecast value or outcome was read)
  Real-time starts 23 Feb - 2 Mar 2026 with a reforecast of the same month-day:
    ECMWF 02-23/25/27, 03-01 (model year 2026); ECCC 02-23/26, 03-02 (2026);
    HMCR 02-26 (2026); KMA 02-25, 03-01 (2026); NCEP 02-25/28, 03-01; JMA 02-25.
    CMA (real time 02-23/26, 03-02; model year 2025 reforecasts 02-20/24/27, 03-03),
    CNRM (real time 02-26; reforecasts 02-25, 03-01) and CNR-ISAC (02-26; 02-20/25,
    03-02) have none, and by the rule above are not used.
  Observations: one series per variable, WeatherBench 2 where it exists, then
    ARCO-ERA5 (2020-24 and the --late file for 2024-26), each ARCO part minus the
    single overlap-mean offset of s2s_heldout_test.splice (same overlap, same
    0.99 correlation rule); the 2023-24 values are therefore those of that test.
  Null: the reforecast files through MM.Centre (leave-one-year-out anomalies), the
    P' calendar-window and quorum rules; 4 March 2026 is added to the onsets that
    pseudo-onsets must avoid.
  Pooled four-event test: the three 2023-24 events with the systems of
    s2s_heldout_test (cnrm, ecmwf2025, cma2025) and 2026 with the primary systems;
    pooled mean rank against the jointly drawn null.

Output: results/current/6_predictability/s2s_heldout2026_test.json
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
import s2s_heldout_test as H                         # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

NAME = "s2s_heldout2026_test"
N_NULL = 10000
EVENT = pd.Timestamp("2026-03-04")
PRIMARY = ["ecmwf", "eccc", "hmcr", "kma", "ncep"]
SECONDARY = ["cma", "cnrm", "jma", "cnr_isac"]
# reforecast file tag per system: the model version of the 2026 real-time model
REF = {"ecmwf": "ecmwf2026", "eccc": "eccc2026", "hmcr": "hmcr2026", "kma": "kma",
       "ncep": "ncep", "cma": "cma2025", "cnrm": "cnrm", "jma": "jma", "cnr_isac": "cnr_isac"}
LATE_PSL = ING / "era5_psl_cap_00utc_arco_late.parquet"
LATE_T2M = ING / "era5_t2m_regions_daily_arco_late.parquet"
SPL_PSL = ING / "era5_psl_cap_00utc_spliced_late.parquet"
SPL_T2M = ING / "era5_t2m_regions_daily_spliced_late.parquet"


def splice():
    """H.splice (WB2 + ARCO 2020-24, overlap offset removed), then the 2024-26
    ARCO values with the same offsets."""
    stats, ok = H.splice()
    p = pd.read_parquet(H.SPL_PSL); p["time"] = pd.to_datetime(p["time"])
    lp = pd.read_parquet(LATE_PSL); lp["time"] = pd.to_datetime(lp["time"])
    lp["psl_cap_N"] -= stats["psl_cap_N"]["offset"]
    p = pd.concat([p, lp[lp["time"] > p["time"].max()]]).sort_values("time")
    if not p["time"].is_unique:
        raise ValueError("repeated times in the spliced polar cap")
    p.reset_index(drop=True).to_parquet(SPL_PSL)
    t = pd.read_parquet(H.SPL_T2M); t["date"] = pd.to_datetime(t["date"])
    lt = pd.read_parquet(LATE_T2M); lt["date"] = pd.to_datetime(lt["date"])
    for r in SR.REGIONS:
        lt[r] -= stats[r]["offset"]
    t = pd.concat([t, lt[lt["date"] > t["date"].max()]]).sort_values("date")
    if not t["date"].is_unique:
        raise ValueError("repeated dates in the spliced temperature")
    t.reset_index(drop=True).to_parquet(SPL_T2M)
    return stats, ok


def ref_file(c, outcome):
    v = "psl_cap" if outcome == "psl" else "t2m_regions"
    tg = "ecmf" if REF[c] == "ecmwf" else REF[c]
    return ING / f"s2s_{tg}_{v}.parquet"


def rt_file(c, outcome):
    return ING / f"s2s_rt_{c}_{'psl_cap' if outcome == 'psl' else 't2m_regions'}.parquet"


def setup(outcome):
    T.OBS = SPL_PSL
    if outcome == "psl":
        T.WIN, T.load = (8, 25), T.__dict__["_orig_load"]
    else:
        SR.OBS = SPL_T2M
        T.WIN, T.load = SR.WIN_T, SR.loader(outcome)


class RealTime:
    """2026 real-time ensemble of one system scored like MM.Centre.at: anomaly =
    real time minus the mean of the reforecasts (same file as the null) over all
    their hindcast years at the same start month-day and lead; observed anomaly
    against the same hindcast years (T.obs_anom); rank with the real-time members."""
    def __init__(self, c, outcome, ref):
        col = "psl_cap_N" if outcome == "psl" else outcome
        rf = pd.read_parquet(ref_file(c, outcome))
        rf["md"] = pd.to_datetime(rf["model_date"]).dt.strftime("%m-%d")
        clim = rf.groupby(["md", "lead_day"])[col].mean()
        rt = pd.read_parquet(rt_file(c, outcome))
        rt["init"] = pd.to_datetime(rt["init"])
        rt["md"] = rt["init"].dt.strftime("%m-%d")
        a = rt[col].values - clim.reindex(pd.MultiIndex.from_frame(rt[["md", "lead_day"]])).values
        # T.start_stats takes A = -anom; for temperature SR.loader stores -anom so
        # that A is the temperature anomaly: same convention here
        rt["anom"] = a if outcome == "psl" else -a
        self.ob, self.hyears, self.min_other = ref.ob, ref.hyears, ref.min_other
        self.by_init = {pd.Timestamp(i): g for i, g in rt.groupby("init")}
        self.n_mem = {i: int(g["member"].nunique()) for i, g in self.by_init.items()}
        self.matched = sorted(str(i.date()) for i, g in self.by_init.items() if g["anom"].notna().all())

    def at(self, o):
        ss = []
        for i in sorted(self.by_init):
            k = (o - i).days
            if MM.K_SHORT[0] <= k <= MM.K_SHORT[1]:
                s = T.start_stats(self.by_init, self.ob, i, k, self.hyears, self.n_mem[i], self.min_other)
                if s is not None:
                    ss.append(s)
        if not ss:
            return None
        return {"A_ens": float(np.mean([s["A_ens"] for s in ss])), "A_obs": ss[0]["A_obs"],
                "pit": float(np.mean([s["pit"] for s in ss])), "n_starts": len(ss),
                "n_members": [int(len(s["A_members"])) for s in ss]}


def event_and_null(systems, outcome, cat, rng):
    """Mean rank at 4 March 2026 over the systems with matched starts, and the
    calendar-window quorum null drawn from their reforecasts."""
    ref = {c: MM.Centre(c, ref_file(c, outcome)) for c in systems if ref_file(c, outcome).exists()
           and rt_file(c, outcome).exists()}
    rts = {c: RealTime(c, outcome, ref[c]) for c in ref}
    at = {c: rts[c].at(EVENT) for c in rts}
    cover = [c for c in rts if at[c] is not None]
    out = {"systems_with_files": list(ref), "matched_starts": {c: rts[c].matched for c in rts},
           "systems_scored": cover}
    if not cover:
        return out, None
    def zone_free(p):
        return bool(np.all(np.abs((cat - p).days) > T.ZONE_SEP))
    years = sorted({y for c in cover for y in ref[c].hyears} | {y + 1 for c in cover for y in ref[c].hyears})
    cands = [p for y in years for dd in range(-MM.NULL_WINDOW, MM.NULL_WINDOW + 1)
             for p in [EVENT + pd.DateOffset(years=y - EVENT.year) + pd.Timedelta(days=dd)] if zone_free(p)]
    need = int(np.ceil(len(cover) / 2))
    pool = [p for p in cands if sum(ref[c].at(p) is not None for c in cover) >= need]
    obs = float(np.mean([at[c]["pit"] for c in cover]))
    null = np.array([np.mean([ref[c].at(p)["pit"] for c in cover if ref[c].at(p) is not None])
                     for p in (pool[rng.integers(len(pool))] for _ in range(N_NULL))])
    out.update({"mean_rank": round(obs, 4), "null_mean": round(float(null.mean()), 4),
                "null_q025_q975": [round(float(q), 4) for q in np.quantile(null, [0.025, 0.975])],
                "p": round(float(np.mean(null <= obs)), 4), "n_pseudo_candidates": len(pool),
                "per_system": {c: {"rank": round(at[c]["pit"], 4), "n_starts": at[c]["n_starts"],
                                   "n_members": at[c]["n_members"], "A_ens": round(at[c]["A_ens"], 3)}
                               for c in cover},
                "A_obs": round(at[cover[0]]["A_obs"], 3)})
    return out, null


def main():
    T.__dict__["_orig_load"] = T.load
    res = {"plan_approved": "2026-09-30", "registered_commit": "a511fbc",
           "disclosure_commit": "a04a4a8", "n_null": N_NULL, "seed_rule": "crc32(NAME|outcome|set)",
           "event": str(EVENT.date()), "primary": PRIMARY, "secondary": SECONDARY, "reforecasts": REF}
    stats, ok = splice()
    res["splice"] = stats
    if not ok:
        res["status"] = "splice failed (overlap correlation below 0.99); test not run"
    else:
        cat = load_catalogue("primary").append(pd.DatetimeIndex([EVENT]))
        for outcome, lab in (("psl", "polar_cap"), ("NEURASIA", "NEURASIA")):
            setup(outcome)
            H.setup(outcome)
            T.OBS, SR.OBS = SPL_PSL, SPL_T2M          # H.setup points at its own splice
            sd = lambda s: np.random.default_rng(zlib.crc32(f"{NAME}|{outcome}|{s}".encode()))
            r, null = event_and_null(PRIMARY, outcome, cat, sd("primary"))
            res[("HO2026_" if outcome == "psl" else "HO2026_secondary_") + lab] = r
            print(lab, "primary", json.dumps(r), flush=True)
            r2, _ = event_and_null(PRIMARY + SECONDARY, outcome, cat, sd("with_secondary"))
            res[f"with_secondary_{lab}"] = r2
            print(lab, "with secondary", json.dumps(r2), flush=True)
            if null is None:
                continue
            # four held-out events: 2023-24 with the s2s_heldout_test systems
            H.N_NULL = N_NULL
            cen = H.centres(H.files(outcome))
            rng = sd("four_events")
            a = H.test(cen, H.HELDOUT, cat, rng)
            if a:
                na = H.test_null(cen, list(map(pd.Timestamp, a["events"])), cat, rng)
                n = a["n_events"] + 1
                obs = (a["mean_rank"] * a["n_events"] + r["mean_rank"]) / n
                pooled = (na * a["n_events"] + null) / n
                res[f"four_events_{lab}"] = {
                    "n_events": n, "mean_rank": round(float(obs), 4),
                    "null_mean": round(float(pooled.mean()), 4),
                    "null_q025_q975": [round(float(q), 4) for q in np.quantile(pooled, [0.025, 0.975])],
                    "p": round(float(np.mean(pooled <= obs)), 4),
                    "events_2023_24": {k: v["mean_rank"] for k, v in a["events"].items()},
                    "event_2026": r["mean_rank"]}
                print(lab, "four events", json.dumps(res[f"four_events_{lab}"]), flush=True)
    (RESULTS / "s2s_heldout2026_test.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> s2s_heldout2026_test.json")


if __name__ == "__main__":
    main()
