#!/usr/bin/env python3
"""
s2s_postonset_test.py
=====================
IS THE UNDER-FORECAST SHIFT A MISSED SSW, OR A COUPLING DEFICIT?

Plan approved 2026-09-29. In P' (and the regional test) observed outcomes after
SSWs lie low in ensembles started 2-9 days BEFORE onset. That mixes forecasts that
did not predict the warming with forecasts that did. Written before any of the
short-lead or 10 hPa data were analysed.

DATA
  Forecasts: the P' systems and model versions; polar-cap msl (leads 1-9 from
  s2s_<tag>_psl_cap_short.parquet + 10.. from the P' files) and N-Eurasian 2 m
  temperature (windows 1-9 + 10..); zonal-mean u at 10 hPa, 60N, leads 1-15
  (s2s_<tag>_u10_60N.parquet). Observations and anomalies as in P'.

TESTS (multi-model mean of the nine confirmatory systems; outcome windows as in P'
and the regional test; calendar-window null with the same start offsets)
  T1 POST-ONSET: starts 0-7 days AFTER onset (the SSW is in the initial state).
     H1-post: mean rank lower than on event-free dates; p = P(null <= observed).
  T2 HITS vs MISSES (starts 2-9 d before onset): a start is a HIT if at least half
     of its members have u(10 hPa, 60N) < 0 on some lead within +-3 days of the
     observed onset, else a MISS. Mean rank over hit starts and over miss starts
     (event level: mean over that system's starts of the class; multi-model mean
     over systems), each against the P' null of all pre-onset starts; and the
     hit-minus-miss difference with a 95% interval from 10,000 bootstrap resamples
     of events (events having both classes).
  Reading, fixed in advance:
    T1 and T2-hit low     -> a deficit in the response GIVEN the SSW (coupling or
                             amplitude): a new result, stated as such;
    T1 and T2-hit ~ null, T2-miss low -> the deficit is failing to forecast the SSW
                             itself; the paper says so.

Output: results/current/6_predictability/s2s_postonset_test.json
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

NAME = "s2s_postonset_test"
N_NULL = 10000
N_BOOT = 10000
POST = (-7, 0)                 # k = onset - init: starts 0..7 days AFTER onset
PRE = MM.K_SHORT               # 2..9 days before
HIT_TOL = 3


def tag(c):
    return "ecmf" if c == "ecmwf" else c


def load_all(c, outcome):
    """Forecast table with leads 1.. (short + main files) and anomalies, as T.load."""
    if outcome == "psl":
        main = MM.FILES[c]; short = ING / f"s2s_{tag(c)}_psl_cap_short.parquet"; col = "psl_cap_N"
    else:
        main = SR.FILES[c]; short = ING / f"s2s_{tag(c)}_t2m_regions_short.parquet"; col = outcome
    fc = pd.concat([pd.read_parquet(main), pd.read_parquet(short)], ignore_index=True)
    fc["init"] = pd.to_datetime(fc["init"]); fc["hyear"] = fc["init"].dt.year
    fc["md"] = pd.to_datetime(fc["model_date"]).dt.strftime("%m-%d")
    g = fc.groupby(["md", "lead_day", "hyear"])[col].agg(["sum", "count"])
    tot = g.groupby(level=[0, 1]).sum()
    loo = tot.reindex(g.index.droplevel(2)).values - g.values
    clim = pd.Series(loo[:, 0] / loo[:, 1], index=g.index)
    anom = fc[col].values - clim.reindex(pd.MultiIndex.from_frame(fc[["md", "lead_day", "hyear"]])).values
    fc["anom"] = anom if outcome == "psl" else -anom          # A = -anom (psl) / +anom (T)
    return fc


def observed(outcome):
    if outcome == "psl":
        ob = pd.read_parquet(T.OBS); ob = ob[pd.to_datetime(ob["time"]).dt.hour == 0]
        return ob.set_index("time")["psl_cap_N"]
    et = pd.read_parquet(SR.OBS)
    return pd.Series(-et[outcome].values, index=pd.to_datetime(et["date"]))


class Sys:
    def __init__(self, c, outcome):
        fc = load_all(c, outcome)
        self.ob = observed(outcome)
        self.n_mem = int(fc["member"].nunique())
        self.hyears = sorted(int(y) for y in fc["hyear"].unique())
        self.min_other = max(8, int(0.75 * (len(self.hyears) - 1)))
        self.inits = sorted(pd.Timestamp(i) for i in fc["init"].unique())
        self.by_init = {pd.Timestamp(i): g for i, g in fc.groupby("init")}
        u = pd.read_parquet(ING / f"s2s_{tag(c)}_u10_60N.parquet")
        u["init"] = pd.to_datetime(u["init"])
        self.u = {pd.Timestamp(i): g.pivot(index="member", columns="lead_day", values="u10_60N")
                  for i, g in u.groupby("init")}
        self._c = {}

    def starts(self, p, kr):
        key = (pd.Timestamp(p), kr)
        if key in self._c:
            return self._c[key]
        out = []
        for i in self.inits:
            k = (key[0] - i).days
            if kr[0] <= k <= kr[1]:
                s = T.start_stats(self.by_init, self.ob, i, k, self.hyears, self.n_mem, self.min_other)
                if s is not None:
                    hit = None
                    if k >= 1 and i in self.u:
                        lead = [L for L in range(k - HIT_TOL, k + HIT_TOL + 1) if L in self.u[i].columns]
                        if lead:
                            hit = bool(((self.u[i][lead] < 0).any(axis=1)).mean() >= 0.5)
                    out.append({**s, "hit": hit})
        self._c[key] = out
        return out


def run(outcome):
    rng = np.random.default_rng(zlib.crc32(f"{NAME}|{outcome}".encode()))
    T.WIN = (8, 25) if outcome == "psl" else SR.WIN_T
    sysn = [c for c in MM.CONFIRM]
    S = {c: Sys(c, outcome) for c in sysn}
    cat = load_catalogue("primary")

    def zone_free(p):
        return bool(np.all(np.abs((cat - p).days) > T.ZONE_SEP))

    def ev_rank(p, kr, cls=None):
        vals = []
        for c in sysn:
            st = S[c].starts(p, kr)
            if cls is not None:
                st = [s for s in st if s["hit"] is cls]
            if st:
                vals.append(np.mean([s["pit"] for s in st]))
        return float(np.mean(vals)) if vals else None

    out = {}
    years = sorted({y for s in S.values() for y in s.hyears} | {y + 1 for s in S.values() for y in s.hyears})
    for label, kr in (("T1_post_onset", POST), ("pre_onset_all", PRE)):
        ev = [pd.Timestamp(o) for o in cat if ev_rank(o, kr) is not None]
        obs = float(np.mean([ev_rank(o, kr) for o in ev]))
        cands = {o: [p for y in years for dd in range(-MM.NULL_WINDOW, MM.NULL_WINDOW + 1)
                     for p in [o + pd.DateOffset(years=y - o.year) + pd.Timedelta(days=dd)]
                     if zone_free(p)] for o in ev}
        pools = {o: [p for p in cands[o] if ev_rank(p, kr) is not None] for o in ev}
        null = np.array([np.mean([ev_rank(pools[o][rng.integers(len(pools[o]))], kr) for o in ev])
                         for _ in range(N_NULL)])
        out[label] = {"n_events": len(ev), "mean_rank": round(obs, 4), "null_mean": round(float(null.mean()), 4),
                      "null_q025_q975": [round(float(q), 4) for q in np.quantile(null, [0.025, 0.975])],
                      "p": round(float(np.mean(null <= obs)), 4)}
        if label == "pre_onset_all":
            pre_null = null
            for cls, nm in ((True, "T2_hits"), (False, "T2_misses")):
                evc = [o for o in ev if ev_rank(o, kr, cls) is not None]
                v = float(np.mean([ev_rank(o, kr, cls) for o in evc])) if evc else np.nan
                out[nm] = {"n_events": len(evc), "mean_rank": round(v, 4),
                           "p_vs_all_start_null": round(float(np.mean(pre_null <= v)), 4) if evc else None,
                           "n_starts": int(sum(1 for c in sysn for o in evc for s in S[c].starts(o, kr) if s["hit"] is cls))}
            both = [o for o in ev if ev_rank(o, kr, True) is not None and ev_rank(o, kr, False) is not None]
            d = np.array([ev_rank(o, kr, True) - ev_rank(o, kr, False) for o in both])
            if len(d) > 2:
                bs = [rng.choice(d, len(d)).mean() for _ in range(N_BOOT)]
                out["T2_hit_minus_miss"] = {"n_events": len(d), "mean": round(float(d.mean()), 4),
                                            "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])]}
    return out


def main():
    res = {"plan_approved": "2026-09-29", "n_null": N_NULL, "n_boot": N_BOOT, "post_k": list(POST),
           "pre_k": list(PRE), "hit_rule": f">=50% of members u10(60N)<0 within +-{HIT_TOL} d of onset",
           "outcomes": {}}
    for oc in ("psl", "NEURASIA"):
        res["outcomes"][oc] = run(oc)
        for k, v in res["outcomes"][oc].items():
            print(oc, k, v, flush=True)
    (RESULTS / "s2s_postonset_test.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> s2s_postonset_test.json")


if __name__ == "__main__":
    main()
