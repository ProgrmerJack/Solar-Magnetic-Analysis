#!/usr/bin/env python3
"""
s2s_multimodel_test.py
======================
MULTI-MODEL S2S TEST: EVENT DISCRIMINATION (D) AND THE SIZE OF THE SHIFT (H1, H2)

Plan approved 2026-09-25 (item A). Hypotheses H1 and H2 were written into this
docstring before any non-ECMWF forecast was retrieved; ECMWF (result P) is the
DISCOVERY set, every other centre is the CONFIRMATORY set.

DATA
  s2s_<origin>_psl_cap.parquet (acquire_s2s_reforecasts.py <origin>): ECCC, CMA,
  HMCR, KMA, CNRM, JMA, CNR-ISAC, NCEP, CPTEC; one model version each; Dec-Mar
  hindcast starts; leads 10-34 d; polar-cap (60-90N) msl at 00 UTC reduced as ERA5.
  ERA5 00 UTC polar cap as the observation. Events: load_catalogue("primary").

UNITS AND ANOMALIES (as result P, s2s_forecast_test.py, whose functions are used)
  forecast anomaly = leave-one-year-out mean of the centre's other hindcast years at
  the same start key and lead; observed anomaly = ERA5 minus its leave-one-year-out
  mean over the same hindcast years at the same calendar dates (at least
  max(8, 0.75 x (n-1)) other years). A = -anomaly over days +8..+25 after onset.
  SHORT LEAD ONLY: starts 2-9 d before onset. An event's value for a centre is the
  mean over that centre's starts in that window.

TESTS (per centre; multi-model = mean over centres covering the event)
  D   discrimination: r(ensemble-mean A, observed A) across events; conditional
      calendar-window null as in P (pseudo-onsets within +-21 calendar days of the
      event date, > 135 d from every catalogued onset; draws whose across-event
      variance of observed A is within a factor 1.25 of the events').
  H1  (confirmatory) observed outcomes after SSWs lie on the DOWNWARD side of the
      forecast ensembles: mean PIT over events is LOWER than on pseudo-events.
      p = P(null mean PIT <= observed).
  H2  (confirmatory) correcting only the size of the shift -- adding to every
      member the leave-one-event-out mean of (observed - ensemble-mean) A --
      lowers the ensemble CRPS at SSWs more than the same correction does on
      pseudo-events. p = P(null gain >= observed gain).
  Primary: the multi-model mean over the CONFIRMATORY centres. Per-centre and
  "all centres incl. ECMWF" are reported as secondary.

Output: results/current/6_predictability/s2s_multimodel_test.json
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
from build_catalogue import load_catalogue          # noqa: E402

CONFIRM = ["eccc", "cma", "hmcr", "kma", "cnrm", "jma", "cnr_isac", "ncep", "cptec"]
DISCOVERY = "ecmwf"
FILES = {"ecmwf": ING / "s2s_ecmf_psl_cap.parquet",
         **{c: ING / f"s2s_{c}_psl_cap.parquet" for c in CONFIRM}}
K_SHORT = (2, 9)
N_NULL = 10000
NULL_WINDOW = 21
NAME = "s2s_multimodel_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)


def crps_ens(x, y):
    x = np.asarray(x, float)
    return float(np.mean(np.abs(x - y)) - 0.5 * np.mean(np.abs(x[:, None] - x[None, :])))


class Centre:
    def __init__(self, name, path):
        T.FC = path
        fc, ob = T.load()
        self.name, self.ob = name, ob
        self.n_mem = int(fc["member"].nunique())
        self.hyears = sorted(int(y) for y in fc["hyear"].unique())
        self.min_other = max(8, int(0.75 * (len(self.hyears) - 1)))
        self.inits = sorted(pd.Timestamp(i) for i in fc["init"].unique())
        self.by_init = {pd.Timestamp(i): g for i, g in fc.groupby("init")}
        self._cache = {}

    def at(self, p):
        """Event-level statistics for a (pseudo-)onset p, or None."""
        key = pd.Timestamp(p)
        if key in self._cache:
            return self._cache[key]
        ss = []
        for i in self.inits:
            k = (key - i).days
            if K_SHORT[0] <= k <= K_SHORT[1]:
                s = T.start_stats(self.by_init, self.ob, i, k, self.hyears,
                                  self.n_mem, self.min_other)
                if s is not None:
                    ss.append(s)
        out = None
        if ss:
            out = {"A_ens": float(np.mean([s["A_ens"] for s in ss])),
                   "A_obs": ss[0]["A_obs"], "pit": float(np.mean([s["pit"] for s in ss])),
                   "members": [s["A_members"] for s in ss], "n_starts": len(ss)}
        self._cache[key] = out
        return out


def h2_gain(rows):
    """Total CRPS gain from adding the leave-one-event-out mean (obs - ens) shift."""
    d = np.array([r["A_obs"] - r["A_ens"] for r in rows])
    g = 0.0
    for e, r in enumerate(rows):
        delta = (d.sum() - d[e]) / (len(d) - 1)
        raw = np.mean([crps_ens(m, r["A_obs"]) for m in r["members"]])
        cor = np.mean([crps_ens(np.asarray(m) + delta, r["A_obs"]) for m in r["members"]])
        g += raw - cor
    return float(g)


def summarise(rows):
    a = np.array([r["A_ens"] for r in rows]); o = np.array([r["A_obs"] for r in rows])
    return {"r": float(np.corrcoef(a, o)[0, 1]) if len(rows) > 2 else np.nan,
            "var_obs": float(np.var(o)), "pit": float(np.mean([r["pit"] for r in rows])),
            "gain": h2_gain(rows) if len(rows) > 2 else np.nan,
            "A_ens": float(a.mean()), "A_obs": float(o.mean())}


def mm_rows(per):
    """Combine centres at event level: means of A_ens, A_obs and PIT; members pooled
    per centre so each centre's CRPS enters with equal weight."""
    return {"A_ens": float(np.mean([p["A_ens"] for p in per])),
            "A_obs": float(np.mean([p["A_obs"] for p in per])),
            "pit": float(np.mean([p["pit"] for p in per])),
            "members": [m for p in per for m in p["members"]]}


def main():
    rng = np.random.default_rng(SEED)
    cat = load_catalogue("primary")
    cen = {c: Centre(c, f) for c, f in FILES.items() if f.exists()}
    missing = [c for c in FILES if c not in cen]
    print(f"centres: {list(cen)}; missing: {missing}", flush=True)

    def zone_free(p):
        return bool(np.all(np.abs((cat - p).days) > T.ZONE_SEP))

    res = {"plan_approved": "2026-09-25", "k_short": list(K_SHORT), "n_null": N_NULL,
           "seed": SEED, "null_window_days": NULL_WINDOW, "confirmatory": CONFIRM,
           "discovery": DISCOVERY, "missing_centres": missing, "centres": {}, "multimodel": {}}

    # real events per centre
    real = {c: {} for c in cen}
    for c, C in cen.items():
        for o in cat:
            r = C.at(o)
            if r is not None:
                real[c][o] = r
    events = sorted({o for c in cen for o in real[c]})
    # candidate pseudo-onsets per event: calendar window, zone-free, any year
    years = sorted({y for C in cen.values() for y in C.hyears} | {y + 1 for C in cen.values() for y in C.hyears})
    cands = {}
    for o in events:
        lst = []
        for y in years:
            base = o + pd.DateOffset(years=y - o.year)
            for dd in range(-NULL_WINDOW, NULL_WINDOW + 1):
                p = base + pd.Timedelta(days=dd)
                if zone_free(p):
                    lst.append(p)
        cands[o] = lst

    def per_centre(c, C):
        ev = [o for o in events if o in real[c]]
        rows = [real[c][o] for o in ev]
        obs = summarise(rows)
        pools = {o: [p for p in cands[o] if C.at(p) is not None] for o in ev}
        nr, npit, ngain, nvar = [], [], [], []
        if all(pools[o] for o in ev) and len(ev) > 2:
            for _ in range(N_NULL):
                rr = [C.at(pools[o][rng.integers(len(pools[o]))]) for o in ev]
                s = summarise(rr); nr.append(s["r"]); npit.append(s["pit"])
                ngain.append(s["gain"]); nvar.append(s["var_obs"])
        return ev, obs, np.array(nr), np.array(npit), np.array(ngain), np.array(nvar)

    def report(label, ev, obs, nr, npit, ngain, nvar):
        if len(nr) == 0:
            return {"n_events": len(ev), "note": "null not computable"}
        mv = np.abs(np.log(nvar / obs["var_obs"])) < np.log(1.25)
        out = {"n_events": len(ev), "events": [str(o.date()) for o in ev],
               "D_r": round(obs["r"], 4), "D_r_null_mean": round(float(np.nanmean(nr)), 4),
               "D_p_conditional": round(float(np.nanmean(nr[mv] >= obs["r"])), 4) if mv.sum() else None,
               "D_n_conditional": int(mv.sum()),
               "H1_mean_pit": round(obs["pit"], 4), "H1_null_mean_pit": round(float(np.mean(npit)), 4),
               "H1_p": round(float(np.mean(npit <= obs["pit"])), 4),
               "H2_crps_gain": round(obs["gain"], 3), "H2_null_mean_gain": round(float(np.nanmean(ngain)), 3),
               "H2_p": round(float(np.nanmean(ngain >= obs["gain"])), 4),
               "A_ens_mean": round(obs["A_ens"], 1), "A_obs_mean": round(obs["A_obs"], 1)}
        print(f"{label:22s} n={len(ev):2d}  D r={out['D_r']:+.2f} (null {out['D_r_null_mean']:+.2f}) "
              f"p_cond={out['D_p_conditional']} | H1 PIT {out['H1_mean_pit']:.3f} "
              f"(null {out['H1_null_mean_pit']:.3f}) p={out['H1_p']:.4f} | H2 gain {out['H2_crps_gain']:+.1f} "
              f"(null {out['H2_null_mean_gain']:+.1f}) p={out['H2_p']:.4f} | A ens/obs "
              f"{out['A_ens_mean']:+.0f}/{out['A_obs_mean']:+.0f} Pa", flush=True)
        return out

    for c, C in cen.items():
        ev, obs, nr, npit, ngain, nvar = per_centre(c, C)
        res["centres"][c] = {"n_members": C.n_mem, "hindcast_years": [C.hyears[0], C.hyears[-1]],
                             **report(c, ev, obs, nr, npit, ngain, nvar)}

    # multi-model: event-level mean over the centres covering the event; each null
    # draw moves every event to one pseudo-onset valid for ALL of those centres
    for label, group in (("confirmatory", [c for c in CONFIRM if c in cen]),
                         ("all_incl_ecmwf", list(cen))):
        cover = {o: [c for c in group if o in real[c]] for o in events}
        ev = [o for o in events if cover[o]]
        rows = [mm_rows([real[c][o] for c in cover[o]]) for o in ev]
        obs = summarise(rows)
        pools = {o: [p for p in cands[o] if all(cen[c].at(p) is not None for c in cover[o])]
                 for o in ev}
        nr, npit, ngain, nvar = [], [], [], []
        if ev and all(pools[o] for o in ev):
            for _ in range(N_NULL):
                rr = []
                for o in ev:
                    p = pools[o][rng.integers(len(pools[o]))]
                    rr.append(mm_rows([cen[c].at(p) for c in cover[o]]))
                s = summarise(rr); nr.append(s["r"]); npit.append(s["pit"])
                ngain.append(s["gain"]); nvar.append(s["var_obs"])
        out = report(f"MULTI-MODEL {label}", ev, obs, np.array(nr), np.array(npit),
                     np.array(ngain), np.array(nvar))
        out["centres_per_event"] = {str(o.date()): cover[o] for o in ev}
        out["null_candidates_per_event"] = {str(o.date()): len(pools[o]) for o in ev}
        res["multimodel"][label] = out

    (RESULTS / "s2s_multimodel_test.json").write_text(json.dumps(res, indent=2, default=str),
                                                      encoding="utf8", newline="\n")
    print("Saved -> s2s_multimodel_test.json")


if __name__ == "__main__":
    main()
