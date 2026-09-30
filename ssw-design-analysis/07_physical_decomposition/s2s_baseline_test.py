#!/usr/bin/env python3
"""
s2s_baseline_test.py
====================
IS THE POST-SSW RANK DEFICIT SPECIFIC TO SSWs, OR A GENERAL VORTEX-STATE ERROR?

Plan approved 2026-09-29 (referee concern: the event-free baseline of P' is drawn
from SSW-free winters, which are enriched in strong-vortex states, so a deficit of
predictable-signal amplitude that is symmetric in vortex state would put ranks
below 0.5 after SSWs AND above 0.5 in strong-vortex winters, and 0.40 against 0.54
would count it twice). Written before these tests were run.

DATA
  The forecasts, observations, anomalies, systems, start window (2-9 d before the
  date) and event statistics of result P' (s2s_multimodel_test.Centre) and of the
  regional test (s2s_regional_test.loader, NEURASIA). Vortex state: NCEP 60-90N
  polar-cap mean zonal wind at 10 hPa on the (pseudo-)onset date
  (data/processed/atmospheric/ncep_stratosphere_1958.parquet, uwnd_ms_10hPa), a
  monotone index of vortex strength.

TESTS (multi-model mean of the nine confirmatory systems, as in P'; each outcome:
polar-cap NAM proxy, and NEURASIA 2 m temperature)
  B0  mean rank after SSWs against 0.5 (calibration in absolute terms): 95%
      interval from 10,000 bootstrap resamples of events.
  B1  ALL-WINTER BASELINE: pseudo-onsets from every December-March date of the
      hindcast years more than 30 days from every catalogued onset (SSW winters
      included, SSW windows excluded), within +-21 calendar days of each event;
      10,000 draws; p = P(null mean rank <= observed).
  B2  VORTEX-STATE REGRESSION: over all dates of B1, rank = a + b u10 (OLS on the
      multi-model date values). Residual of the SSWs = mean over events of
      rank - (a + b u10 at onset); null: the same residual on B1 draws; p =
      P(null <= observed). If B2 is not significant while B1 is, the deficit is a
      general error of amplitude with vortex state, reached at SSWs because they
      are the weakest vortex states; if B2 is significant, it is specific to SSWs.
  Both answers are reported; the manuscript's wording follows B1 and B2.

Output: results/current/6_predictability/s2s_baseline_test.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
import s2s_forecast_test as T                        # noqa: E402
import s2s_multimodel_test as MM                     # noqa: E402
import s2s_regional_test as SR                       # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

NAME = "s2s_baseline_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_NULL = 10000
N_BOOT = 10000
SEP_ALLWINTER = 30
U10 = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere_1958.parquet"


def run(outcome):
    """outcome: 'psl' (P' setup) or 'NEURASIA' (regional setup)."""
    rng = np.random.default_rng(zlib.crc32(f"{NAME}|{outcome}".encode()))
    if outcome == "psl":
        T.WIN, T.load = (8, 25), T.__dict__["_orig_load"]
        files = {c: f for c, f in MM.FILES.items() if c in MM.CONFIRM}
    else:
        T.WIN, T.load = SR.WIN_T, SR.loader(outcome)
        files = {c: f for c, f in SR.FILES.items() if c in MM.CONFIRM}
    cen = {c: MM.Centre(c, f) for c, f in files.items() if f.exists()}
    cat = load_catalogue("primary")
    u = pd.read_parquet(U10)["uwnd_ms_10hPa"]
    u.index = pd.to_datetime(u.index).tz_localize(None).normalize()

    years = sorted({y for C in cen.values() for y in C.hyears} | {y + 1 for C in cen.values() for y in C.hyears})
    days = [d for y in years for d in pd.date_range(f"{y - 1}-12-01", f"{y}-03-31")]
    allw = [d for d in days if np.all(np.abs((cat - d).days) > SEP_ALLWINTER)]

    def mm_at(p):
        vals = [C.at(p) for C in cen.values()]
        vals = [v for v in vals if v is not None]
        return float(np.mean([v["pit"] for v in vals])) if vals else None

    ev = []
    for o in cat:
        if any(C.at(o) is not None for C in cen.values()):
            ev.append(pd.Timestamp(o))
    rank_ev = np.array([mm_at(o) for o in ev])
    u_ev = np.array([u.get(o, np.nan) for o in ev])
    pool = {}
    for d in allw:
        r = mm_at(d)
        if r is not None and np.isfinite(u.get(d, np.nan)):
            pool[d] = (r, float(u[d]))
    pdates = np.array(list(pool))
    prank = np.array([pool[d][0] for d in pdates]); pu = np.array([pool[d][1] for d in pdates])
    b, a = np.polyfit(pu, prank, 1)
    # B0
    bs = np.array([rng.choice(rank_ev, len(rank_ev)).mean() for _ in range(N_BOOT)])
    # candidates per event within +-21 calendar days
    doy = np.array([d.dayofyear for d in pdates])
    cand = {}
    for o in ev:
        dd = np.abs(doy - o.dayofyear); dd = np.minimum(dd, 365 - dd)
        cand[o] = np.where(dd <= T.NULL_WINDOW)[0]
    null_rank, null_res = np.empty(N_NULL), np.empty(N_NULL)
    for k in range(N_NULL):
        idx = np.array([rng.choice(cand[o]) for o in ev])
        null_rank[k] = prank[idx].mean()
        null_res[k] = (prank[idx] - (a + b * pu[idx])).mean()
    ok = np.isfinite(u_ev)
    res_ev = float(np.mean(rank_ev[ok] - (a + b * u_ev[ok])))
    return {
        "seed": int(zlib.crc32(f"{NAME}|{outcome}".encode())),
        "n_events": int(len(ev)), "n_allwinter_dates": int(len(pdates)),
        "systems": list(cen),
        "B0_mean_rank": round(float(rank_ev.mean()), 4),
        "B0_ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
        "B1_null_mean": round(float(null_rank.mean()), 4),
        "B1_null_q025_q975": [round(float(q), 4) for q in np.quantile(null_rank, [0.025, 0.975])],
        "B1_p": round(float(np.mean(null_rank <= rank_ev.mean())), 4),
        "B2_slope_per_ms": round(float(b), 5), "B2_intercept": round(float(a), 4),
        "B2_u10_mean_events": round(float(np.nanmean(u_ev)), 2),
        "B2_u10_range_allwinter": [round(float(pu.min()), 1), round(float(pu.max()), 1)],
        "B2_residual_events": round(res_ev, 4),
        "B2_null_q025_q975": [round(float(q), 4) for q in np.quantile(null_res, [0.025, 0.975])],
        "B2_p": round(float(np.mean(null_res <= res_ev)), 4),
        "B2_rank_predicted_at_event_u10": round(float(np.mean(a + b * u_ev[ok])), 4),
    }


def main():
    T.__dict__["_orig_load"] = T.load
    res = {"plan_approved": "2026-09-29", "n_null": N_NULL, "n_boot": N_BOOT,
           "seed_rule": "crc32(NAME|outcome), recorded per outcome",
           "sep_allwinter_days": SEP_ALLWINTER, "vortex_index": "NCEP 60-90N mean u, 10 hPa",
           "inputs": [str(U10.relative_to(ROOT))], "outcomes": {}}
    for oc in ("psl", "NEURASIA"):
        r = run(oc)
        res["outcomes"][oc] = r
        print(f"{oc:9s} n={r['n_events']} rank {r['B0_mean_rank']} {r['B0_ci95']} | all-winter null "
              f"{r['B1_null_mean']} p={r['B1_p']} | slope {r['B2_slope_per_ms']}/(m/s) predicted at "
              f"event u10 {r['B2_rank_predicted_at_event_u10']} residual {r['B2_residual_events']} "
              f"p={r['B2_p']}", flush=True)
    (RESULTS / "s2s_baseline_test.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> s2s_baseline_test.json")


if __name__ == "__main__":
    main()
