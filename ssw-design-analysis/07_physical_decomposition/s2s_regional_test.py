#!/usr/bin/env python3
"""
s2s_regional_test.py
====================
DO OPERATIONAL FORECASTS UNDER-FORECAST REGIONAL COLD AFTER SSWs?

Plan approved 2026-09-26 (regional cold-risk plan). Written before any regional
forecast or observed temperature was analysed. The operational counterpart of
snapsi_regional_test.py, and the regional counterpart of result P' (the
multi-model polar-cap test, s2s_multimodel_test.py), whose design, null and
known-truth calibration it reuses unchanged.

DATA
  s2s_<origin>_t2m_regions.parquet (acquire_s2s_reforecasts.py --var t2m): daily
  mean 2 m temperature over four regions, the same ten systems and model
  versions as P', starts 2-9 d before onset. Observation:
  era5_t2m_regions_daily.parquet (acquire_era5_t2m_regions.py), same regions.
  Regions: NEURASIA 50-65N 10-130E (PRIMARY; Kretschmer et al. 2018), HI_EUROPE
  55-70N 0-60E, MID_EASIA 35-55N 90-150E, MID_NAMER 35-55N 120-60W (Huang et al.
  2021).

OUTCOME
  A_T = regional temperature anomaly averaged over the daily means of post-onset
  days +8..+24, anomalies as in P' (forecast: leave-one-year-out mean of the other
  hindcast years at the same start key and lead; observed: the other years at the
  same calendar dates). Negative = colder than normal.

TESTS (s2s_multimodel_test machinery; statistic = mean over the confirmatory
systems covering each event; 10,000 calendar-window pseudo-onset draws)
  H1-T (primary, NEURASIA)  observed outcomes after SSWs lie on the COLD side of the
        forecast ensembles more than on event-free dates: mean PIT lower than the
        null; p = P(null mean PIT <= observed). Falsifier: p >= 0.05 -> no evidence
        that the polar-cap under-forecast (P') reaches regional temperature.
  Secondary: H1-T in the other three regions (reported without correction, and
        with Holm's correction over the four); D (discrimination) and H2 (shift
        correction) as in P'; the system-centred and leave-one-event-out variants.

Output: results/current/6_predictability/s2s_regional_test.json
"""
import contextlib
import io
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
import s2s_forecast_test as T                        # noqa: E402
import s2s_multimodel_test as MM                     # noqa: E402

NAME = "s2s_regional_test"
REGIONS = ["NEURASIA", "HI_EUROPE", "MID_EASIA", "MID_NAMER"]
PRIMARY = "NEURASIA"
FILES = {c: ING / (f"s2s_ecmf_t2m_regions.parquet" if c == "ecmwf"
                   else f"s2s_{c}_t2m_regions.parquet") for c in MM.FILES}
OBS = ING / "era5_t2m_regions_daily.parquet"
WIN_T = (8, 24)


def loader(region):
    """Replacement for T.load: (forecast with 'anom', observed series), both with
    the sign REVERSED, because T.start_stats takes A = -anomaly (the NAM sign for
    pressure). With the reversal A is the temperature anomaly itself, so a low PIT
    means colder than forecast, as a low PIT meant 'more downward' in P'."""
    def load():
        fc = pd.read_parquet(T.FC)
        fc["init"] = pd.to_datetime(fc["init"])
        fc["hyear"] = fc["init"].dt.year
        fc["md"] = pd.to_datetime(fc["model_date"]).dt.strftime("%m-%d")
        g = fc.groupby(["md", "lead_day", "hyear"])[region].agg(["sum", "count"])
        tot = g.groupby(level=[0, 1]).sum()
        loo = tot.reindex(g.index.droplevel(2)).values - g.values
        clim = pd.Series(loo[:, 0] / loo[:, 1], index=g.index)
        anom = fc[region].values - clim.reindex(
            pd.MultiIndex.from_frame(fc[["md", "lead_day", "hyear"]])).values
        fc["anom"] = -anom
        ob = pd.read_parquet(OBS)
        ob = pd.Series(-ob[region].values, index=pd.to_datetime(ob["date"]))
        if not ob.index.is_unique:
            raise ValueError("repeated dates in the ERA5 regional series")
        return fc, ob
    return load


def holm(ps):
    order = np.argsort(ps)
    adj, run = np.empty(len(ps)), 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - rank) * ps[i]))
        adj[i] = run
    return adj


def main():
    missing = [c for c, f in FILES.items() if not f.exists()]
    if missing or not OBS.exists():
        sys.exit(f"missing inputs: {missing} {'' if OBS.exists() else OBS.name}")
    T.WIN = WIN_T
    MM.FILES = FILES
    res = {"plan_approved": "2026-09-26", "primary": PRIMARY, "window_days": list(WIN_T),
           "n_null": MM.N_NULL, "seed": MM.SEED, "inputs": [str(f.relative_to(ROOT)) for f in FILES.values()]
           + [str(OBS.relative_to(ROOT))], "regions": {}}
    td = Path(tempfile.mkdtemp())
    MM.RESULTS = td
    for r in REGIONS:
        T.load = loader(r)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            MM.main()
        tmp = td / f"{MM.NAME}.json"          # MM.main writes here (RESULTS redirected)
        out = json.loads(tmp.read_text())
        tmp.unlink()
        res["regions"][r] = {"multimodel": out["multimodel"], "centres": out["centres"]}
        c = out["multimodel"]["confirmatory"]
        print(f"{r:10s} H1-T PIT {c.get('H1_mean_pit')} (null {c.get('H1_null_mean_pit')}) "
              f"p={c.get('H1_p')}  centred p={c.get('H1c_p')}  LOO max p={c.get('H1_p_loo_max')}  "
              f"D r={c.get('D_r')} p={c.get('D_p_conditional')}  H2 p={c.get('H2_p')}  "
              f"A ens/obs {c.get('A_ens_mean')}/{c.get('A_obs_mean')} K", flush=True)
    ps = np.array([res["regions"][r]["multimodel"]["confirmatory"]["H1_p"] for r in REGIONS], float)
    res["H1_holm_over_regions"] = dict(zip(REGIONS, [round(float(x), 4) for x in holm(ps)]))
    res["falsifier_primary"] = {"rule": "H1-T p >= 0.05 in NEURASIA",
                                "met": bool(res["regions"][PRIMARY]["multimodel"]["confirmatory"]["H1_p"] >= 0.05)}
    print("Holm:", res["H1_holm_over_regions"], " falsifier met:", res["falsifier_primary"]["met"])
    (RESULTS / "s2s_regional_test.json").write_text(json.dumps(res, indent=2, default=str),
                                          encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")


if __name__ == "__main__":
    main()
