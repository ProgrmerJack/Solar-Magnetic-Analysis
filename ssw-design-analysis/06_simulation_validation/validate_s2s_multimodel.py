#!/usr/bin/env python3
"""
validate_s2s_multimodel.py
==========================
Known-truth calibration of s2s_multimodel_test.py on every centre's real layout
(starts, members, leads) and the real ERA5 outcome.

  equal_skill  forecast = outcome + member noise + noise shared by the members of
               each (start, lead), at every date     -> D, H1, H2 reject ~5%
  noise        forecast carries no information       -> D rejects ~5%. H1 is NOT
               null here: a forecast that knows nothing misses the shift after
               SSWs, which is what it detects (fired in 100/100). H2 did not fire
               (0/100; mean p 0.19): its power is low. (The first draft of this
               docstring listed noise as a null for H1/H2.)
  damped       forecast = 0.6 x outcome + noise, at every date. The shift is
               under-forecast but NOT specifically at SSWs; H1 and H2 are expected
               to fire here, which is what they are designed to detect (an
               under-forecast shift after SSWs), whatever its cause.
  ssw_skill    (added 2026-09-26, after the real result, because D's power in the
               multi-model design had not been measured) as equal_skill for starts
               2-9 d before a catalogued SSW; elsewhere the ERA5 outcome of another
               year at the same calendar dates plus the same noise -> power of D

Rates are for the primary statistic, the confirmatory multi-model mean.
Output: results/current/3_calibration/validate_s2s_multimodel.json

--var t2m (added 2026-09-26): the same cases for s2s_regional_test.py, with the
real ERA5 northern-Eurasia temperature (days +8..+24 daily means) as the outcome
and every system's real t2m layout; all noise amplitudes are scaled by the ratio
of the winter daily s.d. of the temperature series to that of polar-cap pressure,
so each case has the same signal-to-noise as in the pressure calibration.
Output: results/current/3_calibration/validate_s2s_multimodel_t2m.json
"""
import contextlib
import functools
import io
import json
import multiprocessing as mp
import sys
import tempfile
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "3_calibration"
sys.path.insert(0, str(HERE.parents[0] / "07_physical_decomposition"))
import s2s_multimodel_test as MM                     # noqa: E402
import s2s_forecast_test as T                        # noqa: E402
import s2s_regional_test as SR                       # noqa: E402
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
from build_catalogue import load_catalogue          # noqa: E402

N_SEEDS = 100
N_NULL = 1000
WORKERS = 16
NAME = "validate_s2s_multimodel"
CASES = ("equal_skill", "noise", "damped", "ssw_skill")
VAR = "msl"                  # "t2m" with --var t2m
_G = {}


def _winter_sd(series):
    m = pd.DatetimeIndex(series.index).month
    return float(series[np.isin(m, (12, 1, 2, 3))].std())


def _setup():
    ob = pd.read_parquet(T.OBS)
    ob = ob[pd.to_datetime(ob["time"]).dt.hour == 0].set_index("time")["psl_cap_N"]
    scale, col = 1.0, "psl_cap_N"
    if VAR == "t2m":
        et = pd.read_parquet(SR.OBS)
        et = pd.Series(et[SR.PRIMARY].values, index=pd.to_datetime(et["date"]))
        scale, col = _winter_sd(et) / _winter_sd(ob), SR.PRIMARY
        ob = et
        MM.FILES = SR.FILES
        T.WIN = SR.WIN_T
        T.load = SR.loader(SR.PRIMARY)
    lay = {}
    for c, f in MM.FILES.items():
        if f.exists():
            b = pd.read_parquet(f)[["model_date", "init", "member", "lead_day", "valid"]]
            truth = pd.to_datetime(b["valid"]).map(ob).fillna(float(ob.mean())).values
            key, _ = pd.factorize(b["init"].astype(str) + "|" + b["lead_day"].astype(str))
            valid = pd.to_datetime(b["valid"])
            other = pd.Series(pd.DatetimeIndex(valid) + pd.DateOffset(years=1)).map(ob)
            other = other.fillna(pd.Series(pd.DatetimeIndex(valid) - pd.DateOffset(years=1)).map(ob))
            other = other.fillna(float(ob.mean())).values
            ini = pd.to_datetime(b["init"])
            ssw = np.zeros(len(b), bool)
            for o in load_catalogue("primary"):
                k = (o - ini).dt.days
                ssw |= ((k >= MM.K_SHORT[0]) & (k <= MM.K_SHORT[1])).values
            lay[c] = (b, truth, key, other, ssw)
    orig = T.obs_anom
    # the loader hands T.start_stats the SIGN-REVERSED temperature (see
    # s2s_regional_test.loader); the cached observed anomaly must match it
    ob_used = -ob if VAR == "t2m" else ob

    @functools.lru_cache(maxsize=None)
    def oa(a, offs, hy, mo):
        return orig(ob_used, [a], list(offs), list(hy), mo)[a]
    T.obs_anom = lambda o, anchors, offs, hy, mo=15: {
        a: oa(a, tuple(int(x) for x in offs), tuple(hy), mo) for a in anchors}
    _G.update(lay=lay, td=Path(tempfile.mkdtemp()), scale=scale, col=col,
              mean=float(ob.mean()))


def _one(task):
    case, seed = task
    if not _G:
        _setup()
    files = {}
    for c, (b, truth, key, other, ssw) in _G["lay"].items():
        rng = np.random.default_rng(zlib.crc32(f"{NAME}|{case}|{c}|{seed}".encode()))
        n = len(b)
        sc, mu = _G["scale"], (101300 if VAR == "msl" else _G["mean"])
        shared = rng.normal(0, 2500 * sc, key.max() + 1)[key]
        noise = rng.normal(0, 2500 * sc, n) + shared
        if case == "equal_skill":
            v = truth + noise
        elif case == "ssw_skill":
            v = np.where(ssw, truth, other) + noise
        elif case == "damped":
            v = mu + 0.6 * (truth - mu) + noise
        else:
            v = mu + rng.normal(0, 900 * sc, n)
        f = _G["td"] / f"{case}_{seed}_{c}.parquet"
        b.assign(**{_G["col"]: v}).to_parquet(f)
        files[c] = f
    MM.FILES, MM.N_NULL, MM.RESULTS = files, N_NULL, _G["td"]
    MM.ROOT = Path("/")
    with contextlib.redirect_stdout(io.StringIO()):
        MM.main()
    out = _G["td"] / "s2s_multimodel_test.json"
    d = json.loads(out.read_text())["multimodel"]["confirmatory"]
    # /tmp is RAM-backed here: 400 runs x 10 systems of synthetic forecasts filled
    # it (16 GB) and killed the first four-case run
    for f in [*files.values(), out]:
        f.unlink()
    return case, seed, {k: d.get(k) for k in ("D_p_conditional", "H1_p", "H2_p", "H1c_p")}


def main():
    ctx = mp.get_context("fork")
    tasks = [(c, s) for c in CASES for s in range(N_SEEDS)]
    with ProcessPoolExecutor(WORKERS, mp_context=ctx) as ex:
        got = list(ex.map(_one, tasks))
    if VAR == "t2m":
        MM.FILES = SR.FILES                           # the layouts to be used
    res = {"var": VAR, "n_seeds": N_SEEDS, "n_null": N_NULL, "centres": [c for c in MM.FILES],
           "binomial_se_at_0.05": round(float(np.sqrt(0.05 * 0.95 / N_SEEDS)), 3), "cases": {}}
    for case in CASES:
        rows = [r for c, s, r in got if c == case]
        out = {}
        for k in ("D_p_conditional", "H1_p", "H2_p", "H1c_p"):
            p = np.array([r[k] for r in rows if r[k] is not None], float)
            out[k] = {"rate_below_0.05": round(float(np.mean(p < 0.05)), 3),
                      "mean": round(float(np.mean(p)), 3), "n": int(len(p))}
        res["cases"][case] = out
        print(f"{case:12s} " + "  ".join(f"{k} {v['rate_below_0.05']:.2f} (mean p {v['mean']:.2f})"
                                          for k, v in out.items()), flush=True)
    out = (RESULTS / "validate_s2s_multimodel.json" if VAR == "msl"
           else RESULTS / "validate_s2s_multimodel_t2m.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(f"Saved -> {out.name}")


if __name__ == "__main__":
    if "--var" in sys.argv:
        VAR = sys.argv[sys.argv.index("--var") + 1]
        assert VAR in ("msl", "t2m"), VAR
    main()
