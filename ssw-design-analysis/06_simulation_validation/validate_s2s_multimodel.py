#!/usr/bin/env python3
"""
validate_s2s_multimodel.py
==========================
Known-truth calibration of s2s_multimodel_test.py on every centre's real layout
(starts, members, leads) and the real ERA5 outcome.

  equal_skill  forecast = outcome + member noise + noise shared by the members of
               each (start, lead), at every date     -> D, H1, H2 reject ~5%
  noise        forecast carries no information       -> D, H1, H2 reject ~5%
  damped       forecast = 0.6 x outcome + noise, at every date. The shift is
               under-forecast but NOT specifically at SSWs; H1 and H2 are expected
               to fire here, which is what they are designed to detect (an
               under-forecast shift after SSWs), whatever its cause.

Rates are for the primary statistic, the confirmatory multi-model mean.
Output: results/current/3_calibration/validate_s2s_multimodel.json
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

N_SEEDS = 24
N_NULL = 1000
WORKERS = 12
NAME = "validate_s2s_multimodel"
CASES = ("equal_skill", "noise", "damped")
_G = {}


def _setup():
    ob = pd.read_parquet(T.OBS)
    ob = ob[pd.to_datetime(ob["time"]).dt.hour == 0].set_index("time")["psl_cap_N"]
    lay = {}
    for c, f in MM.FILES.items():
        if f.exists():
            b = pd.read_parquet(f)[["model_date", "init", "member", "lead_day", "valid"]]
            truth = pd.to_datetime(b["valid"]).map(ob).fillna(float(ob.mean())).values
            key, _ = pd.factorize(b["init"].astype(str) + "|" + b["lead_day"].astype(str))
            lay[c] = (b, truth, key)
    orig = T.obs_anom

    @functools.lru_cache(maxsize=None)
    def oa(a, offs, hy, mo):
        return orig(ob, [a], list(offs), list(hy), mo)[a]
    T.obs_anom = lambda o, anchors, offs, hy, mo=15: {
        a: oa(a, tuple(int(x) for x in offs), tuple(hy), mo) for a in anchors}
    _G.update(lay=lay, td=Path(tempfile.mkdtemp()))


def _one(task):
    case, seed = task
    if not _G:
        _setup()
    files = {}
    for c, (b, truth, key) in _G["lay"].items():
        rng = np.random.default_rng(zlib.crc32(f"{NAME}|{case}|{c}|{seed}".encode()))
        n = len(b)
        shared = rng.normal(0, 2500, key.max() + 1)[key]
        noise = rng.normal(0, 2500, n) + shared
        if case == "equal_skill":
            v = truth + noise
        elif case == "damped":
            v = 101300 + 0.6 * (truth - 101300) + noise
        else:
            v = 101300 + rng.normal(0, 900, n)
        f = _G["td"] / f"{case}_{seed}_{c}.parquet"
        b.assign(psl_cap_N=v).to_parquet(f)
        files[c] = f
    MM.FILES, MM.N_NULL, MM.RESULTS = files, N_NULL, _G["td"]
    MM.ROOT = Path("/")
    with contextlib.redirect_stdout(io.StringIO()):
        MM.main()
    d = json.loads((_G["td"] / "s2s_multimodel_test.json").read_text())["multimodel"]["confirmatory"]
    return case, seed, {k: d.get(k) for k in ("D_p_conditional", "H1_p", "H2_p")}


def main():
    ctx = mp.get_context("fork")
    tasks = [(c, s) for c in CASES for s in range(N_SEEDS)]
    with ProcessPoolExecutor(WORKERS, mp_context=ctx) as ex:
        got = list(ex.map(_one, tasks))
    res = {"n_seeds": N_SEEDS, "n_null": N_NULL, "centres": [c for c in MM.FILES],
           "binomial_se_at_0.05": round(float(np.sqrt(0.05 * 0.95 / N_SEEDS)), 3), "cases": {}}
    for case in CASES:
        rows = [r for c, s, r in got if c == case]
        out = {}
        for k in ("D_p_conditional", "H1_p", "H2_p"):
            p = np.array([r[k] for r in rows if r[k] is not None], float)
            out[k] = {"rate_below_0.05": round(float(np.mean(p < 0.05)), 3),
                      "mean": round(float(np.mean(p)), 3), "n": int(len(p))}
        res["cases"][case] = out
        print(f"{case:12s} " + "  ".join(f"{k} {v['rate_below_0.05']:.2f} (mean p {v['mean']:.2f})"
                                          for k, v in out.items()), flush=True)
    (RESULTS / "validate_s2s_multimodel.json").write_text(json.dumps(res, indent=2),
                                                          encoding="utf8", newline="\n")
    print("Saved -> validate_s2s_multimodel.json")


if __name__ == "__main__":
    main()
