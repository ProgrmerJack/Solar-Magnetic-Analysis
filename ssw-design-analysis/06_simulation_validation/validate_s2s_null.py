#!/usr/bin/env python3
"""
validate_s2s_null.py
====================
Is the discrimination test in s2s_forecast_test.py calibrated, and how much power
does it have? Known truth, on the real forecast layout and the real ERA5 outcomes.

Synthetic forecasts replace the S2S values (same model dates, starts, members,
leads) and the test is run exactly as on the real data:

  noise        forecasts carry no information              -> p ~ U(0,1)
  equal_skill  forecasts = outcome + noise everywhere       -> p ~ U(0,1): skill that
               is not SSW-specific must not be flagged
  ssw_skill    outcome in the forecast ONLY for starts 2-17 d before a catalogued
               SSW, matched-variance noise elsewhere         -> power

and, for the record, the plan's first null design (NULL_MODE "year_shift") on
noise, which this check rejected (20% false positives at the short lead when
first run in scratch on 2026-09-25).

Output: results/current/3_calibration/validate_s2s_null.json
"""
import contextlib
import functools
import io
import json
import sys
import tempfile
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "3_calibration"
sys.path.insert(0, str(HERE.parents[0] / "07_physical_decomposition"))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
import s2s_forecast_test as T                        # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

N_SEEDS = 40
N_NULL = 1000
NAME = "validate_s2s_null"
CASES = [("calendar_window", "noise"), ("calendar_window", "equal_skill"),
         ("calendar_window", "ssw_skill"), ("year_shift", "noise")]
KEYS = ("r_ensmean_vs_obs", "p_r", "p_r_conditional", "BSS_P_dw", "p_BSS",
        "p_BSS_conditional")


WORKERS = 16
_G = {}


def _setup():
    """Per-process state: the layout, the outcome series and a cached obs_anom
    (the observed outcome does not depend on the forecast)."""
    base = pd.read_parquet(T.FC)[["model_date", "init", "member", "lead_day", "valid"]]
    ob = pd.read_parquet(T.OBS)
    ob = ob[pd.to_datetime(ob["time"]).dt.hour == 0].set_index("time")["psl_cap_N"]
    valid = pd.to_datetime(base["valid"])
    truth = valid.map(ob).fillna(float(ob.mean())).values
    # the real outcome from ANOTHER year at the same calendar dates: same
    # variance and persistence as the truth, no information about it
    other = pd.Series(pd.DatetimeIndex(valid) + pd.DateOffset(years=1)).map(ob)
    other = other.fillna(pd.Series(pd.DatetimeIndex(valid) - pd.DateOffset(years=1)).map(ob))
    other = other.fillna(float(ob.mean())).values
    # one error shared by the 11 members of each (start, lead), keyed explicitly:
    # the file is sorted by init, member, lead, so reshaping by 11 rows was wrong
    key, _ = pd.factorize(base["init"].astype(str) + "|" + base["lead_day"].astype(str))
    ini = pd.to_datetime(base["init"])
    ssw_start = np.zeros(len(base), bool)
    for o in load_catalogue("primary"):
        k = (o - ini).dt.days
        ssw_start |= ((k >= 2) & (k <= 17)).values
    orig = T.obs_anom

    @functools.lru_cache(maxsize=None)
    def oa(a, offs, hy):
        return orig(ob, [a], list(offs), list(hy))[a]
    T.obs_anom = lambda o, anchors, offs, hy: {
        a: oa(a, tuple(int(x) for x in offs), tuple(hy)) for a in anchors}
    T.N_NULL, T.ROOT = N_NULL, Path("/")
    _G.update(base=base, truth=truth, other=other, key=key, ssw_start=ssw_start,
              td=Path(tempfile.mkdtemp()))


def _one(task):
    mode, h0, seed = task
    if not _G:
        _setup()
    base, truth, other, key = _G["base"], _G["truth"], _G["other"], _G["key"]
    n = len(base)
    rng = np.random.default_rng(zlib.crc32(f"{NAME}|{mode}|{h0}|{seed}".encode()))
    shared = rng.normal(0, 2500, key.max() + 1)[key]
    if h0 == "noise":
        v = 101300 + rng.normal(0, 900, n)
    elif h0 == "equal_skill":
        v = truth + rng.normal(0, 2500, n) + shared
    else:
        noise = rng.normal(0, 2500, n) + shared
        v = np.where(_G["ssw_start"], truth + noise, other + noise)
    tag = f"{mode}_{h0}_{seed}"
    f = _G["td"] / f"{tag}.parquet"
    base.assign(psl_cap_N=v).to_parquet(f)
    T.FC, T.NAME, T.RESULTS, T.NULL_MODE = f, tag, _G["td"], mode
    with contextlib.redirect_stdout(io.StringIO()):
        T.main()
    d = json.loads((_G["td"] / f"{tag}.json").read_text())
    return (mode, h0, seed), {b: {k: d["bins"][b]["discrimination"][k] for k in KEYS}
                              for b in d["bins"]}


def main():
    from concurrent.futures import ProcessPoolExecutor
    fc_path = T.FC
    res = {"n_seeds": N_SEEDS, "n_null": N_NULL, "forecast_layout": str(fc_path.relative_to(ROOT)),
           "synthetic": {"noise": "N(0, 900 Pa) around 101300 Pa",
                         "equal_skill": "outcome + N(0, 2500) per member + N(0, 2500) shared by the "
                                        "members of each (start, lead)",
                         "ssw_skill": "as equal_skill at SSW starts; elsewhere the ERA5 outcome of "
                                      "another year at the same calendar dates, plus the same noise"},
           "workers": WORKERS, "cases": {}}
    tasks = [(m, h, s) for m, h in CASES for s in range(N_SEEDS)]
    with ProcessPoolExecutor(WORKERS) as ex:
        got = dict(ex.map(_one, tasks))
    for mode, h0 in CASES:
        rows = [got[(mode, h0, s)] for s in range(N_SEEDS)]
        out = {}
        for b in ("short", "mid"):
            for k in ("p_r", "p_r_conditional", "p_BSS", "p_BSS_conditional"):
                p = np.array([r[b][k] for r in rows], float)
                out[f"{b}.{k}"] = {"rate_below_0.05": round(float(np.mean(p < 0.05)), 3),
                                   "mean": round(float(np.nanmean(p)), 3)}
        res["cases"][f"{mode}|{h0}"] = out
        print(f"{mode:16s} {h0:12s} " + "  ".join(
            f"{b} r-cond {out[f'{b}.p_r_conditional']['rate_below_0.05']:.2f} "
            f"(uncond {out[f'{b}.p_r']['rate_below_0.05']:.2f}) BSS-cond "
            f"{out[f'{b}.p_BSS_conditional']['rate_below_0.05']:.2f}"
            for b in ("short", "mid")), flush=True)
    res["binomial_se_at_0.05"] = round(float(np.sqrt(0.05 * 0.95 / N_SEEDS)), 3)
    (RESULTS / "validate_s2s_null.json").write_text(json.dumps(res, indent=2),
                                                    encoding="utf8", newline="\n")
    print("Saved -> validate_s2s_null.json")


if __name__ == "__main__":
    main()
