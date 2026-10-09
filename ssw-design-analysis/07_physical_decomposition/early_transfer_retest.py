#!/usr/bin/env python3
"""
early_transfer_retest.py
========================
DOES THE OUT-OF-SAMPLE FAILURE OF THE SHIFT RULE IN 1941-1958 SURVIVE AN INDEPENDENT
DETECTOR AND AN INDEPENDENT OUTCOME?

Registered 2026-10-09, committed with its inputs' acquisition scripts BEFORE any JRA-3Q wind
or GHCN-Daily temperature value was read. The failure being re-tested (shift_rule_forecast
R6, revision 3): after 12 SSWs detected in ERA5 with onsets 1941-1958, 0 northern-Eurasian
(50-65 N, 10-130 E) days 8-24 windows fell below the period's event-free 10th percentile,
against 0.32 predicted by the SNAPSI shift rule (binomial p = 0.012). Two suspects: ERA5's
sparsely constrained early stratosphere (which events are detected) and its early surface
temperature (the outcome).
INPUTS
  detector  JRA-3Q daily-mean zonal-mean u(10 hPa, 60.0 N) (acquire_jra3q_u10.py), from
            November 1947; the project's Charlton-Polvani detector
            (ensemble_precursor.detect_ssw), unchanged.
  outcome   GHCN-Daily station anomaly over the same box (acquire_ghcnd_neurasia.py),
            base 1941-1958 for the tests.
  ERA5      the R6 onsets and the R6 outcome series (shift_rule_forecast inputs), unchanged.
  rule      the SNAPSI probabilities stored by shift_rule_forecast (0.24, 0.32, 0.45).
GATES (fixed now; a failed gate stops that arm and is reported)
  D  detector: on November 1958 - April 1980, JRA-3Q onsets matched within +-3 days to at
     least 80% of the primary catalogue's onsets in that span; extra onsets reported.
  O  outcome: correlation of station-based (base 1959-2022) and ERA5 WeatherBench 2
     (shift_rule_forecast.t_anom_region, base 1959-2022) days 8-24 window means starting on
     every November-March day of 1959-2022 >= 0.90. The same correlation against the early
     ERA5 series (base 1941-1958) on 1941-1958 is reported, not gated.
TESTS (verification exactly as R6: shift_rule_forecast.obs_verify, event-free November-March
dates of 1941-1958 more than 135 days from every onset, q = 0.05, 0.10, 0.20)
  P1 PRIMARY  JRA-3Q onsets 1 November 1947 - 31 March 1958, station outcome, q = 0.10.
     Exclusion set: JRA-3Q onsets up to 1 June 1959, plus ERA5 onsets before September 1947.
  S1  ERA5 R6 onsets (12), station outcome (isolates the outcome source); ERA5 exclusion set.
  S2  JRA-3Q onsets, ERA5 early outcome (isolates the detector); P1 exclusion set.
  S3  onsets found by both detectors (within +-3 days; ERA5 date), station outcome.
READING (fixed now; P1 at q = 0.10)
  binomial p under the rule < 0.05 with fewer cold windows than predicted: "the failure
    replicates with an independent detector and outcome";
  p >= 0.05: "the failure does not replicate", reported with the probability that this
    test would reject the rule if the true rate were climatology (q), i.e. its power against
    the observed alternative.
  S1-S3 are descriptive: which source carries the failure.
Output: results/current/6_predictability/early_transfer_retest.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RES = ROOT / "results" / "current" / "6_predictability"
OUT = RES / "early_transfer_retest.json"
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
import ensemble_precursor as EP                     # noqa: E402
import shift_rule_forecast as SRF                   # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

NAME = "early_transfer_retest"
SEED = zlib.crc32(NAME.encode()) % (2 ** 32)
QS = (0.05, 0.10, 0.20)
P1_SPAN = (pd.Timestamp("1947-11-01"), pd.Timestamp("1958-03-31"))
VAL_SPAN = (pd.Timestamp("1958-11-01"), pd.Timestamp("1980-04-30"))
GATE_D, GATE_O = 0.80, 0.90


def wmeans(ta, starts, win=SRF.T_WIN):
    return np.array([SRF.window(ta, d, win) for d in starts])


def power_vs_climatology(n, p_rule, q):
    """P(two-sided binomial test of the rule rejects at 0.05) when the true rate is q."""
    return float(sum(stats.binom.pmf(k, n, q) for k in range(n + 1) if stats.binomtest(k, n, p_rule).pvalue < 0.05))


def main():
    rng = np.random.default_rng(SEED)
    res = {"registered": "2026-10-09, committed before any JRA-3Q or GHCN-Daily value was read", "seed": int(SEED)}
    srf = json.loads((RES / "shift_rule_forecast.json").read_text())
    r6 = srf["revision3"]["R6_out_of_sample_1940_1958"]
    pr = {str(q): r6["primary"][str(q)]["p_rule"] for q in QS}
    era5_on = pd.DatetimeIndex(r6["era5_onsets_1940_1958"])
    _, era5_all = SRF.detect_era5()
    # detector
    j = pd.read_parquet(ING / "jra3q_u10_60N_daily.parquet")
    u = pd.Series(j.u10_60N.values, index=pd.to_datetime(j.date))
    jon = pd.DatetimeIndex(EP.detect_ssw(u.values, u.index))
    cat = load_catalogue("primary"); cv = cat[(cat >= VAL_SPAN[0]) & (cat <= VAL_SPAN[1])]
    jv = jon[(jon >= VAL_SPAN[0]) & (jon <= VAL_SPAN[1])]
    hit = [bool(len(jv) and np.min(np.abs((jv - c).days)) <= 3) for c in cv]
    extra = [str(o.date()) for o in jv if np.min(np.abs((cv - o).days)) > 3]
    res["gate_D"] = {"n_catalogue": int(len(cv)), "matched_within_3d": int(sum(hit)), "share": round(float(np.mean(hit)), 4),
                     "missed": [str(c.date()) for c, h in zip(cv, hit) if not h], "extra_jra3q": extra,
                     "passes": bool(np.mean(hit) >= GATE_D)}
    # outcome
    g = pd.read_parquet(ING / "ghcnd_neurasia_daily.parquet"); g.index = pd.to_datetime(g.date)
    st_early, st_mod = g["anom_b1941_1958"], g["anom_b1959_2022"]
    era_mod = SRF.t_anom_region(ING / "era5_t2m_regions_daily.parquet", SRF.REG)
    era_early = SRF.t_anom_region(ING / "era5_t2m_regions_daily_arco_early.parquet", SRF.REG, base=(1941, 1958))
    def corr(a, b, y0, y1):
        days = [d for d in pd.date_range(f"{y0}-01-01", f"{y1}-12-31") if d.month in (11, 12, 1, 2, 3)]
        x, y = wmeans(a, days), wmeans(b, days); k = np.isfinite(x) & np.isfinite(y)
        return round(float(np.corrcoef(x[k], y[k])[0, 1]), 4), int(k.sum())
    r_mod, n_mod = corr(st_mod, era_mod, 1959, 2022); r_early, n_early = corr(st_early, era_early, 1941, 1958)
    res["gate_O"] = {"r_1959_2022": r_mod, "n_windows": n_mod, "passes": bool(r_mod >= GATE_O),
                     "r_1941_1958_vs_early_ERA5_not_gated": r_early, "n_windows_early": n_early,
                     "stations_used_1941_1958": int(g["n_stations_b1941_1958"].max()), "cells_max": int(g["n_cells_b1941_1958"].max())}
    print(json.dumps({"gate_D": res["gate_D"], "gate_O": res["gate_O"]}, indent=1), flush=True)
    p1_on = jon[(jon >= P1_SPAN[0]) & (jon <= P1_SPAN[1])]
    excl_j = jon[jon <= pd.Timestamp("1959-06-01")].append(era5_all[era5_all < pd.Timestamp("1947-09-01")]).sort_values()
    excl_e = era5_all[era5_all <= pd.Timestamp("1959-06-01")]
    both = [e for e in era5_on if len(p1_on) and np.min(np.abs((p1_on - e).days)) <= 3]
    res["onsets"] = {"jra3q_P1": [str(o.date()) for o in p1_on], "era5_R6": [str(o.date()) for o in era5_on],
                     "both_within_3d_era5_dates": [str(o.date()) for o in both]}
    # added after the run (descriptive, inputs only, no outcome): reanalysis agreement and event-set overlap
    eu = pd.read_parquet(ING / "era5_u10_60N_daily_cds.parquet"); eu = pd.Series(eu["u10_60N"].values, index=pd.to_datetime(eu["date"]))
    x = pd.concat([u, eu], axis=1, keys=["j", "e"], sort=True).dropna(); cmp_ = {}
    for a, b in (("1947-11-01", "1958-04-30"), ("1958-11-01", "1980-04-30")):
        y = x[a:b]
        cmp_[f"{a[:4]}-{b[:4]}"] = {"n_days": int(len(y)), "r": round(float(y.j.corr(y.e)), 4), "bias_jra_minus_era5": round(float((y.j - y.e).mean()), 3),
                                     "rms": round(float(np.sqrt(((y.j - y.e) ** 2).mean())), 3), "sign_disagreement": round(float(((y.j < 0) != (y.e < 0)).mean()), 4)}
    e47 = era5_on[era5_on >= pd.Timestamp("1947-11-01")]
    cmp_["era5_onsets_from_nov1947"] = int(len(e47))
    cmp_["of_which_jra3q_within_3d"] = int(sum(np.min(np.abs((p1_on - e).days)) <= 3 for e in e47))
    cmp_["jra3q_onsets_not_in_era5_within_3d"] = int(sum(np.min(np.abs((e47 - o).days)) > 3 for o in p1_on))
    cmp_["missed_catalogue_nearest_jra3q_days"] = {m: int(np.min(np.abs((jv - pd.Timestamp(m)).days))) for m in res["gate_D"]["missed"]}
    res["input_comparison_post_run"] = cmp_
    arms = {}
    if res["gate_D"]["passes"] and res["gate_O"]["passes"]:
        arms["P1"] = SRF.obs_verify(st_early, list(p1_on), excl_j, (1941, 1958), SRF.T_WIN, pr, rng, "P1")
        arms["S3"] = SRF.obs_verify(st_early, both, excl_j, (1941, 1958), SRF.T_WIN, pr, rng, "S3")
    if res["gate_O"]["passes"]:
        arms["S1"] = SRF.obs_verify(st_early, list(era5_on), excl_e, (1941, 1958), SRF.T_WIN, pr, rng, "S1")
    if res["gate_D"]["passes"]:
        arms["S2"] = SRF.obs_verify(era_early, list(p1_on), excl_j, (1941, 1958), SRF.T_WIN, pr, rng, "S2")
    for a in arms.values():
        for q in QS:
            r = a[str(q)]; r["power_vs_climatology"] = round(power_vs_climatology(r["n"], r["p_rule"], q), 4)
    res["arms"] = arms
    if "P1" in arms:
        r = arms["P1"]["0.1"]
        res["reading_P1"] = ("the failure replicates with an independent detector and outcome"
                             if r["binom_p_under_rule"] < 0.05 and r["f"] < r["p_rule"] else "the failure does not replicate")
    else:
        res["reading_P1"] = "not run: gate failed"
    for k, a in arms.items():
        print(k, {q: (a[str(q)]["count"], a[str(q)]["n"], a[str(q)]["binom_p_under_rule"], a[str(q)]["binom_p_under_climatology"],
                      a[str(q)]["power_vs_climatology"]) for q in QS}, flush=True)
    print(res["reading_P1"])
    OUT.write_text(json.dumps(res, indent=1, default=str), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
