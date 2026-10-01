#!/usr/bin/env python3
"""
shift_rule_forecast.py
======================
A COLD-RISK FORECAST FROM THE SHIFT ALONE: DOES A MODEL-DERIVED PROBABILITY VERIFY
AGAINST OBSERVED SSWs, AND WHAT IS IT WORTH TO A DECISION-MAKER?

Plan approved 2026-10-01. Written before the verification below was run.

QUESTION
  The paper's practical claim is that after an SSW the product is the shifted
  distribution, not a class. Test it as a forecast: fix, from the models alone, the
  probability that the northern-Eurasian (50-65N, 10-130E) fortnight (days 8-24
  after onset) is colder than the q-quantile of the no-SSW distribution, for
  q = 0.05, 0.10, 0.20 -- the "shift rule" -- and verify it, untuned, against the
  observed SSWs. Also for a negative annular mode over days 8-52.

MODEL PROBABILITIES (no observational input)
  Cold: SNAPSI, every centre x initialisation pair of the regional test
    (snapsi_regional_test.build_pairs layout): share of nudged members below the
    control ensemble's q-quantile, mean over pairs; centre-bootstrap 95% interval
    (10,000); also without ECCC.
  Negative annular mode: CMIP6, share of the 1,517 events with a negative days
    8-52 annular-mode anomaly (forecast_value_test data).

OBSERVED VERIFICATION
  ERA5 regional temperature anomalies (era5_regional_test definitions: 31-day
  smoothed 1959-2022 day-of-year climatology; event-free = November-March days more
  than 135 d from every catalogued onset; q-quantiles from the event-free days).
  Verification set V1 (primary): the 39 events of the regional test (WeatherBench 2
  series). V2 (secondary): V1 plus the 2023-2024 events with the spliced ARCO-ERA5
  series of s2s_heldout_test.py. Annular mode: CPC daily AO, days 8-52, 43 events,
  negative if the window mean < 0; climatological base rate from event-free days.
  PROSPECTIVE (run only with --prospective, after s2s_heldout2026_test.py): the
  4 March 2026 SSW, whose 2 m temperature outcome had not been examined when this
  was written (the CPC AO for 30-31 March 2026 had been seen; disclosed there).

METRICS (per q)
  observed frequency f with Wilson 95% interval; two-sided binomial p of the count
  under the shift rule and under climatology (q); Brier skill of the rule against
  climatology over the events; risk ratio f/q against the rule's p/q; relative
  economic value for cost/loss ratios alpha (a user acting when the forecast
  probability exceeds alpha), with event-bootstrap intervals at alpha = q, 2q.
READING (fixed now): the rule is calibrated if the binomial p under it is > 0.05,
  and adds information if the count rejects climatology; value is reported as a
  curve. Only the predictor is out of sample: the observed frequencies at q = 0.10
  had been computed before (era5_regional_test).

Output: results/current/6_predictability/shift_rule_forecast.json
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
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
import snapsi_regional_test as SRT                  # noqa: E402
import era5_regional_test as ERT                     # noqa: E402
import within_model_check as W                       # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

NAME = "shift_rule_forecast"
SEED = zlib.crc32(NAME.encode()) % (2 ** 32)
QS = (0.05, 0.10, 0.20)
N_BOOT = 10000
T_WIN = (8, 24)
AO_WIN = (8, 52)
ZONE_SEP = 135
REG = "NEURASIA"
ALPHAS = np.round(np.arange(0.02, 0.81, 0.02), 2)


def wilson(k, n):
    if n == 0:
        return [None, None]
    z = 1.96; p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [round(c - h, 4), round(c + h, 4)]


def snapsi_probs(rng):
    """Share of nudged members below the control's q-quantile, per pair."""
    centres = sorted({p.name.split("_")[0] for p in SRT.TAS.glob("*.npz")})
    rows = []
    for c in centres:
        bad, _ = SRT.L.corruption_guard(c)
        if bad:
            continue
        for init in SRT.L.ONSET:
            tn, tc = SRT.region_series(c, "nudged", init), SRT.region_series(c, "control", init)
            if len(tn) < 10 or len(tc) < 10:
                continue
            row = {"centre": c, "init": init}
            for q in QS:
                row[f"p{q}"] = float(np.mean(tn[REG] < np.quantile(tc[REG], q)))
            rows.append(row)
    df = pd.DataFrame(rows)
    out = {}
    for sub, d in (("all", df), ("excl_ECCC", df[df.centre != "ECCC"])):
        o = {}
        cs = d["centre"].unique()
        for q in QS:
            col = f"p{q}"
            by = {c: d.loc[d.centre == c, col].values for c in cs}
            bs = [np.mean(np.concatenate([by[c] for c in rng.choice(cs, len(cs))])) for _ in range(N_BOOT)]
            o[str(q)] = {"p": round(float(d[col].mean()), 4),
                         "ci95": [round(float(x), 4) for x in np.quantile(bs, [0.025, 0.975])]}
        out[sub] = {"n_pairs": int(len(d)), "n_centres": int(len(cs)), "q": o}
    return out


def cmip6_p_negative():
    ys = []
    for f in sorted(W.RAW.glob("*_zm.nc")):
        c = W.prepare_member(f)
        if c is None:
            continue
        ys.append(W.C6.anom(c["am"], c["on"], c["cl"], W.OUT_WIN))
    y = np.concatenate(ys)
    y = y[np.isfinite(y)]
    return {"n_events": int(len(y)), "p_negative": round(float(np.mean(y < 0)), 4)}


def t_anomalies(spliced):
    f = ING / ("era5_t2m_regions_daily_spliced.parquet" if spliced else "era5_t2m_regions_daily.parquet")
    t = pd.read_parquet(f).set_index("date")[[REG]]
    t.index = pd.to_datetime(t.index)
    base = t[(t.index.year >= 1959) & (t.index.year <= 2022)]
    doy = base.groupby(base.index.dayofyear).mean().reindex(range(1, 367))
    doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    doy.index = range(1, 367)
    return (t - doy.reindex(t.index.dayofyear).values)[REG]


def window(s, a, w):
    d = s.reindex(pd.date_range(a + pd.Timedelta(days=w[0]), a + pd.Timedelta(days=w[1])))
    return float(d.mean()) if d.notna().all() else np.nan


def rev_curve(p_rule, q, f):
    """Relative economic value of acting when p > alpha, on SSW occasions with true
    frequency f, against a user who knows only climatology q."""
    out = []
    for a in ALPHAS:
        e_clim = a if a < q else f            # expense in units of L
        e_rule = a if a < p_rule else f
        e_perf = f * a
        den = e_clim - e_perf
        out.append(round(float((e_clim - e_rule) / den), 4) if abs(den) > 1e-12 else None)
    return out


def verify(events, outcomes, p_rule, q, rng, label):
    y = np.asarray(outcomes, float)
    k, n = int(y.sum()), int(len(y))
    bs_rule = np.mean((p_rule - y) ** 2); bs_clim = np.mean((q - y) ** 2)
    boot_f = np.array([rng.choice(y, n).mean() for _ in range(N_BOOT)])
    rev_at = {}
    for a in (q, 2 * q):
        vals = [rev_curve(p_rule, q, fb)[int(np.argmin(np.abs(ALPHAS - a)))] for fb in boot_f[:2000]]
        vals = np.array([v for v in vals if v is not None], float)
        rev_at[str(round(a, 2))] = {"rev": rev_curve(p_rule, q, k / n)[int(np.argmin(np.abs(ALPHAS - a)))],
                                    "ci95": [round(float(x), 4) for x in np.quantile(vals, [0.025, 0.975])] if len(vals) else None}
    return {"set": label, "n": n, "count": k, "f": round(k / n, 4), "f_wilson95": wilson(k, n),
            "p_rule": p_rule, "q": q,
            "binom_p_under_rule": round(float(stats.binomtest(k, n, p_rule).pvalue), 4),
            "binom_p_under_climatology": round(float(stats.binomtest(k, n, q).pvalue), 4),
            "brier_rule": round(float(bs_rule), 4), "brier_climatology": round(float(bs_clim), 4),
            "bss_vs_climatology": round(float(1 - bs_rule / bs_clim), 4),
            "risk_ratio_observed": round((k / n) / q, 3), "risk_ratio_rule": round(p_rule / q, 3),
            "rev_curve": dict(zip([str(a) for a in ALPHAS], rev_curve(p_rule, q, k / n))),
            "rev_at": rev_at, "events": [str(pd.Timestamp(e).date()) for e in events]}


def main():
    rng = np.random.default_rng(SEED)
    prospective = "--prospective" in sys.argv
    res = {"plan_approved": "2026-10-01", "seed": SEED, "n_boot": N_BOOT, "qs": list(QS),
           "region": REG, "t_window": list(T_WIN), "ao_window": list(AO_WIN)}
    sp = snapsi_probs(rng)
    res["model_probabilities_snapsi"] = sp
    res["model_p_negative_am_cmip6"] = cmip6_p_negative()
    print("model", json.dumps(sp), res["model_p_negative_am_cmip6"], flush=True)
    cat = load_catalogue("primary")
    res["verification"] = {}
    for label, spliced, last in (("V1", False, pd.Timestamp("2023-01-10")), ("V2", True, pd.Timestamp("2024-04-30"))):
        ta = t_anomalies(spliced)
        ta = ta[ta.index <= last]
        free_days = [d for d in pd.date_range("1959-01-01", "2022-12-31")
                     if d.month in (11, 12, 1, 2, 3) and np.all(np.abs((cat - d).days) > ZONE_SEP)]
        ta_v1 = t_anomalies(False)
        free = np.array([window(ta_v1, d, T_WIN) for d in free_days]); free = free[np.isfinite(free)]
        ev = [o for o in cat if o <= last - pd.Timedelta(days=T_WIN[1])]
        vals = np.array([window(ta, o, T_WIN) for o in ev])
        ok = np.isfinite(vals)
        ev = [e for e, k in zip(ev, ok) if k]; vals = vals[ok]
        res["verification"][label] = {}
        for q in QS:
            thr = float(np.quantile(free, q))
            p_rule = sp["all"]["q"][str(q)]["p"]
            r = verify(ev, vals < thr, p_rule, q, rng, label)
            r["threshold_K"] = round(thr, 3)
            r["p_rule_excl_ECCC"] = sp["excl_ECCC"]["q"][str(q)]["p"]
            r["binom_p_under_rule_excl_ECCC"] = round(float(stats.binomtest(r["count"], r["n"], r["p_rule_excl_ECCC"]).pvalue), 4)
            res["verification"][label][str(q)] = r
            print(label, q, {k: r[k] for k in ("n", "count", "f", "f_wilson95", "p_rule", "binom_p_under_rule",
                                                "binom_p_under_climatology", "bss_vs_climatology", "rev_at")}, flush=True)
    # annular mode: CPC AO days 8-52
    import re
    recs = [re.match(r"\s*(\d{4})\s+(\d+)\s+(\d+)\s*(-?\d+\.\d+)", ln)
            for ln in (ROOT / "data/processed/atmospheric/ao_daily_cpc.txt").read_text().splitlines()]
    recs = [m.groups() for m in recs if m]
    ao = pd.Series([float(v) for *_, v in recs], index=pd.to_datetime([f"{y}-{m}-{d}" for y, m, d, _ in recs]))
    ao = ao.where(ao > -90)
    ev = [o for o in cat if o <= pd.Timestamp("2025-01-01")]
    yv = np.array([window(ao, o, AO_WIN) for o in ev]); okk = np.isfinite(yv)
    free_days = [d for d in pd.date_range("1958-01-01", "2024-12-31")
                 if d.month in (11, 12, 1, 2, 3) and np.all(np.abs((cat - d).days) > ZONE_SEP)]
    fr = np.array([window(ao, d, AO_WIN) for d in free_days]); fr = fr[np.isfinite(fr)]
    base = float(np.mean(fr < 0))
    res["verification"]["AO_negative"] = verify([e for e, k in zip(ev, okk) if k], yv[okk] < 0,
                                                res["model_p_negative_am_cmip6"]["p_negative"], base, rng, "AO")
    res["verification"]["AO_negative"]["climatological_base_rate"] = round(base, 4)
    print("AO", {k: res["verification"]["AO_negative"][k] for k in ("n", "count", "f", "p_rule", "q",
                                                                      "binom_p_under_rule", "binom_p_under_climatology",
                                                                      "bss_vs_climatology")}, flush=True)
    if prospective:
        late = pd.read_parquet(ING / "era5_t2m_regions_daily_arco_late.parquet").set_index("date")[[REG]]
        late.index = pd.to_datetime(late.index)
        spl = json.loads((RESULTS / "s2s_heldout_test.json").read_text())["splice"][REG]["offset"]
        t_all = pd.concat([pd.read_parquet(ING / "era5_t2m_regions_daily_spliced.parquet").set_index("date")[[REG]],
                           late[late.index > pd.Timestamp("2024-04-30")] - spl])
        t_all.index = pd.to_datetime(t_all.index)
        base_ = t_all[(t_all.index.year >= 1959) & (t_all.index.year <= 2022)]
        doy = base_.groupby(base_.index.dayofyear).mean().reindex(range(1, 367))
        doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
        doy.index = range(1, 367)
        ta = (t_all - doy.reindex(t_all.index.dayofyear).values)[REG]
        o26 = pd.Timestamp("2026-03-04")
        v = window(ta, o26, T_WIN)
        pr = {}
        free = np.array([window(t_anomalies(False), d, T_WIN) for d in free_days if d.year >= 1959])
        free = free[np.isfinite(free)]
        for q in QS:
            pr[str(q)] = {"cold": bool(v < float(np.quantile(free, q))), "p_rule": sp["all"]["q"][str(q)]["p"]}
        res["prospective_2026"] = {"anomaly_K": round(v, 3), "by_q": pr,
                                   "ao_window_mean": None if not np.isfinite(window(ao, o26, AO_WIN)) else round(window(ao, o26, AO_WIN), 3)}
        print("prospective 2026", res["prospective_2026"], flush=True)
    (RESULTS / "shift_rule_forecast.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> shift_rule_forecast.json")


if __name__ == "__main__":
    main()
