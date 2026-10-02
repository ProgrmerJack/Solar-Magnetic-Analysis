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

REVISION-3 ADDENDUM (approved 2026-10-01: "Do all!"; committed BEFORE any ERA5
surface value before 1959 was retrieved, before any observed outcome over days
15-42 after an SSW or in the three other regions was computed in this script, and
before any SNAPSI probability other than northern Eurasia days 8-24 was computed).
Written to the same JSON under "revision3"; the registered numbers above are not
changed.
  R1 STRICT INDEPENDENCE: V1 without the two observed SSWs that SNAPSI imposes
     (12 February 2018, 2 January 2019): 37 events.
  R2 SEVERITY: q = 0.025, 0.05, 0.10, 0.20, 0.33 (SNAPSI p, observed count, Wilson
     interval, binomial p under the rule and under climatology, risk ratios); the
     shift reading predicts risk ratios that rise as q falls.
  R3 REGIONS: R2 for HI_EUROPE, MID_EASIA and MID_NAMER as well. Reading fixed now:
     the regional test found residuals beyond the circulation in all three (most in
     North America), so the rule may fail there; a rejection under the rule is
     reported as a regional failure of the circulation shift, not hidden.
  R4 WEEKS 3-6: windows days 15-28 and 29-42 after onset (q = 0.05, 0.10, 0.20,
     northern Eurasia primary, all regions secondary). SNAPSI p from every pair whose
     members cover the window (region_series' full-coverage rule; weeks 3-4: both
     later initialisations at all nine centres; weeks 5-6: mainly 8 January 2019),
     with the number of pairs and centres; computed only with at least four
     centres. Observed: the same windows after the V1 events, quantiles from
     event-free dates with the same window.
  R5 A SECOND, INDEPENDENT PREDICTOR (CMIP6 x ERA5 event-free): p2(q) = mean over
     the CMIP6 events of Phi((thr_q - a - b N_i) / s_e), where N_i is each CMIP6
     event's days 8-24 annular-mode anomaly (ensemble_precursor index, daily-s.d.
     units) and a, b, s_e the ERA5 regression of northern-Eurasian days 8-24
     temperature anomaly on the days 8-24 NAM (era5_nam_daily nam_1000, divided by
     its event-free November-March daily s.d.) over EVENT-FREE dates only. No
     post-SSW observation enters. Verified on V1 as for the SNAPSI rule.
  R6 OUT OF SAMPLE, 1940-1958: SSWs detected by the project's detector
     (ensemble_precursor.detect_ssw, the CP07 rule validated on NCEP) on the ERA5
     u(10 hPa, 60N) daily mean from the CDS (acquire_era5_u10_cds.py), onsets
     1 November 1940 - 31 March 1958 (before the catalogue begins; none of their
     surface outcomes has been examined). Outcome: ERA5 (ARCO-ERA5, reduced to the
     WeatherBench 2 grid and regions as acquire_era5_arco_extension.py) daily-mean
     2 m temperature. Because of the warming trend, anomalies and thresholds are
     taken within the period: anomaly from the 31-day-smoothed 1941-1958
     day-of-year mean; q-quantiles from that period's event-free November-March
     dates (> 135 d from every detected onset). Primary: northern Eurasia days
     8-24, q = 0.05, 0.10, 0.20, SNAPSI rule; counts with Wilson intervals,
     binomial p under the rule and under climatology; also pooled with V1.
     Secondary: negative polar-cap NAM proxy (ERA5 polar-cap mean sea-level
     pressure, 00 UTC, NAM sign) days 8-52 against the CMIP6 0.74. Sensitivity:
     onsets from 1946 only (ERA5's stratosphere has a cold bias and sparse
     upper-air data before 1946; Soci et al. 2024). Reading fixed now: with about
     ten events the test can only reject gross miscalibration; consistency with the
     rule and the direction against climatology are reported as such.

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


# ------------------------------------------------------------------ revision 3
R3_QS = (0.025, 0.05, 0.10, 0.20, 0.33)
R4_WINS = ((15, 28), (29, 42))
REGIONS = ("NEURASIA", "HI_EUROPE", "MID_EASIA", "MID_NAMER")
SNAPSI_EVENTS = (pd.Timestamp("2018-02-12"), pd.Timestamp("2019-01-02"))
EARLY = (pd.Timestamp("1940-11-01"), pd.Timestamp("1958-03-31"))


def snapsi_probs_gen(rng, win, qs, min_centres=4):
    """SNAPSI p(q) per region for a post-onset window (full-coverage pairs only)."""
    old = SRT.TAS_DAYS
    SRT.TAS_DAYS = win
    try:
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
                for r in REGIONS:
                    for q in qs:
                        row[f"{r}|{q}"] = float(np.mean(tn[r] < np.quantile(tc[r], q)))
                rows.append(row)
    finally:
        SRT.TAS_DAYS = old
    df = pd.DataFrame(rows)
    out = {"n_pairs": int(len(df)), "n_centres": int(df["centre"].nunique()) if len(df) else 0,
           "pairs": sorted(f"{c}/{i}" for c, i in zip(df.get("centre", []), df.get("init", []))), "p": {}}
    if out["n_centres"] < min_centres:
        out["status"] = f"fewer than {min_centres} centres cover the window"
        return out
    cs = df["centre"].unique()
    for r in REGIONS:
        out["p"][r] = {}
        for q in qs:
            col = f"{r}|{q}"
            by = {c: df.loc[df.centre == c, col].values for c in cs}
            bs = [np.mean(np.concatenate([by[c] for c in rng.choice(cs, len(cs))])) for _ in range(N_BOOT)]
            out["p"][r][str(q)] = {"p": round(float(df[col].mean()), 4),
                                   "ci95": [round(float(x), 4) for x in np.quantile(bs, [0.025, 0.975])]}
    return out


def t_anom_region(path, region, base=(1959, 2022)):
    t = pd.read_parquet(path).set_index("date")[[region]]
    t.index = pd.to_datetime(t.index)
    b = t[(t.index.year >= base[0]) & (t.index.year <= base[1])]
    doy = b.groupby(b.index.dayofyear).mean().reindex(range(1, 367))
    doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    doy.index = range(1, 367)
    return (t - doy.reindex(t.index.dayofyear).values)[region]


def obs_verify(ta, events, onsets_all, free_years, win, p_by_q, rng, label):
    free_days = [d for d in pd.date_range(f"{free_years[0]}-01-01", f"{free_years[1]}-12-31")
                 if d.month in (11, 12, 1, 2, 3) and np.all(np.abs((onsets_all - d).days) > ZONE_SEP)]
    free = np.array([window(ta, d, win) for d in free_days]); free = free[np.isfinite(free)]
    vals = np.array([window(ta, o, win) for o in events]); ok = np.isfinite(vals)
    ev = [e for e, k in zip(events, ok) if k]; vals = vals[ok]
    out = {"n_free_windows": int(len(free))}
    for q, pr in p_by_q.items():
        thr = float(np.quantile(free, float(q)))
        r = verify(ev, vals < thr, pr, float(q), rng, label)
        r.pop("rev_curve"); r["threshold_K"] = round(thr, 3)
        out[str(q)] = r
    return out


def detect_era5():
    import ensemble_precursor as EP
    u = pd.read_parquet(ING / "era5_u10_60N_daily_cds.parquet")
    u = pd.Series(u["u10_60N"].values, index=pd.to_datetime(u["date"]))
    on = pd.DatetimeIndex(EP.detect_ssw(u.values, u.index))
    return u, on


def revision3():
    rng = np.random.default_rng(zlib.crc32(f"{NAME}|revision3".encode()))
    res = json.loads((RESULTS / "shift_rule_forecast.json").read_text())
    r3 = {"label": "revision-3 addendum, registered in commit add7ba9 before these outcomes were retrieved"}
    cat = load_catalogue("primary")
    wb2 = ING / "era5_t2m_regions_daily.parquet"
    full = ING / "era5_t2m_regions_daily_full.parquet"
    v1_ev = [o for o in cat if o <= pd.Timestamp("2023-01-10") - pd.Timedelta(days=T_WIN[1])]
    # SNAPSI p for every region, days 8-24, all q
    sp = snapsi_probs_gen(rng, T_WIN, R3_QS)
    r3["snapsi_days8_24"] = sp
    # R1 strict independence (NEURASIA, registered q)
    indep = [o for o in v1_ev if np.all(np.abs((pd.DatetimeIndex(SNAPSI_EVENTS) - o).days) > 3)]
    ta = t_anom_region(wb2, REG)
    r3["R1_strict_independence"] = obs_verify(ta, indep, cat, (1959, 2022), T_WIN,
                                              {str(q): sp["p"][REG][str(q)]["p"] for q in QS}, rng, "R1")
    r3["R1_strict_independence"]["dropped"] = [str(o.date()) for o in v1_ev if o not in indep]
    # R2/R3 severity x regions, days 8-24
    r3["R2_R3_days8_24"] = {}
    for r in REGIONS:
        tr = t_anom_region(wb2, r)
        r3["R2_R3_days8_24"][r] = obs_verify(tr, v1_ev, cat, (1959, 2022), T_WIN,
                                             {str(q): sp["p"][r][str(q)]["p"] for q in R3_QS}, rng, f"R2R3-{r}")
    # R4 weeks 3-6
    r3["R4_weeks"] = {}
    for w in R4_WINS:
        spw = snapsi_probs_gen(rng, w, QS)
        k = f"days{w[0]}_{w[1]}"
        r3["R4_weeks"][k] = {"snapsi": spw}
        if spw.get("status"):
            continue
        for r in REGIONS:
            tr = t_anom_region(full, r)
            r3["R4_weeks"][k][r] = obs_verify(tr, v1_ev, cat, (1959, 2022), w,
                                              {str(q): spw["p"][r][str(q)]["p"] for q in QS}, rng, f"R4-{k}-{r}")
    # R5 second predictor: CMIP6 annular-mode shift through the ERA5 event-free relation
    from scipy.stats import norm as _norm
    nam = pd.read_parquet(ING / "era5_nam_daily.parquet")
    nam.index = pd.to_datetime(nam.index)
    if np.corrcoef(nam["nam_1000"], nam["z1000_m"])[0, 1] > 0:
        raise ValueError("nam_1000 is not in the NAM sign (positive = low polar height)")
    free_days = [d for d in pd.date_range("1959-01-01", "2022-12-31")
                 if d.month in (11, 12, 1, 2, 3) and np.all(np.abs((cat - d).days) > ZONE_SEP)]
    fd = pd.DatetimeIndex(free_days)
    sd_free = float(nam.loc[nam.index.isin(fd), "nam_1000"].std())
    nb = nam["nam_1000"] / sd_free
    clim = nb[nb.index.isin(fd)].groupby(nb[nb.index.isin(fd)].index.dayofyear).mean()
    nanom = nb - clim.reindex(nb.index.dayofyear).values
    X = np.array([window(nanom, d, T_WIN) for d in free_days]); Y = np.array([window(ta, d, T_WIN) for d in free_days])
    ok = np.isfinite(X) & np.isfinite(Y)
    b, a_ = np.polyfit(X[ok], Y[ok], 1); se = float(np.std(Y[ok] - (a_ + b * X[ok]), ddof=2))
    Ni, sds = [], []
    for f in sorted(W.RAW.glob("*_zm.nc")):
        c = W.prepare_member(f)
        if c is None:
            continue
        # same units as the ERA5 regressor: the member's event-free November-March s.d.
        # (the index is standardised by its all-year s.d.; unit fix 2026-10-02 after review)
        am = c["am"]; msk = W.C6.influence_mask(am.index, c["on"])
        sdm = float(am[np.isin(am.index.month, (11, 12, 1, 2, 3)) & ~msk].std())
        Ni.append(W.C6.anom(c["am"], c["on"], c["cl"], T_WIN) / sdm); sds.append(sdm)
    Ni = np.concatenate(Ni); Ni = Ni[np.isfinite(Ni)]
    Yall = Y[np.isfinite(Y)]                      # the windows obs_verify uses for its thresholds
    p2 = {}
    for q in QS:
        thr = float(np.quantile(Yall, q))
        p2[str(q)] = round(float(np.mean(_norm.cdf((thr - a_ - b * Ni) / se))), 4)
    r3["R5_second_predictor"] = {"era5_free_regression": {"slope_K_per_sd": round(float(b), 4), "intercept": round(float(a_), 4),
                                                         "resid_sd": round(se, 4), "n_free": int(ok.sum())},
                                 "cmip6_n_events": int(len(Ni)), "cmip6_mean_nam_days8_24": round(float(Ni.mean()), 4),
                                 "cmip6_member_ndjfm_sd_range": [round(min(sds), 3), round(max(sds), 3)],
                                 "p_rule2": p2,
                                 "verification": obs_verify(ta, v1_ev, cat, (1959, 2022), T_WIN, p2, rng, "R5")}
    # R6 out of sample 1940-1958 (needs the ERA5 wind series)
    if not (ING / "era5_u10_60N_daily_cds.parquet").exists():
        r3["R6_out_of_sample_1940_1958"] = {"status": "pending: ERA5 u(10 hPa) series not yet retrieved"}
        res["revision3"] = r3
        (RESULTS / "shift_rule_forecast.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
        print(json.dumps({k: v for k, v in r3.items() if k != "snapsi_days8_24"}, default=str)[:8000])
        print("Saved -> shift_rule_forecast.json (revision3, R6 pending)")
        return
    u, on = detect_era5()
    early_on = on[(on >= EARLY[0]) & (on <= EARLY[1])]
    cat_match = {str(o.date()): (int(np.min(np.abs((on - o).days))) if len(on) else None) for o in cat}
    r6 = {"era5_onsets_1940_1958": [str(o.date()) for o in early_on],
          "detector_vs_catalogue": {"n_catalogue": int(len(cat)),
                                    "matched_within_3d": int(sum(v is not None and v <= 3 for v in cat_match.values())),
                                    "era5_onsets_1958_2024_not_in_catalogue": [str(o.date()) for o in on
                                        if pd.Timestamp("1958-01-01") <= o <= pd.Timestamp("2024-04-30")
                                        and np.min(np.abs((cat - o).days)) > 3]}}
    ea = ING / "era5_t2m_regions_daily_arco_early.parquet"
    te = t_anom_region(ea, REG, base=(1941, 1958))
    early_free_on = on[on <= pd.Timestamp("1959-06-01")]
    pr = {str(q): sp["p"][REG][str(q)]["p"] for q in QS}
    r6["primary"] = obs_verify(te, list(early_on), early_free_on, (1941, 1958), T_WIN, pr, rng, "R6")
    r6["from_1946"] = obs_verify(te, [o for o in early_on if o >= pd.Timestamp("1946-01-01")], early_free_on,
                                 (1941, 1958), T_WIN, pr, rng, "R6-1946")
    # pooled with V1 (counts and binomial under the rule / climatology)
    r6["pooled_with_V1"] = {}
    for q in QS:
        a1, b1 = res["verification"]["V1"][str(q)], r6["primary"][str(q)]
        k, n = a1["count"] + b1["count"], a1["n"] + b1["n"]
        r6["pooled_with_V1"][str(q)] = {"count": k, "n": n, "f": round(k / n, 4), "f_wilson95": wilson(k, n),
                                        "binom_p_under_rule": round(float(stats.binomtest(k, n, pr[str(q)]).pvalue), 4),
                                        "binom_p_under_climatology": round(float(stats.binomtest(k, n, q).pvalue), 4)}
    # secondary: negative polar-cap NAM proxy days 8-52
    pe = pd.read_parquet(ING / "era5_psl_cap_00utc_arco_early.parquet")
    ps = pd.Series(pe["psl_cap_N"].values, index=pd.to_datetime(pe["time"]).dt.floor("D"))
    b_ = ps[(ps.index.year >= 1941) & (ps.index.year <= 1958)]
    dclim = b_.groupby(b_.index.dayofyear).mean().reindex(range(1, 367))
    dclim = pd.concat([dclim.iloc[-15:], dclim, dclim.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    dclim.index = range(1, 367)
    pan = ps - dclim.reindex(ps.index.dayofyear).values          # > 0 = negative NAM
    yv = np.array([window(pan, o, AO_WIN) for o in early_on]); okk = np.isfinite(yv)
    fdays = [d for d in pd.date_range("1941-01-01", "1958-12-31")
             if d.month in (11, 12, 1, 2, 3) and np.all(np.abs((early_free_on - d).days) > ZONE_SEP)]
    fr = np.array([window(pan, d, AO_WIN) for d in fdays]); fr = fr[np.isfinite(fr)]
    base = float(np.mean(fr > 0))
    r6["secondary_negative_nam"] = verify([o for o, k in zip(early_on, okk) if k], yv[okk] > 0,
                                          res["model_p_negative_am_cmip6"]["p_negative"], base, rng, "R6-NAM")
    r6["secondary_negative_nam"].pop("rev_curve"); r6["secondary_negative_nam"]["climatological_base_rate"] = round(base, 4)
    r3["R6_out_of_sample_1940_1958"] = r6
    res["revision3"] = r3
    (RESULTS / "shift_rule_forecast.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(json.dumps({k: v for k, v in r3.items() if k != "snapsi_days8_24"}, default=str)[:6000])
    print("Saved -> shift_rule_forecast.json (revision3)")


def r6_diagnostics():
    """EXPLORATORY, added after the R6 result was seen (labelled so in the JSON):
    how the 1941-58 SSWs differ from those of 1959-2022."""
    res = json.loads((RESULTS / "shift_rule_forecast.json").read_text())
    u, on = detect_era5()
    early = [o for o in on if EARLY[0] <= o <= EARLY[1]]
    cat = load_catalogue("primary")
    later = [o for o in cat if pd.Timestamp("1959-01-01") <= o <= pd.Timestamp("2022-12-15")]
    umin = lambda o: float(u[o:o + pd.Timedelta(days=20)].min())
    ue, ul = np.array([umin(o) for o in early]), np.array([umin(o) for o in later])
    te = t_anom_region(ING / "era5_t2m_regions_daily_arco_early.parquet", REG, base=(1941, 1958))
    tl = t_anom_region(ING / "era5_t2m_regions_daily.parquet", REG)
    Te = np.array([window(te, o, T_WIN) for o in early]); Tl = np.array([window(tl, o, T_WIN) for o in later])
    pe = pd.read_parquet(ING / "era5_psl_cap_00utc_arco_early.parquet")
    ps = pd.Series(pe["psl_cap_N"].values, index=pd.to_datetime(pe["time"]).dt.floor("D"))
    b_ = ps[(ps.index.year >= 1941) & (ps.index.year <= 1958)]
    c = b_.groupby(b_.index.dayofyear).mean().reindex(range(1, 367))
    c = pd.concat([c.iloc[-15:], c, c.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]; c.index = range(1, 367)
    pan = ps - c.reindex(ps.index.dayofyear).values
    sd = float(pan[pan.index.month.isin([11, 12, 1, 2, 3])].std())
    nam_e = np.array([-window(pan, o, (8, 25)) / sd for o in early])
    free = [d for d in pd.date_range("1941-01-01", "1958-12-31") if d.month in (11, 12, 1, 2, 3)
            and np.all(np.abs((on[on <= pd.Timestamp("1959-06-01")] - d).days) > ZONE_SEP)]
    nam_f = np.array([-window(pan, d, (8, 25)) / sd for d in free]); nam_f = nam_f[np.isfinite(nam_f)]
    vt = json.loads((ROOT / "results" / "current" / "5_mechanism" / "vortex_threshold_continuity.json").read_text())
    dose = vt["era5_1940_2025"]["all"]["Y3_T_NEURASIA"]["E1_E2"]["dose_slope_per_10ms"]
    dmed = float(np.median(ue) - np.median(ul))
    out = {"label": "exploratory, added after the R6 result was seen",
           "median_umin_early": round(float(np.median(ue)), 2), "median_umin_1959_2022": round(float(np.median(ul)), 2),
           "share_umin_above_minus4_early": round(float(np.mean(ue > -4)), 3),
           "share_umin_above_minus4_1959_2022": round(float(np.mean(ul > -4)), 3),
           "mean_T_anomaly_early_K": round(float(np.nanmean(Te)), 3), "mean_T_anomaly_1959_2022_K": round(float(np.nanmean(Tl)), 3),
           "nam_proxy_days8_25_early_minus_free_sd": round(float(np.nanmean(nam_e) - nam_f.mean()), 3),
           "era5_dose_K_per_10ms": dose,
           "cooling_difference_explained_by_dose_K": round(float(-dose * dmed / 10.0), 3),
           "n_early": len(early), "n_1959_2022": len(later)}
    res["revision3"]["R6_out_of_sample_1940_1958"]["exploratory_diagnostics"] = out
    (RESULTS / "shift_rule_forecast.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(json.dumps(out))


if __name__ == "__main__":
    if "--r6-diagnostics" in sys.argv:
        r6_diagnostics()
        sys.exit(0)
    if "--revision3" in sys.argv:
        revision3()
    else:
        main()
