#!/usr/bin/env python3
"""
state_forecast.py
=================
AN ISSUE-TIME FORECAST CONDITIONED ON THE CONTINUOUS STATE, EVALUATED AGAINST
CLIMATOLOGY, THE FIXED SSW RULES AND OPERATIONAL ENSEMBLES.

Plan approved 2026-10-02 (user chose "Add a new forecast" after the external audit,
section 7). Committed BEFORE it was run. All observed SSW outcomes used here have
been examined before (V, R6, T-g); the evaluation is therefore RETROSPECTIVE
cross-validation with whole-winter folds, not prospective evidence. The models
were specified now, after those results were known.

QUESTION
  After an SSW, does a probabilistic forecast that uses the continuous state known
  at issue time (stratospheric wind and its fall, the surface annular mode and
  regional temperature) predict regional cold better than (a) climatology, (b) the
  fixed SSW rules of the paper, and (c) raw and calibrated operational ensembles?
  And does knowing that the state is an "SSW" add anything beyond that state?

TARGET (primary, fixed): northern-Eurasian (50-65N, 10-130E) 2 m temperature
  anomaly, mean over days 8-24 after issue; CRPS. Secondary: P(below the
  event-free 10th percentile), Brier score.
OBSERVATIONS: ERA5 daily means 1940-2026 (ARCO-ERA5 1940-58 and 2023-26,
  WeatherBench 2 1959-2023, spliced as in vortex_threshold_continuity.era5_outcomes);
  u(10 hPa, 60N) daily mean from the CDS. Anomalies: 31-day-smoothed day-of-year
  climatology and a linear trend, BOTH FITTED ON THE TRAINING WINTERS OF EACH FOLD.
CASES
  E1 (primary): every ERA5 SSW onset (project CP07 detector on the CDS wind),
    November 1940 - March 2026; issue date = onset day.
  E2 (secondary): the catalogued SSWs of 1998-2021 with S2S reforecast starts 0-7
    days after onset (ten systems); issue date = start date; target window days
    8-24 after onset.
INFORMATION AT ISSUE (day d; nothing after d)
  u10(d); u10(d) - u10(d-10); min u10 over d-10..d; polar-cap NAM proxy mean over
  d-5..d; regional temperature anomaly mean over d-5..d; sin/cos day of year.
FORECASTS (all fitted within the fold)
  F0 climatology: event-free training days within +-15 days of the issue day of
     year (> 135 d from every onset); empirical distribution.
  F1 fixed SNAPSI rule: F0 displaced by -0.826 s.d. of F0 (the Gaussian shift that
     gives the SNAPSI probability 0.324 below the 10th percentile).
  F2 constant CMIP6 x ERA5 rule (as revision-3 R5): mixture over the 1,517 CMIP6
     events of N(a + b N_i, s_e), with a, b, s_e from the training event-free days.
  F3 STATE MODEL: Gaussian, mean linear in the issue-time predictors, constant
     residual s.d.; trained on ALL training winter days (Nov-Mar), SSW or not. It
     is told nothing about SSWs.
  F4 STATE + SSW: F3 plus an indicator of an SSW onset on day d.
  F3* ORACLE (reported separately, not a forecast): F3 plus the realised minimum
     u10 over d..d+20.
  E2 only: RAW ensemble (members' window-mean anomaly, each system against its
     leave-one-year-out reforecast climatology) and CALIBRATED ensemble (mean bias
     and spread factor from that system's starts in the training winters, same
     lead window); scored per system and averaged per event.
SCORES AND TESTS
  Paired CRPS differences by event; intervals from 2,000 resamples of winters;
  skill = 1 - CRPS_model / CRPS_reference. Reported by era (onsets before 1979,
  from 1979) and the leave-one-winter-out range.
DECISION RULE (fixed now)
  A forecasting headline is kept only if F3's CRPS skill over BOTH F0 and F1 is
  >= 0.05 with its 95% interval excluding zero, positive in both eras, and
  positive with any single winter left out. Otherwise the result is reported as a
  limitation. "SSW adds beyond the state" requires F4 to beat F3 with an interval
  excluding zero.

Output: results/current/6_predictability/state_forecast.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
NAME = "state_forecast"
SEED = zlib.crc32(NAME.encode()) % (2 ** 32)
N_BOOT = 2000
WIN = (8, 24)
SNAPSI_SHIFT = -0.826
ZONE = 135


def gcrps(mu, sd, y):
    z = (y - mu) / sd
    return float(sd * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi)))


def ecrps(x, y):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    xs = np.sort(x); n = len(xs)
    # O(n log n) form of mean|x-y| - 0.5 mean|x-x'|
    i = np.arange(1, n + 1)
    return float(np.mean(np.abs(xs - y)) - np.sum((2 * i - n - 1) * xs) / n ** 2)


def mix_crps_sample(mus, sd, y, rng, n=4000):
    s = rng.choice(mus, n) + rng.normal(0, sd, n)
    return ecrps(s, y), s


def winter(t):
    return t.year if t.month <= 6 else t.year + 1


def raw_series():
    """Daily raw polar-cap msl (00 UTC) and N-Eurasian T (daily mean), 1940-2026,
    spliced exactly as vortex_threshold_continuity.era5_outcomes (no anomalies)."""
    import s2s_heldout_test as H
    st, _ = H.splice()
    op, ot = st["psl_cap_N"]["offset"], st["NEURASIA"]["offset"]

    def psl(f):
        d = pd.read_parquet(ING / f); t = pd.to_datetime(d["time"]); d = d[t.dt.hour == 0]
        return pd.Series(d["psl_cap_N"].values, index=pd.to_datetime(d["time"]).dt.floor("D"))
    wb = psl("era5_psl_cap_6h_full.parquet")
    ar = pd.concat([psl("era5_psl_cap_00utc_arco.parquet"), psl("era5_psl_cap_00utc_arco_late.parquet")]) - op
    ea = psl("era5_psl_cap_00utc_arco_early.parquet") - op
    p = pd.concat([ea[ea.index < wb.index.min()], wb, ar[ar.index > wb.index.max()]]).sort_index()

    def t2(f):
        d = pd.read_parquet(ING / f)
        return pd.Series(d["NEURASIA"].values, index=pd.to_datetime(d["date"]))
    tw = t2("era5_t2m_regions_daily_full.parquet")
    tsp = t2("era5_t2m_regions_daily_spliced.parquet")
    tl = t2("era5_t2m_regions_daily_arco_late.parquet") - ot
    te = t2("era5_t2m_regions_daily_arco_early.parquet") - ot
    t = pd.concat([te[te.index < tw.index.min()], tw, tsp[tsp.index > tw.index.max()], tl[tl.index > tsp.index.max()]]).sort_index()
    u = pd.read_parquet(ING / "era5_u10_60N_daily_cds.parquet"); u = pd.Series(u["u10_60N"].values, index=pd.to_datetime(u["date"]))
    idx = pd.date_range("1940-01-01", "2026-05-31", freq="D")
    return (p[~p.index.duplicated()].reindex(idx), t[~t.index.duplicated()].reindex(idx), u.reindex(idx))


def anomalies(x, train_mask, trend=True, sd_months=(11, 12, 1, 2, 3)):
    """Smoothed day-of-year climatology (+ linear trend) fitted on training days."""
    tr = x[train_mask & x.notna()]
    c = tr.groupby(tr.index.dayofyear).mean().reindex(range(1, 367))
    c = pd.concat([c.iloc[-15:], c, c.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    c.index = range(1, 367)
    a = x - c.reindex(x.index.dayofyear).values
    if trend:
        yr = (a.index - pd.Timestamp("1980-01-01")).days.values / 3652.5
        m = (train_mask & a.notna()).values
        b = np.polyfit(yr[m], a.values[m], 1)
        a = a - np.polyval(b, yr)
    return a


def rolling_mean(s, a, b):
    """Mean of s over d+a..d+b for every d (NaN unless >= 90% present)."""
    v = s.values; n = len(v); out = np.full(n, np.nan)
    ok = np.isfinite(v); cs = np.concatenate([[0], np.cumsum(np.where(ok, v, 0))]); cn = np.concatenate([[0], np.cumsum(ok)])
    L = b - a + 1
    for d in range(n):
        lo, hi = d + a, d + b + 1
        if lo < 0 or hi > n:
            continue
        k = cn[hi] - cn[lo]
        if k >= 0.9 * L:
            out[d] = (cs[hi] - cs[lo]) / k
    return pd.Series(out, index=s.index)


def features(p, t, u, train_mask, a=WIN[0], b=WIN[1]):
    pa = -anomalies(p, train_mask, trend=False)
    pa = pa / pa[train_mask & pa.index.month.isin([11, 12, 1, 2, 3])].std()
    ta = anomalies(t, train_mask, trend=True)
    X = pd.DataFrame({"u0": u, "du": u - u.shift(10), "umin": u.rolling(11).min(),
                      "nam5": rolling_mean(pa, -5, 0), "t5": rolling_mean(ta, -5, 0),
                      "s": np.sin(2 * np.pi * u.index.dayofyear / 365.25), "c": np.cos(2 * np.pi * u.index.dayofyear / 365.25)})
    y = rolling_mean(ta, a, b)
    nam_w = rolling_mean(pa, a, b)
    oracle = u[::-1].rolling(21, min_periods=18).min()[::-1]       # min over d..d+20
    return X, y, nam_w, oracle


SYSTEMS = ["ecmwf", "eccc", "cma", "hmcr", "kma", "cnrm", "jma", "cnr_isac", "ncep", "cptec"]


def load_t(c):
    tg = "ecmf" if c == "ecmwf" else c
    d = pd.concat([pd.read_parquet(ING / f"s2s_{tg}_t2m_regions_short.parquet"), pd.read_parquet(ING / f"s2s_{tg}_t2m_regions.parquet")])
    d["init"] = pd.to_datetime(d["init"]); d = d.drop_duplicates(["init", "member", "lead_day"])
    w = d.pivot_table(index=["init", "member"], columns="lead_day", values="NEURASIA").reindex(columns=range(1, 34))
    yr = pd.to_datetime(w.index.get_level_values("init")).year.values
    md = pd.to_datetime(w.index.get_level_values("init")).strftime("%m-%d").values
    an = w.copy()
    for m in np.unique(md):                                    # leave-one-year-out model climatology
        sel = md == m
        sub = w[sel]; ys = yr[sel]
        tot = sub.groupby(ys).sum(min_count=1); cnt = sub.groupby(ys).count()
        T, N = tot.sum(), cnt.sum()
        for y_ in np.unique(ys):
            clim = (T - tot.loc[y_]) / (N - cnt.loc[y_])
            idx = np.flatnonzero(sel)[ys == y_]
            an.iloc[idx] = w.iloc[idx].values - clim.values
    return an


def e2(p, t, u, rng, skill, far):
    """Operational comparison: catalogued SSWs 1998-2021, starts 0-7 d after onset."""
    from build_catalogue import load_catalogue
    cat = load_catalogue("primary")
    ev = [o for o in cat if pd.Timestamp("1998-11-01") <= o <= pd.Timestamp("2021-04-30")]
    A = {c: load_t(c) for c in SYSTEMS}
    win = np.array([winter(d) for d in u.index]); ndjfm = u.index.month.isin([11, 12, 1, 2, 3])
    cols = ["u0", "du", "umin", "nam5", "t5", "s", "c"]
    rows = []
    for w in sorted({winter(o) for o in ev}):
        train = pd.Series(win != w, index=u.index)
        Xf, _, _, _ = features(p, t, u, train)
        ta = anomalies(t, train, trend=True)
        cache = {}
        for o in [o for o in ev if winter(o) == w]:
            yo = rolling_mean(ta, 8, 24).loc[o]
            if not np.isfinite(yo):
                continue
            rec = {"onset": str(o.date()), "winter": w, "raw": [], "cal": [], "F3": [], "F0": [], "lag": []}
            for c, an in A.items():
                inits = pd.to_datetime(an.index.get_level_values("init")).unique()
                for s0 in [s_ for s_ in inits if 0 <= (s_ - o).days <= 7]:
                    k = (s0 - o).days; a_, b_ = 8 - k, 24 - k
                    mem = an.loc[s0].loc[:, a_:b_].mean(axis=1).values
                    mem = mem[np.isfinite(mem)]
                    if len(mem) < 3 or not Xf[cols].loc[s0].notna().all():
                        continue                               # paired: every method scored on this start or none
                    rec["raw"].append(ecrps(mem, yo)); rec["lag"].append(k)
                    key = (c, a_, b_)
                    if key not in cache:                       # calibration on training-winter starts
                        ii = pd.to_datetime(an.index.get_level_values("init"))
                        wtr = np.array([winter(x) for x in ii]) != w
                        em = an[wtr].loc[:, a_:b_].mean(axis=1).groupby(level="init").agg(["mean", "var"])
                        yo_tr = rolling_mean(ta, a_, b_).reindex(em.index).values
                        okc = np.isfinite(yo_tr) & np.isfinite(em["mean"].values)
                        bias = float(np.mean(yo_tr[okc] - em["mean"].values[okc]))
                        mse = float(np.mean((yo_tr[okc] - em["mean"].values[okc] - bias) ** 2))
                        phi = float(np.sqrt(mse / np.nanmean(em["var"].values[okc])))
                        cache[key] = (bias, phi)
                    bias, phi = cache[key]
                    rec["cal"].append(gcrps(mem.mean() + bias, max(phi * mem.std(ddof=1), 1e-3), yo))
                    kk = ("F3", a_, b_)
                    if kk not in cache:
                        _, yk, _, _ = features(p, t, u, train, a_, b_)
                        trd = train.values & ndjfm & Xf[cols].notna().all(1).values & yk.notna().values
                        Ak = np.column_stack([np.ones(trd.sum()), Xf[cols].values[trd]])
                        bk = np.linalg.lstsq(Ak, yk.values[trd], rcond=None)[0]
                        cache[kk] = (bk, float(np.std(yk.values[trd] - Ak @ bk, ddof=Ak.shape[1])), yk)
                    bk, sk, yk = cache[kk]
                    # the state model is REISSUED at this start date s0 (ERA5 state up to s0) and verified on the
                    # same window, days 8-24 after onset = days a_..b_ after s0: same issue date, valid time, case
                    xi = Xf[cols].loc[s0]
                    rec["F3"].append(gcrps(float(np.r_[1.0, xi.values] @ bk), sk, yo))
                    dd = np.abs(u.index.dayofyear.values - s0.dayofyear); dd = np.minimum(dd, 366 - dd)
                    # event-free, as the design and E1 (fix 2026-10-02: the first run omitted far)
                    clim = yk.values[train.values & ndjfm & far & (dd <= 15) & yk.notna().values]
                    rec["F0"].append(ecrps(clim, yo))
            if rec["raw"]:
                rows.append({"onset": rec["onset"], "winter": w, **{k: float(np.mean(rec[k])) for k in ("raw", "cal", "F3", "F0") if rec[k]},
                             "n_forecasts": len(rec["raw"]), "mean_issue_lag_days": float(np.mean(rec["lag"])),
                             "n_scored_per_method": {k: len(rec[k]) for k in ("raw", "cal", "F3", "F0")}})
    d = pd.DataFrame(rows)
    out = {"pairing": "each start scored by all four methods or none; F3 and F0 reissued at the start date, same valid window",
           "n_events": int(len(d)), "mean_crps": {k: round(float(d[k].mean()), 4) for k in ("raw", "cal", "F3", "F0")},
           "events": d.round(4).to_dict(orient="records")}
    for a_, b_ in (("F3", "raw"), ("F3", "cal"), ("cal", "raw"), ("F3", "F0"), ("cal", "F0")):
        out[f"{a_}_vs_{b_}"] = skill(a_, b_, d.dropna(subset=[a_, b_]))
    return out


def main():
    rng = np.random.default_rng(SEED)
    import ensemble_precursor as EP
    from build_catalogue import load_catalogue
    import within_model_check as Wm
    p, t, u = raw_series()
    on = pd.DatetimeIndex(EP.detect_ssw(u.dropna().values, u.dropna().index))
    E1 = [o for o in on if pd.Timestamp("1940-11-01") <= o <= pd.Timestamp("2026-03-31")]
    win = np.array([winter(d) for d in u.index])
    ndjfm = u.index.month.isin([11, 12, 1, 2, 3])
    far = np.ones(len(u), bool)
    for o in on:
        far &= np.abs((u.index - o).days) > ZONE
    # CMIP6 N_i for F2 (days 8-24 annular-mode anomalies)
    Ni = []
    for f in sorted(Wm.RAW.glob("*_zm.nc")):
        c = Wm.prepare_member(f)
        if c is not None:
            # same units as the observed regressor (event-free Nov-Mar s.d.; fix 2026-10-02)
            am = c["am"]; msk = Wm.C6.influence_mask(am.index, c["on"])
            sdm = float(am[np.isin(am.index.month, (11, 12, 1, 2, 3)) & ~msk].std())
            Ni.append(Wm.C6.anom(c["am"], c["on"], c["cl"], WIN) / sdm)
    Ni = np.concatenate(Ni); Ni = Ni[np.isfinite(Ni)]
    cols = ["u0", "du", "umin", "nam5", "t5", "s", "c"]
    rows = []
    for w in sorted({winter(o) for o in E1}):
        train = pd.Series(win != w, index=u.index)
        X, y, namw, orc = features(p, t, u, train)
        trd = train.values & ndjfm & X[cols].notna().all(1).values & y.notna().values
        A = np.column_stack([np.ones(trd.sum()), X[cols].values[trd]])
        b3 = np.linalg.lstsq(A, y.values[trd], rcond=None)[0]; s3 = float(np.std(y.values[trd] - A @ b3, ddof=A.shape[1]))
        ssw_day = np.isin(u.index, on).astype(float)
        A4 = np.column_stack([A, ssw_day[trd]])
        b4 = np.linalg.lstsq(A4, y.values[trd], rcond=None)[0]; s4 = float(np.std(y.values[trd] - A4 @ b4, ddof=A4.shape[1]))
        Ao = np.column_stack([A, orc.values[trd]]); okk = np.isfinite(Ao).all(1)
        bo = np.linalg.lstsq(Ao[okk], y.values[trd][okk], rcond=None)[0]; so = float(np.std(y.values[trd][okk] - Ao[okk] @ bo, ddof=Ao.shape[1]))
        fr = train.values & ndjfm & far & y.notna().values
        fm = fr & namw.notna().values
        b2, a2 = np.polyfit(namw.values[fm], y.values[fm], 1); s2 = float(np.std(y.values[fm] - (a2 + b2 * namw.values[fm]), ddof=2))
        for o in [o for o in E1 if winter(o) == w]:
            i = u.index.get_loc(o)
            if not (np.isfinite(y.iloc[i]) and X[cols].iloc[i].notna().all()):
                continue
            yo = float(y.iloc[i]); doy = o.dayofyear
            dd = np.abs(u.index.dayofyear.values - doy); dd = np.minimum(dd, 366 - dd)
            clim = y.values[fr & (dd <= 15)]
            thr = float(np.quantile(clim, 0.10))
            x = np.r_[1.0, X[cols].iloc[i].values]
            mu3, mu4 = float(x @ b3), float(np.r_[x, 1.0] @ b4)
            r = {"onset": str(o.date()), "winter": w, "y": yo, "cold": int(yo < thr)}
            r["F0"] = ecrps(clim, yo); r["F1"] = ecrps(clim + SNAPSI_SHIFT * clim.std(), yo)
            c2, smp = mix_crps_sample(a2 + b2 * Ni, s2, yo, rng); r["F2"] = c2
            r["F3"] = gcrps(mu3, s3, yo); r["F4"] = gcrps(mu4, s4, yo)
            if np.isfinite(orc.iloc[i]):
                r["F3_oracle"] = gcrps(float(np.r_[x, orc.iloc[i]] @ bo), so, yo)
            r["P0"] = float(np.mean(clim < thr)); r["P1"] = float(np.mean(clim + SNAPSI_SHIFT * clim.std() < thr))
            r["P2"] = float(np.mean(smp < thr)); r["P3"] = float(norm.cdf((thr - mu3) / s3)); r["P4"] = float(norm.cdf((thr - mu4) / s4))
            r["mu3"], r["sd3"] = mu3, s3
            rows.append(r)
        print(f"winter {w}: {sum(1 for r in rows if r['winter'] == w)} events", flush=True)
    df = pd.DataFrame(rows)
    res = {"plan_approved": "2026-10-02", "seed": SEED, "n_boot": N_BOOT, "label": "retrospective whole-winter cross-validation",
           "n_events": int(len(df)), "n_winters": int(df.winter.nunique()), "coef_note": "per fold; see events",
           "events": df.round(4).to_dict(orient="records")}

    def skill(a, b, d):
        uw = d.winter.unique()
        est = 1 - d[a].mean() / d[b].mean()
        bs = []
        for _ in range(N_BOOT):
            s_ = pd.concat([d[d.winter == x] for x in rng.choice(uw, len(uw))])
            bs.append(1 - s_[a].mean() / s_[b].mean())
        lowo = [1 - d[d.winter != x][a].mean() / d[d.winter != x][b].mean() for x in uw]
        return {"skill": round(float(est), 4), "ci95": [round(float(q), 4) for q in np.percentile(bs, [2.5, 97.5])],
                "lowo_min": round(float(np.min(lowo)), 4), "lowo_max": round(float(np.max(lowo)), 4)}
    comp = {}
    for a, b in (("F3", "F0"), ("F3", "F1"), ("F3", "F2"), ("F1", "F0"), ("F2", "F0"), ("F4", "F3"), ("F3_oracle", "F3")):
        d = df.dropna(subset=[a, b])
        comp[f"{a}_vs_{b}"] = {"all": skill(a, b, d),
                               "before_1979": skill(a, b, d[d.winter < 1979]) if (d.winter < 1979).sum() >= 5 else None,
                               "from_1979": skill(a, b, d[d.winter >= 1979]) if (d.winter >= 1979).sum() >= 5 else None}
    res["E1_crps_skill"] = comp
    res["E1_mean_crps"] = {k: round(float(df[k].mean()), 4) for k in ("F0", "F1", "F2", "F3", "F4") }
    res["E1_brier"] = {k: round(float(np.mean((df[f"P{k[1]}"] - df["cold"]) ** 2)), 4) for k in ("F0", "F1", "F2", "F3", "F4")}
    res["E1_cold_count"] = {"observed": int(df.cold.sum()), "n": int(len(df)),
                            "expected": {k: round(float(df[f"P{k[1]}"].sum()), 2) for k in ("F0", "F1", "F2", "F3", "F4")}}
    a3, a1, a0 = res["E1_crps_skill"]["F3_vs_F1"], res["E1_crps_skill"]["F3_vs_F0"], None
    ok = all(c is not None and c["skill"] >= 0.05 and c["ci95"][0] > 0 for c in (a3["all"], res["E1_crps_skill"]["F3_vs_F0"]["all"])) \
        and all(x is not None and x["skill"] > 0 for x in (a3["before_1979"], a3["from_1979"], a1["before_1979"], a1["from_1979"])) \
        and a3["all"]["lowo_min"] > 0 and a1["all"]["lowo_min"] > 0
    res["decision"] = "forecasting headline kept" if ok else "forecasting reported as a limitation (decision rule not met)"
    f4 = res["E1_crps_skill"]["F4_vs_F3"]["all"]
    res["ssw_adds_beyond_state"] = bool(f4["ci95"][0] > 0)
    res["E2"] = e2(p, t, u, rng, skill, far)
    for k in ("E1_mean_crps", "E1_crps_skill", "E1_brier", "E1_cold_count", "decision", "ssw_adds_beyond_state", "E2"):
        print(k, json.dumps(res[k]), flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "state_forecast.json").write_text(json.dumps(res, indent=2, default=float), encoding="utf8", newline="\n")
    print("Saved -> state_forecast.json")


if __name__ == "__main__":
    main()
