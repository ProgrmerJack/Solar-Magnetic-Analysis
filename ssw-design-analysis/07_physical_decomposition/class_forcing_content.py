#!/usr/bin/env python3
"""
class_forcing_content.py
========================
HOW MUCH OF A CLASS CONTRAST IS A DIFFERENCE IN FORCING, WHEN FORCING DIFFERS
BETWEEN THE UNITS THAT ARE CLASSIFIED?

Registered 2026-10-09, committed BEFORE any statistic below was computed. The
inputs (SNAPSI nudged/control and free members, the vortex-weakening unit tables,
IBTrACS/ERSST storms, HadGEM3 DAMIP boxes, CMIP6 zonal means) were analysed before
for other statistics; none of the quantities defined here had been computed.

WHY
  Within one nudged ensemble forcing is shared, so its class contrast can only be
  selection on internal variability. The field applies the label ACROSS events,
  whose forcing differs, and reads "downward events were more strongly forced".
  Write a unit's outcome X = F + E (F its forced shift, E internal). Classes cut
  from X differ in mean outcome by D and in mean forcing by D_F. If F and E are
  jointly Gaussian, E[F | X] is linear in X and
        D_F = h2 * D,   h2 = Var(F) / Var(X),
  the breeder's equation of quantitative genetics (response = heritability x
  selection differential; Lush 1937; Falconer 1965): a class contrast reveals a
  forcing difference only in proportion to the share of outcome variance the
  forcing explains. Under the paper's law (one location-scale population per
  forcing level) the share R = D_F / D can be computed exactly for any level
  structure (R_LS below); h2 is its Gaussian value (binary forcing, shift 1.36 s.d.,
  half the units forced: R 0.335 against h2 0.316; computed 2026-10-09 by
  quadrature, before any data statistic).

ESTIMANDS (every set)
  D     mean outcome of class minus mean outcome of the rest
  D_F   the same difference of the units' forced shift F
  R_obs = D_F / D (observed share);  R_h2 = Var(F) / Var(X);
  R_LS  = the share when each forced level is the reference level mapped
          x -> m_f + k (x - m_r) (the law), pooled with the observed level shares
          and classified identically (sets with known levels only).

P1 PRIMARY, known binary forcing (SNAPSI zonal-mean nudging). In each of the 36
  NH centre x initialisation pairs (NRL excluded by the corruption guard, as
  throughout), nudged and control members are pooled; outcome, standardisation and
  classification (conditions 1-2 on the 6-hourly NAM proxy, days 8-25) exactly as
  in snapsi_selection_test / snapsi_free_natural_experiment; F = Delta (nudged
  mean minus control mean) for nudged members, 0 for control. Statistic: mean over
  pairs of R_obs - R_LS; interval from 10,000 resamples of centres.
P2 PRIMARY, continuous forcing across real events (observations). The 114
  vortex-weakening episodes of NCEP-NCAR 1958-2025 (unit table of
  vortex_threshold_continuity.json): outcome the CPC AO over days 8-52, class the
  downward label (conditions 1-2 on the AO), forcing proxy the episode's minimum
  10 hPa, 60N wind u_min. F_hat = b (u_min - mean), b the least-squares slope of
  the AO on u_min. Statistic R_obs - R_h2 (R_h2 = r^2); interval from 10,000
  resamples of winters. Also reported: R_obs itself.
READINGS (fixed now; tolerance 0.05, a twentieth of the contrast)
  interval of the statistic inside +-0.05  -> "the forced share of the contrast is
     what the variance share predicts (the rule holds)"
  interval excluding 0 with |mean| > 0.05   -> "the label selects forcing beyond the
     variance share (the rule fails)"
  otherwise                                 -> inconclusive.
  P2 headline reading on R_obs: interval entirely below 0.5 -> "most of the observed
  downward-non-downward contrast is not a difference in vortex dose"; entirely
  above 0.5 -> "most of it is"; otherwise undetermined.
  Holm over P1 and P2 for any p value reported (two-sided bootstrap p of 0).
VALIDATION GATE (run with --validate, before P1/P2 are computed; results stored)
  P1: 200 synthetic data sets with the real pair layout, nudged members drawn
      with replacement from the pair's control members and shifted by the pair's
      Delta (the law exactly); 1,000 centre resamples each. Gate: one-sided 95%
      upper bound of the rate at which the interval excludes 0 below 7%.
      Power (reported, no gate): nudged = resampled control + B * Delta / 0.5,
      B ~ Bernoulli(0.5) ("couples or does not", same mean shift).
  P2: 1,000 synthetic data sets keeping the real AO and labels, with
      u_min* = mean + beta (AO - mean) + Gaussian noise matched to the real
      regression of u_min on the AO (E[u_min | AO] linear, so the rule holds
      exactly); 1,000 winter resamples each. Gate as P1.
  If a gate fails, that primary is not interpreted; the failure is reported.
SECONDARY (each reported with its own interval, not in the Holm family)
  S1 SNAPSI free natural experiment: members that develop the SSW vs not
     (snapsi_free_natural_experiment classification), R_obs - R_LS.
  S2 ERA5 1940-2025 episodes (144; unit table era5_1940_2025): polar-cap NAM
     proxy, label on the proxy, u_min; as P2.
  S3 tropical cyclones (tc_class_law inputs): levels = global genesis-SST
     terciles; classes major (LMI >= 96 kt, outcome LMI) and RI (24-h increase
     >= 30 kt, outcome that increase); F = the storm's level mean; R_obs - R_LS,
     2,000 season resamples. Repeated with basin-relative terciles (cuts within
     each basin, SA excluded), which removes basin composition from the levels.
  S4 human forcing, HadGEM3-GC31-LL historical vs hist-nat 1985-2014 (the
     damip_class_law targets and classes): per year pool both arms, F = Z * shift;
     R_obs - R_LS averaged over years; 10,000 year resamples.
  S5 CMIP6 1,517 events (response_shape_dose inputs): outcome the days 8-52
     annular mode, class outcome < 0 (condition 1; daily series not used), F_hat
     the least-squares fit on the realised 100 hPa wind dose (DOSE columns) with
     member fixed effects; R_obs - R_h2; 2,000 member resamples. The dose is
     observed after onset and partly reflects the troposphere.
  S6 derived, no test: R_max = s2 / (1 + s2) for the forced between-event
     variance s2 of CMIP6 (point and upper bounds) and observations (upper
     bounds) from the forced-variance-ceiling result, and for the ICON
     event-conditioned variance (0.41-0.43 between, ~0.8 of within); also at the
     strongest damping the experiment allows (upper s_f 0.42 CMIP6, 0.78 obs).
  T0 theory, no data: the largest R / h2 over two-point forced distributions
     (two forced types, weights 0.01-0.99, either sign) of the same variance, for
     class shares 0.5 and 0.65 -- how far the Gaussian rule can understate R
     (computed 2026-10-09 before registration: 1.06, 1.09, 1.16 for forced
     variances 0.019, 0.117, 0.46).
AMENDMENT 2026-10-09, before any statistic or validation was run: the P1 gate drew
  the interval's two limits from two separate bootstraps; it now uses one.
OUTPUT results/current/8_experiment/class_forcing_content.json
"""
import json
import sys
import warnings
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
OUT = RESULTS / "class_forcing_content.json"
sys.path.insert(0, str(HERE))
import snapsi_selection_test as L                   # noqa: E402
import snapsi_free_natural_experiment as NE         # noqa: E402

warnings.filterwarnings("ignore", message=".*(variant|p-value (capped|floored)).*")
NAME = "class_forcing_content"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_BOOT, TOL, GATE = 10000, 0.05, 0.07
CONT = ROOT / "results" / "current" / "5_mechanism" / "vortex_threshold_continuity.json"
TCS = HERE.parents[0] / "03_data_ingestion" / "tc_storms.parquet"
RI_KT, MAJOR_KT = 30.0, 96.0


# ------------------------------------------------------------------ core
def share(x, c, F, w=None):
    """D, D_F, R = D_F / D for outcome x, class mask c, forcing F, optional weights."""
    w = np.ones(len(x)) if w is None else w
    if c.sum() < 3 or (~c).sum() < 3:
        return np.nan, np.nan, np.nan
    mx = lambda v, m: float(np.sum(w[m] * v[m]) / np.sum(w[m]))
    D = mx(x, c) - mx(x, ~c); DF = mx(F, c) - mx(F, ~c)
    return D, DF, DF / D


def h2(x, F):
    return float(np.var(F) / np.var(x))


def pair_share(S, N):
    """S, N: lists of 6-hourly member series (forced, reference). Pooled shares."""
    ms, dws = NE.classify(S); mn, dwn = NE.classify(N)
    m_s, m_n = float(ms.mean()), float(mn.mean()); k = float(ms.std(ddof=1) / mn.std(ddof=1))
    delta = m_s - m_n
    x = np.concatenate([ms, mn]); c = np.concatenate([dws, dwn])
    F = np.concatenate([np.full(len(ms), delta), np.zeros(len(mn))])
    D, DF, R = share(x, c, F)
    mapped = [m_s + k * (s - m_n) for s in N]
    mm, dwm = NE.classify(mapped)
    p = len(S) / (len(S) + len(N))
    xl = np.concatenate([mm, mn]); cl = np.concatenate([dwm, dwn])
    Fl = np.concatenate([np.full(len(mm), delta), np.zeros(len(mn))])
    wl = np.concatenate([np.full(len(mm), p / len(mm)), np.full(len(mn), (1 - p) / len(mn))])
    _, _, R_LS = share(xl, cl, Fl, wl)
    return {"D": D, "D_F": DF, "R_obs": R, "R_LS": R_LS, "R_h2": h2(x, F), "delta": delta, "sd_ratio": k,
            "p_forced": p, "stat": R - R_LS}


def centre_summary(df, col, rng, n_boot=N_BOOT):
    d = df.dropna(subset=[col])
    cs = d.centre.unique(); by = {c: d.loc[d.centre == c, col].values for c in cs}
    bs = np.array([np.mean(np.concatenate([by[c] for c in rng.choice(cs, len(cs))])) for _ in range(n_boot)])
    return {"mean": round(float(d[col].mean()), 4), "ci95": [round(float(q), 4) for q in np.quantile(bs, [0.025, 0.975])],
            "n_units": int(len(d)), "n_centres": int(len(cs)),
            "p_boot_two_sided": round(float(min(1.0, 2 * min((bs <= 0).mean(), (bs >= 0).mean()))), 4)}


def reading(r):
    lo, hi = r["ci95"]
    if -TOL <= lo and hi <= TOL:
        return "the forced share of the contrast is what the variance share predicts (the rule holds)"
    if (lo > 0 or hi < 0) and abs(r["mean"]) > TOL:
        return "the label selects forcing beyond the variance share (the rule fails)"
    return "inconclusive"


def upper95(k, n):
    from scipy.stats import beta
    return float(beta.ppf(0.95, k + 1, n - k)) if k < n else 1.0


# ------------------------------------------------------------------ P1 / S1 data
def snapsi_pairs():
    centres = sorted({p.name.split("_")[0] for p in L.RED.glob("*.parquet")})
    out = []
    for c in centres:
        bad, _ = L.corruption_guard(c)
        if bad:
            continue
        for init in L.ONSET:
            con, nud = L.load(c, "control", init), L.load(c, "nudged", init)
            if con is None or nud is None or not (L.spans_window(con, init) and L.spans_window(nud, init)):
                continue
            wm = L.window_means(con, init); base, sd = float(wm.mean()), float(wm.std(ddof=1))
            N = list(NE.member_series(con, init, base, sd).values())
            S = list(NE.member_series(nud, init, base, sd).values())
            out.append((c, init, S, N))
    return out


def free_pairs():
    centres = sorted({p.name.split("_")[0] for p in NE.FREE_RED.glob("*.parquet")} - {"NRL"})
    out = []
    for c in centres:
        for init in NE.NH_INITS:
            p, u = NE.load_psl(c, init), NE.load_u10(c, init)
            if p is None or u is None or not L.spans_window(p, init):
                continue
            lab = NE.classify_members(u, init)
            non = [m for m, v in lab.items() if v == "non"]; ssw = [m for m, v in lab.items() if v == "ssw"]
            if len(ssw) < NE.MIN_EACH or len(non) < NE.MIN_EACH:
                continue
            pm = p[p.member.isin(non)]; wm = L.window_means(pm, init)
            base, sd = float(wm.mean()), float(wm.std(ddof=1))
            S = list(NE.member_series(p[p.member.isin(ssw)], init, base, sd).values())
            N = list(NE.member_series(pm, init, base, sd).values())
            out.append((c, init, S, N))
    return out


def pairs_frame(pairs):
    rows = []
    for c, init, S, N in pairs:
        r = pair_share(S, N); r.update(centre=c, init=init, n_forced=len(S), n_ref=len(N)); rows.append(r)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ P2 / S2 data
def episodes(which):
    d = json.loads(CONT.read_text())
    t = pd.DataFrame(d["unit_table"] if which == "ncep" else d["era5_1940_2025"]["unit_table"])
    y = "Y1_AO" if which == "ncep" else "Y2E_NAMproxy"
    t = t.dropna(subset=["X", y, "Y4_downward"]).copy()
    if "winter" not in t:
        dt = pd.to_datetime(t.date); t["winter"] = np.where(dt.dt.month >= 7, dt.dt.year + 1, dt.dt.year)
    return t["X"].values.astype(float), t[y].values.astype(float), t["Y4_downward"].values > 0.5, t["winter"].values


def episode_stats(X, Y, c):
    b = np.cov(X, Y, ddof=0)[0, 1] / np.var(X)
    F = b * (X - X.mean())
    D, DF, R = share(Y, c, F)
    return {"D": D, "D_F": DF, "R_obs": R, "R_h2": h2(Y, F), "stat": R - h2(Y, F), "slope_per_10ms": 10 * b}


def winter_boot(X, Y, c, win, rng, n_boot=N_BOOT):
    ws = np.unique(win); idx = {w: np.where(win == w)[0] for w in ws}
    keys = ("R_obs", "R_h2", "stat"); bs = {k: [] for k in keys}
    for _ in range(n_boot):
        ii = np.concatenate([idx[w] for w in rng.choice(ws, len(ws))])
        s = episode_stats(X[ii], Y[ii], c[ii])
        for k in keys:
            bs[k].append(s[k])
    return {k: np.array(v) for k, v in bs.items()}


def summ(v, bs):
    b = bs[np.isfinite(bs)]
    return {"mean": round(float(v), 4), "ci95": [round(float(q), 4) for q in np.quantile(b, [0.025, 0.975])],
            "p_boot_two_sided": round(float(min(1.0, 2 * min((b <= 0).mean(), (b >= 0).mean()))), 4)}


def headline(r):
    lo, hi = r["ci95"]
    return ("most of the observed downward-non-downward contrast is not a difference in vortex dose" if hi < 0.5 else
            "most of it is a difference in vortex dose" if lo > 0.5 else "undetermined")


# ------------------------------------------------------------------ validation
def validate():
    rng = np.random.default_rng(SEED + 11)
    pairs = snapsi_pairs()
    deltas = {(c, i): float(np.mean([s.mean() for s in S]) - np.mean([s.mean() for s in N])) for c, i, S, N in pairs}
    res = {}
    for label, mix in (("law", False), ("couples_or_not", True)):
        hits, n = 0, 200
        for _ in range(n):
            syn = []
            for c, i, S, N in pairs:
                d = deltas[(c, i)]
                draw = [N[j] for j in rng.integers(len(N), size=len(S))]
                syn.append((c, i, [s + (d / 0.5) * (rng.random() < 0.5) if mix else s + d for s in draw], N))
            df = pairs_frame(syn)
            lo, hi = centre_summary(df, "stat", rng, 1000)["ci95"]
            hits += int(not (lo <= 0 <= hi))
        res[f"P1_{label}"] = {"rate_interval_excludes_0": hits / n, "upper95": round(upper95(hits, n), 4), "n_sim": n}
    X, Y, c, win = episodes("ncep")
    beta = np.cov(X, Y, ddof=0)[0, 1] / np.var(Y); resid = X - X.mean() - beta * (Y - Y.mean()); s = resid.std()
    hits, n = 0, 1000
    for _ in range(n):
        Xs = X.mean() + beta * (Y - Y.mean()) + rng.normal(0, s, len(Y))
        bs = winter_boot(Xs, Y, c, win, rng, 1000)["stat"]; bs = bs[np.isfinite(bs)]
        lo, hi = np.quantile(bs, [0.025, 0.975]); hits += int(lo > 0 or hi < 0)
    res["P2_law"] = {"rate_interval_excludes_0": hits / n, "upper95": round(upper95(hits, n), 4), "n_sim": n}
    res["P1_passes"] = res["P1_law"]["upper95"] < GATE
    res["P2_passes"] = res["P2_law"]["upper95"] < GATE
    out = json.loads(OUT.read_text()) if OUT.exists() else {}
    out["validation"] = res; out["validation_fixed"] = "2026-10-09, before P1/P2 were computed"
    OUT.write_text(json.dumps(out, indent=1), encoding="utf8", newline="\n")
    print(json.dumps(res, indent=1)); return 0


# ------------------------------------------------------------------ secondary sets
def tc_set(rng, basin_relative):
    d = pd.read_parquet(TCS)
    d = d[(d.season >= 1982) & (d.season <= 2024) & np.isfinite(d.sst_g)].reset_index(drop=True)
    if basin_relative:
        d = d[d.basin != "SA"].reset_index(drop=True)
        lev = np.zeros(len(d), int)
        for b, g in d.groupby("basin"):
            cuts = np.quantile(g.sst_g, [1 / 3, 2 / 3]); lev[g.index] = np.digitize(g.sst_g.values, cuts)
    else:
        lev = np.digitize(d.sst_g.values, np.quantile(d.sst_g, [1 / 3, 2 / 3]))
    seasons = d.season.values; us = np.unique(seasons); idx = {s: np.where(seasons == s)[0] for s in us}
    out = {}
    for name, col, thr in (("MAJOR", "y_lmi", MAJOR_KT), ("RI", "y_ri", RI_KT)):
        ok = np.isfinite(d[col].values)

        def stats(ii):
            x = d[col].values[ii]; lv = lev[ii]; m = ok[ii]; x, lv = x[m], lv[m]
            mu = np.array([x[lv == j].mean() for j in range(3)]); sd = np.array([x[lv == j].std(ddof=1) for j in range(3)])
            F = mu[lv]; c = x >= thr
            D, DF, R = share(x, c, F)
            ref = x[lv == 0]; pj = np.array([(lv == j).mean() for j in range(3)])
            xs = np.concatenate([mu[j] + sd[j] / sd[0] * (ref - mu[0]) for j in range(3)])
            Fs = np.concatenate([np.full(len(ref), mu[j]) for j in range(3)])
            ws = np.concatenate([np.full(len(ref), pj[j] / len(ref)) for j in range(3)])
            _, _, RLS = share(xs, xs >= thr, Fs, ws)
            return {"D": D, "D_F": DF, "R_obs": R, "R_LS": RLS, "R_h2": h2(x, F), "stat": R - RLS}
        s0 = stats(np.arange(len(d)))
        bs = np.array([stats(np.concatenate([idx[s] for s in rng.choice(us, len(us))]))["stat"] for _ in range(2000)])
        r = summ(s0["stat"], bs); r["reading"] = reading(r)
        out[name] = {k: round(float(v), 4) for k, v in s0.items()}; out[name]["stat_summary"] = r
    return out


def damip_set(rng):
    import damip_class_law as DM
    raw = pd.read_parquet(DM.DATA)
    out = {}
    for tgt, (box, season, side) in DM.TARGETS.items():
        s = DM.seasonal(raw, box, season); q = 2 / 3 if side == "hot" else 1 / 3
        rows = []
        for yr in DM.YEARS:
            xh = s[(s.arm == "historical") & (s.syear == yr)]["mean"].values
            xn = s[(s.arm == "hist-nat") & (s.syear == yr)]["mean"].values
            if len(xh) < 15 or len(xn) < 15:
                continue
            mu, sd = xn.mean(), xn.std(ddof=1); h, n = (xh - mu) / sd, (xn - mu) / sd
            thr = float(np.quantile(n, q)); cls = (lambda v: v > thr) if side == "hot" else (lambda v: v < thr)
            delta = h.mean() - n.mean(); k = h.std(ddof=1) / n.std(ddof=1)
            x = np.concatenate([h, n]); F = np.concatenate([np.full(len(h), delta), np.zeros(len(n))])
            D, DF, R = share(x, cls(x), F)
            mp_ = h.mean() + k * (n - n.mean()); p = len(h) / (len(h) + len(n))
            xl = np.concatenate([mp_, n])
            Fl = np.concatenate([np.full(len(mp_), delta), np.zeros(len(n))])
            wl = np.concatenate([np.full(len(mp_), p / len(mp_)), np.full(len(n), (1 - p) / len(n))])
            _, _, RLS = share(xl, cls(xl), Fl, wl)
            rows.append({"year": yr, "R_obs": R, "R_LS": RLS, "R_h2": h2(x, F), "stat": R - RLS})
        df = pd.DataFrame(rows).dropna()
        v = df.stat.values
        bs = np.array([v[rng.integers(len(v), size=len(v))].mean() for _ in range(N_BOOT)])
        r = summ(v.mean(), bs); r["reading"] = reading(r)
        out[tgt] = {"n_years": int(len(df)), "R_obs": round(float(df.R_obs.mean()), 4), "R_LS": round(float(df.R_LS.mean()), 4),
                    "R_h2": round(float(df.R_h2.mean()), 4), "stat_summary": r}
    return out


def cmip6_set(rng):
    import response_shape_dose as RSD
    W = RSD.W
    Xs, ys, gs = [], [], []
    for f in sorted(W.RAW.glob("*_zm.nc")):
        c = W.prepare_member(f)
        if c is None:
            continue
        Fe = W.H.features_for(c["u"], c["plev"], c["lat"], c["idx"], c["on"])
        Xs.append(Fe[RSD.DOSE].values); ys.append(W.C6.anom(c["am"], c["on"], c["cl"], W.OUT_WIN)); gs.append(np.full(len(ys[-1]), f.stem))
    X = np.concatenate(Xs); y = np.concatenate(ys); g = np.concatenate(gs)
    ok = np.isfinite(y) & np.isfinite(X).all(1); X, y, g = X[ok], y[ok], g[ok]
    mem = np.unique(g); idx = {m: np.where(g == m)[0] for m in mem}

    def stats(ii):
        Xi, yi, gi = X[ii], y[ii], g[ii]
        Xd = Xi.copy(); yd = yi.copy()
        for m in np.unique(gi):
            s = gi == m; Xd[s] -= Xd[s].mean(0); yd[s] -= yd[s].mean()
        beta = np.linalg.lstsq(Xd, yd, rcond=None)[0]; F = Xd @ beta
        D, DF, R = share(yi, yi < 0, F)
        return {"D": D, "D_F": DF, "R_obs": R, "R_h2": h2(yd, F), "stat": R - h2(yd, F)}
    s0 = stats(np.arange(len(y)))
    bs = [stats(np.concatenate([idx[m] for m in rng.choice(mem, len(mem))])) for _ in range(2000)]
    out = {k: round(float(v), 4) for k, v in s0.items()}
    out["n_events"] = int(len(y)); out["n_members"] = int(len(mem))
    for k in ("R_obs", "R_h2", "stat"):
        out[k + "_summary"] = summ(s0[k], np.array([b[k] for b in bs]))
    out["stat_summary"] = out.pop("stat_summary"); out["stat_summary"]["reading"] = reading(out["stat_summary"])
    return out


def two_point_max(s2, qs=(0.5, 0.65)):
    from scipy.optimize import brentq
    from scipy.stats import norm
    best = 0.0
    for q in qs:
        for p in np.linspace(0.01, 0.99, 99):
            for d in (np.sqrt(s2 / (p * (1 - p))), -np.sqrt(s2 / (p * (1 - p)))):
                c = brentq(lambda c: (1 - p) * norm.cdf(c) + p * norm.cdf(c - d) - q, -30, 30)
                P1 = p * norm.cdf(c - d)
                exlo = ((1 - p) * -norm.pdf(c) + p * (d * norm.cdf(c - d) - norm.pdf(c - d))) / q
                exhi = (p * d - q * exlo) / (1 - q)
                R = (d * P1 / q - (p * d - d * P1) / (1 - q)) / (exlo - exhi)
                best = max(best, R / (s2 / (1 + s2)))
    return round(float(best), 4)


def derived():
    fv = ROOT / "results" / "current"
    src = sorted(fv.rglob("*forced_variance*ceiling*.json")) + sorted(fv.rglob("forced_variance_ceiling*.json"))
    out = {"formula": "R_max = s2 / (1 + s2), s2 the forced between-event variance in units of the internal variance"}
    for name, s in (("cmip6_point_var_0.019", 0.019), ("cmip6_upper_sd_0.342", 0.342 ** 2), ("cmip6_upper_sd_0.269", 0.269 ** 2),
                    ("obs_upper_sd_0.596", 0.596 ** 2), ("obs_upper_sd_0.681", 0.681 ** 2),
                    ("cmip6_damped_upper_sd_0.42", 0.42 ** 2), ("obs_damped_upper_sd_0.78", 0.78 ** 2)):
        out[name] = {"R_max_gaussian": round(s / (1 + s), 4), "two_point_inflation": two_point_max(s),
                     "R_max_two_point": round(s / (1 + s) * two_point_max(s), 4)}
    out["icon_event_conditioned_between_over_within_0.8"] = {"R_max_gaussian": round(0.8 / 1.8, 4), "two_point_inflation": two_point_max(0.8)}
    out["sources"] = ["Methods 'Bounds on event differences' (upper s_f 0.342/0.269 CMIP6, 0.596/0.681 obs)",
                      "Methods 'Identification audit' (upper s_f at the strongest damping: 0.42 CMIP6, 0.78 obs)",
                      "Supplementary Note 1 (CMIP6 forced variance 0.019)", "icon_event_heterogeneity.json (~80% of within)"]
    out["source_files_found"] = [str(p.relative_to(ROOT)) for p in src]
    return out


# ------------------------------------------------------------------ main
def main():
    if "--validate" in sys.argv:
        return validate()
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    v = res.get("validation")
    if not v:
        sys.exit("validation gate not run (run --validate first)")
    rng = np.random.default_rng(SEED)
    res.update({"registered": "2026-10-09, committed before any statistic was computed", "seed": SEED, "n_boot": N_BOOT,
                "tolerance_share": TOL, "inputs": ["03_data_ingestion/_snapsi_reduced", "03_data_ingestion/_snapsi_reduced_free",
                                                  "03_data_ingestion/_snapsi_u10_free", str(CONT.relative_to(ROOT)),
                                                  "03_data_ingestion/tc_storms.parquet", "03_data_ingestion/damip_hadgem_boxes.parquet",
                                                  "CMIP6 *_zm.nc via response_shape_dose"]})
    p1 = pairs_frame(snapsi_pairs())
    r1 = centre_summary(p1, "stat", rng)
    r1["reading"] = reading(r1) if v["P1_passes"] else "not interpreted: validation gate failed"
    res["P1_snapsi_pooled_arms"] = {"stat_R_obs_minus_R_LS": r1, **{k: centre_summary(p1, k, rng) for k in ("R_obs", "R_LS", "R_h2")},
                                    "n_pairs": int(len(p1)), "pairs": p1.round(4).to_dict(orient="records")}
    print("P1", json.dumps(r1), flush=True)
    X, Y, c, win = episodes("ncep")
    s2 = episode_stats(X, Y, c); bs = winter_boot(X, Y, c, win, np.random.default_rng(SEED + 2))
    r2 = summ(s2["stat"], bs["stat"]); r2["reading"] = reading(r2) if v["P2_passes"] else "not interpreted: validation gate failed"
    ro = summ(s2["R_obs"], bs["R_obs"]); ro["reading"] = headline(ro)
    res["P2_obs_ncep_episodes"] = {"n_units": int(len(Y)), "n_downward": int(c.sum()), "D": round(s2["D"], 4), "D_F": round(s2["D_F"], 4),
                                   "slope_AO_per_10ms": round(s2["slope_per_10ms"], 4), "stat_R_obs_minus_R_h2": r2,
                                   "R_obs": ro, "R_h2": summ(s2["R_h2"], bs["R_h2"])}
    print("P2", json.dumps(res["P2_obs_ncep_episodes"]), flush=True)
    ps = [r1["p_boot_two_sided"], r2["p_boot_two_sided"]]; o = np.argsort(ps); adj = [0.0, 0.0]; run = 0.0
    for rank, i in enumerate(o):
        run = max(run, min(1.0, (2 - rank) * ps[i])); adj[i] = run
    res["holm_P1_P2"] = [round(a, 4) for a in adj]
    p_s1 = pairs_frame(free_pairs()); r = centre_summary(p_s1, "stat", np.random.default_rng(SEED + 3)); r["reading"] = reading(r)
    res["S1_free_natural_experiment"] = {"stat_R_obs_minus_R_LS": r, **{k: centre_summary(p_s1, k, np.random.default_rng(SEED + 4))
                                                                       for k in ("R_obs", "R_LS", "R_h2")}, "n_ensembles": int(len(p_s1))}
    X, Y, c, win = episodes("era5")
    s = episode_stats(X, Y, c); bs = winter_boot(X, Y, c, win, np.random.default_rng(SEED + 5))
    r = summ(s["stat"], bs["stat"]); r["reading"] = reading(r)
    res["S2_obs_era5_episodes"] = {"n_units": int(len(Y)), "stat_R_obs_minus_R_h2": r, "R_obs": summ(s["R_obs"], bs["R_obs"]),
                                   "R_h2": summ(s["R_h2"], bs["R_h2"]), "D": round(s["D"], 4)}
    res["S3_tropical_cyclones"] = {"global_terciles": tc_set(np.random.default_rng(SEED + 6), False),
                                   "basin_relative_terciles": tc_set(np.random.default_rng(SEED + 7), True)}
    res["S4_human_forcing_hadgem3"] = damip_set(np.random.default_rng(SEED + 8))
    res["S5_cmip6_events_dose"] = cmip6_set(np.random.default_rng(SEED + 9))
    res["S6_derived_bounds"] = derived()
    for k in ("S1_free_natural_experiment", "S2_obs_era5_episodes", "S3_tropical_cyclones", "S4_human_forcing_hadgem3",
              "S5_cmip6_events_dose", "S6_derived_bounds"):
        print(k, json.dumps(res[k], default=float)[:600], flush=True)
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
