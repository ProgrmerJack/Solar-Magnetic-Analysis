#!/usr/bin/env python3
"""
enso_extreme_class.py
=====================
IS "EXTREME EL NINO" A FIXED REGIME WHOSE FREQUENCY GREENHOUSE FORCING CHANGES, OR A
THRESHOLD ON ONE POPULATION THAT FORCING SHIFTS AND RESCALES?

Registered 2026-10-09 (plan approved by the user the same day), committed together with
its validation gate BEFORE any Nino3 value was retrieved (03_data_ingestion/
acquire_cmip6_nino3.py is committed in the same commit and had not been run).

LITERATURE (opened 2026-10-09)
  Cai et al. 2014, Nat. Clim. Change 4, 111 (doi:10.1038/nclimate2100): extreme El Nino =
    December-February (DJF) Nino3 (5S-5N, 150W-90W) total rainfall > 5 mm/day; rainfall
    "increases nonlinearly with Nino3 SST"; frequency doubles under greenhouse warming in
    models selected by Nino3 rainfall skewness > 1 and rainfall exceeding 5 mm/day.
  Takahashi & Dewitte 2016 (Clim. Dyn.) and Takahashi, Karamperidou & Dewitte 2019
    (Clim. Dyn. 52, 7477): strong and moderate El Nino are two regimes (bimodal peak SST),
    from an SST threshold for deep convection in the eastern Pacific.
  So both readings are live: a convective regime (fixed types) or one population pushed
  across a threshold.

DATA  CMIP6 Pangeo store, Amon pr (mm/day) and ts (K; over the all-ocean Nino3 box this is
  SST on the atmospheric grid), area-weighted mean over grid cells with centres in 5S-5N,
  210-270E. Members present in both historical and ssp585. DJF means labelled by the
  January year, all three months required.
  LEVELS  L0 DJF 1951-1980 (historical), L1 DJF 2021-2050 (ssp585), L2 DJF 2070-2099 (ssp585).
  A member enters if it has all 90 DJF seasons; a model (CanESM5 physics p1 and p2 analysed
  separately, as in class_law_general) enters if it has >= 15 members.
  CANDIDATES  CanESM5, MIROC6, ACCESS-ESM1-5, EC-Earth3 (pr); CanESM5, MIROC6, ACCESS-ESM1-5 (ts).
STEP 0, MODEL SELECTION (Cai et al.'s rule, applied before any test statistic): pooled DJF
  Nino3 rainfall over the three levels has skewness > 1 and a maximum > 5 mm/day. Models
  failing it are reported and enter only the unselected sensitivity.
TESTS (unit: member-season; inference: two-way bootstrap of mixture_moment_test, whole
  member trajectories x 5-year blocks of seasons, B = 1999; members are paired across
  levels because ssp585 continues historical)
  T1 PRIMARY  fixed two-type mixing on DJF Nino3 rainfall, mixture_moment_test.tests over
     L0, L1, L2, per selected model; Holm over models. Category per model:
     "incompatible with fixed two-type mixing" (Holm p < 0.05); "compatible" (95% interval
     of rho inside [-0.25, +0.25]); otherwise "unresolved".
     READING (fixed now): "extreme El Nino is not a fixed regime whose frequency forcing
     changes" if at least one selected model is incompatible and none is compatible;
     "consistent with a fixed regime" if a majority of selected models is compatible;
     otherwise "not resolved". If no model passes STEP 0, T1 is not run and that is reported.
  T2 SECONDARY  the same on DJF Nino3 ts.
  T3 SECONDARY  the class law for the extreme class (rainfall > 5 mm/day, absolute), from
     L0 to L1 and from L0 to L2: in units of the reference level's s.d., the forced level's
     contrast (class mean minus rest) minus the contrast of the reference level mapped
     x -> m_f + k (x - m_r), classified at 5 mm/day; seasons pooled within a level; interval
     from the same two-way bootstrap; read against +-0.25 as everywhere ("holds" if inside,
     "fails" if excluding 0 with |mean| > 0.25, else "inconclusive"); undefined if either
     class has < 10 member-seasons.
  T4 SECONDARY  forced share of the extreme-vs-rest contrast, levels pooled, F = the level's
     mean rainfall: R_obs, R_LS (each level the reference mapped by its own mean and s.d.,
     pooled with the observed level shares), h2 = Var(F)/Var(X) (class_forcing_content);
     statistic R_obs - R_LS read against +-0.05.
  T5 DESCRIPTIVE  shape: quantile differences of the level-standardised L2 against L0
     rainfall (5th-95th percentiles); bimodality of L0 DJF Nino3 ts: Gaussian mixture with 1
     vs 2 components (BIC) and the number of modes of the fitted two-component density.
  Also reported, unselected: T1 for every candidate model.
VALIDATION GATE (--validate; before any data). 1,000 simulations per null case, B = 999,
  three levels of 30 seasons x 50 members (ids shared across levels), year effects s.d. 0.1,
  no autocorrelation (DJF Nino3 has weak year-to-year memory):
  E1 two gamma types, moderate (mean 1.5, s.d. 1.0) and extreme (mean 8, s.d. 2.5) mm/day,
     extreme weight 0.05/0.10/0.20      E2 the same, weights 0.10/0.20/0.35
  E3 lognormal types of the same means and s.d., weights 0.05/0.10/0.20
  E4 as E1 with 36 members per level
  GATE: one-sided 95% Clopper-Pearson upper bound of the T false-rejection rate at 0.05
  below 0.07 in every null case. Power (reported, no gate): PA one Gaussian SST population,
  mean 0/0.5/1.0, through the fixed convex map rain = 1 + 2 exp(0.8 sst); PB one gamma
  population with constant coefficient of variation, mean 2/2.6/3.4 mm/day.
  A failed gate stops T1-T3 and is reported.
Output: results/current/8_experiment/enso_extreme_class.json
Run:    enso_extreme_class.py --validate ; enso_extreme_class.py
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, skew

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "8_experiment"
OUT = RESULTS / "enso_extreme_class.json"
DATA = HERE.parents[0] / "03_data_ingestion" / "cmip6_nino3.parquet"
sys.path.insert(0, str(HERE))
import mixture_moment_test as MM                    # noqa: E402

NAME = "enso_extreme_class"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
B, B_VAL, N_SIM, ALPHA, GATE, TOL, TOL_R = 1999, 999, 1000, 0.05, 0.07, 0.25, 0.05
EXTREME_MM = 5.0
LEVELS = {"L0": ("historical", 1951, 1980), "L1": ("ssp585", 2021, 2050), "L2": ("ssp585", 2070, 2099)}
MIN_MEMBERS, MIN_CLASS = 15, 10
QS = np.round(np.arange(0.05, 0.951, 0.05), 2)


# ------------------------------------------------------------------ data
def djf(d):
    """model-variant x var x member x Jan-year DJF means (all three months present)."""
    d = d.copy()
    d["model"] = np.where(d.model == "CanESM5", "CanESM5-" + d.member.str.extract(r"i1(p\d)")[0], d.model)
    d = d[d.month.isin([12, 1, 2])]
    d["jyear"] = np.where(d.month == 12, d.year + 1, d.year)
    g = d.groupby(["model", "var", "experiment", "member", "jyear"]).NINO3.agg(["mean", "count"]).reset_index()
    return g[g["count"] == 3].drop(columns="count")


def arms_for(g, model, var):
    """three arms (ids, X[members x 30]) with only members complete in every level."""
    s = g[(g.model == model) & (g["var"] == var)]
    mats = {}
    for L, (exp, y0, y1) in LEVELS.items():
        w = s[(s.experiment == exp) & s.jyear.between(y0, y1)].pivot(index="member", columns="jyear", values="mean")
        mats[L] = w.reindex(columns=range(y0, y1 + 1))
    common = sorted(set.intersection(*[set(m.dropna().index) for m in mats.values()]))
    return [(np.array(common), mats[L].loc[common].values.astype(float)) for L in LEVELS], common


# ------------------------------------------------------------------ statistics
def contrast_w(x, w, thr):
    c = x > thr
    if w[c].sum() <= 0 or w[~c].sum() <= 0:
        return np.nan
    return float(np.sum(w[c] * x[c]) / w[c].sum() - np.sum(w[~c] * x[~c]) / w[~c].sum())


def law_residual(ref, frc, thr, wr=None, wf=None):
    wr = np.ones(ref.size) if wr is None else wr; wf = np.ones(frc.size) if wf is None else wf
    if (wf * (frc > thr)).sum() < MIN_CLASS or (wf * (frc <= thr)).sum() < MIN_CLASS:
        return np.nan
    mr = np.sum(wr * ref) / wr.sum(); sr = np.sqrt(np.sum(wr * (ref - mr) ** 2) / wr.sum())
    mf = np.sum(wf * frc) / wf.sum(); sf = np.sqrt(np.sum(wf * (frc - mf) ** 2) / wf.sum())
    mapped = mf + sf / sr * (ref - mr)
    return (contrast_w(frc, wf, thr) - contrast_w(mapped, wr, thr)) / sr


def boot_weights(n_ids, nY, rng, B_):
    return MM.weights(rng, n_ids, nY, B_)


def law_test(arms, a, b, thr, rng):
    (ids, Xr), (_, Xf) = arms[a], arms[b]
    r0 = law_residual(Xr.ravel(), Xf.ravel(), thr)
    Wm, Wy = boot_weights(len(ids), Xr.shape[1], rng, B)
    bs = []
    for k in range(B):
        w = np.outer(Wm[k], Wy[k]).ravel()
        bs.append(law_residual(Xr.ravel(), Xf.ravel(), thr, w, w))
    bs = np.array(bs); bs = bs[np.isfinite(bs)]
    if not np.isfinite(r0) or len(bs) < 0.9 * B:
        return {"residual": None, "reading": "undefined (class below 10 member-seasons)"}
    lo, hi = np.quantile(bs, [0.025, 0.975])
    rd = "holds" if (-TOL <= lo and hi <= TOL) else "fails" if ((lo > 0 or hi < 0) and abs(r0) > TOL) else "inconclusive"
    n_ext = [int((Xr > thr).sum()), int((Xf > thr).sum())]
    return {"residual": round(float(r0), 4), "ci95": [round(float(lo), 4), round(float(hi), 4)], "reading": rd,
            "n_extreme_ref_forced": n_ext, "rate_ref_forced": [round(n_ext[0] / Xr.size, 4), round(n_ext[1] / Xf.size, 4)]}


def forced_share(arms, thr, rng):
    def stats(Xs, ws=None):
        x = np.concatenate([X.ravel() for X in Xs])
        lev = np.concatenate([np.full(X.size, j) for j, X in enumerate(Xs)])
        w = np.ones(x.size) if ws is None else np.concatenate(ws)
        mu = np.array([np.sum(w[lev == j] * x[lev == j]) / w[lev == j].sum() for j in range(3)])
        sd = np.array([np.sqrt(np.sum(w[lev == j] * (x[lev == j] - mu[j]) ** 2) / w[lev == j].sum()) for j in range(3)])
        F = mu[lev]; c = x > thr
        if (w * c).sum() < MIN_CLASS or (w * ~c).sum() < MIN_CLASS:
            return None
        mx = lambda v, m: np.sum(w[m] * v[m]) / w[m].sum()
        D = mx(x, c) - mx(x, ~c); R = (mx(F, c) - mx(F, ~c)) / D
        Fm = np.sum(w * F) / w.sum(); xm = np.sum(w * x) / w.sum()
        h2 = np.sum(w * (F - Fm) ** 2) / np.sum(w * (x - xm) ** 2)
        ref = x[lev == 0]; wr = w[lev == 0]; pj = np.array([w[lev == j].sum() for j in range(3)]) / w.sum()
        xs = np.concatenate([mu[j] + sd[j] / sd[0] * (ref - mu[0]) for j in range(3)])
        Fs = np.concatenate([np.full(ref.size, mu[j]) for j in range(3)])
        wsy = np.concatenate([pj[j] * wr / wr.sum() for j in range(3)])
        cs = xs > thr
        D_ = np.sum(wsy[cs] * xs[cs]) / wsy[cs].sum() - np.sum(wsy[~cs] * xs[~cs]) / wsy[~cs].sum()
        RLS = (np.sum(wsy[cs] * Fs[cs]) / wsy[cs].sum() - np.sum(wsy[~cs] * Fs[~cs]) / wsy[~cs].sum()) / D_
        return {"R_obs": R, "R_LS": RLS, "R_h2": h2, "stat": R - RLS}
    Xs = [a[1] for a in arms]; s0 = stats(Xs)
    if s0 is None:
        return {"reading": "undefined"}
    Wm, Wy = boot_weights(len(arms[0][0]), Xs[0].shape[1], rng, B)
    bs = []
    for k in range(B):
        w = np.outer(Wm[k], Wy[k]).ravel()
        r = stats(Xs, [w, w, w])
        if r:
            bs.append(r["stat"])
    lo, hi = np.quantile(bs, [0.025, 0.975])
    rd = ("the forced share is what the law predicts" if (-TOL_R <= lo and hi <= TOL_R) else
          "the class selects forcing beyond the law" if ((lo > 0 or hi < 0) and abs(s0["stat"]) > TOL_R) else "inconclusive")
    return {**{k: round(float(v), 4) for k, v in s0.items()}, "stat_ci95": [round(float(lo), 4), round(float(hi), 4)], "reading": rd}


def shape_and_modes(arms, ts_arms):
    out = {}
    z0 = arms[0][1].ravel(); z2 = arms[2][1].ravel()
    z0 = (z0 - z0.mean()) / z0.std(); z2 = (z2 - z2.mean()) / z2.std()
    out["rain_std_quantile_diff_L2_minus_L0"] = dict(zip([str(q) for q in QS], np.round(np.quantile(z2, QS) - np.quantile(z0, QS), 3).tolist()))
    if ts_arms is not None:
        from sklearn.mixture import GaussianMixture
        x = ts_arms[0][1].ravel().reshape(-1, 1)
        g1 = GaussianMixture(1, random_state=0).fit(x); g2 = GaussianMixture(2, n_init=10, random_state=0).fit(x)
        grid = np.linspace(x.min(), x.max(), 4001).reshape(-1, 1)
        f = np.exp(g2.score_samples(grid))
        modes = int(np.sum((f[1:-1] > f[:-2]) & (f[1:-1] > f[2:])))
        out["L0_ts_bic_1_minus_2"] = round(float(g1.bic(x) - g2.bic(x)), 2)
        out["L0_ts_modes_of_2component_fit"] = modes
        out["L0_ts_skewness"] = round(float(skew(x.ravel())), 3)
    return out


def category(r, p_holm):
    ci = r["rho_ci95"]
    return ("incompatible with fixed two-type mixing" if p_holm < 0.05 else
            "compatible" if ci and -TOL <= ci[0] and ci[1] <= TOL else "unresolved")


def holm(ps):
    o = np.argsort(ps); adj = np.empty(len(ps)); run = 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (len(ps) - r) * ps[i])); adj[i] = run
    return adj


# ------------------------------------------------------------------ validation
def sim(kind, rng, nY=30):
    n = 36 if kind == "E4" else 50
    from scipy.stats import gamma
    def g(mean, sd, size):
        k = (mean / sd) ** 2; return gamma.rvs(k, scale=mean / k, size=size, random_state=rng)
    def ln(mean, sd, size):
        s2 = np.log(1 + (sd / mean) ** 2); return rng.lognormal(np.log(mean) - s2 / 2, np.sqrt(s2), size)
    arms = []
    ids = np.array([f"m{i}" for i in range(n)])
    for j in range(3):
        ye = 0.1 * rng.normal(size=nY)[None, :]
        if kind in ("E1", "E2", "E3", "E4"):
            p = {"E1": (0.05, 0.10, 0.20), "E2": (0.10, 0.20, 0.35), "E3": (0.05, 0.10, 0.20), "E4": (0.05, 0.10, 0.20)}[kind][j]
            draw = ln if kind == "E3" else g
            typ = rng.random((n, nY)) < p
            X = np.where(typ, draw(8.0, 2.5, (n, nY)), draw(1.5, 1.0, (n, nY)))
        elif kind == "PA":
            X = 1 + 2 * np.exp(0.8 * (rng.normal(size=(n, nY)) + (0.0, 0.5, 1.0)[j]))
        else:
            m = (2.0, 2.6, 3.4)[j]; X = g(m, 0.6 * m, (n, nY))
        arms.append((ids, X + ye))
    return arms


def _val_job(args):
    kind, seed, nsim = args
    rng = np.random.default_rng(seed); rej = 0
    for _ in range(nsim):
        rej += MM.tests(sim(kind, rng), B_VAL, rng)["p_T"] < ALPHA
    return kind, rej


def validate():
    kinds = ["E1", "E2", "E3", "E4", "PA", "PB"]; chunks = 10
    jobs = [(k, SEED + 1000 * i + c, N_SIM // chunks) for i, k in enumerate(kinds) for c in range(chunks)]
    import multiprocessing as mp
    with mp.get_context("fork").Pool(8) as pool:
        outs = pool.map(_val_job, jobs)
    agg = {}
    for k, r in outs:
        agg[k] = agg.get(k, 0) + r
    ub = lambda r: float(beta.ppf(0.95, r + 1, N_SIM - r)) if r < N_SIM else 1.0
    val = {"n_sim": N_SIM, "B": B_VAL, "cases": {k: {"T_reject": v / N_SIM, "T_upper95": round(ub(v), 4)} for k, v in agg.items()}}
    val["passes"] = all(val["cases"][k]["T_upper95"] < GATE for k in ("E1", "E2", "E3", "E4"))
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res["validation"] = val; res["validation_fixed"] = "2026-10-09, before any Nino3 value was retrieved"
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(json.dumps(val, indent=1)); return 0


# ------------------------------------------------------------------ main
def main():
    if "--validate" in sys.argv:
        return validate()
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    if not res.get("validation", {}).get("passes"):
        sys.exit("validation gate not passed (run --validate)")
    import platform
    import scipy
    import sklearn
    g = djf(pd.read_parquet(DATA))
    res.update({"registered": "2026-10-09, committed before any Nino3 value was retrieved", "seed": SEED, "B": B,
                "inputs": [str(DATA.relative_to(ROOT))], "command": " ".join(sys.argv),
                "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                             "scipy": scipy.__version__, "sklearn": sklearn.__version__}})
    models = sorted(g.model.unique())
    sel, info, arms_pr, arms_ts = [], {}, {}, {}
    for m in models:
        a, common = arms_for(g, m, "pr")
        x = np.concatenate([X.ravel() for _, X in a]) if common else np.array([])
        info[m] = {"n_members_pr": len(common)}
        if len(common) >= MIN_MEMBERS:
            arms_pr[m] = a
            info[m].update({"pooled_skewness": round(float(skew(x)), 3), "max_mm_day": round(float(x.max()), 3),
                            "mean_L0_L1_L2": [round(float(X.mean()), 3) for _, X in a],
                            "sd_L0_L1_L2": [round(float(X.std()), 3) for _, X in a],
                            "extreme_rate_L0_L1_L2": [round(float((X > EXTREME_MM).mean()), 4) for _, X in a]})
            info[m]["passes_cai_rule"] = bool(info[m]["pooled_skewness"] > 1 and info[m]["max_mm_day"] > EXTREME_MM)
            if info[m]["passes_cai_rule"]:
                sel.append(m)
        t, ct = arms_for(g, m, "ts")
        info[m]["n_members_ts"] = len(ct)
        if len(ct) >= MIN_MEMBERS:
            arms_ts[m] = t
    res["step0_models"] = info; res["selected"] = sel
    print("STEP0", json.dumps(info), flush=True)
    t1 = {m: MM.tests(arms_pr[m], B, np.random.default_rng(SEED + 1 + i)) for i, m in enumerate(sorted(arms_pr))}
    for fam, keys in (("selected", sel), ("all", sorted(arms_pr))):
        if keys:
            adj = holm([t1[m]["p_T"] for m in keys])
            for m, a in zip(keys, adj):
                t1[m][f"p_T_holm_{fam}"] = round(float(a), 4); t1[m][f"category_{fam}"] = category(t1[m], a)
    res["T1_rain_mixing"] = t1
    cats = [t1[m]["category_selected"] for m in sel]
    inc = sum(c.startswith("incompatible") for c in cats); com = sum(c == "compatible" for c in cats)
    res["T1_reading"] = ("not run: no model passes Cai et al.'s rule" if not sel else
                         "extreme El Nino is not a fixed regime whose frequency forcing changes" if inc >= 1 and com == 0 else
                         "consistent with a fixed regime" if com > len(sel) / 2 else "not resolved")
    print("T1", res["T1_reading"], json.dumps({m: (t1[m]["rho"], t1[m]["rho_ci95"], t1[m]["p_T"]) for m in t1}), flush=True)
    t2 = {m: MM.tests(arms_ts[m], B, np.random.default_rng(SEED + 20 + i)) for i, m in enumerate(sorted(arms_ts))}
    if t2:
        adj = holm([t2[m]["p_T"] for m in sorted(t2)])
        for m, a in zip(sorted(t2), adj):
            t2[m]["p_T_holm"] = round(float(a), 4); t2[m]["category"] = category(t2[m], a)
    res["T2_ts_mixing"] = t2
    res["T3_law_extreme"] = {m: {"L0_to_L1": law_test(arms_pr[m], 0, 1, EXTREME_MM, np.random.default_rng(SEED + 40 + i)),
                                 "L0_to_L2": law_test(arms_pr[m], 0, 2, EXTREME_MM, np.random.default_rng(SEED + 60 + i))}
                             for i, m in enumerate(sorted(arms_pr))}
    res["T4_forced_share"] = {m: forced_share(arms_pr[m], EXTREME_MM, np.random.default_rng(SEED + 80 + i)) for i, m in enumerate(sorted(arms_pr))}
    res["T5_shape_modes"] = {m: shape_and_modes(arms_pr[m], arms_ts.get(m)) for m in sorted(arms_pr)}
    for k in ("T2_ts_mixing", "T3_law_extreme", "T4_forced_share", "T5_shape_modes"):
        print(k, json.dumps(res[k], default=float)[:1500], flush=True)
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
