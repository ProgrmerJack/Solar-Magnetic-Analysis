#!/usr/bin/env python3
"""
s2s_member_experiment.py
========================
A NATURAL EXPERIMENT INSIDE OPERATIONAL ENSEMBLES: MEMBERS THAT DO AND DO NOT
REVERSE THE VORTEX FROM THE SAME INITIAL STATE.

Plan approved 2026-10-02 (user: "Please find a way to address the items 1-6!";
standing approval of new analyses). Committed BEFORE any member-level split by
reversal was computed. The per-member reforecast files on disk had been used only
for event-based tests (P', regional, lead-lag), never split by member reversal.

WHY
  The decisive causal evidence (SNAPSI) imposes two observed Northern Hemisphere
  SSWs. Every reforecast start of the ten S2S systems (4,500 starts, 1991-2024
  hindcast years, December-March) shares one initial state across its members;
  some members spontaneously reverse the 10 hPa, 60N wind and others do not. That
  is a matched comparison like SNAPSI's -- same initial conditions, SSW or no SSW
  -- over hundreds of winters, in ten other models, with SSWs that arise on their
  own rather than being imposed.

DATA (03_data_ingestion, per system: s2s_<c>_u10_60N[_long], _psl_cap[_short],
_t2m_regions[_short]; the P' confirmatory systems plus ECMWF; model versions as P')
  u(10 hPa, 60N) leads 1-34; polar-cap mean sea-level pressure leads 1-34;
  northern-Eurasian 2 m temperature leads 1-33.

MEMBERS
  Reversing member: westerly at leads 1-2, first westerly-to-easterly change of
  u(10 hPa, 60N) at a lead k_m in 3..20. Non-reversing: westerly at leads 1-2 and
  westerly throughout leads 3..20. Others (easterly at start) are dropped.
  Mixed start: at least two reversing and two non-reversing members.
  Anchor k* = median onset lead of the reversing members (rounded down).
OUTCOMES (window means, leads k*+8 .. k*+25 for the NAM proxy, k*+8 .. k*+24 for
  temperature; a start enters only if the window lies within the data, i.e.
  k* <= 9; secondary windows k*+8 .. k*+20 with k* <= 14)
  A: polar-cap NAM proxy = -(psl_cap - its mean over all starts of the system with
     the same start month-day and lead); T: northern-Eurasian temperature anomaly,
     same construction. Units: the system's s.d. of the window mean over all
     members of all starts (sigma_sys).
TESTS (pooled over systems and mixed starts; 2,000 bootstrap resamples of
  hindcast winters, clustering all starts and systems of a winter)
  N1 SHIFT: mean over starts of mean(reversing) - mean(non-reversing).
  N2 TRANSLATION: pooled within-start variance of reversing members about their
     own start mean over that of non-reversing members (each group's sum of
     squares over its degrees of freedom).
  N3 THRESHOLD: downward label (window mean < 0 and > 50% of days < 0) in each
     member; rate in each group; class contrast (DW minus NDW mean, pooled within
     group after removing start means) in each group; difference of contrasts.
  N4 STEP AT REVERSAL, INITIAL STATE FIXED: all members of all starts (westerly at
     leads 1-2); running variable = the member's minimum u over leads 3-20; key
     lead = lead of that minimum, <= 9; outcome window key+8 .. key+25; start fixed
     effects (outcome and running variable demeaned within start) ; dose slope per
     10 m/s and the coefficient of 1(u_min < 0) with separate slopes, as the
     observational test; covariate: the member's NAM proxy over leads 1-2.
  N5 BALANCE: NAM proxy over leads 1..k*-1 (before onset), reversing minus
     non-reversing -- how far members that will reverse already differ.
  N6 TEMPERATURE: N1-N3 for northern-Eurasian temperature.
  Descriptive: quantile shift (q 0.05..0.95) of reversing minus non-reversing,
  standardised within start; skewness of each group.
READING (fixed now): one shifted population predicts N1 < 0, N2 within about
  [0.9, 1.1] (no added spread), N3 contrast difference ~ 0, N4 dose > 0 with a step
  consistent with 0. Two kinds of response predict N2 > 1 (a mixture widens), a
  contrast that differs between the groups, or a step at reversal. A large N5
  means the comparison is not clean (members that reverse were already
  different) and is reported as such; N1 is then also reported net of N5.

Output: results/current/5_mechanism/s2s_member_experiment.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import skew

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RESULTS = ROOT / "results" / "current" / "5_mechanism"
NAME = "s2s_member_experiment"
SEED = zlib.crc32(NAME.encode()) % (2 ** 32)
N_BOOT = 2000
SYSTEMS = ["ecmwf", "eccc", "cma", "hmcr", "kma", "cnrm", "jma", "cnr_isac", "ncep", "cptec"]
LEADS = np.arange(1, 35)
QS = np.round(np.arange(0.05, 0.951, 0.05), 2)


def tag(c):
    return "ecmf" if c == "ecmwf" else c


def piv(df, col, maxlead):
    df = df[df["lead_day"] <= maxlead]
    w = df.pivot_table(index=["init", "member"], columns="lead_day", values=col)
    return w.reindex(columns=range(1, maxlead + 1))


def anomaly(w):
    md = pd.to_datetime(w.index.get_level_values("init")).strftime("%m-%d")
    clim = w.groupby(md).transform("mean")
    return w - clim


def load(c):
    t = tag(c)
    u = pd.concat([pd.read_parquet(ING / f"s2s_{t}_u10_60N.parquet"), pd.read_parquet(ING / f"s2s_{t}_u10_60N_long.parquet")])
    p = pd.concat([pd.read_parquet(ING / f"s2s_{t}_psl_cap_short.parquet"), pd.read_parquet(ING / f"s2s_{t}_psl_cap.parquet")])
    tt = pd.concat([pd.read_parquet(ING / f"s2s_{t}_t2m_regions_short.parquet"), pd.read_parquet(ING / f"s2s_{t}_t2m_regions.parquet")])
    for d in (u, p, tt):
        d["init"] = pd.to_datetime(d["init"])
        d.drop_duplicates(["init", "member", "lead_day"], inplace=True)
    U = piv(u, "u10_60N", 34)
    A = -anomaly(piv(p, "psl_cap_N", 34))              # NAM proxy (Pa, NAM sign)
    T = anomaly(piv(tt, "NEURASIA", 33))
    idx = U.index.intersection(A.index).intersection(T.index)
    return U.loc[idx], A.loc[idx], T.loc[idx]


def classify(u):
    """-> (status, onset lead): status 'rev', 'non' or None."""
    if not (np.all(u[:2] > 0) and np.all(np.isfinite(u[:20]))):
        return None, None
    for k in range(3, 21):
        if u[k - 1] < 0 and u[k - 2] > 0:
            return "rev", k
    return ("non", None) if np.all(u[2:20] > 0) else (None, None)


def winter(t):
    return t.year if t.month <= 6 else t.year + 1


def build(c):
    U, A, T = load(c)
    starts, members = [], []
    # sigma of window means by window start k (all members, all starts)
    sigA = {k: float(np.nanstd(A.loc[:, k + 8:k + 25].mean(axis=1))) for k in range(3, 15)}
    sigA20 = {k: float(np.nanstd(A.loc[:, k + 8:k + 20].mean(axis=1))) for k in range(3, 15)}
    sigT = {k: float(np.nanstd(T.loc[:, k + 8:k + 24].mean(axis=1))) for k in range(3, 10)}
    for init, g in U.groupby(level="init"):
        u = g.values
        st = [classify(x) for x in u]
        rev = [i for i, (s_, _) in enumerate(st) if s_ == "rev"]
        non = [i for i, (s_, _) in enumerate(st) if s_ == "non"]
        a, t = A.loc[init].values, T.loc[init].values
        # N4 members: westerly at leads 1-2
        for i, x in enumerate(u):
            if np.all(x[:2] > 0) and np.all(np.isfinite(x[2:20])):
                kmin = int(np.argmin(x[2:20])) + 3
                if kmin <= 9:
                    y = np.nanmean(a[i, kmin + 7:kmin + 25]) / sigA[kmin]
                    pre = np.nanmean(a[i, :2]) / sigA[kmin]
                    members.append({"sys": c, "init": init, "winter": winter(init), "umin": float(x[kmin - 1]),
                                    "y": y, "pre": pre})
        if len(rev) < 2 or len(non) < 2:
            continue
        ks = int(np.floor(np.median([st[i][1] for i in rev])))
        row = {"sys": c, "init": init, "winter": winter(init), "k": ks, "n_rev": len(rev), "n_non": len(non)}
        if ks <= 14:
            w20 = a[:, ks + 7:ks + 20].mean(axis=1) / sigA20[ks]
            row["A20_rev"], row["A20_non"] = w20[rev].tolist(), w20[non].tolist()
        if ks <= 9:
            wa = a[:, ks + 7:ks + 25]
            ma = wa.mean(axis=1) / sigA[ks]
            dw = (wa.mean(axis=1) < 0) & ((wa < 0).mean(axis=1) > 0.5)
            wt = t[:, ks + 7:ks + 24].mean(axis=1) / sigT[ks]
            dwt = wt < 0
            pre = a[:, :max(ks - 1, 1)].mean(axis=1) / sigA[ks]
            row.update({"A_rev": ma[rev].tolist(), "A_non": ma[non].tolist(),
                        "DW_rev": dw[rev].astype(int).tolist(), "DW_non": dw[non].astype(int).tolist(),
                        "T_rev": wt[rev].tolist(), "T_non": wt[non].tolist(),
                        "pre_rev": pre[rev].tolist(), "pre_non": pre[non].tolist()})
        starts.append(row)
    return starts, members


def stats_starts(S, key="A"):
    """N1 shift, N2 variance ratio, N3 rates/contrasts from a list of start rows."""
    d, ssr, dfr, ssn, dfn = [], 0.0, 0, 0.0, 0
    for r in S:
        a, b = np.array(r[f"{key}_rev"]), np.array(r[f"{key}_non"])
        d.append(a.mean() - b.mean())
        ssr += ((a - a.mean()) ** 2).sum(); dfr += len(a) - 1
        ssn += ((b - b.mean()) ** 2).sum(); dfn += len(b) - 1
    return float(np.mean(d)), float((ssr / dfr) / (ssn / dfn))


def contrasts(S):
    out = {}
    for g in ("rev", "non"):
        v, l = [], []
        for r in S:
            a = np.array(r[f"A_{g}"]); v.append(a - a.mean()); l.append(np.array(r[f"DW_{g}"]))
        v, l = np.concatenate(v), np.concatenate(l).astype(bool)
        out[g] = {"rate": float(l.mean()), "contrast": float(v[l].mean() - v[~l].mean()) if l.any() and (~l).any() else np.nan}
    out["contrast_diff"] = out["rev"]["contrast"] - out["non"]["contrast"]
    return out


def boot(units, fn, rng):
    w = np.array([u["winter"] for u in units]); uw = np.unique(w)
    by = {x: [u for u in units if u["winter"] == x] for x in uw}
    out = []
    for _ in range(N_BOOT):
        pick = rng.choice(uw, len(uw))
        out.append(fn([u for x in pick for u in by[x]]))
    return np.array(out, dtype=float)


def ci(a):
    return [round(float(x), 4) for x in np.nanpercentile(a, [2.5, 97.5])]


def n4(M, rng):
    df = pd.DataFrame(M)
    df["D"] = (df["umin"] < 0).astype(float)
    def fit(d):
        g = d.groupby("init")
        cols = ["y", "umin", "pre", "D"]
        dm = d[cols] - g[cols].transform("mean")
        dm["Dx"] = (d["D"] * d["umin"]) - (d["D"] * d["umin"]).groupby(d["init"]).transform("mean")
        X1 = dm[["umin", "pre"]].values; b1 = np.linalg.lstsq(X1, dm["y"].values, rcond=None)[0]
        X2 = dm[["umin", "pre", "D", "Dx"]].values; b2 = np.linalg.lstsq(X2, dm["y"].values, rcond=None)[0]
        return 10 * b1[0], b2[2]
    est = fit(df)
    uw = df["winter"].unique(); byw = {w: df[df.winter == w] for w in uw}
    bs = np.array([fit(pd.concat([byw[w] for w in rng.choice(uw, len(uw))])) for _ in range(N_BOOT)])
    return {"n_members": int(len(df)), "n_below": int(df.D.sum()), "n_starts": int(df.init.nunique()),
            "dose_per_10ms_sigma": round(float(est[0]), 4), "dose_ci95": ci(bs[:, 0]),
            "step_sigma": round(float(est[1]), 4), "step_ci95": ci(bs[:, 1]),
            "step_ci90": [round(float(x), 4) for x in np.nanpercentile(bs[:, 1], [5, 95])],
            "step_p_two_sided": round(float(min(1, 2 * min(np.mean(bs[:, 1] <= 0), np.mean(bs[:, 1] >= 0)))), 4)}


def main():
    rng = np.random.default_rng(SEED)
    S, M = [], []
    per = {}
    for c in SYSTEMS:
        s_, m_ = build(c)
        S += s_; M += m_
        per[c] = {"mixed_starts": len(s_), "mixed_starts_k_le_9": sum(1 for r in s_ if "A_rev" in r),
                  "n4_members": len(m_)}
        print(c, per[c], flush=True)
    S9 = [r for r in S if "A_rev" in r]
    S14 = [r for r in S if "A20_rev" in r]
    res = {"plan_approved": "2026-10-02", "registered_commit": "be804a7", "seed": SEED, "n_boot": N_BOOT,
           "systems": SYSTEMS, "per_system": per, "n_mixed_starts_k_le_9": len(S9),
           "n_winters": int(len({r["winter"] for r in S9})),
           "n_members_rev": int(sum(r["n_rev"] for r in S9)), "n_members_non": int(sum(r["n_non"] for r in S9))}
    sh, vr = stats_starts(S9, "A")
    b = boot(S9, lambda u: stats_starts(u, "A"), rng)
    res["N1_shift_sigma"] = {"est": round(sh, 4), "ci95": ci(b[:, 0])}
    res["N2_variance_ratio"] = {"est": round(vr, 4), "ci95": ci(b[:, 1])}
    c3 = contrasts(S9)
    bc = boot(S9, lambda u: (lambda o: (o["rev"]["rate"], o["non"]["rate"], o["contrast_diff"]))(contrasts(u)), rng)
    res["N3_threshold"] = {"rate_rev": round(c3["rev"]["rate"], 4), "rate_non": round(c3["non"]["rate"], 4),
                           "contrast_rev": round(c3["rev"]["contrast"], 4), "contrast_non": round(c3["non"]["contrast"], 4),
                           "contrast_diff": round(c3["contrast_diff"], 4), "contrast_diff_ci95": ci(bc[:, 2]),
                           "rate_rev_ci95": ci(bc[:, 0]), "rate_non_ci95": ci(bc[:, 1])}
    res["N4_step_initial_state_fixed"] = n4(M, rng)
    shp, _ = stats_starts(S9, "pre")
    bp = boot(S9, lambda u: stats_starts(u, "pre"), rng)
    res["N5_balance_pre_onset_sigma"] = {"est": round(shp, 4), "ci95": ci(bp[:, 0]),
                                         "N1_net_of_N5": round(sh - shp, 4)}
    sht, vrt = stats_starts(S9, "T")
    bt = boot(S9, lambda u: stats_starts(u, "T"), rng)
    res["N6_temperature"] = {"shift_sigma": round(sht, 4), "shift_ci95": ci(bt[:, 0]),
                             "variance_ratio": round(vrt, 4), "variance_ratio_ci95": ci(bt[:, 1])}
    s20, v20 = stats_starts(S14, "A20")
    b20 = boot(S14, lambda u: stats_starts(u, "A20"), rng)
    res["secondary_window_8_20"] = {"n_starts": len(S14), "shift": round(s20, 4), "shift_ci95": ci(b20[:, 0]),
                                    "variance_ratio": round(v20, 4), "variance_ratio_ci95": ci(b20[:, 1])}
    # descriptive: quantile shift and skewness, standardised within start by the non-reversing mean
    rv = np.concatenate([np.array(r["A_rev"]) - np.mean(r["A_non"]) for r in S9])
    nv = np.concatenate([np.array(r["A_non"]) - np.mean(r["A_non"]) for r in S9])
    rd = np.concatenate([np.array(r["A_rev"]) - np.mean(r["A_rev"]) for r in S9])
    res["descriptive"] = {"quantile_shift": dict(zip([str(q) for q in QS], [round(float(x), 4) for x in np.quantile(rv, QS) - np.quantile(nv, QS)])),
                          "skew_rev_within": round(float(skew(rd)), 4), "skew_non_within": round(float(skew(nv)), 4)}
    for k in ("N1_shift_sigma", "N2_variance_ratio", "N3_threshold", "N4_step_initial_state_fixed",
              "N5_balance_pre_onset_sigma", "N6_temperature", "secondary_window_8_20", "descriptive"):
        print(k, json.dumps(res[k]), flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "s2s_member_experiment.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> s2s_member_experiment.json")


if __name__ == "__main__":
    main()
