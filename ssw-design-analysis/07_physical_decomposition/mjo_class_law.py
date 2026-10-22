#!/usr/bin/env python3
"""
mjo_class_law.py
================
DO MJO EVENTS THAT CROSS THE MARITIME CONTINENT AND EVENTS IT BLOCKS FORM TWO FIXED TYPES WHOSE
FREQUENCY THE QBO CHANGES, OR IS "CROSSING" A THRESHOLD ON ONE POPULATION OF PROPAGATION EXTENTS
THAT THE QBO SHIFTS?

Registered 2026-10-22 (the user asked for the MJO test, the optional third analysis of the approved
plan), committed with its validation gate BEFORE any RMM or QBO value was retrieved
(03_data_ingestion/acquire_indices.py gains the two sources in the same commit).

LITERATURE (retrieved 2026-10-22; abstracts and repository copies)
  Zhang & Ling 2017, J. Climate (barrier effect of the Maritime Continent): MJO events are sorted into
    those that propagate through the Maritime Continent and those that are blocked; about half of the
    Indian Ocean events get through. Ling et al. tie the blocking to the diurnal cycle of land convection,
    a mechanism specific to the blocked class.
  Son et al. 2017, J. Climate 30, 1909: boreal-winter MJO activity around the Maritime Continent is
    stronger and more organised in easterly QBO winters. Zhang & Zhang 2018: more MJO days in easterly
    QBO winters, from more Indian Ocean initiations and a weaker barrier effect.
  Wheeler & Hendon 2004, Mon. Wea. Rev.: the RMM index (BoM file; ENSO signal removed before 2014,
    only the 120-day mean afterwards).

DATA  BoM RMM daily (year, month, day, RMM1, RMM2, phase, amplitude), 1979-03 onward (the OLR record
  has a gap in 1978); CPC monthly 50 hPa QBO zonal wind index.
EVENTS  start = a day with amplitude >= 1 in phase 2 or 3 (Indian Ocean) whose previous day was not,
  with no event started in the 30 preceding days; the event runs while the amplitude stays >= 1, ending
  at the last such day before 3 consecutive days below 1, or after 60 days. Phase position from the
  angle theta = atan2(RMM2, RMM1): pos = ((theta + 180 deg) mod 360) / 45 + 1 (phase k spans [k, k+1));
  checked against the file's phase column (>= 99% agreement required, else the run stops). Unwrapped
  along the event, so it can pass 9.
OUTCOME  Y = the largest unwrapped phase position reached during the event (continuous).
CLASS  "crossing" = Y >= 6 (reaches the western Pacific beyond the Maritime Continent); secondary
  class Y >= 7.
SAMPLE  events starting November-March (boreal winter; winter labelled by its January year), winters
  1979/80-2024/25. FORCING LEVELS terciles over events of the QBO u50 of the start month: level 0
  easterly (lowest), 2 westerly. Unit: event; uncertainty: 1,999 bootstrap resamples of whole winters.
TESTS (machinery of tc_class_law, validated there for clustered events)
  M1 PRIMARY  fixed two-type mixing of Y across the three QBO levels (raw-moment collinearity T, Wald
     statistic with t on winters - 1 d.f.). Category: "incompatible with fixed two-type mixing"
     (p < 0.05); "compatible" (95% interval of rho inside +-0.25); else "unresolved".
     READING (fixed now): "crossing and blocked are not fixed types whose frequency the QBO changes" if
     incompatible; "consistent with fixed types" if compatible; "not resolved" otherwise.
  M2 SECONDARY  class law for "crossing", easterly to westerly QBO (tc_class_law.law_residual),
     residual read against +-0.25 s.d.; also for Y >= 7.
  M3 SECONDARY  forced share of the crossing-versus-blocked contrast: R = (mean of the level mean F in
     the class - rest) / outcome contrast, against the location-scale prediction R_LS and h2 =
     Var(F)/Var(Y); statistic R - R_LS read against +-0.05.
  M4 DESCRIPTIVE  crossing rate per level; the QBO effect on Y (level 0 minus level 2 mean).
  S1 SENSITIVITY  M1 and M2 on winters 1979/80-2012/13 only (before the RMM processing change).
VALIDATION GATE (--validate; before any data): 1,000 simulations of 46 winters with Poisson(3) events,
  a winter-constant covariate whose terciles define the levels, winter effects s.d. 0.3; fixed types
  blocked = 2.5 + Gamma(2, 0.5), crossing = 6 + Gamma(3, 0.5), crossing weight 1/(1 + exp(-1.2 cov));
  false rejection of M1 at 0.05: one-sided 95% Clopper-Pearson upper bound < 0.07, else M1 is not
  interpreted. Power reported against one population whose location rises with cov.
Output: results/current/8_experiment/mjo_class_law.json
Run:    mjo_class_law.py --validate ; mjo_class_law.py
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "results" / "current" / "8_experiment" / "mjo_class_law.json"
RAW = HERE.parents[0] / "03_data_ingestion" / "raw"
sys.path.insert(0, str(HERE))
import tc_class_law as TC                            # noqa: E402

NAME = "mjo_class_law"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
N_SIM, GATE, TOL, TOL_R = 1000, 0.07, 0.25, 0.05
CROSS, CROSS2 = 6.0, 7.0
WIN = (1980, 2025)                                   # January years of the winters used
INC, COM = TC.INC, TC.COM


# ------------------------------------------------------------------ data
def read_rmm():
    rows = []
    for line in (RAW / "rmm.txt").read_text(errors="replace").splitlines():
        p = line.split()
        if len(p) < 7:
            continue
        try:
            y, m, d = int(p[0]), int(p[1]), int(p[2]); r1, r2, ph, amp = float(p[3]), float(p[4]), int(p[5]), float(p[6])
        except ValueError:
            continue
        if abs(r1) > 900 or abs(r2) > 900 or amp > 900 or ph < 1 or ph > 8:
            continue
        rows.append((y, m, d, r1, r2, ph, amp))
    d = pd.DataFrame(rows, columns=["year", "month", "day", "rmm1", "rmm2", "phase", "amp"])
    d["date"] = pd.to_datetime(d[["year", "month", "day"]])
    d = d[d.date >= "1979-03-01"].set_index("date").asfreq("D")
    return d


def read_qbo():
    vals = {}
    for line in (RAW / "qbo_u50.txt").read_text(errors="replace").splitlines():
        p = line.split()
        if len(p) == 13 and p[0].isdigit() and len(p[0]) == 4:
            for k, v in enumerate(p[1:], start=1):
                try:
                    x = float(v)
                except ValueError:
                    continue
                if abs(x) < 900:
                    vals[(int(p[0]), k)] = x
    return vals


def position(r1, r2):
    th = np.degrees(np.arctan2(r2, r1))
    return np.mod(th + 180.0, 360.0) / 45.0 + 1.0


def events(d):
    amp, pos = d.amp.values, position(d.rmm1.values, d.rmm2.values)
    act = (amp >= 1) & np.isin(np.floor(pos), (2, 3))
    out, last = [], -10 ** 6
    for i in range(1, len(d)):
        if not (act[i] and not act[i - 1]) or i - last < 30 or not np.isfinite(amp[i]):
            continue
        p0 = pos[i]; up = p0; ymax = p0; below = 0; end = i
        for j in range(i + 1, min(i + 61, len(d))):
            if not np.isfinite(amp[j]):
                break
            step = pos[j] - pos[j - 1]
            step = step - 8.0 * np.round(step / 8.0)          # unwrap on the 8-phase circle
            up += step
            if amp[j] >= 1:
                below = 0; end = j; ymax = max(ymax, up)
            else:
                below += 1
                if below >= 3:
                    break
        out.append({"start": d.index[i], "start_pos": p0, "Y": ymax, "days": int(end - i + 1)})
        last = i
    return pd.DataFrame(out)


# ------------------------------------------------------------------ statistics
def forced_share(y, lev, win, n_win, thr, rng, B=1999):
    def stat(w):
        ww = w[win]
        F = np.array([np.sum(ww[lev == k] * y[lev == k]) / ww[lev == k].sum() for k in range(3)])[lev]
        c = y >= thr; mx = lambda v, m: np.sum(ww[m] * v[m]) / ww[m].sum()
        D = mx(y, c) - mx(y, ~c); R = (mx(F, c) - mx(F, ~c)) / D
        Fm, ym = np.sum(ww * F) / ww.sum(), np.sum(ww * y) / ww.sum()
        h2 = np.sum(ww * (F - Fm) ** 2) / np.sum(ww * (y - ym) ** 2)
        # location-scale prediction: each level's events replaced by level 0 mapped by that level's mean and s.d.
        mu = [np.sum(ww[lev == k] * y[lev == k]) / ww[lev == k].sum() for k in range(3)]
        sd = [np.sqrt(np.sum(ww[lev == k] * (y[lev == k] - mu[k]) ** 2) / ww[lev == k].sum()) for k in range(3)]
        ref, wr = y[lev == 0], ww[lev == 0]; pj = [ww[lev == k].sum() / ww.sum() for k in range(3)]
        xs = np.concatenate([mu[k] + sd[k] / sd[0] * (ref - mu[0]) for k in range(3)])
        Fs = np.concatenate([np.full(ref.size, mu[k]) for k in range(3)]); wsy = np.concatenate([pj[k] * wr / wr.sum() for k in range(3)])
        cs = xs >= thr
        Ds = np.sum(wsy[cs] * xs[cs]) / wsy[cs].sum() - np.sum(wsy[~cs] * xs[~cs]) / wsy[~cs].sum()
        RLS = (np.sum(wsy[cs] * Fs[cs]) / wsy[cs].sum() - np.sum(wsy[~cs] * Fs[~cs]) / wsy[~cs].sum()) / Ds
        return R, RLS, h2
    R, RLS, h2 = stat(np.ones(n_win))
    bs = np.array([stat(np.bincount(rng.integers(n_win, size=n_win), minlength=n_win).astype(float)) for _ in range(B)])
    lo, hi = np.nanquantile(bs[:, 0] - bs[:, 1], [0.025, 0.975])
    rd = ("the forced share is what the law predicts" if -TOL_R <= lo and hi <= TOL_R else
          "the class selects forcing beyond the law" if (lo > 0 or hi < 0) and abs(R - RLS) > TOL_R else "inconclusive")
    return {"R_obs": round(float(R), 4), "R_LS": round(float(RLS), 4), "R_h2": round(float(h2), 4),
            "stat_ci95": [round(float(lo), 4), round(float(hi), 4)], "reading": rd}


def category(r):
    ci = r["rho_ci95"]
    return INC if r["p_T"] < 0.05 else COM if ci and -TOL <= ci[0] and ci[1] <= TOL else "unresolved"


def analyse(ev, rng):
    win = ev.win_idx.values; n_win = int(win.max() + 1); y = ev.Y.values; lev = ev.level.values
    m1 = TC.moment_test(y, lev, win, n_win, rng); m1["category"] = category(m1)
    out = {"n_events": int(len(ev)), "n_winters": n_win, "M1_mixing": m1,
           "M2_law_crossing": TC.law_residual(y, lev, win, n_win, CROSS, rng),
           "M2_law_Y7": TC.law_residual(y, lev, win, n_win, CROSS2, rng),
           "M3_forced_share": forced_share(y, lev, win, n_win, CROSS, rng),
           "M4": {"crossing_rate_by_level": [round(float((y[lev == k] >= CROSS).mean()), 4) for k in range(3)],
                  "mean_Y_by_level": [round(float(y[lev == k].mean()), 4) for k in range(3)],
                  "n_by_level": [int((lev == k).sum()) for k in range(3)]}}
    out["M1_reading"] = ("crossing and blocked are not fixed types whose frequency the QBO changes" if m1["category"] == INC else
                         "consistent with fixed types" if m1["category"] == COM else "not resolved")
    return out


# ------------------------------------------------------------------ validation
def simulate(rng, mixing=True, n_win=46):
    cnt = rng.poisson(3, n_win); win = np.repeat(np.arange(n_win), cnt)
    cw = rng.normal(size=n_win); cov = cw[win]
    lev = np.digitize(cov, np.quantile(cov, [1 / 3, 2 / 3]))
    weff = 0.3 * rng.normal(size=n_win)[win]
    if mixing:
        typ = rng.random(win.size) < 1 / (1 + np.exp(-1.2 * cov))
        y = np.where(typ, 6.0 + rng.gamma(3, 0.5, win.size), 2.5 + rng.gamma(2, 0.5, win.size)) + weff
    else:
        y = 4.5 + 0.8 * cov + rng.gamma(3, 0.6, win.size) + weff
    return y, lev, win, n_win


def validate():
    rng = np.random.default_rng(SEED)
    rej = sum(TC.moment_test(*simulate(rng, True), rng=rng)["p_T"] < 0.05 for _ in range(N_SIM))
    pw = sum(TC.moment_test(*simulate(rng, False), rng=rng)["p_T"] < 0.05 for _ in range(N_SIM // 5))
    ub = float(beta.ppf(0.95, rej + 1, N_SIM - rej))
    val = {"n_sim": N_SIM, "false_rejection_rate": rej / N_SIM, "upper95": round(ub, 4), "gate": GATE, "passes": ub < GATE,
           "power_vs_one_population": round(pw / (N_SIM // 5), 4)}
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res["validation"] = val; res["validation_fixed"] = "2026-10-22, before any RMM or QBO value was retrieved"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1), encoding="utf8", newline="\n")
    print(json.dumps(val, indent=1)); return 0


# ------------------------------------------------------------------ main
def main():
    if "--validate" in sys.argv:
        return validate()
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    import platform
    import scipy
    res.update({"registered": "2026-10-22, committed before any RMM or QBO value was retrieved", "seed": SEED, "B": TC.B,
                "inputs": ["03_data_ingestion/raw/rmm.txt", "03_data_ingestion/raw/qbo_u50.txt"], "command": " ".join(sys.argv),
                "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__}})
    d = read_rmm()
    ok = d.amp.notna()
    agree = float((np.floor(position(d.rmm1[ok].values, d.rmm2[ok].values)).astype(int) == d.phase[ok].astype(int).values).mean())
    res["phase_check_agreement"] = round(agree, 4)
    if agree < 0.99:
        res["status"] = "STOPPED: phase convention does not match the file"
        OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n"); print(res["status"], agree); return 1
    ev = events(d)
    ev["month"] = ev.start.dt.month; ev["year"] = ev.start.dt.year
    ev["jyear"] = np.where(ev.month >= 11, ev.year + 1, ev.year)
    qbo = read_qbo(); ev["qbo"] = [qbo.get((y, m), np.nan) for y, m in zip(ev.year, ev.month)]
    res["n_events_all_seasons"] = int(len(ev))
    w = ev[ev.month.isin([11, 12, 1, 2, 3]) & ev.jyear.between(*WIN) & ev.qbo.notna()].copy()
    w["level"] = np.digitize(w.qbo, np.quantile(w.qbo, [1 / 3, 2 / 3]))
    w["win_idx"] = pd.factorize(w.jyear, sort=True)[0]
    res["qbo_tercile_bounds_ms"] = [round(float(x), 2) for x in np.quantile(w.qbo, [1 / 3, 2 / 3])]
    gate_ok = res.get("validation", {}).get("passes")
    rng = np.random.default_rng(SEED + 1)
    res["primary"] = analyse(w, rng)
    if not gate_ok:
        res["primary"]["M1_reading"] = "not interpreted: the validation gate did not pass"
    s = w[w.jyear <= 2013].copy(); s["win_idx"] = pd.factorize(s.jyear, sort=True)[0]
    res["S1_before_2014"] = analyse(s, np.random.default_rng(SEED + 2))
    print(json.dumps({k: res[k] for k in ("phase_check_agreement", "n_events_all_seasons", "qbo_tercile_bounds_ms")}), flush=True)
    print(json.dumps(res["primary"], default=float), flush=True)
    print(json.dumps(res["S1_before_2014"], default=float), flush=True)
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf8", newline="\n")
    print(f"Saved -> {OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
