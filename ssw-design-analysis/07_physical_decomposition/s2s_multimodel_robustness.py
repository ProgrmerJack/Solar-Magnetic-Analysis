#!/usr/bin/env python3
"""
s2s_multimodel_robustness.py
============================
VERIFICATION DETAILS A REFEREE WILL ASK FOR, ON THE P' DATA AND DESIGN.

Plan approved 2026-09-29. Written before these were computed. Same systems,
anomalies, starts (2-9 d before onset), events and calendar-window null as
s2s_multimodel_test.py (P'); outcomes: polar-cap NAM proxy ("psl") and northern
Eurasian 2 m temperature ("NEURASIA", s2s_regional_test.loader).

  R1 RELIABILITY OF THE PROBABILITY FORECAST the paper recommends:
     p = P(outcome < 0) = share of members with a negative window mean, averaged
     over a system's starts, then over the nine confirmatory systems; y = 1 if the
     observed outcome is negative. Brier score and its Murphy decomposition
     (reliability, resolution, uncertainty; bins of width 0.2) after SSWs and on
     every zone-free calendar-window date (the P' candidate pool); mean forecast
     probability against observed frequency in both sets.
  R2 DISCRIMINATION, QUANTIFIED: delta r = r(SSWs) - mean r of the conditional
     null (P' D); 95% interval of r(SSWs) from 10,000 bootstrap resamples of events;
     the minimum detectable delta r at alpha 0.05 with 80% power, taking the
     conditional null's spread as the sampling spread (MDE = q95 - mean + 0.84 sd);
     power at delta r = 0.2 and 0.3 under a shift of the null.
  R3 EVENTS AS THE UNIT: a linear mixed model of the per-(system, date) rank,
     rank ~ SSW + system (fixed) + (1 | date), on the events and up to 2,000
     randomly drawn candidate dates (all of them when fewer exist: 971 for the
     polar cap, 1,495 for temperature); the SSW coefficient with its 95% interval.
  R4 WELL-COVERED EVENTS ONLY: H1 (mean rank below the null, P' design) on the
     events covered by at least five confirmatory systems.

Output: results/current/6_predictability/s2s_multimodel_robustness.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
import s2s_forecast_test as T                        # noqa: E402
import s2s_multimodel_test as MM                     # noqa: E402
import s2s_regional_test as SR                       # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

NAME = "s2s_multimodel_robustness"
N_NULL = 10000
N_BOOT = 10000
N_MIXED = 2000
MIN_SYSTEMS = 5


def brier_decomp(p, y, bins=np.linspace(0, 1, 6)):
    p, y = np.asarray(p, float), np.asarray(y, float)
    bs = float(np.mean((p - y) ** 2)); ybar = float(y.mean())
    k = np.clip(np.digitize(p, bins[1:-1]), 0, len(bins) - 2)
    rel = res = 0.0; table = []
    for b in range(len(bins) - 1):
        i = k == b
        if i.any():
            pk, ok = p[i].mean(), y[i].mean()
            rel += i.sum() * (pk - ok) ** 2; res += i.sum() * (ok - ybar) ** 2
            table.append({"bin": [round(bins[b], 1), round(bins[b + 1], 1)], "n": int(i.sum()),
                          "mean_p": round(float(pk), 3), "obs_freq": round(float(ok), 3)})
    n = len(p)
    return {"n": n, "brier": round(bs, 4), "reliability": round(rel / n, 4),
            "resolution": round(res / n, 4), "uncertainty": round(ybar * (1 - ybar), 4),
            "mean_p": round(float(p.mean()), 4), "obs_freq": round(ybar, 4), "bins": table}


def run(outcome):
    rng = np.random.default_rng(zlib.crc32(f"{NAME}|{outcome}".encode()))
    if outcome == "psl":
        T.WIN, T.load = (8, 25), T.__dict__["_orig_load"]
        files = {c: f for c, f in MM.FILES.items() if c in MM.CONFIRM}
    else:
        T.WIN, T.load = SR.WIN_T, SR.loader(outcome)
        files = {c: f for c, f in SR.FILES.items() if c in MM.CONFIRM}
    cen = {c: MM.Centre(c, f) for c, f in files.items() if f.exists()}
    cat = load_catalogue("primary")

    def zone_free(p):
        return bool(np.all(np.abs((cat - p).days) > T.ZONE_SEP))

    real = {c: {o: C.at(o) for o in cat if C.at(o) is not None} for c, C in cen.items()}
    events = sorted({o for c in cen for o in real[c]})
    cover = {o: [c for c in cen if o in real[c]] for o in events}
    years = sorted({y for C in cen.values() for y in C.hyears} | {y + 1 for C in cen.values() for y in C.hyears})
    cands = {o: [p for y in years for dd in range(-MM.NULL_WINDOW, MM.NULL_WINDOW + 1)
                 for p in [o + pd.DateOffset(years=y - o.year) + pd.Timedelta(days=dd)] if zone_free(p)]
             for o in events}
    need = {o: int(np.ceil(len(cover[o]) / 2)) for o in events}
    pools = {o: [p for p in cands[o] if sum(cen[c].at(p) is not None for c in cover[o]) >= need[o]]
             for o in events}

    def prob(r):
        return float(np.mean([np.mean(np.asarray(m) < 0) for m in r["members"]]))

    def mm(o_or_p, systems):
        v = [cen[c].at(o_or_p) for c in systems]
        v = [x for x in v if x is not None]
        if not v:
            return None
        return {"pit": float(np.mean([x["pit"] for x in v])), "A_ens": float(np.mean([x["A_ens"] for x in v])),
                "A_obs": float(np.mean([x["A_obs"] for x in v])), "p": float(np.mean([prob(x) for x in v])),
                "per": {c: x["pit"] for c, x in zip([c for c in systems if cen[c].at(o_or_p) is not None], v)}}

    ev = {o: mm(o, cover[o]) for o in events}
    # --- R1 reliability
    all_c = sorted({p for o in events for p in pools[o]})
    free = {p: mm(p, list(cen)) for p in all_c}
    free = {p: v for p, v in free.items() if v is not None}
    r1 = {"after_SSWs": brier_decomp([ev[o]["p"] for o in events], [ev[o]["A_obs"] < 0 for o in events]),
          "event_free_dates": brier_decomp([v["p"] for v in free.values()], [v["A_obs"] < 0 for v in free.values()])}
    # --- R2 discrimination
    a = np.array([ev[o]["A_ens"] for o in events]); ob = np.array([ev[o]["A_obs"] for o in events])
    r_real = float(np.corrcoef(a, ob)[0, 1]); var_real = float(np.var(ob))
    nr, nv = [], []
    for _ in range(N_NULL):
        rr = [mm(pools[o][rng.integers(len(pools[o]))], cover[o]) for o in events]
        aa = np.array([x["A_ens"] for x in rr]); oo = np.array([x["A_obs"] for x in rr])
        nr.append(np.corrcoef(aa, oo)[0, 1]); nv.append(np.var(oo))
    nr, nv = np.array(nr), np.array(nv)
    mv = np.abs(np.log(nv / var_real)) < np.log(1.25)
    rc = nr[mv]
    boot = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(events), len(events))
        boot.append(np.corrcoef(a[i], ob[i])[0, 1])
    sd0, mu0, q95 = float(np.std(rc)), float(np.mean(rc)), float(np.quantile(rc, 0.95))
    r2 = {"r_ssw": round(r_real, 4), "r_ssw_ci95_event_bootstrap": [round(float(q), 4) for q in np.nanquantile(boot, [0.025, 0.975])],
          "null_cond_mean": round(mu0, 4), "null_cond_sd": round(sd0, 4), "n_null_cond": int(mv.sum()),
          "delta_r": round(r_real - mu0, 4),
          "delta_r_ci95": [round(float(q) - mu0, 4) for q in np.nanquantile(boot, [0.025, 0.975])],
          "mde_alpha05_power80": round(q95 - mu0 + 0.84 * sd0, 4),
          "power_at_delta_0.2": round(float(np.mean(rc + 0.2 >= q95)), 3),
          "power_at_delta_0.3": round(float(np.mean(rc + 0.3 >= q95)), 3)}
    # --- R3 mixed model
    import statsmodels.formula.api as smf
    rows = []
    for o in events:
        for c, v in ev[o]["per"].items():
            rows.append({"date": f"e{o.date()}", "system": c, "ssw": 1, "rank": v})
    pick = rng.choice(len(all_c), min(N_MIXED, len(all_c)), replace=False)
    for k in pick:
        p = all_c[k]
        if p in free:
            for c, v in free[p]["per"].items():
                rows.append({"date": f"p{p.date()}", "system": c, "ssw": 0, "rank": v})
    df = pd.DataFrame(rows)
    fit = smf.mixedlm("rank ~ ssw + C(system)", df, groups=df["date"]).fit(method="lbfgs")
    ci = fit.conf_int().loc["ssw"].values
    r3 = {"n_rows": int(len(df)), "n_event_dates": len(events), "n_free_dates": int(df[df.ssw == 0].date.nunique()),
          "ssw_coef": round(float(fit.params["ssw"]), 4), "ssw_ci95": [round(float(x), 4) for x in ci],
          "ssw_p": round(float(fit.pvalues["ssw"]), 5)}
    # --- R4 well-covered events
    wc = [o for o in events if len(cover[o]) >= MIN_SYSTEMS]
    obs_pit = float(np.mean([ev[o]["pit"] for o in wc]))
    npit = np.array([np.mean([mm(pools[o][rng.integers(len(pools[o]))], cover[o])["pit"] for o in wc])
                     for _ in range(N_NULL)])
    r4 = {"min_systems": MIN_SYSTEMS, "n_events": len(wc), "mean_rank": round(obs_pit, 4),
          "null_mean": round(float(npit.mean()), 4), "p": round(float(np.mean(npit <= obs_pit)), 4)}
    return {"n_events": len(events), "systems": list(cen), "R1_reliability": r1,
            "R2_discrimination": r2, "R3_mixed_model": r3, "R4_well_covered": r4}


def main():
    T.__dict__["_orig_load"] = T.load
    res = {"plan_approved": "2026-09-29", "n_null": N_NULL, "n_boot": N_BOOT, "n_mixed_dates": N_MIXED,
           "seed_base": NAME, "outcomes": {}}
    for oc in ("psl", "NEURASIA"):
        r = run(oc)
        res["outcomes"][oc] = r
        a, f = r["R1_reliability"]["after_SSWs"], r["R1_reliability"]["event_free_dates"]
        print(f"{oc:9s} R1 after SSWs: mean p {a['mean_p']} obs {a['obs_freq']} Brier {a['brier']} rel {a['reliability']} | "
              f"event-free: mean p {f['mean_p']} obs {f['obs_freq']} Brier {f['brier']} rel {f['reliability']}", flush=True)
        print(f"{'':9s} R2 {r['R2_discrimination']}", flush=True)
        print(f"{'':9s} R3 {r['R3_mixed_model']}", flush=True)
        print(f"{'':9s} R4 {r['R4_well_covered']}", flush=True)
    (RESULTS / "s2s_multimodel_robustness.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> s2s_multimodel_robustness.json")


if __name__ == "__main__":
    main()
