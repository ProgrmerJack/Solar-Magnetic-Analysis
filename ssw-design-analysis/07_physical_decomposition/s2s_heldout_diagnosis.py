#!/usr/bin/env python3
"""
s2s_heldout_diagnosis.py
========================
WHY DID THE HELD-OUT REPLICATION FAIL? (exploratory diagnosis, post hoc)

Requested 2026-09-30 after s2s_heldout_test.py (registered in 4c8b3f2) found the
three SSWs of 2023-24 on the upward side of the ensembles (mean rank 0.72 against
0.53) where the 17 SSWs of 1998-2021 had been on the downward side (0.40 against
0.54). Everything here is EXPLORATORY and POST HOC: it can explain the failure, it
cannot rescue the registered test, and it is reported as such. Questions and
readings written before this script was first run:

  D1 DEFECT CHECK (conflict protocol: new work first). With the held-out code path,
     recompute CNRM's P' result on the 1998-2021 events (i) with the original ERA5
     series and (ii) with the spliced series; (i) must reproduce P' exactly.
  D2 VERSIONS OR EVENTS? The P' design (starts 2-9 d before onset, calendar-window
     null, quorum) applied to the 1998-2021 events covered by ECMWF model year 2025
     (hindcasts 2005-2021 of 2005-2024) and CMA model year 2025 (2010-2021), and to
     the same events with ECMWF model year 2022. Reading: new versions below the
     null on the old events -> the change is in the events, not the versions; new
     versions near or above the null -> a version change is a candidate cause.
     Also listed: the P' per-system result against each system's model-version year.
  D2b (added after review of D2): paired event bootstrap (10,000) of the ECMWF
     2025-minus-2022 rank difference on the same events.
  D3 WINTER EFFECT. For CNRM, ECMWF-2025 and CMA-2025, the mean rank of every
     December-March date of each hindcast winter (starts 2-9 d before the date),
     (a) all dates and (b) dates more than 30 d from every catalogued onset.
     Reading: winters 2022/23 and 2023/24 unusually high on dates without SSWs ->
     a winter-level forecast bias (e.g. an end-of-period trend in the leave-one-
     year-out climatology) rather than an SSW-specific change.
  D4 SAMPLING. Event table for all 20 events (multi-system rank, observed and
     ensemble-mean response, observed sign). Mean rank by observed sign. Probability
     that three events drawn from the 17 earlier events' ranks average >= 0.72
     (10,000 draws), and the same under the event-free null.
  D3b/D4b (added after review): correlation of winter-mean ranks (dates > 30 d
     from SSWs) between the three systems over common winters; mean rank by
     observed sign separately for 1998-2021 and the held-out events.
  D5 SIGNAL STRENGTH. Across events, observed response against the multi-system
     ensemble mean (OLS slope, event bootstrap); whether the held-out events lie
     inside the 95% prediction interval of the 17-event regression.

  D6 (added after D1-D5 were seen) EVIDENCE CARRIED BY THE HELD-OUT SIGNS. For each
     event, the forecast probability of a negative outcome (share of members with
     A < 0, mean over starts and systems) and the observed sign. The likelihood of
     the three held-out signs (a) under the forecast probabilities and (b) under
     those probabilities recalibrated by the 1998-2021 under-forecast (a constant
     logit shift fitted by maximum likelihood to the 17 earlier events); their
     ratio says how strongly three events discriminate the two.

Output: results/current/6_predictability/s2s_heldout_diagnosis.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ING = HERE.parents[0] / "03_data_ingestion"
RESULTS = ROOT / "results" / "current" / "6_predictability"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
import s2s_forecast_test as T                        # noqa: E402
import s2s_multimodel_test as MM                     # noqa: E402
import s2s_heldout_test as H                         # noqa: E402
from build_catalogue import load_catalogue          # noqa: E402

NAME = "s2s_heldout_diagnosis"
N_DRAW = 10000
ORIG_OBS = T.OBS
VERSION_YEAR = {"ecmwf": 2022, "eccc": 2025, "cma": 2022, "hmcr": 2025, "kma": 2026,
                "cnrm": 2025, "jma": 2022, "cnr_isac": 2023, "ncep": 2011, "cptec": 2023}


def centre(tag, path):
    return MM.Centre(tag, path)


def main():
    rng = np.random.default_rng(zlib.crc32(NAME.encode()))
    T.__dict__["_orig_load"] = T.load
    T.WIN = (8, 25)
    cat = load_catalogue("primary")
    old = [o for o in cat if o < pd.Timestamp("2022-06-01")]
    res = {"label": "exploratory, post hoc; requested 2026-09-30", "n_draw": N_DRAW}

    # D1 defect check
    T.OBS = ORIG_OBS
    c_orig = centre("cnrm", MM.FILES["cnrm"])
    pm = json.loads((RESULTS / "s2s_multimodel_test.json").read_text())
    ev = [o for o in old if c_orig.at(o) is not None]
    mean_orig = float(np.mean([c_orig.at(o)["pit"] for o in ev]))
    T.OBS = H.SPL_PSL
    c_spl = centre("cnrm", MM.FILES["cnrm"])
    ev2 = [o for o in old if c_spl.at(o) is not None]
    mean_spl = float(np.mean([c_spl.at(o)["pit"] for o in ev2]))
    res["D1"] = {"cnrm_Pprime_H1_mean_pit": pm["centres"]["cnrm"]["H1_mean_pit"],
                 "cnrm_recomputed_original_obs": round(mean_orig, 4), "n_events_original": len(ev),
                 "cnrm_recomputed_spliced_obs": round(mean_spl, 4), "n_events_spliced": len(ev2),
                 "reproduces": bool(abs(mean_orig - pm["centres"]["cnrm"]["H1_mean_pit"]) < 1e-3)}
    print("D1", res["D1"], flush=True)

    # D2 versions or events (spliced obs: the 2025 versions' hindcast years reach 2024)
    sysd = {"ecmwf2025": ING / "s2s_ecmwf2025_psl_cap.parquet", "cma2025": ING / "s2s_cma2025_psl_cap.parquet",
            "ecmwf2022": MM.FILES["ecmwf"], "cnrm": MM.FILES["cnrm"]}
    cen = {k: centre(k, v) for k, v in sysd.items()}
    res["D2"] = {}
    for k, C in cen.items():
        r = H.test({k: C}, old, cat, np.random.default_rng(zlib.crc32(f"{NAME}|D2|{k}".encode())))
        res["D2"][k] = {kk: r[kk] for kk in ("n_events", "mean_rank", "null_mean", "null_q025_q975", "p")} if r else None
        print("D2", k, res["D2"][k], flush=True)
    both = [o for o in old if cen["ecmwf2025"].at(o) is not None and cen["ecmwf2022"].at(o) is not None]
    res["D2"]["ecmwf_same_events"] = {
        "n_events": len(both),
        "mean_rank_2022": round(float(np.mean([cen["ecmwf2022"].at(o)["pit"] for o in both])), 4),
        "mean_rank_2025": round(float(np.mean([cen["ecmwf2025"].at(o)["pit"] for o in both])), 4),
        "per_event": {str(o.date()): [round(cen["ecmwf2022"].at(o)["pit"], 3), round(cen["ecmwf2025"].at(o)["pit"], 3)]
                      for o in both}}
    dd = np.array([cen["ecmwf2025"].at(o)["pit"] - cen["ecmwf2022"].at(o)["pit"] for o in both])
    bsd = np.array([rng.choice(dd, len(dd)).mean() for _ in range(N_DRAW)])
    res["D2"]["ecmwf_same_events"]["diff_2025_minus_2022"] = {
        "mean": round(float(dd.mean()), 4), "ci95": [round(float(q), 4) for q in np.percentile(bsd, [2.5, 97.5])],
        "label": "D2b, added after review"}
    res["D2"]["Pprime_by_version_year"] = {c: {"version_year": VERSION_YEAR[c], "H1_mean_pit": pm["centres"][c].get("H1_mean_pit"),
                                               "null_mean": pm["centres"][c].get("H1_null_mean_pit"), "n_events": pm["centres"][c].get("n_events")}
                                           for c in VERSION_YEAR if c in pm["centres"]}
    print("D2 same events", res["D2"]["ecmwf_same_events"], flush=True)

    # D3 winter effect
    res["D3"] = {}
    for k in ("cnrm", "ecmwf2025", "cma2025"):
        C = cen[k]
        rows = []
        for y in C.hyears:
            for d in pd.date_range(f"{y - 1}-12-01", f"{y}-03-31"):
                r = C.at(d)
                if r is not None:
                    rows.append({"winter": y, "pit": r["pit"],
                                 "far": bool(np.all(np.abs((cat - d).days) > 30))})
        df = pd.DataFrame(rows)
        g_all = df.groupby("winter")["pit"].mean()
        g_far = df[df["far"]].groupby("winter")["pit"].agg(["mean", "count"])
        res["D3"][k] = {"all_dates": {int(w): round(float(v), 3) for w, v in g_all.items()},
                        "dates_30d_from_onsets": {int(w): [round(float(v["mean"]), 3), int(v["count"])] for w, v in g_far.iterrows()},
                        "all_winters_mean_far": round(float(df[df["far"]]["pit"].mean()), 4)}
        print("D3", k, res["D3"][k], flush=True)

    ks = list(res["D3"])
    cor = {}
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            a_, b_ = res["D3"][ks[i]]["dates_30d_from_onsets"], res["D3"][ks[j]]["dates_30d_from_onsets"]
            w = sorted(set(a_) & set(b_))
            cor[f"{ks[i]}~{ks[j]}"] = {"n_winters": len(w), "r": round(float(np.corrcoef(
                [a_[x][0] for x in w], [b_[x][0] for x in w])[0, 1]), 3)}
    res["D3"]["winter_mean_correlations"] = cor
    print("D3 correlations", cor, flush=True)
    # D4 sampling: 20-event table with the multi-system means of the registered tests
    conf = {c: centre(c, MM.FILES[c]) for c in MM.CONFIRM}
    hsys = {k: cen[k] for k in ("cnrm", "ecmwf2025", "cma2025")}
    table = []
    for group, events in ((conf, old), (hsys, H.HELDOUT)):
        for o in events:
            cov = [c for c in group if group[c].at(o) is not None]
            if not cov:
                continue
            table.append({"event": str(o.date()), "set": "1998-2021" if group is conf else "held-out",
                          "n_systems": len(cov),
                          "mean_rank": round(float(np.mean([group[c].at(o)["pit"] for c in cov])), 4),
                          "A_ens": round(float(np.mean([group[c].at(o)["A_ens"] for c in cov])), 1),
                          "A_obs": round(float(group[cov[0]].at(o)["A_obs"]), 1)})
    tb = pd.DataFrame(table)
    early = tb[tb.set == "1998-2021"]
    draws = np.array([rng.choice(early["mean_rank"].values, 3, replace=False).mean() for _ in range(N_DRAW)])
    res["D4"] = {"table": table,
                 "mean_rank_by_observed_sign": {s: {"n": int(len(g)), "mean_rank": round(float(g["mean_rank"].mean()), 4)}
                                                for s, g in (("negative_NAM (A_obs<0)", tb[tb.A_obs < 0]),
                                                             ("positive_NAM (A_obs>=0)", tb[tb.A_obs >= 0]))},
                 "mean_rank_by_sign_and_period": {f"{per}|{sg}": {"n": int(len(g)), "mean_rank": round(float(g["mean_rank"].mean()), 4)}
                                                  for per in ("1998-2021", "held-out")
                                                  for sg, g in (("negative", tb[(tb.set == per) & (tb.A_obs < 0)]),
                                                                ("positive", tb[(tb.set == per) & (tb.A_obs >= 0)]))
                                                  if len(g)},
                 "earlier_share_negative": round(float((early.A_obs < 0).mean()), 3),
                 "heldout_share_negative": round(float((tb[tb.set == 'held-out'].A_obs < 0).mean()), 3),
                 "P_3_of_earlier_ge_0.72": round(float(np.mean(draws >= 0.7167)), 4),
                 "earlier_rank_q025_q975": [round(float(q), 3) for q in np.quantile(early["mean_rank"], [0.025, 0.975])]}
    print("D4", {k: v for k, v in res["D4"].items() if k != "table"}, flush=True)
    for row in table:
        print("   ", row, flush=True)

    # D5 signal strength
    x, yv = early["A_ens"].values, early["A_obs"].values
    b, a = np.polyfit(x, yv, 1)
    resid_sd = float(np.std(yv - (a + b * x), ddof=2))
    bs = []
    for _ in range(N_DRAW):
        i = rng.integers(0, len(x), len(x))
        if np.std(x[i]) > 0:
            bs.append(np.polyfit(x[i], yv[i], 1)[0])
    ho = tb[tb.set == "held-out"]
    pred = []
    for _, r in ho.iterrows():
        yhat = a + b * r["A_ens"]
        se = resid_sd * np.sqrt(1 + 1 / len(x) + (r["A_ens"] - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
        pred.append({"event": r["event"], "A_ens": r["A_ens"], "A_obs": r["A_obs"], "predicted": round(float(yhat), 1),
                     "pi95": [round(float(yhat - 2.16 * se), 1), round(float(yhat + 2.16 * se), 1)],
                     "inside": bool(abs(r["A_obs"] - yhat) <= 2.16 * se)})
    res["D5"] = {"slope_obs_on_ens_1998_2021": round(float(b), 3), "intercept_Pa": round(float(a), 1),
                 "slope_ci95": [round(float(q), 3) for q in np.percentile(bs, [2.5, 97.5])],
                 "residual_sd_Pa": round(resid_sd, 1), "heldout_vs_prediction": pred}
    print("D5", res["D5"], flush=True)
    # D6 evidence carried by the held-out signs
    from scipy.optimize import minimize_scalar
    from scipy.special import expit, logit

    def pneg(group, o):
        cov = [c for c in group if group[c].at(o) is not None]
        return float(np.mean([np.mean([np.mean(np.asarray(m) < 0) for m in group[c].at(o)["members"]]) for c in cov]))
    pe = np.clip(np.array([pneg(conf, pd.Timestamp(e)) for e in early["event"]]), 0.01, 0.99)
    ye = (early["A_obs"].values < 0).astype(float)
    nll = lambda d: -np.sum(ye * np.log(expit(logit(pe) + d)) + (1 - ye) * np.log(1 - expit(logit(pe) + d)))
    dlt = float(minimize_scalar(nll, bounds=(-5, 5), method="bounded").x)
    ph = np.clip(np.array([pneg(hsys, o) for o in H.HELDOUT]), 0.01, 0.99)
    yh = (ho["A_obs"].values < 0).astype(float)
    lik = lambda p: float(np.prod(np.where(yh == 1, p, 1 - p)))
    la, lb = lik(ph), lik(expit(logit(ph) + dlt))
    res["D6"] = {"label": "added after D1-D5 were seen",
                 "earlier_mean_forecast_p_negative": round(float(pe.mean()), 3),
                 "earlier_observed_share_negative": round(float(ye.mean()), 3),
                 "logit_shift_fitted_on_earlier": round(dlt, 3),
                 "heldout_forecast_p_negative": [round(float(v), 3) for v in ph],
                 "heldout_recalibrated_p_negative": [round(float(v), 3) for v in expit(logit(ph) + dlt)],
                 "heldout_observed_negative": [bool(v) for v in yh],
                 "likelihood_forecast_as_issued": round(la, 4), "likelihood_recalibrated": round(lb, 4),
                 "ratio_issued_over_recalibrated": round(la / lb, 2)}
    print("D6", res["D6"], flush=True)
    (RESULTS / "s2s_heldout_diagnosis.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> s2s_heldout_diagnosis.json")


if __name__ == "__main__":
    main()
