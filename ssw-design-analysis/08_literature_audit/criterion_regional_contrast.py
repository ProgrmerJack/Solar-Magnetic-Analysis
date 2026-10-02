#!/usr/bin/env python3
"""
criterion_regional_contrast.py
==============================
DOES THE MATCHED NULL REPRODUCE A CLASS DIFFERENCE IN A VARIABLE NOT USED TO
DEFINE THE CLASSES? NORTHERN-EURASIAN TEMPERATURE UNDER THE NAM-BASED DW/NDW LABEL.

Plan approved 2026-10-02 (user: "address all those 5 issues", external review
item 2). Committed BEFORE it was run: no DW-minus-NDW temperature contrast, observed
or null, has been computed by this project. (era5_regional_test.py regressed the
temperature RESIDUAL beyond the circulation on the label; it did not compare the
class contrast with a matched null. The 39 events' temperature anomalies and labels
have been seen before, so this is RETROSPECTIVE.)

QUESTION
  After observed SSWs, events labelled downward (DW) by the published surface
  conditions are colder over northern Eurasia than non-downward (NDW) ones. Is that
  regional contrast what the same classification produces on matched event-free
  dates displaced by the SSW circulation shift -- i.e. selection on the NAM plus
  the ordinary circulation-temperature relation -- or does an SSW-specific regional
  difference between the classes remain?

DATA
  Events: the 39 catalogued SSWs of criterion_transport_cv (ERA5 NAM days 8-52
  inside the record), restricted to those with a complete ERA5 northern-Eurasian
  temperature window (era5_t2m_regions_daily.parquet, to January 2023).
  Classification (primary): published surface conditions, 1000 hPa NAM over days
  8-52 after onset, window mean < 0 and more than half of days < 0 (Karpechko et
  al. 2017 conditions 1-2). Secondary: with condition 3 (150 hPa NAM < 0 on more
  than 70% of days).
  Target (primary, fixed now): northern Eurasia (50-65N, 10-130E) 2 m temperature
  anomaly, mean over days 8-24 after onset (the paper's regional window), anomaly
  from the 1959-2022 31-day-smoothed day-of-year climatology (era5_regional_test).
  Secondary targets: the same over days 8-52; the other three regions, days 8-24.

STATISTIC
  C = mean T(DW) - mean T(NDW) over the events (K).

NULL (primary: whole-winter cross-validation, as criterion_transport_cv)
  For each held-out winter, from the other winters only: event-free candidate days
  (zone-free rule) and the shift of the 1000 and 150 hPa NAM (training SSWs minus
  200 day-of-year-matched event-free sets, criterion_sweep's estimator). Each event
  gets candidates = training-winter event-free days within +-10 days of its day of
  year. One null replicate draws one candidate per event, classifies its NAM
  displaced by the training shift with the identical rule, and takes the same
  date's own temperature anomaly (joint circulation-temperature variability is
  kept; temperature is not shifted, because a constant shift cancels in C).
  4,000 replicates; replicates with fewer than 3 events in either class are
  dropped and counted.
  Sensitivity: the same with the shift estimated from all winters.
TESTS (fixed now)
  p = 2 min(P(C_null <= C_obs), P(C_null >= C_obs)).
  Share reproduced = mean(C_null) / C_obs, with a 4,000-draw bootstrap of events
  (resampled with replacement, labels kept) for the interval of C_obs.
  Power: in each of 1,000 planted replays, one null replicate plays the events, its
  DW members are made colder by 0.5 s.d. of the target on November-March
  event-free dates, and it is tested
  against the null; detection rate at p < 0.05 reported.
READING (fixed now; primary target and classification)
  p >= 0.05: the regional class contrast is reproduced by selection plus the ordinary
    circulation-temperature relation (consequence of the critique beyond the
    defining variable; how much is the share).
  p < 0.05 and C_obs more negative than the null: an SSW-specific regional
    difference between the classes remains (delimits the critique).
  p < 0.05 and C_obs less negative: the classes differ less regionally than the
    null implies.
  Either outcome is reported. A non-significant p with low power is reported as
  inconclusive, not as reproduction.

Output: results/current/9_literature/criterion_regional_contrast.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import criterion_sweep as CS                         # noqa: E402

G, load_catalogue = CS.G, CS.load_catalogue
NAME = "criterion_regional_contrast"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
RESULTS = CS.RESULTS
ING = HERE.parents[0] / "03_data_ingestion"
T_FILE = ING / "era5_t2m_regions_daily.parquet"
REGIONS = ["NEURASIA", "HI_EUROPE", "MID_EASIA", "MID_NAMER"]
N_SHIFT_SETS, N_REP, N_BOOT, N_PLANT, HALFWIN, MIN_CLASS = 200, 4000, 4000, 1000, 10, 3
PLANT_SD = 0.5


def winter(t):
    return t.year if t.month <= 6 else t.year + 1


def t_anomalies():
    t = pd.read_parquet(T_FILE).set_index("date")[REGIONS]
    t.index = pd.to_datetime(t.index)
    base = t[(t.index.year >= 1959) & (t.index.year <= 2022)]
    doy = base.groupby(base.index.dayofyear).mean().reindex(range(1, 367))
    doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    doy.index = range(1, 367)
    return t - doy.reindex(t.index.dayofyear).values


def t_window(tan, dates, a, b):
    """Mean of each region over days a..b after each date (NaN unless complete)."""
    out = np.full((len(dates), len(REGIONS)), np.nan)
    pos = tan.index.get_indexer(pd.DatetimeIndex(dates).normalize())
    v = tan.values
    for i, p in enumerate(pos):
        if p < 0 or p + b >= len(v):
            continue
        seg = v[p + a:p + b + 1]
        if np.isfinite(seg).all():
            out[i] = seg.mean(0)
    return out


def contrast(lab, y):
    ok = np.isfinite(y)
    l, y = lab[ok].astype(bool), y[ok]
    if l.sum() < MIN_CLASS or (~l).sum() < MIN_CLASS:
        return np.nan
    return float(y[l].mean() - y[~l].mean())


def main():
    rng = np.random.default_rng(SEED)
    e = pd.read_parquet(CS.ERA5)
    tan = t_anomalies()
    allev = load_catalogue("primary")
    real = allev[(allev >= e.index.min() + pd.Timedelta(days=70)) & (allev <= e.index.max() - pd.Timedelta(days=70))]
    clean = G.zone_free_index(e.index, allev)
    levels = ("nam_1000", "nam_150")
    Mr = {c: CS.window_matrix(e[c], real) for c in levels}
    Mc = {c: CS.window_matrix(e[c], clean) for c in levels}
    targets = {"d8_24": (8, 24), "d8_52": (8, 52)}
    Tr = {k: t_window(tan, real, *w) for k, w in targets.items()}
    Tc = {k: t_window(tan, clean, *w) for k, w in targets.items()}
    keep = np.isfinite(Tr["d8_24"][:, 0])
    print(f"events with NAM windows {len(real)}; with complete N-Eurasia days 8-24 {int(keep.sum())}", flush=True)
    wr = np.array([winter(t) for t in real]); wc = np.array([winter(t) for t in clean])
    dr = np.array([t.dayofyear for t in real]); dc = np.array([t.dayofyear for t in clean])
    pos = {t: i for i, t in enumerate(clean)}
    clf = {"primary_surface": False, "secondary_with_150hPa": True}

    def shifts(train_ev, train_clean):
        by = G.build_by_doy(clean[train_clean])
        sets = [G.draw_clean(clean[train_clean], dr[train_ev], rng, by) for _ in range(N_SHIFT_SETS)]
        out = {}
        for c in levels:
            pool = np.concatenate([CS.wmean(Mc[c][[pos[t] for t in s_]], 52) for s_ in sets])
            out[c] = float(np.nanmean(CS.wmean(Mr[c][train_ev], 52)) - np.nanmean(pool))
        return out

    # per-event shift and candidate set: whole-winter CV (primary) and all winters (sensitivity)
    ev_idx = np.flatnonzero(keep)
    plan = {"cv": {}, "full": {}}
    full_sh = shifts(np.arange(len(real)), np.arange(len(clean)))
    for w in sorted(set(wr[ev_idx])):
        tr = np.flatnonzero(wr != w); ctr = np.flatnonzero(wc != w)
        sh = shifts(tr, ctr)
        for i in ev_idx[wr[ev_idx] == w]:
            dd = np.abs(dc[ctr] - dr[i]); dd = np.minimum(dd, 366 - dd)
            plan["cv"][i] = (sh, ctr[dd <= HALFWIN])
            dd2 = np.abs(dc - dr[i]); dd2 = np.minimum(dd2, 366 - dd2)
            plan["full"][i] = (full_sh, np.flatnonzero(dd2 <= HALFWIN))
    res = {"plan_approved": "2026-10-02", "seed": SEED, "n_rep": N_REP, "n_boot": N_BOOT, "n_plant": N_PLANT,
           "n_shift_sets": N_SHIFT_SETS, "halfwin_days": HALFWIN, "min_per_class": MIN_CLASS,
           "label": "retrospective (event outcomes and labels seen before; contrast not computed before)",
           "inputs": [str(CS.ERA5.relative_to(CS.ROOT)), str(T_FILE.relative_to(CS.ROOT))],
           "n_events": int(len(ev_idx)), "events": [str(real[i].date()) for i in ev_idx],
           "shift_all_winters": {k: round(v, 4) for k, v in full_sh.items()}}
    ndjfm = np.array([t.month in (11, 12, 1, 2, 3) for t in clean])
    sd_free = {k: np.nanstd(Tc[k][ndjfm], axis=0, ddof=1) for k in targets}

    for cname, c3 in clf.items():
        lab_obs, ok_obs = CS.classify(Mr["nam_1000"][ev_idx], Mr["nam_150"][ev_idx], 52, 0.0, 0.5, c3)
        lab_obs = np.where(ok_obs, lab_obs, False)
        block = {"n_dw": int(lab_obs.sum()), "n_ndw": int((~lab_obs).sum())}
        for variant in ("cv", "full"):
            # null replicates: one candidate per event, classified with the shifted NAM
            draws = np.array([[plan[variant][i][1][rng.integers(len(plan[variant][i][1]))] for i in ev_idx]
                              for _ in range(N_REP)])
            labs = np.zeros(draws.shape, bool)
            for j, i in enumerate(ev_idx):
                sh = plan[variant][i][0]
                l_, o_ = CS.classify(Mc["nam_1000"][draws[:, j]], Mc["nam_150"][draws[:, j]], 52, 0.0, 0.5, c3,
                                     sh["nam_1000"], sh["nam_150"] if c3 else 0.0)
                labs[:, j] = np.where(o_, l_, False)
            out = {}
            for tk in targets:
                for ri, reg in enumerate(REGIONS):
                    if tk == "d8_52" and reg != "NEURASIA":
                        continue
                    y_obs = Tr[tk][ev_idx, ri]
                    c_obs = contrast(lab_obs, y_obs)
                    cn = np.array([contrast(labs[r], Tc[tk][draws[r], ri]) for r in range(N_REP)])
                    dropped = int(np.isnan(cn).sum()); cn = cn[np.isfinite(cn)]
                    p = float(min(1.0, 2 * min(np.mean(cn <= c_obs), np.mean(cn >= c_obs))))
                    rec = {"C_obs_K": round(c_obs, 4), "C_null_mean_K": round(float(cn.mean()), 4),
                           "C_null_ci95_K": [round(float(x), 4) for x in np.percentile(cn, [2.5, 97.5])],
                           "p_two_sided": round(p, 4), "share_reproduced": round(float(cn.mean() / c_obs), 4) if c_obs else None,
                           "n_rep_used": int(len(cn)), "n_rep_dropped_small_class": dropped,
                           "event_free_sd_K": round(float(sd_free[tk][ri]), 4)}
                    if variant == "cv" and tk == "d8_24" and reg == "NEURASIA":
                        # bootstrap of the observed contrast (events resampled, labels kept)
                        bs = []
                        for _ in range(N_BOOT):
                            ix = rng.integers(0, len(y_obs), len(y_obs))
                            bs.append(contrast(lab_obs[ix], y_obs[ix]))
                        bs = np.array(bs)
                        rec["C_obs_ci95_K"] = [round(float(x), 4) for x in np.nanpercentile(bs, [2.5, 97.5])]
                        # power: a null replicate plays the events, DW made colder by PLANT_SD event-free s.d.
                        delta = PLANT_SD * sd_free[tk][ri]
                        det = 0
                        for _ in range(N_PLANT):
                            r = rng.integers(N_REP)
                            yp = Tc[tk][draws[r], ri] - delta * labs[r]
                            cp = contrast(labs[r], yp)
                            if np.isfinite(cp) and min(1.0, 2 * min(np.mean(cn <= cp), np.mean(cn >= cp))) < 0.05:
                                det += 1
                        rec["power_planted"] = {"delta_K": round(float(delta), 4), "detection_rate": round(det / N_PLANT, 4)}
                        if p >= 0.05:
                            rec["reading"] = ("reproduced by selection plus the ordinary circulation-temperature relation"
                                              if det / N_PLANT >= 0.5 else "inconclusive (not rejected, power below 0.5)")
                        else:
                            rec["reading"] = ("SSW-specific regional class difference remains" if c_obs < cn.mean()
                                              else "classes differ less than the null implies")
                    out[f"{reg}_{tk}"] = rec
            block[variant] = out
            print(cname, variant, json.dumps(out.get("NEURASIA_d8_24")), flush=True)
        res[cname] = block
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")


if __name__ == "__main__":
    main()
