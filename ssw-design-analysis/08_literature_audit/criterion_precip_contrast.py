#!/usr/bin/env python3
"""
criterion_precip_contrast.py
============================
A NEW MEASURED ENDPOINT: STATION-BASED PRECIPITATION AFTER OBSERVED SSWs. DOES THE
MATCHED NULL REPRODUCE THE DOWNWARD/NON-DOWNWARD DIFFERENCE IN RAINFALL, AND DOES THE
SHIFTED NULL PREDICT HELD-OUT WINTERS' RAINFALL?

Registered 2026-10-08, committed BEFORE any E-OBS value was fetched
(03_data_ingestion/acquire_eobs_precip.py). No precipitation outcome of any SSW has
been examined by this project; the NAM labels of the 39 events have (RETROSPECTIVE in
that sense, as criterion_regional_contrast).

EVENTS, LABEL, NULL: exactly criterion_regional_contrast (39 catalogued SSWs with ERA5
  NAM windows; published surface conditions 1-2 on the 1000 hPa NAM, days 8-52; the
  whole-winter cross-validated matched null: event-free candidates within +-10 days of
  each event's day of year from the other winters, NAM displaced by the training shift
  and classified identically, the candidate date's OWN precipitation taken; 4,000
  replicates), restricted to events whose target window is complete in E-OBS.
TARGETS (boxes and window from Dai et al. 2025, fixed now)
  PRIMARY    IBERIA 36-44N, 10-1W, mean daily precipitation anomaly over days 1-25
             after onset (mm/day; anomaly from the 1959-2022 31-day-smoothed
             day-of-year climatology)
  SECONDARY  SCANDINAVIA 55-70N, 4-20E and BRITISH_ISLES 49-61N, 10W-1E, days 1-25;
             all three over days 8-52; PERSISTENCE: longest run of wet days (box mean
             >= 1 mm/day) within days 1-25, IBERIA and SCANDINAVIA.
STATISTIC  C = mean(DW) - mean(NDW).
TESTS (as criterion_regional_contrast)  p = 2 min(P(C_null <= C_obs), P(C_null >= C_obs));
  share reproduced = mean(C_null) / C_obs; bootstrap interval of C_obs (4,000);
  power (primary): DW members of a null replicate made WETTER by 0.5 event-free s.d.
  (November-March), detection rate at p < 0.05 over 1,000 replays.
READING (primary, fixed now)
  p >= 0.05 and power >= 0.5: "the rainfall class difference is reproduced by selection
    plus the ordinary circulation-precipitation relation";
  p >= 0.05 and power < 0.5: "inconclusive";
  p < 0.05: "an SSW-specific rainfall difference between the classes remains" with its
    sign relative to the null.
FORECAST (secondary, fixed now): whole-winter CV. For each event, the shifted-null
  forecast of IBERIA (and SCANDINAVIA) days 1-25 precipitation = its candidates' values
  plus the training shift (training SSWs' mean minus the mean of 200 day-of-year-matched
  event-free sets, as the NAM shift); climatology = candidates unshifted. CRPS gain
  (climatology minus shifted), mean over events, 95% interval from resampling winters.
  Reading: interval above zero "improves on climatology"; below "worse"; else "not
  resolved".

Output: results/current/9_literature/criterion_precip_contrast.json
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
import criterion_regional_contrast as CRC            # noqa: E402  winter, contrast, constants

G, load_catalogue = CS.G, CS.load_catalogue
NAME = "criterion_precip_contrast"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)
RESULTS = CS.RESULTS
P_FILE = HERE.parents[0] / "03_data_ingestion" / "eobs_precip_regions_daily.parquet"
NULL_NPZ = HERE.parents[0] / "03_data_ingestion" / "_null_replicates" / "criterion_precip_contrast.npz"   # read by equivalence_observational
REGIONS = ["IBERIA", "SCANDINAVIA", "BRITISH_ISLES"]
N_SHIFT_SETS, N_REP, N_BOOT, N_PLANT, HALFWIN = CRC.N_SHIFT_SETS, CRC.N_REP, CRC.N_BOOT, CRC.N_PLANT, CRC.HALFWIN
PLANT_SD, WET = 0.5, 1.0


def load():
    p = pd.read_parquet(P_FILE).set_index("date")[REGIONS]
    p.index = pd.to_datetime(p.index)
    base = p[(p.index.year >= 1959) & (p.index.year <= 2022)]
    doy = base.groupby(base.index.dayofyear).mean().reindex(range(1, 367))
    doy = pd.concat([doy.iloc[-15:], doy, doy.iloc[:15]]).rolling(31, center=True, min_periods=15).mean().iloc[15:-15]
    doy.index = range(1, 367)
    return p, p - doy.reindex(p.index.dayofyear).values


def window(frame, dates, a, b, stat="mean"):
    out = np.full((len(dates), len(REGIONS)), np.nan)
    for i, d in enumerate(pd.DatetimeIndex(dates).normalize()):
        seg = frame.reindex(pd.date_range(d + pd.Timedelta(days=a), d + pd.Timedelta(days=b))).values
        if np.isfinite(seg).all():
            if stat == "mean":
                out[i] = seg.mean(0)
            else:                                         # longest run of wet days
                for r in range(seg.shape[1]):
                    best = run = 0
                    for v in seg[:, r] >= WET:
                        run = run + 1 if v else 0; best = max(best, run)
                    out[i, r] = best
    return out


def crps_ens(x, y):
    x = np.asarray(x, float)
    return float(np.mean(np.abs(x - y)) - 0.5 * np.mean(np.abs(x[:, None] - x[None, :])))


def main():
    rng = np.random.default_rng(SEED)
    e = pd.read_parquet(CS.ERA5)
    raw, pan = load()
    allev = load_catalogue("primary")
    real = allev[(allev >= e.index.min() + pd.Timedelta(days=70)) & (allev <= e.index.max() - pd.Timedelta(days=70))]
    clean = G.zone_free_index(e.index, allev)
    levels = ("nam_1000", "nam_150")
    Mr = {c: CS.window_matrix(e[c], real) for c in levels}
    Mc = {c: CS.window_matrix(e[c], clean) for c in levels}
    targets = {"d1_25": (1, 25), "d8_52": (8, 52)}
    Yr = {k: window(pan, real, *w) for k, w in targets.items()}
    Yc = {k: window(pan, clean, *w) for k, w in targets.items()}
    Yr["wetrun_d1_25"] = window(raw, real, 1, 25, "run"); Yc["wetrun_d1_25"] = window(raw, clean, 1, 25, "run")
    keep = np.isfinite(Yr["d1_25"][:, 0])
    print(f"events with NAM windows {len(real)}; with complete Iberian days 1-25 {int(keep.sum())}", flush=True)
    wr = np.array([CRC.winter(t) for t in real]); wc = np.array([CRC.winter(t) for t in clean])
    dr = np.array([t.dayofyear for t in real]); dc = np.array([t.dayofyear for t in clean])
    pos = {t: i for i, t in enumerate(clean)}

    def shifts(train_ev, train_clean):
        by = G.build_by_doy(clean[train_clean])
        sets = [G.draw_clean(clean[train_clean], dr[train_ev], rng, by) for _ in range(N_SHIFT_SETS)]
        out = {}
        for c in levels:
            pool = np.concatenate([CS.wmean(Mc[c][[pos[t] for t in s_]], 52) for s_ in sets])
            out[c] = float(np.nanmean(CS.wmean(Mr[c][train_ev], 52)) - np.nanmean(pool))
        for ri, reg in enumerate(REGIONS):                    # precipitation shift (forecast only)
            pool = np.concatenate([Yc["d1_25"][[pos[t] for t in s_], ri] for s_ in sets])
            out[f"pr_{reg}"] = float(np.nanmean(Yr["d1_25"][train_ev, ri]) - np.nanmean(pool))
        return out

    ev_idx = np.flatnonzero(keep)
    plan = {}
    for w in sorted(set(wr[ev_idx])):
        tr = np.flatnonzero((wr != w) & keep); ctr = np.flatnonzero(wc != w)
        sh = shifts(tr, ctr)
        for i in ev_idx[wr[ev_idx] == w]:
            dd = np.abs(dc[ctr] - dr[i]); dd = np.minimum(dd, 366 - dd)
            plan[i] = (sh, ctr[dd <= HALFWIN])
    res = {"registered": "2026-10-08, committed before any E-OBS value was fetched", "seed": SEED, "n_rep": N_REP,
           "n_boot": N_BOOT, "n_plant": N_PLANT, "n_shift_sets": N_SHIFT_SETS, "halfwin_days": HALFWIN,
           "inputs": [str(CS.ERA5.relative_to(CS.ROOT)), str(P_FILE.relative_to(CS.ROOT)) if P_FILE.is_relative_to(CS.ROOT) else str(P_FILE)],
           "n_events": int(len(ev_idx)), "events": [str(real[i].date()) for i in ev_idx]}
    lab_obs, ok_obs = CS.classify(Mr["nam_1000"][ev_idx], Mr["nam_150"][ev_idx], 52, 0.0, 0.5, False)
    lab_obs = np.where(ok_obs, lab_obs, False)
    res["n_dw"], res["n_ndw"] = int(lab_obs.sum()), int((~lab_obs).sum())
    draws = np.array([[plan[i][1][rng.integers(len(plan[i][1]))] for i in ev_idx] for _ in range(N_REP)])
    labs = np.zeros(draws.shape, bool)
    for j, i in enumerate(ev_idx):
        l_, o_ = CS.classify(Mc["nam_1000"][draws[:, j]], Mc["nam_150"][draws[:, j]], 52, 0.0, 0.5, False,
                             plan[i][0]["nam_1000"], 0.0)
        labs[:, j] = np.where(o_, l_, False)
    ndjfm = np.array([t.month in (11, 12, 1, 2, 3) for t in clean])
    out = {}; store = {}
    for tk in ("d1_25", "d8_52", "wetrun_d1_25"):
        for ri, reg in enumerate(REGIONS):
            if tk == "wetrun_d1_25" and reg == "BRITISH_ISLES":
                continue
            y_obs = Yr[tk][ev_idx, ri]
            c_obs = CRC.contrast(lab_obs, y_obs)
            cn = np.array([CRC.contrast(labs[r], Yc[tk][draws[r], ri]) for r in range(N_REP)])
            dropped = int(np.isnan(cn).sum()); cn = cn[np.isfinite(cn)]
            store[f"{reg}_{tk}"] = cn
            p = float(min(1.0, 2 * min(np.mean(cn <= c_obs), np.mean(cn >= c_obs))))
            sd = float(np.nanstd(Yc[tk][ndjfm, ri], ddof=1))
            rec = {"C_obs": round(c_obs, 4), "C_null_mean": round(float(cn.mean()), 4),
                   "C_null_ci95": [round(float(x), 4) for x in np.percentile(cn, [2.5, 97.5])],
                   "p_two_sided": round(p, 4), "share_reproduced": round(float(cn.mean() / c_obs), 4) if c_obs else None,
                   "n_rep_used": int(len(cn)), "n_rep_dropped_small_class": dropped, "event_free_sd": round(sd, 4)}
            if tk == "d1_25" and reg == "IBERIA":
                bs = []
                for _ in range(N_BOOT):
                    ix = rng.integers(0, len(y_obs), len(y_obs)); bs.append(CRC.contrast(lab_obs[ix], y_obs[ix]))
                rec["C_obs_ci95"] = [round(float(x), 4) for x in np.nanpercentile(bs, [2.5, 97.5])]
                delta = PLANT_SD * sd; det = 0
                for _ in range(N_PLANT):
                    r = rng.integers(N_REP)
                    cp = CRC.contrast(labs[r], Yc[tk][draws[r], ri] + delta * labs[r])
                    if np.isfinite(cp) and min(1.0, 2 * min(np.mean(cn <= cp), np.mean(cn >= cp))) < 0.05:
                        det += 1
                rec["power_planted"] = {"delta": round(delta, 4), "detection_rate": round(det / N_PLANT, 4)}
                rec["reading"] = ((("the rainfall class difference is reproduced by selection plus the ordinary "
                                    "circulation-precipitation relation") if det / N_PLANT >= 0.5 else "inconclusive")
                                  if p >= 0.05 else
                                  ("an SSW-specific rainfall difference between the classes remains ("
                                   + ("DW wetter" if c_obs > cn.mean() else "DW drier") + " than the null)"))
            out[f"{reg}_{tk}"] = rec
            print(reg, tk, rec, flush=True)
    res["contrast"] = out
    fc = {}
    for reg in ("IBERIA", "SCANDINAVIA"):
        ri = REGIONS.index(reg); gains, wins = [], []
        for i in ev_idx:
            y = Yr["d1_25"][i, ri]; cand = Yc["d1_25"][plan[i][1], ri]; cand = cand[np.isfinite(cand)]
            if not np.isfinite(y) or len(cand) < 20:
                continue
            gains.append(crps_ens(cand, y) - crps_ens(cand + plan[i][0][f"pr_{reg}"], y)); wins.append(wr[i])
        gains, wins = np.array(gains), np.array(wins); uw = np.unique(wins); bs = []
        for _ in range(N_BOOT):
            pick = rng.choice(uw, len(uw)); bs.append(np.concatenate([gains[wins == w] for w in pick]).mean())
        ci = [round(float(x), 4) for x in np.percentile(bs, [2.5, 97.5])]
        fc[reg] = {"mean_crps_gain": round(float(gains.mean()), 4), "ci95": ci, "n_events": int(len(gains)),
                   "reading": "improves on climatology" if ci[0] > 0 else "worse than climatology" if ci[1] < 0 else "not resolved"}
        print("forecast", reg, fc[reg], flush=True)
    res["forecast_d1_25"] = fc
    NULL_NPZ.parent.mkdir(parents=True, exist_ok=True); np.savez_compressed(NULL_NPZ, **store)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"{NAME}.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print(f"Saved -> {NAME}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
