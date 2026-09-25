#!/usr/bin/env python3
"""
s2s_forecast_test.py
====================
DO OPERATIONAL FORECASTS STARTED BEFORE AN SSW KNOW WHICH EVENTS WILL COUPLE?

Plan approved 2026-09-25 (before any data were downloaded).

QUESTION
  ECMWF S2S reforecasts (model year 2022, hindcasts 2002-2021, 11 members) started
  2-9 days ("short") and 10-17 days ("mid") before an observed SSW:
  (a) RELIABILITY -- does the observed polar-cap outcome over days +8..+25 after
      onset fall inside the forecast ensemble as often as it should, i.e. is the
      rank histogram flat, as at ordinary winter starts?
  (b) DISCRIMINATION -- across events, does the forecast say which events end up
      downward-propagating (Karpechko et al. 2017 conditions 1-2) and how strongly,
      better than the same forecasts do for event-free winter dates?
  The one-shifted-population reading predicts (a) yes and (b) no.

INPUTS
  03_data_ingestion/s2s_ecmf_psl_cap.parquet   (acquire_s2s_reforecasts.py)
  03_data_ingestion/era5_psl_cap_6h.parquet    (00 UTC values; same cap reduction)
  events from load_catalogue("primary") only.

METHOD
  Anomalies: forecast minus the leave-one-year-out mean of the other 19 hindcast
  years (all members) at the same model date and lead; observation minus the
  leave-one-year-out mean of ERA5 at the same calendar start and lead. Both in Pa,
  so model bias and drift cancel. NAM sign: A = -anomaly (A < 0 = downward).
  Window: valid days onset+8 .. onset+25 (18 daily values, 00 UTC).
  Classification: DW = window mean A < 0 and more than half of days A < 0.
  Event level: an event's starts within a lead bin are averaged (the observed
  window is the same for all of them).
  Discrimination: r(ensemble-mean A, observed A) across events; Brier skill of
  P(DW) = member fraction against the constant forecast mean(P).
  Null: N_NULL draws; each event is replaced by a pseudo-onset on a zone-free
  day (more than 135 d from every catalogued onset) within +-21 calendar days
  of the event's date in any hindcast year, with that pseudo-onset's own starts
  in the same lead bin, gathered exactly as for the event; the statistics are
  recomputed. p = P(null >= observed). The plan's first design moved each event
  to another year at the same calendar date; with only 8 event-free winters it
  failed the known-truth calibration (NULL_MODE switch kept for that check).
  Reliability: rank (PIT) of the observed window mean among the 11 members;
  flatness of the 12-bin rank histogram (Monte-Carlo chi-square), outer-bin share
  (U shape, exact binomial against 2/12) and spread-error ratio, at SSW starts and
  at zone-free control starts.

WHAT WOULD SHOW THE THESIS WRONG
  discrimination at SSW starts above the null at p < 0.05, or a U-shaped rank
  histogram at SSW starts (outcomes in both tails: a second population the
  forecast misses).

Output: results/current/6_predictability/s2s_forecast_test.json
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
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
from build_catalogue import load_catalogue          # noqa: E402

FC = ING / "s2s_ecmf_psl_cap.parquet"
OBS = ING / "era5_psl_cap_6h.parquet"
WIN = (8, 25)
BINS = {"short": (2, 9), "mid": (10, 17)}
ZONE_SEP = 135                       # |pseudo - real| > 75 - (-60), as zone_free_index
N_NULL = 10000
NULL_MODE = "calendar_window"        # or "year_shift"; chosen by the calibration check
NULL_WINDOW = 21                     # calendar_window: pseudo-onsets within +-21 d of the event's date
NAME = "s2s_forecast_test"
SEED = zlib.crc32(NAME.encode()) % (2 ** 31)


def load():
    fc = pd.read_parquet(FC)
    fc["init"] = pd.to_datetime(fc["init"])
    fc["hyear"] = fc["init"].dt.year
    fc["md"] = pd.to_datetime(fc["model_date"]).dt.strftime("%m-%d")
    # leave-one-year-out forecast climatology at (model date, lead)
    g = fc.groupby(["md", "lead_day", "hyear"])["psl_cap_N"].agg(["sum", "count"])
    tot = g.groupby(level=[0, 1]).sum()
    loo = (tot.reindex(g.index.droplevel(2)).values - g.values)
    clim = pd.Series(loo[:, 0] / loo[:, 1], index=g.index)
    fc["anom"] = fc["psl_cap_N"].values - clim.reindex(
        pd.MultiIndex.from_frame(fc[["md", "lead_day", "hyear"]])).values
    ob = pd.read_parquet(OBS)
    ob["time"] = pd.to_datetime(ob["time"])
    ob = ob[ob["time"].dt.hour == 0].set_index("time")["psl_cap_N"]
    if not ob.index.is_unique:           # .get would return several values
        raise ValueError(f"{OBS.name}: repeated timestamps")
    return fc, ob


def obs_anom(ob, anchors, offsets, hyears):
    """Observed anomaly at anchor + offset days, against the leave-one-year-out
    ERA5 mean at the same calendar dates in the other hindcast years (the whole
    window shifted by whole years). Anchored on the (pseudo-)onset, so every start
    of one event sees the same observed outcome; anchoring on the start date made
    it depend on which start was listed first (review 2026-09-25)."""
    out = {}
    for a in anchors:
        vals = []
        for d in offsets:
            v = ob.get(a + pd.Timedelta(days=int(d)), np.nan)
            others = [ob.get(a + pd.DateOffset(years=int(y) - a.year)
                             + pd.Timedelta(days=int(d)), np.nan)
                      for y in hyears if int(y) != a.year]
            ok = np.isfinite(others)
            vals.append(v - np.mean(np.asarray(others)[ok]) if ok.sum() >= 15 else np.nan)
        out[a] = np.array(vals)
    return out


def classify(a):
    """a: (..., days) NAM-sign anomaly -> (window mean, DW)."""
    m = np.nanmean(a, axis=-1)
    return m, (m < 0) & (np.nanmean(a < 0, axis=-1) > 0.5)


def start_stats(fc_by_init, ob, init, k, hyears):
    """Window statistics for one start whose (pseudo-)onset is init + k days."""
    leads = np.arange(k + WIN[0], k + WIN[1] + 1)
    f = fc_by_init.get(init)
    if f is None:
        return None
    piv = f[f["lead_day"].isin(leads)].pivot(index="member", columns="lead_day",
                                              values="anom")
    if piv.shape != (11, len(leads)) or piv.isna().any().any():
        return None
    onset = init + pd.Timedelta(days=int(k))
    o = obs_anom(ob, [onset], range(WIN[0], WIN[1] + 1), hyears)[onset]
    if np.isnan(o).any():
        return None
    Am, DWm = classify(-piv.values)
    Ao, DWo = classify(-o)
    return {"A_members": Am, "A_ens": float(Am.mean()), "P_dw": float(DWm.mean()),
            "A_obs": float(Ao), "DW_obs": bool(DWo),
            "pit": float((np.sum(Am < Ao) + 0.5 * np.sum(Am == Ao) + 0.5) / 12.0),
            "spread": float(Am.std(ddof=1))}


def discrimination(ev):
    """ev: list of per-event dicts (A_ens, P_dw, A_obs, DW_obs) -> r, BSS."""
    a = np.array([e["A_ens"] for e in ev]); o = np.array([e["A_obs"] for e in ev])
    p = np.array([e["P_dw"] for e in ev]); d = np.array([e["DW_obs"] for e in ev], float)
    r = float(np.corrcoef(a, o)[0, 1]) if a.std() > 0 and o.std() > 0 else np.nan
    bs = np.mean((p - d) ** 2); bs_ref = np.mean((p.mean() - d) ** 2)
    bss = float(1 - bs / bs_ref) if bs_ref > 0 else np.nan
    return r, bss


def main():
    rng = np.random.default_rng(SEED)
    fc, ob = load()
    hyears = sorted(fc["hyear"].unique())
    inits = sorted(fc["init"].unique())
    by_init = {i: g for i, g in fc.groupby("init")}
    cat = load_catalogue("primary")
    lo, hi = pd.Timestamp(min(inits)), pd.Timestamp(max(inits)) + pd.Timedelta(days=42)
    events = [o for o in cat if lo + pd.Timedelta(days=17) <= o <= hi - pd.Timedelta(days=WIN[1])]
    print(f"{len(inits)} starts {pd.Timestamp(inits[0]).date()}..{pd.Timestamp(inits[-1]).date()}, "
          f"hindcast years {hyears[0]}-{hyears[-1]}; {len(events)} catalogued events in range")

    def zone_free(p):
        return bool(np.all(np.abs((cat - p).days) > ZONE_SEP))

    res = {"plan_approved": "2026-09-25", "model_year": 2022, "origin": "ecmwf",
           "window_post_onset_days": list(WIN), "lead_bins": BINS,
           "n_null": N_NULL, "seed": SEED, "zone_sep_days": ZONE_SEP,
           "null_mode": NULL_MODE, "null_window_days": NULL_WINDOW,
           "inputs": [str(FC.relative_to(ROOT)), str(OBS.relative_to(ROOT))],
           "n_starts_total": len(inits), "bins": {}}
    for bname, (k0, k1) in BINS.items():
        # SSW starts
        ev_rows, starts = [], []
        for o in events:
            ss = []
            for init in inits:
                k = (o - pd.Timestamp(init)).days
                if k0 <= k <= k1:
                    s = start_stats(by_init, ob, pd.Timestamp(init), k, hyears)
                    if s is not None:
                        ss.append((pd.Timestamp(init), k, s))
            if not ss:
                continue
            starts += [dict(event=str(o.date()), init=str(i.date()), k=k, **{x: s[x] for x in
                        ("A_ens", "P_dw", "A_obs", "DW_obs", "pit", "spread")}) for i, k, s in ss]
            ev_rows.append({"event": str(o.date()),
                            "A_ens": float(np.mean([s["A_ens"] for *_, s in ss])),
                            "P_dw": float(np.mean([s["P_dw"] for *_, s in ss])),
                            "A_obs": ss[0][2]["A_obs"], "DW_obs": ss[0][2]["DW_obs"],
                            "starts": [(str(i.date()), k) for i, k, _ in ss]})
        r_obs, bss_obs = discrimination(ev_rows)

        # control pool: zone-free starts, by (month of init, k)
        pool = {}
        for init in inits:
            init = pd.Timestamp(init)
            for k in range(k0, k1 + 1):
                if zone_free(init + pd.Timedelta(days=k)):
                    pool.setdefault((init.month, k), []).append(init)
        cache = {}

        def ctl(init, k):
            key = (init, k)
            if key not in cache:
                cache[key] = start_stats(by_init, ob, init, k, hyears)
            return cache[key]

        # all control starts once (reliability), at every k of the bin
        ctl_all = [s for (m, k), lst in pool.items() for i in lst
                   if (s := ctl(i, k)) is not None]
        # null pool per event: lists of start-statistics for one pseudo-onset each
        alt = []
        for e in ev_rows:
            o = pd.Timestamp(e["event"]); ok = []
            for y in hyears:
                d = pd.DateOffset(years=int(y) - o.year)
                if NULL_MODE == "year_shift":
                    if int(y) == o.year or not zone_free(o + d):
                        continue
                    ss = [ctl(pd.Timestamp(i) + d, k) for i, k in e["starts"]]
                    if all(x is not None for x in ss):
                        ok.append(ss)
                else:                                    # calendar_window
                    for dd in range(-NULL_WINDOW, NULL_WINDOW + 1):
                        p = o + d + pd.Timedelta(days=dd)
                        if not zone_free(p):
                            continue
                        ss = [ctl(pd.Timestamp(i), (p - pd.Timestamp(i)).days)
                              for i in inits if k0 <= (p - pd.Timestamp(i)).days <= k1]
                        ss = [x for x in ss if x is not None]
                        if ss:
                            ok.append(ss)
            alt.append(ok)
        n_alt = [len(a) for a in alt]
        null_r, null_bss, null_var, null_ndw = [], [], [], []
        if min(n_alt) > 0:
            for _ in range(N_NULL):
                ev = []
                for a in alt:
                    ss = a[rng.integers(len(a))]
                    ev.append({"A_ens": np.mean([x["A_ens"] for x in ss]),
                               "P_dw": np.mean([x["P_dw"] for x in ss]),
                               "A_obs": ss[0]["A_obs"], "DW_obs": ss[0]["DW_obs"]})
                rr, bb = discrimination(ev)
                null_r.append(rr); null_bss.append(bb)
                null_var.append(np.var([x["A_obs"] for x in ev]))
                null_ndw.append(int(sum(x["DW_obs"] for x in ev)))
        null_r, null_bss = np.array(null_r), np.array(null_bss)
        # CONDITIONAL null. r shrinks when the predictand varies less, and these
        # events' outcomes vary less across events than random date sets do; the
        # known-truth calibration showed the unconditional test is then
        # conservative (mean p 0.67 with equal skill everywhere) -- biased toward
        # "no event-specific skill". Compare only against draws whose across-event
        # variance of the observed outcome is within +-25% of the events' (r), and
        # whose number of observed DW events is within +-1 (Brier skill).
        # Band: 0.80 <= var_null / var_obs <= 1.25 (|log ratio| < log 1.25).
        var_obs = np.var([e["A_obs"] for e in ev_rows])
        ndw_obs = int(sum(e["DW_obs"] for e in ev_rows))
        mv = np.abs(np.log(np.array(null_var) / var_obs)) < np.log(1.25) if len(null_var) else np.array([], bool)
        mb = np.abs(np.array(null_ndw) - ndw_obs) <= 1 if len(null_ndw) else np.array([], bool)

        pit_s = np.array([s["pit"] for s in starts]); pit_c = np.array([s["pit"] for s in ctl_all])

        def rank_flatness(pit):
            """Monte-Carlo chi-square p for a flat 12-bin rank histogram. The PIT
            takes only 12 values, so a KS test against a continuous uniform
            rejects a perfect ensemble at large n (review 2026-09-25)."""
            r = np.clip(np.round(pit * 12 - 0.5).astype(int), 0, 11)
            c = np.bincount(r, minlength=12); n = c.sum(); e = n / 12
            chi = ((c - e) ** 2 / e).sum()
            sim = rng.multinomial(n, [1 / 12] * 12, size=20000)
            p = float(np.mean(((sim - e) ** 2 / e).sum(axis=1) >= chi))
            k = int(c[0] + c[11])
            return {"counts": c.tolist(), "chi2": round(float(chi), 2), "p_flat": round(p, 4),
                    "outer_share": round(k / n, 3),
                    "p_outer_binomial": round(float(stats.binomtest(k, n, 2 / 12).pvalue), 4)}
        outer = lambda p: float(np.mean((p < 1.5 / 12) | (p > 10.5 / 12)))
        def spread_error(ss):
            err = np.sqrt(np.mean([(s["A_ens"] - s["A_obs"]) ** 2 for s in ss]))
            spr = np.sqrt(np.mean([s["spread"] ** 2 for s in ss])) * np.sqrt(12 / 11)
            return float(spr / err)
        out = {
            "n_events": len(ev_rows), "n_ssw_starts": len(starts),
            "n_control_starts": len(ctl_all), "n_null_valid": int(len(null_r)),
            "null_years_per_event": dict(zip([e["event"] for e in ev_rows], n_alt)),
            "events": ev_rows, "ssw_starts": starts,
            "discrimination": {
                "r_ensmean_vs_obs": round(r_obs, 4),
                "r_null_mean": round(float(np.nanmean(null_r)), 4),
                "r_null_q95": round(float(np.nanpercentile(null_r, 95)), 4),
                "p_r": round(float(np.nanmean(null_r >= r_obs)), 4),
                "BSS_P_dw": round(bss_obs, 4),
                "BSS_null_mean": round(float(np.nanmean(null_bss)), 4),
                "BSS_null_q95": round(float(np.nanpercentile(null_bss, 95)), 4),
                "p_BSS": round(float(np.nanmean(null_bss >= bss_obs)), 4),
                "p_r_conditional": round(float(np.nanmean(null_r[mv] >= r_obs)), 4) if mv.sum() else None,
                "n_null_conditional_r": int(mv.sum()),
                "r_null_conditional_mean": round(float(np.nanmean(null_r[mv])), 4) if mv.sum() else None,
                "p_BSS_conditional": round(float(np.nanmean(null_bss[mb] >= bss_obs)), 4) if mb.sum() else None,
                "n_null_conditional_BSS": int(mb.sum()),
                "var_A_obs_events": round(float(var_obs), 1),
                "var_A_obs_null_median": round(float(np.median(null_var)), 1) if len(null_var) else None,
                "n_DW_obs_events": ndw_obs},
            "shift": {"A_ens_mean_ssw": round(float(np.mean([e["A_ens"] for e in ev_rows])), 1),
                      "A_ens_mean_control": round(float(np.mean([s["A_ens"] for s in ctl_all])), 1),
                      "A_obs_mean_ssw": round(float(np.mean([e["A_obs"] for e in ev_rows])), 1),
                      "P_dw_mean_ssw": round(float(np.mean([e["P_dw"] for e in ev_rows])), 3),
                      "DW_obs_rate_ssw": round(float(np.mean([e["DW_obs"] for e in ev_rows])), 3),
                      "P_dw_mean_control": round(float(np.mean([s["P_dw"] for s in ctl_all])), 3),
                      "DW_obs_rate_control": round(float(np.mean([s["DW_obs"] for s in ctl_all])), 3)},
            "reliability": {
                "pit_ssw_mean": round(float(pit_s.mean()), 3),
                "pit_control_mean": round(float(pit_c.mean()), 3),
                "outer_bin_share_ssw": round(outer(pit_s), 3),
                "outer_bin_share_control": round(outer(pit_c), 3),
                "outer_bin_share_expected": round(2 / 12, 3),
                "rank_hist_ssw": rank_flatness(pit_s),
                "rank_hist_control": rank_flatness(pit_c),
                "spread_error_ratio_ssw": round(spread_error(starts), 3),
                "spread_error_ratio_control": round(spread_error(ctl_all), 3),
                "note": "starts of one event share an observed window, and control starts at "
                        "neighbouring k share windows: tests treat starts as independent and are "
                        "therefore liberal; the event, not the start, is the independent unit"}}
        res["bins"][bname] = out
        d, rl = out["discrimination"], out["reliability"]
        print(f"\n=== {bname} ({k0}-{k1} d before onset): {out['n_events']} events, "
              f"{out['n_ssw_starts']} SSW starts, {out['n_control_starts']} control starts")
        print(f"  r(ens mean, obs) = {d['r_ensmean_vs_obs']:+.3f}  null mean {d['r_null_mean']:+.3f}, "
              f"95th {d['r_null_q95']:+.3f}, p = {d['p_r']:.3f}")
        print(f"  conditional on outcome spread: p_r = {d['p_r_conditional']} "
              f"(n = {d['n_null_conditional_r']}); on DW count: p_BSS = {d['p_BSS_conditional']} "
              f"(n = {d['n_null_conditional_BSS']}); spread events/null "
              f"{d['var_A_obs_events']:.0f}/{d['var_A_obs_null_median']:.0f} Pa^2")
        print(f"  BSS of P(DW)     = {d['BSS_P_dw']:+.3f}  null mean {d['BSS_null_mean']:+.3f}, "
              f"95th {d['BSS_null_q95']:+.3f}, p = {d['p_BSS']:.3f}")
        print(f"  shift: A_ens SSW {out['shift']['A_ens_mean_ssw']:+.0f} Pa vs control "
              f"{out['shift']['A_ens_mean_control']:+.0f}; P(DW) {out['shift']['P_dw_mean_ssw']:.2f} "
              f"vs {out['shift']['P_dw_mean_control']:.2f}; observed DW rate "
              f"{out['shift']['DW_obs_rate_ssw']:.2f} vs {out['shift']['DW_obs_rate_control']:.2f}")
        print(f"  PIT mean {rl['pit_ssw_mean']:.2f} (control {rl['pit_control_mean']:.2f}); outer-bin share "
              f"{rl['outer_bin_share_ssw']:.2f} vs {rl['outer_bin_share_control']:.2f} (expected 0.17; "
              f"binomial p SSW {rl['rank_hist_ssw']['p_outer_binomial']:.3f}); flat-histogram p "
              f"SSW {rl['rank_hist_ssw']['p_flat']:.3f}, control {rl['rank_hist_control']['p_flat']:.3f}; spread/error "
              f"{rl['spread_error_ratio_ssw']:.2f} vs {rl['spread_error_ratio_control']:.2f}")
    out = "s2s_forecast_test.json" if NAME == "s2s_forecast_test" else f"{NAME}.json"
    (RESULTS / out).write_text(json.dumps(res, indent=2, default=str),
                               encoding="utf8", newline="\n")
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main()
