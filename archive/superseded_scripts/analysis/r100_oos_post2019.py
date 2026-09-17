#!/usr/bin/env python3
"""
r100_oos_post2019.py
====================
Referee fix M3 - frozen-pipeline out-of-sample (OOS) test on post-2019 SSW events.

Generalises scripts/analysis/35_prospective_2021_test.py from the single Jan-2021
event to ALL post-study SSW events, applying the SAME matched-control design used
in the primary analysis, with the reference climatology frozen to data the primary
analysis could/would use.

IMPORTANT, HONEST LIMITATIONS (documented, not hidden):
  1. The post-study era (2019/20 - 2024/25) contains only ~2 MAJOR mid-winter SSWs
     (5 Jan 2021, 16 Feb 2023). This is a hard physical limit: n_OOS ~ 2, so NO
     post-2019 OOS test can be statistically powered. It is supportive context.
  2. The PRIMARY endpoint (natural dry-slab COUNTS) ends 2019-05 in all local data.
     A count-channel OOS therefore requires a new SLF data request; the function
     `oos_counts()` is written and will run the instant that file is provided.
  3. Available OOS channels are therefore: (a) accidents (SLF, to 2025) - the WEAK
     consequence channel; (b) danger ratings (EnviDat, to 2024) - the bulletin arm.

SSW identification is frozen: events are verified from the NCEP 10 hPa 60N
zonal-mean zonal-wind reversal (u < 0), the same criterion family as the catalog.

Output: data/results/r100_oos_post2019.json
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r100_oos_post2019.json"

WINDOW = 15                      # +/- days, matches primary design
POST2019_START = pd.Timestamp("2019-06-01", tz="UTC")
# Frozen candidate post-study major mid-winter SSW onsets (documented major SSWs;
# verified below via 10 hPa wind reversal). March final-warming reversals (e.g.
# 2022-03-19) are excluded by design, matching the DJF mid-winter catalog scope.
CANDIDATE_SSW = ["2021-01-05", "2023-02-16", "2024-01-16"]
# Post-2019 winters and whether each contains a major SSW (for DOY-matched control).
POST_WINTERS = {
    "2019-20": ("2019-12-01", "2020-04-30", False),
    "2020-21": ("2020-12-01", "2021-04-30", True),   # Jan 2021 SSW
    "2021-22": ("2021-12-01", "2022-04-30", False),
    "2022-23": ("2022-12-01", "2023-04-30", True),   # Feb 2023 SSW
    "2023-24": ("2023-12-01", "2024-04-30", True),    # Jan 2024 SSW
    "2024-25": ("2024-12-01", "2025-04-30", False),
}


def to_utc(ts):
    ts = pd.Timestamp(ts)
    return ts if ts.tzinfo is not None else ts.tz_localize("UTC")


def verify_ssw(ncep_strat):
    """Confirm each candidate event is a genuine SSW via 10 hPa wind reversal,
    and scan all post-2019 winters for any reversal we might have missed."""
    u10 = ncep_strat["uwnd_ms_10hPa"].dropna()
    u10 = u10[u10.index >= POST2019_START]
    verified = []
    for d in CANDIDATE_SSW:
        onset = to_utc(d)
        w = u10[(u10.index >= onset - pd.Timedelta(days=20)) &
                (u10.index <= onset + pd.Timedelta(days=20))]
        reversed_ = bool((w < 0).any())
        verified.append({
            "onset": d, "u10_min_ms": float(w.min()) if len(w) else None,
            "reversal_confirmed": reversed_,
            "n_reversal_days_pm20": int((w < 0).sum()),
        })
    # independent scan: any sustained (>=1 day) DJFM reversal post-2019
    djfm = u10[u10.index.month.isin([12, 1, 2, 3])]
    neg = djfm[djfm < 0]
    scan_onsets = []
    if len(neg):
        prev = None
        for d in neg.index:
            if prev is None or (d - prev).days > 20:
                scan_onsets.append(str(d.date()))
            prev = d
    return verified, scan_onsets


def doy_window_mask(index, onset, window=WINDOW):
    doy = onset.dayofyear
    lo, hi = doy - window, doy + window
    idoy = index.dayofyear
    if lo >= 1 and hi <= 366:
        return (idoy >= lo) & (idoy <= hi)
    lo = (lo - 1) % 366 + 1
    hi = (hi - 1) % 366 + 1
    return (idoy >= lo) | (idoy <= hi)


def daily_rate_in_window(daily, onset, win_start, win_end, window=WINDOW):
    """Mean daily value within +/- window of onset, restricted to [win_start,win_end]."""
    onset = to_utc(onset)
    a = daily[(daily.index >= onset - pd.Timedelta(days=window)) &
              (daily.index <= onset + pd.Timedelta(days=window)) &
              (daily.index >= to_utc(win_start)) & (daily.index <= to_utc(win_end))]
    return a.mean(), len(a)


def matched_expectation(daily, onset, control_winters):
    """DOY-matched expectation from post-2019 NON-SSW winters (controls the secular
    trend level while remaining out-of-sample w.r.t. the pre-2019 training)."""
    onset = to_utc(onset)
    vals = []
    for (cs, ce) in control_winters:
        seg = daily[(daily.index >= to_utc(cs)) & (daily.index <= to_utc(ce))]
        if len(seg) == 0:
            continue
        m = doy_window_mask(seg.index, onset)
        sv = seg[m]
        if len(sv):
            vals.append(sv.mean())
    return (np.mean(vals) if vals else np.nan), len(vals)


def oos_channel(daily, label, predict_direction):
    """Generic OOS matched test for a daily Swiss-mean series (accidents or danger)."""
    control_winters = [(s, e) for k, (s, e, is_ssw) in POST_WINTERS.items() if not is_ssw]
    rows = []
    for d in CANDIDATE_SSW:
        onset = to_utc(d)
        # find which winter this event sits in
        host = next((k for k, (s, e, _) in POST_WINTERS.items()
                     if to_utc(s) <= onset <= to_utc(e)), None)
        ws, we, _ = POST_WINTERS[host]
        obs, n_obs = daily_rate_in_window(daily, onset, ws, we)
        exp, n_ctrl = matched_expectation(daily, onset, control_winters)
        rr = float(obs / exp) if (exp and exp > 0 and np.isfinite(obs)) else None
        rows.append({
            "onset": d, "host_winter": host,
            "observed_mean": float(obs) if np.isfinite(obs) else None,
            "expected_doy_matched": float(exp) if np.isfinite(exp) else None,
            "rate_ratio": rr, "n_obs_days": int(n_obs), "n_control_winters": n_ctrl,
            "direction": (None if rr is None else
                          ("DECREASE" if rr < 1 else "INCREASE" if rr > 1 else "FLAT")),
        })
    rrs = [r["rate_ratio"] for r in rows if r["rate_ratio"] is not None]
    n = len(rrs)
    if predict_direction == "decrease":
        hits = sum(rr < 1 for rr in rrs)
    else:
        hits = sum(rr > 1 for rr in rrs)
    gm = float(np.exp(np.mean(np.log(rrs)))) if rrs else None
    sign_p = stats.binomtest(hits, n, 0.5).pvalue if n else None
    return {
        "channel": label, "predicted_direction": predict_direction,
        "per_event": rows, "n_events": n, "gm_rate_ratio": gm,
        "hits_in_predicted_direction": f"{hits}/{n}",
        "sign_test_p": float(sign_p) if sign_p is not None else None,
        "note": ("UNDERPOWERED: only ~2 post-study major SSWs exist; "
                 "this is supportive context, not a powered confirmation."),
    }


def load_accidents_daily():
    acc = pd.read_parquet(ROOT / "data/processed/cryosphere/slf_accidents.parquet")
    if acc.index.tz is None:
        acc.index = acc.index.tz_localize("UTC")
    daily = acc.assign(_n=1).groupby(acc.index.normalize())["_n"].sum()
    # reindex to continuous daily calendar (post-2019) so zero-accident days count
    full = pd.date_range(POST2019_START, acc.index.max(), freq="D", tz="UTC")
    return daily.reindex(full, fill_value=0).astype(float)


def load_danger_daily():
    """Swiss-mean daily danger level from EnviDat post-2019 files (best-effort)."""
    frames = []
    f1 = ROOT / "data/cryosphere/envidat/swiss_danger_2020_2023.csv"
    f2 = ROOT / "data/cryosphere/envidat/swiss_danger_2023_2024.csv"
    if f1.exists():
        d = pd.read_csv(f1)
        d = d.rename(columns={"validDate": "date", "dangerLevel": "danger"})
        frames.append(d[["date", "danger"]])
    if f2.exists():
        d = pd.read_csv(f2)
        d = d.rename(columns={"dangerlevel_tidy": "danger"})
        dc = "date" if "date" in d.columns else d.columns[0]
        frames.append(d[[dc, "danger"]].rename(columns={dc: "date"}))
    if not frames:
        return None
    alld = pd.concat(frames, ignore_index=True)
    alld["date"] = pd.to_datetime(alld["date"], errors="coerce", utc=True)
    alld["danger"] = pd.to_numeric(alld["danger"], errors="coerce")
    alld = alld.dropna(subset=["date", "danger"])
    return alld.groupby(alld["date"].dt.normalize())["danger"].mean()


def oos_counts():
    """COUNT-CHANNEL OOS (primary endpoint). Runs IFF post-2019 natural dry-slab
    counts are provided at the path below (NOT in local data; SLF request needed)."""
    p = ROOT / "data/processed/cryosphere/slf_activity_post2019.parquet"
    if not p.exists():
        return {"status": "DATA_UNAVAILABLE",
                "needed_file": str(p),
                "request": ("Daily SLF natural-trigger dry-slab activity index "
                            "(aai_dry_natural / dry_natural_size_1234) for winters "
                            "2019/20-2024/25 from WSL/SLF."),
                "note": "Function is implemented; will execute on file arrival."}
    slf = pd.read_parquet(p)
    if slf.index.tz is None:
        slf.index = slf.index.tz_localize("UTC")
    col = "aai_dry_natural" if "aai_dry_natural" in slf.columns else slf.columns[0]
    return oos_channel(slf[col].astype(float), f"natural_dry_slab_counts ({col})", "decrease")


def main():
    ncep_strat = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet")
    if ncep_strat.index.tz is None:
        ncep_strat.index = ncep_strat.index.tz_localize("UTC")

    verified, scan = verify_ssw(ncep_strat)
    print("SSW verification (NCEP 10 hPa wind reversal):")
    for v in verified:
        print(f"  {v['onset']}: u10_min={v['u10_min_ms']:.1f} m/s, "
              f"reversal={v['reversal_confirmed']} ({v['n_reversal_days_pm20']} neg days)")
    print(f"Independent DJFM reversal scan onsets (post-2019): {scan}")

    acc_daily = load_accidents_daily()
    acc_res = oos_channel(acc_daily, "fatal_accidents", "increase")  # loaded-gun: accidents up

    danger_daily = load_danger_daily()
    if danger_daily is not None and len(danger_daily):
        dang_res = oos_channel(danger_daily, "bulletin_danger_level", "increase")  # danger up
    else:
        dang_res = {"channel": "bulletin_danger_level", "status": "unavailable"}

    counts_res = oos_counts()

    print(f"\n=== OOS RESULTS (post-2019, n={len(CANDIDATE_SSW)} verified SSWs) ===")
    for res in (acc_res, dang_res):
        if "per_event" in res:
            print(f"\n{res['channel']} (predict {res['predicted_direction']}):")
            for r in res["per_event"]:
                print(f"  {r['onset']} [{r['host_winter']}]: obs={r['observed_mean']}, "
                      f"exp={r['expected_doy_matched']}, RR={r['rate_ratio']} {r['direction']}")
            print(f"  -> {res['hits_in_predicted_direction']} in predicted direction; "
                  f"gmRR={res['gm_rate_ratio']}; sign P={res['sign_test_p']}")
    print(f"\ncount channel (PRIMARY endpoint): {counts_res.get('status', 'ran')}")

    out = {
        "description": "M3 frozen-pipeline OOS test on post-2019 SSW events",
        "design": ("matched-control DOY design; reference = post-2019 NON-SSW winters "
                   "(controls secular trend, out-of-sample w.r.t. pre-2019 training)"),
        "hard_limitation": ("only ~2 major mid-winter SSWs exist post-study "
                            "(2021-01-05, 2023-02-16); n_OOS~2 -> underpowered by construction"),
        "ssw_verification": verified,
        "post2019_reversal_scan": scan,
        "accident_channel": acc_res,
        "danger_channel": dang_res,
        "count_channel_primary": counts_res,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
