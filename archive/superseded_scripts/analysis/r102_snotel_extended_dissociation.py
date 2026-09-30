#!/usr/bin/env python3
"""
r102_snotel_extended_dissociation.py
====================================
BREAKTHROUGH analysis: test the loaded-gun MECHANISM (loading maintained while
natural release triggers are suppressed) at HIGH POWER, independently of the
Swiss/ERA5/SNOWPACK chain.

Why this raises power past the n=16 ceiling: SSWs are hemispheric, so the SAME
stratospheric events that the Swiss avalanche record could only sample 16 times
(1998/99-2018/19) are present back to 1979 in the NCEP record. The continental-US
SNOTEL network (945 stations, 1980-2026) lets us observe the snowpack response to
~37 major SSWs - more than double the Alpine event sample - with a physically
independent observing system.

Design (mirrors the paper's matched-control approach, event = unit of inference):
  1. Detect major mid-winter SSWs 1979-2024 from NCEP 10 hPa 60N zonal wind using
     the Charlton-Polvani algorithm (20-day westerly separation; final-warming
     exclusion requiring >=10 westerly days to return before 30 Apr).
  2. For each SSW in the SNOTEL era, compute national-mean snowpack/forcing
     metrics in the [0,+30] d surface-impact window vs a DOY-matched non-SSW
     climatology.
  3. LOADING arm  : SWE level, SWE accumulation rate, precipitation.
     TRIGGER arm   : melt-day fraction (Tmax>0), rain-on-snow fraction.
  4. Event-level sign / t tests across ~37 events.

Output: data/results/r102_snotel_extended_dissociation.json
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r102_snotel_extended_dissociation.json"
WIN = 30           # surface-impact window [0, +WIN] days
DOY_BW = 5         # +/- day-of-year bandwidth for matched climatology


def detect_ssw(u):
    """Charlton-Polvani major mid-winter SSW central dates from 10 hPa 60N wind."""
    u = u.dropna().sort_index()
    dates, vals = u.index, u.values
    events, west_run = [], 0
    for k in range(len(u)):
        d, v = dates[k], vals[k]
        if v > 0:
            west_run += 1
            continue
        # easterly day
        in_season = d.month in (11, 12, 1, 2, 3, 4)
        if in_season and west_run >= 20:
            season_end_year = d.year + 1 if d.month >= 11 else d.year
            apr30 = pd.Timestamp(year=season_end_year, month=4, day=30, tz="UTC")
            fut = u[(u.index > d) & (u.index <= apr30)]
            maxrun = cur = 0
            for vv in fut.values:
                if vv > 0:
                    cur += 1; maxrun = max(maxrun, cur)
                else:
                    cur = 0
            if maxrun >= 10:               # vortex recovers -> mid-winter SSW (not final warming)
                events.append(d)
        west_run = 0
    return pd.DatetimeIndex(events)


# Snow-climate classification by SNOTEL state code (paper's continental hypothesis)
CONTINENTAL = {"CO", "UT", "WY", "MT", "ID", "NV", "NM", "AZ", "SD"}
MARITIME = {"WA", "OR", "CA", "AK"}


def _state(sid):
    parts = str(sid).split(":")
    return parts[1] if len(parts) > 1 else "??"


def build_snotel_daily_by_climate():
    """Aggregate 10.6M station-days to daily metrics, split by snow climate."""
    sn = pd.read_parquet(ROOT / "data/processed/cryosphere/snotel_daily.parquet",
                         columns=["station_id", "wteq_mm", "prec_mm", "tavg_c", "tmax_c"])
    if sn.index.tz is None:
        sn.index = sn.index.tz_localize("UTC")
    sn = sn.rename_axis("date").reset_index()
    sn = sn.sort_values(["station_id", "date"])
    sn["dswe"] = sn.groupby("station_id")["wteq_mm"].diff()
    sn["melt"] = (sn["tmax_c"] > 0).astype(float)
    sn["ros"] = ((sn["prec_mm"] > 5) & (sn["tmax_c"] > 1)).astype(float)
    sn["accum"] = ((sn["prec_mm"] > 2) & (sn["tavg_c"] < 0)).astype(float)
    st = sn["station_id"].map(_state)
    sn["climate"] = np.where(st.isin(CONTINENTAL), "continental",
                             np.where(st.isin(MARITIME), "maritime", "other"))

    def agg(df):
        d = df.groupby("date").agg(
            swe=("wteq_mm", "mean"), dswe=("dswe", "mean"), prec=("prec_mm", "mean"),
            melt_frac=("melt", "mean"), ros_frac=("ros", "mean"),
            accum_frac=("accum", "mean"), n_stations=("station_id", "nunique"))
        return d[d["n_stations"] >= 15]

    out = {"all": agg(sn), "continental": agg(sn[sn["climate"] == "continental"]),
           "maritime": agg(sn[sn["climate"] == "maritime"])}
    for k, v in out.items():
        ns = sn[sn["climate"] == k]["station_id"].nunique() if k != "all" else sn["station_id"].nunique()
        print(f"  {k}: {ns} stations, {len(v)} daily rows")
    return out


METRICS = {
    "swe":        ("LOADING",  "up"),
    "dswe":       ("LOADING",  "up"),
    "prec":       ("LOADING",  "up"),
    "accum_frac": ("LOADING",  "up"),
    "melt_frac":  ("TRIGGER",  "down"),
    "ros_frac":   ("TRIGGER",  "down"),
}


def doy_matched_expectation(series, onset, ssw_years):
    """Mean over the [doy, doy+WIN] calendar window across NON-SSW years."""
    doys = [(onset.dayofyear + k - 1) % 366 + 1 for k in range(WIN + 1)]
    lo, hi = min(doys), max(doys)
    idx = series.index
    if hi - lo <= WIN + 2:
        m = (idx.dayofyear >= lo) & (idx.dayofyear <= hi)
    else:  # wrapped around year end
        m = (idx.dayofyear >= lo) | (idx.dayofyear <= hi)
    seg = series[m & (~idx.year.isin(ssw_years))]
    return seg.mean()


def main():
    strat = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet")
    if strat.index.tz is None:
        strat.index = strat.index.tz_localize("UTC")
    ssw = detect_ssw(strat["uwnd_ms_10hPa"])
    print(f"Detected {len(ssw)} major mid-winter SSWs 1979-2024 (Charlton-Polvani).")
    print("  by decade:", ssw.to_series().groupby(ssw.year // 10 * 10).count().to_dict())

    by_climate = build_snotel_daily_by_climate()

    def analyze(daily, label):
        ev = ssw[(ssw >= daily.index.min()) & (ssw <= daily.index.max() - pd.Timedelta(days=WIN))]
        ssw_years = set(ev.year) | set((ev + pd.Timedelta(days=WIN)).year)
        per_event = []
        for onset in ev:
            win = daily[(daily.index >= onset) & (daily.index <= onset + pd.Timedelta(days=WIN))]
            if len(win) < WIN * 0.5:
                continue
            row = {"onset": str(onset.date())}
            for col in METRICS:
                obs = win[col].mean()
                exp = doy_matched_expectation(daily[col], onset, ssw_years)
                row[col] = float(obs - exp) if np.isfinite(obs) and np.isfinite(exp) else np.nan
            per_event.append(row)
        pe = pd.DataFrame(per_event)
        n = len(pe)
        out = {"n_events_tested": n, "metrics": {}}
        print(f"\n[{label}] dissociation across n={n} SSW events:")
        for col, (arm, direction) in METRICS.items():
            x = pe[col].dropna().values
            if len(x) < 5:
                continue
            n_pos = int((x > 0).sum())
            hits = n_pos if direction == "up" else len(x) - n_pos
            sign_p = stats.binomtest(hits, len(x), 0.5, alternative="greater").pvalue
            t, t_p = stats.ttest_1samp(x, 0)
            d_eff = float(np.mean(x) / np.std(x, ddof=1)) if np.std(x, ddof=1) > 0 else np.nan
            out["metrics"][col] = {
                "arm": arm, "predicted_direction": direction, "n": len(x),
                "mean_anomaly": float(np.mean(x)),
                "hits_in_predicted_direction": f"{hits}/{len(x)}",
                "sign_test_p": float(sign_p), "t_p_two_sided": float(t_p), "cohens_d": d_eff}
            print(f"  {col:11s} [{arm:7s} {direction:4s}]: {hits}/{len(x)} dir, "
                  f"d={d_eff:+.2f}, sign P={sign_p:.4f}, t-P={t_p:.4f}")
        trig = [c for c in METRICS if METRICS[c][0] == "TRIGGER" and c in out["metrics"]]
        load = [c for c in METRICS if METRICS[c][0] == "LOADING" and c in out["metrics"]]
        out["triggers_all_down"] = bool(all(out["metrics"][c]["mean_anomaly"] < 0 for c in trig))
        out["loading_preservation_up"] = bool(
            out["metrics"].get("accum_frac", {}).get("mean_anomaly", -1) > 0 and
            out["metrics"].get("dswe", {}).get("mean_anomaly", -1) > 0)
        return out

    results = {"n_ssw_detected": len(ssw), "ssw_dates": [str(d.date()) for d in ssw],
               "window_days": WIN, "by_climate": {}}
    for label, daily in by_climate.items():
        results["by_climate"][label] = analyze(daily, label)

    cont = results["by_climate"]["continental"]
    mar = results["by_climate"]["maritime"]
    results["geographic_test"] = {
        "hypothesis": "loaded-gun dissociation is continental-specific (paper's claim)",
        "continental_triggers_down": cont["triggers_all_down"],
        "continental_loading_up": cont["loading_preservation_up"],
        "maritime_triggers_down": mar["triggers_all_down"],
        "confirms_continental_specificity": bool(
            cont["triggers_all_down"] and cont["loading_preservation_up"]),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print("\n=== GEOGRAPHIC TEST (paper's continental-specificity claim) ===")
    print(f"  continental: triggers_down={cont['triggers_all_down']}, "
          f"loading_up={cont['loading_preservation_up']}")
    print(f"  maritime:    triggers_down={mar['triggers_all_down']}")
    print(f"Saved -> {OUT}")


if __name__ == "__main__":
    main()
