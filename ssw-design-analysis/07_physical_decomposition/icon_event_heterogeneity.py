#!/usr/bin/env python3
"""
icon_event_heterogeneity.py
===========================
HOW MUCH DO EVENT-CONDITIONED MEAN RESPONSES DIFFER? A SUMMARY-LEVEL COMPARISON
WITH THE 18 ICON SPIN-OFF ENSEMBLES OF LOEFFEL ET AL. (2026).

Plan approved 2026-10-02 (user: summary-level comparison only). Committed BEFORE
any value in the files was read: only the file structure had been inspected
(variables, dimensions, attributes).

DATA (public archive, data.ub.uni-muenchen.de/808, retrieved through the Internet
Archive copy because the archive host timed out; raw/loeffel2026/)
  18 spin-off ensembles (40 members each; Loeffel et al. 2026, WCD 7, 895-913):
  ensemble mean and s.d. of standardised polar-cap (60-90N) geopotential height
  anomalies at 52 levels, and of u(10 hPa, 60N), 6-hourly, days 0-60 from the
  spin-off start. No member data, no cross-time covariance.

WHAT THIS CAN AND CANNOT ESTIMATE
  These are initial-condition ensembles conditioned on a model SSW and on its
  tropospheric state, not SSW-versus-no-SSW interventions; their mean differences
  are event-conditioned responses, not pure stratospheric causal effects. The 18
  events were selected by the authors, so their spread is not a population
  variance. With only daily spreads, the s.d. of a window mean is bounded, not
  known: it lies between 0 and the mean daily s.d.

ESTIMATES (fixed now)
  Onset of each ensemble: first time the ENSEMBLE-MEAN u(10 hPa, 60N) is negative
  (ensembles without one are reported and excluded). Outcome: NAM sign (minus the
  standardised anomaly) at the level nearest 1000 hPa, mean over onset+8 ..
  onset+25 days (primary; the paper's window); secondary onset+26 .. onset+42
  where the record allows.
  M1 between-event s.d. of the ensemble means (raw).
  M2 between-event variance net of finite-ensemble noise, V_b = Var(means) -
     mean(s_w^2)/40, for two bounds on the window-mean s.d. s_w: s_w = mean daily
     s.d. (upper bound on noise -> lower bound on V_b) and s_w = 0 (upper bound).
  M3 V_b as a share of the within-ensemble variance (mean daily s.d. squared):
     the analogue of the forced share of response variance; compared with the
     CMIP6 assumption-dependent estimate (0.019 [-0.033, 0.072] in daily-s.d.
     units, days 8-52), stating the difference in windows and units.
  M4 across events, correlation of the outcome mean with the week-2 (onset+8 ..
     +14) mean 100 hPa standardised anomaly (Loeffel et al.'s between-event
     relation), with a bootstrap over events.
READING (from the external audit, fixed now)
  Substantial V_b (share well above the CMIP6 upper bound) with overlapping
  distributions -> retract "events differ, but by little"; keep "different odds,
  not two realised classes". V_b comparable to CMIP6 -> report agreement with its
  limits. Member-level shape and label rates cannot be assessed from these files.

Output: results/current/5_mechanism/icon_event_heterogeneity.json
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = HERE.parents[0] / "03_data_ingestion" / "raw" / "loeffel2026" / "unpacked" / "spin-offs"
RESULTS = ROOT / "results" / "current" / "5_mechanism"
NAME = "icon_event_heterogeneity"
SEED = zlib.crc32(NAME.encode()) % (2 ** 32)
N_BOOT, N_MEM = 2000, 40
CMIP6 = {"est": 0.019, "ci95": [-0.033, 0.072], "window": "days 8-52", "units": "member daily s.d."}


def window(ds, var, lev, onset, a, b):
    t = ds["time_6h"].values
    sel = (t >= onset + a) & (t <= onset + b)
    if t.max() < onset + b:
        return None
    x = ds[var].sel(plev=lev, method="nearest").values[sel]
    return x


def main():
    res = {"primary_registered": run("first_negative"), "sensitivity_crossing_only": run("crossing"),
           "note": ("Loeffel et al. (2026, WCD 7, 895-913, Sect. 2) initialise every spin-off 'at the SSW onset or the "
                    "day before', onset being the first day the daily-mean u(10 hPa, 60N) is negative. Six ensembles "
                    "whose ensemble-mean wind is already negative at t = 0 therefore start on their onset day, and the "
                    "registered rule (first negative ensemble-mean wind) anchors all 18 within 0-1.25 days: it is the "
                    "primary estimate. A code review (2026-10-02) had read t = 0 as 'no onset in record' and kept only "
                    "the 12 westerly-to-easterly crossings; that reading ignored the published design and was reverted "
                    "the same day. The 12-event subset is kept as a sensitivity.")}
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "events"} if isinstance(v, dict) else v for k, v in res.items()}, indent=1))
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "icon_event_heterogeneity.json").write_text(json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> icon_event_heterogeneity.json")


def run(rule):
    rng = np.random.default_rng(SEED)
    rows, excl = [], {}
    for f in sorted(RAW.glob("indices_*.nc")):
        ds = xr.open_dataset(f)
        u = ds["u10_60_mean"].values; t = ds["time_6h"].values
        if rule == "crossing":
            cr = np.flatnonzero((u[1:] < 0) & (u[:-1] >= 0)) + 1
            if u[0] < 0:
                excl[f.stem] = f"ensemble-mean wind already negative at the start ({u[0]:.1f} m/s): started on its onset day; excluded in this sensitivity only"
                continue
            neg = cr
        else:
            neg = np.flatnonzero(u < 0)
        if len(neg) == 0:
            excl[f.stem] = "ensemble-mean u(10 hPa, 60N) never negative"
            continue
        on = float(t[neg[0]])
        lev_s = float(ds["plev"].sel(plev=100000.0, method="nearest"))
        lev_100 = float(ds["plev"].sel(plev=10000.0, method="nearest"))
        m1 = window(ds, "GPHano_mean", lev_s, on, 8, 25)
        s1 = window(ds, "GPHano_stdv", lev_s, on, 8, 25)
        w2 = window(ds, "GPHano_mean", lev_100, on, 8, 14)
        m2 = window(ds, "GPHano_mean", lev_s, on, 26, 42)
        s2 = window(ds, "GPHano_stdv", lev_s, on, 26, 42)
        if m1 is None:
            excl[f.stem] = f"record ends before onset+25 (onset day {on:.2f})"
            continue
        rows.append({"id": f.stem, "onset_day": on, "y": -float(np.mean(m1)), "sd_daily": float(np.mean(s1)),
                     "z100_wk2": -float(np.mean(w2)),
                     "y2": None if m2 is None else -float(np.mean(m2)), "sd2": None if s2 is None else float(np.mean(s2))})
    y = np.array([r["y"] for r in rows]); sd = np.array([r["sd_daily"] for r in rows]); x = np.array([r["z100_wk2"] for r in rows])

    def vb(yv, sdv):
        v = np.var(yv, ddof=1)
        return v - np.mean(sdv ** 2) / N_MEM, v
    lo_vb, raw_v = vb(y, sd)
    within = float(np.mean(sd ** 2))
    bs = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        a, b = vb(y[i], sd[i]); r_ = np.corrcoef(y[i], x[i])[0, 1] if np.std(x[i]) > 0 else np.nan
        bs.append((a, b, r_))
    bs = np.array(bs)
    res = {"plan_approved": "2026-10-02", "seed": SEED, "n_boot": N_BOOT, "n_members_per_ensemble": N_MEM,
           "n_events": len(rows), "excluded": excl, "events": rows,
           "M1_between_event_sd_of_means": round(float(np.sqrt(raw_v)), 4),
           "M2_between_event_variance": {"lower_bound_noise_sd_equals_daily_sd": round(float(lo_vb), 4),
                                         "upper_bound_noise_zero": round(float(raw_v), 4),
                                         "ci95_lower_bound": [round(float(q), 4) for q in np.percentile(bs[:, 0], [2.5, 97.5])],
                                         "ci95_upper_bound": [round(float(q), 4) for q in np.percentile(bs[:, 1], [2.5, 97.5])]},
           "M3_share_of_within_ensemble_variance": {"within_daily_variance": round(within, 4),
                                                    "share_lower": round(float(lo_vb / within), 4),
                                                    "share_upper": round(float(raw_v / within), 4),
                                                    "cmip6_reference": CMIP6},
           "M4_corr_week2_100hPa_vs_outcome": {"r": round(float(np.corrcoef(y, x)[0, 1]), 4),
                                               "ci95": [round(float(q), 4) for q in np.nanpercentile(bs[:, 2], [2.5, 97.5])]},
           "mean_response_days8_25": round(float(y.mean()), 4)}
    y2 = np.array([r["y2"] for r in rows if r["y2"] is not None]); s2 = np.array([r["sd2"] for r in rows if r["y2"] is not None])
    if len(y2) >= 5:
        l2, v2 = vb(y2, s2)
        res["secondary_days26_42"] = {"n_events": int(len(y2)), "between_variance_lower": round(float(l2), 4),
                                      "between_variance_upper": round(float(v2), 4),
                                      "share_lower": round(float(l2 / np.mean(s2 ** 2)), 4),
                                      "share_upper": round(float(v2 / np.mean(s2 ** 2)), 4)}
    return res


if __name__ == "__main__":
    main()
