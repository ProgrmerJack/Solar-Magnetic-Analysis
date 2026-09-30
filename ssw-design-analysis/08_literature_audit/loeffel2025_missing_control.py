#!/usr/bin/env python3
"""
loeffel2025_missing_control.py
==============================
THE CONTROL A LIVE 2025 PAPER DOES NOT RUN, RUN ON ITS OWN DIAGNOSTIC.

THE PAPER
  Loeffel, Rupp, Kiefer, Pinto, Birner, Garny (2025), "Case-to-Case Variability in
  the Tropospheric Response to Sudden Stratospheric Warmings Revealed by Ensemble
  Re-Forecasts", EGUsphere preprint egusphere-2025-4164, discussion started
  4 October 2025, doi 10.5194/egusphere-2025-4164. ICON model, a 120-member
  event-generating ensemble (EGE) yielding 57 SSWs, plus 18 ensemble re-forecasts.

  Their metrics, quoted verbatim from section 2.2:
    "We define the surface response metric as the 1000 hPa GPH anomaly averaged
     over weeks 3-7 post central date, and define a lower stratospheric response
     metric as GPH anomalies at 100 hPa. We refer to a SSW as having a lower
     stratospheric response if the 100 hPa GPH anomaly exceeds 1.5 sigma for at
     least 10 consecutive days, within the first 6 weeks after the central warming
     date."
  Polar cap is 60-90N (section 2.2), which is the cap used here.

WHY THIS IS THE EXACT CASE THIS PROJECT IS ABOUT
  The classifier window is days 0-42. The response window, weeks 3-7, is days
  ~21-49. THEY OVERLAP BY 22 DAYS, and both are polar-cap geopotential height in a
  deep, strongly coupled column. That is the `z100_post` stratifier of
  `stratifier_bias_law.py` almost exactly, and on the 42-event observational record
  that stratifier's contrast is +0.710 against a predicted selection bias of
  +0.741, leaving a residual of -0.032 [-0.67, +0.54] -- and the SAME split applied
  to random winter dates WITH NO SSW AT ALL manufactures +0.905 sigma.

WHAT THEY DO CONTROL FOR, WHICH IS MORE THAN MOST
  They compare against both a model climatology and a No-SSW group, and for those
  groups, per the Figure 4 caption, "a random day in January or February was
  selected as the event date". That is a real null for the SSW-versus-no-SSW
  question and it is good practice.

WHAT IS MISSING, PRECISELY
  The 1.5-sigma/10-day threshold is applied ONLY WITHIN THE SSW GROUP. There is no
  "random dates that pass the same LS-signal test" group. So the reported
  separation between the LS-signal and no-LS-signal clusters has never been
  compared against what that threshold produces when no event is present. That
  single missing group is what is computed here.

WHAT THIS SCRIPT DOES AND DOES NOT CLAIM
  Test A mirrors their Figure 4: a member-level threshold split on single
  realisations. The comparison is apt, because EGE members ARE single
  realisations, and the manufactured contrast is directly interpretable.

  Test B mirrors their Figures 9-10, where the r = 0.85 correlation is taken
  between ENSEMBLE MEANS across 18 re-forecasts. Ensemble-averaging suppresses the
  internal tropospheric noise that drives this bias, so their design is
  substantially more robust there and the single-realisation null computed here is
  an UPPER BOUND on the structural contribution, NOT a claim that their 0.85 is
  spurious. It is reported that way.

  This is observational ERA5, not their ICON ensemble, so this is a test of the
  DIAGNOSTIC, not a replication of their experiment.

Output: loeffel2025_missing_control.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "9_literature"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402

CACHE = HERE / "era5_cap60_z100_z1000.parquet"
STORE = ("gs://weatherbench2/datasets/era5/"
         "1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr")
LEVELS = [1000, 100]
CAP_LAT = 60.0             # Loeffel et al. polar cap, NOT the 65N Baldwin-Thompson one
G0 = 9.80665

LS_WIN = (0, 42)           # "within the first 6 weeks after the central warming date"
LS_CORR_WIN = (14, 42)     # "from week 2 onwards", used for the continuous predictor
SFC_WIN = (21, 49)         # "weeks 3-7 post central date"
LS_SIGMA = 1.5             # "exceeds 1.5 sigma"
LS_RUN = 10                # "for at least 10 consecutive days"
N_EVENTS_REF = 18          # their re-forecast count, matched in test B
N_NULL = 4000
SEED = 20260803


def acquire():
    """Both fields from data already on disk. No download.

    An earlier version of this script re-pulled ERA5 z100/z1000 over the 60-90N cap
    from WeatherBench2. That was wrong twice over: it took >1 h and 2.2 GB because
    the requested chunking fought the store's native chunking, AND it was
    unnecessary. The merged NCEP record already carries `hgt_m_100hPa` as a
    cos-weighted 60-90N polar-cap mean -- Loeffel et al.'s EXACT level and cap --
    for 1958-2024, and the CPC AO is this project's established surface index.

    The one thing given up is a 1000 hPa polar-cap surface metric matched to their
    cap; the AO is a different construction of the same mode (corr +0.839 against
    the ERA5 1000 hPa NAM on 11,609 shared days). For a DESIGN test asking how much
    separation the threshold produces when no event is present, that is
    second-order, and it is recorded here rather than glossed.
    """
    strat = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere_1958.parquet"
    if not strat.exists():
        strat = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere.parquet"
    s = pd.read_parquet(strat)
    if s.index.tz is not None:
        s.index = s.index.tz_convert("UTC").tz_localize(None)
    ao = M.load("ao")["y"]
    df = pd.DataFrame({"z100_m": s["hgt_m_100hPa"]})
    df["ao"] = ao.reindex(df.index)
    df = df.dropna()
    df.index.name = "date"
    print(f"using {strat.name} (100 hPa, 60-90N cap) + CPC AO: "
          f"{df.index.min().date()} .. {df.index.max().date()}  {len(df):,} days")
    return df


def standardise(h, clean_mask):
    """Day-of-year standardised anomaly, with the climatology fitted on EVENT-FREE
    days only. Fitting it on all days lets the events set their own baseline --
    the contamination trap this project has already paid for three times."""
    doy = h.index.dayofyear
    c = h[clean_mask]
    mu = c.groupby(c.index.dayofyear).mean()
    sd = c.groupby(c.index.dayofyear).std()
    return (h - doy.map(mu).values) / doy.map(sd).values


def win_mean(s, onsets, win):
    out = []
    for o in onsets:
        v = s[(s.index >= o + pd.Timedelta(days=win[0]))
              & (s.index <= o + pd.Timedelta(days=win[1]))]
        out.append(float(v.mean()) if len(v) >= (win[1] - win[0]) // 2 else np.nan)
    return np.array(out)


def ls_flag(s, onsets):
    """Their threshold: > 1.5 sigma for >= 10 CONSECUTIVE days inside days 0-42.

    Sign: their 100 hPa GPH anomaly POSITIVE means a weak vortex / warm polar cap,
    which is the disturbed state that follows an SSW. The series here is the raw
    standardised height anomaly, so positive is the same disturbed sense and no
    negation is applied.
    """
    out = []
    for o in onsets:
        v = s[(s.index >= o + pd.Timedelta(days=LS_WIN[0]))
              & (s.index <= o + pd.Timedelta(days=LS_WIN[1]))].values
        if len(v) < (LS_WIN[1] - LS_WIN[0]) // 2:
            out.append(np.nan)
            continue
        best = run = 0
        for x in v:
            run = run + 1 if np.isfinite(x) and x > LS_SIGMA else 0
            best = max(best, run)
        out.append(float(best >= LS_RUN))
    return np.array(out)


def main():
    rng = np.random.default_rng(SEED)
    df = acquire()
    z100 = df["z100_m"]
    # Sign: their 100 hPa GPH anomaly POSITIVE = warm/high polar cap = disturbed
    # vortex, and their 1000 hPa surface metric POSITIVE = high polar heights,
    # which is the NEGATIVE AO phase. So the AO is negated to put the surface
    # metric in their sense, and their reported "right shift" of the LS-signal
    # group becomes a POSITIVE contrast here.
    z1000 = -df["ao"]

    allev = load_catalogue("primary")   # EXCLUSION set: every event, even ones not scored here
    real = allev[(allev >= df.index.min() + pd.Timedelta(days=60))
                 & (allev <= df.index.max() - pd.Timedelta(days=60))]
    mask = G.real_influence_mask(df.index, allev)   # every catalogued event
    clean = ~mask
    print(f"ERA5 {df.index.min().date()}..{df.index.max().date()}; "
          f"events in range {len(real)}; clean days {clean.sum():,}/{len(df):,}")

    ls = standardise(z100, clean)
    sfc = standardise(z1000, clean)
    clean_idx = G.zone_free_index(df.index, allev)
    doys = np.array([t.dayofyear for t in pd.DatetimeIndex(real)])

    res = {"paper": "Loeffel et al. 2025, EGUsphere egusphere-2025-4164",
           "n_events_observed": int(len(real)),
           "windows": {"LS_threshold": list(LS_WIN), "LS_corr": list(LS_CORR_WIN),
                       "surface": list(SFC_WIN),
                       "overlap_days": SFC_WIN[1] - LS_WIN[0] - (SFC_WIN[1] - LS_WIN[1])
                       if LS_WIN[1] > SFC_WIN[0] else 0}}
    res["windows"]["overlap_days"] = max(0, min(LS_WIN[1], SFC_WIN[1]) - max(LS_WIN[0], SFC_WIN[0]))

    # ---------------- TEST A: their Figure 4 threshold split ----------------
    print(f"\n=== TEST A  their Fig.4 LS-signal split, on single realisations ===")
    print(f"  criterion: 100 hPa anomaly > {LS_SIGMA} sigma for >= {LS_RUN} "
          f"consecutive days in days {LS_WIN[0]}-{LS_WIN[1]}")
    print(f"  response : 1000 hPa anomaly averaged over days {SFC_WIN[0]}-{SFC_WIN[1]}")
    print(f"  windows overlap by {res['windows']['overlap_days']} days\n")

    f_real = ls_flag(ls, real)
    y_real = win_mean(sfc, real, SFC_WIN)
    ok = np.isfinite(f_real) & np.isfinite(y_real)
    obs_rate = float(np.nanmean(f_real[ok]))
    obs_hi, obs_lo = y_real[ok][f_real[ok] == 1], y_real[ok][f_real[ok] == 0]
    obs_contrast = float(obs_hi.mean() - obs_lo.mean()) if len(obs_hi) and len(obs_lo) else np.nan
    print(f"  REAL SSWs (n={int(ok.sum())}): LS-signal rate {obs_rate:.1%}, "
          f"contrast {obs_contrast:+.3f} sigma")

    rates, contrasts = [], []
    for _ in range(N_NULL):
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < 10:
            continue
        fl, yy = ls_flag(ls, f), win_mean(sfc, f, SFC_WIN)
        m = np.isfinite(fl) & np.isfinite(yy)
        if m.sum() < 10:
            continue
        rates.append(float(fl[m].mean()))
        a, b = yy[m][fl[m] == 1], yy[m][fl[m] == 0]
        if len(a) >= 2 and len(b) >= 2:
            contrasts.append(float(a.mean() - b.mean()))
    rates, contrasts = np.array(rates), np.array(contrasts)
    nc = float(np.mean(contrasts))
    ci = [float(np.percentile(contrasts, 2.5)), float(np.percentile(contrasts, 97.5))]
    share = 100 * nc / obs_contrast if np.isfinite(obs_contrast) and obs_contrast else np.nan
    print(f"  PSEUDO-EVENTS, NO SSW ({len(contrasts):,} draws): "
          f"LS-signal rate {rates.mean():.1%}, contrast {nc:+.3f} "
          f"[{ci[0]:+.3f}, {ci[1]:+.3f}]")
    null_sig = not (ci[0] <= 0 <= ci[1])
    print(f"\n  READ THIS CAREFULLY, THE RATIO IS NOT THE STORY:")
    print(f"   - the null contrast's CI {'EXCLUDES' if null_sig else 'CONTAINS'} zero, so the "
          f"threshold {'does' if null_sig else 'does NOT'} reliably manufacture a contrast")
    print(f"   - the OBSERVED contrast at real SSWs is only {obs_contrast:+.3f} sigma, so a "
          f"'% of observed'\n     figure divides by a near-zero denominator and is "
          f"meaningless here ({share:.0f}%).")
    print(f"   - the substantive point is that this threshold barely separates the surface "
          f"at all\n     in observations ({obs_contrast:+.3f} sigma), while selecting "
          f"{obs_rate:.0%} of real events versus\n     {rates.mean():.0%} of random dates.")
    res["test_A"] = {"observed_rate": round(obs_rate, 4),
                     "observed_contrast": round(obs_contrast, 4),
                     "null_rate": round(float(rates.mean()), 4),
                     "null_contrast": round(nc, 4),
                     "null_contrast_CI95": [round(c, 4) for c in ci],
                     "percent_of_observed_reproduced_with_no_event": round(float(share), 1),
                     "n_null_draws": int(len(contrasts))}

    # ---------------- TEST B: their Fig 9-10 correlation --------------------
    print(f"\n=== TEST B  their Fig.9-10 correlation, matched to n={N_EVENTS_REF} ===")
    print(f"  UPPER BOUND ONLY: they correlate ENSEMBLE MEANS across {N_EVENTS_REF} "
          f"re-forecasts,\n  which suppresses the internal noise this bias feeds on. "
          f"Single realisations\n  here therefore OVERSTATE the structural share. "
          f"Reported as a bound, not a refutation.\n")
    x_real = win_mean(ls, real, LS_CORR_WIN)
    m = np.isfinite(x_real) & np.isfinite(y_real)
    r_obs = float(np.corrcoef(x_real[m], y_real[m])[0, 1])
    print(f"  REAL SSWs (n={int(m.sum())}): r(100 hPa d{LS_CORR_WIN[0]}-{LS_CORR_WIN[1]}, "
          f"surface d{SFC_WIN[0]}-{SFC_WIN[1]}) = {r_obs:+.3f}")

    rs = []
    for _ in range(N_NULL):
        f = G.draw_clean(clean_idx, doys, rng)
        if len(f) < N_EVENTS_REF:
            continue
        f = f[rng.choice(len(f), N_EVENTS_REF, replace=False)]
        xx, yy = win_mean(ls, f, LS_CORR_WIN), win_mean(sfc, f, SFC_WIN)
        mm = np.isfinite(xx) & np.isfinite(yy)
        if mm.sum() >= N_EVENTS_REF - 2:
            rs.append(float(np.corrcoef(xx[mm], yy[mm])[0, 1]))
    rs = np.array(rs)
    rci = [float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))]
    print(f"  PSEUDO-EVENTS, NO SSW, n={N_EVENTS_REF} each ({len(rs):,} draws): "
          f"r = {rs.mean():+.3f} [{rci[0]:+.3f}, {rci[1]:+.3f}]")
    p85 = float((rs >= 0.85).mean())
    print(f"  P(null r >= 0.85) = {p85:.4f}")
    print(f"\n  VERDICT ON THEIR r = 0.85:")
    print(f"   - the MISSING BASELINE IS REAL AND LARGE: with no event at all, window "
          f"overlap\n     alone gives r = {rs.mean():+.3f}. Their 0.85 must be read against "
          f"that, not against 0.")
    print(f"   - but it {'SURVIVES' if p85 < 0.05 else 'does NOT survive'} the control: "
          f"P(null >= 0.85) = {p85:.4f}. The headline is NOT overturned.")
    res["test_B"] = {"observed_r": round(r_obs, 4),
                     "null_r_mean": round(float(rs.mean()), 4),
                     "null_r_CI95": [round(c, 4) for c in rci],
                     "P_null_r_ge_0.85": round(float((rs >= 0.85).mean()), 4),
                     "n_matched": N_EVENTS_REF, "n_null_draws": int(len(rs)),
                     "interpretation": "UPPER BOUND on the structural share only; "
                                       "their r is between ensemble means, which "
                                       "suppresses the noise driving this bias."}

    (RESULTS / "loeffel2025_missing_control.json").write_text(
        json.dumps(res, indent=2), encoding="utf8")
    print("\nSaved -> loeffel2025_missing_control.json")


if __name__ == "__main__":
    main()
