#!/usr/bin/env python3
"""
predictability_wave_driving.py
==============================
CLOSES THE LAST GAP IN `predictability_ceiling.py`: NO WAVE-DRIVING PREDICTOR.

THE GAP
  The headline result -- out-of-sample R^2 = 0.100 (pre-onset) against ~0.63
  implied by the DW/NDW split -- was built from window means of zonal-mean u at
  three levels and two latitude bands. The CMIP6 archive here holds zonal-mean u
  and psl ONLY (21 *_zm.nc files, no 3-D fields), so there was no eddy heat flux
  and therefore no wave-driving predictor. Wave driving is the single most
  physically obvious thing that could raise 0.100, and leaving it out would be the
  first thing a referee asked about.

  Computing a true v'T' in CMIP6 would need daily 3-D va and ta for 21 members --
  hundreds of GB. That is not worth it, and it is not necessary, because the gap
  can be closed from two directions instead.

DIRECTION 1: WAVE-DRIVING PROXIES FROM ZONAL-MEAN u (CMIP6, n=1888)
  The TEM zonal momentum equation is
        du/dt  =  f v*  +  (1 / rho0 a cos(phi)) div F  +  X
  so the wave forcing div F appears directly in the zonal-mean wind TENDENCY. This
  is not a stretch: the Charlton-Polvani wind-tendency SSW definition is built on
  exactly that identity, and vortex deceleration is the field's standard
  wave-driving diagnostic when eddy fluxes are unavailable. Added here:
    du/dt        day-to-day tendency of u, windowed  -- wave-driven deceleration
    shear        u(10 hPa) - u(100 hPa) at 60N       -- where the forcing acts
    geometry     u(75N) - u(45N) at 10 hPa           -- vortex width/position,
                                                        which sets wave refraction
  None of these is recoverable from the window MEANS already used, so this is new
  information rather than a relabelling.

DIRECTION 2: THE REAL EDDY HEAT FLUX (OBSERVATIONS, n=42)
  Observations do have v'T' -- `heatflux_daily.parquet`, 100 hPa, 45-75N,
  1958-2024. So the question "does genuine wave driving add skill beyond the
  vortex state?" is answerable directly, just with far less power.

  PREDICTED BEFORE RUNNING, and recorded here so the outcome cannot be
  reinterpreted afterwards: at n=42 the permutation null for CV R^2 sits around
  0.10-0.15, so this arm most likely CANNOT resolve a true effect of 0.10. The
  expected deliverable is a BOUND, not a detection. It is run with a deliberately
  small feature set and a nested A/B comparison rather than 38 features, because
  that is the only version with any power at all.

Output: predictability_wave_driving.json
"""
import json
import sys
import warnings
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "6_predictability"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
SIM = ROOT / "ssw-design-analysis" / "06_simulation_validation"
sys.path.insert(0, str(SIM))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "05_corrected_estimators"))
from build_catalogue import load_catalogue          # noqa: E402
import gate3_clean_null as G                        # noqa: E402
import multi_index_event_study as M                 # noqa: E402
import ensemble_precursor as EP                     # noqa: E402
import stratifier_law_cmip6 as C6                   # noqa: E402
import predictability_ceiling as P                  # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
STRAT = ROOT / "data" / "processed" / "atmospheric" / "ncep_stratosphere_1958.parquet"
HF = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "heatflux_daily.parquet"
OUT_WIN = (8, 52)
SEASON = (11, 12, 1, 2, 3, 4)
SEED = 20260804


def doy_std(s):
    d = s.index.dayofyear
    return (s - s.groupby(d).transform("mean")) / s.groupby(d).transform("std")


def build_cmip6_plus():
    """Original features PLUS wave-driving proxies."""
    rows, Ys, groups = [], [], []
    for f in sorted(RAW.glob("*_zm.nc")):
        try:
            m = EP.load_member(f)
        except Exception:
            continue
        m2 = m[np.isin(m.index.month, SEASON)].dropna()
        if len(m2) < 2000:
            continue
        on = EP.detect_ssw(m["u10"].values, m.index)   # full daily series: CP07 needs contiguous days
        if len(on) < 15:
            continue
        am = m2["am"]
        msk = C6.influence_mask(am.index, on)
        cl = am[~msk].groupby(am[~msk].index.dayofyear).mean()
        Y = C6.anom(am, on, cl, OUT_WIN)

        ds = xr.open_dataset(f)
        lat, plev = ds["lat"].values, ds["plev"].values
        u = ds["u_zm"].values
        ds.close()
        idx = m.index
        if len(idx) != u.shape[0]:
            continue
        i60 = int(np.argmin(np.abs(lat - 60.0)))
        i45 = int(np.argmin(np.abs(lat - 45.0)))
        i75 = int(np.argmin(np.abs(lat - 75.0)))
        cap = lat >= 65.0
        wcap = np.cos(np.deg2rad(lat[cap]))

        feat = {}
        lev_series = {}
        for lev in P.LEVELS_PA:
            j = int(np.argmin(np.abs(plev - lev)))
            if not np.isclose(plev[j], lev, rtol=0.1):
                continue
            s60 = pd.Series(u[:, j, i60], index=idx)
            scap = pd.Series(np.average(u[:, j, cap], weights=wcap, axis=1), index=idx)
            lev_series[lev] = (s60, scap, j)
            for src, tag in ((s60, "60N"), (scap, "cap")):
                z = doy_std(src)
                for wname, win, _ in P.WINDOWS:
                    feat[f"{P.LEVEL_LAB[lev]}_{tag}_{wname}"] = P._wm(z, on, win)
                # WAVE-DRIVING PROXY: zonal-mean wind TENDENCY. div F enters the
                # TEM momentum budget through du/dt, so this is the archive's
                # best available stand-in for the eddy heat flux.
                zt = doy_std(src.diff())
                for wname, win, _ in P.WINDOWS:
                    feat[f"dudt_{P.LEVEL_LAB[lev]}_{tag}_{wname}"] = P._wm(zt, on, win)

        # vertical shear 10-100 hPa at 60N: where the wave forcing is deposited
        if 1000.0 in lev_series and 10000.0 in lev_series:
            sh = doy_std(lev_series[1000.0][0] - lev_series[10000.0][0])
            for wname, win, _ in P.WINDOWS:
                feat[f"shear10_100_{wname}"] = P._wm(sh, on, win)
        # vortex meridional geometry at 10 hPa: sets wave refraction
        if 1000.0 in lev_series:
            j = lev_series[1000.0][2]
            gm = doy_std(pd.Series(u[:, j, i75] - u[:, j, i45], index=idx))
            for wname, win, _ in P.WINDOWS:
                feat[f"geom75_45_{wname}"] = P._wm(gm, on, win)

        feat["doy_sin"] = np.sin(2 * np.pi * pd.DatetimeIndex(on).dayofyear / 365.25)
        feat["doy_cos"] = np.cos(2 * np.pi * pd.DatetimeIndex(on).dayofyear / 365.25)
        rows.append(pd.DataFrame(feat))
        Ys.append(Y)
        groups.append(np.full(len(Y), f.stem))
    if not rows:
        return None, None, None
    return pd.concat(rows, ignore_index=True), np.concatenate(Ys), np.concatenate(groups)


def build_obs():
    """Observational arm WITH the real eddy heat flux. Deliberately few features."""
    ao = M.load("ao")["y"]
    s = pd.read_parquet(STRAT)
    if s.index.tz is not None:
        s.index = s.index.tz_convert("UTC").tz_localize(None)
    hf = pd.read_parquet(HF).set_index("date")["vT_100hPa_45_75N"]

    real = load_catalogue("primary")
    lo = max(s.index.min(), hf.index.min(), ao.index.min()) + pd.Timedelta(days=60)
    hi = min(s.index.max(), hf.index.max(), ao.index.max()) - pd.Timedelta(days=60)
    real = real[(real >= lo) & (real <= hi)]

    mask = G.real_influence_mask(ao.index, real)
    clim = ao[~mask].groupby(ao[~mask].index.dayofyear).mean()
    Y = np.array([np.nanmean(
        (ao[(ao.index >= o + pd.Timedelta(days=OUT_WIN[0]))
            & (ao.index <= o + pd.Timedelta(days=OUT_WIN[1]))]
         - clim.reindex(ao[(ao.index >= o + pd.Timedelta(days=OUT_WIN[0]))
                           & (ao.index <= o + pd.Timedelta(days=OUT_WIN[1]))]
                        .index.dayofyear).values)) for o in real])

    u10 = doy_std(s["uwnd_ms_10hPa"])
    z100 = -doy_std(s["hgt_m_100hPa"])
    vt = doy_std(hf)
    dudt = doy_std(s["uwnd_ms_10hPa"].diff())

    W = [("m45_31", (-45, -31)), ("m30_16", (-30, -16)), ("m15_01", (-15, -1))]
    base, wave = {}, {}
    for name, win in W:
        base[f"u10_{name}"] = P._wm(u10, real, win)
        base[f"z100_{name}"] = P._wm(z100, real, win)
        wave[f"vT_{name}"] = P._wm(vt, real, win)
        wave[f"dudt_{name}"] = P._wm(dudt, real, win)
    winters = np.where(pd.DatetimeIndex(real).month >= 11,
                       pd.DatetimeIndex(real).year + 1, pd.DatetimeIndex(real).year)
    return pd.DataFrame(base), pd.DataFrame(wave), Y, winters


def main():
    rng = np.random.default_rng(SEED)
    res = {}

    print("=== CMIP6 WITH WAVE-DRIVING PROXIES ===")
    X, y, g = build_cmip6_plus()
    if X is not None:
        print(f"n={np.isfinite(y).sum()} events, {X.shape[1]} features "
              f"(was 38), {len(np.unique(g))} groups")
        res["cmip6"] = P.run("CMIP6+wave", X, y, g, rng, {})
        # the same pipeline WITHOUT the wave proxies, on the same events, computed
        # here: the typed-in 0.0998/0.1121/0.4245 used until 2026-09-25 came from
        # an older event set and made the "gain" a cross-sample difference
        Xp, yp, gp, _ = P.build_cmip6()
        assert np.array_equal(np.isfinite(yp), np.isfinite(y)) and \
            np.allclose(yp[np.isfinite(yp)], y[np.isfinite(y)]), "event sets differ"
        prev = {}
        for tier, lab in ((1, "P1 pre-onset"), (2, "P2 at-onset"),
                          (3, "P3 + post-onset stratosphere")):
            cols = P.tier_cols(Xp.columns, tier) + ["doy_sin", "doy_cos"]
            prev[lab] = round(P.cv_r2(Xp[cols].values, yp, gp, "ridge"), 4)
        print(f"\n  {'tier':<32s} {'was':>8s} {'now':>8s} {'gain':>8s}")
        print("  " + "-" * 58)
        for k, v in prev.items():
            now = res["cmip6"][k]["ridge"]["cv_r2"]
            print(f"  {k:<32s} {v:+8.4f} {now:+8.4f} {now - v:+8.4f}")
        res["cmip6_previous_without_wave"] = prev

    print("\n\n=== OBSERVATIONS WITH THE REAL EDDY HEAT FLUX ===")
    Xb, Xw, yo, wint = build_obs()
    print(f"n={np.isfinite(yo).sum()} events, {len(np.unique(wint))} winters")
    print("PREDICTION MADE BEFORE RUNNING: at this n the null p95 is expected")
    print("around 0.10-0.15, so this arm should yield a BOUND, not a detection.\n")
    ob = {}
    for lab, Xt in (("vortex only (u10, z100)", Xb),
                    ("vortex + WAVE (v'T', du/dt)", pd.concat([Xb, Xw], axis=1))):
        r = P.cv_r2(Xt.values, yo, wint, "ridge", n_splits=5)
        null = P.permutation_null(Xt.values, yo, wint, "ridge", rng, 400)
        p95 = float(np.percentile(null, 95))
        p = float((null >= r).mean())
        # With 42 events in 35 winters most groups tie in size, so GroupKFold's
        # fold assignment is an arbitrary tie-break, and it moved between
        # library versions (cv_r2 -0.030 -> -0.126 on identical data). Report
        # the spread over 200 relabelings of the winters, which vary only the
        # tie-breaking, so the verdict does not rest on one arbitrary split.
        frng = np.random.default_rng(zlib.crc32(f"folds|{lab}".encode()) % (2 ** 32))
        uw = np.unique(wint)
        spread = []
        for _ in range(200):
            relab = dict(zip(uw, frng.permutation(len(uw))))
            spread.append(P.cv_r2(Xt.values, yo, np.array([relab[w] for w in wint]),
                                  "ridge", n_splits=5))
        spread = np.array(spread)
        ob[lab] = {"n_features": Xt.shape[1], "cv_r2": round(r, 4),
                   "null_p95": round(p95, 4), "p_value": round(p, 4),
                   # skill = beats the null AND beats the mean (R^2 > 0)
                   "significant": bool(p < 0.05 and r > 0),
                   "fold_assignment_spread": {
                       "n": 200, "median": round(float(np.median(spread)), 4),
                       "p05_p95": [round(float(np.percentile(spread, 5)), 4),
                                   round(float(np.percentile(spread, 95)), 4)],
                       "max": round(float(spread.max()), 4),
                       "fraction_positive": round(float((spread > 0).mean()), 3)}}
        fs = ob[lab]["fold_assignment_spread"]
        print(f"  {lab:<30s} CV R^2 = {r:+.4f}   null p95 = {p95:+.4f}   "
              f"p = {p:.3f}   {'SKILL' if (p < 0.05 and r > 0) else 'no skill'}"
              f"   | 200 fold assignments: median {fs['median']:+.4f}, "
              f"max {fs['max']:+.4f}, positive {fs['fraction_positive']:.0%}")
    res["observations"] = ob
    gain = ob["vortex + WAVE (v'T', du/dt)"]["cv_r2"] - ob["vortex only (u10, z100)"]["cv_r2"]
    print(f"\n  wave-driving gain in observations: {gain:+.4f}")
    print(f"  null p95 here is {ob['vortex only (u10, z100)']['null_p95']:+.4f}, so any")
    print(f"  true effect below that is undetectable at n={np.isfinite(yo).sum()} -- as predicted.")
    res["observations_wave_gain"] = round(float(gain), 4)

    (RESULTS / "predictability_wave_driving.json").write_text(
        json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("\nSaved -> predictability_wave_driving.json")


if __name__ == "__main__":
    main()
