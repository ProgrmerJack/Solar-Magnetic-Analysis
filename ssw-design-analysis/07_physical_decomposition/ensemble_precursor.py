#!/usr/bin/env python3
"""
ensemble_precursor.py
=====================
Settles the question 43 observed events cannot: is the pre-onset surface anomaly
before a sudden stratospheric warming REAL?

THE OBSERVATIONAL DEAD END
  On the frozen 43-event catalogue the AO at -30..-16 d is -0.82 (p=0.015). It
  loses significance under every restriction -- isolated events -0.55 (p=0.117),
  satellite era -0.39 (p=0.277), both -0.28 (p=0.496) -- but EVERY ONE of those
  intervals still contains -0.82. Neither mechanism test separates (era p=0.115,
  clustering p=0.234). The intervals widen as n falls; nothing is resolved.
  The binding constraint is 47 winters of one realisation of history.

THE FIX
  CanESM5 historical, daily, from the CMIP6 Google Cloud zarr archive with no
  authentication. Each member is 165 years; ten members give ~1,650 winters
  against 47 observed, and several hundred SSWs against 43.

  This is a model, and a model is not the atmosphere. What it can do is decide
  whether an event-study design applied to a system with KNOWN, abundant
  realisations recovers a pre-onset anomaly -- and if so, how large. If hundreds
  of model SSWs show no pre-onset anomaly while the estimator is unbiased, the
  observed -0.82 is very likely sampling noise. If they show one, it is real and
  its magnitude is finally pinned down.

SSW DETECTION -- Charlton & Polvani (2007), applied identically in every member
  (detect_ssw; validated against the published NCEP-NCAR compendium dates by
  06_simulation_validation/validate_ssw_detector.py).

SURFACE INDEX
  Annular mode as the leading EOF of zonal-mean SLP over 20-90N (Nov-Mar),
  computed per member and standardised (see load_member). Positive = strong
  vortex / positive AO, the same sign convention as the CPC AO index used
  throughout this project.

ESTIMATOR
  Identical to the observational canonical study: mutually exclusive lead/lag
  bins, nearest-onset assignment, winter fixed effects, 3 annual harmonics,
  baseline = winter days >75 d from any onset, winter-block bootstrap.

Output: ensemble_precursor_<model>.json, one per model
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "7_ensemble"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[0] / "06_simulation_validation"))
import calibrate_seasonality as C                   # noqa: E402

RAW = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "raw" / "cmip6"
BINS = [(-60, -46), (-45, -31), (-30, -16), (-15, -1),
        (0, 14), (15, 29), (30, 44), (45, 60)]
LABELS = [f"{a:+d}..{b:+d}" for a, b in BINS]
BASELINE_GAP = 75
N_BOOT = 1000
SEASON = (11, 12, 1, 2, 3, 4)


CF_UNITS = "days since 1850-01-01 12:00:00.000000"
CF_CALENDAR = "noleap"


def load_member(f):
    import xarray as xr
    ds = xr.open_dataset(f)
    # Files written before the calendar fix store time as raw integer day counts
    # (0..60224) with the CF units/calendar dropped. Decode here, at the single
    # point every analysis reads a member, rather than patching files one by one
    # -- an undecoded file otherwise collapses 165 years into one "winter" and
    # yields 1 detected SSW instead of ~90.
    if not np.issubdtype(ds["time"].dtype, np.datetime64):
        import cftime
        d = cftime.num2date(np.asarray(ds["time"].values), CF_UNITS,
                            calendar=CF_CALENDAR)
        ds = ds.assign_coords(time=pd.to_datetime(
            [f"{x.year:04d}-{x.month:02d}-{x.day:02d}" for x in d]))
    lat = ds["lat"].values
    plev = ds["plev"].values
    # plev is in Pa, so 10 hPa = 1000 Pa. argmin alone would silently return
    # 50 hPa (or worse) for any model whose daily archive omits 10 hPa, and the
    # SSW detection would then be running on the wrong level without complaint.
    j10 = int(np.argmin(np.abs(plev - 1000.0)))
    if not np.isclose(plev[j10], 1000.0, rtol=0.05):
        raise ValueError(f"{f.name}: nearest level to 10 hPa is "
                         f"{plev[j10]/100:.1f} hPa -- 10 hPa not archived, "
                         f"Charlton-Polvani detection not valid here")
    i60 = int(np.argmin(np.abs(lat - 60.0)))
    if not np.isclose(lat[i60], 60.0, atol=2.0):
        raise ValueError(f"{f.name}: nearest latitude to 60N is {lat[i60]:.1f}")
    u = ds["u_zm"].values[:, j10, i60]
    psl = ds["psl_zm"].values
    t = pd.to_datetime(ds["time"].values.astype("datetime64[ns]"))
    ds.close()

    # Annular mode as the leading EOF of zonal-mean SLP over 20-90N, computed
    # per model rather than imposed as fixed latitude bands.
    #
    # The fixed-band index (35-55N minus 65-90N) works for CanESM5 but gives
    # MIROC6 a vortex-surface correlation of only +0.073 against CanESM5's +0.337,
    # while MIROC6's raw post-onset composite is a healthy -0.503. A band index
    # assumes every model puts its annular-mode nodes in the same place; an EOF
    # lets each model define its own structure. Sign is fixed so that positive =
    # low polar pressure = positive AO, matching the CPC convention used for the
    # observations throughout this project.
    nh = lat >= 20
    # clip before sqrt: some grids report |lat| marginally >90, making cos negative
    # and the weights NaN, which silently emptied MIROC-ES2L (7 usable events)
    w = np.sqrt(np.clip(np.cos(np.deg2rad(lat[nh])), 0.0, None))
    X = psl[:, nh] * w
    win = np.isin(t.month, (11, 12, 1, 2, 3))
    Xw = X[win]
    Xw = Xw - np.nanmean(Xw, axis=0)
    good = np.isfinite(Xw).all(axis=0)
    if good.sum() < 5:
        am = np.full(len(t), np.nan)
    else:
        _, _, Vt = np.linalg.svd(Xw[:, good], full_matrices=False)
        e = Vt[0]
        pol = lat[nh][good] >= 65
        if np.sign(np.mean(e[pol])) > 0:              # want negative loading at pole
            e = -e
        am = (X[:, good] - np.nanmean(X[win][:, good], axis=0)) @ e
    return pd.DataFrame({"u10": u, "am": (am - np.nanmean(am)) / np.nanstd(am)},
                        index=t)


def detect_ssw(u, idx):
    """Charlton & Polvani (2007) central dates, as stated by Butler et al.
    (2017, ESSD 9, 63): "the central date ... occurs when the daily-mean
    zonal-mean zonal winds at 10 hPa and 60N first change from westerly to
    easterly between November and March. The winds must return to westerly
    for 20 consecutive days between events ... If the winds do not return to
    westerly for at least 10 consecutive days before 30 April, the warming is
    a final warming and is not included."

    The first version flagged ANY easterly day and skipped 20 steps, so 254 of
    1,888 CMIP6 onsets were not westerly-to-easterly changes and 84 re-counted
    one easterly spell (found by the methods audit, 2026-09-25). Days must be
    consecutive for a change or a run to count: a gap of more than 2 days
    (missing data; a noleap calendar's absent 29 Feb is allowed) breaks a run,
    and NaN counts as neither westerly nor easterly.
    """
    u = np.asarray(u, dtype=float)
    idx = pd.DatetimeIndex(idx)
    n = len(u)
    step = np.r_[np.inf, np.diff(idx.values).astype("timedelta64[D]").astype(float)]
    contiguous = step <= 2
    west = u > 0
    east = u < 0
    ndjfm = np.isin(idx.month, (11, 12, 1, 2, 3))
    onsets = []
    separated = True          # no event yet, so the first needs no separation
    run = 0                   # consecutive westerly days since the last event
    for i in range(1, n):
        if step[i] > 60:
            # A season-filtered series jumps from 30 April to 1 November; the
            # summer in between always separates two winters' events.
            separated = True
        run = (run + 1 if (west[i] and contiguous[i]) else (1 if west[i] else 0))
        if not separated and run >= 20:
            separated = True
        if not (ndjfm[i] and east[i] and west[i - 1] and contiguous[i]):
            continue
        if not separated:
            continue
        # final-warming test: 10 consecutive westerly days before 30 April
        o = idx[i]
        yr = o.year + 1 if o.month >= 11 else o.year
        end = pd.Timestamp(year=yr, month=4, day=30)
        k, r, ok = i + 1, 0, False
        while k < n and idx[k] <= end:
            r = (r + 1 if (west[k] and contiguous[k]) else (1 if west[k] else 0))
            if r >= 10:
                ok = True
                break
            k += 1
        if ok:
            onsets.append(o)
            separated = False
            run = 0
    return pd.DatetimeIndex(onsets)


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def build(d, onsets):
    on = np.sort(np.array([np.datetime64(pd.Timestamp(o), "D") for o in onsets]))
    dd = d.index.values.astype("datetime64[D]")
    pos = np.searchsorted(on, dd)
    lo = np.clip(pos - 1, 0, len(on) - 1)
    hi = np.clip(pos, 0, len(on) - 1)
    dlo = (dd - on[lo]).astype(int)
    dhi = (dd - on[hi]).astype(int)
    lag = np.where(np.abs(dlo) <= np.abs(dhi), dlo, dhi)
    code = np.full(len(dd), -1)
    for i, (a, b) in enumerate(BINS):
        code[(lag >= a) & (lag <= b)] = i
    o = d.copy()
    o["bin"] = code
    o["is_baseline"] = np.abs(lag) > BASELINE_GAP
    return o[(o["bin"] >= 0) | o["is_baseline"]].copy()


def demean(A, codes):
    _, c = np.unique(codes, return_inverse=True)
    n = c.max() + 1
    cnt = np.bincount(c, minlength=n).astype(float)
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        s = np.bincount(c, weights=A[:, j], minlength=n)
        out[:, j] = A[:, j] - (s / cnt)[c]
    return out


def fit(dd):
    y = dd["am"].values.astype(float)
    S = C.harmonics(dd["doy"].values, 3)
    D = np.column_stack([(dd["bin"].values == i).astype(float)
                         for i in range(len(BINS))])
    Ad = demean(np.column_stack([y, D, S]), dd["wid"].values)
    b, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
    return b[:len(BINS)]


def main():
    model = sys.argv[1] if len(sys.argv) > 1 else "CanESM5"
    files = sorted(RAW.glob(f"{model}_*_zm.nc"))
    if not files:
        print(f"no {model} member files in {RAW}; run acquire_cmip6_ensemble.py first")
        return 1
    print(f"{model}: {len(files)} members on disk")

    frames, n_ssw = [], 0
    for f in files:
        m = load_member(f)
        full = m
        m = m[np.isin(m.index.month, SEASON)]
        on = detect_ssw(full["u10"].values, full.index)   # full daily series: CP07 needs contiguous days
        n_ssw += len(on)
        mem = f.stem.split("_")[1]
        d = m.copy()
        d["doy"] = d.index.dayofyear
        d["wid"] = [f"{mem}_{w}" for w in winter_of(d.index)]
        dd = build(d, on)
        frames.append(dd)
        print(f"  {mem}: {len(on)} SSWs in {d['wid'].nunique()} winters "
              f"({len(on)/d['wid'].nunique():.2f}/winter)", flush=True)

    all_d = pd.concat(frames)
    nw = all_d["wid"].nunique()
    print(f"\nPOOLED: {n_ssw:,} SSWs, {nw:,} winters, {len(all_d):,} days")
    print(f"  observations for comparison: 43 SSWs, 36 winters")
    print(f"  -> {n_ssw/43:.0f}x the events, {nw/36:.0f}x the winters")

    est = fit(all_d)
    rng = np.random.default_rng(20260731)
    wids = np.array(sorted(all_d["wid"].unique()))
    idxw = {w: np.flatnonzero((all_d["wid"] == w).values) for w in wids}
    boot = []
    for _ in range(N_BOOT):
        pick = rng.choice(wids, len(wids), replace=True)
        sub = all_d.iloc[np.concatenate([idxw[w] for w in pick])].copy()
        sub["wid"] = np.concatenate(
            [np.full(len(idxw[w]), j) for j, w in enumerate(pick)])
        try:
            b = fit(sub)
            if np.all(np.isfinite(b)):
                boot.append(b)
        except Exception:
            pass
    boot = np.array(boot)

    print(f"\n{'bin':10s} {'effect':>9s} {'95% CI':>22s} {'p':>9s}")
    print("-" * 54)
    prof = {}
    for i, l in enumerate(LABELS):
        c = boot[:, i]
        p = float(2 * min((c >= 0).mean(), (c <= 0).mean()))
        ci = [float(np.percentile(c, 2.5)), float(np.percentile(c, 97.5))]
        prof[l] = {"effect": round(float(est[i]), 4),
                   "CI95": [round(ci[0], 4), round(ci[1], 4)], "p": p}
        print(f"{l:10s} {est[i]:+9.3f} [{ci[0]:+8.3f},{ci[1]:+8.3f}] {p:9.4f}"
              f"{'  *' if p < 0.05 else ''}")

    pre = prof["-30..-16"]
    print(f"\n=== THE PRE-ONSET QUESTION ===")
    print(f"  observed (43 events):  -0.820  CI [-1.413, -0.202]  p=0.015")
    print(f"  ensemble ({n_ssw:,} events): {pre['effect']:+.3f}  "
          f"CI [{pre['CI95'][0]:+.3f}, {pre['CI95'][1]:+.3f}]  p={pre['p']:.4f}")
    w = pre["CI95"][1] - pre["CI95"][0]
    print(f"  interval width: {w:.3f} vs 1.211 observed"
          + (f" -> {1.211/w:.1f}x tighter" if w > 1e-6 else "  (degenerate)"))
    if pre["CI95"][0] <= -0.820 <= pre["CI95"][1]:
        print("  the observed -0.820 IS consistent with the ensemble")
    else:
        print("  the observed -0.820 is NOT consistent with the ensemble")

    # ---- is the observed -0.820 what 43 events would give if truth were the
    # ensemble value? Model-vs-observation disagreement has two readings -- the
    # observation is noisy, or the model is wrong. Subsampling separates them:
    # draw 36-winter blocks from the ensemble, refit, and see how often an
    # estimate as extreme as -0.820 appears when the truth is known.
    print(f"\n=== IS THE OBSERVED -0.820 JUST 43-EVENT SAMPLING NOISE? ===")
    sub_est = []
    rng2 = np.random.default_rng(99)
    ssw_wids = np.array(sorted(set(all_d.loc[all_d["bin"] == 4, "wid"])))
    for _ in range(600):
        pick = rng2.choice(ssw_wids, min(36, len(ssw_wids)), replace=False)
        s = all_d[all_d["wid"].isin(pick)]
        try:
            b = fit(s)
            if np.all(np.isfinite(b)):
                sub_est.append(b[2])          # the -30..-16 bin
        except Exception:
            pass
    sub_est = np.array(sub_est)
    frac = float((sub_est <= -0.820).mean())
    print(f"  36-winter subsamples of the ensemble: mean {sub_est.mean():+.3f}, "
          f"sd {sub_est.std():.3f}")
    print(f"  5th-95th pct: [{np.percentile(sub_est,5):+.3f}, "
          f"{np.percentile(sub_est,95):+.3f}]")
    print(f"  P(estimate <= -0.820 | truth = {pre['effect']:+.3f}) = {frac:.3f}")
    print("  -> the observed value is ordinary sampling noise" if frac > 0.05
          else "  -> the observed value is extreme even allowing for 43-event noise")

    res = {"n_members": len(files), "n_ssw": int(n_ssw), "n_winters": int(nw),
           "n_days": int(len(all_d)), "n_boot": N_BOOT,
           "bins": LABELS, "profile": prof,
           "observed_pre_onset": {"effect": -0.820, "CI95": [-1.413, -0.202],
                                  "p": 0.015, "n_events": 43},
           "subsample_test": {
               "n_subsamples": int(len(sub_est)), "winters_per_subsample": 36,
               "mean": round(float(sub_est.mean()), 4),
               "sd": round(float(sub_est.std()), 4),
               "pct5_95": [round(float(np.percentile(sub_est, 5)), 4),
                           round(float(np.percentile(sub_est, 95)), 4)],
               "P_le_observed": frac}}
    (RESULTS / f"ensemble_precursor_{model}.json").write_text(json.dumps(res, indent=2),
                                                  encoding="utf8")
    print(f"\nSaved -> ensemble_precursor_{model}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
