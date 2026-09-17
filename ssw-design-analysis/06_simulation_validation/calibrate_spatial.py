#!/usr/bin/env python3
"""
calibrate_spatial.py
====================
GATE 3 for SPATIAL outcomes. Gate 3 was closed only for index-type series
(AO, NAO). Gridded fields and station networks are a different problem, because
seasonality differs by unit -- elevation, latitude, continentality -- so a single
shared seasonal curve cannot be assumed adequate. Nothing spatial may be
analysed until this passes.

TEST CASE
  SNOTEL daily mean air temperature, aggregated to (state, day). Temperature is
  deliberately chosen: it has the strongest and most unit-dependent seasonal
  cycle of anything in the archive, so it is the hardest case. If a seasonal
  specification calibrates here it will calibrate on smoother fields.

  Units are US states (continental + maritime), 1980-2026. State x winter are
  the strata; bootstrap resamples whole winters, taking every state of that
  winter together, because states share weather within a winter.

SPECIFICATIONS COMPARED
  shared_harm3      one 3-harmonic curve for all units  <- the plan calls this
                    inadmissible for spatial data; tested to show whether that
                    is true here rather than assumed
  unit_harm3        3 harmonics interacted with unit    <- unit-specific shape
  unit_harm2        2 harmonics interacted with unit    <- fewer parameters,
                    since Gate 3 found LESS flexibility calibrated better

NULL CONSTRUCTION (as fixed in Gate 3)
  Real-event influence is removed from the DATA first ([-60, +75] d of any real
  onset), then pseudo-onsets are drawn among the surviving days preserving the
  observed day-of-year distribution. Both bins and baseline are therefore clean.

BUDGET
  Null means at high resolution (cheap: one fit per draw). Coverage and FPR at
  lower resolution because each draw needs 1,000 bootstrap replicates on ~10^5
  rows; the resulting Clopper-Pearson intervals are correspondingly wide and are
  reported as such rather than being quoted as precise.

Output: spatial_calibration.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import calibrate_seasonality as C
import gate3_clean_null as G

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "02_event_catalogues"))
from build_catalogue import load_catalogue      # noqa: E402

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "3_calibration"
RESULTS.mkdir(parents=True, exist_ok=True)
OUT = RESULTS / "spatial_calibration.json"

SEASON = (11, 12, 1, 2, 3, 4)
BINS = C.BINS
TEST_BIN = C.TEST_BIN
BASELINE_GAP = C.BASELINE_GAP
N_MEAN = 1500
N_COVER = 60           # low: each draw costs 1,000 fits on ~10^5 rows
COVER_BEST_ONLY = True  # coverage only for the spec with least bias
N_BOOT = 1000          # Gate 3 standing rule
STATES = {"CO", "UT", "WY", "MT", "ID", "NV", "WA", "OR", "CA", "AK"}


def clopper_pearson(k, n, a=0.05):
    lo = stats.beta.ppf(a / 2, k, n - k + 1) if k else 0.0
    hi = stats.beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return round(float(lo), 3), round(float(hi), 3)


def load_spatial():
    sn = pd.read_parquet(
        C.ROOT / "data/processed/cryosphere/snotel_daily.parquet",
        columns=["station_id", "tavg_c"])
    if not isinstance(sn.index, pd.DatetimeIndex):
        sn = sn.set_index("date")
    if sn.index.tz is not None:
        sn.index = sn.index.tz_convert("UTC").tz_localize(None)
    sn = sn[np.isin(sn.index.month, SEASON)]
    st = sn["station_id"].map(lambda x: str(x).split(":")[1] if ":" in str(x) else "??")
    sn = sn[st.isin(STATES)].assign(unit=st[st.isin(STATES)].values)
    sn = sn.rename_axis("date").reset_index()
    agg = (sn.groupby(["unit", "date"], as_index=False)
             .agg(y=("tavg_c", "mean"), n=("station_id", "nunique")))
    agg = agg[agg["n"] >= 5].dropna(subset=["y"])
    agg["winter"] = np.where(agg["date"].dt.month >= 11,
                             agg["date"].dt.year + 1, agg["date"].dt.year)
    agg["doy"] = agg["date"].dt.dayofyear
    return agg


def seasonal_design(d, spec):
    doy = d["doy"].values
    if spec == "shared_harm3":
        return C.harmonics(doy, 3)
    k = 3 if spec == "unit_harm3" else 2
    H = C.harmonics(doy, k)
    units = pd.Categorical(d["unit"])
    cols = []
    for i in range(len(units.categories)):
        sel = (units.codes == i).astype(float)[:, None]
        cols.append(H * sel)
    return np.hstack(cols)


def assign(d, onsets):
    on = np.sort(np.array([np.datetime64(pd.Timestamp(o), "D") for o in onsets]))
    dd = pd.DatetimeIndex(d["date"]).values.astype("datetime64[D]")
    pos = np.searchsorted(on, dd)
    lo = np.clip(pos - 1, 0, len(on) - 1)
    hi = np.clip(pos, 0, len(on) - 1)
    dlo = (dd - on[lo]).astype(int)
    dhi = (dd - on[hi]).astype(int)
    lag = np.where(np.abs(dlo) <= np.abs(dhi), dlo, dhi)
    code = np.full(len(dd), -1)
    for i, (a, b) in enumerate(BINS):
        code[(lag >= a) & (lag <= b)] = i
    return code, lag


def demean(A, codes):
    _, c = np.unique(codes, return_inverse=True)
    n = c.max() + 1
    cnt = np.bincount(c, minlength=n).astype(float)
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        s = np.bincount(c, weights=A[:, j], minlength=n)
        out[:, j] = A[:, j] - (s / cnt)[c]
    return out


def fit(d, onsets, spec):
    code, lag = assign(d, onsets)
    keep = (code >= 0) | (np.abs(lag) > BASELINE_GAP)
    dd = d[keep]
    code = code[keep]
    if len(dd) < 500:
        return None
    D = np.column_stack([(code == i).astype(float) for i in range(len(BINS))])
    S = seasonal_design(dd, spec)
    A = np.column_stack([dd["y"].values.astype(float), D, S])
    strata = pd.Categorical(dd["unit"].astype(str) + "|"
                            + dd["winter"].astype(str)).codes
    Ad = demean(A, strata)
    try:
        beta, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
        return beta[:len(BINS)]
    except Exception:
        return None


def boot_ci(d, onsets, spec, n_boot=N_BOOT, seed=0):
    rng = np.random.default_rng(seed)
    winters = np.array(sorted(d["winter"].unique()))
    idx = {w: np.flatnonzero((d["winter"] == w).values) for w in winters}
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(winters, len(winters), replace=True)
        rows = np.concatenate([idx[w] for w in pick])
        sub = d.iloc[rows].copy()
        sub["winter"] = np.concatenate(
            [np.full(len(idx[w]), j) for j, w in enumerate(pick)])
        b = fit(sub, onsets, spec)
        if b is not None and np.isfinite(b[TEST_BIN]):
            vals.append(b[TEST_BIN])
    if len(vals) < 50:
        return None
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def draw_clean(clean_dates, doys, rng):
    by = {}
    for t in clean_dates:
        by.setdefault(t.dayofyear, []).append(t)
    out = []
    for doy in rng.permutation(doys):
        c = by.get(int(doy))
        if c:
            out.append(c[rng.integers(len(c))])
    return pd.DatetimeIndex(out)


def main():
    d = load_spatial()
    onsets = load_catalogue("primary")
    dmin, dmax = d["date"].min(), d["date"].max()
    onsets = onsets[(onsets >= dmin) & (onsets <= dmax)]

    infl = G.real_influence_mask(pd.DatetimeIndex(d["date"]), onsets)
    clean = d[~infl].copy()
    uniq = pd.DatetimeIndex(sorted(set(clean["date"])))
    doys = pd.DatetimeIndex(onsets).dayofyear.values
    print(f"SNOTEL state-days {len(d):,} across {d['unit'].nunique()} units, "
          f"{d['winter'].nunique()} winters, {dmin.date()}..{dmax.date()}")
    print(f"events in range {len(onsets)}; removed {infl.sum():,} influenced "
          f"({100*infl.mean():.0f}%); clean {len(clean):,} rows, "
          f"{len(uniq):,} distinct clean dates")

    res = {"n_mean_draws": N_MEAN, "n_cover_draws": N_COVER, "n_boot": N_BOOT,
           "n_rows": int(len(d)), "n_units": int(d["unit"].nunique()),
           "n_winters": int(d["winter"].nunique()), "n_events": int(len(onsets)),
           "specs": {}}

    # ---- phase 1: null means for every spec (cheap, one fit per draw) -------
    means = {}
    for spec in ("shared_harm3", "unit_harm2", "unit_harm3"):
        rng = np.random.default_rng(11)
        ests = []
        for _ in range(N_MEAN):
            fk = draw_clean(uniq, doys, rng)
            if len(fk) < 5:
                continue
            b = fit(clean, fk, spec)
            if b is not None and np.isfinite(b[TEST_BIN]):
                ests.append(b[TEST_BIN])
        ests = np.array(ests)

        means[spec] = ests
        e0 = {"null_mean_degC": round(float(ests.mean()), 4),
              "null_sd": round(float(ests.std()), 4),
              "n_draws": int(len(ests)),
              "null_mean_in_sd_units": round(float(ests.mean() / clean["y"].std()), 4)}
        res["specs"][spec] = e0
        print(f"  [means] {spec:14s} null={e0['null_mean_degC']:+.4f} degC "
              f"({e0['null_mean_in_sd_units']:+.3f} sd)  sd={e0['null_sd']:.3f}  "
              f"n={e0['n_draws']}")

    # ---- phase 2: coverage only for the least-biased spec -------------------
    best = min(means, key=lambda s: abs(means[s].mean()))
    print(f"\n  least-biased spec: {best} -> computing coverage "
          f"({N_COVER} draws x {N_BOOT} reps)")
    for spec in ([best] if COVER_BEST_ONLY else list(means)):
        ests = means[spec]
        rng = np.random.default_rng(23)
        cov = rej = 0
        for _ in range(N_COVER):
            fk = draw_clean(uniq, doys, rng)
            ci = boot_ci(clean, fk, spec, seed=int(rng.integers(1e6)))
            if ci:
                if ci[0] <= 0 <= ci[1]:
                    cov += 1
                else:
                    rej += 1
        n = cov + rej
        res["specs"][spec].update({
            "mean_CI95": [round(float(ests.mean() - 1.96 * ests.std() / np.sqrt(len(ests))), 4),
                          round(float(ests.mean() + 1.96 * ests.std() / np.sqrt(len(ests))), 4)],
            "coverage": round(cov / n, 3) if n else None,
            "coverage_CI95": clopper_pearson(cov, n) if n else None,
            "fpr": round(rej / n, 3) if n else None,
            "n_cover": int(n)})
        e = res["specs"][spec]
        print(f"  [coverage] {spec:14s} cov={e['coverage']} {e['coverage_CI95']}  "
              f"FPR={e['fpr']}  (n={e['n_cover']}; interval wide by design)")
    res["coverage_spec"] = best

    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
