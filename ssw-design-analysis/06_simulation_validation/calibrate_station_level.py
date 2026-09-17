#!/usr/bin/env python3
"""
calibrate_station_level.py
==========================
GATE 3 for STATION-LEVEL outcomes -- the case the aggregated test could not
speak to.

WHY THE AGGREGATED RESULT DOES NOT TRANSFER
  Calibrating on state-aggregated SNOTEL temperature found shared and
  unit-specific seasonality indistinguishable (+0.157 vs +0.157 degC). But
  averaging ~90 stations into a state removes precisely the heterogeneity the
  gate exists to test. The acquired metadata shows how much is removed:
  **elevation spans 3 to 3,560 m and latitude 33.7 to 71.3 degrees**. A single
  seasonal curve shared between a sea-level Alaskan site and a 3,500 m Colorado
  site is a far stronger assumption than one shared between two state means.

WHAT IS TESTED
  shared_harm3   one seasonal curve for every station
  band_harm3     one curve per elevation band (quartiles) -- a middle option
  station_clim   each station's OWN leave-one-winter-out day-of-year climatology
                 subtracted before fitting; the winter under test never
                 contributes to its own climatology

  A heterogeneity diagnostic is reported first: if per-station seasonal cycles
  were in fact similar, the comparison would be uninformative and that has to be
  visible rather than assumed.

SAMPLING
  Bootstrapping 1,000 replicates over ~5 million station-days is not tractable.
  A stratified subsample of stations is drawn across elevation deciles, which
  preserves the heterogeneity under test while making the fits affordable. The
  subsample is fixed by seed and reported.

NULL CONSTRUCTION  as fixed in Gate 3: real-event influence removed from the
data first, pseudo-onsets drawn among surviving days preserving the observed
day-of-year distribution, winter-block bootstrap resampling whole winters with
all stations of that winter together.

Output: station_calibration.json
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
ROOT = HERE.parents[1]
OUT = RESULTS / "station_calibration.json"
META = ROOT / "ssw-design-analysis" / "03_data_ingestion" / "snotel_stations.csv"

SEASON = (11, 12, 1, 2, 3, 4)
BINS = C.BINS
TEST_BIN = C.TEST_BIN
BASELINE_GAP = C.BASELINE_GAP
N_STATIONS = 120          # stratified across elevation deciles
N_MEAN = 400
N_COVER = 40              # each draw costs 1,000 fits
N_BOOT = 1000             # Gate 3 floor
SEED_SAMPLE = 20260729


def clopper_pearson(k, n, a=0.05):
    lo = stats.beta.ppf(a / 2, k, n - k + 1) if k else 0.0
    hi = stats.beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return round(float(lo), 3), round(float(hi), 3)


def load_stations():
    meta = pd.read_csv(META)
    meta = meta.dropna(subset=["elevation_m", "latitude"])
    rng = np.random.default_rng(SEED_SAMPLE)
    meta["band"] = pd.qcut(meta["elevation_m"], 10, labels=False, duplicates="drop")
    per = max(1, N_STATIONS // meta["band"].nunique())
    keep = (meta.groupby("band", group_keys=False)
                .apply(lambda g: g.sample(min(per, len(g)), random_state=SEED_SAMPLE)))
    sel = set(keep["station_id"])

    sn = pd.read_parquet(ROOT / "data/processed/cryosphere/snotel_daily.parquet",
                         columns=["station_id", "tavg_c"])
    if not isinstance(sn.index, pd.DatetimeIndex):
        sn = sn.set_index("date")
    if sn.index.tz is not None:
        sn.index = sn.index.tz_convert("UTC").tz_localize(None)
    sn = sn[np.isin(sn.index.month, SEASON)]
    sn = sn[sn["station_id"].isin(sel)].dropna(subset=["tavg_c"])
    sn = sn.rename_axis("date").reset_index().rename(columns={"tavg_c": "y"})
    sn["winter"] = np.where(sn["date"].dt.month >= 11,
                            sn["date"].dt.year + 1, sn["date"].dt.year)
    sn["doy"] = sn["date"].dt.dayofyear
    sn = sn.merge(keep[["station_id", "elevation_m", "latitude", "band"]],
                  on="station_id", how="left")
    return sn, keep


def heterogeneity(d):
    """Do per-station seasonal cycles actually differ? Report, do not assume."""
    clim = (d.groupby(["station_id", "doy"])["y"].mean()
              .reset_index())
    amp = (clim.groupby("station_id")["y"]
               .agg(lambda s: s.max() - s.min()))
    mean_lvl = clim.groupby("station_id")["y"].mean()
    el = d.groupby("station_id")["elevation_m"].first()
    j = pd.concat([amp.rename("amp"), mean_lvl.rename("lvl"), el], axis=1).dropna()
    r_amp = stats.pearsonr(j["elevation_m"], j["amp"])
    r_lvl = stats.pearsonr(j["elevation_m"], j["lvl"])
    return {
        "n_stations": int(len(j)),
        "seasonal_amplitude_degC": {
            "min": round(float(j["amp"].min()), 2),
            "median": round(float(j["amp"].median()), 2),
            "max": round(float(j["amp"].max()), 2)},
        "corr_elevation_vs_amplitude": [round(float(r_amp[0]), 3), float(r_amp[1])],
        "corr_elevation_vs_mean_level": [round(float(r_lvl[0]), 3), float(r_lvl[1])],
    }


def station_climatology(d):
    """Leave-one-winter-out per-station day-of-year climatology."""
    tot = d.groupby(["station_id", "doy"])["y"].agg(["sum", "count"])
    out = np.empty(len(d))
    for w, idx in d.groupby("winter").groups.items():
        sub = d.loc[idx]
        own = sub.groupby(["station_id", "doy"])["y"].agg(["sum", "count"])
        rest = tot.subtract(own, fill_value=0)
        clim = (rest["sum"] / rest["count"].replace(0, np.nan))
        key = pd.MultiIndex.from_arrays([sub["station_id"], sub["doy"]])
        vals = clim.reindex(key).values
        # a station-day with no other winter falls back to the station mean
        fallback = sub.groupby("station_id")["y"].transform("mean").values
        out[d.index.get_indexer(idx)] = np.where(np.isnan(vals), fallback, vals)
    return out


def seasonal_design(d, spec):
    if spec == "shared_harm3":
        return C.harmonics(d["doy"].values, 3)
    if spec == "band_harm3":
        H = C.harmonics(d["doy"].values, 3)
        b = pd.Categorical(d["band"] // 3)          # quartile-ish groups
        cols = [H * (b.codes == i).astype(float)[:, None]
                for i in range(len(b.categories))]
        return np.hstack(cols)
    return np.zeros((len(d), 0))                     # station_clim: y residualised


def assign(dates, onsets):
    on = np.sort(np.array([np.datetime64(pd.Timestamp(o), "D") for o in onsets]))
    dd = pd.DatetimeIndex(dates).values.astype("datetime64[D]")
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
    o = np.empty_like(A)
    for j in range(A.shape[1]):
        s = np.bincount(c, weights=A[:, j], minlength=n)
        o[:, j] = A[:, j] - (s / cnt)[c]
    return o


def fit(d, onsets, spec, yresid=None):
    code, lag = assign(d["date"], onsets)
    keep = (code >= 0) | (np.abs(lag) > BASELINE_GAP)
    dd = d[keep]
    code = code[keep]
    if len(dd) < 2000:
        return None
    y = (dd["y"].values if yresid is None else yresid[keep]).astype(float)
    D = np.column_stack([(code == i).astype(float) for i in range(len(BINS))])
    S = seasonal_design(dd, spec)
    A = np.column_stack([y, D, S]) if S.size else np.column_stack([y, D])
    strata = pd.Categorical(dd["station_id"].astype(str) + "|"
                            + dd["winter"].astype(str)).codes
    Ad = demean(A, strata)
    try:
        beta, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
        return beta[:len(BINS)]
    except Exception:
        return None


def draw_clean(dates, doys, rng):
    by = {}
    for t in dates:
        by.setdefault(t.dayofyear, []).append(t)
    out = []
    for doy in rng.permutation(doys):
        c = by.get(int(doy))
        if c:
            out.append(c[rng.integers(len(c))])
    return pd.DatetimeIndex(out)


def main():
    d, meta = load_stations()
    print(f"station-level sample: {d['station_id'].nunique()} stations, "
          f"{len(d):,} station-days, {d['winter'].nunique()} winters")
    e = meta["elevation_m"]
    print(f"  elevation {e.min():.0f}-{e.max():.0f} m, "
          f"latitude {meta['latitude'].min():.1f}-{meta['latitude'].max():.1f} deg")

    het = heterogeneity(d)
    print(f"\nheterogeneity: seasonal amplitude "
          f"{het['seasonal_amplitude_degC']['min']}-"
          f"{het['seasonal_amplitude_degC']['max']} degC "
          f"(median {het['seasonal_amplitude_degC']['median']}); "
          f"corr(elev, amplitude) = {het['corr_elevation_vs_amplitude'][0]}, "
          f"corr(elev, level) = {het['corr_elevation_vs_mean_level'][0]}")

    onsets = load_catalogue("primary")
    onsets = onsets[(onsets >= d["date"].min()) & (onsets <= d["date"].max())]
    infl = G.real_influence_mask(pd.DatetimeIndex(d["date"]), onsets)
    clean = d[~infl].reset_index(drop=True)
    uniq = pd.DatetimeIndex(sorted(set(clean["date"])))
    doys = pd.DatetimeIndex(onsets).dayofyear.values
    print(f"events {len(onsets)}; removed {infl.sum():,} influenced "
          f"({100*infl.mean():.0f}%); clean {len(clean):,} rows")

    resid = clean["y"].values - station_climatology(clean)
    res = {"n_stations": int(clean["station_id"].nunique()),
           "n_rows": int(len(clean)), "n_winters": int(clean["winter"].nunique()),
           "n_events": int(len(onsets)), "heterogeneity": het,
           "n_mean_draws": N_MEAN, "n_cover_draws": N_COVER, "n_boot": N_BOOT,
           "specs": {}}

    sd = float(clean["y"].std())
    means = {}
    for spec in ("shared_harm3", "band_harm3", "station_clim"):
        yr = resid if spec == "station_clim" else None
        rng = np.random.default_rng(11)
        ests = []
        for _ in range(N_MEAN):
            fk = draw_clean(uniq, doys, rng)
            b = fit(clean, fk, spec, yresid=yr)
            if b is not None and np.isfinite(b[TEST_BIN]):
                ests.append(b[TEST_BIN])
        ests = np.array(ests)
        means[spec] = ests
        res["specs"][spec] = {
            "null_mean_degC": round(float(ests.mean()), 4),
            "null_mean_in_sd_units": round(float(ests.mean() / sd), 4),
            "null_sd": round(float(ests.std()), 4), "n_draws": int(len(ests))}
        print(f"  [means] {spec:14s} null={ests.mean():+.4f} degC "
              f"({ests.mean()/sd:+.3f} sd)  sd={ests.std():.3f}  n={len(ests)}")

    best = min(means, key=lambda s: abs(means[s].mean()))
    print(f"\n  least-biased: {best} -> coverage ({N_COVER} draws x {N_BOOT} reps)")
    yr = resid if best == "station_clim" else None
    rng = np.random.default_rng(23)
    cov = rej = 0
    for _ in range(N_COVER):
        fk = draw_clean(uniq, doys, rng)
        vals = []
        r2 = np.random.default_rng(int(rng.integers(1e6)))
        winters = np.array(sorted(clean["winter"].unique()))
        idx = {w: np.flatnonzero((clean["winter"] == w).values) for w in winters}
        for _ in range(N_BOOT):
            pick = r2.choice(winters, len(winters), replace=True)
            rows = np.concatenate([idx[w] for w in pick])
            sub = clean.iloc[rows].copy()
            sub["winter"] = np.concatenate(
                [np.full(len(idx[w]), j) for j, w in enumerate(pick)])
            b = fit(sub, fk, best, yresid=None if yr is None else yr[rows])
            if b is not None and np.isfinite(b[TEST_BIN]):
                vals.append(b[TEST_BIN])
        if len(vals) < 50:
            continue
        lo, hi = np.percentile(vals, 2.5), np.percentile(vals, 97.5)
        if lo <= 0 <= hi:
            cov += 1
        else:
            rej += 1
    n = cov + rej
    res["coverage_spec"] = best
    res["specs"][best].update({
        "coverage": round(cov / n, 3) if n else None,
        "coverage_CI95": clopper_pearson(cov, n) if n else None,
        "fpr": round(rej / n, 3) if n else None, "n_cover": int(n)})
    print(f"  [coverage] {best}: cov={res['specs'][best]['coverage']} "
          f"{res['specs'][best]['coverage_CI95']}  "
          f"FPR={res['specs'][best]['fpr']}  (n={n})")

    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
