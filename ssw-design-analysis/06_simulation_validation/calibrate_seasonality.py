#!/usr/bin/env python3
"""
calibrate_seasonality.py
========================
GO/NO-GO CONDITION 3. Everything downstream is blocked until this passes.

THE PROBLEM
  In the superseded work, pseudo-onset randomisation returned a null distribution
  centred at -0.258 rather than 0, meaning the estimator reported an "effect"
  where none existed. Two candidate causes, which must be separated:

    (a) the seasonal adjustment (three annual harmonics) is too rigid to remove
        the true seasonal cycle, leaving structure that the lag bins absorb;

    (b) the pseudo-onset draw itself was wrong. The old test drew fake onsets
        uniformly from 1 Dec + 0..100 days, but real SSWs cluster in Jan-Feb.
        A null built from a different calendar distribution than the real events
        is not a null for those events -- it measures the seasonal cycle.

  (b) is a defect in the TEST, not the estimator, and would make the estimator
  look broken when it is not. It is checked first.

DESIGN
  Pseudo-onsets are drawn by PERMUTING the observed onset day-of-year values
  across winters. This reproduces the observed calendar distribution exactly, so
  any residual is attributable to the seasonal model rather than to calendar
  mismatch. Two draws are compared:
      matched    permute observed DOYs across winters   <- the correct null
      uniform    the old uniform Dec-Mar draw           <- reproduces the defect

SEASONAL SPECIFICATIONS COMPARED
  harm3      3 annual harmonics                        (the superseded default)
  harm6      6 annual harmonics
  cyclic_spline  cyclic cubic spline in day-of-year, knots fixed in advance
  lowo_clim  leave-one-winter-out day-of-year climatology, subtracted before
             fitting, so the seasonal estimate never sees the winter being tested

PASS CRITERIA (all four, else the estimator is not usable)
  1. mean estimate under the matched null within +/-0.05 of zero
  2. 95% interval coverage in 0.93-0.97
  3. false-positive rate in 0.03-0.07
  4. no material dependence on spline complexity

Output: ssw-design-analysis/06_simulation_validation/calibration_results.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "3_calibration"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]  # Solar-Magnetic-Analysis, which holds data/
OUT = RESULTS / "calibration_results.json"

SEASON = (11, 12, 1, 2, 3, 4)
BINS = [(-45, -31), (-30, -16), (-15, -1), (0, 14), (15, 29), (30, 44)]
LABELS = [f"{a:+d}..{b:+d}" for a, b in BINS]
BASELINE_GAP = 75
N_DRAWS = 400
TEST_BIN = 4          # +15..+29, the bin the old work reported as the response


# ----------------------------------------------------------------- data ------
def read_cpc(path):
    d = pd.read_csv(path, sep=r"\s+", header=None,
                    names=["year", "month", "day", "value"])
    d["date"] = pd.to_datetime(dict(year=d.year, month=d.month, day=d.day),
                               errors="coerce")
    return d.dropna(subset=["date"]).set_index("date")["value"].astype(float).sort_index()


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


# ------------------------------------------------- seasonal design blocks ----
def harmonics(doy, k):
    cols = []
    for j in range(1, k + 1):
        cols.append(np.sin(2 * np.pi * j * doy / 365.25))
        cols.append(np.cos(2 * np.pi * j * doy / 365.25))
    return np.column_stack(cols)


def cyclic_spline(doy, n_knots=8, period=365.25):
    """Cyclic cubic spline basis on day-of-year, knots fixed a priori."""
    knots = np.linspace(0, period, n_knots, endpoint=False)
    x = np.asarray(doy, float)
    cols = []
    for k in knots:
        d = np.abs((x - k + period / 2) % period - period / 2) / (period / n_knots)
        # compact cubic (Wendland-type) bump: smooth, local, cyclic by construction
        cols.append(np.where(d < 1, (1 - d) ** 3 * (3 * d + 1), 0.0))
    B = np.column_stack(cols)
    return B[:, 1:]          # drop one column for identifiability


def lowo_climatology(y, doy, winter):
    """Leave-one-winter-out day-of-year climatology, smoothed over +/-7 days."""
    df = pd.DataFrame({"y": y, "doy": doy, "w": winter})
    tot = df.groupby("doy")["y"].agg(["sum", "count"])
    out = np.empty(len(df))
    for w, idx in df.groupby("w").groups.items():
        sub = df.loc[idx]
        own = sub.groupby("doy")["y"].agg(["sum", "count"])
        rest = tot.subtract(own, fill_value=0)
        clim = (rest["sum"] / rest["count"].replace(0, np.nan))
        clim = clim.reindex(range(1, 367)).interpolate(limit_direction="both")
        clim = clim.rolling(15, center=True, min_periods=1).mean()
        out[df.index.get_indexer(idx)] = clim.reindex(sub["doy"]).values
    return out


SPECS = {
    "harm3": lambda d: harmonics(d["doy"].values, 3),
    "harm6": lambda d: harmonics(d["doy"].values, 6),
    "cyclic_spline": lambda d: cyclic_spline(d["doy"].values, n_knots=8),
    "cyclic_spline12": lambda d: cyclic_spline(d["doy"].values, n_knots=12),
    "lowo_clim": None,      # handled by residualising y, no design columns
}


# ------------------------------------------------------------- estimator -----
def assign_bins(dates, onsets):
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
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        s = np.bincount(c, weights=A[:, j], minlength=n)
        out[:, j] = A[:, j] - (s / cnt)[c]
    return out


def fit_once(d, onsets, spec):
    code, lag = assign_bins(d.index, onsets)
    keep = (code >= 0) | (np.abs(lag) > BASELINE_GAP)
    dd = d[keep]
    code = code[keep]
    y = dd["y"].values.astype(float)
    if spec == "lowo_clim":
        y = y - lowo_climatology(y, dd["doy"].values, dd["winter"].values)
        S = np.zeros((len(dd), 0))
    else:
        S = SPECS[spec](dd)
    D = np.column_stack([(code == i).astype(float) for i in range(len(BINS))])
    A = np.column_stack([y, D, S]) if S.size else np.column_stack([y, D])
    Ad = demean(A, dd["winter"].values)
    try:
        beta, *_ = np.linalg.lstsq(Ad[:, 1:], Ad[:, 0], rcond=None)
        return beta[:len(BINS)]
    except Exception:
        return None


def bootstrap_ci(d, onsets, spec, n_boot=200, seed=0):
    """Winter-block bootstrap interval for the tested bin."""
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
        b = fit_once(sub, onsets, spec)
        if b is not None and np.isfinite(b[TEST_BIN]):
            vals.append(b[TEST_BIN])
    if len(vals) < 30:
        return None
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# ---------------------------------------------------------- pseudo-onsets ----
def pseudo_matched(real_onsets, winters_avail, rng):
    """Permute the OBSERVED onset day-of-years across winters.

    Preserves the observed calendar distribution exactly, so a non-zero result
    cannot be blamed on calendar mismatch.
    """
    doys = pd.DatetimeIndex(real_onsets).dayofyear.values
    ws = rng.choice(winters_avail, len(doys), replace=False) \
        if len(winters_avail) >= len(doys) else \
        rng.choice(winters_avail, len(doys), replace=True)
    out = []
    for w, doy in zip(ws, rng.permutation(doys)):
        year = int(w) - 1 if doy > 250 else int(w)
        try:
            out.append(pd.Timestamp(year=year, month=1, day=1)
                       + pd.Timedelta(days=int(doy) - 1))
        except Exception:
            continue
    return pd.DatetimeIndex(out)


def pseudo_uniform(n, winters_avail, rng):
    """The superseded draw: uniform 1 Dec + 0..100 d. Reproduces the old defect."""
    ws = rng.choice(winters_avail, n, replace=True)
    return pd.DatetimeIndex([
        pd.Timestamp(year=int(w) - 1, month=12, day=1)
        + pd.Timedelta(days=int(rng.integers(0, 100))) for w in ws])


def main():
    ao = read_cpc(ROOT / "data/processed/atmospheric/ao_daily_cpc.txt")
    d = pd.DataFrame({"y": ao.values}, index=ao.index)
    d = d[np.isin(d.index.month, SEASON)]
    d["winter"] = winter_of(d.index)
    d["doy"] = d.index.dayofyear

    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    onsets = pd.DatetimeIndex(pd.to_datetime(can["date"]).dropna())
    onsets = onsets[(onsets >= d.index.min()) & (onsets <= d.index.max())]
    winters = np.array(sorted(d["winter"].unique()))
    print(f"AO {d.index.min().date()}..{d.index.max().date()}: {len(d):,} winter days, "
          f"{len(winters)} winters, {len(onsets)} real events")
    print(f"observed onset DOY: median {np.median(pd.DatetimeIndex(onsets).dayofyear):.0f}, "
          f"IQR {np.percentile(pd.DatetimeIndex(onsets).dayofyear,[25,75])}")

    res = {"n_draws": N_DRAWS, "tested_bin": LABELS[TEST_BIN], "specs": {}}
    for spec in SPECS:
        obs = fit_once(d, onsets, spec)
        entry = {"observed_profile": {l: round(float(v), 4)
                                      for l, v in zip(LABELS, obs)}}
        for draw_name, draw in (("matched", pseudo_matched),
                                ("uniform", pseudo_uniform)):
            rng = np.random.default_rng(7)
            ests, covers, rejects = [], 0, 0
            for t in range(N_DRAWS):
                fake = (draw(onsets, winters, rng) if draw_name == "matched"
                        else draw(len(onsets), winters, rng))
                if len(fake) < 5:
                    continue
                b = fit_once(d, fake, spec)
                if b is None or not np.isfinite(b[TEST_BIN]):
                    continue
                ests.append(b[TEST_BIN])
                if t < 60:      # coverage on a subsample: bootstrap is expensive
                    ci = bootstrap_ci(d, fake, spec, n_boot=120, seed=t)
                    if ci:
                        covers += (ci[0] <= 0 <= ci[1])
                        rejects += not (ci[0] <= 0 <= ci[1])
            ests = np.array(ests)
            ncov = covers + rejects
            entry[draw_name] = {
                "null_mean": round(float(ests.mean()), 4),
                "null_sd": round(float(ests.std()), 4),
                "abs_mean": round(float(abs(ests.mean())), 4),
                "coverage": round(covers / ncov, 3) if ncov else None,
                "false_positive_rate": round(rejects / ncov, 3) if ncov else None,
                "n": int(len(ests)), "n_coverage": int(ncov),
            }
            print(f"  {spec:16s} {draw_name:8s} null mean={entry[draw_name]['null_mean']:+.4f} "
                  f"sd={entry[draw_name]['null_sd']:.3f} "
                  f"cov={entry[draw_name]['coverage']} "
                  f"FPR={entry[draw_name]['false_positive_rate']}")
        m = entry["matched"]
        entry["passes"] = bool(
            m["abs_mean"] <= 0.05
            and m["coverage"] is not None and 0.93 <= m["coverage"] <= 0.97
            and m["false_positive_rate"] is not None
            and 0.03 <= m["false_positive_rate"] <= 0.07)
        res["specs"][spec] = entry

    ok = [s for s, v in res["specs"].items() if v["passes"]]
    res["passing_specs"] = ok
    print("\n" + "=" * 66)
    print("PASS (matched null centred, coverage 0.93-0.97, FPR 0.03-0.07):",
          ok or "NONE")
    OUT.write_text(json.dumps(res, indent=2), encoding="utf8")
    print(f"Saved -> {OUT}")


if __name__ == "__main__":
    main()
