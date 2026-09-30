#!/usr/bin/env python3
"""
r114_avalanche_event_study.py
=============================
Recompute the avalanche results (R106-R111) in EVENT-STUDY form.

WHY: R106-R111 each fitted ONE lag window at a time, so the control group for
any window was "all other days in that winter" -- which includes the post-onset
response period. R113 showed on the AO that this specification is not sound: a
window estimated that way is a contrast against a contaminated control group,
and its "placebo" is not a placebo. Every avalanche null therefore has to be
re-derived before it can be relied on.

DESIGN (identical to r113)
  Mutually exclusive lead/lag bins fitted SIMULTANEOUSLY against an explicit
  omitted baseline of winter days more than BASELINE_GAP days from ANY onset.
  Days are assigned to the bin of the NEAREST central date. Counts and binary
  indicators use conditional Poisson (stratum intercepts profiled out);
  inference is a winter-block bootstrap resampling whole winters -- and, in
  multi-region/multi-station analyses, all units of a winter together.

ARMS
  1. CH-Davos natural dry-slab counts        winter strata,        17 events
  2. Pooled Alpine (Davos + Austria LAWIS)   region x winter,      canonical cat
  3. SNOTEL continental instrumental physics state x winter, offset,  31 events
     (kinetic growth TG>10 K/m, melt, new-snow loading, rain-on-snow)

The question is narrow and specific: with a correctly specified event study,
is there ANY post-onset avalanche or snowpack response, and does it differ from
the pre-onset bins?

Output: data/results/r114_avalanche_event_study.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import condpois

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r114_avalanche_event_study.json"
SEASON = (11, 12, 1, 2, 3, 4)
N_HARM = 3
N_BOOT = 1000
BASELINE_GAP = 75

BINS = [(-45, -31), (-30, -16), (-15, -1), (0, 14), (15, 29), (30, 44)]
LABELS = [f"{a:+d}..{b:+d}" for a, b in BINS]

CONTINENTAL = {"CO", "UT", "WY", "MT", "ID", "NV", "NM", "AZ", "SD"}


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def add_bins(d, onsets, datecol="date"):
    """Nearest-onset lag -> mutually exclusive bin codes + baseline flag."""
    on = np.sort(np.array([np.datetime64(pd.Timestamp(o).tz_localize(None), "D")
                           for o in onsets]))
    dd = pd.DatetimeIndex(d[datecol]).tz_localize(None).values.astype("datetime64[D]")
    pos = np.searchsorted(on, dd)
    lo = np.clip(pos - 1, 0, len(on) - 1)
    hi = np.clip(pos, 0, len(on) - 1)
    dlo = (dd - on[lo]).astype(int)
    dhi = (dd - on[hi]).astype(int)
    lag = np.where(np.abs(dlo) <= np.abs(dhi), dlo, dhi)
    code = np.full(len(d), -1)
    for i, (a, b) in enumerate(BINS):
        code[(lag >= a) & (lag <= b)] = i
    d = d.copy()
    d["lag"] = lag
    d["bin"] = code
    d["is_baseline"] = np.abs(lag) > BASELINE_GAP
    return d[(d["bin"] >= 0) | d["is_baseline"]].copy()


def harmonics(d, datecol="date"):
    doy = pd.DatetimeIndex(d[datecol]).dayofyear.values
    for k in range(1, N_HARM + 1):
        d[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        d[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
    return d


def design(d, ycol, unitcol):
    y = d[ycol].values.astype(float)
    X = np.column_stack(
        [(d["bin"].values == i).astype(float) for i in range(len(BINS))]
        + [d[f"{t}{k}"].values for k in range(1, N_HARM + 1) for t in ("s", "c")])
    strata = pd.Categorical(
        d[unitcol].astype(str) + "|" + d["winter"].astype(str)).codes.astype(np.int64)
    return y, X, strata


def run(d, ycol, unitcol, label, offcol=None, seed=0):
    y, X, strata = design(d, ycol, unitcol)
    off = np.log(np.maximum(d[offcol].values.astype(float), 1e-9)) if offcol else None
    if offcol:
        keep = d[offcol].values > 0
        y, X, strata, off = y[keep], X[keep], strata[keep], off[keep]
        d = d[keep]
    b = condpois.fit(y, X, strata, offset=off)
    if b is None:
        return {"arm": label, "status": "unidentified"}

    rng = np.random.default_rng(seed)
    unit = pd.Categorical(d[unitcol]).codes.astype(np.int64)
    n_u = int(unit.max()) + 1
    winters = np.array(sorted(d["winter"].unique()))
    idx_by_w = [np.flatnonzero((d["winter"] == w).values) for w in winters]
    boot = []
    for _ in range(N_BOOT):
        pick = rng.integers(0, len(winters), len(winters))
        rows = np.concatenate([idx_by_w[k] for k in pick])
        st = np.concatenate([unit[idx_by_w[k]] + j * n_u for j, k in enumerate(pick)])
        v = condpois.fit(y[rows], X[rows], st,
                         offset=None if off is None else off[rows])
        if v is not None and np.all(np.isfinite(v[:len(BINS)])):
            boot.append(v[:len(BINS)])
    boot = np.array(boot)

    out = {"arm": label, "n_rows": int(len(d)),
           "n_winters": int(d["winter"].nunique()),
           "n_strata": int(d.groupby([unitcol, "winter"]).ngroups),
           "n_baseline_rows": int(d["is_baseline"].sum()), "bins": {}}
    for i, lab in enumerate(LABELS):
        col = boot[:, i] if len(boot) else np.array([])
        out["bins"][lab] = {
            "RR": round(float(np.exp(b[i])), 4),
            "CI95": [round(float(np.exp(np.percentile(col, 2.5))), 4),
                     round(float(np.exp(np.percentile(col, 97.5))), 4)] if len(col) else None,
            "p_two_sided": float(2 * min((col >= 0).mean(), (col <= 0).mean())) if len(col) else None,
        }
    print(f"\n[{label}]  {out['n_rows']:,} rows, {out['n_winters']} winters, "
          f"{out['n_strata']} strata, {out['n_baseline_rows']:,} baseline")
    for lab in LABELS:
        v = out["bins"][lab]
        flag = " *" if v["p_two_sided"] is not None and v["p_two_sided"] < 0.05 else ""
        print(f"   {lab:>9s}  RR={v['RR']:5.2f}  "
              f"[{v['CI95'][0]:.2f},{v['CI95'][1]:.2f}]  P={v['p_two_sided']:.4f}{flag}")
    return out


def davos():
    p = pd.read_parquet(ROOT / "data/processed/analysis_panel.parquet")
    if p.index.tz is not None:
        p.index = p.index.tz_convert("UTC").tz_localize(None)
    s = p["dry_natural_size_1234"].dropna().astype(float)
    s = s[np.isin(s.index.month, SEASON)]
    d = pd.DataFrame({"date": s.index, "count": s.values})
    d["winter"] = winter_of(pd.DatetimeIndex(d["date"]))
    d["unit"] = "davos"
    return d


def alpine():
    d1 = davos()
    lw = ROOT / "data/cryosphere/lawis_full/incidents_raw.json"
    parts = [d1]
    if lw.exists():
        rec = json.loads(lw.read_text(encoding="utf8"))
        dt = [pd.to_datetime(x["date"], utc=True, errors="coerce")
              for x in rec if x.get("date")
              and (x.get("location") or {}).get("country", {}).get("code") == "AT"]
        dt = pd.DatetimeIndex([x for x in dt if pd.notna(x)]).tz_convert("UTC").tz_localize(None).normalize()
        full = pd.date_range("1992-11-01", "2024-12-31", freq="D")
        full = full[np.isin(full.month, SEASON)]
        c = pd.Series(0.0, index=full)
        v = dt.value_counts()
        common = v.index.intersection(full)
        c.loc[common] = v.loc[common].astype(float)
        d2 = pd.DataFrame({"date": c.index, "count": c.values})
        d2["winter"] = winter_of(pd.DatetimeIndex(d2["date"]))
        d2["unit"] = "austria"
        parts.append(d2)
    return pd.concat(parts, ignore_index=True)


def snotel():
    sn = pd.read_parquet(ROOT / "data/processed/cryosphere/snotel_daily.parquet",
                         columns=["station_id", "snwd_mm", "tavg_c", "tmax_c", "prec_mm"])
    if not isinstance(sn.index, pd.DatetimeIndex):
        sn = sn.set_index("date")
    if sn.index.tz is not None:
        sn.index = sn.index.tz_convert("UTC").tz_localize(None)
    sn = sn[np.isin(sn.index.month, SEASON)]
    st = sn["station_id"].map(lambda x: str(x).split(":")[1] if ":" in str(x) else "??")
    sn = sn[st.isin(CONTINENTAL)]
    st = st[st.isin(CONTINENTAL)]
    ok = (sn["snwd_mm"] >= 200) & sn["tavg_c"].notna()
    tg = np.where(ok, (0.0 - sn["tavg_c"]) / (sn["snwd_mm"] / 1000.0).replace(0, np.nan), np.nan)
    sn = sn.assign(state=st.values, valid_tg=ok.astype(int),
                   kinetic=((tg > 10) & ok).astype(int),
                   melt=(sn["tmax_c"] > 0).astype(int),
                   new_snow=((sn["prec_mm"] > 10) & (sn["tavg_c"] < 0)).astype(int),
                   ros=((sn["prec_mm"] > 5) & (sn["tmax_c"] > 1)).astype(int))
    sn = sn.rename_axis("date").reset_index()
    agg = sn.groupby(["state", "date"]).agg(
        n_report=("station_id", "nunique"), n_valid_tg=("valid_tg", "sum"),
        kinetic=("kinetic", "sum"), melt=("melt", "sum"),
        new_snow=("new_snow", "sum"), ros=("ros", "sum")).reset_index()
    agg = agg[agg["n_report"] >= 5]
    agg["winter"] = winter_of(pd.DatetimeIndex(agg["date"]))
    return agg


def main():
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    onsets = pd.DatetimeIndex(pd.to_datetime(can["date"]).dropna())
    res = {"bins": LABELS, "baseline_gap_days": BASELINE_GAP, "arms": {}}

    print("=" * 74)
    d = harmonics(add_bins(davos(), onsets))
    res["arms"]["CH-Davos natural counts"] = run(d, "count", "unit", "CH-Davos natural counts")

    d = harmonics(add_bins(alpine(), onsets))
    res["arms"]["Pooled Alpine"] = run(d, "count", "unit", "Pooled Alpine")

    sn = harmonics(add_bins(snotel(), onsets))
    for col, off in (("kinetic", "n_valid_tg"), ("melt", "n_report"),
                     ("new_snow", "n_report"), ("ros", "n_report")):
        res["arms"][f"SNOTEL {col}"] = run(sn, col, "state", f"SNOTEL {col}", offcol=off)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
