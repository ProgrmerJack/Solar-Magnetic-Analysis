#!/usr/bin/env python3
"""
r110_albina_alpine_problems.py
==============================
Alpine replication of R109: does the manuscript's mechanism appear in
forecaster-assessed avalanche problems in the ALPS -- the domain where its
"planetary-wave Alpine blocking" mechanism is actually claimed to operate?

Data: data/cryosphere/albina_caaml/  (downloaded by
      scripts/download_extra/download_albina_caaml_problems.py)
      ALBINA / EUREGIO CAAML v6 bulletins, 97 regions in Tyrol, South Tyrol and
      Trentino, 2018-12-04 .. 2024-04-30, 117,047 EAWS-standard problems.

WHY THIS MATTERS
  * R108 falsified the manuscript's structural limb using SNOWPACK model output
    (seasonal artifact). ALBINA's `persistent_weak_layers` problem type is a
    forecaster-assessed PWL measure, completely independent of that model.
  * R109 found the trigger dissociation replicates in NORWAY. Norway is not the
    Alps, and the manuscript's mechanism is Alpine-specific, so an Alpine test
    is the one that matters for the manuscript's own claim.

PRE-SPECIFIED MAPPING (from EAWS problem definitions, not chosen on results)
  persistent_weak_layers : the manuscript's structural claim         -> UP
  spontaneous_types      : gliding_snow or wet_snow -- problems that release
                           without human input; the natural-trigger channel -> DOWN
  skier_triggerable      : persistent_weak_layers or wind_slab -- the classic
                           human-triggerable dry-slab problems              -> UP
  poor_stability         : snowpackStability in {poor, very_poor}
  high_frequency         : frequency in {some, many}  -- expected avalanche
                           abundance, the closest ALBINA proxy for activity -> DOWN
  danger_ge3             : max danger >= 3

DESIGN identical to R107/R109: conditional Poisson, region x winter strata,
shared DOY harmonics, winter-block bootstrap resampling whole winters across all
regions. Pre-onset placebo must be null.

HONEST LIMITATION: the archive spans six winters and FIVE canonical SSWs
(2019-01-01, 2021-01-05, 2023-02-16, 2024-01-16, 2024-03-04). Region-day counts
are large but the winter-block bootstrap has only six winters to resample, so
intervals are wide and this cannot carry event-level inference on its own.

Output: data/results/r110_albina_alpine_problems.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import condpois

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r110_albina_alpine_problems.json"
SEASON = (11, 12, 1, 2, 3, 4)
N_BOOT = 1000
N_HARM = 3

WINDOWS = {
    "placebo_pre": (-45, -16),
    "onset":       (-15, 15),
    "early_post":  (0, 14),
    "mid_post":    (15, 29),
    "late_post":   (30, 44),
    "post_15_44":  (15, 44),
}

SPONTANEOUS = {"gliding_snow", "wet_snow"}
SKIER = {"persistent_weak_layers", "wind_slab"}

OUTCOMES = {
    "persistent_weak_layers": "PWL problem cited -> manuscript predicts UP",
    "spontaneous_types": "gliding/wet-snow problem cited (natural channel) -> predicts DOWN",
    "skier_triggerable": "PWL or wind-slab problem cited (human channel) -> predicts UP",
    "poor_stability": "snowpack stability poor/very_poor",
    "high_frequency": "problem frequency some/many (activity proxy) -> predicts DOWN",
    "danger_ge3": "max danger >= 3",
}


def load():
    p = pd.read_csv(ROOT / "data/cryosphere/albina_caaml/albina_problems.csv")
    d = pd.read_csv(ROOT / "data/cryosphere/albina_caaml/albina_danger.csv")
    p["date"] = pd.to_datetime(p["date"], utc=True)
    d["date"] = pd.to_datetime(d["date"], utc=True)

    # collapse the (date, region) x problem rows to one row per region-day
    g = p.groupby(["date", "region"])
    agg = pd.DataFrame({
        "persistent_weak_layers": g["problem_type"].apply(
            lambda s: int((s == "persistent_weak_layers").any())),
        "spontaneous_types": g["problem_type"].apply(
            lambda s: int(s.isin(SPONTANEOUS).any())),
        "skier_triggerable": g["problem_type"].apply(
            lambda s: int(s.isin(SKIER).any())),
        "poor_stability": g["snowpack_stability"].apply(
            lambda s: int(s.isin(["poor", "very_poor"]).any())),
        "high_frequency": g["frequency"].apply(
            lambda s: int(s.isin(["some", "many"]).any())),
    }).reset_index()

    dd = d.groupby(["date", "region"], as_index=False)["danger_max"].max()
    out = agg.merge(dd, on=["date", "region"], how="left")
    out["danger_ge3"] = (pd.to_numeric(out["danger_max"], errors="coerce") >= 3).astype(int)
    out = out[np.isin(out["date"].dt.month, SEASON)].copy()

    out["winter"] = np.where(out["date"].dt.month >= 11,
                             out["date"].dt.year + 1, out["date"].dt.year)
    doy = out["date"].dt.dayofyear.values
    for k in range(1, N_HARM + 1):
        out[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        out[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
    return out


def tag(d, ev, a, b):
    f = np.zeros(len(d), int)
    dt = d["date"].values
    for o in ev:
        m = (dt >= np.datetime64(o + pd.Timedelta(days=a))) & \
            (dt <= np.datetime64(o + pd.Timedelta(days=b)))
        f[m.nonzero()[0]] = 1
    out = d.copy()
    out["W"] = f
    return out.groupby(["region", "winter"]).filter(lambda g: g["W"].nunique() == 2)


def design(p, col):
    y = p[col].values.astype(float)
    strata = pd.Categorical(
        p["region"].astype(str) + "|" + p["winter"].astype(str)).codes.astype(np.int64)
    H = [p[f"{t}{k}"].values for k in range(1, N_HARM + 1) for t in ("s", "c")]
    return y, np.column_stack([p["W"].values.astype(float)] + H), strata


def run(p, col, label, seed=0):
    y, X, strata = design(p, col)
    b = condpois.fit(y, X, strata)
    if b is None:
        return {"window": label, "status": "unidentified"}
    est = float(b[0])
    rng = np.random.default_rng(seed)
    reg = pd.Categorical(p["region"]).codes.astype(np.int64)
    n_reg = int(reg.max()) + 1
    winters = np.array(sorted(p["winter"].unique()))
    idx_by_w = [np.flatnonzero((p["winter"] == w).values) for w in winters]
    bs = []
    for _ in range(N_BOOT):
        pick = rng.integers(0, len(winters), len(winters))
        rows = np.concatenate([idx_by_w[k] for k in pick])
        st = np.concatenate([reg[idx_by_w[k]] + j * n_reg for j, k in enumerate(pick)])
        v = condpois.fit(y[rows], X[rows], st)
        if v is not None and np.isfinite(v[0]):
            bs.append(float(v[0]))
    bs = np.array(bs)
    return {
        "window": label,
        "RR": round(float(np.exp(est)), 4),
        "CI95": [round(float(np.exp(np.percentile(bs, 2.5))), 4),
                 round(float(np.exp(np.percentile(bs, 97.5))), 4)] if len(bs) else None,
        "p_one_sided_up": float((bs <= 0).mean()) if len(bs) else None,
        "p_one_sided_down": float((bs >= 0).mean()) if len(bs) else None,
        "base_rate": round(float(p[col].mean()), 4),
        "n_region_winters": int(p.groupby(["region", "winter"]).ngroups),
        "n_region_days": int(len(p)), "n_boot": int(len(bs)),
    }


def main():
    d = load()
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    ev = pd.to_datetime(can["date"]).dt.tz_localize("UTC")
    ev = pd.DatetimeIndex(ev[(ev >= d["date"].min()) & (ev <= d["date"].max())])
    print(f"ALBINA: {len(d):,} region-days, {d['region'].nunique()} regions, "
          f"{d['date'].min().date()}..{d['date'].max().date()}, "
          f"{d['winter'].nunique()} winters")
    print(f"Canonical SSWs covered: {len(ev)} -> {[str(x.date()) for x in ev]}")
    print("base rates:", {k: round(float(d[k].mean()), 3) for k in OUTCOMES})

    res = {"n_region_days": int(len(d)), "n_regions": int(d["region"].nunique()),
           "n_ssw": len(ev), "ssw_dates": [str(x.date()) for x in ev],
           "outcome_definitions": OUTCOMES, "results": {}}
    for col, desc in OUTCOMES.items():
        print(f"\n=== {col} :: {desc} ===")
        res["results"][col] = {}
        for lab, (a, b) in WINDOWS.items():
            p = tag(d, ev, a, b)
            if p.empty or p["W"].nunique() < 2:
                continue
            r = run(p, col, lab)
            res["results"][col][lab] = r
            if "RR" in r:
                print(f"  {lab:12s} RR={r['RR']:5.2f}  CI[{r['CI95'][0]:.2f},{r['CI95'][1]:.2f}]  "
                      f"P(up)={r['p_one_sided_up']:.3f} P(down)={r['p_one_sided_down']:.3f}  "
                      f"strata={r['n_region_winters']}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
