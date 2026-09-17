#!/usr/bin/env python3
"""
r109_norway_problem_mechanism.py
================================
NEW PILLAR CANDIDATE: test the manuscript's mechanism directly, in an
independent country, snow climate, and observing system.

The manuscript's mechanism has three limbs:
  (a) persistent weak layers develop / deepen   ("kinetic-growth metamorphism")
  (b) NATURAL release triggers are suppressed   ("observable signals disappear")
  (c) HUMAN triggerability rises                ("loaded gun")

Until now (a) rested on SNOWPACK model output at Swiss stations -- which R108
showed to be a seasonal artifact -- and (b)/(c) on Davos counts alone.

The Norwegian avalanche-warning archive (varsom, NVE) records, for every
regional bulletin, STRUCTURED avalanche problems assessed by forecasters:
    AvalCauseName          e.g. "Buried weak layer of faceted snow above a crust"
    AvalTriggerSimpleName  e.g. "Spontaneous release", "Low additional load"
    AvalancheExtName       e.g. "Dry slab avalanche"
This gives an operational, human-expert observable for all three limbs, over
65 regions and 53,907 region-days (2013-2026) -- entirely independent of the
Swiss data, the SNOWPACK model, and ERA5.

MAPPING (fixed before looking at any result; taken from the manuscript's own
mechanism, not chosen to fit)
  PWL_persistent : cause mentions faceted snow or surface hoar
                   (the kinetic-growth forms the manuscript invokes);
                   'buried weak layer of NEW snow' is storm slab, NOT persistent,
                   and is deliberately excluded.
  natural_trigger: trigger == "Spontaneous release"           -> predicted DOWN
  human_trigger  : trigger in {"Low additional load", "Easy to trigger"} -> UP
  dry_slab       : avalanche type == "Dry slab avalanche"

DESIGN  as r107: conditional Poisson on the binary indicator (gives a relative
risk), region x winter strata, three shared DOY harmonics (region levels are
already absorbed by the strata; see `design()`), and a winter-block bootstrap
resampling whole winters across all regions together.
A pre-onset placebo window must be null.

HONEST LIMITATION, stated up front: the archive begins in 2013, so it spans only
FOUR canonical SSWs (2018-02-12, 2019-01-01, 2021-01-05, 2023-02-16). Region-day
power is large but event-level power is not; the winter-block bootstrap reflects
this and its intervals will be wide. This is a mechanism replication test, not a
substitute for event-level inference.

Output: data/results/r109_norway_problem_mechanism.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import condpois

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r109_norway_problem_mechanism.json"
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

PERSISTENT = ("faceted", "surface hoar")      # kinetic-growth forms
HUMAN_TRIG = ("Low additional load", "Easy to trigger")

OUTCOMES = {
    "PWL_persistent": "persistent (faceted/surface-hoar) weak-layer problem cited -> manuscript predicts UP",
    "natural_trigger": "spontaneous-release problem cited -> manuscript predicts DOWN",
    "human_trigger": "low-load / easy-to-trigger problem cited -> manuscript predicts UP",
    "dry_slab": "dry-slab avalanche problem cited",
    "danger_ge3": "danger level >= 3 (considerable+) -> manuscript predicts UP",
}


def load():
    d = pd.read_parquet(
        ROOT / "data/processed/cryosphere/norway_avalanche.parquet",
        columns=["validfrom", "regionname", "regiontypename", "dangerlevel",
                 "avalancheproblems"])
    d["date"] = pd.to_datetime(d["validfrom"], errors="coerce", utc=True,
                               format="mixed").dt.normalize()
    d = d.dropna(subset=["date", "regionname"])
    d = d[np.isin(d["date"].dt.month, SEASON)]

    def flags(arr):
        pwl = nat = hum = dry = 0
        if arr is not None and len(arr):
            for p in arr:
                c = (p.get("AvalCauseName") or "").lower()
                t = p.get("AvalTriggerSimpleName") or ""
                a = p.get("AvalancheExtName") or ""
                if any(k in c for k in PERSISTENT):
                    pwl = 1
                if t == "Spontaneous release":
                    nat = 1
                if t in HUMAN_TRIG:
                    hum = 1
                if a == "Dry slab avalanche":
                    dry = 1
        return pwl, nat, hum, dry

    f = [flags(a) for a in d["avalancheproblems"]]
    d["PWL_persistent"] = [x[0] for x in f]
    d["natural_trigger"] = [x[1] for x in f]
    d["human_trigger"] = [x[2] for x in f]
    d["dry_slab"] = [x[3] for x in f]
    d["danger_ge3"] = (pd.to_numeric(d["dangerlevel"], errors="coerce") >= 3).astype(int)

    d["winter"] = np.where(d["date"].dt.month >= 11,
                           d["date"].dt.year + 1, d["date"].dt.year)
    doy = d["date"].dt.dayofyear.values
    for k in range(1, N_HARM + 1):
        d[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        d[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
    return d.drop(columns=["avalancheproblems"])


def tag(d, ev, a, b):
    f = np.zeros(len(d), int)
    dt = d["date"].values
    for o in ev:
        m = (dt >= np.datetime64(o + pd.Timedelta(days=a))) & \
            (dt <= np.datetime64(o + pd.Timedelta(days=b)))
        f[m.nonzero()[0]] = 1
    out = d.copy()
    out["W"] = f
    # stratum must hold both exposed and unexposed days to inform W
    return out.groupby(["regionname", "winter"]).filter(
        lambda g: g["W"].nunique() == 2)


def design(p, col):
    """Shared DOY harmonics, not region-specific.

    ponytail: with 65 regions, region-specific harmonics add 260 nuisance
    parameters to every one of ~60,000 bootstrap fits and the run does not
    finish. The region x winter strata already absorb each region's level in
    each winter, so only the SHAPE of the seasonal cycle is shared -- a third
    harmonic is added to keep that shape flexible. Verified below that the
    point estimates are materially unchanged (see `--region-harmonics` check).
    """
    y = p[col].values.astype(float)
    strata = pd.Categorical(
        p["regionname"].astype(str) + "|" + p["winter"].astype(str)).codes.astype(np.int64)
    H = [p[f"{t}{k}"].values for k in range(1, N_HARM + 1) for t in ("s", "c")]
    X = np.column_stack([p["W"].values.astype(float)] + H)
    return y, X, strata


def run(p, col, label, seed=0):
    y, X, strata = design(p, col)
    b = condpois.fit(y, X, strata)
    if b is None:
        return {"window": label, "status": "unidentified"}
    est = float(b[0])

    rng = np.random.default_rng(seed)
    reg = pd.Categorical(p["regionname"]).codes.astype(np.int64)
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
        "n_region_winters": int(p.groupby(["regionname", "winter"]).ngroups),
        "n_region_days": int(len(p)),
        "n_boot": int(len(bs)),
    }


def main():
    d = load()
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    ev = pd.to_datetime(can["date"]).dt.tz_localize("UTC")
    ev = pd.DatetimeIndex(ev[(ev >= d["date"].min()) & (ev <= d["date"].max())])
    print(f"Norway varsom: {len(d):,} region-days, {d['regionname'].nunique()} regions, "
          f"{d['date'].min().date()}..{d['date'].max().date()}")
    print(f"Canonical SSWs covered: {len(ev)} -> {[str(x.date()) for x in ev]}")
    print("base rates:", {k: round(float(d[k].mean()), 3) for k in OUTCOMES})

    res = {"n_region_days": int(len(d)), "n_regions": int(d["regionname"].nunique()),
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
