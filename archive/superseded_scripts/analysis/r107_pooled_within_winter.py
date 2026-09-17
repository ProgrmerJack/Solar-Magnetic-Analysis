#!/usr/bin/env python3
"""
r107_pooled_within_winter.py
============================
THE definitive test: is post-SSW avalanche suppression real, once every known
confounder is removed and inference is honest?

DESIGN (fixed in advance; nothing here is selected on the outcome)
------------------------------------------------------------------
Single pooled Poisson model over stacked region-days:

    count ~ W(lag window) + C(region x winter) + region-specific DOY harmonics

  * C(region x winter) absorbs, exactly, every region-winter-level confounder:
    the ~10x growth in backcountry recreation, changes in reporting practice,
    each winter's snowpack, and each region's baseline activity. Each winter in
    each region is therefore its own control -- the causally correct contrast,
    and the one the manuscript's between-winter design does not make.
  * Region-specific DOY harmonics let each region have its own seasonal cycle
    (SSWs and avalanches both peak mid-winter; without this the comparison is
    confounded by season).
  * A shared W coefficient gives the pooled incidence rate ratio.

INFERENCE
  Winter-block bootstrap resampling WHOLE WINTERS ACROSS ALL REGIONS TOGETHER.
  This is essential: Swiss and Austrian records share the same Alpine weather,
  so their arms are NOT independent. Resampling by winter preserves that
  cross-region correlation, which a per-region meta-analysis would ignore.
  Model standard errors are ~an order of magnitude too small here (verified:
  Tirol model P=0.0004 vs bootstrap P=0.11), so they are never used.

DOMAINS
  Primary   = ALPINE. The manuscript's mechanism is explicitly *Alpine blocking*,
              so the Alps are the pre-specified domain of the hypothesis.
  Secondary = US. A different blocking regime; reported as an out-of-domain
              check, not folded into the primary claim.

CATALOG
  data/processed/atmospheric/ssw_canonical.csv -- the authoritative NOAA CSL /
  Butler et al. (2017) compendium. The manuscript's own 16-event list contains
  2012-01-11 (not a major SSW in any reanalysis) and omits MAR 2000 / MAR 2010.

PLACEBO
  A pre-onset window must show no effect: downward coupling cannot act backwards.

Output: data/results/r107_pooled_within_winter.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import condpois

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r107_pooled_within_winter.json"
SEASON = (11, 12, 1, 2, 3, 4)
N_BOOT = 2000
N_HARM = 2

CONTINENTAL = {"CO", "UT", "WY", "MT", "ID", "NV", "NM", "AZ", "SD"}
MARITIME = {"WA", "OR", "CA", "AK"}

# lag windows (days from SSW central date)
WINDOWS = {
    "placebo_pre":  (-45, -16),
    "onset":        (-15, 15),     # the manuscript's window
    "early_post":   (0, 14),
    "mid_post":     (15, 29),
    "late_post":    (30, 44),
    "post_15_44":   (15, 44),      # main downward-coupling window
}


def winter_of(idx):
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def daily(dates, lo, hi):
    """Event timestamps -> daily counts over winter days in [lo, hi].

    .normalize() is essential: LAWIS API timestamps carry a time of day
    ("2018-02-17T10:49:00+01:00"), which would never intersect a midnight-based
    day index and would silently yield an all-zero series.
    """
    full = pd.date_range(lo, hi, freq="D", tz="UTC")
    full = full[np.isin(full.month, SEASON)]
    c = pd.Series(0.0, index=full)
    v = pd.DatetimeIndex(dates).normalize().value_counts()
    common = v.index.intersection(full)
    c.loc[common] = v.loc[common].astype(float)
    return c


def load_regions():
    """region -> daily count series. Heterogeneous observables on purpose."""
    hi = pd.Timestamp("2024-12-31", tz="UTC")
    R = {}

    # --- CH Davos: natural dry-slab counts (not accidents -> no exposure issue)
    p = pd.read_parquet(ROOT / "data/processed/analysis_panel.parquet")
    if p.index.tz is not None:
        p.index = p.index.tz_convert("UTC").tz_localize(None)
    p.index = p.index.tz_localize("UTC")
    s = p["dry_natural_size_1234"].dropna().astype(float)
    R["CH-Davos-natural"] = s[np.isin(s.index.month, SEASON)]

    # --- Austria: full LAWIS incident record (all provinces)
    lw = ROOT / "data/cryosphere/lawis_full/incidents_raw.json"
    if lw.exists():
        d = json.loads(lw.read_text(encoding="utf8"))
        rec = [(pd.to_datetime(x["date"], utc=True, errors="coerce"),
                (x.get("location") or {}).get("country", {}).get("code"))
               for x in d if x.get("date")]
        df = pd.DataFrame(rec, columns=["date", "cc"]).dropna()
        at = df[df["cc"] == "AT"]["date"]
        R["AT-Austria-LAWIS"] = daily(at, pd.Timestamp("1992-11-01", tz="UTC"), hi)
    else:  # fall back to the Tirol-only CSV if the full download has not run
        t = pd.read_csv(ROOT / "data/cryosphere/austria_lawis/incidents_tirol.csv")
        t["date"] = pd.to_datetime(t["date"], errors="coerce", format="mixed", utc=True)
        R["AT-Tirol-LAWIS"] = daily(t["date"].dropna(),
                                    pd.Timestamp("1992-11-01", tz="UTC"), hi)

    # --- US CAIC accidents, split by snow climate
    c = pd.read_excel(ROOT / "data/cryosphere/caic/caic_accident_data.xlsx")
    c["Date"] = pd.to_datetime(c["Date"], errors="coerce", utc=True)
    c = c.dropna(subset=["Date"])
    lo = pd.Timestamp("1979-01-01", tz="UTC")
    c = c[(c["Date"] >= lo) & (c["Date"] <= hi)]
    R["US-continental"] = daily(c[c["State"].isin(CONTINENTAL)]["Date"], lo, hi)
    R["US-maritime"] = daily(c[c["State"].isin(MARITIME)]["Date"], lo, hi)
    return R


def build_panel(regions, ev, a, b):
    """Stack regions into one frame with region x winter strata."""
    parts = []
    for name, s in regions.items():
        df = pd.DataFrame({"count": s.values}, index=s.index)
        df["region"] = name
        df["winter"] = winter_of(df.index)
        doy = df.index.dayofyear.values
        for k in range(1, N_HARM + 1):
            df[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
            df[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
        f = np.zeros(len(df), int)
        for o in ev:
            m = (df.index >= o + pd.Timedelta(days=a)) & (df.index <= o + pd.Timedelta(days=b))
            f[m.nonzero()[0]] = 1
        df["W"] = f
        # a stratum is informative only with both exposed and unexposed days
        df = df.groupby("winter").filter(
            lambda g: g["W"].nunique() == 2 and g["count"].sum() > 0)
        parts.append(df)
    return pd.concat(parts) if parts else pd.DataFrame()


def design(panel):
    """(y, X, strata): region x winter strata; region-specific DOY harmonics.

    Stratum intercepts are NOT columns here -- they are profiled out
    analytically by the conditional Poisson likelihood (see condpois.py).
    """
    y = panel["count"].values.astype(float)
    strata = pd.Categorical(
        panel["region"].astype(str) + "|" + panel["winter"].astype(str)).codes
    regs = pd.Categorical(panel["region"])
    H = []
    for ri in range(len(regs.categories)):
        sel = (regs.codes == ri).astype(float)
        for k in range(1, N_HARM + 1):
            H.append(panel[f"s{k}"].values * sel)
            H.append(panel[f"c{k}"].values * sel)
    X = np.column_stack([panel["W"].values.astype(float)] + H)
    return y, X, strata


def fit(panel):
    if panel.empty or panel["W"].nunique() < 2:
        return None
    y, X, strata = design(panel)
    b = condpois.fit(y, X, strata)
    if b is None or not np.isfinite(b[0]):
        return None
    return float(b[0])


def bootstrap(panel, n_boot=N_BOOT, seed=0):
    """Resample WHOLE WINTERS, taking every region's data for that winter together.

    Swiss and Austrian records share the same Alpine weather, so their arms are
    not independent; resampling by winter preserves that cross-region
    correlation. Operates on prebuilt numpy blocks -- rebuilding DataFrames per
    replicate dominates the runtime otherwise.
    """
    rng = np.random.default_rng(seed)
    y, X, _ = design(panel)
    # int64: Categorical.codes is int8, and code + j*n_reg overflows it
    reg = pd.Categorical(panel["region"]).codes.astype(np.int64)
    winters = np.array(sorted(panel["winter"].unique()))
    idx_by_w = [np.flatnonzero((panel["winter"] == w).values) for w in winters]
    n_reg = int(reg.max()) + 1

    out = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(winters), len(winters))
        rows = np.concatenate([idx_by_w[k] for k in pick])
        # distinct stratum per (region, resampled-winter-slot)
        strata = np.concatenate([reg[idx_by_w[k]] + j * n_reg
                                 for j, k in enumerate(pick)])
        b = condpois.fit(y[rows], X[rows], strata)
        if b is not None and np.isfinite(b[0]):
            out.append(float(b[0]))
    return np.array(out)


def run(regions, ev, label, a, b, n_boot=N_BOOT):
    panel = build_panel(regions, ev, a, b)
    if panel.empty:
        return {"window": label, "status": "EMPTY"}
    est = fit(panel)
    if est is None:
        return {"window": label, "status": "FIT_FAILED"}
    bs = bootstrap(panel, n_boot)
    r = {
        "window": label, "lag_days": [a, b],
        "IRR": round(float(np.exp(est)), 4),
        "CI95": [round(float(np.exp(np.percentile(bs, 2.5))), 4),
                 round(float(np.exp(np.percentile(bs, 97.5))), 4)] if len(bs) else None,
        "p_one_sided_suppression": float((bs >= 0).mean()) if len(bs) else None,
        "n_region_winters": int(panel.groupby(["region", "winter"]).ngroups),
        "n_days": int(len(panel)),
        "n_events_observed": float(panel["count"].sum()),
        "regions": sorted(panel["region"].unique().tolist()),
        "n_boot": int(len(bs)),
    }
    print(f"  {label:14s} IRR={r['IRR']:5.2f}  CI[{r['CI95'][0]:.2f},{r['CI95'][1]:.2f}]  "
          f"P(suppr)={r['p_one_sided_suppression']:.4f}  "
          f"strata={r['n_region_winters']}")
    return r


def main():
    can = pd.read_csv(ROOT / "data/processed/atmospheric/ssw_canonical.csv")
    ev = pd.DatetimeIndex(pd.to_datetime(can["date"]).dt.tz_localize("UTC")).sort_values()
    print(f"Canonical SSW catalog: {len(ev)} events {ev.min().date()}..{ev.max().date()}")

    regions = load_regions()
    for k, v in regions.items():
        print(f"  {k:20s} {len(v):6d} winter days, {v.sum():7.0f} events, "
              f"{v.index.min().date()}..{v.index.max().date()}")
        # guard: a date-parsing/timezone mismatch silently yields ~0 events
        assert v.sum() >= 50, (
            f"{k}: only {v.sum():.0f} events over {len(v)} winter days -- "
            "date alignment is almost certainly broken")

    alpine = {k: v for k, v in regions.items() if k.startswith(("CH", "AT"))}
    us = {k: v for k, v in regions.items() if k.startswith("US")}

    res = {"catalog": "ssw_canonical (NOAA CSL / Butler et al. 2017)",
           "n_ssw": len(ev), "design": "Poisson, region x winter FE, "
           "region-specific DOY harmonics, winter-block bootstrap across regions",
           "domains": {}}

    for dom, regs in [("ALPINE (primary, mechanism domain)", alpine),
                      ("US (secondary, out-of-domain)", us),
                      ("ALL regions", regions)]:
        print(f"\n=== {dom} ===")
        res["domains"][dom] = {lab: run(regs, ev, lab, a, b)
                               for lab, (a, b) in WINDOWS.items()}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
