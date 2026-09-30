#!/usr/bin/env python3
"""
r105_multiregion_meta.py
========================
The venue-deciding test: is the post-SSW avalanche-activity suppression
DECISIVE when pooled across independent regions and observing systems?

Single-region tests are all underpowered (Swiss n=16 events; Tirol n=25 but
heavy within-winter clustering). Pooling is the honest way to gain power,
because SSWs are hemispheric: the SAME event is sampled independently by
every regional observing network.

REGIONS / OBSERVABLES (deliberately heterogeneous - convergent evidence):
  CH-Davos      natural dry-slab avalanche counts   1998-2019  (paper's core)
  AT-Tirol      LAWIS avalanche incidents           1992-2024  (independent Alps)
  US-continental CAIC accidents, continental states 1979-2024  (independent continent)
  US-maritime   CAIC accidents, maritime states     1979-2024  (negative control:
                                                     paper predicts weaker/absent)

METHOD
  1. Per region: Poisson GLM  count ~ SSW_window + C(winter) + DOY harmonics
     -> winter fixed effects absorb exposure/reporting trends and every
        winter-scale confounder; harmonics absorb the seasonal cycle.
  2. Per region SE from a WINTER-BLOCK BOOTSTRAP, not the model - avalanche
     counts cluster massively within winters, so model SEs are anti-conservative
     (demonstrated: Tirol model P=0.0004 -> bootstrap P=0.11).
  3. Pool with DerSimonian-Laird random effects (heterogeneous observables and
     regions => random, not fixed, effects). Report I^2 and Cochran's Q.

Window [+15,+44] d post-onset is taken from the r104 lag profile, which is the
pre-specified downward-coupling timescale (surface impact lags onset by 1-2 wk).
A leading [-30,-1] d window is also fitted as a PLACEBO: a genuine downward
-coupling effect must be absent before onset.

Output: data/results/r105_multiregion_meta.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r105_multiregion_meta.json"

SEASON = (11, 12, 1, 2, 3, 4)
CONTINENTAL = {"CO", "UT", "WY", "MT", "ID", "NV", "NM", "AZ", "SD"}
MARITIME = {"WA", "OR", "CA", "AK"}
N_BOOT = 2000
FORMULA = "count ~ W + s1 + c1 + s2 + c2 + C(winter)"


def detect_ssw(u):
    """Charlton-Polvani major mid-winter SSW central dates (identical to r102/r104)."""
    u = u.dropna().sort_index()
    dates, vals = u.index, u.values
    events, west_run = [], 0
    for k in range(len(u)):
        d, v = dates[k], vals[k]
        if v > 0:
            west_run += 1
            continue
        if d.month in SEASON and west_run >= 20:
            sey = d.year + 1 if d.month >= 11 else d.year
            fut = u[(u.index > d) & (u.index <= pd.Timestamp(year=sey, month=4, day=30, tz="UTC"))]
            maxrun = cur = 0
            for vv in fut.values:
                if vv > 0:
                    cur += 1
                    maxrun = max(maxrun, cur)
                else:
                    cur = 0
            if maxrun >= 10:
                events.append(d)
        west_run = 0
    return pd.DatetimeIndex(events)


def daily_from_dates(dates, lo, hi):
    full = pd.date_range(lo, hi, freq="D", tz="UTC")
    full = full[np.isin(full.month, SEASON)]
    c = pd.Series(0, index=full, dtype=float)
    v = pd.DatetimeIndex(dates).value_counts()
    common = v.index.intersection(full)
    c.loc[common] = v.loc[common].astype(float)
    return c


def make_frame(counts, ssw, a, b):
    df = pd.DataFrame({"count": counts.values}, index=counts.index).dropna()
    df["winter"] = np.where(df.index.month >= 11, df.index.year + 1, df.index.year)
    doy = df.index.dayofyear.values
    for k in (1, 2):
        df[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        df[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
    f = np.zeros(len(df), int)
    for o in ssw:
        m = (df.index >= o + pd.Timedelta(days=a)) & (df.index <= o + pd.Timedelta(days=b))
        f[m.nonzero()[0]] = 1
    df["W"] = f
    return df.groupby("winter").filter(
        lambda g: g["W"].nunique() == 2 and g["count"].sum() > 0)


def logirr(df):
    m = smf.glm(FORMULA, df, family=sm.families.Poisson()).fit()
    return float(m.params["W"])


def _blocks(df):
    """Per-winter (y, X_core) blocks; winter dummies are rebuilt per resample.

    ponytail: prebuilt numpy design instead of re-parsing the patsy formula on
    every bootstrap fit -- same estimator, ~50x faster over N_BOOT x regions.
    """
    core = ["W", "s1", "c1", "s2", "c2"]
    return [(g["count"].values.astype(float), g[core].values.astype(float),
             g["W"].nunique() == 2)
            for _, g in df.groupby("winter")]


def _fit_blocks(blocks):
    """Poisson GLM with winter fixed effects on stacked blocks -> coef on W."""
    ys, Xs, rows = [], [], []
    keep = [b for b in blocks if b[2]]
    if len(keep) < 5:
        return None
    n = sum(len(b[0]) for b in keep)
    D = np.zeros((n, len(keep)))
    i = 0
    for j, (y, X, _) in enumerate(keep):
        ys.append(y)
        Xs.append(X)
        D[i:i + len(y), j] = 1.0   # winter intercepts (no global const)
        i += len(y)
    Xf = np.hstack([np.vstack(Xs), D])
    yf = np.concatenate(ys)
    try:
        m = sm.GLM(yf, Xf, family=sm.families.Poisson()).fit()
        return float(m.params[0])
    except Exception:
        return None


def region_effect(name, counts, ssw, a, b, seed=0):
    """Point estimate + winter-block-bootstrap SE (honest clustered inference)."""
    df = make_frame(counts, ssw, a, b)
    if df.empty or df["W"].nunique() < 2:
        return None
    est = logirr(df)
    rng = np.random.default_rng(seed)
    blocks = _blocks(df)
    # check: the fast numpy design must reproduce the patsy fit exactly
    est_fast = _fit_blocks(blocks)
    assert est_fast is not None and abs(est_fast - est) < 1e-6, \
        f"{name}: fast design mismatch {est_fast} vs {est}"
    nb = len(blocks)
    boot = []
    for _ in range(N_BOOT):
        pick = rng.integers(0, nb, nb)
        v = _fit_blocks([blocks[k] for k in pick])
        if v is not None and np.isfinite(v):
            boot.append(v)
    boot = np.array(boot)
    se = float(np.std(boot, ddof=1))
    return {
        "region": name, "n_winters": int(df["winter"].nunique()),
        "n_events_total": int(df["count"].sum()),
        "IRR": round(float(np.exp(est)), 4),
        "log_irr": est, "se_winter_block": se,
        "CI95": [round(float(np.exp(np.percentile(boot, 2.5))), 4),
                 round(float(np.exp(np.percentile(boot, 97.5))), 4)],
        "p_one_sided_suppression": float((boot >= 0).mean()),
        "n_boot": int(len(boot)),
    }


def dersimonian_laird(est, se):
    """Random-effects pooling (DerSimonian-Laird)."""
    est, se = np.asarray(est), np.asarray(se)
    wf = 1 / se ** 2
    mu_f = (wf * est).sum() / wf.sum()
    Q = (wf * (est - mu_f) ** 2).sum()
    k = len(est)
    C = wf.sum() - (wf ** 2).sum() / wf.sum()
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = 1 / (se ** 2 + tau2)
    mu = (wr * est).sum() / wr.sum()
    se_mu = np.sqrt(1 / wr.sum())
    I2 = max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0.0
    from scipy import stats as st
    return {
        "pooled_IRR": round(float(np.exp(mu)), 4),
        "CI95": [round(float(np.exp(mu - 1.96 * se_mu)), 4),
                 round(float(np.exp(mu + 1.96 * se_mu)), 4)],
        "p_two_sided": float(2 * st.norm.sf(abs(mu / se_mu))),
        "p_one_sided_suppression": float(st.norm.cdf(mu / se_mu)),
        "tau2": round(float(tau2), 5), "Q": round(float(Q), 3),
        "I2_percent": round(float(I2), 1), "k_regions": int(k),
    }


def load_regions(ssw):
    hi = pd.Timestamp("2024-12-31", tz="UTC")
    out = {}

    # CH-Davos: natural dry-slab counts (the paper's core observable)
    p = pd.read_parquet(ROOT / "data/processed/analysis_panel.parquet")
    if p.index.tz is not None:
        p.index = p.index.tz_convert("UTC").tz_localize(None)
    p.index = p.index.tz_localize("UTC")
    s = p["dry_natural_size_1234"].dropna()
    s = s[np.isin(s.index.month, SEASON)]
    out["CH-Davos natural counts"] = s.astype(float)

    # AT-Tirol incidents
    t = pd.read_csv(ROOT / "data/cryosphere/austria_lawis/incidents_tirol.csv")
    t["date"] = pd.to_datetime(t["date"], errors="coerce", format="mixed", utc=True)
    t = t.dropna(subset=["date"])
    out["AT-Tirol incidents"] = daily_from_dates(
        t["date"], pd.Timestamp("1992-11-01", tz="UTC"), hi)

    # US CAIC accidents, split by snow climate
    c = pd.read_excel(ROOT / "data/cryosphere/caic/caic_accident_data.xlsx")
    c["Date"] = pd.to_datetime(c["Date"], errors="coerce", utc=True)
    c = c.dropna(subset=["Date"])
    lo = pd.Timestamp("1979-01-01", tz="UTC")
    c = c[(c["Date"] >= lo) & (c["Date"] <= hi)]
    out["US-continental accidents"] = daily_from_dates(
        c[c["State"].isin(CONTINENTAL)]["Date"], lo, hi)
    out["US-maritime accidents"] = daily_from_dates(
        c[c["State"].isin(MARITIME)]["Date"], lo, hi)
    return out


def run(regions, ssw, a, b, label):
    print(f"\n=== {label}: window [{a:+d},{b:+d}] d ===")
    rows = []
    for i, (name, counts) in enumerate(regions.items()):
        r = region_effect(name, counts, ssw, a, b, seed=i)
        if r is None:
            continue
        rows.append(r)
        print(f"  {name:28s} IRR={r['IRR']:5.2f}  CI[{r['CI95'][0]:.2f},{r['CI95'][1]:.2f}]  "
              f"winters={r['n_winters']:2d}  n={r['n_events_total']:6.0f}  "
              f"P(suppr)={r['p_one_sided_suppression']:.3f}")
    pooled = dersimonian_laird([r["log_irr"] for r in rows],
                               [r["se_winter_block"] for r in rows])
    print(f"  {'POOLED (random effects)':28s} IRR={pooled['pooled_IRR']:5.2f}  "
          f"CI[{pooled['CI95'][0]:.2f},{pooled['CI95'][1]:.2f}]  "
          f"P={pooled['p_two_sided']:.4f}  I2={pooled['I2_percent']}%")
    return {"per_region": rows, "pooled": pooled}


def main():
    strat = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet")
    if strat.index.tz is None:
        strat.index = strat.index.tz_localize("UTC")
    ssw = detect_ssw(strat["uwnd_ms_10hPa"])
    print(f"SSWs (Charlton-Polvani, NCEP 1979-2024): {len(ssw)}")

    regions = load_regions(ssw)
    res = {"n_ssw": len(ssw)}
    res["main_post_onset"] = run(regions, ssw, 15, 44, "MAIN post-onset")
    res["placebo_pre_onset"] = run(regions, ssw, -30, -1, "PLACEBO pre-onset")
    res["early_post_onset"] = run(regions, ssw, 0, 14, "early post-onset")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
