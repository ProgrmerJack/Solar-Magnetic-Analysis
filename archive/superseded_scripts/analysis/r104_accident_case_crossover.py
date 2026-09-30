#!/usr/bin/env python3
"""
r104_accident_case_crossover.py
===============================
Breaks the n=16 ceiling on the paper's HUMAN-TRIGGER arm.

The n=16 limit was never physical - it came from the Swiss daily-count record
(1998/99-2018/19). Avalanche ACCIDENT records run far longer:
  * CAIC (US, 10 states)   1951-2025 -> covers ALL 40 NCEP-era SSWs
  * LAWIS Tirol (Austria)  1992-2024 -> covers 29 SSWs, independent Alpine region

Accidents ARE the human-trigger arm - the arm the paper predicts goes UP while
natural activity goes DOWN - and they are the societally relevant outcome.

DESIGN (case-crossover / indirect standardisation; each winter is its own control):
  count_day ~ ssw_window + C(winter) + DOY harmonics,  Negative Binomial
    - C(winter) absorbs the exposure trend (recreation grew ~10x since 1950s)
      and every winter-scale confounder (snowpack, participation, reporting).
    - DOY harmonics absorb the seasonal cycle (accidents AND SSWs both peak
      mid-winter; without this the comparison is confounded).
    - Inference by winter-block bootstrap - accidents cluster within winters.

Pre-specified from the paper (NOT chosen post hoc):
  window     = [0, +30] d after SSW onset (surface-impact window, as in r102)
  direction  = accident rate UP (one-sided hypothesis stated in the manuscript)
  geography  = continental > maritime (paper's continental-specificity claim,
               same state classification as r102)

Output: data/results/r104_accident_case_crossover.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "results" / "r104_accident_case_crossover.json"

WIN = 30            # [0, +WIN] d post-onset, pre-specified in r102
N_BOOT = 2000
SEASON = (11, 12, 1, 2, 3, 4)

# same classification as r102 (paper's continental-specificity hypothesis)
CONTINENTAL = {"CO", "UT", "WY", "MT", "ID", "NV", "NM", "AZ", "SD"}
MARITIME = {"WA", "OR", "CA", "AK"}


def detect_ssw(u):
    """Charlton-Polvani major mid-winter SSW central dates (identical to r102)."""
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


def winter_id(idx):
    """Winter label: Nov-Dec belong to the following calendar year's season."""
    return np.where(idx.month >= 11, idx.year + 1, idx.year)


def daily_counts(dates, lo, hi):
    """Accident dates -> complete daily count series over winter days in [lo, hi]."""
    full = pd.date_range(lo, hi, freq="D", tz="UTC")
    full = full[np.isin(full.month, SEASON)]
    c = pd.Series(0, index=full, dtype=int)
    v = dates.value_counts()
    v.index = pd.DatetimeIndex(v.index)
    common = v.index.intersection(full)
    c.loc[common] = v.loc[common].astype(int)
    return c


def build_frame(counts, ssw, win=WIN):
    df = pd.DataFrame({"count": counts.values}, index=counts.index)
    df["winter"] = winter_id(df.index)
    doy = df.index.dayofyear.values
    # 2 harmonics capture the winter accident season shape without overfitting
    for k in (1, 2):
        df[f"s{k}"] = np.sin(2 * np.pi * k * doy / 365.25)
        df[f"c{k}"] = np.cos(2 * np.pi * k * doy / 365.25)
    flag = np.zeros(len(df), dtype=int)
    for onset in ssw:
        m = (df.index >= onset) & (df.index <= onset + pd.Timedelta(days=win))
        flag[m.nonzero()[0]] = 1
    df["ssw_window"] = flag
    # a winter contributes information only if it has both exposed and
    # unexposed days AND at least one accident (conditional/case-crossover logic)
    keep = df.groupby("winter").filter(
        lambda g: g["ssw_window"].nunique() == 2 and g["count"].sum() > 0)
    return keep


def fit_nb(df):
    """NB GLM with winter FE + DOY harmonics; alpha via Cameron-Trivedi."""
    f = "count ~ ssw_window + s1 + c1 + s2 + c2 + C(winter)"
    pois = smf.glm(f, df, family=sm.families.Poisson()).fit()
    mu = pois.fittedvalues
    y = pois.model.endog
    aux = ((y - mu) ** 2 - y) / mu
    alpha = float(np.clip(sm.OLS(aux, mu).fit().params.iloc[0], 1e-6, 10))
    m = smf.glm(f, df, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
    return m, alpha, pois


def boot_irr(df, n_boot=N_BOOT, seed=0):
    """Winter-block bootstrap of the SSW-window IRR."""
    rng = np.random.default_rng(seed)
    winters = df["winter"].unique()
    out = []
    for _ in range(n_boot):
        pick = rng.choice(winters, size=len(winters), replace=True)
        parts = []
        for j, w in enumerate(pick):
            g = df[df["winter"] == w].copy()
            g["winter"] = j  # relabel so resampled winters stay distinct strata
            parts.append(g)
        b = pd.concat(parts)
        b = b.groupby("winter").filter(lambda g: g["ssw_window"].nunique() == 2)
        if b["winter"].nunique() < 5:
            continue
        try:
            m = smf.glm("count ~ ssw_window + s1 + c1 + s2 + c2 + C(winter)",
                        b, family=sm.families.Poisson()).fit()
            out.append(float(np.exp(m.params["ssw_window"])))
        except Exception:
            continue
    return np.array(out)


def per_event_oe(counts, ssw, win=WIN):
    """Indirect standardisation: observed vs DOY+winter-expected, per event."""
    df = pd.DataFrame({"count": counts.values}, index=counts.index)
    df["winter"] = winter_id(df.index)
    df["doy"] = df.index.dayofyear
    # climatological DOY share (pooled across all winters), smoothed
    clim = df.groupby("doy")["count"].mean()
    clim = clim.reindex(range(1, 367)).interpolate(limit_direction="both")
    clim = clim.rolling(15, center=True, min_periods=1).mean()

    rows = []
    for onset in ssw:
        w = winter_id(pd.DatetimeIndex([onset]))[0]
        g = df[df["winter"] == w]
        if len(g) < 60 or g["count"].sum() == 0:
            continue
        m = (g.index >= onset) & (g.index <= onset + pd.Timedelta(days=win))
        if m.sum() < win * 0.5 or (~m).sum() < 30:
            continue
        obs = float(g.loc[m, "count"].sum())
        # expected = winter total x window's climatological DOY share
        share = clim.reindex(g.loc[m, "doy"]).sum() / clim.reindex(g["doy"]).sum()
        exp = float(g["count"].sum() * share)
        if exp <= 0:
            continue
        rows.append({"onset": str(onset.date()), "winter": int(w),
                     "observed": obs, "expected": round(exp, 2),
                     "oe_ratio": round(obs / exp, 3)})
    return rows


def analyse(name, dates, ssw, lo, hi):
    counts = daily_counts(dates, lo, hi)
    df = build_frame(counts, ssw)
    if df.empty or df["ssw_window"].sum() == 0:
        return {"dataset": name, "status": "NO_OVERLAP"}
    m, alpha, _ = fit_nb(df)
    irr = float(np.exp(m.params["ssw_window"]))
    bs = boot_irr(df)
    rows = per_event_oe(counts, ssw)
    oe = np.array([r["oe_ratio"] for r in rows])
    n_up = int((oe > 1).sum())
    tot_o = sum(r["observed"] for r in rows)
    tot_e = sum(r["expected"] for r in rows)

    res = {
        "dataset": name,
        "n_accidents_total": int(counts.sum()),
        "n_winters_informative": int(df["winter"].nunique()),
        "n_ssw_events_tested": len(rows),
        "n_accidents_in_ssw_windows": int(tot_o),
        "nb_alpha": round(alpha, 4),
        "IRR_ssw_window": round(irr, 4),
        "IRR_model_CI95": [round(float(np.exp(m.conf_int().loc["ssw_window", 0])), 4),
                           round(float(np.exp(m.conf_int().loc["ssw_window", 1])), 4)],
        "p_model_two_sided": float(m.pvalues["ssw_window"]),
        "IRR_bootstrap_CI95": [round(float(np.percentile(bs, 2.5)), 4),
                               round(float(np.percentile(bs, 97.5)), 4)] if len(bs) else None,
        "p_bootstrap_one_sided_up": float((bs <= 1).mean()) if len(bs) else None,
        "n_boot_ok": int(len(bs)),
        "observed_vs_expected": {
            "total_observed": round(tot_o, 1), "total_expected": round(tot_e, 1),
            "ratio": round(tot_o / tot_e, 4) if tot_e else None,
            "poisson_exact_p_one_sided": float(
                stats.poisson.sf(tot_o - 1, tot_e)) if tot_e else None,
            "events_above_expectation": f"{n_up}/{len(oe)}",
            "sign_test_p_one_sided": float(
                stats.binomtest(n_up, len(oe), 0.5, alternative="greater").pvalue) if len(oe) else None,
        },
        "per_event": rows,
    }
    print(f"\n[{name}] winters={res['n_winters_informative']} events={res['n_ssw_events_tested']} "
          f"accidents={res['n_accidents_total']} (in-window {res['n_accidents_in_ssw_windows']})")
    print(f"  IRR = {irr:.3f}  model CI {res['IRR_model_CI95']}  P={res['p_model_two_sided']:.4g}")
    if res["IRR_bootstrap_CI95"]:
        print(f"  winter-block bootstrap CI {res['IRR_bootstrap_CI95']}  "
              f"one-sided P={res['p_bootstrap_one_sided_up']:.4g}")
    ove = res["observed_vs_expected"]
    print(f"  O/E = {ove['total_observed']}/{ove['total_expected']} = {ove['ratio']}  "
          f"Poisson P={ove['poisson_exact_p_one_sided']:.4g}  "
          f"events up {ove['events_above_expectation']} sign P={ove['sign_test_p_one_sided']:.4g}")
    return res


def main():
    strat = pd.read_parquet(ROOT / "data/processed/atmospheric/ncep_stratosphere.parquet")
    if strat.index.tz is None:
        strat.index = strat.index.tz_localize("UTC")
    ssw = detect_ssw(strat["uwnd_ms_10hPa"])
    print(f"SSWs detected (Charlton-Polvani, NCEP 1979-2024): {len(ssw)}")

    results = {"n_ssw": len(ssw), "window_days": WIN,
               "ssw_dates": [str(d.date()) for d in ssw], "arms": {}}

    # ---- CAIC (US, 1951-2025) --------------------------------------------
    caic = pd.read_excel(ROOT / "data/cryosphere/caic/caic_accident_data.xlsx")
    caic["Date"] = pd.to_datetime(caic["Date"], errors="coerce", utc=True)
    caic = caic.dropna(subset=["Date"])
    lo, hi = pd.Timestamp("1979-01-01", tz="UTC"), pd.Timestamp("2024-12-31", tz="UTC")
    caic = caic[(caic["Date"] >= lo) & (caic["Date"] <= hi)]

    results["arms"]["CAIC_all_US"] = analyse(
        "CAIC all US", caic["Date"], ssw, lo, hi)
    results["arms"]["CAIC_continental"] = analyse(
        "CAIC continental", caic[caic["State"].isin(CONTINENTAL)]["Date"], ssw, lo, hi)
    results["arms"]["CAIC_maritime"] = analyse(
        "CAIC maritime", caic[caic["State"].isin(MARITIME)]["Date"], ssw, lo, hi)

    # ---- LAWIS Tirol (Austria, 1992-2024) --------------------------------
    tirol = pd.read_csv(ROOT / "data/cryosphere/austria_lawis/incidents_tirol.csv")
    tirol["date"] = pd.to_datetime(tirol["date"], errors="coerce", format="mixed", utc=True)
    tirol = tirol.dropna(subset=["date"])
    tlo = pd.Timestamp("1992-11-01", tz="UTC")
    results["arms"]["LAWIS_Tirol"] = analyse(
        "LAWIS Tirol", tirol["date"], ssw, tlo, hi)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
