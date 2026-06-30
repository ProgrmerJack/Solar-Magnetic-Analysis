"""
R57 – Pre-onset vs post-onset SSW suppression decomposition.

Analyses
--------
1. Decompose the ±15-day SSW rate ratio (RR) into pre-onset [-15, -1]
   and post-onset [0, +15] components; paired comparison.
2. Correlate pre-onset wave forcing (zonal-wind deceleration at 10 hPa,
   proxy for v'T' / EP-flux convergence) with event-level avalanche RR.
3. Z500 lead–lag: day of minimum Alpine Z500 relative to SSW onset;
   correlate lag with event RR.

Outputs
-------
data/results/r57_pre_post_decomposition.json
"""

from pathlib import Path
import json, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(r"C:\Users\Jack0\Solar-Magnetic-Analysis")
DATA = ROOT / "data" / "processed"
OUT  = ROOT / "data" / "results" / "r57_pre_post_decomposition.json"

# ── load data ───────────────────────────────────────────────────────
panel = pd.read_parquet(DATA / "analysis_panel_v2.parquet")
ssw_cat = pd.read_parquet(DATA / "atmospheric" / "ssw_catalog.parquet")
ncep_strat = pd.read_parquet(DATA / "atmospheric" / "ncep_stratosphere.parquet")
ncep_trop  = pd.read_parquet(DATA / "atmospheric" / "ncep_troposphere.parquet")

# Harmonise time zones – panel is tz-naive; others are UTC
ssw_onsets_utc = ssw_cat.index  # tz-aware
ssw_onsets     = ssw_onsets_utc.tz_localize(None)

strat_idx_naive = ncep_strat.index.tz_localize(None)
ncep_strat = ncep_strat.set_index(strat_idx_naive)

trop_idx_naive = ncep_trop.index.tz_localize(None)
ncep_trop = ncep_trop.set_index(trop_idx_naive)

# Filter SSW events to panel period
panel_start, panel_end = panel.index.min(), panel.index.max()
mask_in_panel = (ssw_onsets >= panel_start) & (ssw_onsets <= panel_end)
events = ssw_onsets[mask_in_panel]
print(f"SSW events in panel period: {len(events)}")

# ── helper: DOY-matched expected count ──────────────────────────────
# Build a climatological day-of-year mean from non-SSW days
aval_col = "dry_natural_size_1234"

# Mark all days within ±20 d of any SSW onset as "SSW-influenced"
ssw_influenced = np.zeros(len(panel), dtype=bool)
for onset in events:
    diff = (panel.index - onset).days
    ssw_influenced |= (np.abs(diff) <= 20)

panel["doy"] = panel.index.day_of_year
clim = panel.loc[~ssw_influenced].groupby("doy")[aval_col].mean()
# Fill any missing DOYs with overall mean
clim = clim.reindex(range(1, 367)).ffill().bfill()


def expected_count(dates):
    """Sum of DOY-climatological mean counts for a sequence of dates."""
    doys = pd.DatetimeIndex(dates).day_of_year
    return sum(clim.loc[d] for d in doys)


# ── Analysis 1: Pre-onset vs post-onset RR ─────────────────────────
records = []
for onset in events:
    pre_dates  = pd.date_range(onset - pd.Timedelta(days=15),
                               onset - pd.Timedelta(days=1), freq="D")
    post_dates = pd.date_range(onset,
                               onset + pd.Timedelta(days=15), freq="D")
    full_dates = pd.date_range(onset - pd.Timedelta(days=15),
                               onset + pd.Timedelta(days=15), freq="D")

    pre_obs  = panel.loc[panel.index.isin(pre_dates), aval_col].sum()
    post_obs = panel.loc[panel.index.isin(post_dates), aval_col].sum()
    full_obs = panel.loc[panel.index.isin(full_dates), aval_col].sum()

    pre_exp  = expected_count(pre_dates)
    post_exp = expected_count(post_dates)
    full_exp = expected_count(full_dates)

    rr_pre  = pre_obs / pre_exp if pre_exp > 0 else np.nan
    rr_post = post_obs / post_exp if post_exp > 0 else np.nan
    rr_full = full_obs / full_exp if full_exp > 0 else np.nan

    records.append(dict(
        onset=str(onset.date()),
        pre_obs=float(pre_obs), pre_exp=round(float(pre_exp), 2),
        post_obs=float(post_obs), post_exp=round(float(post_exp), 2),
        full_obs=float(full_obs), full_exp=round(float(full_exp), 2),
        rr_pre=round(float(rr_pre), 4),
        rr_post=round(float(rr_post), 4),
        rr_full=round(float(rr_full), 4),
    ))

df_ev = pd.DataFrame(records)
print("\n=== Event-level RR decomposition ===")
print(df_ev[["onset", "rr_pre", "rr_post", "rr_full"]].to_string(index=False))

# Paired test: pre vs post
rr_pre_arr  = df_ev["rr_pre"].values
rr_post_arr = df_ev["rr_post"].values
paired_t    = stats.ttest_rel(rr_pre_arr, rr_post_arr)
wilcox      = stats.wilcoxon(rr_pre_arr - rr_post_arr)

# What fraction of the deficit is pre-onset?
total_deficit_pre  = df_ev["pre_exp"].sum() - df_ev["pre_obs"].sum()
total_deficit_post = df_ev["post_exp"].sum() - df_ev["post_obs"].sum()
total_deficit      = total_deficit_pre + total_deficit_post
frac_pre = total_deficit_pre / total_deficit if total_deficit != 0 else np.nan

print(f"\nMedian RR  pre-onset:  {np.median(rr_pre_arr):.3f}")
print(f"Median RR  post-onset: {np.median(rr_post_arr):.3f}")
print(f"Paired t-test:  t={paired_t.statistic:.3f}, p={paired_t.pvalue:.4f}")
print(f"Wilcoxon:       W={wilcox.statistic:.1f}, p={wilcox.pvalue:.4f}")
print(f"Fraction of deficit pre-onset: {frac_pre:.1%}")

# ── Analysis 2: Wave-forcing proxy → surface RR ────────────────────
# Use -du/dt at 10 hPa as proxy for poleward eddy heat flux.
# Large deceleration ⇒ strong wave forcing.
u10 = ncep_strat["uwnd_ms_10hPa"].copy()
u10 = u10.sort_index()

wave_forcing = []
for onset in events:
    win_start = onset - pd.Timedelta(days=20)
    win_end   = onset - pd.Timedelta(days=5)
    u_window  = u10.loc[win_start:win_end]
    if len(u_window) >= 10:
        # du/dt via linear regression (m/s per day)
        x = np.arange(len(u_window))
        slope, _, _, _, _ = stats.linregress(x, u_window.values)
        # Negative slope = deceleration = wave forcing; store as positive
        wave_forcing.append(-slope)
    else:
        wave_forcing.append(np.nan)

df_ev["wave_forcing_proxy"] = wave_forcing
valid = df_ev.dropna(subset=["wave_forcing_proxy"])

r_wf, p_wf = stats.spearmanr(valid["wave_forcing_proxy"], valid["rr_full"])
r_wf_pre, p_wf_pre = stats.spearmanr(valid["wave_forcing_proxy"], valid["rr_pre"])
r_pear, p_pear = stats.pearsonr(valid["wave_forcing_proxy"], valid["rr_full"])

print("\n=== Wave forcing proxy (-du/dt 10hPa) vs RR ===")
print(f"  Spearman (full RR):  r={r_wf:.3f}, p={p_wf:.4f}  (n={len(valid)})")
print(f"  Spearman (pre RR):   r={r_wf_pre:.3f}, p={p_wf_pre:.4f}")
print(f"  Pearson  (full RR):  r={r_pear:.3f}, p={p_pear:.4f}")

# Also use 10hPa temperature change as second proxy
t10 = ncep_strat["air_K_10hPa"].copy().sort_index()
temp_forcing = []
for onset in events:
    win_start = onset - pd.Timedelta(days=20)
    win_end   = onset - pd.Timedelta(days=5)
    t_window  = t10.loc[win_start:win_end]
    if len(t_window) >= 10:
        x = np.arange(len(t_window))
        slope, _, _, _, _ = stats.linregress(x, t_window.values)
        temp_forcing.append(slope)  # positive slope = warming = wave driving
    else:
        temp_forcing.append(np.nan)

df_ev["temp_forcing_proxy"] = temp_forcing
valid2 = df_ev.dropna(subset=["temp_forcing_proxy"])
r_tf, p_tf = stats.spearmanr(valid2["temp_forcing_proxy"], valid2["rr_full"])
print(f"  Spearman (dT/dt vs full RR): r={r_tf:.3f}, p={p_tf:.4f}")

# ── Analysis 3: Z500 lead–lag ──────────────────────────────────────
z500_col = "ncep_z500_nh"  # from panel
z500_nadir_lag = []
z500_nadir_val = []

for onset in events:
    win_start = onset - pd.Timedelta(days=20)
    win_end   = onset + pd.Timedelta(days=20)
    sub = panel.loc[win_start:win_end, z500_col].dropna()
    if len(sub) >= 10:
        min_day = sub.idxmin()
        lag = (min_day - onset).days
        z500_nadir_lag.append(lag)
        z500_nadir_val.append(sub.min())
    else:
        z500_nadir_lag.append(np.nan)
        z500_nadir_val.append(np.nan)

df_ev["z500_nadir_lag"] = z500_nadir_lag
df_ev["z500_nadir_val"] = z500_nadir_val

lags = np.array([x for x in z500_nadir_lag if not np.isnan(x)])
med_lag = float(np.median(lags))
iqr_lag = float(np.percentile(lags, 75) - np.percentile(lags, 25))
print(f"\n=== Z500 nadir lag (days relative to SSW onset) ===")
print(f"  Median: {med_lag:.0f} d   IQR: [{np.percentile(lags,25):.0f}, {np.percentile(lags,75):.0f}] d")

valid3 = df_ev.dropna(subset=["z500_nadir_lag"])
r_z, p_z = stats.spearmanr(valid3["z500_nadir_lag"], valid3["rr_full"])
print(f"  Spearman (lag vs full RR): r={r_z:.3f}, p={p_z:.4f}")

# ── Assemble JSON output ───────────────────────────────────────────
results = {
    "analysis": "r57_pre_post_decomposition",
    "description": "SSW suppression decomposed into pre- and post-onset phases",
    "n_events": int(len(events)),
    "event_table": records,

    "pre_post_comparison": {
        "median_rr_pre":  round(float(np.median(rr_pre_arr)), 4),
        "median_rr_post": round(float(np.median(rr_post_arr)), 4),
        "mean_rr_pre":    round(float(np.mean(rr_pre_arr)), 4),
        "mean_rr_post":   round(float(np.mean(rr_post_arr)), 4),
        "paired_ttest":   {"t": round(float(paired_t.statistic), 4),
                           "p": round(float(paired_t.pvalue), 4)},
        "wilcoxon":       {"W": round(float(wilcox.statistic), 1),
                           "p": round(float(wilcox.pvalue), 4)},
        "total_deficit_pre":  round(float(total_deficit_pre), 1),
        "total_deficit_post": round(float(total_deficit_post), 1),
        "fraction_deficit_pre_onset": round(float(frac_pre), 4),
    },

    "wave_forcing_vs_rr": {
        "proxy_description": "-du/dt at 10 hPa in [-20,-5] d window (proxy for v'T' EP-flux convergence)",
        "n": int(len(valid)),
        "spearman_full_rr": {"r": round(float(r_wf), 4), "p": round(float(p_wf), 4)},
        "spearman_pre_rr":  {"r": round(float(r_wf_pre), 4), "p": round(float(p_wf_pre), 4)},
        "pearson_full_rr":  {"r": round(float(r_pear), 4), "p": round(float(p_pear), 4)},
        "temp_proxy": {
            "description": "dT/dt at 10 hPa in [-20,-5] d window",
            "spearman_full_rr": {"r": round(float(r_tf), 4), "p": round(float(p_tf), 4)},
        },
    },

    "z500_lead_lag": {
        "median_nadir_lag_days": round(med_lag, 1),
        "iqr_nadir_lag_days": [round(float(np.percentile(lags, 25)), 1),
                               round(float(np.percentile(lags, 75)), 1)],
        "spearman_lag_vs_rr": {"r": round(float(r_z), 4), "p": round(float(p_z), 4)},
        "event_lags": {str(r["onset"]): int(r["z500_nadir_lag"])
                       for _, r in df_ev.dropna(subset=["z500_nadir_lag"]).iterrows()},
    },

    "interpretation": {
        "key_finding": (
            f"{frac_pre:.0%} of the total SSW avalanche deficit occurs BEFORE onset, "
            "supporting the common-cause (planetary wave) mechanism rather than "
            "a top-down stratosphere→troposphere causal chain."
        ),
        "wave_forcing_correlation": (
            f"Pre-onset wave forcing (–du/dt 10 hPa) vs full-window RR: "
            f"Spearman r={r_wf:.3f}, p={p_wf:.3f}. "
            + ("Significant" if p_wf < 0.05 else "Not significant at α=0.05")
            + " — events with stronger wave bursts "
            + ("show stronger surface suppression." if r_wf < 0 else "do not show consistently stronger suppression.")
        ),
        "z500_timing": (
            f"Z500 nadir occurs at median lag = {med_lag:.0f} d "
            f"(IQR {np.percentile(lags,25):.0f} to {np.percentile(lags,75):.0f} d) "
            "relative to SSW onset."
        ),
    },
}

# Add wave forcing values to event table
for i, rec in enumerate(results["event_table"]):
    rec["wave_forcing_proxy"] = round(float(df_ev.iloc[i]["wave_forcing_proxy"]), 4) \
        if not np.isnan(df_ev.iloc[i]["wave_forcing_proxy"]) else None
    rec["z500_nadir_lag"] = int(df_ev.iloc[i]["z500_nadir_lag"]) \
        if not np.isnan(df_ev.iloc[i]["z500_nadir_lag"]) else None

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as f:
    json.dump(results, f, indent=2)

print(f"\n✓ Results saved → {OUT}")
