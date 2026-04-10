"""
05_enhanced_analysis.py — Deep mechanism analysis for Nature Geoscience
======================================================================
Extends the primary result (RR=0.774, p=0.008 for 1-3d fast pathway)
with mechanistic tests:

1. SOC power-law analysis (flare + avalanche size distributions)
2. SSW-specific avalanche coupling (Butler 2017 catalog)
3. Forbush decrease mechanism (OMNI proton flux as GCR proxy)
4. Precipitation mediation test (SNOTEL)
5. Strong-storm dose-response (Kp>=7 / Dst<=-100)
6. Lag sweep (0-30d day-by-day RR)
"""
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, LOG, load_panel

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


# ═════════════════════════════════════════════════════════════════════════════
# Helper: same NB GLM as primary endpoint
# ═════════════════════════════════════════════════════════════════════════════
def nb_glm(df, outcome, exposure, covariates=None):
    """Run NB GLM and return rate ratio + CI + p."""
    df = df.copy()
    df["y"] = np.round(df[outcome]).astype(int)
    df.loc[df["y"] < 0, "y"] = 0

    covs = [exposure]
    base_covs = ["day_of_season", "day_of_season_sq"]
    optional = ["nao_daily", "qbo_u50", "ncep_z500_nh", "ncep_slp_nh"]
    for c in base_covs + optional:
        if c in df.columns:
            covs.append(c)
    if covariates:
        covs.extend([c for c in covariates if c in df.columns])

    clean = df[["y"] + covs].dropna()
    if len(clean) < 50 or clean["y"].sum() == 0:
        return None

    Y = clean["y"]
    X = sm.add_constant(clean[covs])

    try:
        model = sm.GLM(Y, X, family=sm.families.NegativeBinomial(alpha=1.0))
        result = model.fit(maxiter=200, disp=0)
    except Exception:
        try:
            model = sm.GLM(Y, X, family=sm.families.Poisson())
            result = model.fit(maxiter=200, disp=0)
        except Exception:
            return None

    beta = float(result.params[exposure])
    pval = float(result.pvalues[exposure])
    ci = result.conf_int(alpha=0.05).loc[exposure]
    return {
        "rate_ratio": float(np.exp(beta)),
        "rr_ci_lower": float(np.exp(ci.iloc[0])),
        "rr_ci_upper": float(np.exp(ci.iloc[1])),
        "p_value": pval,
        "n_obs": int(result.nobs),
        "beta": beta,
    }


# ═════════════════════════════════════════════════════════════════════════════
# 1. SOC Power-Law Analysis
# ═════════════════════════════════════════════════════════════════════════════
def soc_power_law():
    """Fit power-law to flare energies and avalanche sizes using MLE."""
    print("\n" + "=" * 70)
    print("1. SOC POWER-LAW ANALYSIS")
    print("=" * 70)

    results = {}

    # --- Flare energy proxy: GOES class → peak flux ---
    fl = pd.read_parquet(PROCESSED / "solar" / "flares.parquet")
    class_map = {"A": 1e-8, "B": 1e-7, "C": 1e-6, "M": 1e-5, "X": 1e-4}

    def class_to_flux(ct):
        if pd.isna(ct) or not isinstance(ct, str):
            return np.nan
        ct = ct.strip()
        base = class_map.get(ct[0].upper(), np.nan)
        if np.isnan(base):
            return np.nan
        try:
            return base * float(ct[1:])
        except (ValueError, IndexError):
            return base

    fl["peak_flux"] = fl["classType"].apply(class_to_flux)
    flare_sizes = fl["peak_flux"].dropna().values
    flare_sizes = flare_sizes[flare_sizes > 0]

    # MLE power-law fit (Clauset et al. 2009 method)
    # For data x >= x_min: alpha = 1 + n / sum(ln(x/x_min))
    x_min_flare = np.percentile(flare_sizes, 10)  # Use 10th percentile as x_min
    above = flare_sizes[flare_sizes >= x_min_flare]
    n = len(above)
    alpha_flare = 1 + n / np.sum(np.log(above / x_min_flare))
    alpha_err_flare = (alpha_flare - 1) / np.sqrt(n)  # SE from Clauset eq. 3.2

    print(f"\nFlare peak flux (GOES catalog):")
    print(f"  N = {len(flare_sizes)}, above x_min = {n}")
    print(f"  x_min = {x_min_flare:.2e} W/m^2")
    print(f"  alpha = {alpha_flare:.3f} +/- {alpha_err_flare:.3f}")
    print(f"  (Lu & Hamilton 1991 reference: alpha ~ 1.8)")

    results["flare_power_law"] = {
        "n_total": len(flare_sizes),
        "n_above_xmin": int(n),
        "x_min": float(x_min_flare),
        "alpha": float(alpha_flare),
        "alpha_se": float(alpha_err_flare),
        "reference_alpha": 1.8,
    }

    # --- Avalanche size: area_m2 from SLF ---
    sl = pd.read_parquet(PROCESSED / "cryosphere" / "slf_snow_events.parquet")
    aval_sizes = sl["area_m2"].dropna().values
    aval_sizes = aval_sizes[aval_sizes > 0]

    x_min_aval = np.percentile(aval_sizes, 10)
    above_a = aval_sizes[aval_sizes >= x_min_aval]
    n_a = len(above_a)
    alpha_aval = 1 + n_a / np.sum(np.log(above_a / x_min_aval))
    alpha_err_aval = (alpha_aval - 1) / np.sqrt(n_a)

    print(f"\nAvalanche area (SLF snow events):")
    print(f"  N = {len(aval_sizes)}, above x_min = {n_a}")
    print(f"  x_min = {x_min_aval:.0f} m^2")
    print(f"  alpha = {alpha_aval:.3f} +/- {alpha_err_aval:.3f}")
    print(f"  (Typical avalanche SOC: alpha ~ 1.3-2.0)")

    results["avalanche_power_law"] = {
        "n_total": len(aval_sizes),
        "n_above_xmin": int(n_a),
        "x_min": float(x_min_aval),
        "alpha": float(alpha_aval),
        "alpha_se": float(alpha_err_aval),
    }

    # --- Compare exponents ---
    z_diff = abs(alpha_flare - alpha_aval) / np.sqrt(alpha_err_flare**2 + alpha_err_aval**2)
    from scipy import stats
    p_diff = 2 * (1 - stats.norm.cdf(z_diff))
    print(f"\nExponent comparison:")
    print(f"  alpha_flare - alpha_aval = {alpha_flare - alpha_aval:.3f}")
    print(f"  z = {z_diff:.2f}, p = {p_diff:.4f}")
    print(f"  {'DIFFERENT' if p_diff < 0.05 else 'COMPATIBLE'} exponents")

    results["exponent_comparison"] = {
        "delta_alpha": float(alpha_flare - alpha_aval),
        "z_stat": float(z_diff),
        "p_value": float(p_diff),
        "compatible": p_diff >= 0.05,
    }

    # --- Avalanche SOC during active vs quiet geomagnetic periods ---
    sl_idx = sl.index.tz_localize(None) if sl.index.tz else sl.index
    panel = load_panel(winter_only=True)
    panel.index = panel.index.tz_localize(None) if panel.index.tz else panel.index
    active_days = set(panel[panel["post_event_0_30d"] == 1].index.date)
    quiet_days = set(panel[panel["post_event_0_30d"] == 0].index.date)

    active_mask = np.array([d.date() in active_days for d in sl_idx])
    quiet_mask = np.array([d.date() in quiet_days for d in sl_idx])

    for label, mask in [("active_30d", active_mask), ("quiet", quiet_mask)]:
        subset = aval_sizes[mask[:len(aval_sizes)]]
        subset = subset[subset >= x_min_aval]
        if len(subset) > 20:
            a = 1 + len(subset) / np.sum(np.log(subset / x_min_aval))
            se = (a - 1) / np.sqrt(len(subset))
            print(f"  alpha ({label}): {a:.3f} +/- {se:.3f} (N={len(subset)})")
            results[f"avalanche_soc_{label}"] = {
                "alpha": float(a), "alpha_se": float(se), "n": int(len(subset))
            }

    return results


# ═════════════════════════════════════════════════════════════════════════════
# 2. SSW-Avalanche Coupling (Butler 2017 catalog)
# ═════════════════════════════════════════════════════════════════════════════
def ssw_coupling():
    """Test whether SSW events are followed by elevated avalanche activity."""
    print("\n" + "=" * 70)
    print("2. SSW-AVALANCHE COUPLING")
    print("=" * 70)

    panel = load_panel(winter_only=True)
    panel.index = panel.index.tz_localize(None) if panel.index.tz else panel.index
    ssw = pd.read_parquet(PROCESSED / "atmospheric" / "ssw_catalog.parquet")
    ssw.index = ssw.index.tz_localize(None) if ssw.index.tz else ssw.index

    results = {}

    # Filter SSW events within panel range
    ssw_in_range = ssw[(ssw.index >= panel.index.min()) & (ssw.index <= panel.index.max())]
    print(f"SSW events in panel range: {len(ssw_in_range)} / {len(ssw)}")

    if len(ssw_in_range) < 3:
        print("Too few SSW events for analysis")
        return {"error": "Too few SSW events in panel range"}

    # Create SSW post-event windows (0-15d, 15-30d, 30-60d)
    for window_name, d_start, d_end in [
        ("ssw_0_15d", 0, 15), ("ssw_15_30d", 15, 30),
        ("ssw_30_60d", 30, 60), ("ssw_0_60d", 0, 60)
    ]:
        panel[window_name] = 0
        for sd in ssw_in_range.index:
            m = (panel.index >= sd + pd.Timedelta(days=d_start)) & \
                (panel.index <= sd + pd.Timedelta(days=d_end))
            panel.loc[m, window_name] = 1

    # Test each window
    for outcome in ["aai_all_natural", "natural_size_234"]:
        if outcome not in panel.columns:
            continue
        sub = panel[panel[outcome].notna()].copy()
        print(f"\n  Outcome: {outcome} (N={len(sub)})")

        for window_name in ["ssw_0_15d", "ssw_15_30d", "ssw_30_60d", "ssw_0_60d"]:
            n_exposed = int(sub[window_name].sum())
            res = nb_glm(sub, outcome, window_name)
            if res:
                sig = "**" if res["p_value"] < 0.01 else "*" if res["p_value"] < 0.05 else ""
                print(f"    {window_name}: RR={res['rate_ratio']:.3f} "
                      f"[{res['rr_ci_lower']:.3f}, {res['rr_ci_upper']:.3f}] "
                      f"p={res['p_value']:.4f} {sig} (exposed days={n_exposed})")
                results[f"{outcome}_{window_name}"] = res
                results[f"{outcome}_{window_name}"]["exposed_days"] = n_exposed

    # Superposed epoch around SSW events
    print("\n  Superposed Epoch (mean AAI around SSW onset):")
    aai = panel["aai_all_natural"].dropna()
    ssw_dates = ssw_in_range.index
    for lag_start, lag_end, label in [(-15, -1, "pre"), (0, 7, "d0-7"),
                                       (8, 21, "d8-21"), (22, 45, "d22-45"),
                                       (46, 60, "d46-60")]:
        values = []
        for sd in ssw_dates:
            window = aai[(aai.index >= sd + pd.Timedelta(days=lag_start)) &
                         (aai.index <= sd + pd.Timedelta(days=lag_end))]
            if len(window) > 0:
                values.append(window.mean())
        if values:
            mean_v = np.mean(values)
            se_v = np.std(values) / np.sqrt(len(values))
            print(f"    {label:>8s}: mean AAI={mean_v:.2f} +/- {se_v:.2f} (n_ssw={len(values)})")
            results[f"sea_{label}"] = {"mean": float(mean_v), "se": float(se_v), "n": len(values)}

    return results


# ═════════════════════════════════════════════════════════════════════════════
# 3. Forbush Decrease Mechanism Test
# ═════════════════════════════════════════════════════════════════════════════
def forbush_mechanism():
    """Test whether the fast-pathway effect is mediated by GCR/precipitation."""
    print("\n" + "=" * 70)
    print("3. FORBUSH DECREASE MECHANISM TEST")
    print("=" * 70)

    panel = load_panel(winter_only=True)
    panel.index = panel.index.tz_localize(None) if panel.index.tz else panel.index
    results = {}

    # Load OMNI hourly for proton flux (GCR proxy — higher flux = more GCR)
    omni = pd.read_parquet(PROCESSED / "solar" / "omni_hourly.parquet")
    omni.index = omni.index.tz_localize(None) if omni.index.tz else omni.index

    # Proton flux at >1 MeV and >10 MeV as solar energetic particle proxies
    # During Forbush decreases, GCR decreases but SEP (solar proton) increases
    flux_cols = [c for c in omni.columns if "proton_flux" in c.lower()]
    print(f"Available proton flux columns: {flux_cols[:6]}")

    # Use flow pressure and Dst as storm proxies in daily panel
    # Compute daily OMNI summaries
    omni_daily = omni[["Dst", "Kp", "flow_speed", "Bz_GSM", "flow_pressure"]].resample("D").agg({
        "Dst": "min",
        "Kp": "max",
        "flow_speed": "max",
        "Bz_GSM": "min",  # Most southward
        "flow_pressure": "max",
    })
    omni_daily.columns = ["omni_dst_min", "omni_kp_max", "omni_vmax",
                           "omni_bz_min", "omni_pressure_max"]

    # Merge with panel
    merged = panel.join(omni_daily, how="left")

    # Test: does Bz-south (stronger Forbush) predict stronger avalanche decrease?
    # Create Bz-south storm indicator
    merged["strong_bz_south"] = (merged["omni_bz_min"] < -10).astype(int)
    merged["moderate_event"] = ((merged["kp_max"] >= 5) & (merged["kp_max"] < 7)).astype(int)
    merged["strong_event"] = (merged["kp_max"] >= 7).astype(int)

    # Create post-Bz-south windows
    bz_south_dates = merged.index[merged["strong_bz_south"] == 1]
    merged["post_bz_south_1_3d"] = 0
    for d in bz_south_dates:
        m = (merged.index > d) & (merged.index <= d + pd.Timedelta(days=3))
        merged.loc[m, "post_bz_south_1_3d"] = 1

    sub = merged[merged["aai_all_natural"].notna()].copy()

    # Test Bz-south → avalanche
    print(f"\nBz-south storms (Bz < -10 nT): {int(merged['strong_bz_south'].sum())} days")
    res = nb_glm(sub, "aai_all_natural", "post_bz_south_1_3d")
    if res:
        print(f"  Post-Bz-south (1-3d): RR={res['rate_ratio']:.3f} "
              f"[{res['rr_ci_lower']:.3f}, {res['rr_ci_upper']:.3f}] p={res['p_value']:.4f}")
        results["bz_south_fast"] = res

    # Wet avalanche specific test with Bz-south
    if "aai_all_wet" in sub.columns:
        res_wet = nb_glm(sub[sub["aai_all_wet"].notna()], "aai_all_wet", "post_bz_south_1_3d")
        if res_wet:
            print(f"  Post-Bz-south wet aval (1-3d): RR={res_wet['rate_ratio']:.3f} "
                  f"[{res_wet['rr_ci_lower']:.3f}, {res_wet['rr_ci_upper']:.3f}] "
                  f"p={res_wet['p_value']:.4f}")
            results["bz_south_fast_wet"] = res_wet

    # Dose-response: is Kp strength correlated with RR magnitude?
    print("\nDose-response by storm intensity:")
    for label, threshold_col, threshold in [
        ("Kp>=5 (moderate)", "kp_max", 5),
        ("Kp>=6 (strong)", "kp_max", 6),
        ("Kp>=7 (intense)", "kp_max", 7),
        ("Dst<=-50 (moderate)", "dst_min", None),
        ("Dst<=-100 (intense)", "dst_min", None),
    ]:
        merged_tmp = sub.copy()
        if "Dst" in label:
            thresh = -50 if "50" in label else -100
            storm_dates = merged_tmp.index[merged_tmp["dst_min"] <= thresh]
        else:
            storm_dates = merged_tmp.index[merged_tmp[threshold_col] >= threshold]

        merged_tmp["dose_post_1_3d"] = 0
        for d in storm_dates:
            m = (merged_tmp.index > d) & (merged_tmp.index <= d + pd.Timedelta(days=3))
            merged_tmp.loc[m, "dose_post_1_3d"] = 1

        n_exposed = int(merged_tmp["dose_post_1_3d"].sum())
        res = nb_glm(merged_tmp, "aai_all_natural", "dose_post_1_3d")
        if res:
            sig = "**" if res["p_value"] < 0.01 else "*" if res["p_value"] < 0.05 else ""
            print(f"  {label}: RR={res['rate_ratio']:.3f} "
                  f"[{res['rr_ci_lower']:.3f}, {res['rr_ci_upper']:.3f}] "
                  f"p={res['p_value']:.4f} {sig} (exposed={n_exposed})")
            results[f"dose_{label.split()[0]}"] = res

    return results


# ═════════════════════════════════════════════════════════════════════════════
# 4. Day-by-Day Lag Sweep (0-30d)
# ═════════════════════════════════════════════════════════════════════════════
def lag_sweep():
    """Compute rate ratios for each individual lag day 0-30."""
    print("\n" + "=" * 70)
    print("4. LAG SWEEP (Day-by-day RR, 0-30d)")
    print("=" * 70)

    panel = load_panel(winter_only=True)
    panel.index = panel.index.tz_localize(None) if panel.index.tz else panel.index
    sub = panel[panel["aai_all_natural"].notna()].copy()

    event_dates = sub.index[sub["geo_event"] == 1]
    results = {}

    for lag in range(0, 31):
        sub[f"lag_{lag}d"] = 0
        for ed in event_dates:
            target = ed + pd.Timedelta(days=lag)
            if target in sub.index:
                sub.loc[target, f"lag_{lag}d"] = 1

        res = nb_glm(sub, "aai_all_natural", f"lag_{lag}d")
        if res:
            sig = "**" if res["p_value"] < 0.01 else "*" if res["p_value"] < 0.05 else ""
            print(f"  Lag {lag:2d}d: RR={res['rate_ratio']:.3f} "
                  f"[{res['rr_ci_lower']:.3f}, {res['rr_ci_upper']:.3f}] "
                  f"p={res['p_value']:.4f} {sig}")
            results[f"lag_{lag}d"] = res

        sub.drop(columns=[f"lag_{lag}d"], inplace=True)

    return results


# ═════════════════════════════════════════════════════════════════════════════
# 5. SNOTEL Precipitation Mediation
# ═════════════════════════════════════════════════════════════════════════════
def precipitation_mediation():
    """Test if precipitation mediates the geomagnetic-avalanche link."""
    print("\n" + "=" * 70)
    print("5. PRECIPITATION MEDIATION TEST (SNOTEL)")
    print("=" * 70)

    panel = load_panel(winter_only=True)
    panel.index = panel.index.tz_localize(None) if panel.index.tz else panel.index
    results = {}

    # Load SNOTEL daily — aggregate to mean SWE change across stations
    snotel = pd.read_parquet(PROCESSED / "cryosphere" / "snotel_daily.parquet")
    snotel.index = snotel.index.tz_localize(None) if snotel.index.tz else snotel.index

    # Compute daily mean SWE change (proxy for new snowfall)
    # SNOTEL uses 'wteq_mm' for Snow Water Equivalent
    swe_col = [c for c in snotel.columns if c in ["wteq_mm", "swe_mm"]]
    if not swe_col:
        swe_col = [c for c in snotel.columns if "swe" in c.lower() or "wteq" in c.lower()]
    prec_col = [c for c in snotel.columns if "prec" in c.lower()]
    print(f"SNOTEL columns: {list(snotel.columns)}")
    print(f"SWE column: {swe_col}, Precip column: {prec_col}")

    # Use precipitation directly if available (better proxy for new snow)
    use_col = prec_col[0] if prec_col else (swe_col[0] if swe_col else None)
    if use_col:
        # Daily mean across all stations
        daily_swe = snotel.groupby(snotel.index.date)[use_col].mean()
        daily_swe.index = pd.DatetimeIndex(daily_swe.index, name="date")
        # For SWE, compute change; for precipitation, use directly
        if "prec" in use_col.lower():
            daily_swe_change = daily_swe
            daily_swe_change.name = "snotel_precip"
            mediator_col = "snotel_precip"
        else:
            daily_swe_change = daily_swe.diff()
            daily_swe_change.name = "snotel_swe_change"
            mediator_col = "snotel_swe_change"
        panel = panel.join(daily_swe_change, how="left")

        sub = panel[panel["aai_all_natural"].notna() & panel[mediator_col].notna()].copy()
        print(f"Panel with SNOTEL: {len(sub)} days")

        # Test: does controlling for precipitation diminish the geomagnetic signal?
        # Model A: avalanche ~ post_event_1_3d + standard confounders
        res_a = nb_glm(sub, "aai_all_natural", "post_event_1_3d")
        # Model B: avalanche ~ post_event_1_3d + precip + confounders
        res_b = nb_glm(sub, "aai_all_natural", "post_event_1_3d",
                       covariates=[mediator_col])

        if res_a and res_b:
            print(f"\n  Model A (no precip control): RR={res_a['rate_ratio']:.3f} p={res_a['p_value']:.4f}")
            print(f"  Model B (with precip ctrl):  RR={res_b['rate_ratio']:.3f} p={res_b['p_value']:.4f}")
            if np.log(res_a['rate_ratio']) != 0:
                attenuation = 1 - (np.log(res_b['rate_ratio']) / np.log(res_a['rate_ratio']))
            else:
                attenuation = 0
            print(f"  Attenuation of log(RR): {attenuation:.1%}")
            print(f"  (>50% attenuation suggests precipitation mediates the effect)")

            results["model_a_no_precip"] = res_a
            results["model_b_with_precip"] = res_b
            results["attenuation_pct"] = float(attenuation * 100)

        # Also test: does geomagnetic activity predict precipitation directly?
        print(f"\n  Direct test: geomagnetic event -> {mediator_col}:")
        sub2 = sub.copy()
        from scipy import stats as sp_stats
        exposed = sub2.loc[sub2["post_event_1_3d"] == 1, mediator_col].dropna()
        unexposed = sub2.loc[sub2["post_event_1_3d"] == 0, mediator_col].dropna()
        if len(exposed) > 10 and len(unexposed) > 10:
            t_stat, t_pval = sp_stats.ttest_ind(exposed, unexposed, equal_var=False)
            print(f"    Exposed (1-3d post-event): mean SWE change = {exposed.mean():.2f} mm")
            print(f"    Unexposed:                 mean SWE change = {unexposed.mean():.2f} mm")
            print(f"    Welch t = {t_stat:.3f}, p = {t_pval:.4f}")
            results["swe_change_test"] = {
                "exposed_mean": float(exposed.mean()),
                "unexposed_mean": float(unexposed.mean()),
                "t_stat": float(t_stat),
                "p_value": float(t_pval),
            }

    return results


# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════
def main():
    all_results = {}

    # Run each analysis sequentially to avoid OOM
    all_results["soc_power_law"] = soc_power_law()
    all_results["ssw_coupling"] = ssw_coupling()
    all_results["forbush_mechanism"] = forbush_mechanism()
    all_results["lag_sweep"] = lag_sweep()
    all_results["precipitation_mediation"] = precipitation_mediation()

    # Save comprehensive results
    out = RESULTS / "enhanced_analysis.json"
    out.write_text(json.dumps(all_results, indent=2, default=str), encoding="utf-8")
    print(f"\nAll results saved to {out}")

    # Final synthesis
    print("\n" + "=" * 70)
    print("SYNTHESIS OF ENHANCED ANALYSIS")
    print("=" * 70)

    # SOC
    soc = all_results["soc_power_law"]
    print(f"\nSOC Power Laws:")
    print(f"  Flare alpha = {soc['flare_power_law']['alpha']:.3f} +/- {soc['flare_power_law']['alpha_se']:.3f}")
    print(f"  Avalanche alpha = {soc['avalanche_power_law']['alpha']:.3f} +/- {soc['avalanche_power_law']['alpha_se']:.3f}")
    ec = soc["exponent_comparison"]
    print(f"  Comparison: z={ec['z_stat']:.2f}, p={ec['p_value']:.4f} -> {'Compatible' if ec['compatible'] else 'Different'}")

    # SSW
    ssw_r = all_results["ssw_coupling"]
    print(f"\nSSW-Avalanche Coupling:")
    for key in sorted(ssw_r.keys()):
        if key.startswith("aai"):
            r = ssw_r[key]
            sig = "**" if r["p_value"] < 0.01 else "*" if r["p_value"] < 0.05 else ""
            print(f"  {key}: RR={r['rate_ratio']:.3f} p={r['p_value']:.4f} {sig}")

    # Lag sweep summary
    lag_r = all_results["lag_sweep"]
    sig_lags = [(k, v) for k, v in lag_r.items() if v["p_value"] < 0.05]
    print(f"\nLag Sweep: {len(sig_lags)} significant lags out of 31")
    for k, v in sig_lags:
        print(f"  {k}: RR={v['rate_ratio']:.3f} p={v['p_value']:.4f}")


if __name__ == "__main__":
    main()
