"""
12_rebuild_enhanced_panel.py — Build comprehensive panel with ALL datasets
=========================================================================
Extends the original 76-column panel with:
  - F10.7 solar radio flux (solar cycle proxy)
  - MODIS Alps snow fraction
  - IMS Northern Hemisphere snow
  - MLS HNO3 + N2O (atmospheric chemistry bridge)
  - POES EPP daily aggregate (energetic particle precipitation)
  - CME event indicators
  - SNOTEL aggregated SWE + precipitation
  - Norway avalanche daily counts (independent validation)
  - GOES XRS daily peak flux

Output: data/processed/analysis_panel_v2.parquet
"""
import sys, gc, logging
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, LOG, save_analysis

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def tz_strip(idx):
    """Remove timezone if present."""
    if hasattr(idx, 'tz') and idx.tz:
        return idx.tz_localize(None)
    return idx


def main():
    # ── Load existing panel as base ──────────────────────────────────────
    LOG.info("Loading base panel...")
    panel = pd.read_parquet(PROCESSED / "analysis_panel.parquet")
    panel.index = tz_strip(panel.index)
    LOG.info("Base panel: %d rows x %d cols", len(panel), len(panel.columns))

    # ── 1. F10.7 Solar Radio Flux ────────────────────────────────────────
    LOG.info("Adding F10.7...")
    si = pd.read_parquet(PROCESSED / "solar" / "solar_indices.parquet")
    si.index = tz_strip(si.index)
    si = si.rename(columns={"ngdc_f107_adjusted_daily": "f107"})
    # Solar cycle phase: high = above median, low = below
    f107_median = si["f107"].median()
    si["f107_high"] = (si["f107"] > f107_median).astype(int)
    # Smoothed F10.7 (81-day running mean, standard solar physics smoothing)
    si["f107_smooth"] = si["f107"].rolling(81, center=True, min_periods=40).mean()
    panel = panel.join(si[["f107", "f107_high", "f107_smooth"]], how="left")
    panel["f107"] = panel["f107"].ffill(limit=3)
    panel["f107_high"] = panel["f107_high"].ffill(limit=3)
    panel["f107_smooth"] = panel["f107_smooth"].ffill(limit=3)
    LOG.info("  F10.7: %d non-null", panel["f107"].notna().sum())
    del si; gc.collect()

    # ── 2. MODIS Alps Snow Fraction ──────────────────────────────────────
    LOG.info("Adding MODIS snow...")
    modis = pd.read_parquet(PROCESSED / "cryosphere" / "modis_alps_stats.parquet")
    modis.index = tz_strip(modis.index)
    modis = modis.rename(columns={
        "snow_fraction": "modis_snow_frac",
        "cloud_fraction": "modis_cloud_frac",
    })
    panel = panel.join(modis[["modis_snow_frac", "modis_cloud_frac"]], how="left")
    LOG.info("  MODIS: %d non-null", panel["modis_snow_frac"].notna().sum())
    del modis; gc.collect()

    # ── 3. IMS Snow ──────────────────────────────────────────────────────
    LOG.info("Adding IMS snow...")
    ims = pd.read_parquet(PROCESSED / "cryosphere" / "ims_snow_daily.parquet")
    ims.index = tz_strip(ims.index)
    ims = ims.rename(columns={
        "nh_snow_fraction": "ims_nh_snow",
        "alps_snow_fraction": "ims_alps_snow",
    })
    panel = panel.join(ims[["ims_nh_snow", "ims_alps_snow"]], how="left")
    LOG.info("  IMS: %d non-null", panel["ims_nh_snow"].notna().sum())
    del ims; gc.collect()

    # ── 4. MLS HNO3 + N2O ───────────────────────────────────────────────
    LOG.info("Adding MLS HNO3...")
    hno3 = pd.read_parquet(PROCESSED / "atmospheric" / "mls_hno3_polar.parquet")
    hno3.index = tz_strip(hno3.index)
    hno3.columns = ["mls_hno3_" + c for c in hno3.columns]
    panel = panel.join(hno3, how="left")
    LOG.info("  MLS HNO3: %d non-null", panel["mls_hno3_lev_6p8hpa"].notna().sum())
    del hno3; gc.collect()

    LOG.info("Adding MLS N2O...")
    n2o = pd.read_parquet(PROCESSED / "atmospheric" / "mls_n2o_polar.parquet")
    n2o.index = tz_strip(n2o.index)
    n2o.columns = ["mls_n2o_" + c for c in n2o.columns]
    panel = panel.join(n2o, how="left")
    LOG.info("  MLS N2O: %d non-null", panel["mls_n2o_lev_6p8hpa"].notna().sum())
    del n2o; gc.collect()

    # ── 5. CME Event Indicators ──────────────────────────────────────────
    LOG.info("Adding CME indicators...")
    cme = pd.read_parquet(PROCESSED / "solar" / "cme_catalog.parquet")
    cme.index = tz_strip(cme.index)
    cme_daily = cme.resample("D").agg(
        cme_count=("speed_km_s", "count"),
        cme_max_speed=("speed_km_s", "max"),
    )
    cme_daily["cme_fast"] = (cme_daily["cme_max_speed"] > 500).astype(int)
    panel = panel.join(cme_daily, how="left")
    panel["cme_count"] = panel["cme_count"].fillna(0).astype(int)
    panel["cme_fast"] = panel["cme_fast"].fillna(0).astype(int)
    LOG.info("  CME: %d days with CMEs", (panel["cme_count"] > 0).sum())
    del cme, cme_daily; gc.collect()

    # ── 6. GOES XRS Daily Peak ───────────────────────────────────────────
    LOG.info("Adding GOES XRS daily peak...")
    xrs = pd.read_parquet(PROCESSED / "solar" / "goes_xrs.parquet")
    xrs.index = tz_strip(xrs.index)
    xrs_daily = xrs[["xrsb_flux"]].resample("D").max()
    xrs_daily.columns = ["xrs_peak_flux"]
    # Log-transform (standard in solar physics)
    xrs_daily["xrs_log_flux"] = np.log10(xrs_daily["xrs_peak_flux"].clip(lower=1e-10))
    panel = panel.join(xrs_daily, how="left")
    LOG.info("  XRS: %d non-null", panel["xrs_peak_flux"].notna().sum())
    del xrs, xrs_daily; gc.collect()

    # ── 7. POES EPP Daily Aggregate ──────────────────────────────────────
    LOG.info("Adding POES EPP (aggregating satellite passes to daily)...")
    epp_frames = []
    atm_dir = PROCESSED / "atmospheric"
    poes_files = sorted(atm_dir.glob("poes_*.parquet"))
    for pf in poes_files:
        LOG.info("  Loading %s...", pf.name)
        try:
            p = pd.read_parquet(pf, columns=[
                "ted_ele_tel0_hi_eflux", "ted_ele_tel30_hi_eflux",
                "ted_pro_tel0_hi_eflux", "ted_pro_tel30_hi_eflux",
            ])
            p.index = tz_strip(p.index)
            # Total high-energy particle flux (electrons + protons, both telescopes)
            p["total_epp_flux"] = (
                p["ted_ele_tel0_hi_eflux"].fillna(0) +
                p["ted_ele_tel30_hi_eflux"].fillna(0) +
                p["ted_pro_tel0_hi_eflux"].fillna(0) +
                p["ted_pro_tel30_hi_eflux"].fillna(0)
            )
            daily = p[["total_epp_flux"]].resample("D").agg(["mean", "max"])
            daily.columns = ["epp_mean", "epp_max"]
            epp_frames.append(daily)
            del p; gc.collect()
        except Exception as e:
            LOG.warning("  Skipped %s: %s", pf.name, e)

    if epp_frames:
        epp_all = pd.concat(epp_frames)
        # Average across satellites for same day
        epp_daily = epp_all.groupby(epp_all.index).mean()
        epp_daily["epp_log_mean"] = np.log10(epp_daily["epp_mean"].clip(lower=1e-5))
        panel = panel.join(epp_daily, how="left")
        LOG.info("  EPP: %d non-null days", panel["epp_mean"].notna().sum())
        del epp_frames, epp_all, epp_daily; gc.collect()

    # ── 8. SNOTEL Aggregated ─────────────────────────────────────────────
    LOG.info("Adding SNOTEL aggregated...")
    # SNOTEL is 10M+ rows — read only what we need, aggregate by date
    snotel = pd.read_parquet(PROCESSED / "cryosphere" / "snotel_daily.parquet",
                              columns=["wteq_mm", "prec_mm", "tavg_c"])
    snotel.index = tz_strip(snotel.index)
    snotel_daily = snotel.groupby(snotel.index.date).agg(
        snotel_swe_mean=("wteq_mm", "mean"),
        snotel_prec_mean=("prec_mm", "mean"),
        snotel_temp_mean=("tavg_c", "mean"),
    )
    snotel_daily.index = pd.DatetimeIndex(snotel_daily.index, name="date")
    panel = panel.join(snotel_daily, how="left")
    LOG.info("  SNOTEL: %d non-null", panel["snotel_swe_mean"].notna().sum())
    del snotel, snotel_daily; gc.collect()

    # ── 9. Norway Avalanche Daily Counts ─────────────────────────────────
    LOG.info("Adding Norway avalanche counts...")
    nor = pd.read_parquet(PROCESSED / "cryosphere" / "norway_avalanche.parquet")
    nor.index = tz_strip(nor.index)
    # Count events per day
    nor_daily = nor.groupby(nor.index.date).size().rename("norway_aval_count")
    nor_daily.index = pd.DatetimeIndex(nor_daily.index, name="date")
    panel = panel.join(nor_daily.to_frame(), how="left")
    panel["norway_aval_count"] = panel["norway_aval_count"].fillna(0).astype(int)
    LOG.info("  Norway: %d days with avalanches", (panel["norway_aval_count"] > 0).sum())
    del nor, nor_daily; gc.collect()

    # ── 10. Reconstruct event windows for new data ───────────────────────
    # Post-event exposure for additional lag windows used in analysis
    event_dates = panel.index[panel["geo_event"] == 1]
    panel["post_event_3_8d"] = 0   # CME propagation lag
    panel["post_event_15_30d"] = 0  # SSW lag
    panel["post_event_30_60d"] = 0  # Long-term rebound

    for ed in event_dates:
        m1 = (panel.index > ed + pd.Timedelta(days=2)) & \
             (panel.index <= ed + pd.Timedelta(days=8))
        m2 = (panel.index >= ed + pd.Timedelta(days=15)) & \
             (panel.index <= ed + pd.Timedelta(days=30))
        m3 = (panel.index >= ed + pd.Timedelta(days=30)) & \
             (panel.index <= ed + pd.Timedelta(days=60))
        panel.loc[m1, "post_event_3_8d"] = 1
        panel.loc[m2, "post_event_15_30d"] = 1
        panel.loc[m3, "post_event_30_60d"] = 1

    # ── Summary ──────────────────────────────────────────────────────────
    LOG.info("Enhanced panel: %d rows x %d cols", len(panel), len(panel.columns))

    # Save
    out = PROCESSED / "analysis_panel_v2.parquet"
    save_analysis(panel, out, {
        "title": "Enhanced daily analysis panel with all datasets",
        "created_by": "12_rebuild_enhanced_panel.py",
        "n_data_sources": "21+",
    })

    # Print summary
    winter = panel[panel["is_winter"] == 1]
    new_cols = [c for c in panel.columns if c not in [
        'kp_max', 'ap_max', 'ap_daily', 'dst_min', 'month', 'day_of_year',
        'winter_id', 'day_of_season', 'day_of_season_sq', 'is_winter', 'is_summer',
        'geo_event_raw', 'geo_event', 'post_event_1_3d', 'post_event_5_21d', 'post_event_0_30d',
    ]]
    print("\n" + "=" * 60)
    print("ENHANCED PANEL SUMMARY")
    print("=" * 60)
    print("Total rows: %d, Total cols: %d" % (len(panel), len(panel.columns)))
    print("Winter days: %d" % len(winter))
    print("\nNew columns added:")
    for c in sorted(panel.columns):
        if c.startswith(("f107", "modis", "ims", "mls_hno3", "mls_n2o",
                         "cme", "xrs", "epp", "snotel", "norway", "post_event_3",
                         "post_event_15", "post_event_30")):
            nn = panel[c].notna().sum()
            nn_w = winter[c].notna().sum() if c in winter.columns else 0
            print("  %s: %d total (%d winter)" % (c, nn, nn_w))

    print("\nData coverage in winter NDJFM:")
    for c in ["f107", "modis_snow_frac", "ims_nh_snow", "mls_hno3_lev_6p8hpa",
              "mls_n2o_lev_6p8hpa", "epp_mean", "snotel_swe_mean", "xrs_peak_flux"]:
        if c in winter.columns:
            nn = winter[c].notna().sum()
            pct = 100 * nn / len(winter)
            print("  %s: %d / %d (%.1f%%)" % (c, nn, len(winter), pct))


if __name__ == "__main__":
    main()
