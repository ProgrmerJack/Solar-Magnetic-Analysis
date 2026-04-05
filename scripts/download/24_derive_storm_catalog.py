"""
Derive geomagnetic storm catalog from Kyoto Dst CSV.
Applies standard thresholds: moderate (Dst < -50 nT) and intense (Dst < -100 nT).
Also parses OMNI2 proton flux (>10 MeV) for Solar Energetic Particle (SEP) events.
Outputs:
  data/geomagnetic/dst_ae_index/storm_catalog.csv
  data/geomagnetic/dst_ae_index/sep_events.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from download.utils import get_logger

import pandas as pd
import numpy as np

logger = get_logger("storm_catalog")
DATA = Path(__file__).resolve().parents[2] / "data"
GEO_DIR = DATA / "geomagnetic" / "dst_ae_index"
OMNI_DIR = DATA / "atmospheric" / "omni_solar_wind" / "low_res"

# -------------------------------------------------------------------
# 1. Load Kyoto Dst
# -------------------------------------------------------------------
def load_kyoto_dst():
    dst_path = GEO_DIR / "kyoto_dst_1957_present.csv"
    if not dst_path.exists():
        logger.error(f"Kyoto Dst not found at {dst_path}")
        return None
    logger.info(f"Loading Kyoto Dst from {dst_path.name}")
    df = pd.read_csv(dst_path, parse_dates=["datetime_utc"], low_memory=False)
    df = df.rename(columns={"datetime_utc": "datetime", "dst_nT": "dst"})
    df["dst"] = pd.to_numeric(df["dst"], errors="coerce")
    df = df[df["dst"].notna() & (df["dst"] != 9999)].copy()
    df.sort_values("datetime", inplace=True)
    df.reset_index(drop=True, inplace=True)
    logger.info(f"  {len(df):,} valid hourly Dst records ({df['datetime'].min().year}–{df['datetime'].max().year})")
    return df

# -------------------------------------------------------------------
# 2. Storm identification (classic minimum-method)
# -------------------------------------------------------------------
def identify_storms(df: pd.DataFrame, moderate_thresh=-50, intense_thresh=-100) -> pd.DataFrame:
    """
    Identify geomagnetic storm events using a simple storm-epoch method:
    - Sudden commencement or gradual onset threshold: Dst drops below -30 nT
    - Storm main phase defined by Dst < threshold (moderate: -50, intense: -100)
    - Recovery phase: Dst returns above -20 nT
    Returns catalog with onset, minimum, recovery datetimes, and peak Dst.
    """
    logger.info("Identifying geomagnetic storms ...")
    storms = []
    in_storm = False
    storm_start = None
    peak_dst = 0.0
    peak_time = None

    for _, row in df.iterrows():
        val = row["dst"]
        t = row["datetime"]

        if not in_storm:
            if val < moderate_thresh:
                in_storm = True
                storm_start = t
                peak_dst = val
                peak_time = t
        else:
            if val < peak_dst:
                peak_dst = val
                peak_time = t
            # Recovery: Dst rises back above -20 nT for at least 1 hour
            if val > -20:
                storms.append({
                    "onset": storm_start,
                    "peak_time": peak_time,
                    "peak_dst_nT": peak_dst,
                    "recovery_start": t,
                    "duration_h": (t - storm_start).total_seconds() / 3600,
                    "class": "intense" if peak_dst < intense_thresh else "moderate"
                })
                in_storm = False
                storm_start = None
                peak_dst = 0.0

    catalog = pd.DataFrame(storms)
    logger.info(f"  Found {len(catalog)} storms total")
    logger.info(f"    Moderate (Dst < {moderate_thresh} nT): {(catalog['class']=='moderate').sum()}")
    logger.info(f"    Intense  (Dst < {intense_thresh} nT): {(catalog['class']=='intense').sum()}")
    return catalog

# -------------------------------------------------------------------
# 3. Load OMNI2 low-res and derive SEP events (proton flux >10 MeV)
# -------------------------------------------------------------------
OMNI2_COLS = {
    0: "year", 1: "doy", 2: "hour",
    38: "kp", 40: "dst_omni", 41: "ae",
    42: "proton_flux_gt1MeV", 43: "proton_flux_gt2MeV",
    44: "proton_flux_gt4MeV", 45: "proton_flux_gt10MeV",
    46: "proton_flux_gt30MeV", 47: "proton_flux_gt60MeV",
    49: "f107",
}

def load_omni2_yearly(year: int) -> pd.DataFrame | None:
    fpath = OMNI_DIR / f"omni2_{year}.dat"
    if not fpath.exists():
        return None
    try:
        df = pd.read_csv(fpath, sep=r"\s+", header=None,
                         usecols=list(OMNI2_COLS.keys()),
                         names=[f"c{i}" for i in range(55)])
        df.columns = list(OMNI2_COLS.values())
        # Build datetime
        df["datetime"] = pd.to_datetime(
            df["year"].astype(str) + df["doy"].astype(str).str.zfill(3) + df["hour"].astype(str).str.zfill(2),
            format="%Y%j%H", errors="coerce"
        )
        # Replace fill values
        for col in ["proton_flux_gt10MeV", "proton_flux_gt30MeV", "ae"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
                df.loc[df[col] > 9e5, col] = np.nan  # fill=999999.99
        return df[["datetime", "kp", "dst_omni", "ae",
                   "proton_flux_gt1MeV", "proton_flux_gt10MeV",
                   "proton_flux_gt30MeV", "f107"]]
    except Exception as e:
        logger.warning(f"  OMNI2 {year}: {e}")
        return None

def derive_sep_events(omni_df: pd.DataFrame, threshold_pfu=10.0) -> pd.DataFrame:
    """
    SEP events: sustained proton flux >10 MeV exceeding threshold (default 10 pfu).
    Standard NOAA threshold for S1 radiation storm is 10 pfu (>10 MeV).
    """
    logger.info(f"Deriving SEP events (>10 MeV proton flux > {threshold_pfu} pfu) ...")
    flux = omni_df[["datetime", "proton_flux_gt10MeV"]].dropna()
    above = flux[flux["proton_flux_gt10MeV"] >= threshold_pfu].copy()

    # Cluster events separated by >24h
    events = []
    if len(above) == 0:
        return pd.DataFrame(events)

    above = above.sort_values("datetime").reset_index(drop=True)
    start = above.loc[0, "datetime"]
    peak_flux = above.loc[0, "proton_flux_gt10MeV"]
    peak_t = above.loc[0, "datetime"]
    prev = above.loc[0, "datetime"]

    for i in range(1, len(above)):
        t = above.loc[i, "datetime"]
        fl = above.loc[i, "proton_flux_gt10MeV"]
        if (t - prev).total_seconds() > 86400:  # gap > 24h = new event
            events.append({"onset": start, "peak_time": peak_t, "peak_flux_pfu": peak_flux,
                           "end": prev, "duration_h": (prev - start).total_seconds()/3600})
            start = t
            peak_flux = fl
            peak_t = t
        else:
            if fl > peak_flux:
                peak_flux = fl
                peak_t = t
        prev = t

    events.append({"onset": start, "peak_time": peak_t, "peak_flux_pfu": peak_flux,
                   "end": prev, "duration_h": (prev - start).total_seconds()/3600})

    sep_df = pd.DataFrame(events)
    logger.info(f"  {len(sep_df)} SEP events (>10 MeV, ≥{threshold_pfu} pfu)")

    # NOAA storm classes
    sep_df["noaa_class"] = "S1"
    sep_df.loc[sep_df["peak_flux_pfu"] >= 100, "noaa_class"] = "S2"
    sep_df.loc[sep_df["peak_flux_pfu"] >= 1000, "noaa_class"] = "S3"
    sep_df.loc[sep_df["peak_flux_pfu"] >= 10000, "noaa_class"] = "S4"
    sep_df.loc[sep_df["peak_flux_pfu"] >= 100000, "noaa_class"] = "S5"
    return sep_df

# -------------------------------------------------------------------
# main
# -------------------------------------------------------------------
def main():
    logger.info("=== Deriving Geomagnetic Storm + SEP Event Catalogs ===")

    # Storm catalog from Kyoto Dst
    dst_df = load_kyoto_dst()
    if dst_df is not None:
        storms = identify_storms(dst_df)
        out = GEO_DIR / "storm_catalog.csv"
        storms.to_csv(out, index=False)
        logger.info(f"✓  Storm catalog → {out.name} ({len(storms)} events)")

    # Load OMNI2 (all available years)
    logger.info("\nLoading OMNI2 low-res data for SEP analysis ...")
    frames = []
    for yr in range(1963, 2026):
        f = load_omni2_yearly(yr)
        if f is not None:
            frames.append(f)
    if frames:
        omni = pd.concat(frames, ignore_index=True)
        omni.sort_values("datetime", inplace=True)
        logger.info(f"  OMNI2 loaded: {len(omni):,} hourly records ({omni['datetime'].min().year}–{omni['datetime'].max().year})")

        sep = derive_sep_events(omni)
        if len(sep) > 0:
            out2 = GEO_DIR / "sep_events.csv"
            sep.to_csv(out2, index=False)
            logger.info(f"✓  SEP catalog → {out2.name}")

        # Also save AE+Kp+F10.7 daily summary from OMNI2
        omni["date"] = omni["datetime"].dt.date
        daily = omni.groupby("date").agg(
            kp_max=("kp", "max"),
            ae_max=("ae", "max"),
            ae_mean=("ae", "mean"),
            dst_min=("dst_omni", "min"),
            proton_flux_gt10_max=("proton_flux_gt10MeV", "max"),
            f107=("f107", "first")
        ).reset_index()
        out3 = GEO_DIR / "omni2_daily_summary.csv"
        daily.to_csv(out3, index=False)
        logger.info(f"✓  OMNI2 daily summary → {out3.name} ({len(daily):,} days)")

    logger.info("\n=== Derived catalogs complete ===")

if __name__ == "__main__":
    main()
