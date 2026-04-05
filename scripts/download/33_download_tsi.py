"""
Download Total Solar Irradiance (TSI) composites from LASP LISIRD.

Datasets:
  - NRL2 TSI daily (~1882-present) — the standard TSI composite used in climate/solar research
  - SORCE/TIM daily 2003-2020 — high-precision instrument-grade TSI
  - TSIS-1/TIM daily 2018-present — current NIST-traceable TSI standard
  - NRL2 annual 1610-present — for long multi-century context
"""
import csv
import requests
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger, DATA_DIR

LOG = get_logger("tsi")
OUT = DATA_DIR / "solar" / "solar_indices"
OUT.mkdir(parents=True, exist_ok=True)


def julian_to_dt(jd):
    """Convert Julian Date to datetime."""
    from datetime import datetime, timedelta
    # JD 2440587.5 = 1970-01-01 00:00 UTC
    epoch = datetime(1970, 1, 1)
    return epoch + timedelta(days=jd - 2440587.5)


def days_since_to_dt(days, ref_year=1610):
    """Convert days since Jan 1, ref_year to datetime."""
    from datetime import date, timedelta
    base = date(ref_year, 1, 1)
    return base + timedelta(days=float(days))


def download_and_parse_nrl2_daily():
    """NRL2 TSI daily composite from LASP."""
    url = "https://lasp.colorado.edu/lisird/latis/dap/nrl2_tsi_P1D.csv"
    LOG.info("Downloading NRL2 TSI daily...")
    r = requests.get(url, timeout=60)
    r.raise_for_status()

    lines = [l.strip() for l in r.text.split("\n") if l.strip()]
    header = lines[0]
    LOG.info("NRL2 daily header: %s", header)

    # Detect time format from header
    use_days_since = "days since 1610" in header.lower()

    rows = []
    for line in lines[1:]:
        parts = line.split(",")
        if len(parts) < 2:
            continue
        try:
            if use_days_since:
                dt = days_since_to_dt(float(parts[0]), ref_year=1610)
            else:
                dt = julian_to_dt(float(parts[0]))
            rows.append({
                "date": dt.isoformat() if hasattr(dt, "isoformat") else str(dt),
                "tsi_1au_Wm2": float(parts[1]) if parts[1].strip() else None,
                "uncertainty_Wm2": float(parts[2]) if len(parts) > 2 and parts[2].strip() else None,
            })
        except (ValueError, IndexError):
            continue

    df = pd.DataFrame(rows)
    out = OUT / "nrl2_tsi_daily.csv"
    df.to_csv(out, index=False)
    LOG.info("NRL2 TSI daily: %d records (%s to %s) -> %s", len(df), df['date'].min(), df['date'].max(), out)
    return df


def download_and_parse_nrl2_annual():
    """NRL2 TSI annual 1610-present from LASP."""
    # Already downloaded by _download_tsi.py but let's re-parse properly
    raw = OUT / "nrl2_tsi_annual.csv"
    if not raw.exists():
        url = "https://lasp.colorado.edu/lisird/latis/dap/nrl2_tsi_P1Y.csv"
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        raw.write_text(r.text)

    rows = []
    for line in raw.read_text().split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        if len(parts) < 2:
            continue
        try:
            # time in days since 1610-01-01
            dt = days_since_to_dt(float(parts[0]), ref_year=1610)
            rows.append({
                "year": dt.year,
                "tsi_1au_Wm2": float(parts[1]),
                "uncertainty_Wm2": float(parts[2]) if len(parts) > 2 else None,
            })
        except (ValueError, IndexError):
            continue

    df = pd.DataFrame(rows)
    # Drop header row if it snuck in
    df = df[pd.to_numeric(df['year'], errors='coerce').notna()]
    out = OUT / "nrl2_tsi_annual_parsed.csv"
    df.to_csv(out, index=False)
    LOG.info("NRL2 TSI annual: %d records (%s to %s) -> %s", len(df), df['year'].min(), df['year'].max(), out)
    return df


def download_and_parse_sorce():
    """SORCE/TIM daily TSI 2003-2020."""
    raw = OUT / "sorce_tsi_daily_2003_2020.csv"
    rows = []
    for line in raw.read_text().split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        if len(parts) < 2:
            continue
        try:
            dt = julian_to_dt(float(parts[0]))
            tsi = float(parts[1])
            unc = float(parts[2]) if len(parts) > 2 else None
            rows.append({"date": dt.date().isoformat(), "tsi_1au_Wm2": tsi, "uncertainty_Wm2": unc})
        except (ValueError, IndexError):
            continue

    df = pd.DataFrame(rows)
    out = OUT / "sorce_tsi_daily_parsed.csv"
    df.to_csv(out, index=False)
    LOG.info("SORCE TSI daily: %d records (%s to %s) -> %s", len(df), df['date'].min(), df['date'].max(), out)
    return df


def download_and_parse_tsis():
    """TSIS-1/TIM daily TSI 2018-present."""
    raw = OUT / "tsis1_tsi_daily_2018_present.csv"
    if not raw.exists():
        raw = OUT / "tsis1_tsi_daily_full.csv"
    rows = []
    for line in raw.read_text().split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        if len(parts) < 2:
            continue
        try:
            dt = julian_to_dt(float(parts[0]))
            tsi = float(parts[1])
            if tsi == 0.0:
                continue  # missing data
            unc = float(parts[2]) if len(parts) > 2 else None
            rows.append({"date": dt.date().isoformat(), "tsi_1au_Wm2": tsi, "uncertainty_Wm2": unc})
        except (ValueError, IndexError):
            continue

    df = pd.DataFrame(rows)
    out = OUT / "tsis1_tsi_daily_parsed.csv"
    df.to_csv(out, index=False)
    LOG.info("TSIS-1 TSI daily: %d records (%s to %s) -> %s", len(df), df['date'].min(), df['date'].max(), out)
    return df


def build_merged_tsi():
    """Merge SORCE + TSIS-1 into a continuous 2003-present daily record."""
    sorce = pd.read_csv(OUT / "sorce_tsi_daily_parsed.csv", parse_dates=["date"])
    tsis = pd.read_csv(OUT / "tsis1_tsi_daily_parsed.csv", parse_dates=["date"])

    # Use SORCE through 2020-02, then TSIS-1 after
    cutoff = pd.Timestamp("2020-02-25")
    sorce_part = sorce[sorce["date"] <= cutoff].copy()
    tsis_part = tsis[tsis["date"] > cutoff].copy()

    merged = pd.concat([sorce_part, tsis_part], ignore_index=True)
    merged = merged.sort_values("date").drop_duplicates("date")
    out = OUT / "tsi_daily_2003_present.csv"
    merged.to_csv(out, index=False)
    LOG.info("Merged TSI 2003-present: %d records -> %s", len(merged), out)
    return merged


if __name__ == "__main__":
    LOG.info("=== TSI Download & Parse ===")
    download_and_parse_nrl2_daily()
    download_and_parse_nrl2_annual()
    download_and_parse_sorce()
    download_and_parse_tsis()
    build_merged_tsi()
    LOG.info("TSI done.")
