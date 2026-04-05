"""
Download solar spectral irradiance proxies and activity indices from LASP LISIRD.

Datasets:
  - Bremen Mg II composite (1978-present) — chromospheric UV proxy, best EPP indicator
  - Composite Lyman-alpha (1947-present) — 121.6 nm H-Lyman-alpha, ionospheric driver
  - NRL2 SSI P1D — Solar Spectral Irradiance 115-100,000 nm daily (1882-present)
  - Penticton/DRAO F10.7 daily (updated) — 10.7cm radio flux
  - TIMED/SEE EUV solar spectrum level 3 — if available
  - OMI solar indices
  - NRL2 preliminary TSI daily (most recent months)
"""
import requests
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger, DATA_DIR

LOG = get_logger("solar_proxies")
OUT = DATA_DIR / "solar" / "solar_indices"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://lasp.colorado.edu/lisird/latis/dap"


def fetch_csv(dataset_id, filename, desc, parse_dt_col=None):
    """Download a LASP LISIRD CSV dataset."""
    url = f"{BASE}/{dataset_id}.csv"
    out = OUT / filename
    if out.exists():
        LOG.info("  %s already exists, skipping", filename)
        return pd.read_csv(out)
    LOG.info("Downloading %s ...", desc)
    try:
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        out.write_text(r.text, encoding="utf-8")
        lines = [l for l in r.text.split("\n") if l.strip()]
        LOG.info("  %s: %d rows -> %s", desc, len(lines) - 1, out)
        return pd.read_csv(out)
    except Exception as exc:
        LOG.warning("  %s FAILED: %s", desc, exc)
        return None


def download_mgii():
    """Bremen Mg II composite — best solar UV/chromospheric activity proxy."""
    # mg II index: solar proxy for EUV/FUV irradiance, key for EPP ionization
    df = fetch_csv("bremen_composite_mgii_v4", "mgii_composite_daily.csv", "Bremen Mg II composite v4")
    if df is None:
        df = fetch_csv("bremen_composite_mgii", "mgii_composite_daily.csv", "Bremen Mg II composite")
    if df is not None:
        LOG.info("  Mg II columns: %s", list(df.columns)[:6])
        LOG.info("  Mg II range: %s to %s", df.iloc[1, 0] if len(df) > 1 else "?", df.iloc[-1, 0])


def download_lyman_alpha():
    """Composite Lyman-alpha 121.6 nm — primary thermosphere ionization source."""
    df = fetch_csv("composite_lyman_alpha_v3", "lyman_alpha_daily.csv", "Composite Lyman-alpha v3")
    if df is None:
        df = fetch_csv("composite_lyman_alpha", "lyman_alpha_daily.csv", "Composite Lyman-alpha")
    if df is not None:
        LOG.info("  Lyman-alpha columns: %s", list(df.columns)[:5])


def download_f107():
    """Penticton/DRAO adjusted F10.7 daily — same scale as our existing LASP record, extended."""
    # Try adjusted Penticton record (standard for solar activity)
    df = fetch_csv("penticton_radio_flux_nearest_noon", "f107_penticton_daily.csv", "Penticton F10.7 nearest noon")
    if df is None:
        df = fetch_csv("noaa_radio_flux", "f107_noaa_daily.csv", "NOAA F10.7 radio flux")
    if df is not None:
        LOG.info("  F10.7 columns: %s", list(df.columns)[:5])
        LOG.info("  F10.7 range: %s rows", len(df))


def download_ssi_daily():
    """NRL2 Solar Spectral Irradiance — full solar spectrum 1882-present."""
    # This is large — try preliminary (most recent) first, then skip full if large
    LOG.info("Checking NRL2 SSI daily size ...")
    try:
        # Check headers only for size
        r = requests.head(f"{BASE}/nrl2_ssi_prelim_P1D.csv", timeout=10)
        LOG.info("  NRL2 SSI prelim HEAD: %s, size: %s", r.status_code, r.headers.get("Content-Length", "unknown"))
    except Exception as e:
        LOG.warning("  SSI prelim head failed: %s", e)

    # Download preliminary (most recent solar cycle data)
    df = fetch_csv("nrl2_ssi_prelim_P1D", "nrl2_ssi_prelim_daily.csv", "NRL2 SSI prelim daily")
    if df is not None:
        LOG.info("  NRL2 SSI prelim columns: %s ... (%d total)", list(df.columns)[:4], len(df.columns))


def download_timed_see():
    """TIMED/SEE solar EUV level 3 — direct EUV flux driving thermospheric chemistry."""
    df = fetch_csv("timed_see_ssi_l3", "timed_see_euv_daily.csv", "TIMED/SEE EUV L3 daily")
    if df is not None:
        LOG.info("  TIMED/SEE columns: %s ... (%d total)", list(df.columns)[:4], len(df.columns))


def download_omi_solar():
    """OMI solar indices from Aura/OMI."""
    df = fetch_csv("omi_solar_indices", "omi_solar_indices_daily.csv", "OMI solar indices")
    if df is not None:
        LOG.info("  OMI columns: %s", list(df.columns)[:6])


def download_historical_tsi():
    """Historical TSI composite (1610-present) — may be PMOD/WRC or similar."""
    df = fetch_csv("historical_tsi", "historical_tsi_1610_present.csv", "Historical TSI composite")
    if df is None:
        df = fetch_csv("nrl2_historical_tsi", "nrl2_historical_tsi.csv", "NRL2 historical TSI")
    if df is not None:
        LOG.info("  Historical TSI: %d records", len(df))


def download_nrl2_tsi_prelim():
    """NRL2 TSI preliminary (most recent months, not in main record)."""
    df = fetch_csv("nrl2_tsi_prelim_P1D", "nrl2_tsi_prelim_daily.csv", "NRL2 TSI prelim daily")
    if df is not None:
        LOG.info("  NRL2 TSI prelim: %d records", len(df))


if __name__ == "__main__":
    LOG.info("=== Solar Spectral Proxies & Indices ===")
    download_mgii()
    download_lyman_alpha()
    download_f107()
    download_timed_see()
    download_omi_solar()
    download_historical_tsi()
    download_nrl2_tsi_prelim()
    # Skip SSI daily (potentially very large - many wavelength bins)
    LOG.info("=== Done ===")
