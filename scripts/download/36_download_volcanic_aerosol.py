"""
Download volcanic eruption catalog and stratospheric aerosol forcing data.

Datasets:
  1. Smithsonian GVP eruption list via WFS API (VEI >= 4, 1900-present)
  2. NASA GISS stratospheric aerosol optical depth tau.line (1850-2012)
  3. Simple major eruption list for confound flagging
"""
import requests
import json
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger, DATA_DIR

LOG = get_logger("volcanic")
OUT = DATA_DIR / "atmospheric" / "volcanic_aerosol"
OUT.mkdir(parents=True, exist_ok=True)


def download_gvp_eruptions():
    """Download Smithsonian GVP eruption list (all Holocene eruptions with VEI)."""
    # WFS GetFeature for eruption list
    url = (
        "https://webservices.volcano.si.edu/geoserver/GVP-VOTW/ows"
        "?service=WFS&version=2.0.0&request=GetFeature"
        "&typeName=GVP-VOTW:Smithsonian_VOTW_Eruption_List"
        "&outputFormat=csv"
    )
    LOG.info("Downloading GVP eruption list (CSV)...")
    r = requests.get(url, timeout=60)
    LOG.info("  Status: %s, size: %d B", r.status_code, len(r.text))
    if r.ok:
        (OUT / "gvp_eruption_list_holocene.csv").write_text(r.text)
        df = pd.read_csv(pd.io.common.StringIO(r.text))
        LOG.info("  Columns: %s", list(df.columns)[:10])
        LOG.info("  Total eruptions: %d", len(df))
        # Filter for significant recent eruptions
        if "VEI" in df.columns:
            vei4 = df[pd.to_numeric(df["VEI"], errors="coerce") >= 4].copy()
            LOG.info("  VEI >= 4 eruptions: %d", len(vei4))
            (OUT / "gvp_eruptions_vei4plus.csv").write_text(vei4.to_csv(index=False))
        return df
    return None


def download_giss_tau():
    """Download NASA GISS stratospheric aerosol optical depth (1850-2012)."""
    url = "https://data.giss.nasa.gov/modelforce/strataer/tau.line_2012.12.txt"
    LOG.info("Downloading NASA GISS stratospheric AOD ...")
    r = requests.get(url, timeout=30)
    if r.ok:
        (OUT / "giss_strat_aod_tau_1850_2012.txt").write_text(r.text)
        lines = [l for l in r.text.split("\n") if l.strip() and not l.startswith(" ")]
        # Parse data lines: year.fraction, global, N.hemisphere, S.hemisphere
        data_lines = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    year_frac = float(parts[0])
                    if 1850 <= year_frac <= 2013:
                        data_lines.append({
                            "decimal_year": year_frac,
                            "year": int(year_frac),
                            "month": int(round((year_frac % 1) * 12)) + 1,
                            "aod_global": float(parts[1]),
                            "aod_north": float(parts[2]),
                            "aod_south": float(parts[3]),
                        })
                except ValueError:
                    continue

        df = pd.DataFrame(data_lines)
        (OUT / "giss_strat_aod_tau_1850_2012.csv").write_text(df.to_csv(index=False))
        LOG.info("  GISS AOD: %d records (%s to %s)", len(df),
                 df['decimal_year'].min(), df['decimal_year'].max())
        # Show major eruptions (global AOD > 0.05)
        major = df[df['aod_global'] > 0.05]
        if len(major) > 0:
            LOG.info("  Major eruption periods (global AOD > 0.05): %d months", len(major))
        return df
    LOG.warning("  GISS AOD download failed: %s", r.status_code)
    return None


def create_major_eruption_catalog():
    """
    Create a simplified catalog of major stratospheric eruptions for confound flagging.
    Based on published literature (Robock 2000, Toohey & Sigl 2017).
    """
    major_eruptions = [
        {"name": "Agung", "date": "1963-03-17", "VEI": 5, "notes": "~0.1 global AOD"},
        {"name": "Fuego", "date": "1974-10-14", "VEI": 4, "notes": "weak signal"},
        {"name": "St. Helens", "date": "1980-05-18", "VEI": 5, "notes": "mostly tropospheric"},
        {"name": "El Chichon", "date": "1982-04-04", "VEI": 5, "notes": "~0.1 global AOD; 2yr cooling"},
        {"name": "Nevado del Ruiz", "date": "1985-11-13", "VEI": 3, "notes": "minor"},
        {"name": "Kelut", "date": "1990-02-11", "VEI": 4, "notes": "minor strat signal"},
        {"name": "Pinatubo", "date": "1991-06-15", "VEI": 6, "notes": "~0.15 global AOD; largest 20thC"},
        {"name": "Hudson", "date": "1991-08-12", "VEI": 5, "notes": "S. hemisphere"},
        {"name": "Rabaul", "date": "1994-09-19", "VEI": 4, "notes": "minor"},
        {"name": "Ruapehu", "date": "1996-06-17", "VEI": 3, "notes": "minor"},
        {"name": "Soufriere Hills", "date": "1997-08-04", "VEI": 3, "notes": "minor"},
        {"name": "Shishaldin", "date": "1999-04-19", "VEI": 3, "notes": "minor"},
        {"name": "Tungurahua", "date": "2001-08-16", "VEI": 3, "notes": "minor"},
        {"name": "Reventador", "date": "2002-11-03", "VEI": 4, "notes": "small strat"},
        {"name": "Manam", "date": "2005-01-27", "VEI": 4, "notes": "moderate strat"},
        {"name": "Soufriere Hills2", "date": "2006-05-20", "VEI": 3, "notes": "minor"},
        {"name": "Sarychev", "date": "2009-06-12", "VEI": 4, "notes": "notable strat"},
        {"name": "Merapi", "date": "2010-10-26", "VEI": 4, "notes": "moderate"},
        {"name": "Nabro", "date": "2011-06-12", "VEI": 4, "notes": "notable SO2 stratospheric"},
        {"name": "Puyehue-Cordon Caulle", "date": "2011-06-04", "VEI": 5, "notes": "S hemisphere"},
        {"name": "Kelut 2014", "date": "2014-02-13", "VEI": 4, "notes": "moderate"},
        {"name": "Calbuco", "date": "2015-04-22", "VEI": 4, "notes": "S hemisphere, moderate"},
        {"name": "Sarychev2", "date": "2019-06-21", "VEI": 3, "notes": "minor"},
        {"name": "Taal", "date": "2020-01-12", "VEI": 4, "notes": "phreatomagmatic, minor strat"},
        {"name": "Hunga Tonga", "date": "2022-01-15", "VEI": 5,
         "notes": "exceptional water vapor injection; unusual strat warming"},
        {"name": "Shishaldin2", "date": "2023-10-25", "VEI": 3, "notes": "minor"},
    ]
    df = pd.DataFrame(major_eruptions)
    out = OUT / "major_eruptions_catalog.csv"
    df.to_csv(out, index=False)
    LOG.info("Major eruption catalog: %d events → %s", len(df), out)
    return df


if __name__ == "__main__":
    LOG.info("=== Volcanic Aerosol & Eruption Data ===")
    download_gvp_eruptions()
    download_giss_tau()
    create_major_eruption_catalog()
    LOG.info("=== Done ===")
