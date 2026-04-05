"""
Download and parse the Bremen Mg II composite (V8) from IUP Bremen.
This is the canonical solar chromospheric UV proxy used in EPP-atmosphere studies.
Updated daily: 1978-present.
URL: https://www.iup.uni-bremen.de/gome/solar/MgII_composite.dat
"""
import requests
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger, DATA_DIR

LOG = get_logger("mgii")
OUT = DATA_DIR / "solar" / "solar_indices"
OUT.mkdir(parents=True, exist_ok=True)

RAW = OUT / "bremen_mgii_composite_v8.dat"
CSV = OUT / "bremen_mgii_composite_v8.csv"

URL = "https://www.iup.uni-bremen.de/gome/solar/MgII_composite.dat"


def download_raw():
    LOG.info("Downloading Bremen Mg II V8 composite from IUP Bremen ...")
    r = requests.get(URL, timeout=60)
    r.raise_for_status()
    RAW.write_text(r.text, encoding="utf-8")
    LOG.info("  Downloaded %d B → %s", len(r.text), RAW)
    return r.text


def parse_and_save(text=None):
    if text is None:
        text = RAW.read_text(encoding="utf-8")

    # Find source descriptions from comments
    source_map = {}
    rows = []
    for line in text.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(";"):
            # Try to parse source legend: ";  N  description  count"
            import re
            m = re.match(r";\s+(\d+)\s+(.+?)\s+(\d+)\s*$", line)
            if m:
                source_map[int(m.group(1))] = m.group(2).strip()
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            decimal_year = float(parts[0])
            month = int(parts[1])
            day = int(parts[2])
            mgii = float(parts[3])
            uncertainty = float(parts[4])
            source_id = int(parts[5])
            # Derive year from decimal year
            year = int(decimal_year)
            rows.append({
                "decimal_year": decimal_year,
                "year": year,
                "month": month,
                "day": day,
                "date": f"{year:04d}-{month:02d}-{day:02d}",
                "mgii_index": mgii,
                "uncertainty": uncertainty,
                "source_id": source_id,
                "source": source_map.get(source_id, f"source_{source_id}"),
            })
        except (ValueError, IndexError):
            continue

    df = pd.DataFrame(rows)
    df = df.sort_values("date").drop_duplicates("date")
    df.to_csv(CSV, index=False)
    LOG.info("Bremen Mg II V8: %d records (%s to %s) → %s",
             len(df), df['date'].min(), df['date'].max(), CSV)
    LOG.info("Sources represented: %s", dict(df['source_id'].value_counts().head()))
    return df


if __name__ == "__main__":
    text = download_raw()
    parse_and_save(text)
