"""
Download remaining climate teleconnection indices from NOAA PSL/ESRL.

PDO, AMO, MEI are critical confounds for long-period avalanche/snowpack analysis.
QBO 30 hPa (to complement our existing 50 hPa QBO).
"""
import requests
import re
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger, DATA_DIR

LOG = get_logger("climate_idx")
OUT = DATA_DIR / "atmospheric" / "climate_indices"
OUT.mkdir(parents=True, exist_ok=True)

# NOAA PSL teleconnection indices
INDICES = [
    ("https://psl.noaa.gov/data/correlation/pdo.data", "pdo_monthly.txt",
     "PDO (Pacific Decadal Oscillation) monthly"),
    ("https://psl.noaa.gov/data/correlation/amo.data", "amo_monthly.txt",
     "AMO (Atlantic Multidecadal Oscillation) monthly"),
    ("https://psl.noaa.gov/enso/mei/data/meiv2.data", "mei_v2_bimonthly.txt",
     "MEI v2 (Multivariate ENSO Index) bimonthly"),
    ("https://psl.noaa.gov/data/correlation/censo.data", "censo_monthly.txt",
     "CENSO (Cold/Warm Atlantic/Pacific pattern) monthly"),
    ("https://psl.noaa.gov/data/correlation/wp.data", "wp_monthly.txt",
     "WP (Western Pacific) index monthly"),
    ("https://psl.noaa.gov/data/correlation/ep.data", "ep_monthly.txt",
     "EP (East Pacific) index monthly"),
    # AO daily (CPC)
    ("https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii", "ao_monthly_cpc.txt",
     "AO monthly CPC"),
    ("https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii", "aao_monthly_cpc.txt",
     "AAO monthly CPC (Antarctic Oscillation)"),
    # QBO at various levels
    ("https://psl.noaa.gov/data/correlation/qbo.u30.data", "qbo_u30_monthly.txt",
     "QBO 30 hPa zonal wind monthly"),
    ("https://psl.noaa.gov/data/correlation/qbo.u10.data", "qbo_u10_monthly.txt",
     "QBO 10 hPa zonal wind monthly"),
    # Solar cycle proxy from NOAA
    ("https://psl.noaa.gov/data/correlation/solar.data", "noaa_solar_cycle_index.txt",
     "NOAA solar cycle index monthly"),
    # Polar vortex index (Holton-Tan)
    ("https://psl.noaa.gov/data/correlation/pvortex.data", "polar_vortex_index_monthly.txt",
     "Polar vortex index monthly"),
]

for url, fname, desc in INDICES:
    out = OUT / fname
    try:
        r = requests.get(url, timeout=30)
        LOG.info("%s %s: %d B", r.status_code, desc, len(r.text))
        if r.ok:
            lines = [l for l in r.text.split("\n") if l.strip()]
            LOG.info("  %d rows | first: %s", len(lines), lines[0][:60] if lines else "empty")
            if lines:
                LOG.info("  last: %s", lines[-1][:60])
            out.write_text(r.text)
        else:
            LOG.warning("  FAILED: %s", r.status_code)
    except Exception as exc:
        LOG.warning("  %s ERROR: %s", desc, exc)
