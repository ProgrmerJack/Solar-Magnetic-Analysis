"""
SDO/HMI SHARP Active Region Flux Summary via JSOC/drms.
No registration required — SHARP summary parameters are publicly accessible.

Downloads HMI Active Region Patch (HARP) flux statistics needed for SOC analysis:
  - USFLUX: Total unsigned magnetic flux per active region [Mx]
  - AREA_ACR: Active region area [microhemispheres]
  - TOTUSJH: Total unsigned current helicity
  - MEANGBZ: Mean vertical gradient of Bz
  - HARP number, NOAA AR number (when available), timestamps

These are the per-active-region scalar summaries, NOT the full magnetogram images.
The SHARP summary for all active regions 2010–2025 is ~200 MB total.

Reference: Bobra et al. 2014, Sol. Phys. doi:10.1007/s11207-014-0529-3
JSOC series: hmi.sharp_cea_720s_nrt or hmi.sharp_720s (definitive, 12-min cadence)
"""
import drms
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from utils import get_logger

LOG = get_logger("42_hmi_sharp")
OUT = Path(__file__).parents[2] / "data" / "solar" / "sdo_hmi_aia"
OUT.mkdir(parents=True, exist_ok=True)

# SHARP parameters needed for SOC analysis
# See Bobra et al. 2014 Table 1 for full descriptions
SHARP_PARAMS = [
    "T_REC",          # timestamp
    "HARPNUM",        # HARP number
    "NOAA_AR",        # NOAA active region number
    "LONMIN", "LONMAX", "LATMIN", "LATMAX",  # bounding box
    "USFLUX",         # Total unsigned flux [Mx] — PRIMARY SOC SIZE PROXY
    "MEANGBZ",        # Mean gradient of Bz
    "AREA_ACR",       # Active region area [microhemispheres]
    "TOTUSJH",        # Total unsigned current helicity
    "TOTUSJZ",        # Total unsigned vertical current
    "ABSNJZH",        # Absolute net current helicity
    "SAVNCPP",        # Sum of abs. values of net currents per polarity
    "MEANPOT",        # Mean photospheric magnetic free energy density
    "TOTPOT",         # Total photospheric magnetic free energy
    "MEANJZH",        # Mean current helicity Bz contribution
    "MEANALP",        # Mean twist parameter alpha
    "MEANSHR",        # Mean shear angle
    "SHRGT45",        # Fraction of area with shear > 45 degrees
    "SIZE",           # Area of patch [microhemispheres]
    "NACR",           # Number of pixels in active region mask
]

# Build parameter string for drms query (drop T_REC — already in index)
params_str = ", ".join(p for p in SHARP_PARAMS if p != "T_REC")

LOG.info("=== SDO/HMI SHARP Active Region Flux Download ===")
LOG.info("Connecting to JSOC (no auth required for summary data)...")

try:
    c = drms.Client()
    LOG.info("JSOC connection OK")

    # Solar-max peak years (2014, 2015) exceed single-query limit → split into halves
    SPLIT_YEARS = {2014, 2015}
    # Recent years: definitive series (hmi.sharp_720s) has ~6-month lag; use NRT for latest
    NRT_YEARS = {2024, 2025}

    for yr in range(2010, 2026):
        out_csv = OUT / f"hmi_sharp_flux_{yr}.csv"
        if out_csv.exists() and out_csv.stat().st_size > 5000:
            LOG.info("  skip %d (exists, %d KB)", yr, out_csv.stat().st_size // 1024)
            continue

        series = "hmi.sharp_cea_720s_nrt" if yr in NRT_YEARS else "hmi.sharp_720s"

        if yr in SPLIT_YEARS:
            # Split into H1 (Jan-Jun) and H2 (Jul-Dec) to avoid JSOC timeout
            halves = [
                (f"{yr}.01.01_TAI", f"{yr}.06.30_TAI"),
                (f"{yr}.07.01_TAI", f"{yr}.12.31_TAI"),
            ]
            parts = []
            for start, stop in halves:
                LOG.info("  querying SHARP %d %s-%s ...", yr, start[:7], stop[:7])
                query = f"{series}[1-9999][{start}-{stop}@12h]{{{params_str}}}"
                try:
                    keys = c.query(query, key=params_str, rec_index=True)
                    if keys is not None and len(keys) > 0:
                        parts.append(keys)
                        LOG.info("    → %d records", len(keys))
                except Exception as e:
                    LOG.warning("    half query failed: %s", e)
            if parts:
                import pandas as pd
                df = pd.concat(parts).drop_duplicates()
                df.to_csv(out_csv, index=True)
                LOG.info("  %d: %d records total → %s", yr, len(df), out_csv.name)
            continue

        if yr in NRT_YEARS:
            # NRT series: range+stride syntax fails; use /Nd duration per month
            import pandas as pd
            from datetime import date
            LOG.info("  querying SHARP %d (NRT, monthly chunks) ...", yr)
            parts = []
            for mo in range(1, 13):
                start = date(yr, mo, 1)
                end = date(yr+1, 1, 1) if mo == 12 else date(yr, mo+1, 1)
                if start > date.today():
                    break
                days = (end - start).days
                tai = start.strftime("%Y.%m.%d_TAI")
                ds = f"{series}[][{tai}/{days}d@12h]{{{params_str}}}"
                url = (f"http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_info?"
                       f"ds={requests.utils.quote(ds)}&op=rs_list&key={params_str}&max=200000")
                import requests as rq
                try:
                    resp = rq.get(url, timeout=120)
                    data = resp.json()
                    if "keywords" in data and data["keywords"]:
                        n = len(data["keywords"][0]["values"])
                        records = [{k["name"]: k["values"][i] for k in data["keywords"]} for i in range(n)]
                        parts.append(pd.DataFrame(records))
                        LOG.info("    %d-%02d: %d records", yr, mo, n)
                except Exception as e:
                    LOG.warning("    %d-%02d failed: %s", yr, mo, e)
            if parts:
                df = pd.concat(parts, ignore_index=True).drop_duplicates()
                df.to_csv(out_csv, index=False)
                LOG.info("  %d: %d records → %s", yr, len(df), out_csv.name)
            continue

    LOG.info("=== HMI SHARP download complete ===")
    LOG.info("Output: %s", OUT)

except Exception as e:
    LOG.error("JSOC connection failed: %s", e)
    LOG.info("Falling back to JSOC HTTP API...")

    import requests, json

    def jsoc_export_url(series, segments, start, stop, cadence="12h"):
        """JSOC export via HTTP API."""
        ds = f"{series}[1-9999][{start}-{stop}@{cadence}]{{{','.join(segments)}}}"
        url = (f"http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_info?"
               f"ds={requests.utils.quote(ds)}&op=rs_list&key={','.join(segments)}&max=200000")
        return url

    for yr in range(2010, 2026):
        out_csv = OUT / f"hmi_sharp_flux_{yr}.csv"
        if out_csv.exists() and out_csv.stat().st_size > 5000:
            continue
        LOG.info("  HTTP API %d ...", yr)
        url = jsoc_export_url(
            "hmi.sharp_720s",
            ["HARPNUM","NOAA_AR","USFLUX","AREA_ACR","TOTPOT","MEANGBZ"],
            f"{yr}.01.01_TAI",
            f"{yr}.12.31_TAI",
        )
        try:
            r = requests.get(url, timeout=300)
            data = r.json()
            if "keywords" in data:
                keys_data = data["keywords"]
                records = []
                if len(keys_data) > 0 and "values" in keys_data[0]:
                    n_rows = len(keys_data[0]["values"])
                    for i in range(n_rows):
                        row = {k["name"]: k["values"][i] for k in keys_data}
                        records.append(row)
                    df = pd.DataFrame(records)
                    df.to_csv(out_csv, index=False)
                    LOG.info("  %d: %d records", yr, len(df))
        except Exception as e2:
            LOG.warning("  HTTP API %d failed: %s", yr, e2)
