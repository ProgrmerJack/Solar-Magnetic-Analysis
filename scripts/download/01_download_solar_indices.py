"""
Script 01 — Solar Indices
Downloads:
  • SILSO International Sunspot Number (daily, monthly, yearly) from the Royal
    Observatory of Belgium (1818 / 1749 – present)
  • GFZ Potsdam definitive Kp / ap / Ap / F10.7 combined file (1932 – present)
  • LISIRD F10.7 daily solar flux (1947 – present)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from utils import DATA_DIR, get_logger, download_file

logger = get_logger("01_solar_indices")
OUT = DATA_DIR / "solar" / "solar_indices"
OUT.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------- #
# 1. SILSO Sunspot Number  (Royal Observatory of Belgium)                      #
# --------------------------------------------------------------------------- #
SILSO = "https://www.sidc.be/SILSO/DATA"
silso_files = {
    "SN_d_tot_V2.0.csv":  f"{SILSO}/SN_d_tot_V2.0.csv",   # daily total
    "SN_m_tot_V2.0.csv":  f"{SILSO}/SN_m_tot_V2.0.csv",   # monthly total
    "SN_y_tot_V2.0.csv":  f"{SILSO}/SN_y_tot_V2.0.csv",   # yearly total
    "SN_d_hem_V2.0.csv":  f"{SILSO}/SN_d_hem_V2.0.csv",   # daily hemispheric
    "SN_ms_tot_V2.0.csv": f"{SILSO}/SN_ms_tot_V2.0.csv",  # 13-month smoothed
}

logger.info("=== SILSO Sunspot Number ===")
for fname, url in silso_files.items():
    download_file(url, OUT / fname, desc=f"SILSO {fname}")

# --------------------------------------------------------------------------- #
# 2. GFZ Potsdam Kp + ap + Ap + F10.7 definitive file (1932 – present)        #
# --------------------------------------------------------------------------- #
KP_OUT = DATA_DIR / "geomagnetic" / "kp_index"
KP_OUT.mkdir(parents=True, exist_ok=True)
logger.info("=== GFZ Kp / ap / F10.7 ===")

kp_files = {
    "Kp_ap_Ap_SN_F107_since_1932.txt":
        "https://kp.gfz-potsdam.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt",
    "Kp_ap_nowcast.txt":
        "https://kp.gfz-potsdam.de/app/files/Kp_ap_nowcast.txt",
    "Kp_ap_forecast.txt":
        "https://kp.gfz-potsdam.de/app/files/Kp_ap_forecast.txt",
}
for fname, url in kp_files.items():
    download_file(url, KP_OUT / fname, desc=f"GFZ {fname}")

# --------------------------------------------------------------------------- #
# 3. NOAA NGDC / SWPC monthly Solar Region Summary and Geomagnetic indices    #
# --------------------------------------------------------------------------- #
logger.info("=== NOAA SWPC Geomagnetic Activity ===")
SWPC_BASE = "https://services.swpc.noaa.gov"
swpc_files = {
    "planetary_k_index_1m.json":    f"{SWPC_BASE}/json/planetary_k_index_1m.json",
    "geomag_dst_7day.json":         f"{SWPC_BASE}/products/noaa-estimated-planetary-k-index-1-minute.json",
    "solar_regions.json":           f"{SWPC_BASE}/json/solar_regions.json",
    "xrays_1d.json":                f"{SWPC_BASE}/json/goes/primary/xrays-1-day.json",
    "xrays_7d.json":                f"{SWPC_BASE}/json/goes/primary/xrays-7-day.json",
    "flares_24h.json":              f"{SWPC_BASE}/json/goes/primary/xray-flares-24-hours.json",
    "flares_7d.json":               f"{SWPC_BASE}/json/goes/primary/xray-flares-latest.json",
    "proton_flux_7d.json":          f"{SWPC_BASE}/json/goes/primary/integral-protons-7-day.json",
    "solar_wind_plasma_7d.json":    f"{SWPC_BASE}/products/solar-wind/plasma-7-day.json",
    "solar_wind_mag_7d.json":       f"{SWPC_BASE}/products/solar-wind/mag-7-day.json",
    "geospace_dst_7d.json":         f"{SWPC_BASE}/products/geospace/dst-7-day.json",
}
for fname, url in swpc_files.items():
    download_file(url, OUT / fname, desc=f"SWPC {fname}")

# --------------------------------------------------------------------------- #
# 4. NOAA Sunspot Cycle data (archived cycles summary)                        #
# --------------------------------------------------------------------------- #
logger.info("=== NOAA Solar Cycle Progression ===")
cycle_files = {
    "solar_cycle_progression.json": f"{SWPC_BASE}/json/solar-cycle/observed-solar-cycle-indices.json",
    "solar_cycle_predicted.json":   f"{SWPC_BASE}/json/solar-cycle/predicted-solar-cycle.json",
}
for fname, url in cycle_files.items():
    download_file(url, OUT / fname, desc=f"Solar cycle {fname}")

logger.info("=== Script 01 complete ===")
