"""Master pipeline runner — processes all domains in order (largest → smallest)."""
import subprocess
import sys
import logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _utils import setup_logging, disk_free_gb, LOG, PROCESSED_ROOT

SCRIPTS = [
    "01_process_geomagnetic.py",
    "02_process_solar_catalogs.py",
    "03_process_omni.py",
    "04_process_goes_xrs.py",
    "05_process_goes_particle.py",
    "06_process_poes.py",
    "07_process_mls.py",
    "08_process_era5.py",
    "09_process_modis.py",   # largest — frees ~24 GB
    "10_process_cryosphere.py",
    "11_process_psp_ace.py",
]


def main() -> None:
    setup_logging()
    script_dir = Path(__file__).parent
    LOG.info("=== Solar-Magnetic-Avalanche Processing Pipeline ===")
    LOG.info("Output: %s", PROCESSED_ROOT)
    LOG.info("Disk free: %.1f GB", disk_free_gb())

    failed: list[str] = []

    for script in SCRIPTS:
        path = script_dir / script
        if not path.exists():
            LOG.warning("Script not found, skipping: %s", script)
            failed.append(script)
            continue

        LOG.info("\n--- Running %s ---", script)
        result = subprocess.run(
            [sys.executable, str(path)],
            capture_output=False,
        )
        if result.returncode != 0:
            LOG.error("FAILED: %s (exit %d)", script, result.returncode)
            failed.append(script)
        else:
            LOG.info("OK: %s | Disk free: %.1f GB", script, disk_free_gb())

    LOG.info("\n=== Pipeline complete ===")
    if failed:
        LOG.warning("Failed scripts (%d): %s", len(failed), ", ".join(failed))
        sys.exit(1)
    else:
        LOG.info("All scripts succeeded.")


if __name__ == "__main__":
    main()
