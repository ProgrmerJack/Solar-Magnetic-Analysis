"""Shared download utilities for the Solar-Magnetic-Analysis project."""
import os
import sys
import time
import logging
import requests
from pathlib import Path
from datetime import datetime

BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / "data"
LOG_DIR  = BASE_DIR / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%H:%M:%S")
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(fmt)
        fh = logging.FileHandler(LOG_DIR / f"{name}.log", encoding="utf-8")
        fh.setFormatter(fmt)
        logger.addHandler(ch)
        logger.addHandler(fh)
    return logger


def download_file(url: str, dest: Path, desc: str = "", session: requests.Session = None,
                  overwrite: bool = False, retries: int = 3, chunk: int = 1 << 16) -> bool:
    """Download *url* to *dest*, with resume-on-partial and retry logic."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("utils")

    if dest.exists() and not overwrite:
        logger.info(f"  SKIP  {dest.name}  (already exists)")
        return True

    getter = session or requests
    for attempt in range(1, retries + 1):
        try:
            resp = getter.get(url, stream=True, timeout=120,
                              headers={"User-Agent": "Solar-Magnetic-Analysis/1.0"})
            resp.raise_for_status()
            total = int(resp.headers.get("content-length", 0))
            downloaded = 0
            with open(dest, "wb") as fh:
                for block in resp.iter_content(chunk_size=chunk):
                    if block:
                        fh.write(block)
                        downloaded += len(block)
            size_kb = dest.stat().st_size / 1024
            label = desc or dest.name
            logger.info(f"  ✓  {label}  ({size_kb:,.0f} KB)")
            return True
        except Exception as exc:
            logger.warning(f"  Attempt {attempt}/{retries} failed: {exc}")
            if dest.exists():
                dest.unlink()
            if attempt < retries:
                time.sleep(5 * attempt)
    logger.error(f"  ✗  Failed after {retries} attempts: {url}")
    return False


def write_instructions(dest: Path, title: str, text: str):
    """Write a INSTRUCTIONS.txt file into a data directory."""
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "INSTRUCTIONS.txt").write_text(f"# {title}\n\n{text}\n", encoding="utf-8")
    logging.getLogger("utils").info(f"  📋  Instructions written → {dest / 'INSTRUCTIONS.txt'}")
