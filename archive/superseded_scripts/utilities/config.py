"""Project-level configuration: derive ROOT from this file's location."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # scripts/ → repo root
DATA = ROOT / 'data'
PROCESSED = DATA / 'processed'
RESULTS = DATA / 'results'
