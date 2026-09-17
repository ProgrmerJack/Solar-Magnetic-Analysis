#!/usr/bin/env python3
"""
acquire_enso.py
===============
Acquires the ONI (Oceanic Nino Index, 3-month running Nino3.4 SST anomaly) from
NOAA PSL. Open, no authentication.

WHY THIS IS NEEDED NOW
  The AAO negative control fails: -0.659 at +45..+60 d, p=0.015, stable across
  60 fixed seeds, with the winter-block bootstrap verified calibrated at that
  bin (boot_calibration.py: SE ratio 1.06, coverage 0.945) and the clean
  pseudo-onset null centred on zero (aao_null.py: mean -0.018). The estimator is
  not manufacturing the signal, so something links NH SSW onsets to the Southern
  annular mode.

  ENSO is the obvious common cause. It modulates NH polar vortex variability
  (and therefore SSW likelihood) AND projects onto the SAM. If so, the AAO is
  not a valid negative control for SSW event studies, and Gate 6's premise is
  wrong rather than the estimator.

Outputs (this directory):
  enso_oni.parquet        monthly ONI, tidy
  enso_checksum.txt       sha256 of the raw source
"""
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
URL = "https://psl.noaa.gov/data/correlation/oni.data"
OUT = HERE / "enso_oni.parquet"
SUM = HERE / "enso_checksum.txt"
# The file states its own missing-value sentinel on a bare line after the data
# (currently " -99.9"; the sibling nina34.data uses -99.99). Parse it rather than
# hard-coding: a wrong sentinel silently admits 2026's unfilled months as real
# ENSO values of -99.9, which is what the physical-range assert below caught.
MISSING_FALLBACK = -99.9
MISSING_TOL = 0.05


def main():
    r = requests.get(URL, timeout=120)
    r.raise_for_status()
    raw = r.content
    lines = r.text.splitlines()
    y0, y1 = (int(x) for x in lines[0].split()[:2])

    # the sentinel sits on a bare single-number line in the trailer
    missing = MISSING_FALLBACK
    for ln in lines[1:]:
        p = ln.split()
        if len(p) == 1:
            try:
                v = float(p[0])
            except ValueError:
                continue
            if v < -50:
                missing = v
                break
    print(f"missing-value sentinel read from file: {missing}")

    rows = []
    for ln in lines[1:]:
        p = ln.split()
        if len(p) != 13:
            continue
        try:
            yr = int(p[0])
        except ValueError:
            continue
        if not (y0 <= yr <= y1):
            continue
        for m, v in enumerate(p[1:], start=1):
            val = float(v)
            if abs(val - missing) < MISSING_TOL or val < -50:
                continue
            rows.append({"date": pd.Timestamp(yr, m, 15), "oni": val})

    df = pd.DataFrame(rows).sort_values("date").reset_index(drop=True)
    assert len(df) > 600, f"only {len(df)} monthly ONI values parsed"
    assert df["oni"].abs().max() < 5, "ONI out of physical range"
    df.to_parquet(OUT, index=False)

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    SUM.write_text(f"{hashlib.sha256(raw).hexdigest()}  oni.data\n"
                   f"# {URL}, retrieved {stamp}\n", encoding="utf8")
    print(f"{len(df)} monthly ONI values "
          f"{df['date'].min().date()}..{df['date'].max().date()}")
    print(f"  range {df['oni'].min():+.2f}..{df['oni'].max():+.2f}, "
          f"sd {df['oni'].std():.2f}")
    print(f"-> {OUT.name}, {SUM.name}")


if __name__ == "__main__":
    main()
