#!/usr/bin/env python3
"""
validate_ssw_detector.py
========================
Does `ensemble_precursor.detect_ssw`, the detector every CMIP6 analysis uses,
reproduce the published Charlton & Polvani (2007) central dates when it is run on
the reanalysis those dates were computed from?

Test: NCEP/NCAR R1 daily zonal-mean u at 10 hPa, 60N, 1958-2024, fetched by the
existing NCSS reducer (`03_data_ingestion/extend_ncep_presatellite.series`), run
through `detect_ssw`, and compared with the NCEP-NCAR column of the frozen NOAA
CSL compendium (`build_catalogue.parse_compendium`). The compendium's NCEP dates
are CP07 dates on this reanalysis, so a correct detector matches them day for
day. Detections after the compendium's last season are reported separately,
not scored.

Pass: every published date found within +-1 day, and no detection inside the
compendium's coverage that the compendium does not list.

Output: results/current/1_catalogue/validate_ssw_detector.json
"""
import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "current" / "1_catalogue"
RESULTS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE.parents[0] / "02_event_catalogues"))
sys.path.insert(0, str(HERE.parents[0] / "03_data_ingestion"))
sys.path.insert(0, str(HERE.parents[0] / "07_physical_decomposition"))
import build_catalogue as BC                        # noqa: E402
import extend_ncep_presatellite as EN               # noqa: E402
import ensemble_precursor as EP                     # noqa: E402

YEARS = range(1958, 2025)
LEVEL = 10
COL = f"uwnd_ms_60N_{LEVEL}hPa"
TOL_DAYS = 1


def u10_60n():
    with ThreadPoolExecutor(EN.WORKERS) as ex:
        parts = list(ex.map(lambda y: EN.series("uwnd", y, LEVEL), YEARS))
    s = pd.concat(parts)[COL].sort_index()
    s.index = pd.DatetimeIndex(s.index).tz_localize(None).normalize()
    return s[~s.index.duplicated()]


def main():
    raw = BC.RAW.read_bytes()
    comp = BC.parse_compendium(raw.decode("utf8", errors="replace"))
    pub = pd.DatetimeIndex(comp["NCEP-NCAR"].dropna().sort_values())
    last_season_end = pd.Timestamp(year=BC.winter_of(pub.max()), month=4, day=30)

    u = u10_60n()
    gaps = int((np.diff(u.index.values).astype("timedelta64[D]").astype(int) > 1).sum())
    print(f"NCEP R1 u10 60N: {u.index.min().date()} .. {u.index.max().date()}, "
          f"{len(u):,} days, {gaps} gaps")
    det = EP.detect_ssw(u.values, u.index)

    scored = det[(det >= u.index.min()) & (det <= last_season_end)]
    pub_in = pub[pub >= u.index.min()]
    match, missed = [], []
    for d in pub_in:
        off = (scored - d).days
        k = np.flatnonzero(np.abs(off) <= TOL_DAYS)
        if len(k):
            match.append((d, int(off[k[0]])))
        else:
            missed.append(d)
    matched_det = {d + pd.Timedelta(days=o) for d, o in match}
    extra = [d for d in scored if d not in matched_det]
    after = det[det > last_season_end]

    exact = sum(1 for _, o in match if o == 0)
    print(f"published NCEP-NCAR dates: {len(pub_in)}  "
          f"(compendium through the {BC.winter_of(pub.max())} season)")
    print(f"matched within +-{TOL_DAYS} d: {len(match)}  (exact: {exact})")
    print(f"missed: {[str(d.date()) for d in missed]}")
    print(f"extra detections inside coverage: {[str(d.date()) for d in extra]}")
    print(f"detections after coverage (not scored): {[str(d.date()) for d in after]}")
    ok = not missed and not extra
    print("PASS" if ok else "FAIL")

    res = {"question": "does detect_ssw reproduce the published CP07 NCEP-NCAR "
                       "central dates on NCEP R1?",
           "passes": ok,
           "inputs": {"compendium": str(BC.RAW.relative_to(ROOT)),
                      "compendium_sha256": hashlib.sha256(raw).hexdigest(),
                      "reanalysis": "NCEP/NCAR R1 daily mean u, 10 hPa, zonal mean "
                                    "at 60N, via extend_ncep_presatellite.series",
                      "period": [str(u.index.min().date()), str(u.index.max().date())],
                      "n_days": int(len(u)), "n_gaps": gaps},
           "tolerance_days": TOL_DAYS,
           "coverage_end": str(last_season_end.date()),
           "n_published": int(len(pub_in)), "n_matched": len(match),
           "n_exact": exact,
           "missed": [str(d.date()) for d in missed],
           "extra_in_coverage": [str(d.date()) for d in extra],
           "after_coverage": [str(d.date()) for d in after],
           "offsets_days": {str(d.date()): o for d, o in match if o != 0}}
    (RESULTS / "validate_ssw_detector.json").write_text(
        json.dumps(res, indent=2), encoding="utf8", newline="\n")
    print("Saved -> validate_ssw_detector.json")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
