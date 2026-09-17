#!/usr/bin/env python3
"""
download_ssw_compendium_canonical.py
====================================
Build THE canonical major mid-winter SSW catalog used by every analysis in this
project, from the authoritative published source rather than a local detector.

SOURCE (authoritative):
  NOAA CSL Sudden Stratospheric Warming Compendium
  https://csl.noaa.gov/groups/csl8/sswcompendium/majorevents.html
  Butler et al. (2017), Earth Syst. Sci. Data 9, 63-76, doi:10.5194/essd-9-63-2017

The page tabulates the central date of every major mid-winter SSW separately for
six reanalyses (NCEP-NCAR, ERA40, ERA-Interim, JRA-55, MERRA2, ERA5). '****'
means that reanalysis did not register the event as major; those are kept as
nulls so that reanalysis-dependence is explicit and auditable.

WHY THIS EXISTS
  Two conflicting catalogs were in the repo and neither matches the compendium:
    * data/results/ssw_event_catalog.csv (16 events, used by the manuscript)
        - contains 2012-01-11, which is NOT a major SSW in ANY reanalysis
        - omits MAR 2000 and MAR 2010, which ARE major SSWs inside the Davos era
    * data/processed/atmospheric/ssw_catalog.parquet (Butler2015, ends 2021)
  A locally-implemented Charlton-Polvani detector also over-detected (40 events
  1979-2024 vs 28 authoritative), admitting final warmings and double-counts.
  All downstream results are only as good as this catalog, so it is pinned here.

PRIMARY reanalysis = ERA5 (most complete: 42 events), with NCEP-NCAR retained
because this repo's stratospheric wind series is NCEP-derived.

Outputs:
  data/atmospheric/ssw_catalog/majorevents_raw.html   verbatim source
  data/processed/atmospheric/ssw_canonical.csv        tidy canonical catalog
"""
import io
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "atmospheric" / "ssw_catalog" / "majorevents_raw.html"
OUT = ROOT / "data" / "processed" / "atmospheric" / "ssw_canonical.csv"
URL = "https://csl.noaa.gov/groups/csl8/sswcompendium/majorevents.html"

REANALYSES = ["NCEP-NCAR", "ERA40", "ERA-Interim", "JRA-55", "MERRA2", "ERA5"]


def fetch():
    RAW.parent.mkdir(parents=True, exist_ok=True)
    if RAW.exists():
        print(f"[raw] cached -> {RAW}")
        return RAW.read_text(encoding="utf8", errors="replace")
    r = requests.get(URL, timeout=60, headers={"User-Agent": "Solar-Magnetic-Analysis/1.0"})
    r.raise_for_status()
    RAW.write_text(r.text, encoding="utf8")
    print(f"[raw] downloaded -> {RAW}")
    return r.text


def parse(html):
    df = pd.read_html(io.StringIO(html))[0]
    # drop the trailing summary rows ("Total Events", "Frequency: ...")
    df = df[~df["Event Name"].astype(str).str.contains("Total|Frequency", na=False)].copy()

    yr = df["Event Name"].str.extract(r"(\d{4})")[0].astype(int)

    def d(v):
        s = str(v).strip()
        if s in ("****", "nan", "NaN", ""):
            return pd.NaT
        return pd.to_datetime(s, format="%d-%b-%y", errors="coerce")

    def fix_century(dt, label):
        """%y maps 58 -> 2058; the event-name year is the authority.

        A DEC event can have a central date in the following January, so the
        admissible years are {label, label+1}.
        """
        if pd.isna(dt):
            return dt
        for delta in (0, -100, 100):
            if dt.year + delta in (label, label + 1):
                return dt.replace(year=dt.year + delta)
        return pd.NaT

    for c in REANALYSES:
        df[c] = [fix_century(d(v), y) for v, y in zip(df[c], yr)]
        df[c] = pd.to_datetime(df[c])

    for c in REANALYSES:
        bad = df[c].notna() & ~df[c].dt.year.isin(set(yr) | set(yr + 1))
        if bad.any():
            raise ValueError(f"year mismatch in {c}: {df.loc[bad, 'Event Name'].tolist()}")

    out = pd.DataFrame({
        "event_name": df["Event Name"].str.strip(),
        "year_label": yr,
        "enso": df.get("ENSO"),
        "qbo_50mb": df.get("QBO 50mb"),
    })
    for c in REANALYSES:
        out[c.replace("-", "_").lower()] = df[c].dt.strftime("%Y-%m-%d")

    # canonical date: ERA5 where available (most complete), else JRA-55, else NCEP
    pref = ["ERA5", "JRA-55", "NCEP-NCAR", "ERA-Interim", "MERRA2", "ERA40"]
    canon, src = [], []
    for _, r in df.iterrows():
        for c in pref:
            if pd.notna(r[c]):
                canon.append(r[c].strftime("%Y-%m-%d"))
                src.append(c)
                break
        else:
            canon.append(None)
            src.append(None)
    out["date"] = canon
    out["date_source"] = src
    out["n_reanalyses_detecting"] = df[REANALYSES].notna().sum(axis=1).values
    return out.sort_values("date").reset_index(drop=True)


# The NOAA CSL compendium table ends at FEB 2023. Later major mid-winter SSWs
# are added from the peer-reviewed literature, flagged by `date_source` so the
# compendium-derived and literature-derived rows stay distinguishable.
#   Lee, S. H. et al. (2025), "Two major sudden stratospheric warmings during
#   winter 2023/2024", Weather, doi:10.1002/wea.7656
# The March 2025 event is deliberately EXCLUDED: it is documented as the early
# FINAL warming of the 2024/25 winter, and final warmings are outside the
# Charlton-Polvani major mid-winter definition used throughout this project.
POST_COMPENDIUM = [
    {"event_name": "JAN 2024", "date": "2024-01-16", "date_source": "Lee2025_Weather"},
    {"event_name": "MAR 2024", "date": "2024-03-04", "date_source": "Lee2025_Weather"},
]


def main():
    cat = parse(fetch())
    extra = pd.DataFrame(POST_COMPENDIUM)
    extra["year_label"] = 2024
    extra["n_reanalyses_detecting"] = np.nan
    cat = pd.concat([cat, extra], ignore_index=True).sort_values("date").reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cat.to_csv(OUT, index=False)
    print(f"[catalog] {len(cat)} major mid-winter SSWs -> {OUT}")

    d = pd.to_datetime(cat["date"])
    print(f"  span {d.min().date()} .. {d.max().date()}")
    for lo, hi, lab in [("1979-01-01", "2024-12-31", "NCEP era (1979+)"),
                        ("1998-11-01", "2019-05-31", "Davos record era")]:
        m = (d >= lo) & (d <= hi)
        print(f"  {lab}: {int(m.sum())} events")
    print("\n  events in the Davos record era:")
    m = (d >= "1998-11-01") & (d <= "2019-05-31")
    for _, r in cat[m].iterrows():
        print(f"    {r['date']}  {r['event_name']:10s} "
              f"(detected by {r['n_reanalyses_detecting']}/6 reanalyses)")

    # check: the manuscript catalog must be a subset of the canonical one
    paper = ROOT / "data/results/ssw_event_catalog.csv"
    if paper.exists():
        p = pd.to_datetime(pd.read_csv(paper)["date"])
        canon = pd.to_datetime(cat["date"])
        # allow +/-5 d slack for reanalysis-dependent central dates
        unmatched = [str(x.date()) for x in p
                     if (canon - x).abs().min() > pd.Timedelta(days=5)]
        missing = [str(x.date()) for x in canon[m]
                   if (p - x).abs().min() > pd.Timedelta(days=5)]
        print(f"\n  manuscript catalog: {len(p)} events")
        print(f"    NOT major SSWs in any reanalysis : {unmatched or 'none'}")
        print(f"    major SSWs missing from paper    : {missing or 'none'}")


if __name__ == "__main__":
    main()
