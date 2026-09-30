#!/usr/bin/env python3
"""
build_catalogue.py
==================
Builds the frozen event catalogue from source, exactly as specified in
CATALOGUE_SPECIFICATION.md. Rebuilding from the raw compendium table must
reproduce it byte-for-byte.

This is the ONLY place an SSW event list is constructed. Every analysis imports
`load_catalogue()` from here. Inline date lists are forbidden -- four mutually
incompatible catalogues in the superseded project came from exactly that.

Primary rule: consensus >= 2/3 of the reanalyses COVERING that date (era-
independent); central date = median of the
qualifying dates. Sensitivity sets are built in the same pass so they can never
drift apart.

Outputs (this directory):
  event_catalogue.csv          machine-readable audit table
  catalogue_sets.json          the frozen id -> date-list mapping
  catalogue_checksum.txt       sha256 of the raw source and of the built table

Run `python build_catalogue.py --test` for the §9 invariance tests.
"""
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results" / "current" / "1_catalogue"
RESULTS.mkdir(parents=True, exist_ok=True)
ROOT = HERE.parents[1]
RAW = ROOT / "data" / "atmospheric" / "ssw_catalog" / "majorevents_raw.html"
OUT_CSV = HERE / "event_catalogue.csv"
OUT_SETS = RESULTS / "catalogue_sets.json"
OUT_SUM = HERE / "catalogue_checksum.txt"

REANALYSES = ["NCEP-NCAR", "ERA40", "ERA-Interim", "JRA-55", "MERRA2", "ERA5"]
COL = {r: r.replace("-", "_").lower() for r in REANALYSES}

# Temporal coverage of each reanalysis. This matters: a count out of six is NOT
# comparable across eras, because only four products span 1958-1978. Under an
# absolute ">=4 of 6" rule, an early event must be detected UNANIMOUSLY while a
# 1979-2002 event needs only 67% -- so the rule silently becomes stricter the
# further back you go, and "marginal" turns into a proxy for "pre-satellite".
# Measured on the previous build: mean consensus 3.3 (1958-78) vs 5.4 (1979-2002),
# and 6 of 8 "marginal" events were pre-1979 (Fisher P = 0.011).
# Consensus is therefore a FRACTION of the products that actually cover the date.
COVERAGE = {
    "NCEP-NCAR":   (1948, 2100),
    "ERA40":       (1957, 2002),
    "ERA-Interim": (1979, 2019),
    "JRA-55":      (1958, 2100),
    "MERRA2":      (1980, 2100),
    "ERA5":        (1940, 2100),
}
CONSENSUS_FRACTION = 2 / 3       # >=2/3 of covering products, era-independent
SEASON_START, SEASON_END = (11, 1), (3, 31)      # 1 Nov .. 31 Mar

# Post-compendium events (spec §7). Literature-sourced; no per-reanalysis counts.
POST_COMPENDIUM = [
    {"event_name": "JAN 2024", "date": "2024-01-16", "source": "Lee2025_Weather"},
    {"event_name": "MAR 2024", "date": "2024-03-04", "source": "Lee2025_Weather"},
]
# Documented and deliberately excluded, with the reason recorded (spec §1, §7)
EXCLUDED = [
    {"event_name": "MAR 2025", "date": "2025-03-09", "source": "GMAO_2025",
     "exclusion_reason": "final warming of the 2024/25 winter; §1 excludes final warmings"},
]


def sha(b):
    return hashlib.sha256(b).hexdigest()


def parse_compendium(html):
    df = pd.read_html(io.StringIO(html))[0]
    df = df[~df["Event Name"].astype(str).str.contains("Total|Frequency", na=False)].copy()
    yr = df["Event Name"].str.extract(r"(\d{4})")[0].astype(int)

    def one(v, label):
        s = str(v).strip()
        if s in ("****", "nan", "NaN", ""):
            return pd.NaT
        d = pd.to_datetime(s, format="%d-%b-%y", errors="coerce")
        if pd.isna(d):
            return pd.NaT
        # two-digit years: the event-name year is the authority; a DEC event may
        # have a central date in the following January
        for delta in (0, -100, 100):
            if d.year + delta in (label, label + 1):
                return d.replace(year=d.year + delta)
        return pd.NaT

    for r in REANALYSES:
        df[r] = [one(v, y) for v, y in zip(df[r], yr)]
        df[r] = pd.to_datetime(df[r])
    df["year_label"] = yr
    return df


def in_season(d):
    if pd.isna(d):
        return False
    return (d.month, d.day) >= SEASON_START or (d.month, d.day) <= SEASON_END


def winter_of(d):
    return d.year + 1 if d.month >= 11 else d.year


def build():
    raw = RAW.read_bytes()
    df = parse_compendium(raw.decode("utf8", errors="replace"))

    rows = []
    for _, r in df.iterrows():
        dates = {x: r[x] for x in REANALYSES if pd.notna(r[x])}
        n = len(dates)
        med = None
        if n:
            # median of qualifying dates, taking the earlier day on a tie (spec §2)
            ordered = sorted(dates.values())
            med = ordered[(n - 1) // 2]
        yr = int(r["year_label"])
        n_avail = sum(1 for x in REANALYSES
                      if COVERAGE[x][0] <= yr <= COVERAGE[x][1])
        frac = n / n_avail if n_avail else 0.0
        row = {
            "event_name": str(r["Event Name"]).strip(),
            "winter": winter_of(med) if med is not None else None,
            "primary_central_date": med.strftime("%Y-%m-%d") if med is not None else None,
            "event_type": "not_classified",
            "final_warming_flag": False,
            "consensus_count": n,
            "n_reanalyses_available": n_avail,
            "consensus_fraction": round(frac, 3),
            "source": "NOAA_CSL_compendium",
            "exclusion_reason": "",
        }
        for x in REANALYSES:
            row[COL[x]] = dates[x].strftime("%Y-%m-%d") if x in dates else ""
        row["primary_inclusion"] = bool(frac >= CONSENSUS_FRACTION)
        if not row["primary_inclusion"]:
            row["exclusion_reason"] = (f"consensus {n}/{n_avail} = {frac:.2f} "
                                       f"< {CONSENSUS_FRACTION:.2f}")
        rows.append(row)

    for e in POST_COMPENDIUM:
        d = pd.Timestamp(e["date"])
        row = {"event_name": e["event_name"], "winter": winter_of(d),
               "primary_central_date": e["date"], "event_type": "not_classified",
               "final_warming_flag": False, "consensus_count": None,
               "n_reanalyses_available": None, "consensus_fraction": None,
               "source": e["source"], "exclusion_reason": "",
               "primary_inclusion": True}
        for x in REANALYSES:
            row[COL[x]] = ""
        rows.append(row)

    for e in EXCLUDED:
        d = pd.Timestamp(e["date"])
        row = {"event_name": e["event_name"], "winter": winter_of(d),
               "primary_central_date": e["date"], "event_type": "final_warming",
               "final_warming_flag": True, "consensus_count": None,
               "n_reanalyses_available": None, "consensus_fraction": None,
               "source": e["source"], "exclusion_reason": e["exclusion_reason"],
               "primary_inclusion": False}
        for x in REANALYSES:
            row[COL[x]] = ""
        rows.append(row)

    cat = pd.DataFrame(rows)
    cat = cat[cat["primary_central_date"].notna()]
    cat = cat.sort_values("primary_central_date").reset_index(drop=True)
    cat.insert(0, "event_id", [f"SSW{i+1:03d}" for i in range(len(cat))])
    cat["sha256_of_source"] = sha(raw)[:16]

    sets = {
        "primary": cat.loc[cat["primary_inclusion"] & ~cat["final_warming_flag"],
                           "primary_central_date"].tolist(),
        "primary_compendium_only": cat.loc[
            cat["primary_inclusion"] & ~cat["final_warming_flag"]
            & (cat["source"] == "NOAA_CSL_compendium"), "primary_central_date"].tolist(),
        "consensus_strict": cat.loc[(cat["consensus_fraction"].fillna(0) >= 0.99)
                                    & ~cat["final_warming_flag"],
                                    "primary_central_date"].tolist(),
        "consensus_half": cat.loc[(cat["consensus_fraction"].fillna(0) >= 0.5)
                                  & ~cat["final_warming_flag"],
                                  "primary_central_date"].tolist(),
        "union": cat.loc[(cat["consensus_count"].fillna(1) >= 1)
                         & ~cat["final_warming_flag"], "primary_central_date"].tolist(),
    }
    for x in REANALYSES:
        sets[COL[x]] = sorted(d for d in cat[COL[x]] if d)

    # lineterminator is pinned: pandas defaults to os.linesep, so the same
    # catalogue written on Windows (CRLF) and Linux (LF) differ by exactly one
    # byte per row and the checksum below stops matching. The frozen file was
    # built on Windows; without this the catalogue does not rebuild
    # byte-for-byte on Linux even though every value is identical.
    cat.to_csv(OUT_CSV, index=False, lineterminator="\n")
    OUT_SETS.write_text(json.dumps(sets, indent=1), encoding="utf8")
    OUT_SUM.write_text(
        f"source_html_sha256 {sha(raw)}\n"
        f"event_catalogue_csv_sha256 {sha(OUT_CSV.read_bytes())}\n", encoding="utf8")
    return cat, sets


def load_catalogue(which="primary"):
    """The single entry point every analysis must use."""
    sets = json.loads(OUT_SETS.read_text(encoding="utf8"))
    if which not in sets:
        raise KeyError(f"unknown catalogue '{which}'; have {sorted(sets)}")
    return pd.DatetimeIndex(pd.to_datetime(sets[which]))


# ------------------------------------------------------------------ tests ----
def run_tests(cat, sets):
    """Spec §9 invariance tests. Any failure blocks the freeze."""
    fails = []
    prim = cat[cat["primary_inclusion"] & ~cat["final_warming_flag"]]
    d = pd.to_datetime(prim["primary_central_date"])

    if cat["event_id"].duplicated().any():
        fails.append("duplicate event_id")
    if d.duplicated().any():
        fails.append("duplicate central date in primary")
    gaps = d.sort_values().diff().dt.days.dropna()
    if (gaps < 20).any():
        fails.append(f"minimum separation violated: {int((gaps < 20).sum())} pairs < 20 d")
    if prim["final_warming_flag"].any():
        fails.append("a final warming is included in primary")
    bad = [x.strftime("%Y-%m-%d") for x in d if not in_season(x)]
    if bad:
        fails.append(f"central date outside 1 Nov–31 Mar: {bad}")
    if not set(prim["source"]) <= {"NOAA_CSL_compendium", "Lee2025_Weather"}:
        fails.append("primary contains an unrecognised source")
    if len(sets["primary"]) != len(prim):
        fails.append("catalogue_sets primary length disagrees with the table")
    if not set(sets["primary_compendium_only"]) <= set(sets["primary"]):
        fails.append("compendium-only set is not a subset of primary")
    if not set(sets["consensus_strict"]) <= set(sets["union"]):
        fails.append("consensus_strict is not a subset of union")
    # literature-sourced events have no consensus fraction, so the subset
    # relation is only meaningful on the compendium-derived part
    if not set(sets["primary_compendium_only"]) <= set(sets["consensus_half"]):
        fails.append("compendium part of primary is not a subset of consensus_half")

    print("\n=== invariance tests (spec §9) ===")
    for f in fails:
        print(f"  FAIL  {f}")
    if not fails:
        print("  all passed")
    return fails


def main():
    cat, sets = build()
    prim = cat[cat["primary_inclusion"] & ~cat["final_warming_flag"]]
    d = pd.to_datetime(prim["primary_central_date"])
    print(f"catalogue built: {len(cat)} rows, {len(prim)} in primary")
    print(f"  event period : {d.min().date()} .. {d.max().date()}")
    print(f"  N events     : {len(prim)}")
    print(f"  N unique winters: {prim['winter'].nunique()}")
    multi = prim.groupby("winter").size()
    print(f"  winters with >1 event: {int((multi > 1).sum())} "
          f"({sorted(multi[multi > 1].index.astype(int))})")
    print("\n  set sizes:")
    for k, v in sets.items():
        print(f"    {k:24s} {len(v)}")
    fails = run_tests(cat, sets)
    print(f"\nwrote {OUT_CSV.name}, {OUT_SETS.name}, {OUT_SUM.name}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
