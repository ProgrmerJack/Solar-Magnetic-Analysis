"""
02_define_events.py — Build the geomagnetic disturbance event catalog
=====================================================================
Outputs event catalog + summary statistics for the paper.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _analysis_utils import PROCESSED, RESULTS, LOG, load_panel

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def build_event_catalog():
    panel = load_panel(winter_only=False)

    # All events (including summer, for falsification)
    events_all = panel[panel["geo_event"] == 1].copy()
    events_all["is_winter_event"] = events_all["is_winter"]

    # Event properties
    events_all["event_id"] = range(1, len(events_all) + 1)

    catalog = events_all[[
        "event_id", "kp_max", "dst_min", "winter_id", "month",
        "is_winter_event", "ssw_within_15d",
    ]].copy()

    if "flare_count" in events_all.columns:
        catalog["flare_count"] = events_all["flare_count"]
    if "flare_max_class" in events_all.columns:
        catalog["flare_max_class"] = events_all["flare_max_class"]
    if "sw_bz_min" in events_all.columns:
        catalog["sw_bz_min"] = events_all["sw_bz_min"]

    # Save
    out = RESULTS / "event_catalog.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_parquet(out)
    LOG.info("Event catalog: %d total events, %d winter events",
             len(catalog), catalog["is_winter_event"].sum())

    # Summary
    w = catalog[catalog["is_winter_event"] == 1]
    print(f"\n{'='*60}")
    print("EVENT CATALOG SUMMARY")
    print(f"{'='*60}")
    print(f"Total events (all seasons): {len(catalog)}")
    print(f"Winter (NDJFM) events: {len(w)}")
    print(f"Summer events (for falsification): {len(catalog) - len(w)}")
    print(f"\nWinter event Kp stats:")
    print(f"  Mean Kp_max: {w['kp_max'].mean():.1f}")
    print(f"  Median Kp_max: {w['kp_max'].median():.1f}")
    print(f"  Max Kp_max: {w['kp_max'].max():.1f}")
    print(f"\nWinter event Dst stats:")
    print(f"  Mean Dst_min: {w['dst_min'].mean():.0f} nT")
    print(f"  Median Dst_min: {w['dst_min'].median():.0f} nT")
    print(f"  Min Dst_min: {w['dst_min'].min():.0f} nT")
    print(f"\nEvents by winter:")
    by_winter = w.groupby("winter_id").size()
    for wid, cnt in by_winter.items():
        print(f"  {wid}: {cnt} events")
    print(f"\nEvents near SSW (±15d): {w['ssw_within_15d'].sum()}")

    return catalog


if __name__ == "__main__":
    build_event_catalog()
