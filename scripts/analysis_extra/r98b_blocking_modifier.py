"""R98b: Blocking-as-modifier secondary analysis.

Per the rubber-duck critique, conditioning on the presence of a 500 hPa
blocking episode and asking "given a block, does an SSW modify its
character?" is the only causally honest formulation of the
blocking-as-treatment idea.  This is a HETEROGENEITY analysis, not an
n=16 escape.

Procedure:
1. Detect blocking episodes from NCEP Z500_NH: define the daily
   blocking-strength index as the standardised Z500_NH anomaly relative to
   the day-of-year climatology.  Episodes = consecutive runs of >=5 days
   with z-score >= 1.0 (top quintile).
2. Tag each episode by whether it overlaps the +/-15 d window of any SSW
   (column ssw_within_15d).
3. Compare (a) total natural avalanche count, (b) dry/wet ratio, and (c)
   accident count, between SSW-coincident and SSW-free blocks.
4. Report rate ratios with bootstrap CIs (resample episodes with replacement
   within block category).
5. Explicit caveat: this conditions on a downstream variable; the test
   answers whether the stratospheric event modifies block character given
   that a block has occurred, not whether it produces additional blocks.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "data" / "processed" / "analysis_panel_v2.parquet"
OUT = ROOT / "data" / "results" / "r98b_blocking_modifier.json"

WINTER_MONTHS = (11, 12, 1, 2, 3, 4)
THRESH_Z = 1.0
MIN_DAYS = 5
BOOT = 5000
RNG = np.random.default_rng(11)


def main() -> dict:
    df = pd.read_parquet(PANEL).copy()
    df["date"] = pd.to_datetime(df.index)
    df["month"] = df["date"].dt.month
    df = df[df["month"].isin(WINTER_MONTHS)].sort_values("date").reset_index(drop=True)

    # Day-of-year climatology of Z500_NH
    doy = df["date"].dt.dayofyear
    clim_mean = df.groupby(doy)["ncep_z500_nh"].transform("mean")
    clim_std = df.groupby(doy)["ncep_z500_nh"].transform("std").replace(0, np.nan)
    df["z500_z"] = (df["ncep_z500_nh"] - clim_mean) / clim_std

    df["block_day"] = (df["z500_z"] >= THRESH_Z).astype(int)

    # Identify episodes (consecutive runs)
    runs = []
    current_start = None
    for i, row in df.iterrows():
        if row["block_day"] == 1 and current_start is None:
            current_start = i
        elif row["block_day"] == 0 and current_start is not None:
            length = i - current_start
            if length >= MIN_DAYS:
                runs.append((current_start, i - 1))
            current_start = None
    if current_start is not None:
        length = len(df) - current_start
        if length >= MIN_DAYS:
            runs.append((current_start, len(df) - 1))

    # Tag each episode by whether it overlaps an SSW window
    episodes = []
    for s, e in runs:
        sub = df.iloc[s:e + 1]
        ssw_overlap = int((sub["ssw_within_15d"] == 1).any())
        natural = float(sub["natural_size_234"].sum(skipna=True))
        accidents = float(sub["accident_count"].sum(skipna=True))
        wet = float(sub["wet_natural_size_1234"].sum(skipna=True)) if "wet_natural_size_1234" in sub.columns else np.nan
        dry = float(sub["dry_natural_size_1234"].sum(skipna=True)) if "dry_natural_size_1234" in sub.columns else np.nan
        episodes.append({
            "start": str(sub["date"].iloc[0].date()),
            "end": str(sub["date"].iloc[-1].date()),
            "length_days": int(len(sub)),
            "ssw_overlap": ssw_overlap,
            "natural_count": natural,
            "accident_count": accidents,
            "wet_count": wet,
            "dry_count": dry,
            "winter_id": str(sub["winter_id"].iloc[0]),
        })

    ep_df = pd.DataFrame(episodes)
    n_total = len(ep_df)
    n_ssw = int(ep_df["ssw_overlap"].sum())
    n_free = n_total - n_ssw

    def bootstrap_rate_ratio(metric: str, per_day: bool = True) -> dict:
        a = ep_df[ep_df["ssw_overlap"] == 1]
        b = ep_df[ep_df["ssw_overlap"] == 0]
        if per_day:
            num_a = a[metric].sum() / a["length_days"].sum()
            num_b = b[metric].sum() / b["length_days"].sum()
        else:
            num_a = a[metric].mean()
            num_b = b[metric].mean()
        rr = float(num_a / num_b) if num_b > 0 else np.nan

        boots = []
        for _ in range(BOOT):
            ai = RNG.integers(0, len(a), size=len(a))
            bi = RNG.integers(0, len(b), size=len(b))
            asub = a.iloc[ai]
            bsub = b.iloc[bi]
            if per_day:
                num_a_b = asub[metric].sum() / asub["length_days"].sum()
                num_b_b = bsub[metric].sum() / bsub["length_days"].sum()
            else:
                num_a_b = asub[metric].mean()
                num_b_b = bsub[metric].mean()
            if num_b_b > 0:
                boots.append(num_a_b / num_b_b)
        boots = np.array(boots)
        return {
            "rate_ratio": rr,
            "ci95": [float(np.quantile(boots, 0.025)),
                     float(np.quantile(boots, 0.975))],
            "n_episodes_ssw": int(len(a)),
            "n_episodes_free": int(len(b)),
        }

    natural_rr = bootstrap_rate_ratio("natural_count")
    accident_rr = bootstrap_rate_ratio("accident_count")

    # Dry / wet ratio comparison (only if both present)
    dry_rr = bootstrap_rate_ratio("dry_count") if ep_df["dry_count"].notna().any() else None
    wet_rr = bootstrap_rate_ratio("wet_count") if ep_df["wet_count"].notna().any() else None

    result = {
        "design_caveat": "This is a heterogeneity analysis CONDITIONAL on the occurrence of a persistent Z500 block. It does NOT estimate the unconditional causal effect of an SSW because blocks are themselves a downstream consequence of the stratospheric state (post-treatment). The question answered is: GIVEN that a block occurs, do its surface-hazard signatures differ when an SSW is present vs absent?",
        "block_definition": {
            "z500_threshold_sigma": THRESH_Z,
            "min_consecutive_days": MIN_DAYS,
            "climatology": "Day-of-year mean and standard deviation of NH-mean Z500 from NCEP, computed across the 21-winter panel.",
        },
        "sample": {
            "total_episodes": n_total,
            "ssw_overlap_episodes": n_ssw,
            "ssw_free_episodes": n_free,
            "ssw_overlap_total_days": int(ep_df.loc[ep_df["ssw_overlap"] == 1, "length_days"].sum()),
            "ssw_free_total_days": int(ep_df.loc[ep_df["ssw_overlap"] == 0, "length_days"].sum()),
        },
        "rate_ratios_per_day": {
            "natural_dry_slab_count_size_2plus": natural_rr,
            "accident_count": accident_rr,
            "dry_natural_count_all_sizes": dry_rr,
            "wet_natural_count_all_sizes": wet_rr,
        },
        "interpretation": "If RR < 1 for natural counts: SSW-coincident blocks have suppressed natural avalanching - consistent with the loaded-gun pattern in which the stratospheric event modifies the block character toward cold/dry/dynamic with reduced natural release. The accident-count direction reflects the consequence channel.",
        "framing_note": "Reported as supportive heterogeneity within the blocking population. The unconditional event-catalog analysis (gmRR=0.32) and the continuous-vortex framework (R98) remain the primary identification strategies.",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    main()
