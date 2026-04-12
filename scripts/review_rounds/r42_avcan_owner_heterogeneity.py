# pyright: reportMissingImports=false
# pylint: disable=import-error,no-name-in-module
"""
Quantify owner-level heterogeneity inside the Avalanche Canada flexible panel.

The national area-weighted summary in r41 is useful as a stress test, but it can
hide internal structure because forecast domains are uneven in spatial extent.
This script reuses the saved flexible-region panel and resolves the SSW signal by
forecast owner.
"""

from __future__ import annotations

import json
from datetime import date, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd  # type: ignore[import-untyped]
from scipy import stats  # type: ignore[import-untyped]

from _avcan_flexible_utils import serialize_event_summary

ROOT = Path(__file__).resolve().parents[2]
PANEL_PATH = (
    ROOT
    / "data"
    / "processed"
    / "cryosphere"
    / "avcan_flexible_ssw_panel.parquet"
)
RESULT_OUT = ROOT / "data" / "results" / "r42_avcan_owner_heterogeneity.json"

SSW_EVENTS = [
    date(2023, 2, 16),
    date(2024, 3, 4),
]
WINDOW_DAYS = 15
BOOTSTRAP_ITERATIONS = 5000


def label_ssw_windows(daily: pd.DataFrame) -> pd.DataFrame:
    daily = daily.copy()
    daily["ssw_event"] = pd.NaT
    daily["event_offset_days"] = np.nan

    for event_date in SSW_EVENTS:
        start = pd.Timestamp(event_date - timedelta(days=WINDOW_DAYS))
        end = pd.Timestamp(event_date + timedelta(days=WINDOW_DAYS))
        mask = daily["date"].between(start, end)
        daily.loc[mask, "ssw_event"] = pd.Timestamp(event_date)
        offsets = daily.loc[mask, "date"] - pd.Timestamp(event_date)
        daily.loc[mask, "event_offset_days"] = offsets.dt.days

    daily["is_ssw_window"] = daily["ssw_event"].notna()
    return daily


def summarize_owner_daily(panel: pd.DataFrame) -> pd.DataFrame:
    def summarize(group: pd.DataFrame) -> pd.Series:
        weights = group["area_sqkm"].to_numpy(dtype=float)
        dangers = group["danger_max"].to_numpy(dtype=float)
        high = group["high_danger"].to_numpy(dtype=float)
        positive_weights = np.isfinite(weights) & (weights > 0)

        if positive_weights.any():
            use_weights = weights[positive_weights]
            use_dangers = dangers[positive_weights]
            use_high = high[positive_weights]
            weighted_mean = np.average(use_dangers, weights=use_weights)
            weighted_high = np.average(use_high, weights=use_weights)
            total_area = use_weights.sum()
        else:
            weighted_mean = float(np.nanmean(dangers))
            weighted_high = float(np.nanmean(high))
            total_area = float("nan")

        return pd.Series(
            {
                "n_regions": int(len(group)),
                "weighted_mean_danger": float(weighted_mean),
                "weighted_high_danger_share": float(weighted_high),
                "total_area_sqkm": float(total_area),
            }
        )

    daily = (
        panel.groupby(["owner", "date", "month_day"])[
            ["region_id", "area_sqkm", "danger_max", "high_danger"]
        ]
        .apply(summarize)
        .reset_index()
    )
    return label_ssw_windows(daily)


def owner_area_shares(panel: pd.DataFrame) -> dict[str, float]:
    national_daily_area = panel.groupby("date")["area_sqkm"].sum().rename("national_area")
    owner_daily_area = (
        panel.groupby(["owner", "date"])["area_sqkm"].sum().rename("owner_area")
    ).reset_index()
    owner_daily_area = owner_daily_area.merge(
        national_daily_area.reset_index(),
        on="date",
        how="left",
    )
    owner_daily_area["area_share"] = (
        owner_daily_area["owner_area"] / owner_daily_area["national_area"]
    )
    shares = owner_daily_area.groupby("owner")["area_share"].mean().sort_values(
        ascending=False
    )
    return {str(owner): float(value) for owner, value in shares.items()}


def bootstrap_summary(pair_df: pd.DataFrame) -> dict[str, list[float]]:
    rng = np.random.default_rng(42)
    diff_samples = np.empty(BOOTSTRAP_ITERATIONS)
    high_samples = np.empty(BOOTSTRAP_ITERATIONS)

    for idx in range(BOOTSTRAP_ITERATIONS):
        random_state = int(rng.integers(0, 2**31 - 1))
        sample = pair_df.sample(
            n=len(pair_df),
            replace=True,
            random_state=random_state,
        )
        diff_samples[idx] = sample["diff"].mean()
        high_samples[idx] = sample["high_diff"].mean()

    return {
        "diff_ci": [float(x) for x in np.percentile(diff_samples, [2.5, 97.5])],
        "high_diff_ci": [
            float(x) for x in np.percentile(high_samples, [2.5, 97.5])
        ],
    }


def build_owner_pairs(owner_daily: pd.DataFrame) -> pd.DataFrame:
    ssw_daily = owner_daily[owner_daily["is_ssw_window"]].copy()
    control_daily = owner_daily[~owner_daily["is_ssw_window"]].copy()
    lookup = control_daily.groupby("month_day").agg(
        ctrl_mean=("weighted_mean_danger", "mean"),
        ctrl_high=("weighted_high_danger_share", "mean"),
    )
    pair_df = ssw_daily.merge(
        lookup,
        on="month_day",
        how="left",
    ).dropna()
    if pair_df.empty:
        raise RuntimeError("No matched-control rows found for owner subset.")

    pair_df["diff"] = pair_df["weighted_mean_danger"] - pair_df["ctrl_mean"]
    pair_df["high_diff"] = (
        pair_df["weighted_high_danger_share"] - pair_df["ctrl_high"]
    )
    pair_df["decrease"] = pair_df["diff"] < 0
    return pair_df


def build_owner_event_summary(pair_df: pd.DataFrame) -> pd.DataFrame:
    event_summary = (
        pair_df.groupby("ssw_event").agg(
            n_days=("date", "nunique"),
            obs_mean=("weighted_mean_danger", "mean"),
            exp_mean=("ctrl_mean", "mean"),
            obs_high=("weighted_high_danger_share", "mean"),
            exp_high=("ctrl_high", "mean"),
            negative_days=("decrease", "sum"),
        )
    ).reset_index()
    event_summary["diff"] = event_summary["obs_mean"] - event_summary["exp_mean"]
    event_summary["high_diff"] = event_summary["obs_high"] - event_summary["exp_high"]
    return event_summary


def summarize_owner_pairs(
    pair_df: pd.DataFrame,
    owner_daily: pd.DataFrame,
    area_share: float,
) -> dict[str, Any]:

    n_pairs = int(len(pair_df))
    n_negative = int(pair_df["decrease"].sum())
    sign_p = stats.binomtest(
        n_negative,
        n_pairs,
        0.5,
        alternative="greater",
    ).pvalue
    wilcoxon = stats.wilcoxon(
        pair_df["diff"],
        alternative="less",
        zero_method="wilcox",
    )
    bootstrap = bootstrap_summary(pair_df)
    event_summary = build_owner_event_summary(pair_df)

    return {
        "n_panel_rows": int(len(owner_daily)),
        "n_dates": int(owner_daily["date"].nunique()),
        "n_regions": int(owner_daily["n_regions"].max()),
        "mean_daily_area_share": float(area_share),
        "obs_mean_danger": float(pair_df["weighted_mean_danger"].mean()),
        "exp_mean_danger": float(pair_df["ctrl_mean"].mean()),
        "mean_diff": float(pair_df["diff"].mean()),
        "obs_high_danger_share": float(pair_df["weighted_high_danger_share"].mean()),
        "exp_high_danger_share": float(pair_df["ctrl_high"].mean()),
        "mean_high_diff": float(pair_df["high_diff"].mean()),
        "n_daily_pairs": n_pairs,
        "negative_days": n_negative,
        "negative_day_fraction": float(n_negative / n_pairs),
        "sign_test_p": float(sign_p),
        "wilcoxon_statistic": float(wilcoxon.statistic),
        "wilcoxon_p": float(wilcoxon.pvalue),
        "bootstrap_mean_diff_ci": bootstrap["diff_ci"],
        "bootstrap_high_diff_ci": bootstrap["high_diff_ci"],
        "by_event": serialize_event_summary(event_summary),
    }


def analyze_owner(owner_daily: pd.DataFrame, area_share: float) -> dict[str, Any]:
    pair_df = build_owner_pairs(owner_daily)
    return summarize_owner_pairs(pair_df, owner_daily, area_share)


def main() -> None:
    panel = pd.read_parquet(PANEL_PATH)
    panel["date"] = pd.to_datetime(panel["date"])
    panel["owner"] = panel["owner"].fillna("UNKNOWN")

    daily = summarize_owner_daily(panel)
    shares = owner_area_shares(panel)

    owners = []
    for owner_name, area_share in shares.items():
        owner_daily = daily[daily["owner"] == owner_name].copy()
        owner_summary = analyze_owner(owner_daily, area_share)
        owner_summary["owner"] = owner_name
        owners.append(owner_summary)

    output = {
        "source_panel": str(PANEL_PATH.relative_to(ROOT)).replace("\\", "/"),
        "n_owner_systems": len(owners),
        "owners": owners,
    }

    RESULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    with RESULT_OUT.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)

    print("Avalanche Canada owner heterogeneity")
    print("=" * 60)
    for owner in owners:
        print(
            f"{owner['owner']}: area_share={owner['mean_daily_area_share']:.1%}, "
            f"diff={owner['mean_diff']:+.3f}, "
            f"negative_days={owner['negative_days']}/{owner['n_daily_pairs']}, "
            f"sign_p={owner['sign_test_p']:.4g}"
        )
    print(f"Saved results: {RESULT_OUT}")


if __name__ == "__main__":
    main()