# pyright: reportMissingImports=false
# pylint: disable=too-many-lines,import-error,no-name-in-module
"""
Build a post-2022 Avalanche Canada flexible-forecast panel and test whether
Canadian forecast danger is suppressed during recent SSW windows.

The documented archive host on docs.avalanche.ca currently fails on legacy
requests, but the live Avalanche Canada website uses the production services
host below. The dated metadata endpoint exposes one record per active forecast
region together with the centroid and highest danger level. The dated areas
endpoint exposes the matching polygons, which lets us area-weight national
danger summaries despite the flexible daily region geometry.

This script targets the usable flexible-forecast era and recent SSWs with
coverage on that system.
"""

from __future__ import annotations

import json
import threading
import time
import concurrent.futures
from datetime import date, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd  # type: ignore[import-untyped]
import requests  # type: ignore[import-untyped]
from requests.adapters import HTTPAdapter  # type: ignore[import-untyped]
from scipy import stats  # type: ignore[import-untyped]
from urllib3.util.retry import Retry

from _avcan_flexible_utils import (
    build_area_lookup,
    build_panel_row,
    serialize_event_summary,
    serialize_offset_summary,
)

ROOT = Path(__file__).resolve().parents[2]
BASE_URL = "https://avcan-services-api.prod.avalanche.ca/forecasts/en"

# Recent SSWs with public Avalanche Canada flexible-forecast coverage.
SSW_EVENTS = [
    date(2023, 2, 16),
    date(2024, 3, 4),
]

WINDOW_DAYS = 15
ARCHIVE_VALID_START = date(2022, 11, 1)
ARCHIVE_VALID_END = date(2026, 4, 1)
QUERY_HOUR_Z = "20:00:00.000Z"
HIGH_DANGER_THRESHOLD = 3
MAX_WORKERS = 8
REQUEST_TIMEOUT = 45
BOOTSTRAP_ITERATIONS = 5000

PANEL_OUT = (
    ROOT
    / "data"
    / "processed"
    / "cryosphere"
    / "avcan_flexible_ssw_panel.parquet"
)
RESULT_OUT = ROOT / "data" / "results" / "r41_avcan_flexible_validation.json"
_THREAD_LOCAL = threading.local()


def build_session() -> requests.Session:
    retry = Retry(
        total=4,
        backoff_factor=0.8,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=MAX_WORKERS,
        pool_maxsize=MAX_WORKERS,
    )
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def get_session() -> requests.Session:
    session = getattr(_THREAD_LOCAL, "session", None)
    if session is None:
        session = build_session()
        _THREAD_LOCAL.session = session
    return session


def daterange(start: date, end: date):
    current = start
    while current < end:
        yield current
        current += timedelta(days=1)


def event_month_days() -> set[str]:
    month_days: set[str] = set()
    for event_date in SSW_EVENTS:
        for offset in range(-WINDOW_DAYS, WINDOW_DAYS + 1):
            target_date = event_date + timedelta(days=offset)
            month_days.add(target_date.strftime("%m-%d"))
    return month_days


def build_target_valid_dates() -> list[date]:
    month_days = event_month_days()
    return [
        current
        for current in daterange(ARCHIVE_VALID_START, ARCHIVE_VALID_END)
        if current.strftime("%m-%d") in month_days
    ]


def archive_datetime(target_date: date) -> str:
    return f"{target_date.isoformat()}T{QUERY_HOUR_Z}"


def fetch_json(url: str) -> Any:
    session = get_session()
    response = session.get(url, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    return response.json()


def fetch_daily_rows(target_date: date) -> list[dict[str, Any]]:
    dt = archive_datetime(target_date)
    metadata_url = f"{BASE_URL}/metadata?date={dt}"
    areas_url = f"{BASE_URL}/areas?date={dt}"

    metadata = fetch_json(metadata_url)
    areas = fetch_json(areas_url)
    area_lookup = build_area_lookup(areas)
    return [
        build_panel_row(
            target_date,
            dt,
            entry,
            area_lookup,
            HIGH_DANGER_THRESHOLD,
        )
        for entry in metadata
    ]


def download_panel(valid_dates: list[date]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    started = time.time()

    with concurrent.futures.ThreadPoolExecutor(
        max_workers=MAX_WORKERS,
    ) as executor:
        futures_map = {
            executor.submit(fetch_daily_rows, target_date): target_date
            for target_date in valid_dates
        }
        for idx, future in enumerate(
            concurrent.futures.as_completed(futures_map),
            start=1,
        ):
            target_date = futures_map[future]
            try:
                rows.extend(future.result())
            except (requests.RequestException, ValueError) as exc:
                print(
                    f"Daily fetch failed for {target_date.isoformat()}: {exc}"
                )

            if idx % 25 == 0 or idx == len(futures_map):
                elapsed = time.time() - started
                rate = idx / elapsed if elapsed else 0.0
                print(
                    f"Fetched {idx}/{len(futures_map)} Canada dates "
                    f"({rate:.2f} dates/s)"
                )

    panel = pd.DataFrame(rows)
    if panel.empty:
        raise RuntimeError(
            "No Avalanche Canada flexible-forecast rows were parsed."
        )

    panel["date"] = pd.to_datetime(panel["date"])
    panel = panel.sort_values(["date", "region_id"]).reset_index(drop=True)
    return panel


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


def summarize_daily(panel: pd.DataFrame) -> pd.DataFrame:
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
                "mean_danger": float(np.nanmean(dangers)),
                "high_danger_share": float(np.nanmean(high)),
                "weighted_mean_danger": float(weighted_mean),
                "weighted_high_danger_share": float(weighted_high),
                "total_area_sqkm": float(total_area),
            }
        )

    daily = (
        panel.groupby(["date", "month_day"])[
            ["region_id", "area_sqkm", "danger_max", "high_danger"]
        ]
        .apply(summarize)
        .reset_index()
    )
    return label_ssw_windows(daily)


def bootstrap_summary(pair_df: pd.DataFrame) -> dict[str, list[float]]:
    rng = np.random.default_rng(42)
    diff_samples = np.empty(BOOTSTRAP_ITERATIONS)
    share_samples = np.empty(BOOTSTRAP_ITERATIONS)

    for idx in range(BOOTSTRAP_ITERATIONS):
        random_state = int(rng.integers(0, 2**31 - 1))
        sample = pair_df.sample(
            n=len(pair_df),
            replace=True,
            random_state=random_state,
        )
        diff_samples[idx] = sample["diff"].mean()
        share_samples[idx] = sample["high_diff"].mean()

    return {
        "diff_ci": [
            float(x) for x in np.percentile(diff_samples, [2.5, 97.5])
        ],
        "high_diff_ci": [
            float(x) for x in np.percentile(share_samples, [2.5, 97.5])
        ],
    }


def build_matched_pairs(daily: pd.DataFrame) -> pd.DataFrame:
    ssw_daily = daily[daily["is_ssw_window"]].copy()
    control_daily = daily[~daily["is_ssw_window"]].copy()

    control_lookup = (
        control_daily.groupby("month_day").agg(
            control_mean=("weighted_mean_danger", "mean"),
            control_high=("weighted_high_danger_share", "mean"),
            control_n=("date", "nunique"),
        )
    ).reset_index()

    matched = ssw_daily.merge(control_lookup, on="month_day", how="left")
    matched = matched.dropna(subset=["control_mean", "control_high"]).copy()
    if matched.empty:
        raise RuntimeError(
            "No matched control rows were found for Avalanche Canada "
            "flexible forecasts."
        )

    pair_df = matched[[
        "date",
        "ssw_event",
        "event_offset_days",
        "weighted_mean_danger",
        "weighted_high_danger_share",
        "control_mean",
        "control_high",
        "n_regions",
        "control_n",
    ]].copy()
    pair_df["diff"] = pair_df["weighted_mean_danger"] - pair_df["control_mean"]
    pair_df["pct_change"] = pair_df["diff"] / pair_df["control_mean"]
    pair_df["high_diff"] = (
        pair_df["weighted_high_danger_share"] - pair_df["control_high"]
    )
    pair_df["decrease"] = pair_df["diff"] < 0
    return pair_df


def summarize_events(pair_df: pd.DataFrame) -> pd.DataFrame:
    event_summary = (
        pair_df.groupby("ssw_event").agg(
            n_days=("date", "nunique"),
            obs_mean=("weighted_mean_danger", "mean"),
            exp_mean=("control_mean", "mean"),
            obs_high=("weighted_high_danger_share", "mean"),
            exp_high=("control_high", "mean"),
            negative_days=("decrease", "sum"),
        )
    ).reset_index()
    event_summary["diff"] = (
        event_summary["obs_mean"] - event_summary["exp_mean"]
    )
    event_summary["high_diff"] = (
        event_summary["obs_high"] - event_summary["exp_high"]
    )
    return event_summary


def summarize_offsets(pair_df: pd.DataFrame) -> pd.DataFrame:
    offset_summary = (
        pair_df.groupby("event_offset_days").agg(
            n_rows=("date", "size"),
            obs_mean=("weighted_mean_danger", "mean"),
            exp_mean=("control_mean", "mean"),
            obs_high=("weighted_high_danger_share", "mean"),
            exp_high=("control_high", "mean"),
        )
    ).reset_index()
    offset_summary["diff"] = (
        offset_summary["obs_mean"] - offset_summary["exp_mean"]
    )
    offset_summary["high_diff"] = (
        offset_summary["obs_high"] - offset_summary["exp_high"]
    )
    return offset_summary


def build_overall_summary(
    pair_df: pd.DataFrame,
    bootstrap: dict[str, list[float]],
    sign_p: float,
    wilcoxon: Any,
) -> dict[str, Any]:
    n_pairs = int(len(pair_df))
    n_decrease = int(pair_df["decrease"].sum())

    return {
        "n_events": int(len(SSW_EVENTS)),
        "n_daily_pairs": n_pairs,
        "obs_mean_danger": float(pair_df["weighted_mean_danger"].mean()),
        "exp_mean_danger": float(pair_df["control_mean"].mean()),
        "mean_diff": float(pair_df["diff"].mean()),
        "mean_pct_change": float(
            pair_df["diff"].mean() / pair_df["control_mean"].mean()
        ),
        "obs_high_danger_share": float(
            pair_df["weighted_high_danger_share"].mean()
        ),
        "exp_high_danger_share": float(pair_df["control_high"].mean()),
        "mean_high_diff": float(pair_df["high_diff"].mean()),
        "negative_day_fraction": float(n_decrease / n_pairs),
        "sign_test_p": float(sign_p),
        "wilcoxon_statistic": float(wilcoxon.statistic),
        "wilcoxon_p": float(wilcoxon.pvalue),
        "bootstrap_mean_diff_ci": bootstrap["diff_ci"],
        "bootstrap_high_diff_ci": bootstrap["high_diff_ci"],
    }


def analyze_daily(daily: pd.DataFrame) -> dict[str, Any]:
    pair_df = build_matched_pairs(daily)

    n_pairs = int(len(pair_df))
    n_decrease = int(pair_df["decrease"].sum())
    sign_p = stats.binomtest(
        n_decrease,
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
    event_summary = summarize_events(pair_df)
    offset_summary = summarize_offsets(pair_df)
    overall = build_overall_summary(pair_df, bootstrap, sign_p, wilcoxon)

    return {
        "overall": overall,
        "by_event": serialize_event_summary(event_summary),
        "offset_summary": serialize_offset_summary(offset_summary),
    }


def main() -> None:
    PANEL_OUT.parent.mkdir(parents=True, exist_ok=True)
    RESULT_OUT.parent.mkdir(parents=True, exist_ok=True)

    valid_dates = build_target_valid_dates()
    print(f"Target valid dates: {len(valid_dates)}")
    print(f"Event month-days: {len(event_month_days())}")

    panel = download_panel(valid_dates)
    panel.to_parquet(PANEL_OUT, index=False)
    print(f"Saved panel: {PANEL_OUT}")

    daily = summarize_daily(panel)
    results = analyze_daily(daily)

    with RESULT_OUT.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print(f"Saved results: {RESULT_OUT}")

    overall = results["overall"]
    print("\nAvalanche Canada flexible-forecast validation")
    print("=" * 60)
    print(f"Daily pairs: {overall['n_daily_pairs']}")
    print(f"Observed weighted mean danger: {overall['obs_mean_danger']:.3f}")
    print(f"Expected weighted mean danger: {overall['exp_mean_danger']:.3f}")
    print(
        f"Mean difference: {overall['mean_diff']:+.3f} "
        f"({overall['mean_pct_change']:+.1%})"
    )
    print(
        f"Observed weighted high-danger share: "
        f"{overall['obs_high_danger_share']:.3f}"
    )
    print(
        f"Expected weighted high-danger share: "
        f"{overall['exp_high_danger_share']:.3f}"
    )
    print(f"High-danger share difference: {overall['mean_high_diff']:+.3f}")
    print(f"Negative days: {overall['negative_day_fraction']:.1%}")
    print(f"Sign-test p: {overall['sign_test_p']:.4f}")
    print(f"Wilcoxon p: {overall['wilcoxon_p']:.4f}")
    print(f"Bootstrap mean-diff CI: {overall['bootstrap_mean_diff_ci']}")


if __name__ == "__main__":
    main()
