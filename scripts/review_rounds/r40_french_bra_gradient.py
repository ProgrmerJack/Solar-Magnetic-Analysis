"""
Quantify internal spatial structure in the French BRA SSW replication.

This script reuses the all-massif French BRA panel produced by
scripts/review_rounds/r39_french_bra_replication.py and tests whether the
weaker all-France pooled signal is explained by a north-south gradient.
"""

from __future__ import annotations

import json
import re
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd  # type: ignore[import-untyped]
from scipy import stats  # type: ignore[import-untyped]

ROOT = Path(__file__).resolve().parents[2]
PANEL_PATH = (
    ROOT
    / "data"
    / "processed"
    / "cryosphere"
    / "french_bra_all_ssw_panel.parquet"
)
META_PATH = (
    ROOT
    / "data"
    / "cryosphere"
    / "france"
    / "Metadata_massif_DP_2024_20240215.geojson"
)
RESULT_OUT = ROOT / "data" / "results" / "r40_french_bra_gradient.json"

SSW_EVENTS = [
    date(2018, 2, 12),
    date(2019, 1, 1),
    date(2021, 1, 5),
    date(2023, 2, 16),
]
WINDOW_DAYS = 15
HIGH_DANGER_THRESHOLD = 3
BOOTSTRAP_ITERATIONS = 5000
EXPECTED_ALL_MASSIFS = [
    "ARAVIS",
    "BAUGES",
    "BEAUFORTAIN",
    "BELLEDONNE",
    "CHABLAIS",
    "CHAMPSAUR",
    "CHARTREUSE",
    "DEVOLUY",
    "EMBRUNAIS-PARPAILLON",
    "GRANDES-ROUSSES",
    "HAUTE-MAURIENNE",
    "HAUTE-TARENTAISE",
    "HAUT-VAR-HAUT-VERDON",
    "MAURIENNE",
    "MERCANTOUR",
    "MONT-BLANC",
    "OISANS",
    "PELVOUX",
    "QUEYRAS",
    "THABOR",
    "UBAYE",
    "VANOISE",
    "VERCORS",
]
NAME_ALIASES = {
    "EMBRUNNAIS-PARPAILLON": "EMBRUNAIS-PARPAILLON",
}


def normalize_name(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "-", value.upper()).strip("-")


def load_massif_metadata() -> dict[str, dict[str, object]]:
    payload = json.loads(META_PATH.read_text(encoding="utf-8"))
    metadata: dict[str, dict[str, object]] = {}
    for feature in payload["features"]:
        props = feature["properties"]
        if not props["mountain"].startswith("Alpes"):
            continue
        key = normalize_name(props["title"])
        metadata[key] = {
            "title": props["title"],
            "lat": float(props["lat_center"]),
            "lon": float(props["lon_center"]),
            "mountain": props["mountain"],
        }

    for metadata_name, panel_name in NAME_ALIASES.items():
        if metadata_name in metadata:
            metadata[panel_name] = metadata[metadata_name]

    return metadata


def label_ssw_windows(panel: pd.DataFrame) -> pd.DataFrame:
    panel = panel.copy()
    panel["ssw_event"] = pd.NaT

    for event_date in SSW_EVENTS:
        start = pd.Timestamp(event_date - timedelta(days=WINDOW_DAYS))
        end = pd.Timestamp(event_date + timedelta(days=WINDOW_DAYS))
        mask = panel["valid_date"].between(start, end)
        panel.loc[mask, "ssw_event"] = pd.Timestamp(event_date)

    panel["is_ssw_window"] = panel["ssw_event"].notna()
    return panel


def build_matched_panel(panel: pd.DataFrame) -> pd.DataFrame:
    panel = label_ssw_windows(panel)
    ssw_panel = panel[panel["is_ssw_window"]].copy()
    control_panel = panel[~panel["is_ssw_window"]].copy()

    control_lookup = (
        control_panel.groupby(["massif", "month_day"]).agg(
            control_mean=("risk_max", "mean"),
            control_high=("high_danger", "mean"),
            control_n=("risk_max", "size"),
        )
    ).reset_index()

    matched = ssw_panel.merge(
        control_lookup,
        on=["massif", "month_day"],
        how="left",
    )
    return matched.dropna(subset=["control_mean", "control_high"]).copy()


def build_pair_panel(matched: pd.DataFrame) -> pd.DataFrame:
    pair_df = (
        matched.groupby(["ssw_event", "massif"]).agg(
            n_days=("valid_date", "nunique"),
            obs_mean=("risk_max", "mean"),
            exp_mean=("control_mean", "mean"),
            obs_high=("high_danger", "mean"),
            exp_high=("control_high", "mean"),
        )
    ).reset_index()
    pair_df["diff"] = pair_df["obs_mean"] - pair_df["exp_mean"]
    pair_df["pct_change"] = pair_df["diff"] / pair_df["exp_mean"]
    pair_df["high_rr"] = pair_df["obs_high"] / pair_df["exp_high"]
    pair_df["decrease"] = pair_df["diff"] < 0
    return pair_df


def bootstrap_mean_diff_ci(pair_df: pd.DataFrame) -> list[float]:
    rng = np.random.default_rng(42)
    diff_samples = np.empty(BOOTSTRAP_ITERATIONS)
    for idx in range(BOOTSTRAP_ITERATIONS):
        sample = pair_df.sample(
            n=len(pair_df),
            replace=True,
            random_state=int(rng.integers(0, 2**31 - 1)),
        )
        diff_samples[idx] = sample["diff"].mean()
    return [float(x) for x in np.percentile(diff_samples, [2.5, 97.5])]


def summarize_group(
    matched: pd.DataFrame,
    pair_df: pd.DataFrame,
) -> dict[str, object]:
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
    obs_mean = matched["risk_max"].mean()
    exp_mean = matched["control_mean"].mean()
    return {
        "n_massifs": int(pair_df["massif"].nunique()),
        "n_event_massif_pairs": n_pairs,
        "obs_mean_danger": float(obs_mean),
        "exp_mean_danger": float(exp_mean),
        "mean_diff": float(obs_mean - exp_mean),
        "mean_pct_change": float((obs_mean - exp_mean) / exp_mean),
        "negative_pairs": n_negative,
        "negative_pair_fraction": float(n_negative / n_pairs),
        "sign_test_p": float(sign_p),
        "wilcoxon_p": float(wilcoxon.pvalue),
        "bootstrap_mean_diff_ci": bootstrap_mean_diff_ci(pair_df),
    }


def corr_stats(x: pd.Series, y: pd.Series) -> dict[str, float]:
    pearson = stats.pearsonr(x, y)
    spearman = stats.spearmanr(x, y)
    return {
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_r": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
    }


def attach_massif_metadata(
    df: pd.DataFrame,
    metadata: dict[str, dict[str, object]],
) -> pd.DataFrame:
    df = df.copy()
    for column in ["mountain", "title", "lat", "lon"]:
        df[column] = df["massif"].map(
            lambda massif, key=column: metadata[massif][key]
        )
    return df


def load_panel(metadata: dict[str, dict[str, object]]) -> pd.DataFrame:
    panel = pd.read_parquet(PANEL_PATH)
    panel["issue_date"] = pd.to_datetime(panel["issue_date"])
    panel["valid_date"] = pd.to_datetime(panel["valid_date"])
    panel["month_day"] = panel["valid_date"].dt.strftime("%m-%d")
    panel["high_danger"] = (
        panel["risk_max"] >= HIGH_DANGER_THRESHOLD
    ).astype(int)
    return attach_massif_metadata(panel, metadata)


def build_output(
    panel: pd.DataFrame,
    matched: pd.DataFrame,
    pair_df: pd.DataFrame,
    massif_summary: pd.DataFrame,
    event_group_summary: pd.DataFrame,
    summaries: dict[str, dict[str, object]],
) -> tuple[dict[str, object], list[str], dict[str, float]]:
    missing_expected_massifs = sorted(
        set(EXPECTED_ALL_MASSIFS) - set(panel["massif"].unique())
    )
    gradient_stats = build_gradient_stats(massif_summary)
    massif_records, event_group_records = build_output_records(
        massif_summary,
        event_group_summary,
    )

    output: dict[str, object] = {
        "panel_path": str(PANEL_PATH.relative_to(ROOT)).replace("\\", "/"),
        "missing_expected_massifs": missing_expected_massifs,
        "group_summaries": summaries,
        "latitude_gradient_mean_diff": gradient_stats["latitude_mean_diff"],
        "latitude_gradient_mean_pct_change": gradient_stats[
            "latitude_mean_pct_change"
        ],
        "longitude_gradient_mean_diff": gradient_stats["longitude_mean_diff"],
        "massif_summary": massif_records,
        "event_group_summary": event_group_records,
        "matched_rows": int(len(matched)),
        "event_massif_pairs": int(len(pair_df)),
    }
    return (
        output,
        missing_expected_massifs,
        gradient_stats["latitude_mean_diff"],
    )


def build_gradient_stats(
    massif_summary: pd.DataFrame,
) -> dict[str, dict[str, float]]:
    return {
        "latitude_mean_diff": corr_stats(
            massif_summary["lat"],
            massif_summary["mean_diff"],
        ),
        "latitude_mean_pct_change": corr_stats(
            massif_summary["lat"],
            massif_summary["mean_pct_change"],
        ),
        "longitude_mean_diff": corr_stats(
            massif_summary["lon"],
            massif_summary["mean_diff"],
        ),
    }


def build_output_records(
    massif_summary: pd.DataFrame,
    event_group_summary: pd.DataFrame,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    massif_records = [
        {
            "massif": row["massif"],
            "title": row["title"],
            "mountain": row["mountain"],
            "lat": float(row["lat"]),
            "lon": float(row["lon"]),
            "n_events": int(row["n_events"]),
            "mean_diff": float(row["mean_diff"]),
            "mean_pct_change": float(row["mean_pct_change"]),
            "mean_high_rr": float(row["mean_high_rr"]),
            "negative_events": int(row["negative_events"]),
        }
        for _, row in massif_summary.iterrows()
    ]
    event_group_records = [
        {
            "ssw_event": row["ssw_event"].date().isoformat(),
            "mountain": row["mountain"],
            "n_massifs": int(row["n_massifs"]),
            "mean_diff": float(row["mean_diff"]),
            "mean_pct_change": float(row["mean_pct_change"]),
            "negative_pairs": int(row["negative_pairs"]),
        }
        for _, row in event_group_summary.iterrows()
    ]
    return massif_records, event_group_records


def build_massif_summary(pair_df: pd.DataFrame) -> pd.DataFrame:
    massif_summary = (
        pair_df.groupby(["massif", "title", "mountain", "lat", "lon"]).agg(
            n_events=("ssw_event", "nunique"),
            mean_diff=("diff", "mean"),
            mean_pct_change=("pct_change", "mean"),
            mean_high_rr=("high_rr", "mean"),
            negative_events=("decrease", "sum"),
        )
    ).reset_index()
    return massif_summary.sort_values("lat", ascending=False).reset_index(
        drop=True
    )


def build_group_summaries(
    matched: pd.DataFrame,
    pair_df: pd.DataFrame,
) -> dict[str, dict[str, object]]:
    summaries = {}
    group_masks = {
        "all": pair_df["massif"].notna(),
        "north": pair_df["mountain"] == "Alpes du Nord",
        "south": pair_df["mountain"] == "Alpes du Sud",
    }
    matched_masks = {
        "all": matched["massif"].notna(),
        "north": matched["mountain"] == "Alpes du Nord",
        "south": matched["mountain"] == "Alpes du Sud",
    }
    for group_name, mask in group_masks.items():
        summaries[group_name] = summarize_group(
            matched.loc[matched_masks[group_name]],
            pair_df.loc[mask],
        )
    return summaries


def build_event_group_summary(pair_df: pd.DataFrame) -> pd.DataFrame:
    return (
        pair_df.groupby(["ssw_event", "mountain"]).agg(
            n_massifs=("massif", "nunique"),
            mean_diff=("diff", "mean"),
            mean_pct_change=("pct_change", "mean"),
            negative_pairs=("decrease", "sum"),
        )
    ).reset_index()


def print_summary(
    summaries: dict[str, dict[str, object]],
    lat_gradient: dict[str, float],
    missing_expected_massifs: list[str],
) -> None:
    print("French BRA spatial gradient summary")
    print("=" * 60)
    for group_name in ["all", "north", "south"]:
        summary = summaries[group_name]
        print(
            f"{group_name:>5}: n_massifs={summary['n_massifs']}, "
            f"diff={summary['mean_diff']:+.3f}, "
            f"negative_pairs={summary['negative_pairs']}/"
            f"{summary['n_event_massif_pairs']}, "
            f"sign_p={summary['sign_test_p']:.4f}, "
            f"wilcoxon_p={summary['wilcoxon_p']:.4f}"
        )

    print(
        f"Latitude gradient (massif mean diff): "
        f"r={lat_gradient['pearson_r']:.3f}, "
        f"p={lat_gradient['pearson_p']:.4g}"
    )
    if missing_expected_massifs:
        print("Missing expected massifs:", ", ".join(missing_expected_massifs))


def run_analysis(
    metadata: dict[str, dict[str, object]],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, dict[str, object]],
]:
    panel = load_panel(metadata)
    matched = build_matched_panel(panel)
    pair_df = attach_massif_metadata(build_pair_panel(matched), metadata)
    massif_summary = build_massif_summary(pair_df)
    summaries = build_group_summaries(matched, pair_df)
    event_group_summary = build_event_group_summary(pair_df)
    return (
        panel,
        matched,
        pair_df,
        massif_summary,
        event_group_summary,
        summaries,
    )


def main() -> None:
    metadata = load_massif_metadata()
    panel, matched, pair_df, massif_summary, event_group_summary, summaries = (
        run_analysis(metadata)
    )

    output, missing_expected_massifs, lat_gradient = build_output(
        panel,
        matched,
        pair_df,
        massif_summary,
        event_group_summary,
        summaries,
    )

    RESULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    with RESULT_OUT.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
    print_summary(summaries, lat_gradient, missing_expected_massifs)
    print(f"Saved results: {RESULT_OUT}")


if __name__ == "__main__":
    main()
