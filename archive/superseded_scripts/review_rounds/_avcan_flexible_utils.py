from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd  # type: ignore[import-untyped]

EARTH_RADIUS_KM = 6371.0088
DANGER_MAP = {
    "no_rating": 0,
    "norating": 0,
    "no_snow": 0,
    "earlyseason": 0,
    "low": 1,
    "moderate": 2,
    "considerable": 3,
    "high": 4,
    "extreme": 5,
    "very_high": 5,
}


def ring_area_sqkm(ring: list[list[float]]) -> float:
    if len(ring) < 4:
        return 0.0

    coords = np.asarray(ring, dtype=float)
    lon = np.deg2rad(coords[:, 0])
    lat = np.deg2rad(coords[:, 1])
    lat0 = float(lat.mean())
    x = EARTH_RADIUS_KM * lon * math.cos(lat0)
    y = EARTH_RADIUS_KM * lat
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def geometry_area_sqkm(geometry: Any) -> float:
    if not geometry:
        return 0.0

    geometry_type = geometry.get("type")
    coordinates = geometry.get("coordinates")
    if geometry_type == "Polygon":
        polygons = [coordinates]
    elif geometry_type == "MultiPolygon":
        polygons = coordinates
    else:
        return 0.0

    total_area = 0.0
    for polygon in polygons:
        if not polygon:
            continue
        outer = ring_area_sqkm(polygon[0])
        holes = sum(ring_area_sqkm(ring) for ring in polygon[1:])
        total_area += max(0.0, outer - holes)
    return total_area


def first_nonempty(*values: Any) -> Any:
    for value in values:
        if value:
            return value
    return ""


def parse_danger_value(entry: dict[str, Any]) -> int:
    highest = entry.get("highestDanger") or {}
    value = str(highest.get("value", "")).strip().lower()
    if value in DANGER_MAP:
        return DANGER_MAP[value]

    icons = entry.get("icons") or []
    for icon in icons:
        ratings = icon.get("ratings") or {}
        values = [
            DANGER_MAP.get(str(v).strip().lower(), 0)
            for v in ratings.values()
        ]
        if values:
            return int(max(values))
    return 0


def build_area_lookup(areas: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        feature["id"]: {
            "area_sqkm": geometry_area_sqkm(feature.get("geometry")),
            "bbox": feature.get("bbox"),
        }
        for feature in areas.get("features", [])
    }


def build_panel_row(
    target_date: Any,
    dt: str,
    entry: dict[str, Any],
    area_lookup: dict[str, dict[str, Any]],
    high_danger_threshold: int,
) -> dict[str, Any]:
    area = entry.get("area") or {}
    centroid = entry.get("centroid") or {}
    owner = entry.get("owner") or {}
    region_id = str(area.get("id") or "")
    region_name = str(
        first_nonempty(
            area.get("name"),
            entry.get("product", {}).get("title"),
        )
    )
    area_meta = area_lookup.get(region_id, {})
    danger_value = parse_danger_value(entry)

    return {
        "date": pd.Timestamp(target_date),
        "month_day": target_date.strftime("%m-%d"),
        "archive_datetime": dt,
        "region_id": region_id,
        "region_name": region_name,
        "owner": str(first_nonempty(owner.get("display"), owner.get("value"))),
        "centroid_lat": centroid.get("latitude"),
        "centroid_lon": centroid.get("longitude"),
        "bbox": area_meta.get("bbox"),
        "area_sqkm": float(area_meta.get("area_sqkm") or 0.0),
        "danger_max": int(danger_value),
        "high_danger": int(danger_value >= high_danger_threshold),
    }


def serialize_event_summary(
    event_summary: pd.DataFrame,
) -> list[dict[str, Any]]:
    return [
        {
            "ssw_event": row["ssw_event"].date().isoformat(),
            "n_days": int(row["n_days"]),
            "obs_mean": float(row["obs_mean"]),
            "exp_mean": float(row["exp_mean"]),
            "diff": float(row["diff"]),
            "obs_high": float(row["obs_high"]),
            "exp_high": float(row["exp_high"]),
            "high_diff": float(row["high_diff"]),
            "negative_days": int(row["negative_days"]),
        }
        for _, row in event_summary.iterrows()
    ]


def serialize_offset_summary(
    offset_summary: pd.DataFrame,
) -> list[dict[str, Any]]:
    return [
        {
            "event_offset_days": int(row["event_offset_days"]),
            "n_rows": int(row["n_rows"]),
            "obs_mean": float(row["obs_mean"]),
            "exp_mean": float(row["exp_mean"]),
            "diff": float(row["diff"]),
            "obs_high": float(row["obs_high"]),
            "exp_high": float(row["exp_high"]),
            "high_diff": float(row["high_diff"]),
        }
        for _, row in offset_summary.iterrows()
    ]
