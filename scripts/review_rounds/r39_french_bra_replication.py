"""
Build a French Alpine BRA danger panel for SSW-window month-days and test whether
French avalanche danger forecasts show an independent suppression signal.

The public Meteo-France archive exposes one JSON index per issue date:
    https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA/bra.YYYYMMDD.json

Each JSON entry lists massif names and available timestamps. XML files are then
available at:
    https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA/BRA.{MASSIF}.{TIMESTAMP}.xml

BRA bulletins are next-day forecasts, so issue_date = valid_date - 1 day.
Each bulletin also embeds a 7-day risk history (BSH), which lets one XML file
cover eight consecutive valid dates. This script uses that history to minimize
archive requests while still covering the union of month-days spanned by the
recent SSW windows (+/- 15 d) and matched non-SSW controls.
"""

from __future__ import annotations

import json
import os
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from scipy import stats
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[2]
BASE_URL = "https://donneespubliques.meteofrance.fr/donnees_libres/Pdf/BRA"
GROUP_NAME = os.getenv("BRA_GROUP", "all").strip().lower()

# French Alpine massifs exposed by the public BRA archive.
ALL_ALPINE_MASSIFS = [
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

NORTH_ALPS_MASSIFS = [
    "ARAVIS",
    "BAUGES",
    "BEAUFORTAIN",
    "BELLEDONNE",
    "CHABLAIS",
    "CHARTREUSE",
    "GRANDES-ROUSSES",
    "HAUTE-MAURIENNE",
    "HAUTE-TARENTAISE",
    "MAURIENNE",
    "MONT-BLANC",
    "OISANS",
    "VANOISE",
    "VERCORS",
]

MASSIF_GROUPS = {
    "all": ALL_ALPINE_MASSIFS,
    "north": NORTH_ALPS_MASSIFS,
}
if GROUP_NAME not in MASSIF_GROUPS:
    raise ValueError(f"Unsupported BRA_GROUP '{GROUP_NAME}'. Use one of {sorted(MASSIF_GROUPS)}.")

SELECTED_MASSIFS = MASSIF_GROUPS[GROUP_NAME]
PANEL_OUT = (
    ROOT
    / "data"
    / "processed"
    / "cryosphere"
    / f"french_bra_{GROUP_NAME}_ssw_panel.parquet"
)
RESULT_OUT = ROOT / "data" / "results" / f"r39_french_bra_{GROUP_NAME}_replication.json"

# Existing repo convention for recent Alpine SSW overlap years.
SSW_EVENTS = [
    date(2018, 2, 12),
    date(2019, 1, 1),
    date(2021, 1, 5),
    date(2023, 2, 16),
]

WINDOW_DAYS = 15
ARCHIVE_VALID_START = date(2016, 12, 1)
ARCHIVE_VALID_END = date(2023, 11, 1)
HIGH_DANGER_THRESHOLD = 3
MAX_WORKERS = 6
REQUEST_TIMEOUT = 45
BOOTSTRAP_ITERATIONS = 5000

_THREAD_LOCAL = threading.local()


@dataclass(frozen=True)
class FetchTask:
    issue_date: date
    massif: str
    timestamp: str


def build_session() -> requests.Session:
    retry = Retry(
        total=4,
        backoff_factor=0.8,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=MAX_WORKERS, pool_maxsize=MAX_WORKERS)
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
            month_days.add((event_date + timedelta(days=offset)).strftime("%m-%d"))
    return month_days


def build_target_valid_dates() -> list[date]:
    month_days = event_month_days()
    return [
        current
        for current in daterange(ARCHIVE_VALID_START, ARCHIVE_VALID_END)
        if current.strftime("%m-%d") in month_days
    ]


def issue_date_coverage(issue_date: date, target_valid_dates: set[date]) -> set[date]:
    start = issue_date - timedelta(days=6)
    end = issue_date + timedelta(days=1)
    return {valid_date for valid_date in target_valid_dates if start <= valid_date <= end}


def choose_issue_dates(target_valid_dates: list[date]) -> list[date]:
    uncovered = set(target_valid_dates)
    candidates = list(daterange(min(target_valid_dates) - timedelta(days=1), max(target_valid_dates) + timedelta(days=7)))
    selected: list[date] = []

    while uncovered:
        best_issue_date = max(
            candidates,
            key=lambda issue_date: (len(issue_date_coverage(issue_date, uncovered)), issue_date),
        )
        covered = issue_date_coverage(best_issue_date, uncovered)
        if not covered:
            raise RuntimeError("Could not cover all French BRA target valid dates.")
        selected.append(best_issue_date)
        uncovered -= covered
        candidates.remove(best_issue_date)

    return sorted(selected)


def fetch_json_index(issue_date: date) -> dict[str, list[str]]:
    session = get_session()
    url = f"{BASE_URL}/bra.{issue_date:%Y%m%d}.json"
    response = session.get(url, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    content_type = response.headers.get("content-type", "")
    text = response.content.decode("utf-8-sig")
    if "html" in content_type.lower() or not text.lstrip().startswith("["):
        return {}
    payload = json.loads(text)
    return {entry["massif"]: entry["heures"] for entry in payload}


def choose_bulletin(root: ET.Element) -> ET.Element | None:
    if root.tag == "BULLETINS_NEIGE_AVALANCHE":
        return root

    bulletins = root.findall(".//BULLETINS_NEIGE_AVALANCHE")
    if not bulletins:
        return None

    def sort_key(node: ET.Element) -> tuple[str, str]:
        return (node.attrib.get("DATEDIFFUSION", ""), node.attrib.get("DATEVALIDITE", ""))

    return max(bulletins, key=sort_key)


def build_row(
    task: FetchTask,
    valid_date: date,
    risk_max: int,
    risk_1: int | None,
    risk_2: int | None,
    altitude_split: int | None,
    datediffusion: str | None,
    source: str,
) -> dict[str, object]:
    return {
        "issue_date": task.issue_date.isoformat(),
        "valid_date": valid_date.isoformat(),
        "massif": task.massif,
        "timestamp": task.timestamp,
        "datediffusion": datediffusion,
        "source": source,
        "risk_max": risk_max,
        "risk_1": risk_1,
        "risk_2": risk_2,
        "altitude_split": altitude_split,
    }


def fetch_xml_records(task: FetchTask) -> list[dict[str, object]]:
    session = get_session()
    url = f"{BASE_URL}/BRA.{task.massif}.{task.timestamp}.xml"
    response = session.get(url, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()

    root = ET.fromstring(response.content)
    bulletin = choose_bulletin(root)
    if bulletin is None:
        return []

    risk_node = bulletin.find(".//CARTOUCHERISQUE/RISQUE")
    if risk_node is None or not risk_node.attrib.get("RISQUEMAXI"):
        return []

    datediffusion = bulletin.attrib.get("DATEDIFFUSION")
    valid_date = pd.to_datetime(bulletin.attrib["DATEVALIDITE"]).date()
    rows = [
        build_row(
            task=task,
            valid_date=valid_date,
            risk_max=int(risk_node.attrib["RISQUEMAXI"]),
            risk_1=int(risk_node.attrib["RISQUE1"]) if risk_node.attrib.get("RISQUE1") else None,
            risk_2=int(risk_node.attrib["RISQUE2"]) if risk_node.attrib.get("RISQUE2") else None,
            altitude_split=int(risk_node.attrib["ALTITUDE"]) if risk_node.attrib.get("ALTITUDE") else None,
            datediffusion=datediffusion,
            source="current",
        )
    ]

    for history_node in bulletin.findall(".//BSH//RISQUES//RISQUE"):
        history_date = history_node.attrib.get("DATE")
        history_risk = history_node.attrib.get("RISQUEMAXI")
        if not history_date or not history_risk:
            continue
        rows.append(
            build_row(
                task=task,
                valid_date=pd.to_datetime(history_date).date(),
                risk_max=int(history_risk),
                risk_1=None,
                risk_2=None,
                altitude_split=None,
                datediffusion=datediffusion,
                source="bsh",
            )
        )

    return rows


def fetch_index_task(issue_date: date) -> tuple[date, dict[str, list[str]], str | None]:
    try:
        return issue_date, fetch_json_index(issue_date), None
    except Exception as exc:
        return issue_date, {}, str(exc)


def build_fetch_tasks(issue_dates: list[date]) -> list[FetchTask]:
    tasks: list[FetchTask] = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(fetch_index_task, issue_date): issue_date
            for issue_date in issue_dates
        }
        for idx, future in enumerate(as_completed(futures), start=1):
            issue_date, day_index, error = future.result()
            if error:
                print(f"Index fetch failed for {issue_date}: {error}")
                continue

            for massif in SELECTED_MASSIFS:
                timestamps = day_index.get(massif)
                if not timestamps:
                    continue
                tasks.append(
                    FetchTask(
                        issue_date=issue_date,
                        massif=massif,
                        timestamp=sorted(timestamps)[-1],
                    )
                )

            if idx % 50 == 0 or idx == len(issue_dates):
                print(
                    f"Indexed {idx}/{len(issue_dates)} issue dates; "
                    f"current task count = {len(tasks)}"
                )

    return tasks


def download_panel(tasks: list[FetchTask]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    started = time.time()
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(fetch_xml_records, task): task for task in tasks}
        for idx, future in enumerate(as_completed(futures), start=1):
            task = futures[future]
            try:
                rows.extend(future.result())
            except Exception as exc:
                print(f"XML fetch failed for {task.massif} {task.issue_date} ({task.timestamp}): {exc}")

            if idx % 500 == 0 or idx == len(futures):
                elapsed = time.time() - started
                rate = idx / elapsed if elapsed else 0.0
                print(f"Fetched {idx}/{len(futures)} XML files ({rate:.1f} files/s)")

    panel = pd.DataFrame(rows)
    if panel.empty:
        raise RuntimeError("No French BRA XML rows were parsed.")

    panel["issue_date"] = pd.to_datetime(panel["issue_date"])
    panel["valid_date"] = pd.to_datetime(panel["valid_date"])
    panel = panel.sort_values(["massif", "valid_date", "datediffusion"]).drop_duplicates(
        subset=["massif", "valid_date"],
        keep="last",
    )
    panel["month_day"] = panel["valid_date"].dt.strftime("%m-%d")
    panel["high_danger"] = (panel["risk_max"] >= HIGH_DANGER_THRESHOLD).astype(int)
    return panel.reset_index(drop=True)


def label_ssw_windows(panel: pd.DataFrame) -> pd.DataFrame:
    panel = panel.copy()
    panel["ssw_event"] = pd.NaT
    panel["event_offset_days"] = np.nan

    for event_date in SSW_EVENTS:
        start = pd.Timestamp(event_date - timedelta(days=WINDOW_DAYS))
        end = pd.Timestamp(event_date + timedelta(days=WINDOW_DAYS))
        mask = panel["valid_date"].between(start, end)
        panel.loc[mask, "ssw_event"] = pd.Timestamp(event_date)
        panel.loc[mask, "event_offset_days"] = (panel.loc[mask, "valid_date"] - pd.Timestamp(event_date)).dt.days

    panel["is_ssw_window"] = panel["ssw_event"].notna()
    return panel


def bootstrap_summary(pair_df: pd.DataFrame) -> dict[str, list[float]]:
    rng = np.random.default_rng(42)
    diff_samples = np.empty(BOOTSTRAP_ITERATIONS)
    rr_samples = np.empty(BOOTSTRAP_ITERATIONS)

    for idx in range(BOOTSTRAP_ITERATIONS):
        sample = pair_df.sample(n=len(pair_df), replace=True, random_state=int(rng.integers(0, 2**31 - 1)))
        diff_samples[idx] = sample["diff"].mean()
        rr_samples[idx] = sample["obs_high"].mean() / sample["exp_high"].mean()

    return {
        "diff_ci": [float(x) for x in np.percentile(diff_samples, [2.5, 97.5])],
        "high_rr_ci": [float(x) for x in np.percentile(rr_samples, [2.5, 97.5])],
    }


def analyze_panel(panel: pd.DataFrame) -> dict[str, object]:
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

    matched = ssw_panel.merge(control_lookup, on=["massif", "month_day"], how="left")
    matched = matched.dropna(subset=["control_mean", "control_high"]).copy()
    if matched.empty:
        raise RuntimeError("No matched control rows were found for French BRA SSW windows.")

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
    pair_df["high_diff"] = pair_df["obs_high"] - pair_df["exp_high"]
    pair_df["decrease"] = pair_df["diff"] < 0

    n_pairs = int(len(pair_df))
    n_decrease = int(pair_df["decrease"].sum())
    sign_p = stats.binomtest(n_decrease, n_pairs, 0.5, alternative="greater").pvalue
    wilcoxon = stats.wilcoxon(pair_df["diff"], alternative="less", zero_method="wilcox")
    bootstrap = bootstrap_summary(pair_df)

    event_summary = (
        pair_df.groupby("ssw_event").agg(
            n_massifs=("massif", "nunique"),
            obs_mean=("obs_mean", "mean"),
            exp_mean=("exp_mean", "mean"),
            obs_high=("obs_high", "mean"),
            exp_high=("exp_high", "mean"),
            negative_pairs=("decrease", "sum"),
        )
    ).reset_index()
    event_summary["diff"] = event_summary["obs_mean"] - event_summary["exp_mean"]
    event_summary["high_rr"] = event_summary["obs_high"] / event_summary["exp_high"]

    massif_summary = (
        pair_df.groupby("massif").agg(
            n_events=("ssw_event", "nunique"),
            obs_mean=("obs_mean", "mean"),
            exp_mean=("exp_mean", "mean"),
            negative_events=("decrease", "sum"),
        )
    ).reset_index()
    massif_summary["diff"] = massif_summary["obs_mean"] - massif_summary["exp_mean"]

    overall = {
        "n_events": int(len(SSW_EVENTS)),
        "n_massifs": int(panel["massif"].nunique()),
        "n_panel_rows": int(len(panel)),
        "n_ssw_rows": int(len(matched)),
        "n_event_massif_pairs": n_pairs,
        "coverage_fraction": float(n_pairs / (len(SSW_EVENTS) * len(SELECTED_MASSIFS))),
        "obs_mean_danger": float(matched["risk_max"].mean()),
        "exp_mean_danger": float(matched["control_mean"].mean()),
        "mean_diff": float(matched["risk_max"].mean() - matched["control_mean"].mean()),
        "mean_pct_change": float((matched["risk_max"].mean() - matched["control_mean"].mean()) / matched["control_mean"].mean()),
        "obs_high_danger_rate": float(matched["high_danger"].mean()),
        "exp_high_danger_rate": float(matched["control_high"].mean()),
        "high_danger_rr": float(matched["high_danger"].mean() / matched["control_high"].mean()),
        "negative_pair_fraction": float(n_decrease / n_pairs),
        "sign_test_p": float(sign_p),
        "wilcoxon_statistic": float(wilcoxon.statistic),
        "wilcoxon_p": float(wilcoxon.pvalue),
        "bootstrap_mean_diff_ci": bootstrap["diff_ci"],
        "bootstrap_high_rr_ci": bootstrap["high_rr_ci"],
    }

    return {
        "overall": overall,
        "by_event": [
            {
                "ssw_event": row["ssw_event"].date().isoformat(),
                "n_massifs": int(row["n_massifs"]),
                "obs_mean": float(row["obs_mean"]),
                "exp_mean": float(row["exp_mean"]),
                "diff": float(row["diff"]),
                "obs_high": float(row["obs_high"]),
                "exp_high": float(row["exp_high"]),
                "high_rr": float(row["high_rr"]),
                "negative_pairs": int(row["negative_pairs"]),
            }
            for _, row in event_summary.iterrows()
        ],
        "massif_summary_top_negative": [
            {
                "massif": row["massif"],
                "n_events": int(row["n_events"]),
                "obs_mean": float(row["obs_mean"]),
                "exp_mean": float(row["exp_mean"]),
                "diff": float(row["diff"]),
                "negative_events": int(row["negative_events"]),
            }
            for _, row in massif_summary.sort_values("diff").head(10).iterrows()
        ],
    }


def main() -> None:
    PANEL_OUT.parent.mkdir(parents=True, exist_ok=True)
    RESULT_OUT.parent.mkdir(parents=True, exist_ok=True)

    valid_dates = build_target_valid_dates()
    issue_dates = choose_issue_dates(valid_dates)

    print(f"Target valid dates: {len(valid_dates)}")
    print(f"Target issue dates: {len(issue_dates)}")
    print(f"Event month-days: {len(event_month_days())}")
    print(f"Massif group: {GROUP_NAME} ({len(SELECTED_MASSIFS)} massifs)")

    tasks = build_fetch_tasks(issue_dates)
    print(f"Total XML fetch tasks: {len(tasks)}")

    panel = download_panel(tasks)
    target_valid_timestamps = pd.to_datetime(sorted(valid_dates))
    panel = panel[panel["valid_date"].isin(target_valid_timestamps)].copy()
    panel.to_parquet(PANEL_OUT, index=False)
    print(f"Saved panel: {PANEL_OUT}")

    results = analyze_panel(panel)
    with RESULT_OUT.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print(f"Saved results: {RESULT_OUT}")

    overall = results["overall"]
    print("\nFrench BRA replication summary")
    print("=" * 60)
    print(f"Panel rows: {overall['n_panel_rows']}")
    print(f"Massifs: {overall['n_massifs']}")
    print(f"Event-massif pairs: {overall['n_event_massif_pairs']}")
    print(f"Observed mean danger: {overall['obs_mean_danger']:.3f}")
    print(f"Expected mean danger: {overall['exp_mean_danger']:.3f}")
    print(f"Mean difference: {overall['mean_diff']:+.3f} ({overall['mean_pct_change']:+.1%})")
    print(f"High-danger RR: {overall['high_danger_rr']:.3f}")
    print(f"Negative pairs: {overall['negative_pair_fraction']:.1%}")
    print(f"Sign-test p: {overall['sign_test_p']:.4f}")
    print(f"Wilcoxon p: {overall['wilcoxon_p']:.4f}")


if __name__ == "__main__":
    main()