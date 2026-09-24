#!/usr/bin/env python3
"""
smoke_test.py
=============
The repository's standing checks. Every one of them exists because something
here once went wrong in that exact way -- see FAILURES.md, which this file is
the executable half of.

    .venv/bin/python ssw-design-analysis/tests/smoke_test.py
    .venv/bin/python ssw-design-analysis/tests/smoke_test.py --quick   # skip the rebuild

Exit 0 only if every check passes. A check that cannot run (missing optional
input) reports SKIP and does not fail the suite; a check that runs and finds a
defect reports FAIL.

WHAT IT DOES NOT DO
  It does not re-run the analyses. It verifies the frozen inputs, the shape of
  the repository, and the machine-checkable half of the project's own rules.
  Scientific verification is re-running the script and reading the table.
"""
from __future__ import annotations

import ast
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SSW = ROOT / "ssw-design-analysis"
RESULTS = ROOT / "results" / "current"

rows: list[tuple[str, str, str]] = []


def rec(status: str, name: str, detail: str = "") -> None:
    rows.append((status, name, detail))


def live_scripts() -> list[Path]:
    """Analysis code only: stage directories, excluding archived trees."""
    out = []
    for stage in sorted(SSW.glob("[0-9][0-9]_*")):
        out.extend(sorted(stage.glob("*.py")))
    return out


# ---------------------------------------------------------------- layout ----

def check_layout() -> None:
    stages = sorted(p.name for p in SSW.glob("[0-9][0-9]_*") if p.is_dir())
    nums = [s[:2] for s in stages]
    dupes = {n for n in nums if nums.count(n) > 1}
    if dupes:
        rec("FAIL", "stage numbering", f"duplicate stage numbers: {sorted(dupes)}")
    else:
        rec("PASS", "stage numbering", f"{len(stages)} stages, no collisions")

    strays = [p for stage in SSW.glob("[0-9][0-9]_*") for p in stage.glob("*.log")]
    if strays:
        rec("FAIL", "run logs", f"{len(strays)} .log files still in stage dirs "
                                f"(they belong in run_logs/<stage>/)")
    else:
        rec("PASS", "run logs", f"{len(list((SSW / 'run_logs').rglob('*.log')))} "
                                f"logs, all under run_logs/")

    # The SNAPSI reduced cache is read BY FILENAME by snapsi_selection_test.py
    # and snapsi_distribution_test.py: the centre is p.name.split("_")[0] and
    # members load via {centre}_{exp}_{init}_*.parquet. On 2026-09-17 a
    # variable prefix was added to those keys and both scripts silently parsed
    # "psl" as a centre, analysing 2,026 of 6,161 members without any error.
    red = SSW / "03_data_ingestion" / "_snapsi_reduced"
    if red.is_dir():
        names = [p.name for p in red.glob("*.parquet")]
        known = {"CCCma", "CNR-ISAC", "ECCC", "ECMWF", "KMA", "Meteo-France",
                 "NCAR", "NRL", "SNU", "UKMO"}
        bad = sorted({n.split("_")[0] for n in names} - known)
        rec("PASS" if not bad else "FAIL", "snapsi cache naming",
            f"{len(names)} members, all parse to a known centre" if not bad
            else f"filenames parse to non-centres {bad} -- the analysis scripts "
                 f"read this directory by filename")
    else:
        rec("SKIP", "snapsi cache naming", "cache absent")

    # Lead 0 must mean 00 UTC on the init date for every member, because every
    # analysis converts lead to post-onset day by that assumption. UKMO files
    # start at 06 UTC and sat 6 h early until the caches were rebased
    # (acquire_snapsi_surface.rebase_cache, 2026-09-24). Schema-only read: fast.
    import pyarrow.parquet as pq
    for label, d, pat in (("psl", red, "*.parquet"),
                          ("zg100", SSW / "03_data_ingestion" / "_snapsi_zg100",
                           "zg100_*.parquet")):
        if not d.is_dir():
            rec("SKIP", f"snapsi {label} lead origin", "cache absent")
            continue
        files = sorted(d.glob(pat))
        unmarked = [f.name for f in files
                    if "lead_origin" not in pq.read_schema(f).names]
        rec("PASS" if not unmarked else "FAIL", f"snapsi {label} lead origin",
            f"{len(files)} members rebased to 00 UTC on the init date" if not unmarked
            else f"{len(unmarked)} of {len(files)} not rebased (run the acquisition "
                 f"script with --rebase), e.g. {unmarked[0]}")

    for required in ("FAILURES.md", "CLAUDE.md", "CATALOG.md"):
        rec("PASS" if (ROOT / required).exists() else "FAIL",
            f"root file {required}",
            "present" if (ROOT / required).exists() else "MISSING")


# ----------------------------------------------------------- environment ----

def check_environment() -> None:
    script = SSW / "environment" / "check_environment.py"
    if not script.exists():
        rec("FAIL", "environment", "check_environment.py missing")
        return
    proc = subprocess.run([sys.executable, str(script), "--no-remote"],
                          capture_output=True, text=True, timeout=300)
    last = [l for l in proc.stdout.splitlines() if l.strip()]
    rec("PASS" if proc.returncode == 0 else "FAIL", "environment",
        last[-1] if last else f"exit {proc.returncode}")


# ------------------------------------------------------------- catalogue ----

def check_catalogue(quick: bool) -> None:
    cat_dir = SSW / "02_event_catalogues"
    csv = cat_dir / "event_catalogue.csv"
    sums = cat_dir / "catalogue_checksum.txt"
    if not (csv.exists() and sums.exists()):
        rec("SKIP", "catalogue", "event_catalogue.csv or checksum absent")
        return

    import hashlib
    recorded = dict(l.split(maxsplit=1) for l in sums.read_text().splitlines() if l.strip())
    actual = hashlib.sha256(csv.read_bytes()).hexdigest()
    expect = recorded.get("event_catalogue_csv_sha256", "").strip()
    rec("PASS" if actual == expect else "FAIL", "catalogue checksum",
        "matches catalogue_checksum.txt" if actual == expect
        else f"on disk {actual[:16]} vs recorded {expect[:16]}")

    # CRLF is what broke this checksum on 2026-09-17.
    crlf = csv.read_bytes().count(b"\r\n")
    rec("PASS" if crlf == 0 else "FAIL", "catalogue line endings",
        "LF only" if crlf == 0 else f"{crlf} CRLF lines -- checksum is platform-dependent")

    sys.path.insert(0, str(cat_dir))
    try:
        import build_catalogue as bc
        prim = bc.load_catalogue("primary")
        n_ev, n_wi = len(prim), len({bc.winter_of(d) for d in prim})
        ok = (n_ev == 43 and n_wi == 36)
        rec("PASS" if ok else "FAIL", "catalogue frozen counts",
            f"{n_ev} events / {n_wi} winters" + ("" if ok else " -- expected 43 / 36"))
    except Exception as exc:  # noqa: BLE001
        rec("FAIL", "catalogue load", f"{type(exc).__name__}: {exc}")

    if quick:
        rec("SKIP", "catalogue rebuild", "--quick")
        return
    before = csv.read_bytes()
    proc = subprocess.run([sys.executable, str(cat_dir / "build_catalogue.py")],
                          capture_output=True, text=True, cwd=cat_dir, timeout=600)
    after = csv.read_bytes()
    if proc.returncode != 0:
        rec("FAIL", "catalogue rebuild", (proc.stderr or proc.stdout).strip()[-200:])
    elif before != after:
        rec("FAIL", "catalogue rebuild", "rebuild is NOT byte-identical")
    elif "all passed" not in proc.stdout:
        rec("FAIL", "catalogue invariance", "spec §9 invariance tests did not all pass")
    else:
        rec("PASS", "catalogue rebuild", "byte-identical, §9 invariance tests pass")


# --------------------------------------------------------------- results ----

def check_results() -> None:
    files = sorted(RESULTS.rglob("*.json"))
    if not files:
        rec("SKIP", "results readable", "no results found")
        return
    bad = []
    payloads = {}
    for f in files:
        try:
            payloads[f] = json.loads(f.read_text())
        except Exception as exc:  # noqa: BLE001
            bad.append(f"{f.name}: {type(exc).__name__}")
    rec("PASS" if not bad else "FAIL", "results readable",
        f"{len(files)} JSON files parse" if not bad else "; ".join(bad[:4]))

    # An R^2 above 1 is impossible and is how an unmatched variance denominator
    # announced itself on 2026-08-03. Match only quantities that really are a
    # variance share: `implied_DW_sigma` is an effect size in sigma and is
    # legitimately > 1, so keying on the word "implied" alone is a false alarm.
    r2_key = re.compile(r"(^|[/_])(r2|r\^2|rsq|variance_explained|implied)$"
                        r"|implied_r2|_r2$|cv_r2", re.I)
    offenders = []
    for f, payload in payloads.items():
        for path, val in walk(payload):
            leaf = path.rsplit("/", 1)[-1]
            container = path.rsplit("/", 2)[1] if path.count("/") >= 2 else ""
            is_share = r2_key.search(leaf) or r2_key.search(container)
            if is_share and not leaf.lower().endswith(("_sigma", "_pa", "_days")):
                if isinstance(val, (int, float)) and not isinstance(val, bool) and val > 1.0:
                    offenders.append(f"{f.relative_to(ROOT)}{path} = {val}")
    rec("PASS" if not offenders else "FAIL", "no impossible R^2",
        "none found" if not offenders else "; ".join(offenders[:4]))

    # N_BOOT >= 1000 is a standing rule (FAILURES.md, 2026-07-28).
    thin = []
    for f, payload in payloads.items():
        for path, val in walk(payload):
            if path.lower().endswith(("n_boot", "n_bootstrap")) and isinstance(val, int):
                if val < 1000:
                    thin.append(f"{f.relative_to(ROOT)}{path} = {val}")
    rec("PASS" if not thin else "FAIL", "bootstrap replicate floor",
        "all n_boot >= 1000" if not thin else "; ".join(thin[:4]))

    # Informational, never a failure: every artifact written on Windows is CRLF
    # and a Linux re-run produces LF with identical values (FAILURES.md,
    # 2026-09-17). Reported so nobody reads a byte difference as a result change.
    crlf = [f for f in files if b"\r\n" in f.read_bytes()]
    rec("INFO", "result line endings",
        f"{len(crlf)}/{len(files)} results are CRLF (Windows-written) -- "
        f"compare parsed, not byte-wise")


def walk(obj, prefix=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk(v, f"{prefix}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from walk(v, f"{prefix}[{i}]")
    else:
        yield prefix, obj


# ----------------------------------------------------------- source rules ----

def check_source_rules() -> None:
    scripts = live_scripts()

    # hash(str) seeding is per-process randomised (FAILURES.md, 2026-07-28).
    # Parsed, not grepped: every one of these scripts *documents* the old bug in
    # a comment, so a regex over the source reports the fix as the defect.
    bad_seed = []
    for p in scripts:
        try:
            tree = ast.parse(p.read_text(encoding="utf8", errors="replace"))
        except SyntaxError:
            continue  # reported separately below
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name) and node.func.id == "hash"
                    and node.args
                    and isinstance(node.args[0], ast.Constant)
                    and isinstance(node.args[0].value, str)):
                bad_seed.append(f"{p.relative_to(ROOT).as_posix()}:{node.lineno}")
    rec("PASS" if not bad_seed else "FAIL", "no hash(str) seeding",
        f"{len(scripts)} scripts clean" if not bad_seed else "; ".join(bad_seed[:4]))

    # to_csv without an explicit lineterminator is platform-dependent when the
    # output is checksummed (FAILURES.md, 2026-09-17).
    bad_csv = []
    for p in scripts:
        src = p.read_text(encoding="utf8", errors="replace")
        for m in re.finditer(r"\.to_csv\(([^)]*)\)", src, flags=re.S):
            if "lineterminator" not in m.group(1):
                bad_csv.append(f"{p.relative_to(ROOT).as_posix()}:"
                               f"{src[:m.start()].count(chr(10)) + 1}")
    rec("PASS" if not bad_csv else "WARN", "to_csv pins lineterminator",
        "all pinned" if not bad_csv else "; ".join(bad_csv[:4]))

    # Every live script should parse -- a syntax error in a stage dir breaks
    # the sibling imports that 147 sites depend on.
    broken = []
    for p in scripts:
        try:
            ast.parse(p.read_text(encoding="utf8", errors="replace"))
        except SyntaxError as exc:
            broken.append(f"{p.relative_to(ROOT).as_posix()}:{exc.lineno}")
    rec("PASS" if not broken else "FAIL", "all live scripts parse",
        f"{len(scripts)} scripts" if not broken else "; ".join(broken[:4]))


# --------------------------------------------------------------- catalog ----

def check_catalog_tool() -> None:
    cj = ROOT / "CATALOG.json"
    if not cj.exists():
        rec("SKIP", "catalogue tool", "CATALOG.json absent -- run tools/catalog_repo.py")
        return
    data = json.loads(cj.read_text())
    results = data.get("results", [])
    if not results:
        rec("SKIP", "catalogue staleness", "CATALOG.json lists no results")
        return
    stale = [r["path"] for r in results if r.get("state") == "stale"]
    orphan = [r["path"] for r in results if not r.get("producers")]
    unreadable = [r["path"] for r in results if r.get("readable") is False]
    rec("PASS" if not stale else "FAIL", "no stale results",
        f"{len(results)} catalogued, {len(stale)} stale, {len(orphan)} orphan, "
        f"{len(unreadable)} unreadable"
        + ("" if not stale else f" -- {stale[:3]}"))


def main() -> int:
    quick = "--quick" in sys.argv
    check_layout()
    check_environment()
    check_catalogue(quick)
    check_results()
    check_source_rules()
    check_catalog_tool()

    width = max(len(n) for _, n, _ in rows)
    for status, name, detail in rows:
        print(f"[{status:4s}] {name:<{width}}  {detail}")

    n_fail = sum(1 for s, _, _ in rows if s == "FAIL")
    n_warn = sum(1 for s, _, _ in rows if s == "WARN")
    n_skip = sum(1 for s, _, _ in rows if s == "SKIP")
    print(f"\n{len(rows)} checks: {n_fail} FAIL, {n_warn} WARN, {n_skip} SKIP")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
