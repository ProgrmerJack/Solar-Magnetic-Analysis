#!/usr/bin/env python3
"""
catalog_repo.py
===============
The single cataloguing tool for this repository. Run it after any analysis and
it rebuilds every index from the files themselves -- nothing here is maintained
by hand, so nothing here can silently go stale.

WHY THIS EXISTS
  Nobody -- human or model -- can hold ~250 results and ~165 scripts in their
  head, and the failure mode is real: claims were made in this project from
  documents that had gone stale against their own JSON.

LAYOUT IT ASSUMES (unified 2026-08-04)
  results/current/<theme>/       53 live results, grouped by the question they
                                 answer (1_catalogue ... 9_literature)
  results/superseded/<family>/   193 retired results, grouped by family
  results/_ALL_RESULTS.json      all 53 live payloads in one file
  results/_PROVENANCE.json       producer + sha256 baseline

  Scripts stay in their stage folders: 147 sites import them as modules by
  stage path, and the numbering encodes pipeline order. Only outputs moved.

WHAT IT EMITS
  CATALOG.json               every script and every result, machine-readable
  CATALOG.md                 the human index, grouped and sorted
  results/_ALL_RESULTS.json  the single-file read path for the numbers

HOW PROVENANCE IS RESOLVED
  Not by filename convention -- by grep. Every script is scanned for the literal
  basename of every result file, so a result whose producer was renamed is still
  attributed, and a result nothing writes is flagged ORPHAN rather than assumed.

STALENESS -- BY CONTENT HASH, NOT MTIME
  A result is STALE when its producing script has CHANGED since the result was
  written. The baseline lives in results/_PROVENANCE.json: producer path plus
  the sha256 the producer had when the result was last known good.

  mtime was tried first and does not work. Any bulk edit -- the 2026-08-03
  migration that repointed 50 write sites into results/ -- rewrites every
  script's mtime and marks every result stale, which is both false and useless.
  A content hash is immune to that: moving or reformatting the repo does not
  change a script's bytes, and editing its logic does.

  Re-record the baseline with `--record` AFTER re-running an analysis, never
  before; recording is an assertion that the JSON on disk came from the script
  now on disk.
"""
from __future__ import annotations

import ast
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "CATALOG.json"
OUT_MD = ROOT / "CATALOG.md"
OUT_RESULTS = ROOT / "results" / "_ALL_RESULTS.json"

# Trees that hold analysis code/results. Anything not listed is data, not code.
SCRIPT_ROOTS = ["ssw-design-analysis", "scripts", "tools"]
# ONE results tree for the whole repo (unified 2026-08-04):
#   results/current/<theme>      the live ssw-design-analysis outputs
#   results/superseded/<family>  the retired project, kept for provenance
RESULT_ROOTS = ["results/current"]
SUPERSEDED_ROOT = "results/superseded"

SKIP_DIR_PARTS = {"__pycache__", ".git", ".mypy_cache", "node_modules", ".codacy"}
# raw/ holds downloaded API payloads, not computed results
SKIP_RESULT_PARTS = {"raw"}

TOUCH_TOLERANCE_MIN = 15   # retained: referenced by nothing after the hash switch

VERDICT_KEYS = ("verdict", "summary", "conclusion", "interpretation", "reading",
                "note", "status", "headline_recommendation", "interpretation_note")
COUNT_KEYS = ("n_events", "n_unique_winters", "n_winters", "n_boot", "n_draw",
              "n_draws", "n_seeds", "n_rep", "n_perm", "n_members", "n_null")


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def walk(rel_roots, suffix):
    for rel in rel_roots:
        base = ROOT / rel
        if not base.exists():
            continue
        for p in sorted(base.rglob(f"*{suffix}")):
            if SKIP_DIR_PARTS & set(p.parts):
                continue
            yield p


def stamp(p: Path) -> str:
    return datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).strftime(
        "%Y-%m-%d %H:%M")


# ---------------------------------------------------------------- scripts ---

def script_title(src: str) -> str:
    """First real sentence of the module docstring, ignoring the ==== rule."""
    try:
        doc = ast.get_docstring(ast.parse(src)) or ""
    except SyntaxError:
        return "(unparseable)"
    lines = [ln.strip() for ln in doc.splitlines()]
    lines = [ln for ln in lines if ln and not set(ln) <= {"=", "-", "~"}]
    if not lines:
        return ""
    # drop a leading bare filename line
    if lines[0].endswith(".py") and len(lines) > 1:
        lines = lines[1:]
    out = []
    for ln in lines:
        out.append(ln)
        if ln.endswith("."):
            break
        if len(" ".join(out)) > 200:
            break
    return " ".join(out).strip()


# A script that has been retired in place says so in its own docstring. Those
# markers were invisible to this tool, so `ao_null_perbin.py` -- whose header
# reads "VOID. DO NOT CITE." -- had its result sitting in results/current/ with
# nothing to distinguish it from a live one. Match only the module docstring:
# the body of a *correct* script often discusses the void one by name.
SCRIPT_FLAG_RE = re.compile(
    r"\b(VOID|DO NOT CITE|DO NOT USE|RETRACTED|SUPERSEDED BY|"
    r"FLAWED|DEPRECATED|DO NOT RUN)\b", re.I)


def script_flags(src: str) -> list[str]:
    """Retirement markers declared in a script's own module docstring."""
    try:
        doc = ast.get_docstring(ast.parse(src)) or ""
    except SyntaxError:
        return []
    # Only the opening block: a long docstring may narrate history further down.
    head = "\n".join(doc.splitlines()[:12])
    return sorted({m.group(1).upper() for m in SCRIPT_FLAG_RE.finditer(head)})


EXT = r"(?:json|csv|parquet|txt|nc)"
# a plain literal: "canonical_event_study.json"
WRITE_RE = re.compile(rf"""["']([A-Za-z0-9_.\-]+\.{EXT})["']""")
# an f-string with an interpolated part: f"ensemble_precursor_{model}.json"
# -> record the literal PREFIX so the output can still be attributed.
FSTR_RE = re.compile(rf"""["']([A-Za-z0-9_.\-]{{4,}})\{{[^"']*\.{EXT}["']""")
MIN_PREFIX = 6


SELF = Path(__file__).resolve()


def scan_scripts():
    rows = []
    for p in walk(SCRIPT_ROOTS, ".py"):
        try:
            src = p.read_text(encoding="utf8", errors="replace")
        except OSError:
            continue
        rel = p.relative_to(ROOT).as_posix()
        if p.resolve() == SELF:
            # this tool's own regex source contains filename-shaped literals;
            # attributing results to it would be pure noise
            rows.append({
                "path": rel, "tree": rel.split("/")[0], "stage": "",
                "title": script_title(src), "loc": src.count("\n") + 1,
                "modified": stamp(p), "sha256_16": sha(p),
                "flags": [],
                "filenames_mentioned": [], "filename_prefixes": [],
            })
            continue
        rows.append({
            "path": rel,
            "tree": rel.split("/")[0],
            "stage": rel.split("/")[1] if rel.count("/") > 1 else "",
            "title": script_title(src),
            "loc": src.count("\n") + 1,
            "modified": stamp(p),
            "sha256_16": sha(p),
            "flags": script_flags(src),
            "filenames_mentioned": sorted(set(WRITE_RE.findall(src))),
            "filename_prefixes": sorted({m for m in FSTR_RE.findall(src)
                                         if len(m) >= MIN_PREFIX}),
        })
    return rows


# ---------------------------------------------------------------- results ---

def flatten(o, prefix="", depth=0, out=None, maxdepth=3):
    if out is None:
        out = []
    if depth > maxdepth:
        return out
    if isinstance(o, dict):
        for k, v in o.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            if isinstance(v, (dict, list)):
                flatten(v, key, depth + 1, out, maxdepth)
            else:
                out.append((key, v))
    elif isinstance(o, list):
        for i, v in enumerate(o[:4]):
            flatten(v, f"{prefix}[{i}]", depth + 1, out, maxdepth)
    return out


def summarise_result(payload):
    """Verdicts and the headline counts, without dumping the whole tree."""
    flat = flatten(payload)
    verdicts, counts = [], {}
    for k, v in flat:
        leaf = k.split(".")[-1].lower()
        if isinstance(v, str) and any(t in leaf for t in VERDICT_KEYS) and len(v) > 3:
            verdicts.append({"key": k, "text": v})
        if leaf in COUNT_KEYS and isinstance(v, (int, float)):
            counts.setdefault(leaf, sorted(set()) or [])
            if v not in counts[leaf]:
                counts[leaf].append(v)
    return verdicts, {k: sorted(v) for k, v in counts.items()}


PROV = ROOT / "results" / "_PROVENANCE.json"


def load_provenance():
    if PROV.exists():
        return json.loads(PROV.read_text(encoding="utf8")).get("results", {})
    return {}


def record_provenance(results):
    """Stamp the current producer hashes as the known-good baseline."""
    out = {}
    for r in results:
        if not r.get("readable", True) or not r.get("producers"):
            continue
        out[Path(r["path"]).name] = {
            "producers": r["producers"],
            "producer_sha256": r["producer_sha256"],
            "result_sha256_16": r["sha256_16"],
            "needs_rerun": False,
        }
    PROV.write_text(json.dumps(
        {"generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "note": ("Known-good baseline: each result was produced by the script "
                  "at the recorded sha256. catalog_repo.py flags STALE when a "
                  "producer's hash no longer matches. Re-record with --record "
                  "only AFTER re-running the analysis."),
         "results": out}, indent=1), encoding="utf8")
    return len(out)


READER_ONLY_DIRS = {"09_figures", "10_tables"}


def scan_results(scripts):
    """Attribute each result to whichever script mentions its basename."""
    prov = load_provenance()
    flag_by_script = {s["path"]: s.get("flags", []) for s in scripts}
    # Display-item scripts READ results and never write them. Attributing by
    # mention made every figure script a "producer" of the JSON it plots, so
    # adding a figure marked K, L and predictability_ceiling stale (2026-09-24).
    readers = [s for s in scripts
               if set(Path(s["path"]).parts) & READER_ONLY_DIRS]
    scripts = [s for s in scripts if s not in readers]
    mention = {}
    for s in scripts:
        for fn in s["filenames_mentioned"]:
            mention.setdefault(fn, []).append(s["path"])
    prefixes = [(pre, s["path"]) for s in scripts
                for pre in s.get("filename_prefixes", [])]

    def producers_for(name: str):
        hits = list(mention.get(name, []))
        for pre, path in prefixes:
            if name.startswith(pre) and path not in hits:
                hits.append(path)
        return sorted(set(hits))

    rows, payloads = [], {}
    for p in walk(RESULT_ROOTS, ".json"):
        if SKIP_RESULT_PARTS & set(p.parts):
            continue
        rel = p.relative_to(ROOT).as_posix()
        if p.name.startswith("_"):        # _ALL_RESULTS.json, _PROVENANCE.json
            continue
        try:
            payload = json.loads(p.read_text(encoding="utf8"))
        except Exception as e:
            rows.append({"path": rel, "readable": False, "error": str(e)[:120]})
            continue

        producers = producers_for(p.name)
        rec = prov.get(p.name, {})
        cur_sha = {q: sha(ROOT / q) for q in producers}
        flags = sorted({f for q in producers for f in flag_by_script.get(q, [])})
        if flags:
            # A result whose own producer is marked VOID / FLAWED / SUPERSEDED
            # outranks every other state: nothing downstream should cite it.
            state = "FLAGGED"
            detail = ("producer declares itself "
                      + ", ".join(flags) + f": {', '.join(producers)}")
        elif not producers:
            state, detail = "ORPHAN", "no script writes this"
        elif rec.get("needs_rerun"):
            state, detail = "STALE", rec.get("note", "flagged for re-run")
        elif not rec:
            state, detail = "unverified", "no provenance baseline recorded"
        elif rec.get("producer_sha256") != cur_sha:
            changed = [q for q in producers
                       if rec.get("producer_sha256", {}).get(q) != cur_sha[q]]
            state, detail = "STALE", f"producer changed: {', '.join(changed)}"
        else:
            state, detail = "current", ""

        verdicts, counts = summarise_result(payload)
        rows.append({
            "path": rel,
            "tree": "current",
            "stage": Path(rel).parent.name,
            "readable": True,
            "bytes": p.stat().st_size,
            "modified": stamp(p),
            "sha256_16": sha(p),
            "producers": producers,
            "state": state,
            "detail": detail,
            "flags": flags,
            "producer_sha256": cur_sha,
            "verdicts": verdicts,
            "counts": counts,
            "top_level_keys": list(payload)[:40] if isinstance(payload, dict)
                              else f"[list of {len(payload)}]",
        })
        payloads[Path(rel).name] = payload
    return rows, payloads


# ------------------------------------------------------------------ emit ---

def write_md(scripts, results):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    L = [f"# Repository catalogue",
         "",
         f"Generated by `tools/catalog_repo.py` at {now}. **Do not edit by hand** "
         f"— re-run the tool.",
         ""]

    stale = [r for r in results if r.get("state") == "STALE"]
    orphan = [r for r in results if r.get("state") == "ORPHAN"]
    unread = [r for r in results if not r.get("readable", True)]
    unver = sum(r.get("state") == "unverified" for r in results)

    L += ["## Health", "",
          f"| check | count |", "|---|---|",
          f"| result files catalogued | {len(results)} |",
          f"| scripts catalogued | {len(scripts)} |",
          f"| **stale** (producer changed since the result was written) | "
          f"**{len(stale)}** |",
          f"| **orphan** (no script writes it) | **{len(orphan)}** |",
          f"| unverified (no provenance baseline) | {unver} |",
          f"| unreadable | {len(unread)} |", ""]

    if stale:
        L += ["### Stale results — re-run these", "",
              "| result | producer | why |", "|---|---|---|"]
        for r in sorted(stale, key=lambda r: r["path"]):
            L.append(f"| `{Path(r['path']).name}` | "
                     f"{', '.join(f'`{Path(q).name}`' for q in r['producers'])} | "
                     f"{r.get('detail','')} |")
        L.append("")

    if orphan:
        L += ["### Orphan results — no script in the repo writes these", "",
              "| result | bytes | modified |", "|---|---|---|"]
        for r in sorted(orphan, key=lambda r: r["path"]):
            L.append(f"| `{r['path']}` | {r['bytes']:,} | {r['modified']} |")
        L.append("")

    # results by stage
    L += ["## Results", ""]
    for tree in sorted({r.get("tree", "") for r in results}):
        sub = [r for r in results if r.get("tree") == tree]
        L += [f"### `{tree}/`", ""]
        for stage in sorted({r.get("stage", "") for r in sub}):
            ss = sorted([r for r in sub if r.get("stage") == stage],
                        key=lambda r: r["path"])
            if not ss:
                continue
            L += [f"**{stage or '(root)'}**", "",
                  "| file | state | n | verdict |", "|---|---|---|---|"]
            for r in ss:
                if not r.get("readable", True):
                    L.append(f"| `{Path(r['path']).name}` | UNREADABLE | | {r.get('error','')} |")
                    continue
                n = "; ".join(f"{k}={','.join(str(x) for x in v[:3])}"
                              for k, v in list(r["counts"].items())[:2])
                v = r["verdicts"][0]["text"] if r["verdicts"] else ""
                v = (v[:150] + "…") if len(v) > 150 else v
                L.append(f"| `{Path(r['path']).name}` | {r['state']} | {n} | {v} |")
            L.append("")

    # scripts by stage
    L += ["## Scripts", ""]
    for tree in sorted({s["tree"] for s in scripts}):
        sub = [s for s in scripts if s["tree"] == tree]
        L += [f"### `{tree}/` — {len(sub)} scripts", ""]
        for stage in sorted({s["stage"] for s in sub}):
            ss = sorted([s for s in sub if s["stage"] == stage],
                        key=lambda s: s["path"])
            L += [f"**{stage or '(root)'}** — {len(ss)}", "",
                  "| script | loc | purpose |", "|---|---|---|"]
            for s in ss:
                t = s["title"][:190] + ("…" if len(s["title"]) > 190 else "")
                L.append(f"| `{Path(s['path']).name}` | {s['loc']} | {t} |")
            L.append("")

    OUT_MD.write_text("\n".join(L), encoding="utf8")


def main():
    import sys
    scripts = scan_scripts()
    results, payloads = scan_results(scripts)

    if "--record" in sys.argv:
        n = record_provenance(results)
        print(f"recorded provenance baseline for {n} results -> "
              f"{PROV.relative_to(ROOT)}")
        results, payloads = scan_results(scripts)   # re-evaluate against it

    OUT_JSON.write_text(json.dumps({
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "tools/catalog_repo.py",
        "touch_tolerance_min": TOUCH_TOLERANCE_MIN,
        "counts": {"scripts": len(scripts), "results": len(results)},
        "scripts": scripts,
        "results": results,
    }, indent=1), encoding="utf8")

    OUT_RESULTS.write_text(json.dumps({
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "tools/catalog_repo.py",
        "note": ("Every ssw-design-analysis result in one file. The per-script "
                 "JSONs remain on disk as write-side provenance; this is the "
                 "read path. Regenerate with tools/catalog_repo.py."),
        "results": payloads,
    }, indent=1), encoding="utf8")

    write_md(scripts, results)

    stale = sum(r.get("state") == "STALE" for r in results)
    orphan = sum(r.get("state") == "ORPHAN" for r in results)
    unver = sum(r.get("state") == "unverified" for r in results)
    print(f"scripts  {len(scripts)}")
    print(f"results  {len(results)}   stale {stale}   orphan {orphan}"
          + (f"   unverified {unver}" if unver else ""))
    print(f"-> {OUT_JSON.relative_to(ROOT)}")
    print(f"-> {OUT_MD.relative_to(ROOT)}")
    print(f"-> {OUT_RESULTS.relative_to(ROOT)}  "
          f"({OUT_RESULTS.stat().st_size/1024:.0f} KB, {len(payloads)} results)")


if __name__ == "__main__":
    main()
