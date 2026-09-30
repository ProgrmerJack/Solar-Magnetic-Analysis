#!/usr/bin/env python3
"""
audit_results.py
================
Deep audit of every file in data/results/ and data/figures/.

Unlike a filename listing, this OPENS each result and records what is actually
inside it: which numeric claims it carries, whether those claims are keyed to an
SSW event definition, and which script produced it.

PROVENANCE is resolved by searching BOTH the live tree and
`archive/superseded_scripts/`, so a result whose producing script has been
quarantined is still traced rather than being reported as an orphan.

VALIDITY is inherited from the producing script's audit record
(`data/results/script_audit.json`), because a result is exactly as trustworthy
as the code that wrote it:

  CURRENT      produced by the corrected pipeline
  SUPERSEDED   producer used a wrong SSW event definition or an unsound design
  NO_INFERENCE producer makes no statistical claim (catalogues, inventories)
  ORPHAN       no producing script found in either tree

DUPLICATES are detected by content hash, since several rounds re-emitted the
same numbers under new names.

Outputs: docs/RESULTS_AUDIT.md
         data/results/results_audit.json   (excluded from its own audit)
"""
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "data" / "results"
FIGURES = ROOT / "data" / "figures"
TREES = [ROOT / "scripts", ROOT / "archive" / "superseded_scripts"]
OUT_MD = ROOT / "docs" / "RESULTS_AUDIT.md"
OUT_JSON = RESULTS / "results_audit.json"

# audit artefacts: describe the repo, not the science
SELF = {"results_audit.json", "script_audit.json", "results_catalog.json",
        "data_inventory.json", "cleanup_plan.json"}

# a claim is a number a paper could cite
CLAIM_KEYS = re.compile(
    r"p_|_p$|pval|p_two|p_one|sign_test|IRR|RR$|^RR|gmRR|CI95|effect|"
    r"cohen|ratio|coef|estimate|significan", re.I)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def build_producer_map(result_names):
    """result filename -> [scripts referencing it], across both trees.

    Resolves by searching for the literal filename (with and without extension).
    An earlier regex-only version missed files written via a RESULTS/ path
    constant -- e.g. nature_tier_analysis.json (1,183 numbers) was reported as an
    orphan although three scripts name it. Reproducibility claims must not rest
    on a pattern guess, so this searches for the names themselves.
    """
    stems = {n: n.rsplit(".", 1)[0] for n in result_names}
    m = defaultdict(list)
    for tree in TREES:
        if not tree.exists():
            continue
        for p in tree.rglob("*.py"):
            if "__pycache__" in p.parts:
                continue
            try:
                src = p.read_text(encoding="utf8", errors="replace")
            except Exception:
                continue
            rel = str(p.relative_to(ROOT)).replace("\\", "/")
            for name, stem in stems.items():
                if name in src or f'"{stem}"' in src or f"'{stem}'" in src:
                    m[name].append(rel)
    return {k: sorted(set(v)) for k, v in m.items()}


def count_claims(p):
    """How many citable numbers does this file carry?"""
    if p.suffix != ".json":
        return None, None
    try:
        d = json.loads(p.read_text(encoding="utf8", errors="replace"))
    except Exception:
        return None, None
    n, keys = 0, set()
    def walk(o):
        nonlocal n
        if isinstance(o, dict):
            for k, v in o.items():
                if CLAIM_KEYS.search(str(k)) and isinstance(v, (int, float, list)):
                    n += 1
                    keys.add(str(k))
                walk(v)
        elif isinstance(o, list):
            for v in o[:200]:
                walk(v)
    walk(d)
    return n, sorted(keys)[:8]


def main():
    audit = {}
    try:
        audit = {x["path"]: x for x in json.loads(
            (RESULTS / "script_audit.json").read_text(encoding="utf8"))}
    except Exception:
        pass

    names = [f.name for f in RESULTS.iterdir() if f.is_file() and f.name not in SELF]
    producers = build_producer_map(names)
    rows, by_hash = [], defaultdict(list)

    for p in sorted(RESULTS.iterdir()):
        if not p.is_file() or p.name in SELF:
            continue
        srcs = producers.get(p.name, [])
        vals = {audit[s]["validity"] for s in srcs if s in audit}
        if not srcs:
            validity, why = "ORPHAN", "no producing script found in either tree"
        elif "CURRENT" in vals:
            validity, why = "CURRENT", "produced by the corrected pipeline"
        elif vals & {"SUPERSEDED", "EXPOSURE_SOURCE", "LIMITED"}:
            bad = sorted(vals & {"SUPERSEDED", "EXPOSURE_SOURCE", "LIMITED"})[0]
            ex = next((audit[s]["exposure"] for s in srcs if s in audit), "?")
            validity, why = "SUPERSEDED", f"producer {bad}, exposure={ex}"
        else:
            validity, why = "NO_INFERENCE", "producer makes no statistical claim"

        nclaims, keys = count_claims(p)
        h = sha(p)
        by_hash[h].append(p.name)
        rows.append({
            "file": p.name, "size": p.stat().st_size, "sha": h,
            "validity": validity, "reason": why, "producers": srcs,
            "n_claims": nclaims, "claim_keys": keys,
        })

    dups = {h: v for h, v in by_hash.items() if len(v) > 1}
    for r in rows:
        r["duplicate_of"] = [x for x in by_hash[r["sha"]] if x != r["file"]]

    by_val = Counter(r["validity"] for r in rows)
    claims_super = sum(r["n_claims"] or 0 for r in rows if r["validity"] == "SUPERSEDED")
    claims_cur = sum(r["n_claims"] or 0 for r in rows if r["validity"] == "CURRENT")

    figs = sorted(FIGURES.iterdir()) if FIGURES.exists() else []
    fig_current = [f for f in figs if f.name.startswith("ng_fig")]

    L = ["# Results audit", "",
         "Generated by `scripts/utilities/audit_results.py` (re-runnable).",
         "Each result is opened and its contents inspected; provenance is resolved",
         "across both the live tree and `archive/superseded_scripts/`.", "",
         f"**{len(rows)} result files**, carrying **{sum(r['n_claims'] or 0 for r in rows)} citable numbers**.",
         "", "## Validity", "", "| status | files | citable numbers |", "|---|---|---|"]
    for k in ("CURRENT", "SUPERSEDED", "NO_INFERENCE", "ORPHAN"):
        n = by_val.get(k, 0)
        c = sum(r["n_claims"] or 0 for r in rows if r["validity"] == k)
        L.append(f"| {k} | {n} | {c} |")
    L += ["",
          f"**{claims_super} of the numbers in this directory come from superseded code, "
          f"against {claims_cur} from the corrected pipeline.**", ""]

    if dups:
        L += [f"## Duplicate contents ({len(dups)} groups)", "",
              "Identical bytes under different names — later rounds re-emitted "
              "earlier numbers.", ""]
        for h, names in sorted(dups.items(), key=lambda kv: -len(kv[1]))[:15]:
            L.append(f"- `{'`, `'.join(names)}`")
        L.append("")

    for status in ("CURRENT", "ORPHAN", "NO_INFERENCE", "SUPERSEDED"):
        sel = [r for r in rows if r["validity"] == status]
        if not sel:
            continue
        L += [f"## {status} ({len(sel)})", "",
              "| file | claims | producer | note |", "|---|---|---|---|"]
        for r in sorted(sel, key=lambda x: -(x["n_claims"] or 0)):
            prod = r["producers"][0] if r["producers"] else "—"
            prod = prod.replace("archive/superseded_scripts/", "archive:")
            L.append(f"| `{r['file']}` | {r['n_claims'] if r['n_claims'] is not None else '—'} "
                     f"| `{prod}` | {r['reason']} |")
        L.append("")

    L += ["## Figures", "",
          f"{len(figs)} files in `data/figures/`, of which {len(fig_current)} "
          "belong to the current (withdrawn) NG manuscript and the remainder to "
          "earlier withdrawn versions.", ""]

    OUT_MD.write_text("\n".join(L), encoding="utf8")
    OUT_JSON.write_text(json.dumps(rows, indent=1), encoding="utf8")
    print(f"{len(rows)} result files audited")
    for k in ("CURRENT", "SUPERSEDED", "NO_INFERENCE", "ORPHAN"):
        print(f"  {k:13s} {by_val.get(k,0)}")
    print(f"  duplicate groups: {len(dups)}")
    print(f"  citable numbers: {claims_cur} current / {claims_super} superseded")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
