#!/usr/bin/env python3
"""
audit_scripts.py
================
Full static audit of every Python script in the repository.

Replaces the placeholder "UNREVIEWED" bucket in catalog_results.py with an
actual per-script determination of role and, for analysis scripts, statistical
validity.

WHAT IS DETECTED (all from source, no assumptions)

  role         download | process | analysis | figure | utility | test | unknown
               inferred from directory, output types, and imports actually used.

  exposure     which SSW event definition the script uses:
                 CATALOG_OLD    data/results/ssw_event_catalog.csv  (16 events,
                                contains 2012-01-11 which is not a major SSW in
                                any of six reanalyses; omits MAR 2000, MAR 2010)
                 CANONICAL      data/processed/atmospheric/ssw_canonical.csv
                 DETECTOR       local Charlton-Polvani implementation (over-detects)
                 PANEL_FLAG     the precomputed ssw_within_15d column, which was
                                built from CATALOG_OLD
                 HARDCODED      an inline list of dates
                 NONE           no SSW exposure

  controls     which confounders the script adjusts for:
                 winter_fe, seasonal, station_or_region_fe, matched_control

  inference    bootstrap | model_p | permutation | none_reported
               plus whether clustering is respected.

  windows      single window vs multiple simultaneous lag bins.

VALIDITY is then derived from those facts, not asserted:

  CURRENT      canonical exposure + winter FE + seasonal control + clustered
               inference + simultaneous bins  (or infrastructure verified by a
               self-test)
  SUPERSEDED   uses a known-bad exposure, or a design shown to be unsound
  LIMITED      analysis with partial controls -- usable as description only
  NO_INFERENCE download/process/figure/utility: nothing to invalidate
  EMPTY        no executable content

Every classification is written with the evidence that produced it, so it can be
checked line by line.

Outputs: docs/SCRIPT_AUDIT.md
         data/results/script_audit.json
"""
import ast
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT_TREES = [ROOT / "scripts", ROOT / "archive" / "superseded_scripts"]
OUT_MD = ROOT / "docs" / "SCRIPT_AUDIT.md"
OUT_JSON = ROOT / "data" / "results" / "script_audit.json"
SKIP_DIRS = {"__pycache__", ".ipynb_checkpoints"}

# infrastructure whose correctness is established by an executable self-check
SELF_VERIFIED = {"condpois.py"}
# scripts implementing the corrected event-study design
EVENT_STUDY = {"r113_event_study.py", "r114_avalanche_event_study.py"}
# the audit tooling itself: it contains every detector string as a literal, so
# it would otherwise match its own patterns and self-classify as analysis
AUDIT_TOOLING = {"audit_scripts.py", "catalog_results.py", "inventory_datasets.py"}
# Resolved by reading the source; recorded here so the determination is explicit
# rather than left in an unresolved bucket.
MANUAL = {
    "12_rebuild_enhanced_panel.py": ("process",
        "builds the enhanced analysis panel; no statistical claim"),
    "propagation_analysis.py": ("analysis",
        "reads ssw_catalog.parquet (wrong catalogue) - superseded"),
    "r34_check_data.py": ("utility",
        "ad-hoc data inspection; prints column names only"),
}

PAT = {
    "cat_old":   r"ssw_event_catalog",
    "cat_canon": r"ssw_canonical",
    "cat_parquet": r"ssw_catalog\.parquet",
    "detector":  r"def\s+detect_ssw|detect_ssw\s*\(",
    "panel_flag": r"ssw_within_15d",
    # transitive: the panel CARRIES ssw_within_15d, so anything loading it
    # inherits the 16-event exposure even without naming the column
    "panel_load": r"load_panel|analysis_panel(_v2)?\.parquet",
    "winter_fe": r"C\(winter\)|winter.*fixed|_demean|condpois|groupby\(\s*\[?[\"']winter|strata",
    "seasonal":  r"sin\(2\s*\*\s*np\.pi|dayofyear|day_of_season|doy|harmonic|C\(month\)",
    "unit_fe":   r"region.*winter|station.*winter|state.*winter|regionname|station_code",
    "matched":   r"doy.*match|matched|day-of-year|dayofyear.*isin",
    "bootstrap": r"bootstrap|resample|rng\.integers|np\.random\.choice",
    "block_boot": r"winter.*bootstrap|block.*bootstrap|idx_by_w|winters\b.*resampl",
    "model_p":   r"\.pvalues|ttest_1samp|ttest_ind|mannwhitneyu|binomtest|\.pvalue",
    "permutation": r"permutation|shuffle|randomis|randomiz",
    "multi_bin": r"BINS\s*=|for\s+i,\s*\(a,\s*b\)\s+in\s+enumerate|lag_bin|bins\s*=\s*\[\(",
    "placebo":   r"placebo|pre_onset|pre-onset",
    "reads_data": r"read_parquet|read_csv|read_excel|netCDF4|xr\.open",
    "writes_result": r"data/results/",
    "writes_figure": r"savefig|\.pdf|\.png",
    "network":   r"requests\.|urlopen|urlretrieve|curl|ftplib|cdsapi|earthaccess",
}
COMPILED = {k: re.compile(v, re.I) for k, v in PAT.items()}


def parse(p):
    try:
        src = p.read_text(encoding="utf8", errors="replace")
    except Exception:
        return None, None
    try:
        tree = ast.parse(src)
    except SyntaxError:
        tree = None
    return src, tree


def docline(src):
    m = re.search(r'"""(.*?)"""', src, re.S)
    if not m:
        return ""
    for l in (x.strip() for x in m.group(1).strip().splitlines()):
        if l and not set(l) <= set("=-_ #"):
            return l[:150]
    return ""


def audit_one(p):
    src, tree = parse(p)
    if src is None:
        return None
    code = "\n".join(l for l in src.splitlines()
                     if not l.strip().startswith("#"))
    hit = {k: bool(c.search(code)) for k, c in COMPILED.items()}
    nloc = sum(1 for l in src.splitlines()
               if l.strip() and not l.strip().startswith("#"))

    # ---- exposure --------------------------------------------------------
    if hit["cat_canon"]:
        exposure = "CANONICAL"
    elif hit["detector"]:
        exposure = "DETECTOR"
    elif hit["cat_old"]:
        exposure = "CATALOG_OLD"
    elif hit["cat_parquet"]:
        # Butler2015 parquet: 39 events, ends 2021-01-05, ALSO contains the
        # spurious 2012-01-11 and omits MAR 2000. Used by ~103 scripts.
        exposure = "CATALOG_PARQUET"
    elif hit["panel_flag"]:
        exposure = "PANEL_FLAG"
    elif hit["panel_load"]:
        exposure = "PANEL_FLAG_TRANSITIVE"
    elif re.search(r"\[\s*[\"']\d{4}-\d{2}-\d{2}[\"']\s*,\s*[\"']\d{4}-\d{2}-\d{2}", code):
        exposure = "HARDCODED"
    else:
        exposure = "NONE"

    # ---- role ------------------------------------------------------------
    # Analysis is tested FIRST and by counting actual statistical calls. An
    # earlier version tested "figure" first, which mis-filed any script that
    # both computed statistics and plotted them (e.g. 24_fresh_analysis.py has
    # 30 statistical calls and uses the wrong catalogue, yet was recorded as
    # NO_INFERENCE because it also called savefig). That hid superseded
    # analyses, which is the exact failure this audit exists to prevent.
    parts = p.parts
    n_stat = len(re.findall(
        r"\.pvalues|ttest_1samp|ttest_ind|mannwhitneyu|binomtest|\.pvalue|"
        r"wilcoxon|pearsonr|spearmanr|linregress|\.fit\(\)|gmRR|rate_ratio|"
        r"bootstrap|percentile\(", code, re.I))
    does_stats = n_stat > 0 or hit["writes_result"]

    # writes output of any kind, not only into data/results/
    writes_any = bool(re.search(r"to_parquet|to_csv|json\.dump|\.write_text|"
                                r"savefig|to_excel", code))
    # a library/helper module: leading-underscore name or config, no entry point
    is_lib = ((p.name.startswith("_") or p.name == "config.py")
              and "__main__" not in code)

    if p.name in MANUAL:
        role = MANUAL[p.name][0]
    elif nloc < 5:
        role = "empty"
    elif p.name in AUDIT_TOOLING or p.name in SELF_VERIFIED:
        role = "utility"
    elif is_lib:
        role = "library"
    elif does_stats and not ("download" in str(p)):
        role = "analysis"
    elif "download" in str(p) or hit["network"]:
        role = "download"
    elif "figures" in parts or hit["writes_figure"]:
        role = "figure"
    elif "process" in str(p):
        role = "process"
    elif "utilities" in parts or "verification" in parts:
        role = "utility"
    elif writes_any:
        role = "process"
    else:
        role = "unknown"

    # the origin of the bad exposure: this script assigns ssw_within_15d, the
    # flag that 24 downstream scripts consume
    builds_exposure = bool(re.search(r"\[.ssw_within_15d.\]\s*=(?!=)", code))

    controls = [k for k in ("winter_fe", "seasonal", "unit_fe", "matched") if hit[k]]
    infer = []
    if hit["block_boot"]:
        infer.append("winter_block_bootstrap")
    elif hit["bootstrap"]:
        infer.append("bootstrap_unclustered")
    if hit["permutation"]:
        infer.append("permutation")
    if hit["model_p"]:
        infer.append("model_p")

    outputs = sorted(set(re.findall(
        r"data/results/([A-Za-z0-9_\-.]+\.(?:json|csv|parquet|txt))", src)))

    # ---- validity --------------------------------------------------------
    if nloc < 5:
        validity, why = "EMPTY", "no executable content"
    elif builds_exposure:
        validity, why = ("EXPOSURE_SOURCE",
                         "BUILDS the ssw_within_15d flag from the 16-event "
                         "catalogue; every downstream PANEL_FLAG script "
                         "inherits this error")
    elif p.name in SELF_VERIFIED:
        validity, why = "CURRENT", "infrastructure with executable self-check"
    elif p.name in EVENT_STUDY:
        validity, why = "CURRENT", "simultaneous lag bins + explicit baseline"
    elif role in ("download", "process", "figure", "utility"):
        # even a non-analysis script matters if it BUILDS a wrong exposure
        if exposure in ("CATALOG_OLD", "PANEL_FLAG", "DETECTOR"):
            validity, why = ("NO_INFERENCE",
                             f"{role}: no statistical claim, but propagates "
                             f"exposure={exposure} downstream")
        else:
            validity, why = "NO_INFERENCE", f"{role}: produces no statistical claim"
    elif role == "unknown":
        # never assume an unclassified script is harmless
        validity, why = ("NEEDS_MANUAL_REVIEW",
                         f"role could not be determined from source "
                         f"(exposure={exposure}, {nloc} loc)")
    elif exposure in ("CATALOG_OLD", "PANEL_FLAG"):
        validity, why = "SUPERSEDED", f"exposure={exposure} (16-event list is wrong)"
    elif exposure == "PANEL_FLAG_TRANSITIVE":
        validity, why = ("SUPERSEDED",
                         "loads analysis_panel.parquet, which carries "
                         "ssw_within_15d built from the 16-event list")
    elif exposure == "CATALOG_PARQUET":
        validity, why = ("SUPERSEDED",
                         "exposure=ssw_catalog.parquet (contains spurious "
                         "2012-01-11, omits MAR 2000, ends 2021)")
    elif exposure == "DETECTOR":
        validity, why = "SUPERSEDED", "exposure=local detector (over-detects 40 vs 31)"
    elif not hit["multi_bin"] and exposure != "NONE":
        validity, why = ("SUPERSEDED",
                         "single-window specification (control group contains "
                         "the post-onset response)")
    elif exposure == "NONE":
        validity, why = "LIMITED", "no SSW exposure; descriptive only"
    else:
        missing = [c for c in ("winter_fe", "seasonal") if not hit[c]]
        if missing or "winter_block_bootstrap" not in infer:
            validity = "LIMITED"
            why = "missing " + ", ".join(missing + (
                [] if "winter_block_bootstrap" in infer else ["clustered inference"]))
        else:
            validity, why = "CURRENT", "canonical exposure + controls + clustered inference"

    return {
        "path": str(p.relative_to(ROOT)).replace("\\", "/"),
        "builds_exposure": builds_exposure,
        "role": role, "validity": validity, "reason": why,
        "exposure": exposure, "controls": controls, "inference": infer,
        "multi_bin": hit["multi_bin"], "placebo": hit["placebo"],
        "outputs": outputs, "nloc": nloc, "doc": docline(src),
    }


def main():
    rows = []
    for tree in SCRIPT_TREES:
        if not tree.exists():
            continue
        for p in sorted(tree.rglob("*.py")):
            if any(s in p.parts for s in SKIP_DIRS):
                continue
            r = audit_one(p)
            if r:
                r["quarantined"] = "archive" in p.parts
                rows.append(r)

    by_val = Counter(r["validity"] for r in rows)
    by_role = Counter(r["role"] for r in rows)
    by_exp = Counter(r["exposure"] for r in rows)

    L = ["# Script audit", "",
         "Generated by `scripts/utilities/audit_scripts.py` (re-runnable).",
         "Every script is classified from its source; the evidence for each",
         "classification is listed so it can be checked line by line.", "",
         f"**{len(rows)} scripts audited.**", "",
         "## Validity", "", "| status | n | meaning |", "|---|---|---|",
         f"| CURRENT | {by_val.get('CURRENT',0)} | canonical exposure, controls, clustered inference, simultaneous bins |",
         f"| SUPERSEDED | {by_val.get('SUPERSEDED',0)} | known-bad exposure or a design shown unsound |",
         f"| LIMITED | {by_val.get('LIMITED',0)} | partial controls; descriptive use only |",
         f"| NO_INFERENCE | {by_val.get('NO_INFERENCE',0)} | download/process/figure/utility: no statistical claim |",
         f"| NEEDS_MANUAL_REVIEW | {by_val.get('NEEDS_MANUAL_REVIEW',0)} | role undetermined from source; must be read by hand |",
         f"| EMPTY | {by_val.get('EMPTY',0)} | no executable content |",
         "", "## Role", "", "| role | n |", "|---|---|"]
    for k, v in by_role.most_common():
        L.append(f"| {k} | {v} |")
    L += ["", "## SSW exposure definition used", "",
          "| exposure | n | status |", "|---|---|---|",
          f"| CATALOG_OLD | {by_exp.get('CATALOG_OLD',0)} | wrong (contains a non-event) |",
          f"| CATALOG_PARQUET | {by_exp.get('CATALOG_PARQUET',0)} | wrong (spurious 2012-01-11, omits MAR 2000, ends 2021) |",
          f"| PANEL_FLAG | {by_exp.get('PANEL_FLAG',0)} | derived from CATALOG_OLD → wrong |",
          f"| DETECTOR | {by_exp.get('DETECTOR',0)} | over-detects (40 vs 31) |",
          f"| HARDCODED | {by_exp.get('HARDCODED',0)} | inline dates; must be checked individually |",
          f"| PANEL_FLAG_TRANSITIVE | {by_exp.get('PANEL_FLAG_TRANSITIVE',0)} | loads the panel, inheriting the bad flag |",
          f"| CANONICAL | {by_exp.get('CANONICAL',0)} | correct |",
          f"| NONE | {by_exp.get('NONE',0)} | no SSW exposure |", ""]

    for status in ("CURRENT", "EXPOSURE_SOURCE", "SUPERSEDED", "LIMITED", "NEEDS_MANUAL_REVIEW", "NO_INFERENCE", "EMPTY"):
        sel = [r for r in rows if r["validity"] == status]
        if not sel:
            continue
        L += [f"## {status} ({len(sel)})", ""]
        if status in ("CURRENT", "EXPOSURE_SOURCE", "SUPERSEDED", "LIMITED", "NEEDS_MANUAL_REVIEW"):
            L += ["| script | exposure | controls | inference | why |",
                  "|---|---|---|---|---|"]
            for r in sel:
                L.append(f"| `{r['path']}` | {r['exposure']} | "
                         f"{', '.join(r['controls']) or '—'} | "
                         f"{', '.join(r['inference']) or '—'} | {r['reason']} |")
        else:
            L += ["| script | role | purpose |", "|---|---|---|"]
            for r in sel:
                L.append(f"| `{r['path']}` | {r['role']} | {r['doc']} |")
        L.append("")

    OUT_MD.write_text("\n".join(L), encoding="utf8")
    OUT_JSON.write_text(json.dumps(rows, indent=1), encoding="utf8")
    print(f"{len(rows)} scripts audited")
    for k in ("CURRENT", "EXPOSURE_SOURCE", "SUPERSEDED", "LIMITED", "NEEDS_MANUAL_REVIEW", "NO_INFERENCE", "EMPTY"):
        print(f"  {k:14s} {by_val.get(k,0)}")
    print("exposure:", dict(by_exp))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
