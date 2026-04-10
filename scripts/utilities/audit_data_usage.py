"""Audit: trace which processed data files are used by analysis scripts."""
import re
from pathlib import Path

scripts_dir = Path("scripts/analysis")
processed_dir = Path("data/processed")

# Get all processed data files
all_data = set()
for f in processed_dir.rglob("*.parquet"):
    all_data.add(str(f.relative_to(processed_dir)).replace("\\", "/"))
for f in processed_dir.rglob("*.nc"):
    all_data.add(str(f.relative_to(processed_dir)).replace("\\", "/"))

print(f"Total processed data files: {len(all_data)}")

# Also scan 01_build_daily_panel.py for what goes INTO the panel
panel_script = scripts_dir / "01_build_daily_panel.py"
panel_content = panel_script.read_text(encoding="utf-8") if panel_script.exists() else ""

# For each script, find PROCESSED / 'subdir' / 'filename' patterns
used_files = {}
all_referenced = set()

for script in sorted(scripts_dir.glob("*.py")):
    content = script.read_text(encoding="utf-8")
    refs = set()
    
    # Match PROCESSED / "subdir" / "file.parquet" patterns
    for m in re.finditer(
        r'PROCESSED\s*/\s*"(\w+)"\s*/\s*"([^"]+)"', content
    ):
        refs.add(f"{m.group(1)}/{m.group(2)}")
    
    # Also match single-quote variant
    for m in re.finditer(
        r"PROCESSED\s*/\s*'(\w+)'\s*/\s*'([^']+)'", content
    ):
        refs.add(f"{m.group(1)}/{m.group(2)}")
    
    # Match analysis_panel
    if "analysis_panel" in content:
        refs.add("analysis_panel.parquet")
    
    # Match RESULTS / "event_catalog"
    if "event_catalog" in content:
        refs.add("(results) event_catalog.parquet")
    
    if refs:
        used_files[script.name] = sorted(refs)
        all_referenced.update(refs)

# Print script-by-script
print("\n" + "=" * 70)
print("SCRIPT-BY-SCRIPT DATA USAGE:")
print("=" * 70)
for sn in sorted(used_files):
    print(f"\n  {sn}:")
    for r in used_files[sn]:
        print(f"    -> {r}")

# Find unreferenced files
print("\n" + "=" * 70)
print("DATA FILES NEVER REFERENCED IN ANY ANALYSIS SCRIPT:")
print("=" * 70)
not_used = []
for f in sorted(all_data):
    found = False
    for ref in all_referenced:
        if ref == f or f.endswith(ref) or ref.endswith(f):
            found = True
            break
        # Check partial match
        fname = f.split("/")[-1]
        rname = ref.split("/")[-1] if "/" in ref else ref
        if fname == rname:
            found = True
            break
    if not found:
        not_used.append(f)
        print(f"  NOT USED: {f}")

print(f"\nSummary: {len(all_data)} total files, {len(not_used)} not used in any analysis script")

# Check what's in the panel vs what could have been
print("\n" + "=" * 70)
print("PANEL COLUMN SOURCES (from 01_build_daily_panel.py):")
print("=" * 70)
panel_refs = set()
for m in re.finditer(r'PROCESSED\s*/\s*"(\w+)"\s*/\s*"([^"]+)"', panel_content):
    panel_refs.add(f"{m.group(1)}/{m.group(2)}")
for m in re.finditer(r"PROCESSED\s*/\s*'(\w+)'\s*/\s*'([^']+)'", panel_content):
    panel_refs.add(f"{m.group(1)}/{m.group(2)}")
for pr in sorted(panel_refs):
    print(f"  IN PANEL: {pr}")
