#!/bin/sh
# One step of the registered prospective test (prospective_next_ssw.py):
# extend ERA5T to six days before today (values frozen as first retrieved),
# detect onsets, score every onset whose outcome window is complete.
# Safe to run repeatedly (e.g. weekly from cron); appends to the run log.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
SSW="$HERE/.."
PY="$SSW/../.venv/bin/python"
LOG="$SSW/run_logs/08_literature_audit/run_prospective.log"
LAST=$("$PY" -c "import pandas as pd; print((pd.read_parquet('$SSW/03_data_ingestion/era5_realtime_cds.parquet').index.max() + pd.Timedelta(days=1)).date())" 2>/dev/null)
END=$(date -d "6 days ago" +%F)
{
  echo "=== $(date -u +%FT%TZ)  data $LAST .. $END"
  if [ "$LAST" \< "$END" ] || [ "$LAST" = "$END" ]; then
    "$PY" "$SSW/03_data_ingestion/acquire_era5_realtime_cds.py" --start "$LAST" --end "$END" 2>&1 | grep -v -i "runtimewarning\|frozen site\|%|" | tail -2
  fi
  "$PY" "$HERE/prospective_next_ssw.py" --detect 2>&1 | grep -v -i "runtimewarning\|frozen site"
  "$PY" "$HERE/prospective_next_ssw.py" --score 2>&1 | grep -v -i "runtimewarning\|frozen site" | tail -12
} >> "$LOG" 2>&1
tail -20 "$LOG"
