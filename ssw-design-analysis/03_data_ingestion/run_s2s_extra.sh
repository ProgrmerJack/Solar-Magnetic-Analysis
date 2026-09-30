#!/usr/bin/env bash
# Retrieval chain for the post-onset and SSW-hit tests (plan approved 2026-09-29).
# One ECDS request runs at a time per account, so the variables go in sequence.
# Second chain (plans approved 2026-09-30, designs committed in 4c8b3f2 before
# retrieval): held-out model-year-2025 systems first, then the lead-lag variables.
# Usage: run_s2s_extra.sh [chain2]
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
SYS="ecmwf eccc cma hmcr kma cnrm jma cnr_isac ncep cptec"
run_var() {   # $1 variable, rest: origins
  local v=$1; shift
  local L=ssw-design-analysis/run_logs/03_data_ingestion/acquire_s2s_reforecasts_$v.log
  echo "=== $(date) packed $v $*" >> $L
  $PY -u -W ignore ssw-design-analysis/03_data_ingestion/acquire_s2s_reforecasts.py --var $v --packed "$@" >> $L 2>&1
  echo "packed exit=$?" >> $L
  for o in "$@"; do
    if [ "$v" = msl ]; then
      $PY -W ignore ssw-design-analysis/03_data_ingestion/acquire_s2s_reforecasts.py $o >> $L 2>&1
    else
      $PY -W ignore ssw-design-analysis/03_data_ingestion/acquire_s2s_reforecasts.py --var $v $o >> $L 2>&1
    fi
    echo "$o assemble exit=$?" >> $L
  done
}
if [ "${1:-}" = chain2 ]; then
  run_var msl ecmwf2025 cma2025
  run_var t2m ecmwf2025 cma2025
  run_var gh100 $SYS
  run_var u10_long $SYS
  echo "chain2 done $(date)"
  exit 0
fi
for v in u10 msl_short t2m_short; do
  L=ssw-design-analysis/run_logs/03_data_ingestion/acquire_s2s_reforecasts_$v.log
  $PY -u -W ignore ssw-design-analysis/03_data_ingestion/acquire_s2s_reforecasts.py --var $v --packed $SYS > $L 2>&1
  echo "packed exit=$?" >> $L
  for o in $SYS; do
    $PY -W ignore ssw-design-analysis/03_data_ingestion/acquire_s2s_reforecasts.py --var $v $o >> $L 2>&1
    echo "$o assemble exit=$?" >> $L
  done
done
echo "chain done $(date)"
