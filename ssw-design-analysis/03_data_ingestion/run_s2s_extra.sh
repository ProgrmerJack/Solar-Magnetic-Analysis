#!/usr/bin/env bash
# Retrieval chain for the post-onset and SSW-hit tests (plan approved 2026-09-29).
# One ECDS request runs at a time per account, so the variables go in sequence.
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
SYS="ecmwf eccc cma hmcr kma cnrm jma cnr_isac ncep cptec"
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
