#!/bin/sh
# One member from each cheap CMIP6 model with daily ua+psl.
# No grep in the pipeline: it buffers and hid all progress on the first attempt.
for M in CESM2-FV2 MPI-ESM1-2-LR IITM-ESM MPI-ESM-1-2-HAM INM-CM5-0; do
  echo "=== $M ==="
  python -u acquire_cmip6_ensemble.py 1 "$M"
done
