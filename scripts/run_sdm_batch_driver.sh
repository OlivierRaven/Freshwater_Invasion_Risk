#!/bin/bash
# Restarts scripts/run_sdm_batch.R repeatedly until every species in
# bulk_target_species has a "success" row in sdm_batch_results.csv.
#
# Why this exists (2026-10-29): running all 353 species in one long-lived R
# session leaked memory somewhere in the terra/GDAL/biomod2 stack -- process
# RSS grew from ~111MB to ~960MB over 38 species, filled swap on this 8GB
# machine, and per-species runtime degraded from ~5 min to 169 min as it
# started thrashing. run_sdm_batch.R now processes at most SDM_BATCH_LIMIT
# species (default 15) per invocation and exits cleanly; this driver just
# keeps restarting it with a fresh process (and a short pause to let the OS
# reclaim memory) until the results CSV shows every species done.
#
# Usage: nohup scripts/run_sdm_batch_driver.sh > /tmp/sdm_batch_driver.log 2>&1 &

cd "/Users/olivierraven/Documents/Education:work/Freshwater_Invasion_Risk" || exit 1

TARGET_LIST="data/derived/bulk_353_species_target_list.csv"
RESULTS="data/derived/sdm_batch_results.csv"

total=$(($(wc -l < "$TARGET_LIST") - 1))

while true; do
  done_count=0
  if [ -f "$RESULTS" ]; then
    done_count=$(awk -F, 'NR>1 && $2=="success"' "$RESULTS" | wc -l | tr -d ' ')
  fi

  if [ "$done_count" -ge "$total" ]; then
    echo "[driver] All $total species done. Exiting."
    break
  fi

  echo "[driver] $done_count/$total done. Launching next batch ($(date))..."
  Rscript scripts/run_sdm_batch.R
  echo "[driver] Batch invocation exited. Pausing 15s before next restart..."
  sleep 15
done
