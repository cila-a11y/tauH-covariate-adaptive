#!/bin/bash
set -euo pipefail

: "${TAUH_PROJECT_DIR:=$PWD}"
: "${TAUH_ACCOUNT:=f202500010hpcvlabuminhox}"
: "${TAUH_PARTITION:=normal-x86}"
: "${TAUH_CORES:=64}"
: "${TAUH_TIME:=48:00:00}"
: "${TAUH_MAX_CONCURRENT:=4}"
: "${TAUH_N_FOCUSED_BOOT_TASKS:=8}"
: "${TAUH_SOURCE_RUN_ID:=paper_blockC_20260531_194657}"
: "${TAUH_SOURCE_RESULTS:=$TAUH_PROJECT_DIR/results_covariate_adaptive_tauH/$TAUH_SOURCE_RUN_ID}"
: "${TAUH_OUTPUT_ROOT:=$TAUH_PROJECT_DIR/results_focused_bootstrap_tauH}"
: "${TAUH_RUN_ID:=focused_bootstrap_$(date +%Y%m%d_%H%M%S)}"
: "${TAUH_FOCUSED_N_BOOT_DATASETS:=300}"
: "${TAUH_FOCUSED_N_BOOT_RESAMPLES:=500}"
: "${TAUH_FOCUSED_SCENARIO_IDS:=20,32,86,110,163,164,175,176}"

export TAUH_PROJECT_DIR TAUH_ACCOUNT TAUH_PARTITION TAUH_CORES TAUH_TIME TAUH_MAX_CONCURRENT
export TAUH_N_FOCUSED_BOOT_TASKS TAUH_SOURCE_RUN_ID TAUH_SOURCE_RESULTS TAUH_OUTPUT_ROOT TAUH_RUN_ID
export TAUH_FOCUSED_N_BOOT_DATASETS TAUH_FOCUSED_N_BOOT_RESAMPLES TAUH_FOCUSED_SCENARIO_IDS
export TAUH_FOCUSED_SCRIPT="$TAUH_PROJECT_DIR/focused_bootstrap_tauH_deucalion.R"
export TAUH_MAIN_SCRIPT="$TAUH_PROJECT_DIR/covariate_adaptive_tauH_simulation_deucalion.R"

SBATCH_SCRIPT="$TAUH_PROJECT_DIR/run_focused_bootstrap_deucalion.sbatch"
mkdir -p "$TAUH_PROJECT_DIR/slurm_logs" "$TAUH_OUTPUT_ROOT"

cat <<INFO
Submitting focused tauH bootstrap
  project_dir       = $TAUH_PROJECT_DIR
  focused_script    = $TAUH_FOCUSED_SCRIPT
  main_script       = $TAUH_MAIN_SCRIPT
  source_results    = $TAUH_SOURCE_RESULTS
  output_root       = $TAUH_OUTPUT_ROOT
  run_id            = $TAUH_RUN_ID
  account           = $TAUH_ACCOUNT
  partition         = $TAUH_PARTITION
  boot tasks        = $TAUH_N_FOCUSED_BOOT_TASKS, max concurrent = $TAUH_MAX_CONCURRENT
  datasets          = $TAUH_FOCUSED_N_BOOT_DATASETS
  resamples         = $TAUH_FOCUSED_N_BOOT_RESAMPLES
  scenario ids      = $TAUH_FOCUSED_SCENARIO_IDS
INFO

boot_job=$(sbatch --parsable \
  --account="$TAUH_ACCOUNT" \
  --partition="$TAUH_PARTITION" \
  --time="$TAUH_TIME" \
  --cpus-per-task="$TAUH_CORES" \
  --array="1-${TAUH_N_FOCUSED_BOOT_TASKS}%${TAUH_MAX_CONCURRENT}" \
  --export=ALL,TAUH_FOCUSED_MODE=bootstrap \
  "$SBATCH_SCRIPT")

echo "focused bootstrap array job: $boot_job"

collect_job=$(sbatch --parsable \
  --account="$TAUH_ACCOUNT" \
  --partition="$TAUH_PARTITION" \
  --time="08:00:00" \
  --cpus-per-task=4 \
  --dependency="afterok:$boot_job" \
  --export=ALL,TAUH_FOCUSED_MODE=collect,TAUH_CORES=4 \
  "$SBATCH_SCRIPT")

echo "focused collect job: $collect_job"
echo "Submitted. Monitor with: squeue --me"
echo "Final results will be under: $TAUH_OUTPUT_ROOT/$TAUH_RUN_ID"
