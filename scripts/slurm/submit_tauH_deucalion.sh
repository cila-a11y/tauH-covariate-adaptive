#!/bin/bash
set -euo pipefail

# ------------------------------------------------------------------------------
# Submission pipeline for Deucalion: truth -> MC array -> bootstrap array -> collect
# It is designed to use at most 4 full x86 nodes concurrently, each with 128 cores.
# ------------------------------------------------------------------------------

PROJECT_DIR="${TAUH_PROJECT_DIR:-$PWD}"
SCRIPT="${TAUH_SCRIPT:-$PROJECT_DIR/covariate_adaptive_tauH_simulation_deucalion.R}"
SBATCH_SCRIPT="${TAUH_SBATCH_SCRIPT:-$PROJECT_DIR/run_tauH_deucalion.sbatch}"
OUTPUT_ROOT="${TAUH_OUTPUT_ROOT:-$PROJECT_DIR/results_covariate_adaptive_tauH}"
RUN_ID="${TAUH_RUN_ID:-tauH_$(date +%Y%m%d_%H%M%S)}"
PROFILE="${TAUH_PROFILE:-paper}"
ACCOUNT="${TAUH_ACCOUNT:-YOUR_ACCOUNT}"
PARTITION="${TAUH_PARTITION:-normal-x86}"
CORES="${TAUH_CORES:-64}"
TIME_LIMIT="${TAUH_TIME:-48:00:00}"

# More tasks than nodes improves load balancing; %MAX_CONCURRENT keeps node usage capped.
N_SCENARIO_TASKS="${TAUH_N_SCENARIO_TASKS:-16}"
N_BOOT_TASKS="${TAUH_N_BOOT_TASKS:-4}"
MAX_CONCURRENT="${TAUH_MAX_CONCURRENT:-4}"
RUN_BOOTSTRAP="${TAUH_RUN_BOOTSTRAP:-1}"

mkdir -p "$PROJECT_DIR/slurm_logs" "$OUTPUT_ROOT"

common_export="ALL,TAUH_PROJECT_DIR=$PROJECT_DIR,TAUH_SCRIPT=$SCRIPT,TAUH_OUTPUT_ROOT=$OUTPUT_ROOT,TAUH_RUN_ID=$RUN_ID,TAUH_PROFILE=$PROFILE,TAUH_CORES=$CORES,TAUH_RUN_BOOTSTRAP=$RUN_BOOTSTRAP"
common_opts=(--account="$ACCOUNT" --partition="$PARTITION" --cpus-per-task="$CORES" --time="$TIME_LIMIT")

echo "Submitting tauH pipeline"
echo "  project_dir       = $PROJECT_DIR"
echo "  script            = $SCRIPT"
echo "  sbatch_script     = $SBATCH_SCRIPT"
echo "  output_root       = $OUTPUT_ROOT"
echo "  run_id            = $RUN_ID"
echo "  account           = $ACCOUNT"
echo "  partition         = $PARTITION"
echo "  scenario tasks    = $N_SCENARIO_TASKS, max concurrent = $MAX_CONCURRENT"
echo "  bootstrap tasks   = $N_BOOT_TASKS, max concurrent = $MAX_CONCURRENT"

jid_truth=$(sbatch --parsable "${common_opts[@]}" \
  --job-name=tauH_truth \
  --export="$common_export,TAUH_MODE=truth" \
  "$SBATCH_SCRIPT")
echo "truth job: $jid_truth"

jid_mc=$(sbatch --parsable "${common_opts[@]}" \
  --dependency=afterok:$jid_truth \
  --array=1-${N_SCENARIO_TASKS}%${MAX_CONCURRENT} \
  --job-name=tauH_mc \
  --export="$common_export,TAUH_MODE=mc,TAUH_N_SCENARIO_TASKS=$N_SCENARIO_TASKS" \
  "$SBATCH_SCRIPT")
echo "MC array job: $jid_mc"

if [[ "$RUN_BOOTSTRAP" == "1" ]]; then
  jid_boot=$(sbatch --parsable "${common_opts[@]}" \
    --dependency=afterok:$jid_mc \
    --array=1-${N_BOOT_TASKS}%${MAX_CONCURRENT} \
    --job-name=tauH_boot \
    --export="$common_export,TAUH_MODE=bootstrap,TAUH_N_BOOT_TASKS=$N_BOOT_TASKS" \
    "$SBATCH_SCRIPT")
  echo "bootstrap array job: $jid_boot"
  dep_collect="afterok:$jid_mc:$jid_boot"
else
  jid_boot=""
  dep_collect="afterok:$jid_mc"
fi

jid_collect=$(sbatch --parsable "${common_opts[@]}" \
  --dependency=$dep_collect \
  --job-name=tauH_collect \
  --export="$common_export,TAUH_MODE=collect" \
  "$SBATCH_SCRIPT")
echo "collect job: $jid_collect"

echo "Submitted. Monitor with: squeue --me"
echo "Final results will be under: $OUTPUT_ROOT/$RUN_ID"
