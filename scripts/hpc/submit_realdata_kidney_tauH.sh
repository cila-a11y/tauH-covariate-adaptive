#!/bin/bash
set -euo pipefail

PROJECT_DIR=${TAUH_PROJECT_DIR:-$(pwd)}
ACCOUNT=${TAUH_ACCOUNT:?Set TAUH_ACCOUNT, e.g. f202500010hpcvlabuminhox}
PARTITION=${TAUH_PARTITION:-normal-x86}
SCRIPT=${TAUH_REAL_SCRIPT:-$PROJECT_DIR/realdata_kidney_tauH_application.R}
SBATCH_SCRIPT=${TAUH_REAL_SBATCH_SCRIPT:-$PROJECT_DIR/run_realdata_kidney_tauH.sbatch}
RUN_ID=${TAUH_REAL_RUN_ID:-kidney_realdata_$(date +%Y%m%d_%H%M%S)}
OUTPUT_ROOT=${TAUH_REAL_OUTPUT_ROOT:-$PROJECT_DIR/results_realdata_kidney_tauH}
CORES=${TAUH_REAL_CORES:-64}
TIME=${TAUH_REAL_TIME:-08:00:00}
N_BOOT=${TAUH_REAL_N_BOOT:-2000}

mkdir -p "$PROJECT_DIR/slurm_logs" "$OUTPUT_ROOT"

export TAUH_PROJECT_DIR="$PROJECT_DIR"
export TAUH_ACCOUNT="$ACCOUNT"
export TAUH_PARTITION="$PARTITION"
export TAUH_REAL_SCRIPT="$SCRIPT"
export TAUH_REAL_RUN_ID="$RUN_ID"
export TAUH_REAL_OUTPUT_ROOT="$OUTPUT_ROOT"
export TAUH_REAL_CORES="$CORES"
export TAUH_REAL_TIME="$TIME"
export TAUH_REAL_N_BOOT="$N_BOOT"

cat <<MSG
Submitting real-data kidney tauH application
  project_dir  = $PROJECT_DIR
  script       = $SCRIPT
  sbatch       = $SBATCH_SCRIPT
  output_root  = $OUTPUT_ROOT
  run_id       = $RUN_ID
  account      = $ACCOUNT
  partition    = $PARTITION
  cores        = $CORES
  time         = $TIME
  n_boot       = $N_BOOT
MSG

jid=$(sbatch --parsable \
  --account="$ACCOUNT" \
  --partition="$PARTITION" \
  --cpus-per-task="$CORES" \
  --time="$TIME" \
  --export=ALL \
  "$SBATCH_SCRIPT")

echo "kidney real-data job: $jid"
echo "Final results will be under: $OUTPUT_ROOT/$RUN_ID"
