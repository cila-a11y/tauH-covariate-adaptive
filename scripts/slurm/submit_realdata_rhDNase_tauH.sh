#!/bin/bash
set -euo pipefail

TAUH_PROJECT_DIR="${TAUH_PROJECT_DIR:-/projects/F202500010HPCVLABUMINHO/cecilia/tauH_simulation}"
TAUH_ACCOUNT="${TAUH_ACCOUNT:-f202500010hpcvlabuminhox}"
TAUH_PARTITION="${TAUH_PARTITION:-normal-x86}"

TAUH_REAL_OUTPUT_ROOT="${TAUH_REAL_OUTPUT_ROOT:-$TAUH_PROJECT_DIR/results_realdata_rhDNase_tauH}"
TAUH_REAL_RUN_ID="${TAUH_REAL_RUN_ID:-rhDNase_real_$(date +%Y%m%d_%H%M%S)}"
TAUH_REAL_CORES="${TAUH_REAL_CORES:-32}"
TAUH_REAL_TIME="${TAUH_REAL_TIME:-02:00:00}"
TAUH_REAL_N_BOOT="${TAUH_REAL_N_BOOT:-200}"

mkdir -p "$TAUH_PROJECT_DIR/slurm_logs"
mkdir -p "$TAUH_REAL_OUTPUT_ROOT"

echo "Submitting rhDNase real-data tauH application"
echo "  project_dir  = $TAUH_PROJECT_DIR"
echo "  script       = $TAUH_PROJECT_DIR/realdata_rhDNase_tauH_application.R"
echo "  sbatch       = $TAUH_PROJECT_DIR/run_realdata_rhDNase_tauH.sbatch"
echo "  output_root  = $TAUH_REAL_OUTPUT_ROOT"
echo "  run_id       = $TAUH_REAL_RUN_ID"
echo "  account      = $TAUH_ACCOUNT"
echo "  partition    = $TAUH_PARTITION"
echo "  cores        = $TAUH_REAL_CORES"
echo "  time         = $TAUH_REAL_TIME"
echo "  n_boot       = $TAUH_REAL_N_BOOT"

jobid=$(sbatch --parsable \
  --job-name=tauH_rhDNase_real \
  --account="$TAUH_ACCOUNT" \
  --partition="$TAUH_PARTITION" \
  --nodes=1 \
  --cpus-per-task="$TAUH_REAL_CORES" \
  --time="$TAUH_REAL_TIME" \
  --output="$TAUH_PROJECT_DIR/slurm_logs/tauH_rhDNase_real_%j.out" \
  --error="$TAUH_PROJECT_DIR/slurm_logs/tauH_rhDNase_real_%j.err" \
  --export=ALL,TAUH_PROJECT_DIR="$TAUH_PROJECT_DIR",TAUH_REAL_OUTPUT_ROOT="$TAUH_REAL_OUTPUT_ROOT",TAUH_REAL_RUN_ID="$TAUH_REAL_RUN_ID",TAUH_REAL_CORES="$TAUH_REAL_CORES",TAUH_REAL_N_BOOT="$TAUH_REAL_N_BOOT" \
  "$TAUH_PROJECT_DIR/run_realdata_rhDNase_tauH.sbatch")

echo "rhDNase real-data job: $jobid"
echo "Final results will be under: $TAUH_REAL_OUTPUT_ROOT/$TAUH_REAL_RUN_ID"
