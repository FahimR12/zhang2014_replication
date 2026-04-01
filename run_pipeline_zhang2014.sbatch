#!/bin/bash
# Run Zhang et al. (2014) TCGA-OV replication pipeline on SLURM.
#
# Submit:
#   sbatch /users/frazai/fastscratch/zhang2014_replication/run_pipeline_zhang2014.sbatch
#
# Optional exports:
#   --export=ALL,FETCH_DATA=1              # run 00_fetch_gdc_data.R first
#   --export=ALL,RUN_STEPS="1"             # run only step 1
#   --export=ALL,RUN_STEPS="3 5"           # run step range 3..5
#   --export=ALL,TCGA_OV_DATA_ROOT=/path/to/TCGA-OV
#   --export=ALL,AUTO_FETCH_IF_MISSING=0   # disable auto-fetch on missing data
#   --export=ALL,HEARTBEAT_SEC=60          # heartbeat interval in seconds

#SBATCH -J zhang2014-ov
#SBATCH -p nodes
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH -D /users/frazai/fastscratch/zhang2014_replication
#SBATCH -o /users/frazai/fastscratch/zhang2014_replication/logs/%x-%j.out
#SBATCH -e /users/frazai/fastscratch/zhang2014_replication/logs/%x-%j.err

set -euo pipefail
mkdir -p /users/frazai/fastscratch/zhang2014_replication/logs

# Activate conda env used for this project.
# Preferred: conda.sh + conda activate
# Fallback: prepend env bin directly when conda.sh is unavailable on compute nodes.
CONDA_ENV_NAME="${CONDA_ENV_NAME:-zhang2014-r}"
CONDA_ENV_PREFIX="${CONDA_ENV_PREFIX:-/users/frazai/fastscratch/conda/envs/${CONDA_ENV_NAME}}"
CONDA_SH_CANDIDATES=(
  "/users/frazai/fastscratch/conda/etc/profile.d/conda.sh"
  "/mnt/fastscratch/users/frazai/conda/etc/profile.d/conda.sh"
  "/users/frazai/.conda/etc/profile.d/conda.sh"
  "/mnt/fastscratch/users/frazai/.conda/etc/profile.d/conda.sh"
)

CONDA_ACTIVATED=0
for conda_sh in "${CONDA_SH_CANDIDATES[@]}"; do
  if [[ -f "${conda_sh}" ]]; then
    # shellcheck disable=SC1090
    source "${conda_sh}"
    conda activate "${CONDA_ENV_NAME}"
    CONDA_ACTIVATED=1
    break
  fi
done

if [[ "${CONDA_ACTIVATED}" -eq 0 ]]; then
  if [[ -x "${CONDA_ENV_PREFIX}/bin/Rscript" ]]; then
    export PATH="${CONDA_ENV_PREFIX}/bin:${PATH}"
    export CONDA_DEFAULT_ENV="${CONDA_ENV_NAME}"
  else
    echo "[zhang2014] ERROR: could not activate conda env."
    echo "[zhang2014] Tried conda.sh at: ${CONDA_SH_CANDIDATES[*]}"
    echo "[zhang2014] Also checked env prefix: ${CONDA_ENV_PREFIX}"
    exit 1
  fi
fi

export OMP_NUM_THREADS=1

echo "[zhang2014] job=${SLURM_JOB_ID}"
echo "[zhang2014] cwd=$(pwd)"
echo "[zhang2014] conda_env=${CONDA_DEFAULT_ENV:-unknown}"
echo "[zhang2014] R=$(which Rscript)"

# Default DATA_ROOT for 01_data_loading.R unless overridden by caller.
if [[ -z "${TCGA_OV_DATA_ROOT:-}" ]]; then
  export TCGA_OV_DATA_ROOT="$(pwd)/data/GDCdata/TCGA-OV"
fi
echo "[zhang2014] TCGA_OV_DATA_ROOT=${TCGA_OV_DATA_ROOT}"

if [[ "${FETCH_DATA:-0}" == "1" ]]; then
  echo "[zhang2014] ── Fetching GDC data ─────────────────────────────────────"
  Rscript 00_fetch_gdc_data.R
fi

# First run convenience: auto-fetch if data root is missing.
if [[ ! -d "${TCGA_OV_DATA_ROOT}" && "${AUTO_FETCH_IF_MISSING:-1}" == "1" ]]; then
  echo "[zhang2014] data root missing, auto-fetching: ${TCGA_OV_DATA_ROOT}"
  Rscript 00_fetch_gdc_data.R
fi

echo "[zhang2014] ── Running pipeline ─────────────────────────────────────────"
STATUS_FILE="$(pwd)/logs/pipeline_status-${SLURM_JOB_ID}.txt"
export PIPELINE_STATUS_FILE="${STATUS_FILE}"
echo "[zhang2014] status_file=${STATUS_FILE}"

HEARTBEAT_SEC="${HEARTBEAT_SEC:-60}"
HEARTBEAT_PID=""
if [[ "${HEARTBEAT_SEC}" =~ ^[0-9]+$ ]] && [[ "${HEARTBEAT_SEC}" -gt 0 ]]; then
  (
    while true; do
      sleep "${HEARTBEAT_SEC}"
      ts="$(date '+%Y-%m-%d %H:%M:%S')"
      current_status="$(tail -n 1 "${STATUS_FILE}" 2>/dev/null || echo "status pending")"
      echo "[zhang2014][heartbeat][${ts}] ${current_status}"
    done
  ) &
  HEARTBEAT_PID=$!
  trap 'if [[ -n "${HEARTBEAT_PID}" ]]; then kill "${HEARTBEAT_PID}" 2>/dev/null || true; fi' EXIT
fi

if [[ -n "${RUN_STEPS:-}" ]]; then
  echo "[zhang2014] steps=${RUN_STEPS}"
  # shellcheck disable=SC2086
  Rscript run_pipeline.R ${RUN_STEPS}
else
  echo "[zhang2014] steps=all"
  Rscript run_pipeline.R
fi

if [[ -n "${HEARTBEAT_PID}" ]]; then
  kill "${HEARTBEAT_PID}" 2>/dev/null || true
fi

echo "[zhang2014] done"
