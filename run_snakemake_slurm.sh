#!/bin/bash
# Run from the pipeline directory. Set PROJECT_DIR and run this script directly (do not use sbatch).
#   PROJECT_DIR=/path/to/project ./run_snakemake_slurm.sh
# The script will submit itself via sbatch with logs in PROJECT_DIR (SLURM does not expand variables in #SBATCH).
#
#SBATCH --job-name=cuttag
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jmurray@mcw.edu # Edit with your email before submitting.
#SBATCH --partition=normal
#SBATCH --account=srrao
#SBATCH --time=3-00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# When run directly (not already in a SLURM job), re-submit with logs in PROJECT_DIR
if [[ -z "${SLURM_JOB_ID:-}" && -n "${PROJECT_DIR:-}" ]]; then
  export PROJECT_DIR
  # Export PIPELINE_DIR so the job gets the real path (SLURM may run the script from a copy, so $0 would be wrong in the job)
  export PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
  exec sbatch --export=ALL \
    -o "${PROJECT_DIR}/cuttag_%j.out" \
    -e "${PROJECT_DIR}/cuttag_%j.err" \
    "$PIPELINE_DIR/$(basename "$0")"
fi

: "${PROJECT_DIR:?Set PROJECT_DIR and run the script (do not use sbatch): PROJECT_DIR=/path/to/project ./run_snakemake_slurm.sh}"
export PROJECT_DIR

# Use PIPELINE_DIR from environment if set (by re-submit); otherwise from script location (e.g. when run inside an interactive job)
if [[ -z "${PIPELINE_DIR:-}" ]]; then
  PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
fi
cd "$PROJECT_DIR"

# Activate the Snakemake environment (uses micromamba; adjust if you use conda)
# ----------------------------------------
# Step 1: Check that micromamba is installed and on PATH.
if ! command -v micromamba &>/dev/null; then
  echo "Error: micromamba not found. Install it first:"
  echo "  \"\${SHELL}\" <(curl -L micro.mamba.pm/install.sh)"
  exit 1
fi

# Step 2: Initialize micromamba for this shell. This sets up the 'activate'
# function and other helpers needed to switch environments.
eval "$(micromamba shell hook -s bash)"

# Step 3: Verify the 'snakemake' environment exists.
if ! micromamba run -n snakemake true 2>/dev/null; then
  echo "Error: micromamba environment 'snakemake' not found."
  echo "Create it with: micromamba create -c conda-forge -c bioconda -n snakemake snakemake snakemake-executor-plugin-slurm"
  exit 1
fi

micromamba activate snakemake


# Use only our profile (unset any cluster/OOD default that might override --profile)
unset SNAKEMAKE_PROFILE
# Give shared FS (e.g. Lustre) time to make job outputs visible before assuming failure (default 5s is often too short)
snakemake -s "$PIPELINE_DIR/Snakefile" --configfile config/config.yaml --profile "$PIPELINE_DIR/profiles/slurm" --latency-wait 60
