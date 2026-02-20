#!/bin/bash
#SBATCH --job-name=cuttag
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jmurray@mcw.edu
#SBATCH --partition=normal
#SBATCH --account=srrao
#SBATCH --time=3-00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Run Snakemake as a SLURM job. You get an email when the pipeline finishes (or fails).
# Usage: cd <pipeline_dir> && sbatch scripts/run_snakemake_slurm.sh <project_directory>
#
# Example:
#   cd /scratch/g/srrao/josiah_ociaml3/cut_tag_pipeline
#   sbatch scripts/run_snakemake_slurm.sh /scratch/g/srrao/josiah_ociaml3/snake_testing
#
# Edit --mail-user above with your email before submitting.

PROJECT_DIR="${1:?Usage: sbatch run_snakemake_slurm.sh <project_directory>}"
PIPELINE_DIR="${SLURM_SUBMIT_DIR}"

cd "$PIPELINE_DIR"

# Activate the Snakemake environment (uses micromamba; adjust if you use conda)
# ----------------------------------------
# Step 1: Check that micromamba is installed and on PATH. If not, the script
# would fail later with a confusing error; we fail early with instructions.
if ! command -v micromamba &>/dev/null; then
  echo "Error: micromamba not found. Install it first:"
  echo "  \"\${SHELL}\" <(curl -L micro.mamba.pm/install.sh)"
  exit 1
fi

# Step 2: Initialize micromamba for this shell. This sets up the 'activate'
# function and other helpers needed to switch environments.
eval "$(micromamba shell hook -s bash)"

# Step 3: Verify the 'snakemake' environment exists. If it doesn't, activating
# would fail; we fail early with the exact command to create it.
if ! micromamba run -n snakemake true 2>/dev/null; then
  echo "Error: micromamba environment 'snakemake' not found."
  echo "Create it with: micromamba create -c conda-forge -c bioconda -n snakemake snakemake snakemake-executor-plugin-slurm"
  exit 1
fi

# Step 4: Activate the environment so 'snakemake' is available below.
micromamba activate snakemake

snakemake --profile profiles/slurm --directory "$PROJECT_DIR"
