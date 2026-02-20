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

# Activate Snakemake (adjust if you use conda instead of micromamba)
eval "$(micromamba shell hook -s bash)"
micromamba activate snakemake

snakemake --profile profiles/slurm --directory "$PROJECT_DIR"
