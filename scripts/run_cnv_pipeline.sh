#!/bin/bash
#
#SBATCH --job-name=cnvkit_purecn_pipeline
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem=4000M
#SBATCH --output=slurm_logs/%x-%j.log

set -euo pipefail

################################################################################
# Usage:
#   sbatch run_cnv_pipeline.sh [MAX_JOBS]
#
#   - MAX_JOBS: Number of Snakemake jobs to run in parallel (default: 50)
#
# Description:
# This script launches the CNVkit-PureCN Snakemake pipeline. It sets up a
# temporary directory for intermediate files and uses a cluster profile
# for submitting individual jobs.
################################################################################

MAX_JOBS=${1:-50}
SNAKEMAKE_FILE="scripts/snakemake/cnv_pipeline.smk"
CONFIG_FILE="scripts/snakemake/config.yaml"

export TMPDIR="$HOME/scratch/tmp"
if [ ! -d "$TMPDIR" ]; then mkdir -p "$TMPDIR"; fi
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

mkdir -p slurm_logs
export SBATCH_DEFAULTS="--output=slurm_logs/%x-%j.log"

echo "INFO: Running CNV pipeline with:"
echo "      Snakefile:    $SNAKEMAKE_FILE"
echo "      Config file:  $CONFIG_FILE"
echo "      Max jobs:     $MAX_JOBS"
echo "      TMPDIR:       $TMPDIR"
date

srun snakemake \
    -s "$SNAKEMAKE_FILE" \
    --use-conda \
    --conda-prefix "conda_envs" \
    --profile=cubi-v1 \
    -j "$MAX_JOBS" \
    --configfile "$CONFIG_FILE" \
    --rerun-incomplete

date
echo "INFO: Pipeline finished."
