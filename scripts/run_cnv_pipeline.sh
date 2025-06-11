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
#   sbatch run_cnv_pipeline.sh [SNAKEFILE] [CONFIG_FILE] [MAX_JOBS]
#
#   - SNAKEFILE: Path to the Snakemake file (default: scripts/snakemake/cnv_pipeline.smk)
#   - CONFIG_FILE: Path to the config file (default: scripts/snakemake/config.yaml)
#   - MAX_JOBS: Number of Snakemake jobs to run in parallel (default: 50)
#
# Description:
# This script launches the CNVkit-PureCN Snakemake pipeline. It sets up a
# temporary directory for intermediate files and uses a cluster profile
# for submitting individual jobs.
################################################################################

SNAKEMAKE_FILE=${1:-"scripts/snakemake/cnv_pipeline.smk"}
CONFIG_FILE=${2:-"scripts/snakemake/config.yaml"}
MAX_JOBS=${3:-50}

# Validate that MAX_JOBS is a number
if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: MAX_JOBS must be a positive integer, got: $MAX_JOBS"
    exit 1
fi

# Validate that input files exist
if [[ ! -f "$SNAKEMAKE_FILE" ]]; then
    echo "ERROR: Snakefile not found: $SNAKEMAKE_FILE"
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

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
