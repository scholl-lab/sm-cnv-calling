#!/bin/bash
#SBATCH --job-name=sm_cnv_calling
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_logs/%x-%j.log
#
# Snakemake 8+ SLURM launcher for sm-cnv-calling — auto-detects BIH, Charite, or local.
#
# Usage:
#   sbatch scripts/run_snakemake.sh [CONFIG_FILE] [EXTRA_ARGS...]
#   sbatch scripts/run_snakemake.sh workflow/Snakefile [CONFIG_FILE] [EXTRA_ARGS...]
#   bash scripts/run_snakemake.sh config/config.yaml --dry-run

set -euo pipefail

detect_cluster() {
    local fqdn
    fqdn=$(hostname -f 2>/dev/null || hostname)
    if [[ -d "/etc/xdg/snakemake/cubi-v1" ]] || [[ "$fqdn" =~ cubi|bihealth ]]; then
        echo "bih"
    elif [[ -f "/etc/profile.d/conda.sh" ]] || [[ "$fqdn" =~ charite|\.sc- ]]; then
        echo "charite"
    else
        echo "local"
    fi
}

CLUSTER=$(detect_cluster)

# Resolve Snakefile, config file, and stage-safe REPO_ROOT
if [[ "${1:-}" == *"Snakefile"* ]] || [[ "${1:-}" == *.smk ]]; then
    SNAKEFILE="$1"
    CONFIG_FILE="${2:-config/config.yaml}"
    shift 2 2>/dev/null || shift 1 2>/dev/null || true
    REPO_ROOT="$(cd "$(dirname "${SNAKEFILE}")/.." && pwd)"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -f "${SCRIPT_DIR}/../workflow/Snakefile" ]]; then
        REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
    else
        REPO_ROOT="$(pwd)"
    fi
    SNAKEFILE="${REPO_ROOT}/workflow/Snakefile"
    CONFIG_FILE="${1:-config/config.yaml}"
    shift 1 2>/dev/null || true
fi

# Conda activation (Charite requires explicit sourcing)
if ! command -v snakemake >/dev/null 2>&1; then
    if [[ "$CLUSTER" == "charite" ]] && [[ -f /etc/profile.d/conda.sh ]]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/conda.sh
    fi
    conda activate "${SNAKEMAKE_CONDA_ENV:-snakemake}" 2>/dev/null || true
fi

# TMPDIR setup (honour existing TMPDIR; default to ~/scratch/tmp on BIH)
if [[ "$CLUSTER" == "bih" ]]; then
    BASE_TMPDIR="${TMPDIR:-${HOME}/scratch/tmp}"
else
    BASE_TMPDIR="${TMPDIR:-/tmp}/snakemake"
fi
mkdir -p "${BASE_TMPDIR}"
TMPDIR=$(mktemp -d "${BASE_TMPDIR}/sm.XXXXXX")
export TMPDIR
trap 'rm -rf "${TMPDIR}"' EXIT

mkdir -p slurm_logs
export SBATCH_DEFAULTS="--output=slurm_logs/%x-%j.log"

echo "=== sm-cnv-calling Launch ==="
echo "  Cluster:    ${CLUSTER}"
echo "  Snakefile:  ${SNAKEFILE}"
echo "  Config:     ${CONFIG_FILE}"
echo "  Profile:    ${REPO_ROOT}/profiles/${CLUSTER}"
echo "  TMPDIR:     ${TMPDIR}"
echo "  Extra args: $*"
echo "  Start:      $(date)"
echo "============================="

snakemake \
    -s "${SNAKEFILE}" \
    --configfile "${CONFIG_FILE}" \
    --workflow-profile "${REPO_ROOT}/profiles/default" \
    --profile "${REPO_ROOT}/profiles/${CLUSTER}" \
    "$@"

echo "=== Finished: $(date) ==="
