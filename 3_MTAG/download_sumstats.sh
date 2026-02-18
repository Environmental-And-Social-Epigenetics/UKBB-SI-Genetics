#!/usr/bin/env bash
set -euo pipefail

# Download MTAG-formatted sumstats from Luria into local 2_GWAS/.
# Override defaults with env vars:
#   REMOTE_USER, REMOTE_HOST, REMOTE_REPO_DIR

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPODIR="$(cd "${SCRIPTDIR}/.." && pwd)"
TARGET_2_GWAS_DIR="${REPODIR}/2_GWAS"

REMOTE_USER="${REMOTE_USER:-mabdel03}"
REMOTE_HOST="${REMOTE_HOST:-luria.mit.edu}"
REMOTE_REPO_DIR="${REMOTE_REPO_DIR:-/home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/UKBB-SI-Genetics}"

DRY_RUN_FLAG=""
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN_FLAG="--dry-run"
fi

if ! command -v rsync >/dev/null 2>&1; then
    echo "ERROR: rsync is required but not found in PATH." >&2
    exit 1
fi

if ! command -v ssh >/dev/null 2>&1; then
    echo "ERROR: ssh is required but not found in PATH." >&2
    exit 1
fi

mkdir -p "${TARGET_2_GWAS_DIR}/mtag_results"
mkdir -p "${TARGET_2_GWAS_DIR}/mtag_results_continuous"

REMOTE_BASE="${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_REPO_DIR}/2_GWAS"

echo "Sync source: ${REMOTE_BASE}"
echo "Sync target: ${TARGET_2_GWAS_DIR}"
if [[ -n "${DRY_RUN_FLAG}" ]]; then
    echo "Mode: dry-run (no files will be written)"
fi
echo

rsync -avz ${DRY_RUN_FLAG} \
    "${REMOTE_BASE}/mtag_results/" \
    "${TARGET_2_GWAS_DIR}/mtag_results/"

rsync -avz ${DRY_RUN_FLAG} \
    "${REMOTE_BASE}/mtag_results_continuous/" \
    "${TARGET_2_GWAS_DIR}/mtag_results_continuous/"

echo
echo "Done. MTAG sumstats are available at:"
echo "  ${TARGET_2_GWAS_DIR}/mtag_results/"
echo "  ${TARGET_2_GWAS_DIR}/mtag_results_continuous/"
