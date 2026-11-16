#!/usr/bin/env bash
#SBATCH --job-name=zam_dp_bins
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --array=1-3

set -euo pipefail

module purge
module load gnu/14.2.0
module load python/3.12.4
module load bcftools/1.17

ROOT="${ROOT:-/N/scratch/xta1/ZambiaElephants}"
INPUT_RUN_DIR="${INPUT_RUN_DIR:?Set INPUT_RUN_DIR when submitting}"
QC_RUN_DIR="${QC_RUN_DIR:?Set QC_RUN_DIR when submitting}"
SCRIPT="${SCRIPT:-$ROOT/slurm/scripts/genotype_dp_bins_zambia.py}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?This script must run as an array task}"

branches=(adapterremoval current_fastp fastp_merge)
if (( TASK_ID < 1 || TASK_ID > ${#branches[@]} )); then
    echo "ERROR: array task ID must be 1-${#branches[@]}: $TASK_ID" >&2
    exit 2
fi
BRANCH="${branches[$((TASK_ID - 1))]}"
BCF="$INPUT_RUN_DIR/$BRANCH/zambia_${BRANCH}.raw_biallelic_snps.bcf"
OUTDIR="$QC_RUN_DIR/$BRANCH"

if [[ ! -s "$SCRIPT" || ! -s "$BCF" || ! -s "$BCF.csi" ]]; then
    echo "ERROR: missing worker script, BCF, or CSI index." >&2
    echo "Script: $SCRIPT" >&2
    echo "BCF: $BCF" >&2
    exit 2
fi
mkdir -p "$OUTDIR"

echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Branch: $BRANCH"
echo "Input BCF: $BCF"
echo "Output directory: $OUTDIR"
echo "Python: $(python3 --version 2>&1)"
echo "BCFtools: $(bcftools --version | head -n 1)"

python3 "$SCRIPT" --bcf "$BCF" --branch "$BRANCH" --outdir "$OUTDIR"
