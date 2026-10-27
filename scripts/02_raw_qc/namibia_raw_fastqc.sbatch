#!/bin/bash
#SBATCH --job-name=namibia_raw_qc
#SBATCH -A c02431
#SBATCH --output=slurm/logs/namibia_raw_qc_%A_%a.out
#SBATCH --error=slurm/logs/namibia_raw_qc_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --array=1-8%8

set -euo pipefail

module purge
module load fastqc/0.12.1

ROOT="/N/scratch/xta1/ZambiaElephants"
LIST="$ROOT/qc/namibia_raw_fastq.list"
OUT="$ROOT/qc/namibia_raw_fastqc"

cd "$ROOT"
mkdir -p "$OUT"

FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${FASTQ:-}" || ! -s "$FASTQ" ]]; then
    echo "ERROR: invalid FASTQ: $FASTQ"
    exit 1
fi

echo "===== INPUT ====="
echo "$FASTQ"
ls -lh "$FASTQ"

fastqc \
    --threads "$SLURM_CPUS_PER_TASK" \
    --outdir "$OUT" \
    "$FASTQ"

echo "===== COMPLETE ====="
date
