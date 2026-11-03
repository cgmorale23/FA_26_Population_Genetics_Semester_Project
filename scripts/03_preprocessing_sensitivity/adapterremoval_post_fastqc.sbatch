#!/bin/bash
#SBATCH --job-name=ar_post_fastqc
#SBATCH -A c02431
#SBATCH --output=slurm/logs/ar_post_fastqc_%A_%a.out
#SBATCH --error=slurm/logs/ar_post_fastqc_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --array=1-162%20

set -euo pipefail

module purge
module load fastqc/0.12.1

ROOT="/N/scratch/xta1/ZambiaElephants"
LIST="$ROOT/qc/adapterremoval_post_fastqc.list"
OUT="$ROOT/qc/adapterremoval_post_fastqc"

mkdir -p "$OUT"

FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${FASTQ:-}" || ! -s "$FASTQ" ]]; then
    echo "ERROR: missing or empty FASTQ: $FASTQ"
    exit 1
fi

echo "===== INPUT ====="
echo "$FASTQ"
ls -lh "$FASTQ"

echo
echo "===== SOFTWARE ====="
fastqc --version

fastqc \
    --threads "$SLURM_CPUS_PER_TASK" \
    --outdir "$OUT" \
    "$FASTQ"

echo
echo "===== COMPLETE ====="
date
