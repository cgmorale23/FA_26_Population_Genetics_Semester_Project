#!/bin/bash
#SBATCH --job-name=public_fastp
#SBATCH -A c02431
#SBATCH --output=slurm/logs/public_fastp_%A_%a.out
#SBATCH --error=slurm/logs/public_fastp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --array=1-42%20

set -euo pipefail

module purge
module load fastp/0.23.4

ROOT="/N/scratch/xta1/ZambiaElephants"

MANIFEST="$ROOT/pilot_trimming/public_pairs.clean.tsv"

OUT="$ROOT/processed/public_fastp"
REPORT="$ROOT/qc/public_fastp_reports"

mkdir -p "$OUT"
mkdir -p "$REPORT"

cd "$ROOT"

LINE=$((SLURM_ARRAY_TASK_ID + 1))

IFS=$'\t' read -r PAIR R1 R2 < <(
    sed -n "${LINE}p" "$MANIFEST"
)

if [[ -z "${PAIR:-}" || ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: invalid input"
    echo "PAIR=$PAIR"
    echo "R1=$R1"
    echo "R2=$R2"
    exit 1
fi

echo "===== SOFTWARE ====="
command -v fastp
fastp --version || true

echo
echo "===== INPUT ====="
echo "PAIR=$PAIR"
echo "R1=$R1"
echo "R2=$R2"

ls -lh "$R1" "$R2"

echo
echo "===== FASTP ====="

fastp \
    --in1 "$R1" \
    --in2 "$R2" \
    --out1 "$OUT/${PAIR}_R1.fastq.gz" \
    --out2 "$OUT/${PAIR}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --poly_g_min_len 10 \
    --disable_quality_filtering \
    --disable_length_filtering \
    --thread "$SLURM_CPUS_PER_TASK" \
    --json "$REPORT/${PAIR}.fastp.json" \
    --html "$REPORT/${PAIR}.fastp.html" \
    --report_title "${PAIR} fastp preprocessing"

echo
echo "===== COMPLETE ====="
echo "$PAIR"
date
