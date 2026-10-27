#!/bin/bash
#SBATCH --job-name=namibia_fastp
#SBATCH -A c02431
#SBATCH --output=slurm/logs/namibia_fastp_%A_%a.out
#SBATCH --error=slurm/logs/namibia_fastp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --array=1-4%4

set -euo pipefail

module purge
module load fastp/0.23.4

ROOT="/N/scratch/xta1/ZambiaElephants"
MANIFEST="$ROOT/qc/namibia_pairs.tsv"
OUT="$ROOT/processed/namibia_fastp"
REPORT="$ROOT/qc/namibia_fastp_reports"

cd "$ROOT"
mkdir -p "$OUT" "$REPORT"

LINE=$((SLURM_ARRAY_TASK_ID + 1))

IFS=$'\t' read -r SAMPLE R1 R2 < <(
    sed -n "${LINE}p" "$MANIFEST"
)

if [[ -z "${SAMPLE:-}" || ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: invalid input"
    echo "SAMPLE=$SAMPLE"
    echo "R1=$R1"
    echo "R2=$R2"
    exit 1
fi

echo "===== INPUT ====="
echo "Sample=$SAMPLE"
echo "R1=$R1"
echo "R2=$R2"
ls -lh "$R1" "$R2"

fastp \
    --in1 "$R1" \
    --in2 "$R2" \
    --out1 "$OUT/${SAMPLE}_R1.fastq.gz" \
    --out2 "$OUT/${SAMPLE}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --poly_g_min_len 10 \
    --disable_quality_filtering \
    --disable_length_filtering \
    --thread "$SLURM_CPUS_PER_TASK" \
    --json "$REPORT/${SAMPLE}.fastp.json" \
    --html "$REPORT/${SAMPLE}.fastp.html" \
    --report_title "${SAMPLE} fastp preprocessing"

echo "===== COMPLETE ====="
date
