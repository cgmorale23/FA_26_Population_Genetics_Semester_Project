#!/bin/bash
#SBATCH --job-name=zambia_fastp_test
#SBATCH -A c02431
#SBATCH --output=slurm/logs/zambia_fastp_test_%A_%a.out
#SBATCH --error=slurm/logs/zambia_fastp_test_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --array=1-8%8

set -euo pipefail

module purge
module load fastp/0.23.4

ROOT="/N/scratch/xta1/ZambiaElephants"
MANIFEST="$ROOT/pilot_trimming/zambia_all/zambia_pairs.tsv"
OUT="$ROOT/pilot_trimming/zambia_all/trimmed"
REPORT="$ROOT/pilot_trimming/zambia_all/fastp_reports"

mkdir -p "$OUT" "$REPORT"

echo "===== SOFTWARE ====="
command -v fastp
fastp --version

# +1 because line 1 is the header
LINE=$((SLURM_ARRAY_TASK_ID + 1))

IFS=$'\t' read -r SAMPLE R1 R2 < <(sed -n "${LINE}p" "$MANIFEST")

if [[ -z "${SAMPLE:-}" || ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: invalid input for task $SLURM_ARRAY_TASK_ID"
    echo "SAMPLE=$SAMPLE"
    echo "R1=$R1"
    echo "R2=$R2"
    exit 1
fi

echo "===== SAMPLE ====="
echo "Sample: $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"
echo "Start: $(date)"

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
    --report_title "${SAMPLE} conservative trimming pilot"

echo "===== COMPLETE ====="
echo "End: $(date)"
