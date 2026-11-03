#!/bin/bash
#SBATCH --job-name=fastp_merge_on
#SBATCH -A c02431
#SBATCH --output=slurm/logs/fastp_merge_on_%A_%a.out
#SBATCH --error=slurm/logs/fastp_merge_on_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --array=1-54%12

set -euo pipefail

module purge
module load fastp/0.23.4

ROOT="/N/scratch/xta1/ZambiaElephants"
MANIFEST="$ROOT/qc/fastp_merge/all_54_raw_pairs.tsv"

OUT="$ROOT/processed/fastp_merge_on"
REPORT="$ROOT/qc/fastp_merge/reports"

cd "$ROOT"

mkdir -p \
    "$OUT/public" \
    "$OUT/namibia" \
    "$OUT/zambia" \
    "$REPORT/public" \
    "$REPORT/namibia" \
    "$REPORT/zambia"

LINE=$((SLURM_ARRAY_TASK_ID + 1))

IFS=$'\t' read -r DATASET UNIT R1 R2 < <(
    sed -n "${LINE}p" "$MANIFEST"
)

if [[ -z "${DATASET:-}" || -z "${UNIT:-}" || ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: invalid input"
    echo "DATASET=$DATASET"
    echo "UNIT=$UNIT"
    echo "R1=$R1"
    echo "R2=$R2"
    exit 1
fi

OUTDIR="$OUT/$DATASET"
REPORTDIR="$REPORT/$DATASET"

echo "===== SOFTWARE ====="
command -v fastp
fastp --version

echo
echo "===== JOB ====="
echo "Job ID:      $SLURM_JOB_ID"
echo "Array task:  $SLURM_ARRAY_TASK_ID"
echo "Dataset:     $DATASET"
echo "Unit:        $UNIT"
echo "Host:        $(hostname)"
echo "Start:       $(date)"

echo
echo "===== INPUT ====="
echo "R1=$R1"
echo "R2=$R2"
ls -lh "$R1" "$R2"

echo
echo "===== FASTP: MERGE ON ====="

fastp \
    --in1 "$R1" \
    --in2 "$R2" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --poly_g_min_len 10 \
    --disable_quality_filtering \
    --disable_length_filtering \
    --merge \
    --merged_out "$OUTDIR/${UNIT}.merged.fastq.gz" \
    --out1 "$OUTDIR/${UNIT}.unmerged_R1.fastq.gz" \
    --out2 "$OUTDIR/${UNIT}.unmerged_R2.fastq.gz" \
    --thread "$SLURM_CPUS_PER_TASK" \
    --json "$REPORTDIR/${UNIT}.fastp.json" \
    --html "$REPORTDIR/${UNIT}.fastp.html" \
    --report_title "${UNIT} fastp merge-on preprocessing"

echo
echo "===== COMPLETE ====="
echo "Dataset: $DATASET"
echo "Unit:    $UNIT"
echo "End:     $(date)"
