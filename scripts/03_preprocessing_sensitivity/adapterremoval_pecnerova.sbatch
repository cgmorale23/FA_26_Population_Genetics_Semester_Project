#!/bin/bash
#SBATCH --job-name=adapterremoval_pec
#SBATCH -A c02431
#SBATCH --output=slurm/logs/adapterremoval_pec_%A_%a.out
#SBATCH --error=slurm/logs/adapterremoval_pec_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --array=1-54%12

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"

AR="/N/u/xta1/Quartz/software/AdapterRemoval-2.3.1/build/AdapterRemoval"

MANIFEST="$ROOT/qc/adapterremoval_pecnerova/all_54_pairs.tsv"

OUT="$ROOT/processed/adapterremoval_pecnerova"
SETTINGS="$ROOT/qc/adapterremoval_pecnerova/settings"

cd "$ROOT"

mkdir -p "$OUT/public" \
         "$OUT/namibia" \
         "$OUT/zambia" \
         "$SETTINGS/public" \
         "$SETTINGS/namibia" \
         "$SETTINGS/zambia"

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

PREFIX="$OUT/$DATASET/$UNIT"
SETTING="$SETTINGS/$DATASET/${UNIT}.settings"

echo "===== SOFTWARE ====="
"$AR" --version

echo
echo "===== JOB ====="
echo "Job ID:      $SLURM_JOB_ID"
echo "Array task:  $SLURM_ARRAY_TASK_ID"
echo "Host:        $(hostname)"
echo "Dataset:     $DATASET"
echo "Unit:        $UNIT"
echo "Start:       $(date)"

echo
echo "===== INPUT ====="
echo "R1=$R1"
echo "R2=$R2"
ls -lh "$R1" "$R2"

echo
echo "===== ADAPTERREMOVAL ====="

"$AR" \
    --file1 "$R1" \
    --file2 "$R2" \
    --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --mm 3 \
    --minlength 1 \
    --collapse \
    --qualitybase 33 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --gzip \
    --basename "$PREFIX" \
    --settings "$SETTING"

echo
echo "===== COMPLETE ====="
echo "Dataset: $DATASET"
echo "Unit:    $UNIT"
echo "End:     $(date)"
