#!/bin/bash
#SBATCH --job-name=zambia_multiqc
#SBATCH -A c02431
#SBATCH --output=slurm/logs/zambia_multiqc_%j.out
#SBATCH --error=slurm/logs/zambia_multiqc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=general

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"

MULTIQC="/N/u/xta1/Quartz/venvs/multiqc/bin/multiqc"

INPUT="$ROOT/qc/zambia_samples_fastqc"
OUT="$ROOT/qc/zambia_samples_multiqc"

cd "$ROOT"

mkdir -p "$OUT"

echo "===== SOFTWARE ====="
"$MULTIQC" --version

echo
echo "===== INPUT CHECK ====="

N_ZIP=$(find "$INPUT" \
    -maxdepth 1 \
    -type f \
    -name "*_fastqc.zip" \
    | wc -l)

N_HTML=$(find "$INPUT" \
    -maxdepth 1 \
    -type f \
    -name "*_fastqc.html" \
    | wc -l)

echo "FastQC ZIP reports:  $N_ZIP"
echo "FastQC HTML reports: $N_HTML"

if [[ "$N_ZIP" -ne 16 ]]; then
    echo "ERROR: expected 16 FastQC ZIP reports, found $N_ZIP"
    exit 1
fi

echo
echo "===== MULTIQC ====="

"$MULTIQC" \
    "$INPUT" \
    -o "$OUT" \
    --filename zambia_samples_multiqc.html \
    --force

echo
echo "===== OUTPUT ====="

find "$OUT" -maxdepth 2 -type f -print

echo
echo "===== COMPLETE ====="
echo "End: $(date)"
