#!/bin/bash
#SBATCH --job-name=multiqc_public
#SBATCH -A c02431
#SBATCH --output=slurm/logs/multiqc_public_%j.out
#SBATCH --error=slurm/logs/multiqc_public_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=general

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
MULTIQC="/N/u/xta1/Quartz/venvs/multiqc/bin/multiqc"

INPUT="$ROOT/qc/fastqc_public_selected"
OUT="$ROOT/qc/multiqc_public_selected"

cd "$ROOT"
mkdir -p "$OUT"

echo "===== SOFTWARE ====="
"$MULTIQC" --version

echo
echo "===== INPUT ====="

N=$(find "$INPUT" -type f -name "*_fastqc.zip" | wc -l)

echo "FastQC reports: $N"

if [[ "$N" -ne 84 ]]; then
    echo "ERROR: expected 84 FastQC reports, found $N"
    exit 1
fi

echo
echo "===== MULTIQC ====="

"$MULTIQC" \
    "$INPUT" \
    -o "$OUT" \
    --force

echo
echo "===== OUTPUT ====="
find "$OUT" -maxdepth 2 -type f -print

echo
echo "Completed: $(date)"
