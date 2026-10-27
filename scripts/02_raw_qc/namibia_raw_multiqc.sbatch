#!/bin/bash
#SBATCH --job-name=namibia_raw_multiqc
#SBATCH -A c02431
#SBATCH --output=slurm/logs/namibia_raw_multiqc_%j.out
#SBATCH --error=slurm/logs/namibia_raw_multiqc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=general

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
MULTIQC="/N/u/xta1/Quartz/venvs/multiqc/bin/multiqc"

INPUT="$ROOT/qc/namibia_raw_fastqc"
OUT="$ROOT/qc/namibia_raw_multiqc"

mkdir -p "$OUT"

"$MULTIQC" \
    "$INPUT" \
    -o "$OUT" \
    --filename namibia_raw_multiqc.html \
    --force
