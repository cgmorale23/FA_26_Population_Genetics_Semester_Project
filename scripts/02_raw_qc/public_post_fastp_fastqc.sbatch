#!/bin/bash
#SBATCH --job-name=public_post_fastp_qc
#SBATCH -A c02431
#SBATCH --output=slurm/logs/public_post_fastp_qc_%A_%a.out
#SBATCH --error=slurm/logs/public_post_fastp_qc_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --array=1-84%20

set -euo pipefail

module purge
module load fastqc/0.12.1

ROOT="/N/scratch/xta1/ZambiaElephants"
LIST="$ROOT/qc/public_post_fastp_fastq.list"
OUT="$ROOT/qc/public_post_fastp_fastqc"

cd "$ROOT"
mkdir -p "$OUT"

echo "===== SOFTWARE ====="
command -v fastqc
fastqc --version

# One post-fastp FASTQ per array task
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [[ -z "${FASTQ:-}" ]]; then
    echo "ERROR: no FASTQ for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

if [[ ! -s "$FASTQ" ]]; then
    echo "ERROR: FASTQ missing or empty:"
    echo "$FASTQ"
    exit 1
fi

FILE=$(basename "$FASTQ")

echo
echo "===== JOB ====="
echo "Job ID:       $SLURM_JOB_ID"
echo "Array task:   $SLURM_ARRAY_TASK_ID"
echo "Host:         $(hostname)"
echo "FASTQ:        $FILE"
echo "Start:        $(date)"

echo
echo "===== INPUT ====="
ls -lh "$FASTQ"

echo
echo "===== FASTQC ====="

fastqc \
    --threads "$SLURM_CPUS_PER_TASK" \
    --outdir "$OUT" \
    "$FASTQ"

echo
echo "===== COMPLETE ====="
echo "FASTQ: $FILE"
echo "End:   $(date)"
