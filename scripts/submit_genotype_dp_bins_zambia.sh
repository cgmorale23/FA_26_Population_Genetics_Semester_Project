#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/N/scratch/xta1/ZambiaElephants}"
INPUT_RUN_DIR="${INPUT_RUN_DIR:-$ROOT/qc/genotypes/zambia_pecnerova_track/20261115_013407}"
SBATCH_SCRIPT="${SBATCH_SCRIPT:-$ROOT/slurm/scripts/genotype_dp_bins_zambia.sbatch}"
WORKER="${WORKER:-$ROOT/slurm/scripts/genotype_dp_bins_zambia.py}"
DP_RUN_ID="${DP_RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
DP_RUN_DIR="$INPUT_RUN_DIR/genotype_qc/dp_bins_$DP_RUN_ID"
LOG_DIR="$ROOT/slurm/logs/genotype_dp_bins/$DP_RUN_ID"

if [[ ! -r "$SBATCH_SCRIPT" || ! -r "$WORKER" ]]; then
    echo "ERROR: missing DP-bin scripts." >&2
    echo "SBATCH script: $SBATCH_SCRIPT" >&2
    echo "Worker: $WORKER" >&2
    exit 2
fi

for branch in adapterremoval current_fastp fastp_merge; do
    bcf="$INPUT_RUN_DIR/$branch/zambia_${branch}.raw_biallelic_snps.bcf"
    if [[ ! -s "$bcf" || ! -s "$bcf.csi" ]]; then
        echo "ERROR: missing BCF or CSI index: $bcf" >&2
        exit 2
    fi
done

mkdir -p "$DP_RUN_DIR" "$LOG_DIR"
echo "Input run: $INPUT_RUN_DIR"
echo "DP-bin output: $DP_RUN_DIR"
echo "Logs: $LOG_DIR"
echo "Array: 3 tasks, one per preprocessing branch"

sbatch \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --export=ALL,ROOT="$ROOT",INPUT_RUN_DIR="$INPUT_RUN_DIR",QC_RUN_DIR="$DP_RUN_DIR",SCRIPT="$WORKER" \
    "$SBATCH_SCRIPT"

echo "After all three tasks complete, run:"
echo "bash $ROOT/slurm/scripts/aggregate_genotype_dp_bins_zambia.sh $DP_RUN_DIR"
