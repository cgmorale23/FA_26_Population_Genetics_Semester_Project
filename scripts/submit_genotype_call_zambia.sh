#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/N/scratch/xta1/ZambiaElephants}"
SCRIPT="${SCRIPT:-$ROOT/slurm/scripts/genotype_call_zambia_array.sbatch}"
BAM_ROOT="${BAM_ROOT:-$ROOT/mapping/sensitivity/duplicate_cleanup_array/20261109_182300/outputs/pecnerova_zambia_dedup_removed}"
REF="${REF:-$ROOT/rawdata/reference/LoxAfr4/LoxAfr4.fa}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="$ROOT/qc/genotypes/zambia_pecnerova_track/$RUN_ID"
LOG_DIR="$ROOT/slurm/logs/genotype_calling"

if [[ ! -r "$SCRIPT" || ! -r "$REF" || ! -r "$REF.fai" ]]; then
    echo "ERROR: missing script, reference, or reference index." >&2
    echo "Script: $SCRIPT" >&2
    echo "Reference: $REF" >&2
    exit 2
fi

for branch in adapterremoval current_fastp fastp_merge; do
    case "$branch" in
        adapterremoval) read_types=(paired collapsed) ;;
        current_fastp) read_types=(paired) ;;
        fastp_merge) read_types=(paired merged) ;;
    esac
    for sample_number in {1..8}; do
        sample="Zambia$(printf '%02d' "$sample_number")"
        for read_type in "${read_types[@]}"; do
            bam="$BAM_ROOT/$branch/dup_cleanup/${sample}.${read_type}.Library1.dedup.bam"
            if [[ ! -s "$bam" || ! -s "$bam.bai" ]]; then
                echo "ERROR: missing BAM or index: $bam" >&2
                exit 2
            fi
        done
    done
done

mkdir -p "$RUN_DIR" "$LOG_DIR"
echo "Run ID: $RUN_ID"
echo "BAM root: $BAM_ROOT"
echo "Reference: $REF"
echo "Run directory: $RUN_DIR"
echo "Array: 81 chromosome tasks (3 branches x 27 autosomes), concurrency 6"

sbatch --export=ALL,ROOT="$ROOT",BAM_ROOT="$BAM_ROOT",REF="$REF",RUN_DIR="$RUN_DIR" "$SCRIPT"

echo "After all 81 tasks complete, aggregate each branch with:"
echo "bash $ROOT/slurm/scripts/aggregate_genotype_call_zambia.sh $RUN_DIR"
