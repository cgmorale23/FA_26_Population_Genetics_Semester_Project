#!/usr/bin/env bash
#SBATCH --job-name=zam_genotype
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-81%6
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/genotype_calling/%x_%A_%a.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/genotype_calling/%x_%A_%a.err

set -euo pipefail

module purge
module load gnu/14.2.0
module load samtools/1.17
module load bcftools/1.17

ROOT="${ROOT:-/N/scratch/xta1/ZambiaElephants}"
BAM_ROOT="${BAM_ROOT:-$ROOT/mapping/sensitivity/duplicate_cleanup_array/20261109_182300/outputs/pecnerova_zambia_dedup_removed}"
REF="${REF:-$ROOT/rawdata/reference/LoxAfr4/LoxAfr4.fa}"
RUN_DIR="${RUN_DIR:?Set RUN_DIR when submitting the array}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?This script must run as a SLURM array task}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

branches=(adapterremoval current_fastp fastp_merge)
branch_index=$(( (TASK_ID - 1) / 27 ))
chromosome_number=$(( (TASK_ID - 1) % 27 + 1 ))
if (( branch_index < 0 || branch_index >= ${#branches[@]} )); then
    echo "ERROR: array task ID is outside 1-81: $TASK_ID" >&2
    exit 2
fi
BRANCH="${branches[$branch_index]}"
CHROMOSOME="chr${chromosome_number}"
BRANCH_DIR="$RUN_DIR/$BRANCH"
TASK_DIR="$BRANCH_DIR/chromosomes"
RG_MAP="$TASK_DIR/${CHROMOSOME}.read_groups.tsv"
OUTPUT="$TASK_DIR/${CHROMOSOME}.snps.bcf"
PARTIAL="$OUTPUT.partial.${SLURM_JOB_ID}_${TASK_ID}"

mkdir -p "$TASK_DIR"
if [[ ! -s "$REF" || ! -s "$REF.fai" ]]; then
    echo "ERROR: missing reference or FASTA index: $REF" >&2
    exit 2
fi

case "$BRANCH" in
    adapterremoval) read_types=(paired collapsed) ;;
    current_fastp) read_types=(paired) ;;
    fastp_merge) read_types=(paired merged) ;;
    *) echo "ERROR: unsupported branch: $BRANCH" >&2; exit 2 ;;
esac

bams=()
: > "$RG_MAP"
for sample_number in {1..8}; do
    sample="Zambia$(printf '%02d' "$sample_number")"
    for read_type in "${read_types[@]}"; do
        bam="$BAM_ROOT/$BRANCH/dup_cleanup/${sample}.${read_type}.Library1.dedup.bam"
        if [[ ! -s "$bam" || ! -s "$bam.bai" ]]; then
            echo "ERROR: missing BAM or index: $bam" >&2
            exit 2
        fi
        bams+=("$bam")
        header=$(samtools view -H "$bam")
        while IFS= read -r rg_line; do
            [[ "$rg_line" == @RG$'\t'* ]] || continue
            rg_id=$(awk -F '\t' '{for (i=2; i<=NF; i++) if ($i ~ /^ID:/) {sub(/^ID:/, "", $i); print $i; exit}}' <<< "$rg_line")
            if [[ -z "$rg_id" ]]; then
                echo "ERROR: read group without ID in $bam: $rg_line" >&2
                exit 2
            fi
            printf '%s\t%s\t%s\n' "$rg_id" "$bam" "$sample" >> "$RG_MAP"
        done <<< "$header"
    done
done

if [[ ! -s "$RG_MAP" ]]; then
    echo "ERROR: no @RG entries found for $BRANCH" >&2
    exit 2
fi

if [[ -e "$OUTPUT" || -e "$OUTPUT.csi" ]]; then
    echo "ERROR: refusing to overwrite existing output: $OUTPUT" >&2
    exit 3
fi
trap 'rm -f -- "$PARTIAL" "$PARTIAL.csi"' EXIT

echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Branch: $BRANCH"
echo "Chromosome: $CHROMOSOME"
echo "Input BAM count: ${#bams[@]}"
echo "Sample count: 8"
echo "Reference: $REF"
echo "BCFtools: $(bcftools --version | head -n 1)"
echo "Samtools: $(samtools --version | head -n 1)"

bcftools mpileup \
    --threads "$THREADS" \
    --fasta-ref "$REF" \
    --regions "$CHROMOSOME" \
    --min-MQ 30 \
    --min-BQ 25 \
    --per-sample-mF \
    --annotate FORMAT/DP,FORMAT/AD \
    --read-groups "$RG_MAP" \
    --output-type u \
    "${bams[@]}" \
| bcftools call --threads "$THREADS" --multiallelic-caller --variants-only --output-type u \
| bcftools view --threads "$THREADS" --types snps --min-alleles 2 --max-alleles 2 --output-type b --output "$PARTIAL"

bcftools index --force "$PARTIAL"
mv -- "$PARTIAL" "$OUTPUT"
mv -- "$PARTIAL.csi" "$OUTPUT.csi"
trap - EXIT

echo "Completed: $OUTPUT"
echo "Records: $(bcftools index --nrecords "$OUTPUT")"
echo "Samples: $(bcftools query --list-samples "$OUTPUT" | paste -sd ',' -)"
