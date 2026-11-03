#!/bin/bash
#SBATCH --job-name=map_cleanup
#SBATCH -A c02431
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/map_cleanup_%A_%a.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/map_cleanup_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=16:00:00
#SBATCH --partition=general

set -euo pipefail

# ============================================================
# ZambiaElephants preprocessing sensitivity mapping
#
# Usage:
#
# CURRENT FASTP:
# sbatch \
#   --array=1-54%10 \
#   --export=ALL,BRANCH=current_fastp \
#   /N/scratch/xta1/ZambiaElephants/slurm/scripts/mapping_align_cleanup.sbatch
#
# FASTP MERGE ON:
# sbatch \
#   --array=1-108%10 \
#   --export=ALL,BRANCH=fastp_merge \
#   /N/scratch/xta1/ZambiaElephants/slurm/scripts/mapping_align_cleanup.sbatch
#
# ADAPTERREMOVAL:
# sbatch \
#   --array=1-108%10 \
#   --export=ALL,BRANCH=adapterremoval \
#   /N/scratch/xta1/ZambiaElephants/slurm/scripts/mapping_align_cleanup.sbatch
#
# For a single-task test:
#
# sbatch \
#   --array=1 \
#   --export=ALL,BRANCH=current_fastp \
#   /N/scratch/xta1/ZambiaElephants/slurm/scripts/mapping_align_cleanup.sbatch
# ============================================================


# ============================================================
# SOFTWARE
# ============================================================

module purge
module load gnu/14.2.0
module load python/3.12.4
module load bwa/0.7.17
module load samtools/1.17


# ============================================================
# PATHS
# ============================================================

ROOT="/N/scratch/xta1/ZambiaElephants"

REF="$ROOT/rawdata/reference/LoxAfr4/LoxAfr4.fa"

PALEOMIX="/N/u/xta1/Quartz/venvs/paleomix_pub2022/bin/paleomix"

THREADS="${SLURM_CPUS_PER_TASK}"

cd "$ROOT"


# ============================================================
# SELECT INPUT MANIFEST
# ============================================================

case "${BRANCH:-}" in

    current_fastp)
        MANIFEST="$ROOT/qc/mapping/current_fastp_mapping.tsv"
        ;;

    fastp_merge)
        MANIFEST="$ROOT/qc/mapping/fastp_merge_mapping.tsv"
        ;;

    adapterremoval)
        MANIFEST="$ROOT/qc/mapping/adapterremoval_mapping.tsv"
        ;;

    *)
        echo "ERROR: BRANCH must be one of:"
        echo "  current_fastp"
        echo "  fastp_merge"
        echo "  adapterremoval"
        exit 1
        ;;
esac


# ============================================================
# CHECK REQUIRED FILES
# ============================================================

if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: manifest does not exist or is empty:"
    echo "$MANIFEST"
    exit 1
fi

if [[ ! -s "$REF" ]]; then
    echo "ERROR: reference FASTA missing:"
    echo "$REF"
    exit 1
fi

if [[ ! -x "$PALEOMIX" ]]; then
    echo "ERROR: PALEOMIX executable missing:"
    echo "$PALEOMIX"
    exit 1
fi


# ============================================================
# READ MANIFEST ROW
#
# Manifest contains one header row, so:
# array task 1 -> line 2
# ============================================================

LINE=$((SLURM_ARRAY_TASK_ID + 1))

IFS=$'\t' read -r \
    BRANCH_IN \
    SOURCE_DATASET \
    ANALYSIS_GROUP \
    SAMPLE \
    LIBRARY \
    PU \
    UNIT \
    READ_TYPE \
    R1 \
    R2 \
    < <(sed -n "${LINE}p" "$MANIFEST")


# ============================================================
# VALIDATE MANIFEST ROW
# ============================================================

if [[ -z "${BRANCH_IN:-}" ]]; then
    echo "ERROR: empty branch in manifest row $LINE"
    exit 1
fi

if [[ -z "${SAMPLE:-}" ]]; then
    echo "ERROR: empty sample in manifest row $LINE"
    exit 1
fi

if [[ -z "${LIBRARY:-}" ]]; then
    echo "ERROR: empty library in manifest row $LINE"
    exit 1
fi

if [[ -z "${PU:-}" ]]; then
    echo "ERROR: empty PU in manifest row $LINE"
    exit 1
fi

if [[ -z "${UNIT:-}" ]]; then
    echo "ERROR: empty unit in manifest row $LINE"
    exit 1
fi

if [[ -z "${READ_TYPE:-}" ]]; then
    echo "ERROR: empty read type in manifest row $LINE"
    exit 1
fi

if [[ -z "${R1:-}" || ! -s "$R1" ]]; then
    echo "ERROR: R1 missing or empty:"
    echo "${R1:-UNDEFINED}"
    exit 1
fi

if [[ "$READ_TYPE" == "paired" ]]; then

    if [[ -z "${R2:-}" || "$R2" == "." || ! -s "$R2" ]]; then
        echo "ERROR: paired record has missing or empty R2:"
        echo "${R2:-UNDEFINED}"
        exit 1
    fi

fi


# ============================================================
# OUTPUT PATHS
# ============================================================

OUTDIR="$ROOT/mapping/sensitivity/$BRANCH_IN/unit_bams/$ANALYSIS_GROUP"

TMPDIR="$ROOT/mapping/sensitivity/$BRANCH_IN/tmp"

mkdir -p "$OUTDIR"
mkdir -p "$TMPDIR"

OUT="$OUTDIR/${UNIT}.${READ_TYPE}.cleanup.bam"

TMP="$TMPDIR/${UNIT}.${READ_TYPE}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"


# ============================================================
# READ GROUP
#
# PALEOMIX hierarchy:
#
# ID = library
# SM = biological sample
# LB = library
# PU = sequencing unit
# PL = ILLUMINA
#
# --rg-id supplies ID separately.
# Each --rg argument therefore uses TAG:VALUE.
# ============================================================

RGID="$LIBRARY"


# ============================================================
# JOB INFORMATION
# ============================================================

echo "============================================================"
echo "MAPPING + PALEOMIX CLEANUP"
echo "============================================================"
echo "Branch:          $BRANCH_IN"
echo "Source dataset:  $SOURCE_DATASET"
echo "Analysis group:  $ANALYSIS_GROUP"
echo "Sample:          $SAMPLE"
echo "Library:         $LIBRARY"
echo "PU:              $PU"
echo "Unit:            $UNIT"
echo "Read type:       $READ_TYPE"
echo "R1:              $R1"
echo "R2:              $R2"
echo "RG ID:           $RGID"
echo "Reference:       $REF"
echo "Output:          $OUT"
echo "Temporary prefix:$TMP"
echo "Host:            $(hostname)"
echo "Start:           $(date)"
echo "============================================================"


# ============================================================
# SOFTWARE VERSIONS
# ============================================================

echo
echo "===== SOFTWARE ====="

echo
echo "BWA:"
bwa 2>&1 | head -3 || true

echo
echo "SAMtools:"
samtools --version | head -2

echo
echo "PALEOMIX:"
"$PALEOMIX" --version


# ============================================================
# REMOVE POSSIBLE INCOMPLETE OUTPUT FROM A PREVIOUS FAILED RUN
# ============================================================

rm -f "$OUT"
rm -f "${OUT}.bai"


# ============================================================
# ALIGNMENT + PALEOMIX CLEANUP
#
# PAIRED:
#
# bwa mem
#   -> samtools fixmate -m
#   -> PALEOMIX cleanup
#   -> coordinate sort
#   -> samtools calmd
#
# MERGED / COLLAPSED / SINGLETON:
#
# bwa mem (single-end)
#   -> PALEOMIX cleanup
#   -> coordinate sort
#   -> samtools calmd
# ============================================================

echo
echo "===== ALIGNMENT ====="


if [[ "$READ_TYPE" == "paired" ]]; then

    echo "Running paired-end BWA-MEM + PALEOMIX cleanup"

    bwa mem \
        -t "$THREADS" \
        "$REF" \
        "$R1" \
        "$R2" \
    | "$PALEOMIX" cleanup \
        --fasta "$REF" \
        --temp-prefix "$TMP" \
        --max-threads "$THREADS" \
        --paired-end \
        --add-mate-score \
        --min-quality 0 \
        --exclude-flags 0x0 \
        --rg-id "$RGID" \
        --rg "SM:${SAMPLE}" \
        --rg "LB:${LIBRARY}" \
        --rg "PU:${PU}" \
        --rg "PL:ILLUMINA" \
        > "$OUT"


elif [[ "$READ_TYPE" == "merged" || \
        "$READ_TYPE" == "collapsed" || \
        "$READ_TYPE" == "singleton" ]]; then

    echo "Running single-end BWA-MEM + PALEOMIX cleanup"

    bwa mem \
        -t "$THREADS" \
        "$REF" \
        "$R1" \
    | "$PALEOMIX" cleanup \
        --fasta "$REF" \
        --temp-prefix "$TMP" \
        --max-threads "$THREADS" \
        --min-quality 0 \
        --exclude-flags 0x0 \
        --rg-id "$RGID" \
        --rg "SM:${SAMPLE}" \
        --rg "LB:${LIBRARY}" \
        --rg "PU:${PU}" \
        --rg "PL:ILLUMINA" \
        > "$OUT"


else

    echo "ERROR: unsupported READ_TYPE:"
    echo "$READ_TYPE"
    exit 1

fi


# ============================================================
# BAM VALIDATION
# ============================================================

echo
echo "===== BAM VALIDATION ====="

if [[ ! -s "$OUT" ]]; then
    echo "ERROR: output BAM is missing or empty:"
    echo "$OUT"
    exit 1
fi

samtools quickcheck -v "$OUT"

echo "samtools quickcheck: PASS"


# ============================================================
# BAM INDEX
# ============================================================

echo
echo "===== BAM INDEX ====="

samtools index \
    -@ "$THREADS" \
    "$OUT"

if [[ ! -s "${OUT}.bai" ]]; then
    echo "ERROR: BAM index was not created:"
    echo "${OUT}.bai"
    exit 1
fi

echo "BAM index: PASS"


# ============================================================
# HEADER / READ GROUP CHECK
# ============================================================

echo
echo "===== BAM HEADER / READ GROUP ====="

samtools view -H "$OUT" \
    | grep -E '^@HD|^@RG' || true


# ============================================================
# FLAGSTAT
# ============================================================

echo
echo "===== FLAGSTAT ====="

samtools flagstat \
    -@ "$THREADS" \
    "$OUT"


# ============================================================
# COMPLETE
# ============================================================

echo
echo "============================================================"
echo "COMPLETE"
echo "Branch:  $BRANCH_IN"
echo "Sample:  $SAMPLE"
echo "Unit:    $UNIT"
echo "Type:    $READ_TYPE"
echo "Output:  $OUT"
echo "End:     $(date)"
echo "============================================================"
