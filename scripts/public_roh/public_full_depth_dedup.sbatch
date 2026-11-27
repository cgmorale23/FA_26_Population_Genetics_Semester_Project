#!/usr/bin/env bash
#SBATCH --job-name=pub_full_dedup
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --array=1-110%6
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/public_full_depth/full_%A_%a.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/public_full_depth/full_%A_%a.err

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
PLAN="$ROOT/qc/duplicate/array_runs/20261109_182300/plans/groups.pecnerova.jsonl"

MARKED_ROOT="$ROOT/mapping/sensitivity/duplicate_cleanup_array/20261109_182300/outputs"
OUT_ROOT="$ROOT/mapping/sensitivity/public_roh/full_depth/deduplicated"

module purge
module load gnu/14.2.0
module load samtools/1.17

TASK_ID="${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array task}"

ROW="$(python3 - "$PLAN" "$TASK_ID" <<'PY'
import json
import sys

plan_path = sys.argv[1]
task_id = int(sys.argv[2])

for line in open(plan_path):
    row = json.loads(line)
    if int(row["task_id"]) == task_id:
        print(
            row["track"],
            row["sample"],
            row["branch"],
            row["bam_type"],
            row["group"],
            sep="\t",
        )
        break
else:
    sys.exit(f"ERROR: task_id={task_id} not found in {plan_path}")
PY
)"

IFS=$'\t' read -r TRACK SAMPLE BRANCH BAM_TYPE GROUP <<< "$ROW"

# This workflow is for public samples; Zambia has its own existing pipeline.
if [[ "$SAMPLE" =~ ^Zambia0[1-8]$ ]]; then
    echo "SKIP task=$TASK_ID sample=$SAMPLE; Zambia sample"
    exit 0
fi

SOURCE="$MARKED_ROOT/$TRACK/$BRANCH/dup_cleanup/${SAMPLE}.${BAM_TYPE}.${GROUP}.dupmark.bam"
DEST_DIR="$OUT_ROOT/$TRACK/$BRANCH/$SAMPLE"
DEST="$DEST_DIR/${SAMPLE}.${BAM_TYPE}.${GROUP}.dedup.bam"
INDEX="$DEST.bai"

is_valid_bam() {
    local bam="$1"
    [[ -s "$bam" && -s "$bam.bai" ]] || return 1
    samtools quickcheck "$bam" >/dev/null 2>&1 || return 1
    samtools idxstats "$bam" >/dev/null 2>&1 || return 1
}

if [[ ! -s "$SOURCE" ]]; then
    echo "ERROR: marked input BAM missing or empty: $SOURCE" >&2
    exit 1
fi

if is_valid_bam "$DEST"; then
    echo "SKIP: valid output already exists: $DEST"
    exit 0
fi

if [[ -e "$DEST" || -e "$INDEX" ]]; then
    echo "ERROR: output exists but did not pass validation; refusing to overwrite: $DEST" >&2
    exit 1
fi

samtools quickcheck "$SOURCE"

mkdir -p "$DEST_DIR"
TMP="$DEST.partial.${SLURM_JOB_ID}"
trap 'rm -f "$TMP" "$TMP.bai"' EXIT

echo "START task=$TASK_ID sample=$SAMPLE branch=$BRANCH type=$BAM_TYPE"
echo "Input: $SOURCE"
echo "Output: $DEST"

# FLAG 0x400 marks duplicate alignments.
samtools view \
    -@ "$SLURM_CPUS_PER_TASK" \
    -b \
    -F 0x400 \
    -o "$TMP" \
    "$SOURCE"

samtools quickcheck "$TMP"
samtools index -@ 2 "$TMP"
samtools idxstats "$TMP" >/dev/null

mv -- "$TMP" "$DEST"
mv -- "$TMP.bai" "$INDEX"
trap - EXIT

samtools quickcheck "$DEST"
samtools idxstats "$DEST" >/dev/null

echo "DONE task=$TASK_ID sample=$SAMPLE output=$DEST"
