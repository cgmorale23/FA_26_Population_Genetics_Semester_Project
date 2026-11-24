#!/usr/bin/env bash
#SBATCH --job-name=fpmerge_resume
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/fastp_merge_exploration/resume_%j.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/fastp_merge_exploration/resume_%j.err

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
RUN="$ROOT/qc/roh/fastp_merge_exploration/20261121_123041"
AUTO="$RUN/autosomes"
OUT="$RUN/pairwise_resume"
ROH_OUT="$RUN/roh/pairwise_shared"
PLINK_BIN="/N/u/xta1/Quartz/tools/plink_1.9_20241022/plink"

mkdir -p "$OUT" "$ROH_OUT"

module purge
module load gnu/14.2.0
module load bcftools/1.17
module load python/3.12.4

if [[ ! -s "$RUN/compare_shared_gt.py" ]]; then
    echo "ERROR: missing comparison helper: $RUN/compare_shared_gt.py" >&2
    exit 1
fi

if [[ ! -x "$PLINK_BIN" ]]; then
    echo "ERROR: PLINK executable not found: $PLINK_BIN" >&2
    exit 1
fi

printf 'comparison\tA_records\tB_records\tshared_records\tA_private\tB_private\n' \
    > "$OUT/pairwise_variant_counts.tsv"

run_roh() {
    local label="$1"
    local source_bcf="$2"
    local dest="$3"
    local prefix="$dest/$label"
    local header="$dest/$label.header"
    local fixed="$dest/$label.header_fixed.bcf"
    local renamed="$dest/$label.autosomes_renamed.bcf"
    local map="$RUN/chr_rename.tsv"

    mkdir -p "$dest"

    if [[ ! -s "$map" ]]; then
        : > "$map"
        for chr in $(seq 1 27); do
            printf 'chr%s\t%s\n' "$chr" "$chr" >> "$map"
        done
    fi

    bcftools view -h "$source_bcf" |
        awk '
            /^##contig=/ {
                if ($0 ~ /ID=chr([1-9]|1[0-9]|2[0-7])[,>]/) print
                next
            }
            { print }
        ' > "$header"

    bcftools reheader -h "$header" -o "$fixed" "$source_bcf"
    bcftools annotate --rename-chrs "$map" \
        -Ob -o "$renamed" "$fixed"

    "$PLINK_BIN" --bcf "$renamed" \
        --chr-set 27 \
        --double-id \
        --memory 12000 \
        --make-bed \
        --out "$prefix"

    "$PLINK_BIN" --bfile "$prefix" \
        --chr-set 27 \
        --memory 12000 \
        --homozyg \
        --homozyg-window-het 3 \
        --homozyg-window-missing 20 \
        --out "$prefix.roh"
}

compare_pair() {
    local a="$1"
    local b="$2"
    local name="${a}_vs_${b}"
    local dir="$OUT/$name"
    local isec="$dir/isec"
    local common_a="$isec/0000.bcf"
    local common_b="$isec/0001.bcf"

    mkdir -p "$isec"

    bcftools isec -p "$isec" -n=2 -c none -Ob \
        "$AUTO/$a.autosomes.bcf" "$AUTO/$b.autosomes.bcf"

    for f in "$common_a" "$common_b"; do
        if [[ ! -s "$f" ]]; then
            echo "ERROR: shared-site BCF not found: $f" >&2
            [[ -f "$isec/README.txt" ]] && cat "$isec/README.txt" >&2
            exit 1
        fi
    done

    local n_a n_b n_shared n_a_private n_b_private
    n_a=$(bcftools view -H "$AUTO/$a.autosomes.bcf" | wc -l)
    n_b=$(bcftools view -H "$AUTO/$b.autosomes.bcf" | wc -l)
    n_shared=$(bcftools view -H "$common_a" | wc -l)
    n_a_private=$((n_a - n_shared))
    n_b_private=$((n_b - n_shared))

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$n_a" "$n_b" "$n_shared" \
        "$n_a_private" "$n_b_private" \
        >> "$OUT/pairwise_variant_counts.tsv"

    python3 "$RUN/compare_shared_gt.py" \
        "$common_a" "$common_b" "$a" "$b" "$dir/$name"

    run_roh "$a" "$common_a" "$ROH_OUT/$name"
    run_roh "$b" "$common_b" "$ROH_OUT/$name"
}

echo "===== Pairwise comparisons ====="
compare_pair current_fastp fastp_merge_paired_only
compare_pair fastp_merge_paired_only fastp_merge
compare_pair current_fastp fastp_merge

echo "===== Collect full-site and shared-site ROH summaries ====="
printf 'dataset\tsample\tROH_segments\ttotal_ROH_kb\tmean_ROH_kb\n' \
    > "$OUT/roh_by_sample.tsv"

find "$RUN/roh/full" "$ROH_OUT" -type f -name '*.roh.hom.indiv' -print | sort |
while read -r file; do
    dataset="${file##*/}"
    dataset="${dataset%.roh.hom.indiv}"
    awk -v d="$dataset" '
        NR == 1 { next }
        NF >= 6 { print d "\t" $2 "\t" $4 "\t" $5 "\t" $6 }
    ' "$file" >> "$OUT/roh_by_sample.tsv"
done

echo "Completed pairwise continuation."
echo "Results: $OUT"
