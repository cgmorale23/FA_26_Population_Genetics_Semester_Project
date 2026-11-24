#!/usr/bin/env bash
#SBATCH --job-name=roh_dp_maf
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/roh_dp_maf_sensitivity/sens_%j.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/roh_dp_maf_sensitivity/sens_%j.err

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
GENO="$ROOT/qc/genotypes"
BASELINE="$ROOT/qc/roh/hardcall_baseline/20261116_222813"
REF_FAI="$ROOT/rawdata/reference/LoxAfr4/LoxAfr4.fa.fai"
RUN_ID="$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual}"
OUT="$ROOT/qc/roh/dp_maf_sensitivity/$RUN_ID"
WORK="$OUT/work"
CHR_MAP="$WORK/chr_rename.tsv"
AUTOSOME_BED="$WORK/autosomes.bed"
METRICS="$OUT/new_scenario_metrics.tsv"
PLINK_BIN="/N/u/xta1/Quartz/tools/plink_1.9_20241022/plink"

mkdir -p "$WORK" "$OUT/roh"

module purge
module load gnu/14.2.0
module load bcftools/1.17
module load python/3.12.4

[[ -x "$PLINK_BIN" ]] || {
    echo "ERROR: PLINK executable not found: $PLINK_BIN" >&2
    exit 1
}
[[ -s "$REF_FAI" ]] || {
    echo "ERROR: reference FAI not found: $REF_FAI" >&2
    exit 1
}

bcftools --version | head -n 1
"$PLINK_BIN" --version | head -n 1

for chr in $(seq 1 27); do
    printf 'chr%s\t%s\n' "$chr" "$chr"
done > "$CHR_MAP"

awk 'BEGIN { OFS="\t" } $1 ~ /^chr([1-9]|1[0-9]|2[0-7])$/ { print $1, 0, $2 }' \
    "$REF_FAI" > "$AUTOSOME_BED"

bed_count=$(wc -l < "$AUTOSOME_BED")
if [[ "$bed_count" -ne 27 ]]; then
    echo "ERROR: expected 27 autosome intervals; found $bed_count" >&2
    cat "$AUTOSOME_BED" >&2
    exit 1
fi

printf 'branch\tdp_mask\tmaf_filter\tsites_retained\tcalled_genotypes\tmissing_genotypes\ttotal_genotype_slots\n' \
    > "$METRICS"

copy_existing_baseline() {
    local branch="$1"
    local source="$BASELINE/${branch}.roh.hom.indiv"
    local destination="$OUT/roh/$branch/none/none/roh.hom.indiv"

    [[ -s "$source" ]] || {
        echo "ERROR: baseline ROH report not found: $source" >&2
        exit 1
    }

    mkdir -p "$(dirname "$destination")"
    cp -- "$source" "$destination"
}

copy_existing_baseline adapterremoval
copy_existing_baseline current_fastp
copy_existing_baseline fastp_merge

validate_autosome_records() {
    local bcf="$1"
    local contigs="$WORK/$(basename "$bcf").record_contigs.tsv"

    # This checks contigs with variant records, not unused header declarations.
    bcftools index -s "$bcf" > "$contigs"

    if ! awk -F '\t' \
        '$1 !~ /^([1-9]|1[0-9]|2[0-7])$/ { print; bad=1 } END { exit bad }' \
        "$contigs"; then
        echo "ERROR: non-autosomal variant records remain in $bcf" >&2
        cat "$contigs" >&2
        exit 1
    fi

    echo "PASS: variant records are restricted to autosomes 1-27 in $bcf"
}

process_branch() {
    local branch="$1"
    local source_bcf="$2"
    local base_bcf="$WORK/${branch}.autosomes.bcf"

    [[ -s "$source_bcf" ]] || {
        echo "ERROR: input BCF not found: $source_bcf" >&2
        return 1
    }

    echo "===== Subset autosomes and rename: $branch ====="
    bcftools view -R "$AUTOSOME_BED" --regions-overlap 0 -Ou "$source_bcf" \
        | bcftools annotate --rename-chrs "$CHR_MAP" -Ob -o "$base_bcf"

    bcftools index -f "$base_bcf"
    validate_autosome_records "$base_bcf"

    for dp_label in none ge5 ge10; do
        local dp_source="$base_bcf"
        local masked_bcf="$WORK/${branch}.${dp_label}.bcf"

        case "$dp_label" in
            none)
                ;;
            ge5)
                echo "===== Mask genotypes with DP < 5: $branch ====="
                bcftools +setGT "$base_bcf" -Ob -o "$masked_bcf" -- \
                    -t q -n . -i 'FMT/DP<5'
                bcftools index -f "$masked_bcf"
                dp_source="$masked_bcf"
                ;;
            ge10)
                echo "===== Mask genotypes with DP < 10: $branch ====="
                bcftools +setGT "$base_bcf" -Ob -o "$masked_bcf" -- \
                    -t q -n . -i 'FMT/DP<10'
                bcftools index -f "$masked_bcf"
                dp_source="$masked_bcf"
                ;;
        esac

        for maf_label in none ge005 ge010 ge0125; do
            # Reuse the completed no-DP/no-MAF baseline instead of rerunning it.
            if [[ "$dp_label" == "none" && "$maf_label" == "none" ]]; then
                continue
            fi

            local maf_value=""
            local scenario_bcf="$dp_source"
            local filtered_bcf="$WORK/${branch}.${dp_label}.${maf_label}.bcf"

            case "$maf_label" in
                none) ;;
                ge005) maf_value="0.05" ;;
                ge010) maf_value="0.10" ;;
                ge0125) maf_value="0.125" ;;
            esac

            if [[ -n "$maf_value" ]]; then
                echo "===== Apply MAF >= $maf_value: $branch / DP $dp_label ====="
                bcftools +fill-tags "$dp_source" -Ou -- -t MAF \
                    | bcftools view -i "INFO/MAF >= $maf_value" \
                        -Ob -o "$filtered_bcf"
                bcftools index -f "$filtered_bcf"
                scenario_bcf="$filtered_bcf"
            fi

            validate_autosome_records "$scenario_bcf"

            local sites stats_file called slots missing
            sites=$(bcftools index -n "$scenario_bcf")
            stats_file="$WORK/${branch}.${dp_label}.${maf_label}.stats.txt"
            bcftools stats -s - "$scenario_bcf" > "$stats_file"

            # PSC columns 3-5 count reference homozygotes, alternate homozygotes,
            # and heterozygotes, respectively.
            called=$(awk -F '\t' \
                '$1 == "PSC" { n += $3 + $4 + $5 } END { print n + 0 }' \
                "$stats_file")
            slots=$((sites * 8))
            missing=$((slots - called))

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$branch" "$dp_label" "$maf_label" \
                "$sites" "$called" "$missing" "$slots" >> "$METRICS"

            local prefix="$OUT/roh/$branch/$dp_label/$maf_label/roh"
            mkdir -p "$(dirname "$prefix")"

            echo "===== PLINK ROH: $branch / DP $dp_label / MAF $maf_label ====="
            "$PLINK_BIN" \
                --bcf "$scenario_bcf" \
                --allow-extra-chr \
                --chr-set 27 \
                --double-id \
                --memory 16000 \
                --homozyg \
                --homozyg-window-het 3 \
                --homozyg-window-missing 20 \
                --out "$prefix"

            [[ -s "$prefix.hom.indiv" ]] || {
                echo "ERROR: missing PLINK individual ROH report: $prefix.hom.indiv" >&2
                return 1
            }

            if [[ "$scenario_bcf" == "$filtered_bcf" ]]; then
                rm -f -- "$filtered_bcf" "$filtered_bcf.csi"
            fi
        done

        if [[ "$dp_source" == "$masked_bcf" ]]; then
            rm -f -- "$masked_bcf" "$masked_bcf.csi"
        fi
    done

    rm -f -- "$base_bcf" "$base_bcf.csi"
}

process_branch adapterremoval \
    "$GENO/zambia_pecnerova_track/20261115_013407/adapterremoval/zambia_adapterremoval.raw_biallelic_snps.bcf"
process_branch current_fastp \
    "$GENO/zambia_pecnerova_track/20261115_013407/current_fastp/zambia_current_fastp.raw_biallelic_snps.bcf"
process_branch fastp_merge \
    "$GENO/zambia_pecnerova_track/20261115_013407/fastp_merge/zambia_fastp_merge.raw_biallelic_snps.bcf"

echo "===== Summarize ROH results, including reused baseline ====="
python3 - "$OUT" <<'PY'
import csv
import os
import statistics
import sys
from collections import defaultdict

out_dir = sys.argv[1]
autosome_kb = 2663685.809
roh_root = os.path.join(out_dir, "roh")
rows = []
groups = defaultdict(list)

for root, _, files in os.walk(roh_root):
    if "roh.hom.indiv" not in files:
        continue

    rel = os.path.relpath(root, roh_root).split(os.sep)
    if len(rel) < 3:
        continue

    branch, dp_mask, maf_filter = rel[:3]
    path = os.path.join(root, "roh.hom.indiv")

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0] == "FID" or len(fields) < 6:
                continue

            row = {
                "branch": branch,
                "dp_mask": dp_mask,
                "maf_filter": maf_filter,
                "sample": fields[1],
                "ROH_segments": int(float(fields[3])),
                "total_ROH_kb": float(fields[4]),
            }
            row["FROH_pct"] = row["total_ROH_kb"] / autosome_kb * 100
            rows.append(row)
            groups[(branch, dp_mask, maf_filter)].append(row)

sample_path = os.path.join(out_dir, "roh_by_sample.tsv")
sample_fields = [
    "branch", "dp_mask", "maf_filter", "sample",
    "ROH_segments", "total_ROH_kb", "FROH_pct",
]
with open(sample_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=sample_fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in sorted(rows, key=lambda r: (
        r["branch"], r["dp_mask"], r["maf_filter"], r["sample"]
    )):
        writer.writerow(row)

cohort_path = os.path.join(out_dir, "roh_cohort_summary.tsv")
cohort_fields = [
    "branch", "dp_mask", "maf_filter", "n_samples",
    "individuals_with_ROH", "mean_ROH_segments",
    "mean_total_ROH_kb", "mean_FROH_pct", "median_FROH_pct",
]
with open(cohort_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=cohort_fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for key, values in sorted(groups.items()):
        froh = [row["FROH_pct"] for row in values]
        writer.writerow({
            "branch": key[0],
            "dp_mask": key[1],
            "maf_filter": key[2],
            "n_samples": len(values),
            "individuals_with_ROH": sum(row["ROH_segments"] > 0 for row in values),
            "mean_ROH_segments": f"{statistics.mean(row['ROH_segments'] for row in values):.4f}",
            "mean_total_ROH_kb": f"{statistics.mean(row['total_ROH_kb'] for row in values):.4f}",
            "mean_FROH_pct": f"{statistics.mean(froh):.6f}",
            "median_FROH_pct": f"{statistics.median(froh):.6f}",
        })

print(f"Per-sample ROH table: {sample_path}")
print(f"Cohort summary: {cohort_path}")
PY

echo "Completed DP/MAF sensitivity run: $OUT"
echo "Retention metrics: $METRICS"
echo "Per-sample ROH: $OUT/roh_by_sample.tsv"
echo "Cohort summary: $OUT/roh_cohort_summary.tsv"
