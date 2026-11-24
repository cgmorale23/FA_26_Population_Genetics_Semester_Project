#!/usr/bin/env bash
#SBATCH --job-name=roh_param_grid
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/roh_parameter_grid/grid_%j.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/roh_parameter_grid/grid_%j.err

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
GENO="$ROOT/qc/genotypes/zambia_pecnerova_track/20261115_013407"
REF_FAI="$ROOT/rawdata/reference/LoxAfr4/LoxAfr4.fa.fai"
PLINK_BIN="/N/u/xta1/Quartz/tools/plink_1.9_20241022/plink"

RUN_ID="$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual}"
OUT="$ROOT/qc/roh/plink_parameter_grid/$RUN_ID"
WORK="$OUT/work"
ROH_DIR="$OUT/roh"
CHR_MAP="$WORK/chr_rename.tsv"
AUTOSOME_BED="$WORK/autosomes.bed"
SUMMARY="$OUT/roh_parameter_grid_summary.tsv"

mkdir -p "$WORK" "$ROH_DIR"

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

validate_autosome_records() {
    local bcf="$1"
    local contigs="$WORK/$(basename "$bcf").record_contigs.tsv"

    bcftools index -s "$bcf" > "$contigs"

    if ! awk -F '\t' \
        '$1 !~ /^([1-9]|1[0-9]|2[0-7])$/ { print; bad=1 } END { exit bad }' \
        "$contigs"; then
        echo "ERROR: non-autosomal variant records remain in $bcf" >&2
        cat "$contigs" >&2
        exit 1
    fi

    echo "PASS: variant records restricted to autosomes 1-27: $bcf"
}

run_parameter_grid() {
    local branch="$1"
    local source_bcf="$2"
    local base_bcf="$WORK/${branch}.autosomes.bcf"

    [[ -s "$source_bcf" ]] || {
        echo "ERROR: input BCF not found: $source_bcf" >&2
        return 1
    }

    echo "===== Prepare autosomes: $branch ====="
    bcftools view -R "$AUTOSOME_BED" --regions-overlap 0 -Ou "$source_bcf" \
        | bcftools annotate --rename-chrs "$CHR_MAP" -Ob -o "$base_bcf"
    bcftools index -f "$base_bcf"
    validate_autosome_records "$base_bcf"

    for dp_label in none ge5; do
        local dp_source="$base_bcf"
        local masked_bcf="$WORK/${branch}.${dp_label}.bcf"

        if [[ "$dp_label" == "ge5" ]]; then
            echo "===== Mask genotype DP < 5: $branch ====="
            bcftools +setGT "$base_bcf" -Ob -o "$masked_bcf" -- \
                -t q -n . -i 'FMT/DP<5'
            bcftools index -f "$masked_bcf"
            validate_autosome_records "$masked_bcf"
            dp_source="$masked_bcf"
        fi

        for maf_label in none ge005; do
            local scenario_bcf="$dp_source"
            local filtered_bcf="$WORK/${branch}.${dp_label}.${maf_label}.bcf"

            if [[ "$maf_label" == "ge005" ]]; then
                echo "===== Apply MAF >= 0.05: $branch / DP $dp_label ====="
                bcftools +fill-tags "$dp_source" -Ou -- -t MAF \
                    | bcftools view -i 'INFO/MAF >= 0.05' -Ob -o "$filtered_bcf"
                bcftools index -f "$filtered_bcf"
                validate_autosome_records "$filtered_bcf"
                scenario_bcf="$filtered_bcf"
            fi

            for min_kb in 500 1000; do
                for min_snp in 50 100; do
                    for window_het in 1 3; do
                        for window_missing in 5 20; do
                            local parameter_label
                            local prefix

                            parameter_label="kb${min_kb}_snp${min_snp}_het${window_het}_miss${window_missing}"
                            prefix="$ROH_DIR/$branch/$dp_label/$maf_label/$parameter_label/roh"
                            mkdir -p "$(dirname "$prefix")"

                            echo "===== ROH: $branch / DP $dp_label / MAF $maf_label / $parameter_label ====="

                            "$PLINK_BIN" \
                                --bcf "$scenario_bcf" \
                                --allow-extra-chr \
                                --chr-set 27 \
                                --double-id \
                                --memory 16000 \
                                --homozyg \
                                --homozyg-kb "$min_kb" \
                                --homozyg-snp "$min_snp" \
                                --homozyg-density 50 \
                                --homozyg-gap 1000 \
                                --homozyg-window-snp 50 \
                                --homozyg-window-threshold 0.05 \
                                --homozyg-window-het "$window_het" \
                                --homozyg-window-missing "$window_missing" \
                                --out "$prefix"

                            [[ -s "$prefix.hom.indiv" ]] || {
                                echo "ERROR: missing PLINK individual report: $prefix.hom.indiv" >&2
                                return 1
                            }
                        done
                    done
                done
            done

            if [[ "$maf_label" == "ge005" ]]; then
                rm -f -- "$filtered_bcf" "$filtered_bcf.csi"
            fi
        done

        if [[ "$dp_label" == "ge5" ]]; then
            rm -f -- "$masked_bcf" "$masked_bcf.csi"
        fi
    done

    rm -f -- "$base_bcf" "$base_bcf.csi"
}

run_parameter_grid adapterremoval \
    "$GENO/adapterremoval/zambia_adapterremoval.raw_biallelic_snps.bcf"
run_parameter_grid current_fastp \
    "$GENO/current_fastp/zambia_current_fastp.raw_biallelic_snps.bcf"
run_parameter_grid fastp_merge \
    "$GENO/fastp_merge/zambia_fastp_merge.raw_biallelic_snps.bcf"

echo "===== Summarize per-sample and cohort ROH ====="

python3 - "$ROH_DIR" "$SUMMARY" "$OUT/roh_parameter_grid_by_sample.tsv" <<'PY'
import csv
import os
import statistics
import sys
from collections import defaultdict

roh_root, cohort_path, sample_path = sys.argv[1:4]
autosome_kb = 2663685.809
groups = defaultdict(list)
sample_rows = []

for root, _, files in os.walk(roh_root):
    if "roh.hom.indiv" not in files:
        continue

    rel = os.path.relpath(root, roh_root).split(os.sep)
    if len(rel) != 4:
        continue

    branch, dp_mask, maf_filter, parameter_label = rel
    parts = parameter_label.split("_")
    min_kb = parts[0].removeprefix("kb")
    min_snp = parts[1].removeprefix("snp")
    window_het = parts[2].removeprefix("het")
    window_missing = parts[3].removeprefix("miss")
    key = (branch, dp_mask, maf_filter, min_kb, min_snp, window_het, window_missing)

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
                "min_ROH_kb": min_kb,
                "min_SNPs": min_snp,
                "window_het": window_het,
                "window_missing": window_missing,
                "sample": fields[1],
                "ROH_segments": int(float(fields[3])),
                "total_ROH_kb": float(fields[4]),
            }
            row["FROH_pct"] = row["total_ROH_kb"] / autosome_kb * 100
            sample_rows.append(row)
            groups[key].append(row)

sample_fields = [
    "branch", "dp_mask", "maf_filter", "min_ROH_kb", "min_SNPs",
    "window_het", "window_missing", "sample", "ROH_segments",
    "total_ROH_kb", "FROH_pct",
]
with open(sample_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=sample_fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(sorted(
        sample_rows,
        key=lambda r: (
            r["branch"], r["dp_mask"], r["maf_filter"],
            int(r["min_ROH_kb"]), int(r["min_SNPs"]),
            int(r["window_het"]), int(r["window_missing"]), r["sample"]
        )
    ))

cohort_fields = [
    "branch", "dp_mask", "maf_filter", "min_ROH_kb", "min_SNPs",
    "window_het", "window_missing", "n_samples", "individuals_with_ROH",
    "mean_ROH_segments", "mean_total_ROH_kb", "mean_FROH_pct", "median_FROH_pct",
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
            "min_ROH_kb": key[3],
            "min_SNPs": key[4],
            "window_het": key[5],
            "window_missing": key[6],
            "n_samples": len(values),
            "individuals_with_ROH": sum(row["ROH_segments"] > 0 for row in values),
            "mean_ROH_segments": f"{statistics.mean(row['ROH_segments'] for row in values):.4f}",
            "mean_total_ROH_kb": f"{statistics.mean(row['total_ROH_kb'] for row in values):.4f}",
            "mean_FROH_pct": f"{statistics.mean(froh):.6f}",
            "median_FROH_pct": f"{statistics.median(froh):.6f}",
        })

print(f"Cohort summary: {cohort_path}")
print(f"Per-sample summary: {sample_path}")
PY

echo "Completed PLINK parameter grid: $OUT"
echo "Cohort summary: $SUMMARY"
echo "Per-sample summary: $OUT/roh_parameter_grid_by_sample.tsv"
