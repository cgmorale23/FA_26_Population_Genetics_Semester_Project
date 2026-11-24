#!/usr/bin/env bash
#SBATCH --job-name=gt_dp_discord
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/fastp_merge_gt_discordance_by_dp/gt_dp_%j.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/fastp_merge_gt_discordance_by_dp/gt_dp_%j.err

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
GENO="$ROOT/qc/genotypes"
EXP="$ROOT/qc/roh/fastp_merge_exploration/20261121_123041"

CURRENT="$GENO/zambia_pecnerova_track/20261115_013407/current_fastp/zambia_current_fastp.raw_biallelic_snps.bcf"
MERGE="$GENO/zambia_pecnerova_track/20261115_013407/fastp_merge/zambia_fastp_merge.raw_biallelic_snps.bcf"
PAIRED_ONLY="$GENO/zambia_preprocessing_control/fastp_merge_paired_only/20261118_093101/fastp_merge_paired_only/zambia_fastp_merge_paired_only.raw_biallelic_snps.bcf"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUT="$EXP/dp_discordance/$RUN_ID"
AUTOSOME_DIR="$OUT/autosomes"

mkdir -p "$AUTOSOME_DIR"

module purge
module load gnu/14.2.0
module load bcftools/1.17
module load python/3.12.4

command -v bcftools
bcftools --version | head -n 1
python3 --version

for f in "$CURRENT" "$MERGE" "$PAIRED_ONLY"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty BCF: $f" >&2
        exit 1
    fi
done

REGIONS="chr1"
for chr in $(seq 2 27); do
    REGIONS="$REGIONS,chr$chr"
done

echo "===== Extract autosomes ====="
bcftools view -r "$REGIONS" -Ob -o "$AUTOSOME_DIR/current_fastp.bcf" "$CURRENT"
bcftools index -f "$AUTOSOME_DIR/current_fastp.bcf"

bcftools view -r "$REGIONS" -Ob -o "$AUTOSOME_DIR/fastp_merge.bcf" "$MERGE"
bcftools index -f "$AUTOSOME_DIR/fastp_merge.bcf"

bcftools view -r "$REGIONS" -Ob -o "$AUTOSOME_DIR/fastp_merge_paired_only.bcf" "$PAIRED_ONLY"
bcftools index -f "$AUTOSOME_DIR/fastp_merge_paired_only.bcf"

echo "===== Create exact shared-site BCF pairs ====="

make_isec() {
    local label="$1"
    local bcf_a="$2"
    local bcf_b="$3"
    local pair_dir="$OUT/pairs/$label"

    mkdir -p "$pair_dir"
    bcftools isec -n=2 -c none -p "$pair_dir/isec" "$bcf_a" "$bcf_b"

    for f in "$pair_dir/isec/0000.bcf" "$pair_dir/isec/0001.bcf"; do
        if [[ ! -s "$f" ]]; then
            echo "ERROR: missing common-site BCF: $f" >&2
            exit 1
        fi
        bcftools index -f "$f"
    done
}

make_isec \
    "current_fastp_vs_fastp_merge" \
    "$AUTOSOME_DIR/current_fastp.bcf" \
    "$AUTOSOME_DIR/fastp_merge.bcf"

make_isec \
    "current_fastp_vs_fastp_merge_paired_only" \
    "$AUTOSOME_DIR/current_fastp.bcf" \
    "$AUTOSOME_DIR/fastp_merge_paired_only.bcf"

make_isec \
    "fastp_merge_paired_only_vs_fastp_merge" \
    "$AUTOSOME_DIR/fastp_merge_paired_only.bcf" \
    "$AUTOSOME_DIR/fastp_merge.bcf"

echo "===== Compare genotypes by DP ====="

python3 - "$OUT" <<'PY'
import csv
import itertools
import os
import subprocess
import sys
from collections import Counter

out_dir = sys.argv[1]

comparisons = [
    ("current_fastp_vs_fastp_merge", "current_fastp", "fastp_merge"),
    (
        "current_fastp_vs_fastp_merge_paired_only",
        "current_fastp",
        "fastp_merge_paired_only",
    ),
    (
        "fastp_merge_paired_only_vs_fastp_merge",
        "fastp_merge_paired_only",
        "fastp_merge",
    ),
]

def get_samples(path):
    result = subprocess.run(
        ["bcftools", "query", "-l", path],
        check=True,
        capture_output=True,
        text=True,
    )
    return [line for line in result.stdout.splitlines() if line]

def query_records(path):
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT:%DP]\\n"
    proc = subprocess.Popen(
        ["bcftools", "query", "-f", fmt, path],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )
    return proc

def parse_gt_dp(value):
    parts = value.split(":", 1)
    gt = parts[0] if parts else "."
    dp_text = parts[1] if len(parts) > 1 else "."

    alleles = gt.replace("|", "/").split("/")
    called = bool(alleles) and all(a != "." for a in alleles)

    if called:
        normalized_gt = "/".join(sorted(alleles, key=int))
    else:
        normalized_gt = None

    try:
        dp = int(dp_text)
    except (TypeError, ValueError):
        dp = None

    return normalized_gt, dp, called

def dp_bin(dp):
    if dp is None:
        return "DP_missing"
    if dp < 5:
        return "DP_lt5"
    if dp < 10:
        return "DP_5_to_9"
    return "DP_ge10"

def write_table(path, header, counter, value_names):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for key in sorted(counter):
            writer.writerow([*key, *counter[key]])

branch_counts = Counter()
joint_counts = Counter()
call_counts = Counter()
site_counts = Counter()

for comparison, branch_a, branch_b in comparisons:
    pair_dir = os.path.join(out_dir, "pairs", comparison, "isec")
    bcf_a = os.path.join(pair_dir, "0000.bcf")
    bcf_b = os.path.join(pair_dir, "0001.bcf")

    samples_a = get_samples(bcf_a)
    samples_b = get_samples(bcf_b)
    if samples_a != samples_b:
        raise RuntimeError(
            f"Sample list/order differs for {comparison}: "
            f"{samples_a} versus {samples_b}"
        )

    proc_a = query_records(bcf_a)
    proc_b = query_records(bcf_b)
    sentinel = object()

    for line_a, line_b in itertools.zip_longest(
        proc_a.stdout, proc_b.stdout, fillvalue=sentinel
    ):
        if line_a is sentinel or line_b is sentinel:
            raise RuntimeError(
                f"Different record counts in shared BCF pair: {comparison}"
            )

        fields_a = line_a.rstrip("\n").split("\t")
        fields_b = line_b.rstrip("\n").split("\t")

        if fields_a[:4] != fields_b[:4]:
            raise RuntimeError(
                f"Site mismatch in {comparison}: "
                f"{fields_a[:4]} versus {fields_b[:4]}"
            )

        site_counts[(comparison,)] += 1

        expected_fields = 4 + len(samples_a)
        if len(fields_a) != expected_fields or len(fields_b) != expected_fields:
            raise RuntimeError(
                f"Unexpected field count in {comparison} at {fields_a[:2]}"
            )

        for idx, sample in enumerate(samples_a, start=4):
            gt_a, dp_a, called_a = parse_gt_dp(fields_a[idx])
            gt_b, dp_b, called_b = parse_gt_dp(fields_b[idx])

            if called_a and called_b:
                is_discordant = gt_a != gt_b
                bin_a = dp_bin(dp_a)
                bin_b = dp_bin(dp_b)

                for sample_label in (sample, "ALL_SAMPLES"):
                    branch_counts[
                        (comparison, sample_label, branch_a, bin_a)
                    ][0] += 1
                    branch_counts[
                        (comparison, sample_label, branch_a, bin_a)
                    ][1] += int(is_discordant)

                    branch_counts[
                        (comparison, sample_label, branch_b, bin_b)
                    ][0] += 1
                    branch_counts[
                        (comparison, sample_label, branch_b, bin_b)
                    ][1] += int(is_discordant)

                    joint_counts[
                        (comparison, sample_label, bin_a, bin_b)
                    ][0] += 1
                    joint_counts[
                        (comparison, sample_label, bin_a, bin_b)
                    ][1] += int(is_discordant)

                status = "both_called_discordant" if is_discordant else "both_called_concordant"
            elif called_a and not called_b:
                status = "A_called_B_missing"
            elif not called_a and called_b:
                status = "A_missing_B_called"
            else:
                status = "both_missing"

            for sample_label in (sample, "ALL_SAMPLES"):
                call_counts[(comparison, sample_label, status)] += 1

    proc_a.stdout.close()
    proc_b.stdout.close()
    rc_a = proc_a.wait()
    rc_b = proc_b.wait()
    if rc_a != 0 or rc_b != 0:
        raise RuntimeError(
            f"bcftools query failed for {comparison}: exit codes {rc_a}, {rc_b}"
        )

branch_path = os.path.join(out_dir, "discordance_rate_by_branch_dp.tsv")
with open(branch_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "comparison", "sample", "branch", "DP_bin",
        "both_called_genotypes", "discordant_GT", "discordance_pct",
    ])
    for key in sorted(branch_counts):
        denominator, numerator = branch_counts[key]
        pct = 100.0 * numerator / denominator if denominator else 0.0
        writer.writerow([*key, denominator, numerator, f"{pct:.6f}"])

joint_path = os.path.join(out_dir, "discordance_rate_by_joint_dp.tsv")
with open(joint_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "comparison", "sample", "branch_A_DP_bin", "branch_B_DP_bin",
        "both_called_genotypes", "discordant_GT", "discordance_pct",
    ])
    for key in sorted(joint_counts):
        denominator, numerator = joint_counts[key]
        pct = 100.0 * numerator / denominator if denominator else 0.0
        writer.writerow([*key, denominator, numerator, f"{pct:.6f}"])

call_path = os.path.join(out_dir, "genotype_callability_counts.tsv")
with open(call_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["comparison", "sample", "call_status", "genotype_count"])
    for key in sorted(call_counts):
        writer.writerow([*key, call_counts[key]])

site_path = os.path.join(out_dir, "shared_site_counts.tsv")
with open(site_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["comparison", "shared_autosomal_variant_sites"])
    for key in sorted(site_counts):
        writer.writerow([key[0], site_counts[key]])

print(f"Results written under: {out_dir}")
print(f"Branch DP rates: {branch_path}")
print(f"Joint DP rates: {joint_path}")
print(f"Callability counts: {call_path}")
print(f"Shared site counts: {site_path}")
PY

echo "Completed DP-stratified genotype discordance analysis: $OUT"
