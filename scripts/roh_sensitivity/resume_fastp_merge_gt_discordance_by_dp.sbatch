#!/usr/bin/env bash
#SBATCH --job-name=gt_dp_resume
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/fastp_merge_gt_discordance_by_dp/resume_%j.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/fastp_merge_gt_discordance_by_dp/resume_%j.err

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
EXP="$ROOT/qc/roh/fastp_merge_exploration/20261121_123041"
OUT="$EXP/dp_discordance/20261121_222411"
AUTOSOME_DIR="$OUT/autosomes"
PAIR_DIR="$OUT/pairs"
RESULTS="$OUT/results"

module purge
module load gnu/14.2.0
module load bcftools/1.17
module load python/3.12.4

mkdir -p "$RESULTS"

for f in \
    "$AUTOSOME_DIR/current_fastp.bcf" \
    "$AUTOSOME_DIR/fastp_merge.bcf" \
    "$AUTOSOME_DIR/fastp_merge_paired_only.bcf"
do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing autosomal BCF: $f" >&2
        exit 1
    fi
done

echo "===== Convert the first pair's existing VCFs to BCF ====="

FIRST="$PAIR_DIR/current_fastp_vs_fastp_merge/isec"
for side in 0000 0001; do
    VCF="$FIRST/$side.vcf"
    BCF="$FIRST/$side.bcf"

    if [[ ! -s "$BCF" ]]; then
        if [[ ! -s "$VCF" ]]; then
            echo "ERROR: missing both VCF and BCF for $side under $FIRST" >&2
            exit 1
        fi
        bcftools view -Ob -o "$BCF.partial" "$VCF"
        mv -- "$BCF.partial" "$BCF"
    fi
    bcftools index -f "$BCF"
done

make_isec() {
    local label="$1"
    local bcf_a="$2"
    local bcf_b="$3"
    local isec_dir="$PAIR_DIR/$label/isec"

    if [[ -s "$isec_dir/0000.bcf" && -s "$isec_dir/0001.bcf" ]]; then
        echo "Using existing BCF intersection: $isec_dir"
    else
        mkdir -p "$PAIR_DIR/$label"
        rm -rf -- "$isec_dir"
        bcftools isec \
            -n=2 \
            -c none \
            -Ob \
            -p "$isec_dir" \
            "$bcf_a" "$bcf_b"
    fi

    for side in 0000 0001; do
        if [[ ! -s "$isec_dir/$side.bcf" ]]; then
            echo "ERROR: expected intersection BCF missing: $isec_dir/$side.bcf" >&2
            find "$isec_dir" -maxdepth 1 -type f -printf '%f\n' >&2 || true
            exit 1
        fi
        bcftools index -f "$isec_dir/$side.bcf"
    done
}

echo "===== Generate remaining exact shared-site BCF pairs ====="

make_isec \
    "current_fastp_vs_fastp_merge_paired_only" \
    "$AUTOSOME_DIR/current_fastp.bcf" \
    "$AUTOSOME_DIR/fastp_merge_paired_only.bcf"

make_isec \
    "fastp_merge_paired_only_vs_fastp_merge" \
    "$AUTOSOME_DIR/fastp_merge_paired_only.bcf" \
    "$AUTOSOME_DIR/fastp_merge.bcf"

echo "===== Calculate DP-stratified discordance rates ====="

python3 - "$PAIR_DIR" "$RESULTS" <<'PY'
import csv
import itertools
import os
import subprocess
import sys
from collections import Counter, defaultdict

pair_root, results_dir = sys.argv[1:3]
os.makedirs(results_dir, exist_ok=True)

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
    return [x for x in result.stdout.splitlines() if x]

def query_records(path):
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT:%DP]\\n"
    return subprocess.Popen(
        ["bcftools", "query", "-f", fmt, path],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )

def parse_gt_dp(value):
    parts = value.split(":", 1)
    gt = parts[0] if parts else "."
    dp_text = parts[1] if len(parts) > 1 else "."

    alleles = gt.replace("|", "/").split("/")
    called = bool(alleles) and all(allele != "." for allele in alleles)
    normalized_gt = "/".join(sorted(alleles, key=int)) if called else None

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

# Each DP stratum stores [both-called denominator, discordant numerator].
branch_counts = defaultdict(lambda: [0, 0])
joint_counts = defaultdict(lambda: [0, 0])
call_counts = Counter()
site_counts = Counter()

for comparison, branch_a, branch_b in comparisons:
    isec_dir = os.path.join(pair_root, comparison, "isec")
    bcf_a = os.path.join(isec_dir, "0000.bcf")
    bcf_b = os.path.join(isec_dir, "0001.bcf")

    for path in (bcf_a, bcf_b):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise FileNotFoundError(path)

    samples_a = get_samples(bcf_a)
    samples_b = get_samples(bcf_b)
    if samples_a != samples_b:
        raise RuntimeError(f"Sample list/order differs in {comparison}")

    proc_a = query_records(bcf_a)
    proc_b = query_records(bcf_b)

    for line_a, line_b in itertools.zip_longest(
        proc_a.stdout, proc_b.stdout, fillvalue=None
    ):
        if line_a is None or line_b is None:
            raise RuntimeError(f"Record counts differ in {comparison}")

        fields_a = line_a.rstrip("\n").split("\t")
        fields_b = line_b.rstrip("\n").split("\t")

        if fields_a[:4] != fields_b[:4]:
            raise RuntimeError(
                f"Variant order/key mismatch in {comparison}: "
                f"{fields_a[:4]} vs {fields_b[:4]}"
            )

        if len(fields_a) != 4 + len(samples_a) or len(fields_b) != 4 + len(samples_b):
            raise RuntimeError(
                f"Unexpected query fields in {comparison} at {fields_a[:2]}"
            )

        site_counts[comparison] += 1

        for column, sample in enumerate(samples_a, start=4):
            gt_a, dp_a, called_a = parse_gt_dp(fields_a[column])
            gt_b, dp_b, called_b = parse_gt_dp(fields_b[column])

            if called_a and called_b:
                discordant = gt_a != gt_b
                bin_a = dp_bin(dp_a)
                bin_b = dp_bin(dp_b)

                for sample_label in (sample, "ALL_SAMPLES"):
                    key_a = (comparison, sample_label, branch_a, bin_a)
                    key_b = (comparison, sample_label, branch_b, bin_b)
                    key_joint = (comparison, sample_label, bin_a, bin_b)

                    branch_counts[key_a][0] += 1
                    branch_counts[key_a][1] += int(discordant)
                    branch_counts[key_b][0] += 1
                    branch_counts[key_b][1] += int(discordant)

                    joint_counts[key_joint][0] += 1
                    joint_counts[key_joint][1] += int(discordant)

                status = (
                    "both_called_discordant"
                    if discordant
                    else "both_called_concordant"
                )
            elif called_a:
                status = "A_called_B_missing"
            elif called_b:
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
            f"bcftools query failed for {comparison}: {rc_a}, {rc_b}"
        )

with open(
    os.path.join(results_dir, "discordance_rate_by_branch_dp.tsv"),
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "comparison", "sample", "branch", "DP_bin",
        "both_called_genotypes", "discordant_GT", "discordance_pct",
    ])
    for key in sorted(branch_counts):
        denominator, numerator = branch_counts[key]
        pct = 100.0 * numerator / denominator if denominator else 0.0
        writer.writerow([*key, denominator, numerator, f"{pct:.6f}"])

with open(
    os.path.join(results_dir, "discordance_rate_by_joint_dp.tsv"),
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "comparison", "sample", "branch_A_DP_bin", "branch_B_DP_bin",
        "both_called_genotypes", "discordant_GT", "discordance_pct",
    ])
    for key in sorted(joint_counts):
        denominator, numerator = joint_counts[key]
        pct = 100.0 * numerator / denominator if denominator else 0.0
        writer.writerow([*key, denominator, numerator, f"{pct:.6f}"])

with open(
    os.path.join(results_dir, "genotype_callability_counts.tsv"),
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["comparison", "sample", "call_status", "genotype_count"])
    for key in sorted(call_counts):
        writer.writerow([*key, call_counts[key]])

with open(
    os.path.join(results_dir, "shared_site_counts.tsv"),
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["comparison", "shared_autosomal_variant_sites"])
    for comparison in sorted(site_counts):
        writer.writerow([comparison, site_counts[comparison]])

print(f"Results written under: {results_dir}")
PY

echo "Completed DP-stratified genotype discordance analysis: $OUT"
