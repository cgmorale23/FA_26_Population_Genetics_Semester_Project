#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="${1:?Usage: aggregate_genotype_call_zambia.sh RUN_DIR}"
ROOT="${ROOT:-/N/scratch/xta1/ZambiaElephants}"
REF="${REF:-$ROOT/rawdata/reference/LoxAfr4/LoxAfr4.fa}"
module purge
module load gnu/14.2.0
module load bcftools/1.17

if [[ ! -s "$REF.fai" ]]; then
    echo "ERROR: missing reference index: $REF.fai" >&2
    exit 2
fi
mapfile -t autosomes < <(awk -F '\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-7])$/ {print $1}' "$REF.fai")
if [[ "${#autosomes[@]}" -ne 27 ]]; then
    echo "ERROR: expected 27 autosomes in $REF.fai, found ${#autosomes[@]}" >&2
    exit 2
fi

expected_samples=()
for sample_number in {1..8}; do
    expected_samples+=("Zambia$(printf '%02d' "$sample_number")")
done

for branch in adapterremoval current_fastp fastp_merge; do
    branch_dir="$RUN_DIR/$branch"
    list="$branch_dir/chromosomes.list"
    output="$branch_dir/zambia_${branch}.raw_biallelic_snps.bcf"
    partial="$output.partial.$$"
    stats="$branch_dir/zambia_${branch}.raw_biallelic_snps.stats.txt"

    : > "$list"
    for chromosome in "${autosomes[@]}"; do
        input="$branch_dir/chromosomes/${chromosome}.snps.bcf"
        if [[ ! -s "$input" || ! -s "$input.csi" ]]; then
            echo "ERROR: missing chromosome BCF or index: $input" >&2
            exit 2
        fi
        printf '%s\n' "$input" >> "$list"
    done

    if [[ -e "$output" || -e "$output.csi" ]]; then
        echo "ERROR: refusing to overwrite existing output: $output" >&2
        exit 3
    fi

    bcftools concat --file-list "$list" --output-type b --output "$partial"
    bcftools index --force "$partial"
    mv -- "$partial" "$output"
    mv -- "$partial.csi" "$output.csi"

    mapfile -t observed_samples < <(bcftools query --list-samples "$output")
    if [[ "${observed_samples[*]}" != "${expected_samples[*]}" ]]; then
        echo "ERROR: unexpected sample list in $output" >&2
        printf 'Expected: %s\nObserved: %s\n' "${expected_samples[*]}" "${observed_samples[*]}" >&2
        exit 4
    fi

    bcftools stats "$output" > "$stats"
    printf 'Branch: %s\nRecords: %s\nSamples: %s\nBCF: %s\nStats: %s\n\n' \
        "$branch" "$(bcftools index --nrecords "$output")" \
        "${observed_samples[*]}" "$output" "$stats"
done

cat > "$RUN_DIR/README.md" <<'EOF'
# Zambia genotype calls by preprocessing branch

Each branch is jointly genotyped across the eight Zambia individuals and kept separate from the other preprocessing branches. Within a branch, paired reads and their complementary collapsed/merged read type are assigned to the same individual using the BAM read-group table; BAMs are streamed directly into `bcftools mpileup` without materializing pooled BAMs.

The calls are restricted to chr1-chr27 and biallelic SNPs. `bcftools mpileup` used minimum mapping quality 30, minimum base quality 25, and `--per-sample-mF`; `bcftools call` used the multiallelic caller. Per-sample DP and AD are annotated. No MAF filter, genotype-depth mask, imputation, or reference-reliability mask is applied in this initial call set. These BCFs are discovery/call inputs for downstream filtering and sensitivity analyses, not final ROH-ready genotypes.

Pečnerová et al. used BCFtools v1.17 with the same core calling options for their high-coverage calls, then applied their reference-reliability filters and removed indels/multiallelic sites. Their imputation was a distinct analysis based on ANGSD genotype likelihoods and BEAGLE across 232 elephants; the Zambia-only eight-sample calls here do not reproduce that imputation panel.
EOF

echo "Wrote run summary: $RUN_DIR/README.md"
