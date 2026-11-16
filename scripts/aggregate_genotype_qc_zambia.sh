#!/usr/bin/env bash
set -euo pipefail

QC_RUN_DIR="${1:?Usage: aggregate_genotype_qc_zambia.sh QC_RUN_DIR}"
branches=(adapterremoval current_fastp fastp_merge)
tables=(genotype_qc_by_sample.tsv genotype_qc_site_callability.tsv genotype_qc_maf_sensitivity.tsv)

for table in "${tables[@]}"; do
    output="$QC_RUN_DIR/$table"
    : > "$output"
    first=1
    for branch in "${branches[@]}"; do
        input="$QC_RUN_DIR/$branch/$table"
        if [[ ! -s "$input" ]]; then
            echo "ERROR: missing branch QC table: $input" >&2
            exit 2
        fi
        if (( first )); then
            cat "$input" > "$output"
            first=0
        else
            tail -n +2 "$input" >> "$output"
        fi
    done
done

cat > "$QC_RUN_DIR/README.md" <<'EOF'
# Zambia raw genotype QC

These summaries describe only biallelic variant records present in each branch-level variant-only BCF. The per-sample call rate is therefore the fraction of those variant sites with a called GT, not a genome-wide callable-genotype rate. DP summaries are calculated at those variant sites. Heterozygous AD balance is summarized as ALT/(REF+ALT); values outside 0.20-0.80 are a diagnostic flag, not an automatic filter.

The MAF table uses hard-called GT allele counts among called genotypes. With eight diploid samples, the smallest nonzero MAF is at least 1/16 (0.0625), so a 0.05 threshold is expected to retain every polymorphic site in a complete set of called genotypes. Inspect callability and the table before choosing filtering thresholds.

No genotypes or source BCF files are modified by this QC step.
EOF

echo "Wrote combined QC tables under $QC_RUN_DIR"
