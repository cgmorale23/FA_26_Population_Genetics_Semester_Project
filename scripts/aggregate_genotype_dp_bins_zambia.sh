#!/usr/bin/env bash
set -euo pipefail

DP_RUN_DIR="${1:?Usage: aggregate_genotype_dp_bins_zambia.sh DP_RUN_DIR}"
branches=(adapterremoval current_fastp fastp_merge)
tables=(genotype_dp_bins_by_sample.tsv genotype_dp_thresholds_by_sample.tsv)

for table in "${tables[@]}"; do
    output="$DP_RUN_DIR/${table%.tsv}.all_branches.tsv"
    : > "$output"
    first=1
    for branch in "${branches[@]}"; do
        input="$DP_RUN_DIR/$branch/$table"
        if [[ ! -s "$input" ]]; then
            echo "ERROR: missing branch result: $input" >&2
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

cat > "$DP_RUN_DIR/README.md" <<'EOF'
# Zambia genotype DP bins

DP is summarized over variant records present in each branch-level variant-only BCF. The bin table separates called and missing GTs and reports DP=0, 1-4, 5-9, 10-19, 20+, and missing DP. The threshold table counts low-DP called genotypes at DP<5 and DP<10. These are diagnostic summaries; no genotype is filtered or changed.
EOF

echo "Wrote combined DP-bin tables under $DP_RUN_DIR"
