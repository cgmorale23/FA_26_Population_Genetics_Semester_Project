#!/usr/bin/env bash
#SBATCH --job-name=roh_hardcall
#SBATCH --account=r01340
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/N/scratch/xta1/ZambiaElephants/slurm/logs/roh_hardcall_%j.out
#SBATCH --error=/N/scratch/xta1/ZambiaElephants/slurm/logs/roh_hardcall_%j.err

set -euo pipefail

ROOT="/N/scratch/xta1/ZambiaElephants"
GENO_DIR="$ROOT/qc/genotypes/zambia_pecnerova_track/20261115_013407"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="$ROOT/qc/roh/hardcall_baseline/$RUN_ID"

mkdir -p "$RUN_DIR"

module purge
module load gnu/14.2.0
module load bcftools/1.17

PLINK_BIN="/N/u/xta1/Quartz/tools/plink_1.9_20241022/plink"

if [[ ! -x "$PLINK_BIN" ]]; then
    echo "ERROR: PLINK executable missing: $PLINK_BIN" >&2
    exit 127
fi

"$PLINK_BIN" --version
command -v bcftools

CHR_MAP="$RUN_DIR/chr_rename.tsv"
AUTOSOME_TARGETS=""
: > "$CHR_MAP"

for chr in $(seq 1 27); do
    printf 'chr%s\t%s\n' "$chr" "$chr" >> "$CHR_MAP"
    if [[ -n "$AUTOSOME_TARGETS" ]]; then
        AUTOSOME_TARGETS+=","
    fi
    AUTOSOME_TARGETS+="chr${chr}"
done

run_branch() {
    local branch="$1"
    local bcf="$2"
    local prefix="$RUN_DIR/$branch"
    local renamed="$RUN_DIR/$branch.autosomes.bcf"

    if [[ ! -s "$bcf" ]]; then
        echo "ERROR: missing input BCF: $bcf" >&2
        return 1
    fi

    echo "===== $branch: subset autosomes and rename ====="
    TIMEFORMAT='Autosome subset and rename elapsed: %3R seconds'
    time {
        bcftools view -t "$AUTOSOME_TARGETS" -Ou "$bcf" |
            bcftools annotate --rename-chrs "$CHR_MAP" \
                -Ob -o "$renamed" -
    }

    echo "===== $branch: convert BCF to PLINK binary ====="
    /usr/bin/time -f "PLINK conversion elapsed: %e seconds" \
        "$PLINK_BIN" --bcf "$renamed" \
        --allow-extra-chr \
        --chr-set 27 \
        --double-id \
        --memory 12000 \
        --make-bed \
        --out "$prefix"

    if [[ ! -s "$prefix.bim" ]]; then
        echo "ERROR: PLINK BIM file is missing or empty: $prefix.bim" >&2
        return 1
    fi

    if ! awk '
        $1 !~ /^([1-9]|1[0-9]|2[0-7])$/ {
            print "ERROR: non-autosomal chromosome in BIM: " $1
            bad = 1
        }
        END { exit bad }
    ' "$prefix.bim"; then
        echo "ERROR: stopping before ROH because BIM contains non-autosomal records." >&2
        return 1
    fi

    echo "===== $branch: run baseline ROH ====="
    /usr/bin/time -f "PLINK ROH elapsed: %e seconds" \
        "$PLINK_BIN" --bfile "$prefix" \
        --chr-set 27 \
        --memory 12000 \
        --homozyg \
        --homozyg-window-het 3 \
        --homozyg-window-missing 20 \
        --out "$prefix.roh"

    rm -f "$renamed"
}

run_branch adapterremoval \
    "$GENO_DIR/adapterremoval/zambia_adapterremoval.raw_biallelic_snps.bcf"

run_branch current_fastp \
    "$GENO_DIR/current_fastp/zambia_current_fastp.raw_biallelic_snps.bcf"

run_branch fastp_merge \
    "$GENO_DIR/fastp_merge/zambia_fastp_merge.raw_biallelic_snps.bcf"

echo "Completed baseline ROH run: $RUN_DIR"
