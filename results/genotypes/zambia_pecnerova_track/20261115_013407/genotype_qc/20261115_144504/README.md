# Zambia raw genotype QC

These summaries describe only biallelic variant records present in each branch-level variant-only BCF. The per-sample call rate is therefore the fraction of those variant sites with a called GT, not a genome-wide callable-genotype rate. DP summaries are calculated at those variant sites. Heterozygous AD balance is summarized as ALT/(REF+ALT); values outside 0.20-0.80 are a diagnostic flag, not an automatic filter.

The MAF table uses hard-called GT allele counts among called genotypes. With eight diploid samples, the smallest nonzero MAF is at least 1/16 (0.0625), so a 0.05 threshold is expected to retain every polymorphic site in a complete set of called genotypes. Inspect callability and the table before choosing filtering thresholds.

No genotypes or source BCF files are modified by this QC step.
