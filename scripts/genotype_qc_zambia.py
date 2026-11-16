#!/usr/bin/env python3
"""Summarize per-sample and per-site QC from branch-level Zambia BCF files."""

import argparse
import csv
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path


EXPECTED_SAMPLES = ["Zambia{:02d}".format(i) for i in range(1, 9)]
MAF_THRESHOLDS = (0.0, 0.05, 0.10, 0.125)


def percentile(histogram, n, fraction):
    if n == 0:
        return "NA"
    target = max(1, math.ceil(n * fraction))
    cumulative = 0
    for value in sorted(histogram):
        cumulative += histogram[value]
        if cumulative >= target:
            return value
    return max(histogram)


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bcf", required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    bcf = Path(args.bcf)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not bcf.is_file() or bcf.stat().st_size == 0:
        raise FileNotFoundError("Missing or empty BCF: {}".format(bcf))
    if not Path(str(bcf) + ".csi").is_file():
        raise FileNotFoundError("Missing CSI index: {}.csi".format(bcf))

    sample_result = subprocess.run(
        ["bcftools", "query", "--list-samples", str(bcf)],
        check=True, text=True, stdout=subprocess.PIPE,
    )
    samples = sample_result.stdout.splitlines()
    if samples != EXPECTED_SAMPLES:
        raise ValueError("Unexpected sample order/list: {}".format(",".join(samples)))

    metrics = {
        sample: {
            "sites": 0,
            "called": 0,
            "missing": 0,
            "hom_ref": 0,
            "heterozygous": 0,
            "hom_alt": 0,
            "dp_all_hist": Counter(),
            "dp_called_hist": Counter(),
            "dp_all_n": 0,
            "dp_called_n": 0,
            "het_ad_missing": 0,
            "het_ad_n": 0,
            "het_alt_sum": 0.0,
            "het_alt_hist": Counter(),
            "het_imbalance": 0,
        }
        for sample in samples
    }
    site_call_counts = Counter()
    maf_counts = Counter()
    n_sites = 0
    query_format = "%CHROM\t%POS[\t%SAMPLE\t%GT\t%DP\t%AD]\n"
    process = subprocess.Popen(
        ["bcftools", "query", "--format", query_format, str(bcf)],
        text=True, stdout=subprocess.PIPE,
    )
    assert process.stdout is not None
    for line_number, line in enumerate(process.stdout, 1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 2 + 4 * len(samples):
            process.kill()
            raise ValueError("Unexpected query field count at line {}: {}".format(line_number, len(fields)))
        n_sites += 1
        site_called = 0
        alt_allele_count = 0
        allele_count = 0
        for offset, sample in enumerate(samples):
            base = 2 + 4 * offset
            observed_sample, gt, dp_text, ad_text = fields[base:base + 4]
            if observed_sample != sample:
                process.kill()
                raise ValueError("Sample order mismatch at site {}: {} != {}".format(
                    line_number, observed_sample, sample,
                ))
            item = metrics[sample]
            item["sites"] += 1
            gt_alleles = gt.replace("|", "/").split("/")
            called_gt = bool(gt_alleles) and all(a != "." for a in gt_alleles)
            if called_gt:
                allele_values = [int(a) for a in gt_alleles]
                item["called"] += 1
                site_called += 1
                allele_count += len(allele_values)
                alt_allele_count += sum(allele_values)
                if all(a == 0 for a in allele_values):
                    item["hom_ref"] += 1
                elif all(a == 1 for a in allele_values):
                    item["hom_alt"] += 1
                else:
                    item["heterozygous"] += 1
            else:
                item["missing"] += 1

            if dp_text not in (".", ""):
                depth = int(dp_text)
                item["dp_all_hist"][depth] += 1
                item["dp_all_n"] += 1
                if called_gt:
                    item["dp_called_hist"][depth] += 1
                    item["dp_called_n"] += 1

            if called_gt and len(gt_alleles) == 2 and gt_alleles[0] != gt_alleles[1]:
                if ad_text in (".", ""):
                    item["het_ad_missing"] += 1
                else:
                    ad = ad_text.split(",")
                    if len(ad) >= 2 and ad[0] != "." and ad[1] != ".":
                        ref_count, alt_count = int(ad[0]), int(ad[1])
                        total = ref_count + alt_count
                        if total > 0:
                            alt_fraction = alt_count / total
                            bin_id = min(1000, max(0, int(round(alt_fraction * 1000))))
                            item["het_ad_n"] += 1
                            item["het_alt_sum"] += alt_fraction
                            item["het_alt_hist"][bin_id] += 1
                            if alt_fraction < 0.20 or alt_fraction > 0.80:
                                item["het_imbalance"] += 1
                        else:
                            item["het_ad_missing"] += 1
                    else:
                        item["het_ad_missing"] += 1

        site_call_counts[site_called] += 1
        if allele_count > 0:
            minor_count = min(alt_allele_count, allele_count - alt_allele_count)
            maf = minor_count / allele_count
            for threshold in MAF_THRESHOLDS:
                if maf >= threshold:
                    maf_counts[threshold] += 1

    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError("bcftools query failed with exit code {}".format(return_code))
    if n_sites == 0:
        raise ValueError("BCF contains no variant records: {}".format(bcf))

    sample_fields = [
        "branch", "sample", "variant_sites_in_bcf", "called_genotypes", "missing_genotypes",
        "call_rate_pct_at_variant_sites", "hom_ref_genotypes", "heterozygous_genotypes",
        "hom_alt_genotypes", "dp_observations_all_genotypes", "mean_dp_all_genotypes",
        "median_dp_all_genotypes", "p05_dp_all_genotypes", "p95_dp_all_genotypes",
        "dp_observations_called_genotypes", "mean_dp_called_genotypes", "median_dp_called_genotypes",
        "p05_dp_called_genotypes", "p95_dp_called_genotypes", "heterozygotes_with_AD",
        "heterozygotes_missing_AD", "mean_het_alt_fraction", "median_het_alt_fraction",
        "het_alt_fraction_outside_0.20_0.80_pct",
    ]
    sample_rows = []
    for sample in samples:
        item = metrics[sample]
        dp_all_sum = sum(depth * count for depth, count in item["dp_all_hist"].items())
        dp_called_sum = sum(depth * count for depth, count in item["dp_called_hist"].items())
        sample_rows.append({
            "branch": args.branch,
            "sample": sample,
            "variant_sites_in_bcf": item["sites"],
            "called_genotypes": item["called"],
            "missing_genotypes": item["missing"],
            "call_rate_pct_at_variant_sites": "{:.6f}".format(100 * item["called"] / item["sites"]),
            "hom_ref_genotypes": item["hom_ref"],
            "heterozygous_genotypes": item["heterozygous"],
            "hom_alt_genotypes": item["hom_alt"],
            "dp_observations_all_genotypes": item["dp_all_n"],
            "mean_dp_all_genotypes": "{:.5f}".format(dp_all_sum / item["dp_all_n"]) if item["dp_all_n"] else "NA",
            "median_dp_all_genotypes": percentile(item["dp_all_hist"], item["dp_all_n"], 0.50),
            "p05_dp_all_genotypes": percentile(item["dp_all_hist"], item["dp_all_n"], 0.05),
            "p95_dp_all_genotypes": percentile(item["dp_all_hist"], item["dp_all_n"], 0.95),
            "dp_observations_called_genotypes": item["dp_called_n"],
            "mean_dp_called_genotypes": "{:.5f}".format(dp_called_sum / item["dp_called_n"]) if item["dp_called_n"] else "NA",
            "median_dp_called_genotypes": percentile(item["dp_called_hist"], item["dp_called_n"], 0.50),
            "p05_dp_called_genotypes": percentile(item["dp_called_hist"], item["dp_called_n"], 0.05),
            "p95_dp_called_genotypes": percentile(item["dp_called_hist"], item["dp_called_n"], 0.95),
            "heterozygotes_with_AD": item["het_ad_n"],
            "heterozygotes_missing_AD": item["het_ad_missing"],
            "mean_het_alt_fraction": "{:.5f}".format(item["het_alt_sum"] / item["het_ad_n"]) if item["het_ad_n"] else "NA",
            "median_het_alt_fraction": "{:.3f}".format(percentile(item["het_alt_hist"], item["het_ad_n"], 0.50) / 1000) if item["het_ad_n"] else "NA",
            "het_alt_fraction_outside_0.20_0.80_pct": "{:.5f}".format(100 * item["het_imbalance"] / item["het_ad_n"]) if item["het_ad_n"] else "NA",
        })
    write_tsv(outdir / "genotype_qc_by_sample.tsv", sample_fields, sample_rows)

    site_rows = []
    for called_samples in range(len(samples) + 1):
        count = site_call_counts[called_samples]
        site_rows.append({
            "branch": args.branch,
            "called_samples_at_site": called_samples,
            "missing_samples_at_site": len(samples) - called_samples,
            "site_count": count,
            "pct_of_variant_sites": "{:.6f}".format(100 * count / n_sites),
        })
    write_tsv(
        outdir / "genotype_qc_site_callability.tsv",
        ["branch", "called_samples_at_site", "missing_samples_at_site", "site_count", "pct_of_variant_sites"],
        site_rows,
    )

    maf_rows = []
    for threshold in MAF_THRESHOLDS:
        retained = maf_counts[threshold]
        maf_rows.append({
            "branch": args.branch,
            "minimum_maf": "{:.3f}".format(threshold),
            "sites_retained": retained,
            "pct_of_all_variant_sites": "{:.6f}".format(100 * retained / n_sites),
            "sites_removed": n_sites - retained,
        })
    write_tsv(
        outdir / "genotype_qc_maf_sensitivity.tsv",
        ["branch", "minimum_maf", "sites_retained", "pct_of_all_variant_sites", "sites_removed"],
        maf_rows,
    )

    print("Branch: {}".format(args.branch))
    print("Variant records scanned: {}".format(n_sites))
    print("Wrote {}".format(outdir / "genotype_qc_by_sample.tsv"))
    print("Wrote {}".format(outdir / "genotype_qc_site_callability.tsv"))
    print("Wrote {}".format(outdir / "genotype_qc_maf_sensitivity.tsv"))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(1)
