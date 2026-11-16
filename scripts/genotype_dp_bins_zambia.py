#!/usr/bin/env python3
"""Count per-sample genotype depths in bins from a branch-level BCF."""

import argparse
import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path


EXPECTED_SAMPLES = ["Zambia{:02d}".format(i) for i in range(1, 9)]
DP_BINS = ("DP=0", "DP=1-4", "DP=5-9", "DP=10-19", "DP>=20", "DP_missing")


def bin_dp(dp_text):
    if dp_text in (".", ""):
        return "DP_missing"
    dp = int(dp_text)
    if dp == 0:
        return "DP=0"
    if dp < 5:
        return "DP=1-4"
    if dp < 10:
        return "DP=5-9"
    if dp < 20:
        return "DP=10-19"
    return "DP>=20"


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
        raise FileNotFoundError("Missing BCF index: {}.csi".format(bcf))

    sample_result = subprocess.run(
        ["bcftools", "query", "--list-samples", str(bcf)],
        check=True, text=True, stdout=subprocess.PIPE,
    )
    samples = sample_result.stdout.splitlines()
    if samples != EXPECTED_SAMPLES:
        raise ValueError("Unexpected sample list/order: {}".format(",".join(samples)))

    counts = {
        sample: {status: Counter() for status in ("called_GT", "missing_GT")}
        for sample in samples
    }
    total_by_status = {
        sample: Counter() for sample in samples
    }
    threshold_counts = {
        sample: Counter() for sample in samples
    }
    n_sites = 0
    query_format = "%CHROM\t%POS[\t%SAMPLE\t%GT\t%DP]\n"
    process = subprocess.Popen(
        ["bcftools", "query", "--format", query_format, str(bcf)],
        text=True, stdout=subprocess.PIPE,
    )
    assert process.stdout is not None

    for line_number, line in enumerate(process.stdout, 1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 2 + 3 * len(samples):
            process.kill()
            raise ValueError("Unexpected field count at query line {}: {}".format(line_number, len(fields)))
        n_sites += 1
        for offset, sample in enumerate(samples):
            base = 2 + 3 * offset
            observed_sample, gt, dp_text = fields[base:base + 3]
            if observed_sample != sample:
                process.kill()
                raise ValueError("Sample order mismatch at query line {}".format(line_number))

            status = "missing_GT" if "." in gt else "called_GT"
            dp_bin = bin_dp(dp_text)
            counts[sample][status][dp_bin] += 1
            total_by_status[sample][status] += 1

            if status == "called_GT" and dp_bin != "DP_missing":
                dp = int(dp_text)
                threshold_counts[sample]["called_gt_with_dp"] += 1
                if dp < 5:
                    threshold_counts[sample]["called_gt_dp_lt_5"] += 1
                if dp < 10:
                    threshold_counts[sample]["called_gt_dp_lt_10"] += 1

    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError("bcftools query failed with exit code {}".format(return_code))
    if n_sites == 0:
        raise ValueError("No variant records found in {}".format(bcf))

    bin_rows = []
    for sample in samples:
        for status in ("called_GT", "missing_GT"):
            denominator = total_by_status[sample][status]
            for dp_bin in DP_BINS:
                count = counts[sample][status][dp_bin]
                bin_rows.append({
                    "branch": args.branch,
                    "sample": sample,
                    "genotype_status": status,
                    "dp_bin": dp_bin,
                    "genotype_count": count,
                    "pct_within_status": "{:.6f}".format(100 * count / denominator) if denominator else "NA",
                })
    write_tsv(
        outdir / "genotype_dp_bins_by_sample.tsv",
        ["branch", "sample", "genotype_status", "dp_bin", "genotype_count", "pct_within_status"],
        bin_rows,
    )

    threshold_rows = []
    for sample in samples:
        item = threshold_counts[sample]
        called_with_dp = item["called_gt_with_dp"]
        threshold_rows.append({
            "branch": args.branch,
            "sample": sample,
            "variant_sites_in_bcf": n_sites,
            "called_genotypes": total_by_status[sample]["called_GT"],
            "called_genotypes_with_DP": called_with_dp,
            "called_genotypes_DP_lt_5": item["called_gt_dp_lt_5"],
            "pct_called_genotypes_DP_lt_5": "{:.6f}".format(100 * item["called_gt_dp_lt_5"] / called_with_dp) if called_with_dp else "NA",
            "called_genotypes_DP_lt_10": item["called_gt_dp_lt_10"],
            "pct_called_genotypes_DP_lt_10": "{:.6f}".format(100 * item["called_gt_dp_lt_10"] / called_with_dp) if called_with_dp else "NA",
            "called_genotypes_DP_missing": counts[sample]["called_GT"]["DP_missing"],
        })
    write_tsv(
        outdir / "genotype_dp_thresholds_by_sample.tsv",
        [
            "branch", "sample", "variant_sites_in_bcf", "called_genotypes",
            "called_genotypes_with_DP", "called_genotypes_DP_lt_5",
            "pct_called_genotypes_DP_lt_5", "called_genotypes_DP_lt_10",
            "pct_called_genotypes_DP_lt_10", "called_genotypes_DP_missing",
        ],
        threshold_rows,
    )

    print("Branch: {}".format(args.branch))
    print("Variant records scanned: {}".format(n_sites))
    print("Wrote {}".format(outdir / "genotype_dp_bins_by_sample.tsv"))
    print("Wrote {}".format(outdir / "genotype_dp_thresholds_by_sample.tsv"))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(1)
