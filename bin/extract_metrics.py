#!/usr/bin/env python3
"""
Extract benchmarking metrics from hap.py summary CSV output.

Parses the hap.py summary.csv file and produces a JSON metrics file
with SNP/INDEL breakdown and a weighted F1 score for use with the
stimulus optimization framework.
"""

import argparse
import csv
import json
import sys


def parse_happy_csv(csv_path):
    """Parse hap.py summary CSV and extract PASS-filtered SNP and INDEL metrics."""
    snp = None
    indel = None

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("Filter") != "PASS":
                continue

            variant_type = row.get("Type", "").upper()
            if variant_type not in ("SNP", "INDEL"):
                continue

            truth_total = int(float(row.get("TRUTH.TOTAL", 0)))
            tp = int(float(row.get("TRUTH.TP", 0)))
            fn = int(float(row.get("TRUTH.FN", 0)))
            fp = int(float(row.get("QUERY.FP", 0)))

            recall_str = row.get("METRIC.Recall", "")
            precision_str = row.get("METRIC.Precision", "")
            f1_str = row.get("METRIC.F1_Score", "")

            recall = float(recall_str) if recall_str and recall_str != "." else 0.0
            precision = float(precision_str) if precision_str and precision_str != "." else 0.0
            f1_score = float(f1_str) if f1_str and f1_str != "." else 0.0

            metrics = {
                "truth_total": truth_total,
                "tp": tp,
                "fn": fn,
                "fp": fp,
                "recall": round(recall, 6),
                "precision": round(precision, 6),
                "f1_score": round(f1_score, 6),
            }

            if variant_type == "SNP":
                snp = metrics
            elif variant_type == "INDEL":
                indel = metrics

    return snp, indel


def compute_weighted_f1(snp, indel):
    """Compute weighted F1 score across SNPs and INDELs."""
    snp_f1 = snp["f1_score"] if snp else 0.0
    indel_f1 = indel["f1_score"] if indel else 0.0
    snp_total = snp["truth_total"] if snp else 0
    indel_total = indel["truth_total"] if indel else 0
    total = snp_total + indel_total

    if total == 0:
        return 0.0

    weighted = (snp_f1 * snp_total + indel_f1 * indel_total) / total
    return round(weighted, 6)


def build_metrics(sample, variant_caller, aligner, snp, indel):
    """Build the full metrics JSON structure."""
    snp_f1 = snp["f1_score"] if snp else 0.0
    indel_f1 = indel["f1_score"] if indel else 0.0

    total_errors = 0
    if snp:
        total_errors += snp["fn"] + snp["fp"]
    if indel:
        total_errors += indel["fn"] + indel["fp"]

    default_variant = {
        "truth_total": 0,
        "tp": 0,
        "fn": 0,
        "fp": 0,
        "recall": 0.0,
        "precision": 0.0,
        "f1_score": 0.0,
    }

    return {
        "sample": sample,
        "variant_caller": variant_caller,
        "aligner": aligner,
        "snp": snp if snp else default_variant,
        "indel": indel if indel else default_variant,
        "summary": {
            "weighted_f1": compute_weighted_f1(snp, indel),
            "total_errors": total_errors,
            "snp_f1": round(snp_f1, 6),
            "indel_f1": round(indel_f1, 6),
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract benchmarking metrics from hap.py summary CSV"
    )
    parser.add_argument(
        "--summary-csv", required=True, help="Path to hap.py summary.csv"
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument(
        "--variant-caller", required=True, help="Variant caller name"
    )
    parser.add_argument("--aligner", required=True, help="Aligner name")
    parser.add_argument(
        "--output", default="metrics.json", help="Output JSON file path"
    )
    args = parser.parse_args()

    try:
        snp, indel = parse_happy_csv(args.summary_csv)
    except Exception as e:
        # On parse failure, output a degraded metrics file
        error_metrics = build_metrics(
            args.sample, args.variant_caller, args.aligner, None, None
        )
        error_metrics["error"] = str(e)
        with open(args.output, "w") as f:
            json.dump(error_metrics, f, indent=2)
        print(f"WARNING: Failed to parse CSV: {e}", file=sys.stderr)
        sys.exit(0)

    metrics = build_metrics(
        args.sample, args.variant_caller, args.aligner, snp, indel
    )
    with open(args.output, "w") as f:
        json.dump(metrics, f, indent=2)

    print(
        f"Metrics: weighted_f1={metrics['summary']['weighted_f1']}, "
        f"snp_f1={metrics['summary']['snp_f1']}, "
        f"indel_f1={metrics['summary']['indel_f1']}, "
        f"total_errors={metrics['summary']['total_errors']}"
    )


if __name__ == "__main__":
    main()
