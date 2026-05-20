#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from ptvlib import read_tsv


def write_counts(path: Path, header: tuple[str, str], counts: Counter[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for key, value in sorted(counts.items()):
            writer.writerow([key, value])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True)
    parser.add_argument("--burden", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tier_counts: Counter[str] = Counter()
    gene_counts: dict[str, set[str]] = defaultdict(set)
    cadd_missing: Counter[str] = Counter()
    for row in read_tsv(args.variants):
        tier = row.get("tier", "")
        tier_counts[tier] += 1
        gene = row.get("stable_gene_id_without_version") or row.get("gene_symbol", "")
        if gene:
            gene_counts[tier].add(gene)
        if tier.startswith("gardner_core") and not row.get("cadd_phred"):
            cadd_missing[tier] += 1

    write_counts(out_dir / "qualifying_variant_counts_by_tier.tsv", ("tier", "n_rows"), tier_counts)
    with (out_dir / "qualifying_gene_counts_by_tier.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tier", "n_genes"])
        for tier, genes in sorted(gene_counts.items()):
            writer.writerow([tier, len(genes)])
    write_counts(out_dir / "cadd_missingness_by_tier.tsv", ("tier", "n_missing_cadd"), cadd_missing)

    burden_summary: dict[str, list[float]] = defaultdict(list)
    for row in read_tsv(args.burden):
        try:
            burden_summary[row["tier"]].append(float(row["shet_burden"]))
        except (KeyError, ValueError):
            pass
    with (out_dir / "burden_distribution_by_tier.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tier", "n_samples", "min", "mean", "max"])
        for tier, values in sorted(burden_summary.items()):
            mean = sum(values) / len(values) if values else 0
            writer.writerow([tier, len(values), min(values) if values else 0, mean, max(values) if values else 0])
    print(f"wrote QC reports to {out_dir}")


if __name__ == "__main__":
    main()
