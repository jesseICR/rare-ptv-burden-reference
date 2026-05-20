#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict

from ptvlib import parse_float, read_tsv, write_tsv


def read_samples(path: str) -> list[dict[str, str]]:
    return list(read_tsv(path))


def read_tier_names(path: str) -> list[str]:
    return [row["tier"] for row in csv.DictReader(open(path), delimiter="\t")]


def empirical_percentiles(values: dict[str, float]) -> dict[str, float]:
    if not values:
        return {}
    sorted_values = sorted(values.values())
    n = len(sorted_values)
    out = {}
    for sample, value in values.items():
        le = sum(1 for x in sorted_values if x <= value)
        out[sample] = 100.0 * le / n
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--tiers", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    samples = read_samples(args.samples)
    sample_ids = [s["sample_id"] for s in samples]
    tiers = read_tier_names(args.tiers)
    sample_meta = {s["sample_id"]: s for s in samples}

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.variants):
        grouped[(row["sample_id"], row["tier"])].append(row)

    rows = []
    burden_by_tier: dict[str, dict[str, float]] = defaultdict(dict)
    for tier in tiers:
        for sample_id in sample_ids:
            hits = grouped.get((sample_id, tier), [])
            genes: dict[str, dict[str, str]] = {}
            for hit in hits:
                key = hit.get("stable_gene_id_without_version") or hit.get("gene_id") or hit.get("gene_symbol")
                if key and key not in genes:
                    genes[key] = hit
            shet_values = [parse_float(hit.get("shet")) for hit in genes.values()]
            shet_values = [v for v in shet_values if v is not None]
            prod = 1.0
            for value in shet_values:
                prod *= 1.0 - value
            burden = 1.0 - prod if shet_values else 0.0
            cns_values = [parse_float(hit.get("cns_expression_percentile")) for hit in genes.values()]
            cns_values = [v for v in cns_values if v is not None]
            pli_values = [parse_float(hit.get("pli")) for hit in genes.values()]
            pli_values = [v for v in pli_values if v is not None]
            row = {
                "sample_id": sample_id,
                "population": sample_meta.get(sample_id, {}).get("population", ""),
                "superpopulation": sample_meta.get(sample_id, {}).get("superpopulation", ""),
                "tier": tier,
                "n_qualifying_variants": len({h["variant_id"] for h in hits}),
                "n_qualifying_genes": len(genes),
                "n_genes_with_shet": len(shet_values),
                "n_genes_missing_shet": max(0, len(genes) - len(shet_values)),
                "n_high_pli_genes": sum(1 for v in pli_values if v >= 0.9),
                "sum_shet": sum(shet_values),
                "max_shet": max(shet_values) if shet_values else 0,
                "shet_burden": burden,
                "mean_cns_expression_percentile": (sum(cns_values) / len(cns_values)) if cns_values else "",
                "max_cns_expression_percentile": max(cns_values) if cns_values else "",
            }
            rows.append(row)
            burden_by_tier[tier][sample_id] = burden

    pct_by_tier = {tier: empirical_percentiles(values) for tier, values in burden_by_tier.items()}
    for row in rows:
        row["empirical_percentile_within_1000G_EUR"] = pct_by_tier[row["tier"]].get(row["sample_id"], "")

    fields = [
        "sample_id",
        "population",
        "superpopulation",
        "tier",
        "n_qualifying_variants",
        "n_qualifying_genes",
        "n_genes_with_shet",
        "n_genes_missing_shet",
        "n_high_pli_genes",
        "sum_shet",
        "max_shet",
        "shet_burden",
        "mean_cns_expression_percentile",
        "max_cns_expression_percentile",
        "empirical_percentile_within_1000G_EUR",
    ]
    write_tsv(args.out, rows, fields)
    print(f"wrote {len(rows)} sample-tier burden rows")


if __name__ == "__main__":
    main()
