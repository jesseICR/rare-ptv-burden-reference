#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from ptvlib import open_text, parse_float, read_tsv, write_tsv


def percentile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * q
    lower = int(pos)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = pos - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def summarize_variants(path: str) -> dict[str, dict[str, object]]:
    by_tier: dict[str, dict[str, set[str] | int]] = defaultdict(
        lambda: {
            "rows": 0,
            "samples": set(),
            "variants": set(),
            "genes": set(),
            "genes_with_shet": set(),
            "genes_with_cns_percentile": set(),
        }
    )
    for row in read_tsv(path):
        tier = row["tier"]
        bucket = by_tier[tier]
        bucket["rows"] = int(bucket["rows"]) + 1
        bucket["samples"].add(row["sample_id"])  # type: ignore[union-attr]
        bucket["variants"].add(row["variant_id"])  # type: ignore[union-attr]
        gene = row.get("stable_gene_id_without_version") or row.get("gene_id") or row.get("gene_symbol")
        if gene:
            bucket["genes"].add(gene)  # type: ignore[union-attr]
            if row.get("shet") not in ("", None):
                bucket["genes_with_shet"].add(gene)  # type: ignore[union-attr]
            if row.get("cns_expression_percentile") not in ("", None):
                bucket["genes_with_cns_percentile"].add(gene)  # type: ignore[union-attr]

    out: dict[str, dict[str, object]] = {}
    for tier, bucket in by_tier.items():
        out[tier] = {
            "qualifying_rows": int(bucket["rows"]),
            "carrier_samples": len(bucket["samples"]),
            "unique_variants": len(bucket["variants"]),
            "unique_genes": len(bucket["genes"]),
            "unique_genes_with_shet": len(bucket["genes_with_shet"]),
            "unique_genes_with_cns_percentile": len(bucket["genes_with_cns_percentile"]),
        }
    return out


def summarize_burdens(path: str) -> tuple[dict[str, dict[str, object]], int]:
    by_tier: dict[str, list[dict[str, str]]] = defaultdict(list)
    sample_ids = set()
    for row in read_tsv(path):
        by_tier[row["tier"]].append(row)
        sample_ids.add(row["sample_id"])

    out = {}
    for tier, rows in by_tier.items():
        n_variants = [float(row["n_qualifying_variants"]) for row in rows]
        n_genes = [float(row["n_qualifying_genes"]) for row in rows]
        shet = [parse_float(row.get("shet_burden")) or 0.0 for row in rows]
        sum_shet = [parse_float(row.get("sum_shet")) or 0.0 for row in rows]
        carriers = [row for row in rows if int(row["n_qualifying_variants"]) > 0]
        shet_sorted = sorted(shet)
        variant_sorted = sorted(n_variants)
        gene_sorted = sorted(n_genes)
        out[tier] = {
            "n_samples": len(rows),
            "n_carriers": len(carriers),
            "carrier_fraction": len(carriers) / len(rows) if rows else 0.0,
            "mean_variants": sum(n_variants) / len(n_variants) if n_variants else 0.0,
            "median_variants": percentile(variant_sorted, 0.5),
            "p95_variants": percentile(variant_sorted, 0.95),
            "max_variants": max(n_variants) if n_variants else 0.0,
            "mean_genes": sum(n_genes) / len(n_genes) if n_genes else 0.0,
            "median_genes": percentile(gene_sorted, 0.5),
            "p95_genes": percentile(gene_sorted, 0.95),
            "max_genes": max(n_genes) if n_genes else 0.0,
            "mean_shet_burden": sum(shet) / len(shet) if shet else 0.0,
            "median_shet_burden": percentile(shet_sorted, 0.5),
            "p95_shet_burden": percentile(shet_sorted, 0.95),
            "max_shet_burden": max(shet) if shet else 0.0,
            "mean_sum_shet": sum(sum_shet) / len(sum_shet) if sum_shet else 0.0,
        }
    return out, len(sample_ids)


def write_summary_tsv(path: str, variant_summary: dict[str, dict[str, object]], burden_summary: dict[str, dict[str, object]]) -> None:
    tiers = sorted(set(variant_summary) | set(burden_summary))
    rows = []
    for tier in tiers:
        merged = {"tier": tier}
        merged.update(variant_summary.get(tier, {}))
        merged.update(burden_summary.get(tier, {}))
        rows.append(merged)
    fields = [
        "tier",
        "n_samples",
        "n_carriers",
        "carrier_fraction",
        "qualifying_rows",
        "unique_variants",
        "unique_genes",
        "unique_genes_with_shet",
        "unique_genes_with_cns_percentile",
        "mean_variants",
        "median_variants",
        "p95_variants",
        "max_variants",
        "mean_genes",
        "median_genes",
        "p95_genes",
        "max_genes",
        "mean_shet_burden",
        "median_shet_burden",
        "p95_shet_burden",
        "max_shet_burden",
        "mean_sum_shet",
    ]
    write_tsv(path, rows, fields)


def top_carriers(path: str, out: str, limit: int) -> None:
    rows = list(read_tsv(path))
    rows.sort(
        key=lambda row: (
            row["tier"],
            -int(row["n_qualifying_variants"]),
            -(parse_float(row.get("shet_burden")) or 0.0),
            row["sample_id"],
        )
    )
    selected = []
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        tier = row["tier"]
        if int(row["n_qualifying_variants"]) == 0:
            continue
        if counts[tier] >= limit:
            continue
        selected.append(row)
        counts[tier] += 1
    if not rows:
        Path(out).write_text("")
        return
    with open_text(out, "wt") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize 1000G EUR PTV reference distributions by tier.")
    parser.add_argument("--qualifying", required=True, help="qualifying_variants_by_tier.tsv.gz")
    parser.add_argument("--burden", required=True, help="burden_by_sample.tsv.gz")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--top-n", type=int, default=10)
    args = parser.parse_args()

    variant_summary = summarize_variants(args.qualifying)
    burden_summary, n_samples = summarize_burdens(args.burden)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_tsv = f"{args.out_prefix}.tier_summary.tsv"
    top_tsv = f"{args.out_prefix}.top_carriers.tsv"
    summary_json = f"{args.out_prefix}.summary.json"
    write_summary_tsv(summary_tsv, variant_summary, burden_summary)
    top_carriers(args.burden, top_tsv, args.top_n)
    Path(summary_json).write_text(
        json.dumps(
            {
                "n_samples": n_samples,
                "variant_summary": variant_summary,
                "burden_summary": burden_summary,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"wrote {summary_tsv}")
    print(f"wrote {top_tsv}")
    print(f"wrote {summary_json}")


if __name__ == "__main__":
    main()
