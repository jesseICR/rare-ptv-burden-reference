#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from statistics import mean

from ptvlib import clean_gene_id, open_text, parse_float


def percentile_rank(values: list[float]) -> list[float]:
    order = sorted((v, i) for i, v in enumerate(values))
    out = [0.0] * len(values)
    n = len(values)
    if n <= 1:
        return [100.0] * n
    for rank, (_v, i) in enumerate(order):
        out[i] = 100.0 * rank / (n - 1)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gct", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    with open_text(args.gct) as handle:
        first = handle.readline()
        if first.startswith("#"):
            handle.readline()
        else:
            handle.seek(0)
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        brain_cols = [c for c in fieldnames if "brain" in c.lower()]
        if not brain_cols:
            raise SystemExit("No GTEx brain tissue columns found")
        for row in reader:
            values = [parse_float(row.get(c)) for c in brain_cols]
            values = [v for v in values if v is not None]
            if not values:
                continue
            gene_id = row.get("Name", "")
            rows.append(
                {
                    "gene_id": gene_id,
                    "stable_gene_id_without_version": clean_gene_id(gene_id),
                    "gene_symbol": row.get("Description", ""),
                    "cns_expression_mean_tpm": mean(values),
                    "n_cns_tissues": len(values),
                }
            )

    ranks = percentile_rank([r["cns_expression_mean_tpm"] for r in rows])
    with open(args.out, "w", newline="") as out:
        fields = [
            "gene_id",
            "stable_gene_id_without_version",
            "gene_symbol",
            "cns_expression_mean_tpm",
            "cns_expression_percentile",
            "n_cns_tissues",
            "expression_source",
        ]
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row, pct in zip(rows, ranks):
            row["cns_expression_percentile"] = f"{pct:.6g}"
            row["expression_source"] = "GTEx_v8_gene_median_TPM_brain_columns"
            writer.writerow(row)
    print(f"wrote {len(rows)} GTEx CNS expression rows to {args.out}")


if __name__ == "__main__":
    main()
