#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv

from ptvlib import clean_gene_id, open_text, parse_float


def find_col(fieldnames: list[str], candidates: list[str]) -> str | None:
    normalized = {f.lower().replace(".", "").replace("_", ""): f for f in fieldnames}
    for cand in candidates:
        key = cand.lower().replace(".", "").replace("_", "")
        if key in normalized:
            return normalized[key]
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open_text(args.input) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(handle, dialect=dialect)
        fields = reader.fieldnames or []
        gene_id_col = find_col(fields, ["gene_id", "ensg", "ensembl_gene_id"])
        symbol_col = find_col(fields, ["gene_symbol", "symbol", "gene", "Gene"])
        shet_col = find_col(fields, ["s_het_drift", "shet_drift", "shet", "s_het", "shet_gene", "selection_coefficient"])
        if not shet_col:
            raise SystemExit("Could not identify s_het column")

        with open(args.out, "w", newline="") as out:
            writer = csv.DictWriter(
                out,
                delimiter="\t",
                fieldnames=["gene_id", "stable_gene_id_without_version", "gene_symbol", "shet", "shet_source", "shet_version"],
                lineterminator="\n",
            )
            writer.writeheader()
            n = 0
            for row in reader:
                shet = parse_float(row.get(shet_col))
                if shet is None:
                    continue
                gene_id = row.get(gene_id_col, "") if gene_id_col else ""
                writer.writerow(
                    {
                        "gene_id": gene_id,
                        "stable_gene_id_without_version": clean_gene_id(gene_id),
                        "gene_symbol": row.get(symbol_col, "") if symbol_col else "",
                        "shet": shet,
                        "shet_source": "user_supplied",
                        "shet_version": "unknown",
                    }
                )
                n += 1
    print(f"wrote {n} normalized s_het rows to {args.out}")


if __name__ == "__main__":
    main()
