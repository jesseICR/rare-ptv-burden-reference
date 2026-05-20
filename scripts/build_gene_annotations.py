#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ptvlib import clean_gene_id, open_text


FIELDS = [
    "gene_id",
    "stable_gene_id_without_version",
    "gene_symbol",
    "shet",
    "shet_source",
    "shet_version",
    "pli",
    "loeuf",
    "constraint_source",
    "cns_expression_mean_tpm",
    "cns_expression_percentile",
    "expression_source",
]

IDENTIFIER_FIELDS = {"gene_id", "stable_gene_id_without_version", "gene_symbol"}


def read_constraint(path: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if not Path(path).exists():
        return out
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        symbol_col = "gene" if "gene" in fields else "gene_symbol" if "gene_symbol" in fields else ""
        gene_id_col = "gene_id" if "gene_id" in fields else ""
        pli_col = "pLI" if "pLI" in fields else "pli" if "pli" in fields else ""
        loeuf_col = "oe_lof_upper" if "oe_lof_upper" in fields else "LOEUF" if "LOEUF" in fields else ""
        for row in reader:
            symbol = row.get(symbol_col, "") if symbol_col else ""
            gene_id = row.get(gene_id_col, "") if gene_id_col else ""
            key = clean_gene_id(gene_id) or symbol
            if not key:
                continue
            out[key] = {
                "gene_id": gene_id,
                "stable_gene_id_without_version": clean_gene_id(gene_id),
                "gene_symbol": symbol,
                "pli": row.get(pli_col, "") if pli_col else "",
                "loeuf": row.get(loeuf_col, "") if loeuf_col else "",
                "constraint_source": "gnomAD_v2.1.1_constraint",
            }
    return out


def read_simple(path: str) -> dict[str, dict[str, str]]:
    if not path or not Path(path).exists():
        return {}
    rows: dict[str, dict[str, str]] = {}
    for row in csv.DictReader(open_text(path), delimiter="\t"):
        key = clean_gene_id(row.get("gene_id", "")) or row.get("stable_gene_id_without_version", "") or row.get("gene_symbol", "")
        if key:
            rows[key] = dict(row)
    return rows


def blank(value: object) -> bool:
    return value is None or str(value).strip() == ""


def symbol_key(symbol: str) -> str:
    return (symbol or "").strip().upper()


def empty_record() -> dict[str, str]:
    return {field: "" for field in FIELDS}


def merge_record(target: dict[str, str], source: dict[str, str]) -> None:
    for field, value in source.items():
        if field not in target or blank(value):
            continue
        if field in IDENTIFIER_FIELDS and not blank(target.get(field)):
            continue
        target[field] = str(value)
    if not target["stable_gene_id_without_version"]:
        target["stable_gene_id_without_version"] = clean_gene_id(target["gene_id"])


def add_rows(records: dict[str, dict[str, str]], rows: dict[str, dict[str, str]], prefer_symbol_match: bool = False) -> None:
    symbol_to_key = {
        symbol_key(record.get("gene_symbol", "")): key
        for key, record in records.items()
        if symbol_key(record.get("gene_symbol", ""))
    }

    for fallback_key, row in rows.items():
        stable_gene_id = clean_gene_id(row.get("stable_gene_id_without_version", "")) or clean_gene_id(row.get("gene_id", ""))
        row_symbol_key = symbol_key(row.get("gene_symbol", ""))

        if prefer_symbol_match and row_symbol_key in symbol_to_key:
            key = symbol_to_key[row_symbol_key]
        elif stable_gene_id:
            key = stable_gene_id
        elif row_symbol_key in symbol_to_key:
            key = symbol_to_key[row_symbol_key]
        else:
            key = fallback_key

        record = records.setdefault(key, empty_record())
        merge_record(record, row)

        updated_symbol_key = symbol_key(record.get("gene_symbol", ""))
        if updated_symbol_key:
            symbol_to_key[updated_symbol_key] = key


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--constraint", required=True)
    parser.add_argument("--shet", required=True)
    parser.add_argument("--gtex", default="")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    constraint = read_constraint(args.constraint)
    shet = read_simple(args.shet)
    gtex = read_simple(args.gtex)

    records: dict[str, dict[str, str]] = {}
    add_rows(records, constraint)
    add_rows(records, gtex)
    add_rows(records, shet, prefer_symbol_match=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for key in sorted(records):
            writer.writerow(records[key])
    print(f"wrote {len(records)} gene annotation rows to {args.out}")


if __name__ == "__main__":
    main()
