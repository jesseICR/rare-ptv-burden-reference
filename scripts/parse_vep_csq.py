#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

from ptvlib import clean_gene_id, is_ptv_consequence, open_text, parse_info, split_consequences, variant_key, write_tsv


def csq_fields_from_header(vcf: str) -> list[str]:
    pattern = re.compile(r'ID=CSQ,.*Format: ([^"]+)')
    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = pattern.search(line)
                if not match:
                    raise SystemExit("Could not parse VEP CSQ header format")
                return match.group(1).split("|")
            if line.startswith("#CHROM"):
                break
    raise SystemExit("VCF is missing VEP CSQ header")


def get(entry: dict[str, str], *names: str) -> str:
    for name in names:
        if name in entry and entry[name] != "":
            return entry[name]
    return ""


def parse_rank(value: str) -> tuple[int | None, int | None]:
    if not value or "/" not in value:
        return None, None
    left, right = value.split("/", 1)
    try:
        return int(left), int(right)
    except ValueError:
        return None, None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    fields = csq_fields_from_header(args.vcf)
    out_fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_id",
        "allele",
        "consequence",
        "gene_id",
        "stable_gene_id_without_version",
        "gene_symbol",
        "transcript_id",
        "canonical",
        "mane_select",
        "biotype",
        "exon",
        "intron",
        "loftee_lof",
        "loftee_filter",
        "loftee_flags",
        "is_stop_gained",
        "is_frameshift",
        "is_splice_acceptor",
        "is_splice_donor",
        "is_candidate_ptv",
        "is_last_exon",
        "is_last_intron",
        "not_last_exon_or_intron",
    ]

    rows = []
    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alts, _qual, _filter, info_s = parts[:8]
            info = parse_info(info_s)
            csq = info.get("CSQ", "")
            if not csq:
                continue
            alt_set = set(alts.split(","))
            for record in csq.split(","):
                values = record.split("|")
                if len(values) < len(fields):
                    values += [""] * (len(fields) - len(values))
                entry = dict(zip(fields, values))
                consequence = get(entry, "Consequence")
                if not is_ptv_consequence(consequence):
                    continue
                allele = get(entry, "Allele")
                alt = allele if allele in alt_set else (next(iter(alt_set)) if len(alt_set) == 1 else allele)
                gene_id = get(entry, "Gene")
                exon = get(entry, "EXON", "Exon")
                intron = get(entry, "INTRON", "Intron")
                exon_rank, exon_total = parse_rank(exon)
                intron_rank, intron_total = parse_rank(intron)
                is_last_exon = bool(exon_rank and exon_total and exon_rank == exon_total)
                is_last_intron = bool(intron_rank and intron_total and intron_rank == intron_total)
                consequences = set(split_consequences(consequence))
                rows.append(
                    {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "variant_id": variant_key(chrom, pos, ref, alt),
                        "allele": allele,
                        "consequence": consequence,
                        "gene_id": gene_id,
                        "stable_gene_id_without_version": clean_gene_id(gene_id),
                        "gene_symbol": get(entry, "SYMBOL", "Gene_symbol"),
                        "transcript_id": get(entry, "Feature"),
                        "canonical": get(entry, "CANONICAL"),
                        "mane_select": get(entry, "MANE_SELECT"),
                        "biotype": get(entry, "BIOTYPE"),
                        "exon": exon,
                        "intron": intron,
                        "loftee_lof": get(entry, "LoF"),
                        "loftee_filter": get(entry, "LoF_filter"),
                        "loftee_flags": get(entry, "LoF_flags"),
                        "is_stop_gained": int("stop_gained" in consequences),
                        "is_frameshift": int("frameshift_variant" in consequences),
                        "is_splice_acceptor": int("splice_acceptor_variant" in consequences),
                        "is_splice_donor": int("splice_donor_variant" in consequences),
                        "is_candidate_ptv": 1,
                        "is_last_exon": int(is_last_exon),
                        "is_last_intron": int(is_last_intron),
                        "not_last_exon_or_intron": int(not (is_last_exon or is_last_intron)),
                    }
                )

    write_tsv(args.out, rows, out_fields)
    print(f"wrote {len(rows)} parsed PTV transcript rows to {args.out}")


if __name__ == "__main__":
    main()
