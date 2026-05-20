#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict

from ptvlib import parse_float, read_tsv, split_consequences, truthy, write_tsv


def read_tiers(path: str) -> list[dict[str, str]]:
    return list(csv.DictReader(open(path), delimiter="\t"))


def read_gene_scores(path: str) -> dict[str, dict[str, str]]:
    scores: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        for key in (row.get("stable_gene_id_without_version", ""), row.get("gene_id", ""), row.get("gene_symbol", "")):
            if key:
                scores[key] = row
    return scores


def passes_tier(row: dict[str, str], tier: dict[str, str], gene: dict[str, str]) -> bool:
    consequences = set(split_consequences(row.get("consequence", "")))
    allowed = set(split_consequences(tier.get("ptv_consequences", "")))
    if allowed and not consequences.intersection(allowed):
        return False

    ac = parse_float(row.get("AC_1000G_EUR_unrelated"))
    max_ac = parse_float(tier.get("kgp_eur_ac_lte"))
    if max_ac is not None and (ac is None or ac > max_ac):
        return False

    af = parse_float(row.get("gnomad_af_selected"))
    max_af = parse_float(tier.get("gnomad_af_lt"))
    if max_af is not None and (af is None or af >= max_af):
        return False

    if truthy(tier.get("require_loftee_hc")) and row.get("loftee_lof") != "HC":
        return False

    cadd_threshold = parse_float(tier.get("require_cadd_gt"))
    if cadd_threshold and cadd_threshold > 0:
        cadd = parse_float(row.get("cadd_phred"))
        if cadd is None or cadd <= cadd_threshold:
            return False

    if truthy(tier.get("require_not_last_exon_or_intron")) and not truthy(row.get("not_last_exon_or_intron")):
        return False

    pli_min = parse_float(tier.get("pli_min"))
    if pli_min is not None:
        pli = parse_float(gene.get("pli"))
        if pli is None:
            return False
        if tier.get("tier", "").startswith("rare_pli30"):
            if pli <= pli_min:
                return False
        elif pli < pli_min:
            return False

    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--carriers", required=True)
    parser.add_argument("--gene-scores", required=True)
    parser.add_argument("--tiers", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    tiers = read_tiers(args.tiers)
    gene_scores = read_gene_scores(args.gene_scores)
    carriers_by_variant: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.carriers):
        carriers_by_variant[row["variant_id"]].append(row)

    out_rows = []
    seen = set()
    for row in read_tsv(args.annotations):
        key = row.get("stable_gene_id_without_version") or row.get("gene_id") or row.get("gene_symbol")
        gene = gene_scores.get(key, {})
        for tier in tiers:
            if not passes_tier(row, tier, gene):
                continue
            for carrier in carriers_by_variant.get(row["variant_id"], []):
                dedupe = (tier["tier"], carrier["sample_id"], row["variant_id"], row.get("stable_gene_id_without_version"), row.get("gene_symbol"))
                if dedupe in seen:
                    continue
                seen.add(dedupe)
                out_rows.append(
                    {
                        "tier": tier["tier"],
                        "sample_id": carrier["sample_id"],
                        "chrom": row.get("chrom", ""),
                        "pos": row.get("pos", ""),
                        "ref": row.get("ref", ""),
                        "alt": row.get("alt", ""),
                        "variant_id": row.get("variant_id", ""),
                        "gene_id": row.get("gene_id", ""),
                        "stable_gene_id_without_version": row.get("stable_gene_id_without_version", ""),
                        "gene_symbol": row.get("gene_symbol", ""),
                        "transcript_id": row.get("transcript_id", ""),
                        "consequence": row.get("consequence", ""),
                        "loftee_lof": row.get("loftee_lof", ""),
                        "cadd_phred": row.get("cadd_phred", ""),
                        "gnomad_af_selected": row.get("gnomad_af_selected", ""),
                        "AC_1000G_EUR_unrelated": row.get("AC_1000G_EUR_unrelated", ""),
                        "pli": gene.get("pli", ""),
                        "shet": gene.get("shet", ""),
                        "cns_expression_percentile": gene.get("cns_expression_percentile", ""),
                    }
                )

    fields = [
        "tier",
        "sample_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_id",
        "gene_id",
        "stable_gene_id_without_version",
        "gene_symbol",
        "transcript_id",
        "consequence",
        "loftee_lof",
        "cadd_phred",
        "gnomad_af_selected",
        "AC_1000G_EUR_unrelated",
        "pli",
        "shet",
        "cns_expression_percentile",
    ]
    write_tsv(args.out, out_rows, fields)
    print(f"wrote {len(out_rows)} qualifying sample-variant-tier rows")


if __name__ == "__main__":
    main()
