#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from ptvlib import open_text, parse_gt, read_tsv, variant_key, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--parsed", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--out-variants", required=True)
    parser.add_argument("--out-carriers", required=True)
    args = parser.parse_args()

    parsed_rows = list(read_tsv(args.parsed))
    wanted = {r["variant_id"] for r in parsed_rows}
    if not wanted:
        write_tsv(args.out_variants, [], list(parsed_rows[0].keys()) if parsed_rows else ["variant_id"])
        write_tsv(args.out_carriers, [], ["variant_id", "sample_id", "genotype", "alt_count"])
        print("no parsed PTV rows")
        return

    requested_samples = [line.strip() for line in open(args.samples) if line.strip()]
    sample_set = set(requested_samples)
    ac_by_variant: dict[str, dict[str, int]] = {}
    carriers = []
    sample_indices: list[tuple[int, str]] = []

    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                sample_indices = [(i, s) for i, s in enumerate(samples) if s in sample_set]
                if not sample_indices:
                    raise SystemExit("No requested samples found in VCF")
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, _id, ref, alts = parts[:5]
            alt = alts.split(",")[0]
            key = variant_key(chrom, pos, ref, alt)
            if key not in wanted:
                continue
            fmt = parts[8].split(":")
            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                continue
            ac = an = n_het = n_hom_alt = 0
            for sample_i, sample_id in sample_indices:
                sample_value = parts[9 + sample_i]
                gt = sample_value.split(":")[gt_idx] if sample_value else ""
                alt_count, allele_count, called = parse_gt(gt)
                if not called:
                    continue
                ac += alt_count
                an += allele_count
                if alt_count == 1:
                    n_het += 1
                    carriers.append({"variant_id": key, "sample_id": sample_id, "genotype": gt, "alt_count": alt_count})
                elif alt_count >= 2:
                    n_hom_alt += 1
            ac_by_variant[key] = {
                "AC_1000G_EUR_unrelated": ac,
                "AN_1000G_EUR_unrelated": an,
                "AF_1000G_EUR_unrelated": (ac / an) if an else "",
                "N_HET_1000G_EUR_unrelated": n_het,
                "N_HOM_ALT_1000G_EUR_unrelated": n_hom_alt,
            }

    out_rows = []
    for row in parsed_rows:
        merged = dict(row)
        merged.update(ac_by_variant.get(row["variant_id"], {
            "AC_1000G_EUR_unrelated": 0,
            "AN_1000G_EUR_unrelated": 0,
            "AF_1000G_EUR_unrelated": "",
            "N_HET_1000G_EUR_unrelated": 0,
            "N_HOM_ALT_1000G_EUR_unrelated": 0,
        }))
        out_rows.append(merged)

    fields = list(out_rows[0].keys()) if out_rows else list(parsed_rows[0].keys())
    write_tsv(args.out_variants, out_rows, fields)
    write_tsv(args.out_carriers, carriers, ["variant_id", "sample_id", "genotype", "alt_count"])
    print(f"wrote AC for {len(out_rows)} annotation rows and {len(carriers)} carriers")


if __name__ == "__main__":
    main()
