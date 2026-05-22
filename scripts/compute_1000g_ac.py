#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os

from ptvlib import open_text, parse_float, parse_gt, read_tsv, truthy, variant_key, write_tsv


def env_float(name: str, default: float | None) -> float | None:
    value = os.environ.get(name)
    if value is None:
        return default
    parsed = parse_float(value)
    return parsed if parsed is not None else default


def env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    return default if value is None else truthy(value)


def format_threshold(value: float | None) -> str:
    return "" if value is None else f"{value:g}"


def sample_field(sample_value: str, fmt_index: dict[str, int], field: str) -> str:
    values = sample_value.split(":") if sample_value else []
    idx = fmt_index.get(field)
    if idx is None or idx >= len(values):
        return ""
    return values[idx]


def numeric_pass(value: str, minimum: float | None) -> bool:
    if minimum is None:
        return True
    parsed = parse_float(value)
    return parsed is not None and parsed >= minimum


def het_ab_pass(value: str, minimum: float | None, maximum: float | None, require_present: bool) -> bool:
    if minimum is None and maximum is None:
        return True
    if value in {"", "."}:
        return not require_present
    parsed = parse_float(value)
    if parsed is None:
        return False
    if minimum is not None and parsed < minimum:
        return False
    if maximum is not None and parsed > maximum:
        return False
    return True


def genotype_qc_pass(
    sample_value: str,
    fmt_index: dict[str, int],
    alt_count: int,
    min_gq: float | None,
    min_dp: float | None,
    het_ab_min: float | None,
    het_ab_max: float | None,
    require_het_ab: bool,
) -> bool:
    if not numeric_pass(sample_field(sample_value, fmt_index, "GQ"), min_gq):
        return False
    if not numeric_pass(sample_field(sample_value, fmt_index, "DP"), min_dp):
        return False
    if alt_count == 1:
        return het_ab_pass(sample_field(sample_value, fmt_index, "AB"), het_ab_min, het_ab_max, require_het_ab)
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--parsed", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--out-variants", required=True)
    parser.add_argument("--out-carriers", required=True)
    parser.add_argument("--require-site-pass", dest="require_site_pass", action="store_true", default=env_bool("KGP_REQUIRE_SITE_PASS", True))
    parser.add_argument("--allow-non-pass-sites", dest="require_site_pass", action="store_false")
    parser.add_argument("--min-gq", type=float, default=env_float("KGP_MIN_GQ", 20.0))
    parser.add_argument("--min-dp", type=float, default=env_float("KGP_MIN_DP", 10.0))
    parser.add_argument("--het-ab-min", type=float, default=env_float("KGP_HET_AB_MIN", 0.2))
    parser.add_argument("--het-ab-max", type=float, default=env_float("KGP_HET_AB_MAX", 0.8))
    parser.add_argument("--require-het-ab", action="store_true", default=env_bool("KGP_REQUIRE_HET_AB", False))
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
            filt = parts[6]
            alt = alts.split(",")[0]
            key = variant_key(chrom, pos, ref, alt)
            if key not in wanted:
                continue
            fmt = parts[8].split(":")
            fmt_index = {field: i for i, field in enumerate(fmt)}
            try:
                gt_idx = fmt_index["GT"]
            except KeyError:
                continue
            site_qc_pass = (not args.require_site_pass) or filt in {"PASS", "."}
            ac = an = n_het = n_hom_alt = 0
            n_genotype_qc_fail = n_carrier_qc_fail = 0
            for sample_i, sample_id in sample_indices:
                sample_value = parts[9 + sample_i]
                gt = sample_field(sample_value, fmt_index, "GT")
                alt_count, allele_count, called = parse_gt(gt)
                if not called:
                    continue
                if not site_qc_pass:
                    continue
                if not genotype_qc_pass(
                    sample_value,
                    fmt_index,
                    alt_count,
                    args.min_gq,
                    args.min_dp,
                    args.het_ab_min,
                    args.het_ab_max,
                    args.require_het_ab,
                ):
                    n_genotype_qc_fail += 1
                    if alt_count > 0:
                        n_carrier_qc_fail += 1
                    continue
                ac += alt_count
                an += allele_count
                if alt_count == 1:
                    n_het += 1
                    carriers.append({
                        "variant_id": key,
                        "sample_id": sample_id,
                        "genotype": gt,
                        "alt_count": alt_count,
                        "genotype_gq": sample_field(sample_value, fmt_index, "GQ"),
                        "genotype_dp": sample_field(sample_value, fmt_index, "DP"),
                        "genotype_ab": sample_field(sample_value, fmt_index, "AB"),
                        "site_filter": filt,
                    })
                elif alt_count >= 2:
                    n_hom_alt += 1
            ac_by_variant[key] = {
                "AC_1000G_EUR_unrelated": ac,
                "AN_1000G_EUR_unrelated": an,
                "AF_1000G_EUR_unrelated": (ac / an) if an else "",
                "N_HET_1000G_EUR_unrelated": n_het,
                "N_HOM_ALT_1000G_EUR_unrelated": n_hom_alt,
                "kgp_site_filter": filt,
                "kgp_site_qc_pass": "1" if site_qc_pass else "0",
                "N_1000G_EUR_unrelated_genotype_qc_fail": n_genotype_qc_fail,
                "N_1000G_EUR_unrelated_carrier_qc_fail": n_carrier_qc_fail,
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
            "kgp_site_filter": "",
            "kgp_site_qc_pass": "0",
            "N_1000G_EUR_unrelated_genotype_qc_fail": 0,
            "N_1000G_EUR_unrelated_carrier_qc_fail": 0,
        }))
        out_rows.append(merged)

    fields = list(out_rows[0].keys()) if out_rows else list(parsed_rows[0].keys())
    write_tsv(args.out_variants, out_rows, fields)
    write_tsv(args.out_carriers, carriers, ["variant_id", "sample_id", "genotype", "alt_count", "genotype_gq", "genotype_dp", "genotype_ab", "site_filter"])
    print(
        "wrote AC for "
        f"{len(out_rows)} annotation rows and {len(carriers)} carriers "
        f"(site PASS required={int(args.require_site_pass)}, min GQ={format_threshold(args.min_gq)}, "
        f"min DP={format_threshold(args.min_dp)}, het AB={format_threshold(args.het_ab_min)}-{format_threshold(args.het_ab_max)}, "
        f"require AB={int(args.require_het_ab)})"
    )


if __name__ == "__main__":
    main()
