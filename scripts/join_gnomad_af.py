#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

from ptvlib import format_float, open_text, parse_float, parse_info, read_tsv, variant_key, with_chr, write_tsv


DEFAULT_FIELDS = ["AF_joint", "AF_joint_nfe", "AF_joint_fin", "fafmax_faf95_max_joint", "AC_joint", "AN_joint"]


def gnomad_path(gnomad_dir: str, chromosome: str) -> str:
    label = with_chr(chromosome)
    path = Path(gnomad_dir) / f"gnomad.joint.v4.1.sites.{label}.vcf.bgz"
    if not path.exists():
        raise SystemExit(f"Missing gnomAD VCF: {path}")
    return str(path)


def header_info_fields(vcf: str) -> set[str]:
    fields = set()
    proc = subprocess.run(["bcftools", "view", "-h", vcf], check=True, text=True, capture_output=True)
    for line in proc.stdout.splitlines():
        if line.startswith("##INFO=<ID="):
            fields.add(line.split("ID=", 1)[1].split(",", 1)[0])
    return fields


def allele_value(value: str, alt_index: int = 0) -> str:
    if not value or value == ".":
        return ""
    parts = value.split(",")
    return parts[alt_index] if alt_index < len(parts) else parts[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True)
    parser.add_argument("--gnomad-dir", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--field-report", required=True)
    args = parser.parse_args()

    rows = list(read_tsv(args.variants))
    if not rows:
        write_tsv(args.out, [], ["variant_id"])
        Path(args.field_report).parent.mkdir(parents=True, exist_ok=True)
        Path(args.field_report).write_text("field\tpresent\tused_for_selected_af\n")
        print("no variants for gnomAD lookup")
        return

    vcf = gnomad_path(args.gnomad_dir, args.chromosome)
    present = header_info_fields(vcf)
    af_fields = [f for f in DEFAULT_FIELDS if f in present and f.startswith(("AF_", "faf"))]

    with open(args.field_report, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["field", "present", "used_for_selected_af"], lineterminator="\n")
        writer.writeheader()
        for field in DEFAULT_FIELDS:
            writer.writerow({"field": field, "present": int(field in present), "used_for_selected_af": int(field in af_fields)})

    unique = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in rows}
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as bed:
        bed_path = bed.name
        for chrom, pos, _ref, _alt in sorted(unique, key=lambda x: (x[0], x[1], x[2], x[3])):
            bed.write(f"{chrom}\t{pos - 1}\t{pos}\n")

    proc = subprocess.run(["bcftools", "view", "-H", "-R", bed_path, vcf], text=True, capture_output=True)
    Path(bed_path).unlink(missing_ok=True)
    if proc.returncode != 0:
        raise SystemExit(proc.stderr)

    found: dict[str, dict[str, str]] = {}
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        chrom, pos, _id, ref, alts, _qual, filt, info_s = parts[:8]
        info = parse_info(info_s)
        for alt_i, alt in enumerate(alts.split(",")):
            key = variant_key(chrom, pos, ref, alt)
            values = {field: allele_value(info.get(field, ""), alt_i) for field in DEFAULT_FIELDS}
            selected_values = [parse_float(values.get(field)) for field in af_fields]
            selected_values = [v for v in selected_values if v is not None]
            values.update(
                {
                    "gnomad_observed": "1",
                    "gnomad_filter": filt,
                    "gnomad_af_selected": format_float(max(selected_values) if selected_values else None),
                    "gnomad_selected_af_fields": ",".join(af_fields),
                }
            )
            found[key] = values

    out_rows = []
    for row in rows:
        key = row["variant_id"]
        merged = dict(row)
        if key in found:
            merged.update(found[key])
        else:
            merged.update({field: "" for field in DEFAULT_FIELDS})
            merged.update(
                {
                    "gnomad_observed": "0",
                    "gnomad_filter": "",
                    "gnomad_af_selected": "0",
                    "gnomad_selected_af_fields": ",".join(af_fields),
                }
            )
        out_rows.append(merged)

    fields = list(rows[0].keys()) + DEFAULT_FIELDS + ["gnomad_observed", "gnomad_filter", "gnomad_af_selected", "gnomad_selected_af_fields"]
    write_tsv(args.out, out_rows, fields)
    print(f"wrote gnomAD annotations for {len(out_rows)} rows")


if __name__ == "__main__":
    main()
