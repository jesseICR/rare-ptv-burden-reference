#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--resources-dir", required=True)
    parser.add_argument("--chromosomes", required=True)
    parser.add_argument("--gnomad-dir", required=True)
    parser.add_argument("--kgp-dir", required=True)
    parser.add_argument("--tiers", default=str(Path(__file__).resolve().parents[1] / "config" / "tiers.tsv"))
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tiers_out = out_dir / "tiers.tsv"
    if Path(args.tiers).resolve() != tiers_out.resolve():
        shutil.copyfile(args.tiers, tiers_out)

    package_files = [
        out_dir / "tool_versions.tsv",
        tiers_out,
        out_dir / "gene_scores_shet_constraint_expression.tsv",
        out_dir / "qualifying_variants_by_tier.tsv.gz",
        out_dir / "burden_by_sample.tsv.gz",
        out_dir / "reference_distribution.tier_summary.tsv",
        out_dir / "reference_distribution.top_carriers.tsv",
        out_dir / "reference_distribution.summary.json",
        out_dir / "README.md",
    ]
    checksums = []
    for path in package_files:
        if path.exists():
            checksums.append({"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size})

    with (out_dir / "checksums.tsv").open("w") as handle:
        handle.write("path\tsha256\tbytes\n")
        for row in checksums:
            handle.write(f"{row['path']}\t{row['sha256']}\t{row['bytes']}\n")

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "genome_build": "GRCh38",
        "chromosomes": args.chromosomes,
        "gnomad_version": "4.1",
        "gnomad_input_name": Path(args.gnomad_dir).name,
        "kgp_input_name": Path(args.kgp_dir).name,
        "kgp_genotype_qc": {
            "require_site_pass": os.environ.get("KGP_REQUIRE_SITE_PASS", "1"),
            "min_gq": os.environ.get("KGP_MIN_GQ", "20"),
            "min_dp": os.environ.get("KGP_MIN_DP", "10"),
            "het_ab_min": os.environ.get("KGP_HET_AB_MIN", "0.2"),
            "het_ab_max": os.environ.get("KGP_HET_AB_MAX", "0.8"),
            "require_het_ab": os.environ.get("KGP_REQUIRE_HET_AB", "0"),
        },
        "cadd_version": "1.7",
        "ptv_tiers": "tiers.tsv",
        "checksums": checksums,
        "note": "Research reference distribution; not a validated clinical score.",
    }
    (out_dir / "reference_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"wrote reference manifest to {out_dir / 'reference_manifest.json'}")


if __name__ == "__main__":
    main()
