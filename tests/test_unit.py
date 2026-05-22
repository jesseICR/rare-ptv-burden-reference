from __future__ import annotations

import csv
import gzip
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(*args: str, cwd: Path | None = None) -> None:
    subprocess.run(args, cwd=cwd or ROOT, check=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_gene_annotation_merge_attaches_symbol_only_shet(tmp_path: Path) -> None:
    constraint = tmp_path / "constraint.tsv"
    constraint.write_text(
        "gene\tgene_id\tpLI\toe_lof_upper\n"
        "GENE1\tENSG000001.2\t0.91\t0.33\n"
    )
    shet = tmp_path / "shet.tsv"
    shet.write_text(
        "gene_id\tstable_gene_id_without_version\tgene_symbol\tshet\tshet_source\tshet_version\n"
        "\t\tGENE1\t0.12\ttest\ttest\n"
    )
    gtex = tmp_path / "gtex.tsv"
    gtex.write_text(
        "gene_id\tstable_gene_id_without_version\tgene_symbol\tcns_expression_mean_tpm\tcns_expression_percentile\tn_cns_tissues\texpression_source\n"
        "ENSG000001.2\tENSG000001\tGENE1\t4.2\t77\t13\ttest\n"
    )
    out = tmp_path / "gene_scores.tsv"
    run(
        "python3",
        "scripts/build_gene_annotations.py",
        "--constraint",
        str(constraint),
        "--shet",
        str(shet),
        "--gtex",
        str(gtex),
        "--out",
        str(out),
    )
    rows = read_tsv(out)
    rows_by_gene_id = {row["stable_gene_id_without_version"]: row for row in rows}
    assert rows_by_gene_id["ENSG000001"]["shet"] == "0.12"
    assert rows_by_gene_id["ENSG000001"]["pli"] == "0.91"
    assert rows_by_gene_id["ENSG000001"]["cns_expression_percentile"] == "77"


def test_parse_ac_tier_and_burden(tmp_path: Path) -> None:
    vcf = tmp_path / "vep.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Gene|SYMBOL|Feature|CANONICAL|MANE_SELECT|BIOTYPE|EXON|INTRON|LoF|LoF_filter|LoF_flags">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
                "chr22\t100\t.\tC\tT\t.\tPASS\tCSQ=T|stop_gained|ENSG000001.2|GENE1|ENST1|YES||protein_coding|2/3||HC||\tGT:GQ:DP:AB\t0/1:99:42:0.5\t0/0:99:40:.",
                "",
            ]
        )
    )
    parsed = tmp_path / "parsed.tsv.gz"
    run("python3", "scripts/parse_vep_csq.py", "--vcf", str(vcf), "--out", str(parsed))
    parsed_rows = read_tsv(parsed)
    assert len(parsed_rows) == 1
    assert parsed_rows[0]["loftee_lof"] == "HC"
    assert parsed_rows[0]["not_last_exon_or_intron"] == "1"

    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\n")
    ac = tmp_path / "ac.tsv.gz"
    carriers = tmp_path / "carriers.tsv.gz"
    run(
        "python3",
        "scripts/compute_1000g_ac.py",
        "--vcf",
        str(vcf),
        "--parsed",
        str(parsed),
        "--samples",
        str(samples),
        "--out-variants",
        str(ac),
        "--out-carriers",
        str(carriers),
    )
    ac_rows = read_tsv(ac)
    assert ac_rows[0]["AC_1000G_EUR_unrelated"] == "1"
    assert len(read_tsv(carriers)) == 1

    annotated = tmp_path / "annotated.tsv.gz"
    rows = ac_rows
    rows[0]["gnomad_af_selected"] = "0.000001"
    rows[0]["cadd_phred"] = "31"
    rows[0]["cadd_version"] = "1.7"
    with gzip.open(annotated, "wt") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    gene_scores = tmp_path / "gene_scores.tsv"
    gene_scores.write_text(
        "gene_id\tstable_gene_id_without_version\tgene_symbol\tshet\tshet_source\tshet_version\tpli\tloeuf\tconstraint_source\tcns_expression_mean_tpm\tcns_expression_percentile\texpression_source\n"
        "ENSG000001.2\tENSG000001\tGENE1\t0.12\ttest\ttest\t0.95\t0.2\ttest\t10\t88\ttest\n"
    )
    qualifying = tmp_path / "qualifying.tsv.gz"
    run(
        "python3",
        "scripts/apply_ptv_tiers.py",
        "--annotations",
        str(annotated),
        "--carriers",
        str(carriers),
        "--gene-scores",
        str(gene_scores),
        "--tiers",
        "config/tiers.tsv",
        "--out",
        str(qualifying),
    )
    qrows = read_tsv(qualifying)
    tiers = {r["tier"] for r in qrows}
    assert "gardner_core_ptv_gnomad1e3_ac1" in tiers
    assert "gardner_core_ptv_gnomad2e6_ac1" in tiers
    assert "rare_pli30_ptv_gnomad1e3_ac1" in tiers

    boundary = tmp_path / "annotated_boundary.tsv.gz"
    rows[0]["gnomad_af_selected"] = "0.000002"
    with gzip.open(boundary, "wt") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    boundary_qualifying = tmp_path / "boundary_qualifying.tsv.gz"
    run(
        "python3",
        "scripts/apply_ptv_tiers.py",
        "--annotations",
        str(boundary),
        "--carriers",
        str(carriers),
        "--gene-scores",
        str(gene_scores),
        "--tiers",
        "config/tiers.tsv",
        "--out",
        str(boundary_qualifying),
    )
    boundary_tiers = {r["tier"] for r in read_tsv(boundary_qualifying)}
    assert "gardner_core_ptv_gnomad2e6_ac1" not in boundary_tiers

    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample_id\tpopulation\tsuperpopulation\tsex\trelatedness_filter\nS1\tGBR\tEUR\t1\ttest\nS2\tGBR\tEUR\t2\ttest\n")
    burden = tmp_path / "burden.tsv.gz"
    run(
        "python3",
        "scripts/compute_burden.py",
        "--variants",
        str(qualifying),
        "--samples",
        str(metadata),
        "--tiers",
        "config/tiers.tsv",
        "--out",
        str(burden),
    )
    brows = read_tsv(burden)
    target = [r for r in brows if r["sample_id"] == "S1" and r["tier"] == "gardner_core_ptv_gnomad1e3_ac1"][0]
    assert target["n_qualifying_genes"] == "1"
    assert abs(float(target["shet_burden"]) - 0.12) < 1e-9


def test_1000g_qc_filters_site_genotype_and_optional_ab(tmp_path: Path) -> None:
    vcf = tmp_path / "qc.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
                '##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
                "chr22\t100\t.\tC\tT\t.\tPASS\t.\tGT:GQ:DP:AB\t0/1:99:30:0.5\t0/0:99:30:.",
                "chr22\t101\t.\tC\tA\t.\tLowQual\t.\tGT:GQ:DP:AB\t0/1:99:30:0.5\t0/0:99:30:.",
                "chr22\t102\t.\tC\tG\t.\tPASS\t.\tGT:GQ:DP:AB\t0/1:10:30:0.5\t0/0:99:30:.",
                "chr22\t103\t.\tC\tA\t.\tPASS\t.\tGT:GQ:DP:AB\t0/1:99:30:.\t0/0:99:30:.",
                "chr22\t104\t.\tC\tG\t.\tPASS\t.\tGT:GQ:DP:AB\t0/1:99:30:0.9\t0/0:99:30:.",
                "",
            ]
        )
    )
    parsed = tmp_path / "parsed.tsv.gz"
    parsed_rows = [
        {"variant_id": "chr22:100:C:T", "is_candidate_ptv": "1"},
        {"variant_id": "chr22:101:C:A", "is_candidate_ptv": "1"},
        {"variant_id": "chr22:102:C:G", "is_candidate_ptv": "1"},
        {"variant_id": "chr22:103:C:A", "is_candidate_ptv": "1"},
        {"variant_id": "chr22:104:C:G", "is_candidate_ptv": "1"},
    ]
    with gzip.open(parsed, "wt") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["variant_id", "is_candidate_ptv"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(parsed_rows)
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\n")
    ac = tmp_path / "ac.tsv.gz"
    carriers = tmp_path / "carriers.tsv.gz"
    run(
        "python3",
        "scripts/compute_1000g_ac.py",
        "--vcf",
        str(vcf),
        "--parsed",
        str(parsed),
        "--samples",
        str(samples),
        "--out-variants",
        str(ac),
        "--out-carriers",
        str(carriers),
    )

    ac_by_variant = {row["variant_id"]: row for row in read_tsv(ac)}
    assert ac_by_variant["chr22:100:C:T"]["AC_1000G_EUR_unrelated"] == "1"
    assert ac_by_variant["chr22:101:C:A"]["kgp_site_qc_pass"] == "0"
    assert ac_by_variant["chr22:101:C:A"]["AC_1000G_EUR_unrelated"] == "0"
    assert ac_by_variant["chr22:102:C:G"]["N_1000G_EUR_unrelated_carrier_qc_fail"] == "1"
    assert ac_by_variant["chr22:104:C:G"]["N_1000G_EUR_unrelated_carrier_qc_fail"] == "1"

    carrier_ids = {row["variant_id"] for row in read_tsv(carriers)}
    assert carrier_ids == {"chr22:100:C:T", "chr22:103:C:A"}
