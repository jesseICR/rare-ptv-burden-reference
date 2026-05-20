#!/usr/bin/env bash
set -euo pipefail
source "${ROOT_DIR}/scripts/lib.sh"

need_cmd python3
need_cmd bcftools
need_cmd tabix
need_cmd bgzip

for chr in ${CHROM_LIST}; do
  label="$(chrom_label "${chr}")"
  gnomad_vcf="${GNOMAD_DIR}/gnomad.joint.v4.1.sites.${label}.vcf.bgz"
  [[ -s "${gnomad_vcf}" ]] || die "Missing gnomAD VCF: ${gnomad_vcf}"
  [[ -s "${gnomad_vcf}.tbi" || -s "${gnomad_vcf}.csi" ]] || die "Missing gnomAD index for ${gnomad_vcf}"
done

if [[ "${STREAM_1000G:-0}" != "1" ]]; then
  [[ -n "${KGP_RAW_DIR:-}" ]] || die "KGP_RAW_DIR is required; pass --kgp-dir to main.sh"
  for chr in ${CHROM_LIST}; do
    raw_vcf="${KGP_RAW_DIR}/chr${chr}.vcf.gz"
    if [[ ! -s "${raw_vcf}" ]]; then
      die "Missing raw 1000G VCF: ${raw_vcf}. Run scripts/download_1000g.sh or pass --download-missing."
    fi
    [[ -s "${raw_vcf}.tbi" || -s "${raw_vcf}.csi" ]] || die "Missing raw 1000G index for ${raw_vcf}"
  done
else
  [[ -n "${KGP_1000G_RAW_BASE:-}" ]] || die "STREAM_1000G=1 requires KGP_1000G_RAW_BASE"
  if ! bcftools view -h "${KGP_1000G_RAW_BASE}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHROM_LIST%%$'\n'*}.recalibrated_variants.vcf.gz" >/dev/null; then
    die "Could not read remote 1000G VCFs through bcftools; install an htslib build with HTTPS support or disable --stream-1000g."
  fi
fi

[[ -s "${KGP_RAW_DIR:-}/20130606_g1k_3202_samples_ped_population.txt" ]] || \
  die "Missing 1000G sample metadata. Run with --download-missing or provide it under ${KGP_RAW_DIR:-<kgp-dir>/raw}."

[[ -s "${RESOURCES_DIR}/reference/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" ]] || \
  die "Missing GENCODE annotation GTF. Run scripts/download_reference.sh or pass --download-missing."
[[ -s "${RESOURCES_DIR}/reference/GRCh38.primary_assembly.genome.fa" ]] || \
  die "Missing GRCh38 FASTA. Run scripts/download_reference.sh or pass --download-missing."
[[ -s "${RESOURCES_DIR}/reference/GRCh38.primary_assembly.genome.fa.fai" ]] || \
  die "Missing GRCh38 FASTA index. Run samtools faidx via scripts/download_reference.sh."
[[ -s "${RESOURCES_DIR}/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz" ]] || \
  die "Missing pLI/constraint resource. Run scripts/download_constraint.sh or pass --download-missing."
[[ -s "${RESOURCES_DIR}/shet/shet.normalized.tsv" ]] || \
  die "Missing normalized s_het resource at ${RESOURCES_DIR}/shet/shet.normalized.tsv. Provide this file before production scoring."
[[ -s "${RESOURCES_DIR}/loftee/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw" ]] || \
  die "Missing LOFTEE GRCh38 GERP bigWig. Run scripts/download_loftee_aux.sh or pass --download-missing."
[[ -s "${RESOURCES_DIR}/loftee/loftee_data/GRCh38/human_ancestor.fa.gz" ]] || \
  die "Missing LOFTEE GRCh38 ancestral FASTA. Run scripts/download_loftee_aux.sh or pass --download-missing."
[[ -s "${RESOURCES_DIR}/loftee/loftee_data/GRCh38/human_ancestor.fa.gz.fai" ]] || \
  die "Missing LOFTEE GRCh38 ancestral FASTA .fai. Run scripts/download_loftee_aux.sh or pass --download-missing."
[[ -s "${RESOURCES_DIR}/loftee/loftee_data/GRCh38/human_ancestor.fa.gz.gzi" ]] || \
  die "Missing LOFTEE GRCh38 ancestral FASTA .gzi. Run scripts/download_loftee_aux.sh or pass --download-missing."

if [[ ! -x "$(command -v vep || true)" ]]; then
  die "VEP executable not found. Install VEP/LOFTEE before running annotation stages."
fi

if [[ "${CHECK_ONLY:-0}" == "1" ]]; then
  echo "check-only validation passed"
else
  echo "input validation passed"
fi
