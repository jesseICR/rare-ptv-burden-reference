#!/usr/bin/env bash
set -euo pipefail

CHR="${1:?chromosome required}"
source "${ROOT_DIR}/scripts/lib.sh"

need_cmd bcftools
need_cmd tabix

mkdir -p "${WORK_DIR}/chr${CHR}"

[[ -n "${KGP_RAW_DIR:-}" ]] || die "KGP_RAW_DIR is required; pass --kgp-dir to main.sh"
raw_vcf="${KGP_RAW_DIR}/chr${CHR}.vcf.gz"
bed="${RESOURCES_DIR}/reference/protein_coding_cds_splice.gencode_v${GENCODE_RELEASE}.bed.gz"
chr_bed="${WORK_DIR}/chr${CHR}/protein_coding_cds_splice.chr${CHR}.bed"
fetch_bed="${WORK_DIR}/chr${CHR}/protein_coding_cds_splice.chr${CHR}.fetch_gap${STREAM_1000G_MERGE_GAP:-10000}.bed"
samples="${RESULTS_DIR}/sample_sets/eur_unrelated.samples.txt"
fasta="${RESOURCES_DIR}/reference/GRCh38.primary_assembly.genome.fa"
out="${WORK_DIR}/chr${CHR}/candidates.chr${CHR}.vcf.gz"
norm="${WORK_DIR}/chr${CHR}/candidates.norm.chr${CHR}.vcf.gz"

if [[ -s "${raw_vcf}" ]]; then
  input_vcf="${raw_vcf}"
elif [[ "${STREAM_1000G:-0}" == "1" ]]; then
  input_vcf="${KGP_1000G_RAW_BASE}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants.vcf.gz"
else
  die "Missing raw VCF ${raw_vcf}"
fi
[[ -s "${bed}" ]] || die "Missing coding/splice BED ${bed}"
[[ -s "${samples}" ]] || die "Missing sample list ${samples}"
[[ -s "${fasta}" ]] || die "Missing FASTA ${fasta}"

if [[ ! -s "${chr_bed}" ]]; then
  zcat "${bed}" | awk -v chr="chr${CHR}" '$1 == chr' > "${chr_bed}"
fi
[[ -s "${chr_bed}" ]] || die "No coding/splice intervals found for chr${CHR} in ${bed}"

if [[ "${STREAM_1000G:-0}" == "1" && "${input_vcf}" =~ ^https?:// ]]; then
  if [[ ! -s "${fetch_bed}" ]]; then
    awk -v gap="${STREAM_1000G_MERGE_GAP:-10000}" '
      NR == 1 { chrom=$1; start=$2; end=$3; next }
      $1 == chrom && $2 <= end + gap { if ($3 > end) end=$3; next }
      { print chrom "\t" start "\t" end; chrom=$1; start=$2; end=$3 }
      END { if (chrom != "") print chrom "\t" start "\t" end }
    ' "${chr_bed}" > "${fetch_bed}"
  fi
  bcftools view \
    --regions-file "${fetch_bed}" \
    --samples-file "${samples}" \
    --output-type u \
    "${input_vcf}" | \
    bcftools view \
      --targets-file "${chr_bed}" \
      --output-type z \
      --output "${out}"
else
  bcftools view \
    --regions-file "${chr_bed}" \
    --samples-file "${samples}" \
    --output-type z \
    --output "${out}" \
    "${input_vcf}"
fi
tabix -f -p vcf "${out}"

bcftools norm -m-any -f "${fasta}" \
  --output-type z \
  --output "${norm}" \
  "${out}"
tabix -f -p vcf "${norm}"
