#!/usr/bin/env bash
set -euo pipefail

CHR="${1:?chromosome required}"
source "${ROOT_DIR}/scripts/lib.sh"

need_cmd vep
need_cmd tabix

mkdir -p "${WORK_DIR}/chr${CHR}"
in_vcf="${WORK_DIR}/chr${CHR}/candidates.norm.chr${CHR}.vcf.gz"
out_vcf="${WORK_DIR}/chr${CHR}/vep_loftee.chr${CHR}.vcf.gz"
loftee_dir="${RESOURCES_DIR}/loftee/loftee"
loftee_data_dir="${RESOURCES_DIR}/loftee/loftee_data/GRCh38"
cache_dir="${RESOURCES_DIR}/vep/cache"
fasta="${RESOURCES_DIR}/reference/GRCh38.primary_assembly.genome.fa"
gerp_bigwig="${loftee_data_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
human_ancestor_fa="${loftee_data_dir}/human_ancestor.fa.gz"
loftee_plugin="LoF,loftee_path:${loftee_dir}/,gerp_bigwig:${gerp_bigwig},human_ancestor_fa:${human_ancestor_fa},conservation_file:false"

[[ -s "${in_vcf}" ]] || die "Missing input VCF ${in_vcf}"
[[ -d "${loftee_dir}" ]] || die "Missing LOFTEE checkout ${loftee_dir}; run scripts/setup_vep_loftee.sh"
[[ -s "${gerp_bigwig}" ]] || die "Missing LOFTEE GERP bigWig ${gerp_bigwig}; run scripts/download_loftee_aux.sh"
[[ -s "${human_ancestor_fa}" ]] || die "Missing LOFTEE ancestral FASTA ${human_ancestor_fa}; run scripts/download_loftee_aux.sh"
[[ -s "${human_ancestor_fa}.fai" ]] || die "Missing LOFTEE ancestral FASTA index ${human_ancestor_fa}.fai; run scripts/download_loftee_aux.sh"
[[ -s "${human_ancestor_fa}.gzi" ]] || die "Missing LOFTEE ancestral FASTA gzip index ${human_ancestor_fa}.gzi; run scripts/download_loftee_aux.sh"
[[ -d "${cache_dir}" ]] || die "Missing VEP cache directory ${cache_dir}; run vep_install for GRCh38"
[[ -s "${fasta}" ]] || die "Missing FASTA ${fasta}"
export PERL5LIB="${loftee_dir}:${PERL5LIB:-}"

rm -f "${out_vcf}" "${out_vcf}.tbi" "${out_vcf}_summary.html" "${out_vcf}_warnings.txt"

vep \
  --input_file "${in_vcf}" \
  --output_file "${out_vcf}" \
  --vcf --compress_output bgzip \
  --cache --offline --assembly GRCh38 \
  --dir_cache "${cache_dir}" \
  --fasta "${fasta}" \
  --cache_version "${VEP_RELEASE}" \
  --species homo_sapiens \
  --fork "${THREADS}" \
  --symbol --canonical --mane --protein --biotype --numbers \
  --dir_plugins "${loftee_dir}" \
  --plugin "${loftee_plugin}" \
  --force_overwrite

tabix -f -p vcf "${out_vcf}"

if ! bcftools view -h "${out_vcf}" | grep -q 'ID=CSQ'; then
  die "VEP output is missing CSQ header"
fi
if ! bcftools view -h "${out_vcf}" | grep -q 'LoF'; then
  die "VEP output header does not mention LOFTEE LoF fields"
fi
