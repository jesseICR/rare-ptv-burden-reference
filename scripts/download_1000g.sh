#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

CHROMOSOMES="1-22"
OUT_DIR=""
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
RAW_DIR="${BASE_URL}/working/20201028_3202_raw_GT_with_annot"
META_URL="${BASE_URL}/20130606_g1k_3202_samples_ped_population.txt"
S3_PREFIX="s3://1000genomes/1000G_2504_high_coverage"
S3_RAW_DIR="${S3_PREFIX}/working/20201028_3202_raw_GT_with_annot"
S3_HTTP_BASE="https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage"
S3_HTTP_RAW_DIR="${S3_HTTP_BASE}/working/20201028_3202_raw_GT_with_annot"
METADATA_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chromosomes) CHROMOSOMES="${2:?}"; shift 2 ;;
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    --metadata-only) METADATA_ONLY=1; shift ;;
    -h|--help) echo "Usage: bash scripts/download_1000g.sh --chromosomes 22 --out-dir /path/to/1000g/raw"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "${OUT_DIR}" ]] || die "--out-dir is required; choose a 1000G raw cache directory on a filesystem with enough space."

mkdir -p "${OUT_DIR}"
metadata="${OUT_DIR}/20130606_g1k_3202_samples_ped_population.txt"
if [[ ! -s "${metadata}" ]]; then
  curl -fLsS --retry 20 --retry-delay 5 --retry-all-errors -C - \
    -o "${metadata}.part" "${S3_HTTP_BASE}/20130606_g1k_3202_samples_ped_population.txt"
  mv "${metadata}.part" "${metadata}"
else
  log "[skip] ${metadata}"
fi

if [[ "${METADATA_ONLY}" == "1" ]]; then
  exit 0
fi

for chr in $(expand_chromosomes "${CHROMOSOMES}"); do
  src="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
  raw="${OUT_DIR}/chr${chr}.vcf.gz"
  index="${raw}.tbi"
  if [[ ! -s "${raw}" || -e "${raw}.lftp-pget-status" ]]; then
    rm -f "${raw}" "${raw}.lftp-pget-status"
    curl -fLsS --retry 20 --retry-delay 10 --retry-all-errors -C - \
      -o "${raw}.part" "${S3_HTTP_RAW_DIR}/${src}"
    mv "${raw}.part" "${raw}"
  else
    log "[skip] ${raw}"
  fi
  if [[ ! -s "${index}" ]]; then
    curl -fLsS --retry 20 --retry-delay 5 --retry-all-errors -C - \
      -o "${index}.part" "${S3_HTTP_RAW_DIR}/${src}.tbi"
    mv "${index}.part" "${index}"
  else
    log "[skip] ${index}"
  fi
done
