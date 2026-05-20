#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

CHROMOSOMES="1-22"
OUT_DIR=""
BASE_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chromosomes) CHROMOSOMES="${2:?}"; shift 2 ;;
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    -h|--help) echo "Usage: bash scripts/download_gnomad.sh --chromosomes 22 --out-dir /path/to/gnomad_v4.1_joint"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "${OUT_DIR}" ]] || die "--out-dir is required; choose a gnomAD v4.1 directory on a filesystem with enough space."

mkdir -p "${OUT_DIR}"
for chr in $(expand_chromosomes "${CHROMOSOMES}"); do
  label="$(chrom_label "${chr}")"
  file="gnomad.joint.v4.1.sites.${label}.vcf.bgz"
  if [[ ! -s "${OUT_DIR}/${file}" ]]; then
    curl -fLsS --retry 20 --retry-delay 10 --retry-all-errors -C - \
      -o "${OUT_DIR}/${file}.part" "${BASE_URL}/${file}"
    mv "${OUT_DIR}/${file}.part" "${OUT_DIR}/${file}"
  else
    log "[skip] ${OUT_DIR}/${file}"
  fi
  if [[ ! -s "${OUT_DIR}/${file}.tbi" ]]; then
    curl -fLsS --retry 20 --retry-delay 5 --retry-all-errors -C - \
      -o "${OUT_DIR}/${file}.tbi.part" "${BASE_URL}/${file}.tbi"
    mv "${OUT_DIR}/${file}.tbi.part" "${OUT_DIR}/${file}.tbi"
  else
    log "[skip] ${OUT_DIR}/${file}.tbi"
  fi
done
