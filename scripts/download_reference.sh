#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

OUT_DIR="resources/reference"
GENCODE_RELEASE="${GENCODE_RELEASE:-47}"
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    --gencode-release) GENCODE_RELEASE="${2:?}"; shift 2 ;;
    -h|--help) echo "Usage: bash scripts/download_reference.sh --out-dir resources/reference --gencode-release 47"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

need_cmd gunzip
need_cmd samtools
mkdir -p "${OUT_DIR}"

gtf="gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"
fa_gz="GRCh38.primary_assembly.genome.fa.gz"

if [[ ! -s "${OUT_DIR}/${gtf}" ]]; then
  curl -fL -C - -o "${OUT_DIR}/${gtf}" "${GENCODE_BASE}/release_${GENCODE_RELEASE}/${gtf}"
else
  log "[skip] ${OUT_DIR}/${gtf}"
fi
if [[ ! -s "${OUT_DIR}/${fa_gz}" && ! -s "${OUT_DIR}/GRCh38.primary_assembly.genome.fa" ]]; then
  curl -fL -C - -o "${OUT_DIR}/${fa_gz}" "${GENCODE_BASE}/release_${GENCODE_RELEASE}/${fa_gz}"
else
  log "[skip] ${OUT_DIR}/${fa_gz}"
fi

if [[ ! -s "${OUT_DIR}/GRCh38.primary_assembly.genome.fa" ]]; then
  gunzip -c "${OUT_DIR}/${fa_gz}" > "${OUT_DIR}/GRCh38.primary_assembly.genome.fa"
fi
if [[ ! -s "${OUT_DIR}/GRCh38.primary_assembly.genome.fa.fai" ]]; then
  samtools faidx "${OUT_DIR}/GRCh38.primary_assembly.genome.fa"
else
  log "[skip] ${OUT_DIR}/GRCh38.primary_assembly.genome.fa.fai"
fi
