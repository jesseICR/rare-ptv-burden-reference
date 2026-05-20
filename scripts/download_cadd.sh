#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

OUT_DIR="resources/cadd/v1.7"
BASE_URL="https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38"
SNV_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    --snv-only) SNV_ONLY=1; shift ;;
    -h|--help) echo "Usage: bash scripts/download_cadd.sh --out-dir resources/cadd/v1.7 [--snv-only]"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

need_cmd curl
mkdir -p "${OUT_DIR}"

download() {
  local name="$1"
  local url="${BASE_URL}/${name}"
  if [[ -s "${OUT_DIR}/${name}" && ! -e "${OUT_DIR}/${name}.lftp-pget-status" ]]; then
    log "[skip] ${OUT_DIR}/${name}"
    return 0
  fi
  log "downloading CADD ${name}"
  if command -v lftp >/dev/null 2>&1; then
    lftp -e "set net:max-retries 5; set net:timeout 30; pget -n 8 -c -O ${OUT_DIR} ${url}; bye" || \
      curl -fL -C - -o "${OUT_DIR}/${name}" "${url}"
  else
    curl -fL -C - -o "${OUT_DIR}/${name}" "${url}"
  fi
}

download "whole_genome_SNVs.tsv.gz"
download "whole_genome_SNVs.tsv.gz.tbi"

if [[ "${SNV_ONLY}" != "1" ]]; then
  download "gnomad.genomes.r4.0.indel.tsv.gz"
  download "gnomad.genomes.r4.0.indel.tsv.gz.tbi"
fi
