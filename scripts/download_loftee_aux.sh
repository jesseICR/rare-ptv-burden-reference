#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

OUT_DIR="resources/loftee/loftee_data/GRCh38"
BASE_URL="https://personal.broadinstitute.org/konradk/loftee_data/GRCh38"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    -h|--help) echo "Usage: bash scripts/download_loftee_aux.sh --out-dir resources/loftee/loftee_data/GRCh38"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

need_cmd curl
mkdir -p "${OUT_DIR}"

download() {
  local name="$1"
  if [[ -s "${OUT_DIR}/${name}" ]]; then
    log "[skip] ${OUT_DIR}/${name}"
    return 0
  fi
  log "downloading ${name}"
  curl -fL -C - -o "${OUT_DIR}/${name}" "${BASE_URL}/${name}"
}

download "gerp_conservation_scores.homo_sapiens.GRCh38.bw"
download "human_ancestor.fa.gz"
download "human_ancestor.fa.gz.fai"
download "human_ancestor.fa.gz.gzi"
