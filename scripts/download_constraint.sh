#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

OUT_DIR="resources/constraint"
URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    -h|--help) echo "Usage: bash scripts/download_constraint.sh --out-dir resources/constraint"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

mkdir -p "${OUT_DIR}"
out="${OUT_DIR}/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
if [[ ! -s "${out}" ]]; then
  curl -fL -C - -o "${out}" "${URL}"
else
  log "[skip] ${out}"
fi
