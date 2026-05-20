#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

OUT_DIR="resources/gtex"
URL="https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    -h|--help) echo "Usage: bash scripts/download_gtex.sh --out-dir resources/gtex"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

mkdir -p "${OUT_DIR}"
gct="${OUT_DIR}/GTEx_Analysis_v8_gene_median_tpm.gct.gz"
out="${OUT_DIR}/gene_cns_expression_percentiles.tsv"
if [[ ! -s "${gct}" ]]; then
  curl -fL -C - -o "${gct}" "${URL}"
else
  log "[skip] ${gct}"
fi
if [[ ! -s "${out}" ]]; then
  python3 "${ROOT_DIR}/scripts/compute_gtex_cns_percentiles.py" \
    --gct "${gct}" \
    --out "${out}"
else
  log "[skip] ${out}"
fi
