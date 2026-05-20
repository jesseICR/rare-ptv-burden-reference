#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

OUT_DIR="resources/shet"
INPUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    --input) INPUT="${2:?}"; shift 2 ;;
    -h|--help) echo "Usage: bash scripts/download_shet.sh --out-dir resources/shet --input path/to/weghorn_shet.tsv"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

mkdir -p "${OUT_DIR}"

if [[ -n "${INPUT}" ]]; then
  python3 "${ROOT_DIR}/scripts/normalize_shet.py" --input "${INPUT}" --out "${OUT_DIR}/shet.normalized.tsv"
  exit 0
fi

if [[ -s "${OUT_DIR}/shet.normalized.tsv" ]]; then
  echo "found ${OUT_DIR}/shet.normalized.tsv"
  exit 0
fi

bundle="${OUT_DIR}/UKBBFertility.rawdata.tar.gz"
curl -L --fail -C - -o "${bundle}" \
  "https://raw.githubusercontent.com/HurlesGroupSanger/UKBBFertility/main/rawdata.tar.gz"
tar -xzf "${bundle}" -C "${OUT_DIR}" rawdata/genelists/shet.weghorn.txt
python3 "${ROOT_DIR}/scripts/normalize_shet.py" \
  --input "${OUT_DIR}/rawdata/genelists/shet.weghorn.txt" \
  --out "${OUT_DIR}/shet.normalized.tsv"
