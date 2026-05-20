#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/lib.sh"

OUT_DIR="resources/vep"
VEP_RELEASE="${VEP_RELEASE:-115}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="${2:?}"; shift 2 ;;
    --vep-release) VEP_RELEASE="${2:?}"; shift 2 ;;
    -h|--help) echo "Usage: bash scripts/setup_vep_loftee.sh --out-dir resources/vep --vep-release 115"; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

need_cmd git
mkdir -p "${OUT_DIR}" "resources/loftee"

if [[ ! -d "${OUT_DIR}/ensembl-vep/.git" ]]; then
  git clone https://github.com/Ensembl/ensembl-vep.git "${OUT_DIR}/ensembl-vep"
fi
(
  cd "${OUT_DIR}/ensembl-vep"
  git fetch --tags origin
  git checkout "release/${VEP_RELEASE}"
)

if [[ ! -d "resources/loftee/loftee/.git" ]]; then
  git clone https://github.com/konradjk/loftee.git "resources/loftee/loftee"
fi
(
  cd "resources/loftee/loftee"
  git fetch origin
  git checkout grch38
  git rev-parse HEAD > "../loftee_commit.txt"
)

cat <<EOF
VEP and LOFTEE source checkouts are present.
Install VEP cache/FASTA with Ensembl's INSTALL.pl as appropriate for this machine:
  cd ${OUT_DIR}/ensembl-vep
  perl INSTALL.pl --AUTO acf --SPECIES homo_sapiens --ASSEMBLY GRCh38 --DESTDIR ../cache
EOF
