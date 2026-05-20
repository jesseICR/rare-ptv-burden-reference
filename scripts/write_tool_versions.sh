#!/usr/bin/env bash
set -euo pipefail

mkdir -p "${RESULTS_DIR}/reference_package"
out="${RESULTS_DIR}/reference_package/tool_versions.tsv"
{
  printf 'tool\tversion\n'
  printf 'pipeline\t%s\n' "${PROJECT_NAME:-kgp-ptv-shet-burden}"
  printf 'bash\t%s\n' "$(bash --version | head -1)"
  printf 'python\t%s\n' "$(python3 --version 2>&1)"
  if command -v bcftools >/dev/null 2>&1; then
    printf 'bcftools\t%s\n' "$(bcftools --version | head -1)"
  fi
  if command -v tabix >/dev/null 2>&1; then
    printf 'tabix\t%s\n' "$(tabix --version 2>&1 | head -1)"
  fi
  if command -v samtools >/dev/null 2>&1; then
    printf 'samtools\t%s\n' "$(samtools --version 2>&1 | head -1)"
  fi
  if command -v plink2 >/dev/null 2>&1; then
    printf 'plink2\t%s\n' "$(plink2 --version 2>&1 | head -1)"
  fi
  if command -v vep >/dev/null 2>&1; then
    printf 'vep\t%s\n' "$(vep --help 2>&1 | head -1)"
  fi
  printf 'cadd\tv%s\n' "${CADD_VERSION:-1.7}"
  printf 'gnomad\tv%s joint\n' "${GNOMAD_VERSION:-4.1}"
} > "${out}"
echo "wrote ${out}"
