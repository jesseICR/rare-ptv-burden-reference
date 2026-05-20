#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/config/default.env"
source "${ROOT_DIR}/scripts/lib.sh"

if [[ -d "${ROOT_DIR}/tools/envs/ptv/bin" ]]; then
  export PATH="${ROOT_DIR}/tools/envs/ptv/bin:${PATH}"
fi

GNOMAD_DIR=""
KGP_DIR="${DEFAULT_KGP_DIR}"
CHROMOSOMES="${DEFAULT_CHROMOSOMES}"
THREADS="${DEFAULT_THREADS}"
DOWNLOAD_MISSING=0
DOWNLOAD_CADD=0
ENABLE_GTEX="${GTEX_ENABLE_DEFAULT}"
DELETE_RAW_AFTER_SUCCESS=0
REFRESH_CADD=0
STREAM_1000G=0

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_chromosomes_sequential.sh --gnomad-dir PATH --kgp-dir PATH [options]

Required:
  --gnomad-dir PATH              Directory containing gnomAD v4.1 joint VCFs.
  --kgp-dir PATH                 Directory for 1000G raw VCF cache/downloads. Raw files live under PATH/raw.

Options:
  --chromosomes LIST             Chromosomes, e.g. 22, 1,2,22, or 1-22. Default: 1-22.
  --threads N                    Threads for tools that support threading. Default: 4.
  --download-missing             Download public non-CADD resources and each missing 1000G chromosome VCF as needed.
  --download-cadd                Download local CADD v1.7 SNV and indel tables.
  --enable-gtex                  Download/use GTEx expression resources.
  --delete-raw-after-success     Delete each raw 1000G chromosome VCF after that chromosome finishes.
  --refresh-cadd                 Recompute CADD joins and tier outputs for chromosomes that already have earlier outputs.
  --stream-1000g                 Read indexed public 1000G VCFs remotely instead of downloading raw chromosome VCFs.
  -h, --help                     Show this help.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gnomad-dir) GNOMAD_DIR="${2:?}"; shift 2 ;;
    --kgp-dir) KGP_DIR="${2:?}"; shift 2 ;;
    --chromosomes) CHROMOSOMES="${2:?}"; shift 2 ;;
    --threads) THREADS="${2:?}"; shift 2 ;;
    --download-missing) DOWNLOAD_MISSING=1; shift ;;
    --download-cadd) DOWNLOAD_CADD=1; shift ;;
    --enable-gtex) ENABLE_GTEX=1; shift ;;
    --delete-raw-after-success) DELETE_RAW_AFTER_SUCCESS=1; shift ;;
    --refresh-cadd) REFRESH_CADD=1; shift ;;
    --stream-1000g) STREAM_1000G=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "${GNOMAD_DIR}" ]] || { usage >&2; die "--gnomad-dir is required"; }
[[ -n "${KGP_DIR}" ]] || { usage >&2; die "--kgp-dir is required"; }
KGP_DIR="${KGP_DIR%/}"
KGP_RAW_DIR="${KGP_DIR}/raw"

export ROOT_DIR PROJECT_NAME GENOME_BUILD GNOMAD_VERSION
export RESOURCES_DIR WORK_DIR RESULTS_DIR LOGS_DIR DOWNLOADS_DIR KGP_DIR KGP_RAW_DIR KGP_1000G_RAW_BASE
export GNOMAD_DIR THREADS FORCE CADD_VERSION VEP_RELEASE GENCODE_RELEASE
export ENABLE_GTEX GNOMAD_ABSENT_POLICY GNOMAD_PRIMARY_AF_FIELDS CADD_THRESHOLD
export MAX_EXTERNAL_LOOKUP_AC STREAM_1000G STREAM_1000G_MERGE_GAP

mkdir -p "${RESOURCES_DIR}" "${WORK_DIR}" "${RESULTS_DIR}" "${LOGS_DIR}" "${DOWNLOADS_DIR}" \
  "${RESULTS_DIR}/qc" "${RESULTS_DIR}/reference_package" "${WORK_DIR}/.done" "${KGP_RAW_DIR}"

CHROM_LIST="$(expand_chromosomes "${CHROMOSOMES}")"

if [[ "${DOWNLOAD_MISSING}" == "1" ]]; then
  bash "${ROOT_DIR}/scripts/download_1000g.sh" \
    --metadata-only --out-dir "${KGP_RAW_DIR}"
  bash "${ROOT_DIR}/scripts/download_reference.sh" \
    --out-dir "${RESOURCES_DIR}/reference" --gencode-release "${GENCODE_RELEASE}"
  bash "${ROOT_DIR}/scripts/download_constraint.sh" \
    --out-dir "${RESOURCES_DIR}/constraint"
  bash "${ROOT_DIR}/scripts/download_shet.sh" \
    --out-dir "${RESOURCES_DIR}/shet"
  bash "${ROOT_DIR}/scripts/download_loftee_aux.sh" \
    --out-dir "${RESOURCES_DIR}/loftee/loftee_data/GRCh38"
  if [[ "${ENABLE_GTEX}" == "1" ]]; then
    bash "${ROOT_DIR}/scripts/download_gtex.sh" \
      --out-dir "${RESOURCES_DIR}/gtex"
  fi
fi

if [[ "${DOWNLOAD_CADD}" == "1" ]]; then
  bash "${ROOT_DIR}/scripts/download_cadd.sh" \
    --out-dir "${RESOURCES_DIR}/cadd/v${CADD_VERSION}"
fi

for CHR in ${CHROM_LIST}; do
  chr_variants="${WORK_DIR}/chr${CHR}/qualifying_variants_by_tier.chr${CHR}.tsv.gz"
  if [[ "${REFRESH_CADD}" != "1" && -s "${chr_variants}" ]]; then
    log "[skip] chr${CHR} already has ${chr_variants}"
    continue
  fi

  raw_vcf="${KGP_RAW_DIR}/chr${CHR}.vcf.gz"
  if [[ "${STREAM_1000G}" != "1" && ! -s "${raw_vcf}" ]]; then
    if [[ "${DOWNLOAD_MISSING}" != "1" ]]; then
      die "Missing raw 1000G VCF ${raw_vcf}; rerun with --download-missing or download it first."
    fi
    bash "${ROOT_DIR}/scripts/download_1000g.sh" \
      --chromosomes "${CHR}" --out-dir "${KGP_RAW_DIR}"
  fi

  if [[ "${REFRESH_CADD}" == "1" ]]; then
    rm -f \
      "${WORK_DIR}/.done/join_cadd_chr${CHR}.done" \
      "${WORK_DIR}/.done/apply_tiers_chr${CHR}.done" \
      "${WORK_DIR}/.done/merge_variant_catalog.done" \
      "${WORK_DIR}/.done/compute_burden.done" \
      "${WORK_DIR}/.done/make_qc_report.done" \
      "${WORK_DIR}/.done/summarize_reference_distribution.done" \
      "${WORK_DIR}/.done/build_reference_package.done"
  fi

  rm -f "${WORK_DIR}/.done/validate_inputs.done"
  args=(--chromosomes "${CHR}" --gnomad-dir "${GNOMAD_DIR}" --kgp-dir "${KGP_DIR}" --threads "${THREADS}")
  if [[ "${ENABLE_GTEX}" == "1" ]]; then
    args+=(--enable-gtex)
  fi
  if [[ "${STREAM_1000G}" == "1" ]]; then
    args+=(--stream-1000g)
  fi
  bash "${ROOT_DIR}/main.sh" "${args[@]}"

  if [[ "${DELETE_RAW_AFTER_SUCCESS}" == "1" && "${STREAM_1000G}" != "1" ]]; then
    rm -f "${raw_vcf}" "${raw_vcf}.tbi" "${raw_vcf}.csi"
  fi
done

inputs=()
for CHR in ${CHROM_LIST}; do
  chr_variants="${WORK_DIR}/chr${CHR}/qualifying_variants_by_tier.chr${CHR}.tsv.gz"
  [[ -s "${chr_variants}" ]] || die "Missing chromosome tier output ${chr_variants}"
  inputs+=("${chr_variants}")
done

python3 "${ROOT_DIR}/scripts/merge_tables.py" \
  --inputs "${inputs[@]}" \
  --out "${RESULTS_DIR}/reference_package/qualifying_variants_by_tier.tsv.gz"
python3 "${ROOT_DIR}/scripts/compute_burden.py" \
  --variants "${RESULTS_DIR}/reference_package/qualifying_variants_by_tier.tsv.gz" \
  --samples "${RESULTS_DIR}/sample_sets/eur_unrelated.metadata.tsv" \
  --tiers "${ROOT_DIR}/config/tiers.tsv" \
  --out "${RESULTS_DIR}/reference_package/burden_by_sample.tsv.gz"
python3 "${ROOT_DIR}/scripts/make_qc_report.py" \
  --variants "${RESULTS_DIR}/reference_package/qualifying_variants_by_tier.tsv.gz" \
  --burden "${RESULTS_DIR}/reference_package/burden_by_sample.tsv.gz" \
  --out-dir "${RESULTS_DIR}/qc"
python3 "${ROOT_DIR}/scripts/summarize_reference_distribution.py" \
  --qualifying "${RESULTS_DIR}/reference_package/qualifying_variants_by_tier.tsv.gz" \
  --burden "${RESULTS_DIR}/reference_package/burden_by_sample.tsv.gz" \
  --out-prefix "${RESULTS_DIR}/reference_package/reference_distribution"
python3 "${ROOT_DIR}/scripts/build_reference_package.py" \
  --results-dir "${RESULTS_DIR}" \
  --resources-dir "${RESOURCES_DIR}" \
  --chromosomes "${CHROMOSOMES}" \
  --gnomad-dir "${GNOMAD_DIR}" \
  --kgp-dir "${KGP_DIR}" \
  --out-dir "${RESULTS_DIR}/reference_package"

log "sequential chromosome run done"
