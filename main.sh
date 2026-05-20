#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT_DIR}/config/default.env"
source "${ROOT_DIR}/scripts/lib.sh"

if [[ -d "${ROOT_DIR}/tools/envs/ptv/bin" ]]; then
  export PATH="${ROOT_DIR}/tools/envs/ptv/bin:${PATH}"
fi

CHROMOSOMES="${DEFAULT_CHROMOSOMES}"
GNOMAD_DIR=""
KGP_DIR="${DEFAULT_KGP_DIR}"
THREADS="${DEFAULT_THREADS}"
FORCE="${FORCE:-0}"
DOWNLOAD_MISSING=0
DOWNLOAD_CADD=0
ENABLE_GTEX="${GTEX_ENABLE_DEFAULT}"
CHECK_ONLY=0
STREAM_1000G=0

usage() {
  cat <<'USAGE'
Usage:
  bash main.sh --gnomad-dir PATH --kgp-dir PATH [options]

Required:
  --gnomad-dir PATH        Directory containing gnomAD v4.1 joint VCFs.
  --kgp-dir PATH           Directory for 1000G raw VCF cache/downloads. Raw files live under PATH/raw.

Options:
  --chromosomes LIST       Chromosomes, e.g. 22, 1,2,22, or 1-22. Default: 1-22.
  --threads N              Threads for tools that support threading. Default: 4.
  --download-missing       Download non-CADD public resources that are absent.
  --download-cadd          Allow CADD v1.7 resource download. This is large.
  --enable-gtex            Download/use GTEx expression resources.
  --stream-1000g           Read indexed public 1000G VCFs remotely instead of requiring local raw VCFs.
  --check-only             Validate arguments/resources and write tool version report only.
  --force                  Re-run stages even when .done sentinels exist.
  -h, --help               Show this help.
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
    --stream-1000g) STREAM_1000G=1; shift ;;
    --check-only) CHECK_ONLY=1; shift ;;
    --force) FORCE=1; shift ;;
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
export CHROMOSOMES GNOMAD_DIR THREADS FORCE CADD_VERSION VEP_RELEASE GENCODE_RELEASE
export ENABLE_GTEX GNOMAD_ABSENT_POLICY GNOMAD_PRIMARY_AF_FIELDS CADD_THRESHOLD
export MAX_EXTERNAL_LOOKUP_AC CHECK_ONLY STREAM_1000G STREAM_1000G_MERGE_GAP

mkdir -p "${RESOURCES_DIR}" "${WORK_DIR}" "${RESULTS_DIR}" "${LOGS_DIR}" "${DOWNLOADS_DIR}" \
  "${RESULTS_DIR}/qc" "${RESULTS_DIR}/reference_package" "${WORK_DIR}/.done" "${KGP_RAW_DIR}"

CHROM_LIST="$(expand_chromosomes "${CHROMOSOMES}")"
export CHROM_LIST

run_stage tool_versions bash "${ROOT_DIR}/scripts/write_tool_versions.sh"

if [[ "${DOWNLOAD_MISSING}" == "1" ]]; then
  if [[ "${STREAM_1000G}" == "1" ]]; then
    run_stage download_1000g_metadata bash "${ROOT_DIR}/scripts/download_1000g.sh" \
      --metadata-only --out-dir "${KGP_RAW_DIR}"
  else
    run_stage download_1000g bash "${ROOT_DIR}/scripts/download_1000g.sh" \
      --chromosomes "${CHROMOSOMES}" --out-dir "${KGP_RAW_DIR}"
  fi
  run_stage download_reference bash "${ROOT_DIR}/scripts/download_reference.sh" \
    --out-dir "${RESOURCES_DIR}/reference" --gencode-release "${GENCODE_RELEASE}"
  run_stage download_constraint bash "${ROOT_DIR}/scripts/download_constraint.sh" \
    --out-dir "${RESOURCES_DIR}/constraint"
  run_stage download_shet bash "${ROOT_DIR}/scripts/download_shet.sh" \
    --out-dir "${RESOURCES_DIR}/shet"
  run_stage download_loftee_aux bash "${ROOT_DIR}/scripts/download_loftee_aux.sh" \
    --out-dir "${RESOURCES_DIR}/loftee/loftee_data/GRCh38"
  if [[ "${ENABLE_GTEX}" == "1" ]]; then
    run_stage download_gtex bash "${ROOT_DIR}/scripts/download_gtex.sh" \
      --out-dir "${RESOURCES_DIR}/gtex"
  fi
else
  log "download-missing disabled; validating already-present resources"
fi

if [[ "${DOWNLOAD_CADD}" == "1" ]]; then
  run_stage download_cadd bash "${ROOT_DIR}/scripts/download_cadd.sh" \
    --out-dir "${RESOURCES_DIR}/cadd/v${CADD_VERSION}"
fi

if [[ "${CHECK_ONLY}" == "1" ]]; then
  bash "${ROOT_DIR}/scripts/validate_inputs.sh" 2>&1 | tee "${LOGS_DIR}/validate_check_only.log"
  log "check-only requested; stopping after validation"
  exit 0
fi

run_stage validate_inputs bash "${ROOT_DIR}/scripts/validate_inputs.sh"

run_stage build_sample_sets python3 "${ROOT_DIR}/scripts/build_sample_sets.py" \
  --metadata "${KGP_RAW_DIR}/20130606_g1k_3202_samples_ped_population.txt" \
  --out-dir "${RESULTS_DIR}/sample_sets"

run_stage build_coding_bed python3 "${ROOT_DIR}/scripts/build_coding_bed.py" \
  --gtf "${RESOURCES_DIR}/reference/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" \
  --out "${RESOURCES_DIR}/reference/protein_coding_cds_splice.gencode_v${GENCODE_RELEASE}.bed.gz"

run_stage build_gene_annotations python3 "${ROOT_DIR}/scripts/build_gene_annotations.py" \
  --constraint "${RESOURCES_DIR}/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz" \
  --shet "${RESOURCES_DIR}/shet/shet.normalized.tsv" \
  --gtex "${RESOURCES_DIR}/gtex/gene_cns_expression_percentiles.tsv" \
  --out "${RESULTS_DIR}/reference_package/gene_scores_shet_constraint_expression.tsv"

for CHR in ${CHROM_LIST}; do
  export CHR
  run_stage "extract_candidates_chr${CHR}" bash "${ROOT_DIR}/scripts/extract_candidates.sh" "${CHR}"
  run_stage "run_vep_loftee_chr${CHR}" bash "${ROOT_DIR}/scripts/run_vep_loftee.sh" "${CHR}"
  run_stage "parse_vep_chr${CHR}" python3 "${ROOT_DIR}/scripts/parse_vep_csq.py" \
    --vcf "${WORK_DIR}/chr${CHR}/vep_loftee.chr${CHR}.vcf.gz" \
    --out "${WORK_DIR}/chr${CHR}/vep_loftee.parsed.chr${CHR}.tsv.gz"
  run_stage "compute_1000g_ac_chr${CHR}" python3 "${ROOT_DIR}/scripts/compute_1000g_ac.py" \
    --vcf "${WORK_DIR}/chr${CHR}/candidates.norm.chr${CHR}.vcf.gz" \
    --parsed "${WORK_DIR}/chr${CHR}/vep_loftee.parsed.chr${CHR}.tsv.gz" \
    --samples "${RESULTS_DIR}/sample_sets/eur_unrelated.samples.txt" \
    --out-variants "${WORK_DIR}/chr${CHR}/ptv_ac.chr${CHR}.tsv.gz" \
    --out-carriers "${WORK_DIR}/chr${CHR}/ptv_carriers.chr${CHR}.tsv.gz"
  run_stage "prefilter_external_chr${CHR}" python3 "${ROOT_DIR}/scripts/prefilter_external_lookup.py" \
    --annotations "${WORK_DIR}/chr${CHR}/ptv_ac.chr${CHR}.tsv.gz" \
    --max-ac "${MAX_EXTERNAL_LOOKUP_AC}" \
    --out "${WORK_DIR}/chr${CHR}/external_lookup_input.chr${CHR}.tsv.gz"
  run_stage "join_gnomad_chr${CHR}" python3 "${ROOT_DIR}/scripts/join_gnomad_af.py" \
    --variants "${WORK_DIR}/chr${CHR}/external_lookup_input.chr${CHR}.tsv.gz" \
    --gnomad-dir "${GNOMAD_DIR}" \
    --chromosome "${CHR}" \
    --out "${WORK_DIR}/chr${CHR}/ptv_gnomad.chr${CHR}.tsv.gz" \
    --field-report "${RESULTS_DIR}/qc/gnomad_info_fields.chr${CHR}.tsv"
  run_stage "prefilter_cadd_chr${CHR}" python3 "${ROOT_DIR}/scripts/prefilter_cadd_lookup.py" \
    --variants "${WORK_DIR}/chr${CHR}/ptv_gnomad.chr${CHR}.tsv.gz" \
    --max-gnomad-af 0.001 \
    --out "${WORK_DIR}/chr${CHR}/cadd_lookup_input.chr${CHR}.tsv.gz"
  run_stage "join_cadd_chr${CHR}" python3 "${ROOT_DIR}/scripts/join_cadd.py" \
    --variants "${WORK_DIR}/chr${CHR}/cadd_lookup_input.chr${CHR}.tsv.gz" \
    --cadd-dir "${RESOURCES_DIR}/cadd/v${CADD_VERSION}" \
    --out "${WORK_DIR}/chr${CHR}/ptv_cadd.chr${CHR}.tsv.gz"
  run_stage "apply_tiers_chr${CHR}" python3 "${ROOT_DIR}/scripts/apply_ptv_tiers.py" \
    --annotations "${WORK_DIR}/chr${CHR}/ptv_cadd.chr${CHR}.tsv.gz" \
    --carriers "${WORK_DIR}/chr${CHR}/ptv_carriers.chr${CHR}.tsv.gz" \
    --gene-scores "${RESULTS_DIR}/reference_package/gene_scores_shet_constraint_expression.tsv" \
    --tiers "${ROOT_DIR}/config/tiers.tsv" \
    --out "${WORK_DIR}/chr${CHR}/qualifying_variants_by_tier.chr${CHR}.tsv.gz"
done

run_stage merge_variant_catalog python3 "${ROOT_DIR}/scripts/merge_tables.py" \
  --inputs $(for CHR in ${CHROM_LIST}; do printf "%s " "${WORK_DIR}/chr${CHR}/qualifying_variants_by_tier.chr${CHR}.tsv.gz"; done) \
  --out "${RESULTS_DIR}/reference_package/qualifying_variants_by_tier.tsv.gz"

run_stage compute_burden python3 "${ROOT_DIR}/scripts/compute_burden.py" \
  --variants "${RESULTS_DIR}/reference_package/qualifying_variants_by_tier.tsv.gz" \
  --samples "${RESULTS_DIR}/sample_sets/eur_unrelated.metadata.tsv" \
  --tiers "${ROOT_DIR}/config/tiers.tsv" \
  --out "${RESULTS_DIR}/reference_package/burden_by_sample.tsv.gz"

run_stage make_qc_report python3 "${ROOT_DIR}/scripts/make_qc_report.py" \
  --variants "${RESULTS_DIR}/reference_package/qualifying_variants_by_tier.tsv.gz" \
  --burden "${RESULTS_DIR}/reference_package/burden_by_sample.tsv.gz" \
  --out-dir "${RESULTS_DIR}/qc"

run_stage summarize_reference_distribution python3 "${ROOT_DIR}/scripts/summarize_reference_distribution.py" \
  --qualifying "${RESULTS_DIR}/reference_package/qualifying_variants_by_tier.tsv.gz" \
  --burden "${RESULTS_DIR}/reference_package/burden_by_sample.tsv.gz" \
  --out-prefix "${RESULTS_DIR}/reference_package/reference_distribution"

run_stage build_reference_package python3 "${ROOT_DIR}/scripts/build_reference_package.py" \
  --results-dir "${RESULTS_DIR}" \
  --resources-dir "${RESOURCES_DIR}" \
  --chromosomes "${CHROMOSOMES}" \
  --gnomad-dir "${GNOMAD_DIR}" \
  --kgp-dir "${KGP_DIR}" \
  --out-dir "${RESULTS_DIR}/reference_package"

log "done"
