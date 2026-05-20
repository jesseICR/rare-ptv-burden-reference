#!/usr/bin/env bash
set -euo pipefail

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

expand_chromosomes() {
  local spec="$1"
  local part start end i out=()
  IFS=',' read -ra parts <<< "${spec}"
  for part in "${parts[@]}"; do
    if [[ "${part}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      start="${BASH_REMATCH[1]}"
      end="${BASH_REMATCH[2]}"
      for ((i=start; i<=end; i++)); do out+=("${i}"); done
    elif [[ "${part}" =~ ^([0-9]+|X|Y|MT|M)$ ]]; then
      out+=("${part}")
    else
      die "Invalid chromosome spec element: ${part}"
    fi
  done
  printf '%s\n' "${out[@]}"
}

run_stage() {
  local stage="$1"
  shift
  local done_file="${WORK_DIR}/.done/${stage}.done"
  local log_file="${LOGS_DIR}/${stage}.log"

  if [[ -f "${done_file}" && "${FORCE:-0}" != "1" ]]; then
    log "[skip] ${stage}"
    return 0
  fi

  log "[run] ${stage}"
  mkdir -p "$(dirname "${done_file}")" "$(dirname "${log_file}")"
  "$@" 2>&1 | tee "${log_file}"
  touch "${done_file}"
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

chrom_label() {
  local chr="$1"
  case "${chr}" in
    chr*) printf '%s\n' "${chr}" ;;
    *) printf 'chr%s\n' "${chr}" ;;
  esac
}
