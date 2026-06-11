#!/usr/bin/env bash
set -euo pipefail

container_image=""
ref_fn=""
fastq_fn=""
out_sam=""
preset="map-ont"
threads="8"
bind_path=""
extra_args=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --container-image) container_image="$2"; shift 2 ;;
    --ref-fn) ref_fn="$2"; shift 2 ;;
    --fastq-fn) fastq_fn="$2"; shift 2 ;;
    --out-sam) out_sam="$2"; shift 2 ;;
    --preset) preset="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --bind-path) bind_path="$2"; shift 2 ;;
    --extra-args) extra_args="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

for required in container_image ref_fn fastq_fn out_sam preset threads; do
  if [[ -z "${!required}" ]]; then
    echo "Missing required argument: $required" >&2
    exit 2
  fi
done

runtime=""
if command -v apptainer >/dev/null 2>&1; then
  runtime="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  runtime="singularity"
else
  echo "Neither apptainer nor singularity is available in PATH" >&2
  exit 1
fi

bind_opt=()
if [[ -n "$bind_path" ]]; then
  bind_opt=(-B "$bind_path")
fi

mkdir -p "$(dirname "$out_sam")"

"$runtime" exec --cleanenv "${bind_opt[@]}" "$container_image" /bin/bash -s -- \
  "$ref_fn" "$fastq_fn" "$out_sam" "$preset" "$threads" "$extra_args" <<'CONTAINER_BASH'
set -euo pipefail

ref_fn="$1"
fastq_fn="$2"
out_sam="$3"
preset="$4"
threads="$5"
extra_args="$6"

cmd=(minimap2 -ax "$preset" -t "$threads")
if [[ -n "$extra_args" ]]; then
  # shellcheck disable=SC2206
  extra_parts=($extra_args)
  cmd+=("${extra_parts[@]}")
fi
cmd+=("$ref_fn" "$fastq_fn")

"${cmd[@]}" > "$out_sam"
CONTAINER_BASH
