#!/usr/bin/env bash
set -euo pipefail

container_image=""
in_sam=""
out_bam=""
out_bai=""
threads="8"
bind_path=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --container-image) container_image="$2"; shift 2 ;;
    --in-sam) in_sam="$2"; shift 2 ;;
    --out-bam) out_bam="$2"; shift 2 ;;
    --out-bai) out_bai="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --bind-path) bind_path="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

for required in container_image in_sam out_bam out_bai threads; do
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

mkdir -p "$(dirname "$out_bam")"

"$runtime" exec --cleanenv "${bind_opt[@]}" "$container_image" /bin/bash -s -- \
  "$in_sam" "$out_bam" "$out_bai" "$threads" <<'CONTAINER_BASH'
set -euo pipefail

in_sam="$1"
out_bam="$2"
out_bai="$3"
threads="$4"

samtools sort -@ "$threads" -o "$out_bam" "$in_sam"
samtools index -@ "$threads" "$out_bam" "$out_bai"
CONTAINER_BASH
