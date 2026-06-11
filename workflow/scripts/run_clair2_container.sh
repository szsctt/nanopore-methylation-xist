#!/usr/bin/env bash
set -euo pipefail

container_image=""
bam_fn=""
ref_fn=""
model_dir=""
out_vcfgz=""
outdir=""
ctg_name=""
sample_name=""
threads=""
threshold=""
min_coverage=""
extra_args=""
bind_path=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --container-image) container_image="$2"; shift 2 ;;
    --bam-fn) bam_fn="$2"; shift 2 ;;
    --ref-fn) ref_fn="$2"; shift 2 ;;
    --model-dir) model_dir="$2"; shift 2 ;;
    --out-vcfgz) out_vcfgz="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --ctg-name) ctg_name="$2"; shift 2 ;;
    --sample-name) sample_name="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --threshold) threshold="$2"; shift 2 ;;
    --min-coverage) min_coverage="$2"; shift 2 ;;
    --extra-args) extra_args="$2"; shift 2 ;;
    --bind-path) bind_path="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

for required in container_image bam_fn ref_fn model_dir out_vcfgz outdir ctg_name sample_name threads threshold min_coverage; do
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

raw_vcf="$outdir/$sample_name.vcf"
mkdir -p "$outdir"

"$runtime" exec --cleanenv "${bind_opt[@]}" "$container_image" /bin/bash -s -- \
  "$bam_fn" "$ref_fn" "$model_dir" "$raw_vcf" "$ctg_name" "$sample_name" "$threads" "$threshold" "$min_coverage" "$extra_args" <<'CONTAINER_BASH'
set -euo pipefail

bam_fn="$1"
ref_fn="$2"
model_dir="$3"
raw_vcf="$4"
ctg_name="$5"
sample_name="$6"
threads="$7"
threshold="$8"
min_coverage="$9"
extra_args="${10}"

clair_bin=$(command -v clair.py || true)
python_bin=$(command -v python3 || command -v python || true)
samtools_bin=$(command -v samtools || true)
pypy_bin=$(command -v pypy3 || command -v pypy || command -v python3 || command -v python || true)

if [[ -z "$clair_bin" || -z "$python_bin" || -z "$samtools_bin" || -z "$pypy_bin" ]]; then
  echo "Missing required binary in Clair2 container." >&2
  echo "clair_bin=$clair_bin" >&2
  echo "python_bin=$python_bin" >&2
  echo "samtools_bin=$samtools_bin" >&2
  echo "pypy_bin=$pypy_bin" >&2
  exit 1
fi

if ! "$pypy_bin" -c 'import intervaltree' >/dev/null 2>&1; then
  echo "[INFO] pypy missing intervaltree; trying PYTHONPATH bridge." >&2
  pypy_py_path=""
  for site_dir in /opt/conda/envs/*/lib/python*/site-packages /usr/local/lib/python*/site-packages /opt/conda/lib/python*/site-packages /usr/lib/python*/site-packages; do
    if [[ -d "$site_dir" ]] && PYTHONPATH="$site_dir" "$pypy_bin" -c 'import intervaltree' >/dev/null 2>&1; then
      pypy_py_path="$site_dir"
      break
    fi
  done

  if [[ -n "$pypy_py_path" ]]; then
    export PYTHONPATH="$pypy_py_path"
    echo "[INFO] using pypy via PYTHONPATH bridge: $pypy_py_path" >&2
  else
    pypy_bin="$python_bin"
    echo "[WARNING] pypy intervaltree unavailable; falling back to python for --pypy." >&2
  fi
else
  echo "[INFO] pypy intervaltree import OK with default environment." >&2
fi

echo "[INFO] --pypy binary set to: $pypy_bin" >&2

cmd=(
  "$python_bin" "$clair_bin" callVarBam
  --chkpnt_fn "$model_dir/model"
  --ref_fn "$ref_fn"
  --bam_fn "$bam_fn"
  --call_fn "$raw_vcf"
  --ctgName "$ctg_name"
  --sampleName "$sample_name"
  --threads "$threads"
  --threshold="$threshold"
  --minCoverage="$min_coverage"
  --samtools "$samtools_bin"
  --pypy "$pypy_bin"
)

if [[ -n "$extra_args" ]]; then
  # shellcheck disable=SC2206
  extra_parts=($extra_args)
  cmd+=("${extra_parts[@]}")
fi

"${cmd[@]}"
CONTAINER_BASH

if [[ ! -s "$raw_vcf" ]]; then
  echo "Expected Clair2 output VCF missing: $raw_vcf" >&2
  find "$outdir" -maxdepth 3 -type f >&2 || true
  exit 1
fi

if command -v bgzip >/dev/null 2>&1; then
  bgzip -f -c "$raw_vcf" > "$out_vcfgz"
else
  gzip -f -c "$raw_vcf" > "$out_vcfgz"
fi
