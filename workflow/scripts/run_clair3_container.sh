#!/usr/bin/env bash
set -euo pipefail

container_image=""
bam_fn=""
ref_fn=""
model_path=""
outdir=""
ctg_name=""
threads=""
bind_path=""
platform="ont"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --container-image) container_image="$2"; shift 2 ;;
    --bam-fn) bam_fn="$2"; shift 2 ;;
    --ref-fn) ref_fn="$2"; shift 2 ;;
    --model-path) model_path="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --ctg-name) ctg_name="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --bind-path) bind_path="$2"; shift 2 ;;
    --platform) platform="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

for required in container_image bam_fn ref_fn model_path outdir ctg_name threads; do
  if [[ -z "${!required}" ]]; then
    echo "Missing required argument: ${required}" >&2
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

bind_args=()
if [[ -n "$bind_path" ]]; then
  bind_args=(-B "$bind_path")
fi

mkdir -p "$outdir"

"$runtime" exec --cleanenv "${bind_args[@]}" "$container_image" /bin/bash -s -- \
  "$bam_fn" "$ref_fn" "$model_path" "$outdir" "$ctg_name" "$threads" "$platform" <<'CONTAINER_BASH'
set -euo pipefail

bam_fn="$1"
ref_fn="$2"
model_path="$3"
outdir="$4"
ctg_name="$5"
threads="$6"
platform="$7"

pick_bin() {
  local default_name="$1"
  shift
  if command -v "$default_name" >/dev/null 2>&1; then
    command -v "$default_name"
    return 0
  fi
  for candidate in "$@"; do
    if [[ -x "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

runner=""
runner_mode="exec"
if command -v run_clair3.sh >/dev/null 2>&1; then
  runner="$(command -v run_clair3.sh)"
elif command -v run_clair3 >/dev/null 2>&1; then
  runner="$(command -v run_clair3)"
else
  for p in /opt/bin/run_clair3.sh /opt/bin/scripts/run_clair3.sh /usr/local/bin/run_clair3.sh /opt/conda/bin/run_clair3.sh; do
    if [[ -x "$p" ]]; then
      runner="$p"
      break
    fi
  done
  if [[ -z "$runner" ]]; then
    for p in /opt/bin/run_clair3.py /opt/bin/scripts/run_clair3.py; do
      if [[ -f "$p" ]]; then
        runner="$p"
        runner_mode="python"
        break
      fi
    done
  fi
fi

if [[ -z "$runner" ]]; then
  echo "Could not locate Clair3 launcher in container." >&2
  echo "PATH: $PATH" >&2
  find /opt -maxdepth 3 -name 'run_clair3*' >&2 || true
  exit 1
fi

python_bin=""
for candidate in \
  /opt/conda/envs/clair3/bin/python3 \
  /opt/conda/envs/clair3/bin/python \
  /opt/conda/bin/python3 \
  /opt/conda/bin/python \
  /usr/local/bin/python3 \
  /usr/bin/python3; do
  if [[ -x "$candidate" ]] && "$candidate" -c 'import tensorflow' >/dev/null 2>&1; then
    python_bin="$candidate"
    break
  fi
done

if [[ -z "$python_bin" ]]; then
  python_bin="$(pick_bin python3 /opt/conda/bin/python3 /usr/local/bin/python3 /usr/bin/python3 /bin/python3 || true)"
fi

samtools_bin="$(pick_bin samtools /opt/conda/bin/samtools /usr/local/bin/samtools /usr/bin/samtools /bin/samtools /opt/bin/samtools || true)"
parallel_bin="$(pick_bin parallel /opt/conda/bin/parallel /usr/local/bin/parallel /usr/bin/parallel /bin/parallel /opt/bin/parallel || true)"
whatshap_bin="$(pick_bin whatshap /opt/conda/bin/whatshap /usr/local/bin/whatshap /usr/bin/whatshap /bin/whatshap /opt/bin/whatshap || true)"
pypy_bin="$(pick_bin pypy3 /opt/conda/bin/pypy3 /usr/local/bin/pypy3 /usr/bin/pypy3 /bin/pypy3 /opt/bin/pypy3 || true)"

if [[ -z "$python_bin" || -z "$samtools_bin" || -z "$parallel_bin" || -z "$whatshap_bin" ]]; then
  echo "Missing required binaries inside container." >&2
  echo "python_bin=$python_bin" >&2
  echo "samtools_bin=$samtools_bin" >&2
  echo "parallel_bin=$parallel_bin" >&2
  echo "whatshap_bin=$whatshap_bin" >&2
  exit 1
fi

if [[ -z "$pypy_bin" ]]; then
  pypy_bin="$python_bin"
  echo "[WARNING] pypy3 not found; using python for --pypy" >&2
fi

echo "[INFO] Clair3 runner: $runner (mode=$runner_mode)" >&2
echo "[INFO] python_bin=$python_bin" >&2
echo "[INFO] pypy_bin=$pypy_bin" >&2
echo "[INFO] samtools_bin=$samtools_bin" >&2
echo "[INFO] parallel_bin=$parallel_bin" >&2
echo "[INFO] whatshap_bin=$whatshap_bin" >&2

if [[ "$runner_mode" == "python" ]]; then
  "$python_bin" "$runner" \
    --bam_fn="$bam_fn" \
    --ref_fn="$ref_fn" \
    --output="$outdir" \
    --threads="$threads" \
    --platform="$platform" \
    --model_path="$model_path" \
    --ctg_name="$ctg_name" \
    --python="$python_bin" \
    --samtools="$samtools_bin" \
    --parallel="$parallel_bin" \
    --whatshap="$whatshap_bin" \
    --pypy="$pypy_bin"
else
  "$runner" \
    --bam_fn="$bam_fn" \
    --ref_fn="$ref_fn" \
    --output="$outdir" \
    --threads="$threads" \
    --platform="$platform" \
    --model_path="$model_path" \
    --ctg_name="$ctg_name" \
    --python="$python_bin" \
    --samtools="$samtools_bin" \
    --parallel="$parallel_bin" \
    --whatshap="$whatshap_bin" \
    --pypy="$pypy_bin"
fi
CONTAINER_BASH
