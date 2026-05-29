import csv

configfile: "config/clair2_benchmark.yaml"

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")


def load_samples(tsv_path):
    samples = {}
    with open(tsv_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {tsv_path}: {', '.join(sorted(missing))}")

        for row in reader:
            sample = row["sample"].strip()
            if not sample:
                continue
            samples[sample] = {k: (v or "").strip() for k, v in row.items()}

    if not samples:
        raise ValueError(f"No samples found in {tsv_path}")

    return samples


SAMPLES = load_samples(config["sample_table"])
ALL_SAMPLE_NAMES = sorted(SAMPLES.keys())
BENCHMARK_SAMPLES = config.get("benchmark_samples", ALL_SAMPLE_NAMES)
XIST = config["xist_region"]
THREADS = config["threads"]
VERSIONS = config["versions"]
VERSION_NAMES = list(VERSIONS)
BENCHMARK_NAME = config.get("benchmark_name", "clair2_xist")
RESULTS_BASE = f"results/{{sample}}/benchmarks/{BENCHMARK_NAME}"
WORK_BASE = f"work/benchmark_{BENCHMARK_NAME}"
LOG_BASE = f"logs/benchmark_{BENCHMARK_NAME}"
MODEL_NAME = config["clair2_model"]["name"]
MODEL_ARCHIVE = f"resources/models/{BENCHMARK_NAME}/{MODEL_NAME}.tar"
MODEL_DIR = f"resources/models/{BENCHMARK_NAME}/{MODEL_NAME}"


def assert_sample_config():
    unknown = sorted(set(BENCHMARK_SAMPLES) - set(ALL_SAMPLE_NAMES))
    if unknown:
        raise ValueError(f"benchmark_samples contains unknown sample(s): {', '.join(unknown)}")


assert_sample_config()


def truth_xist_path(wildcards):
    template = config["truth"]["xist_vcf"]
    return template.format(sample=wildcards.sample)


rule all:
    input:
        expand(
            f"{RESULTS_BASE}/{{version}}/{{sample}}.comparison.tsv",
            sample=BENCHMARK_SAMPLES,
            version=VERSION_NAMES,
        ),
        expand(
            f"{RESULTS_BASE}/{{sample}}.ranking.tsv",
            sample=BENCHMARK_SAMPLES,
        ),


rule create_xist_bed:
    output:
        bed=f"{WORK_BASE}/regions/xist.bed.gz",
    params:
        chrom=XIST["chrom"],
        start=XIST["start"],
        end=XIST["end"],
    shell:
        """
        mkdir -p $(dirname {output.bed})
        # BED is 0-based half-open; config stores 1-based inclusive coordinates.
        # Clair v2 reads BED via `gzip -fdc`, so this must be gzipped.
        printf "%s\t%d\t%d\n" {params.chrom} $(({params.start}-1)) {params.end} | gzip -c > {output.bed}
        """


rule index_reference:
    input:
        ref=config["reference_fasta"],
    output:
        fai=config["reference_fasta"] + ".fai",
    container:
        config["containers"]["samtools"]
    log:
        f"{LOG_BASE}/index_reference.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          samtools faidx {input.ref}
        ) > {log} 2>&1
        """


rule download_clair2_model:
    output:
        archive=MODEL_ARCHIVE,
    params:
        url=config["clair2_model"]["source_url"],
    container:
        config["containers"]["download"]
    log:
        f"{LOG_BASE}/download_model.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.archive})
          model_url="{params.url}"
          if [[ -z "$model_url" ]]; then
            echo "Missing clair2_model.source_url in config/clair2_benchmark.yaml" >&2
            exit 1
          fi

          if command -v curl >/dev/null 2>&1; then
            curl -L --fail --retry 5 --retry-delay 10 -o {output.archive} "$model_url"
          else
            wget -O {output.archive} "$model_url"
          fi
        ) > {log} 2>&1
        """


rule prepare_clair2_model:
    input:
        archive=rules.download_clair2_model.output.archive,
    output:
        model_dir=directory(MODEL_DIR),
    container:
        config["containers"]["download"]
    log:
        f"{LOG_BASE}/prepare_model.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          tmpdir=$(mktemp -d)
          trap 'rm -rf "$tmpdir"' EXIT

          tar -xf {input.archive} -C "$tmpdir"
          model_index=$(find "$tmpdir" -type f -name "model.index" | head -n 1 || true)
          if [[ -z "$model_index" ]]; then
            echo "Could not locate model.index in archive: {input.archive}" >&2
            find "$tmpdir" -maxdepth 3 -type f >&2 || true
            exit 1
          fi

          model_candidate=$(dirname "$model_index")
          rm -rf {output.model_dir}
          mkdir -p {output.model_dir}
          cp -a "$model_candidate"/. {output.model_dir}/

          if [[ ! -f {output.model_dir}/model.index ]]; then
            echo "Prepared model is missing model.index at {output.model_dir}" >&2
            exit 1
          fi
        ) > {log} 2>&1
        """


rule call_clair2_xist:
    input:
        bam=lambda wc: f"results/{wc.sample}/alignment/{wc.sample}.sorted.bam",
        bai=lambda wc: f"results/{wc.sample}/alignment/{wc.sample}.sorted.bam.bai",
        ref=config["reference_fasta"],
        fai=rules.index_reference.output.fai,
        model=rules.prepare_clair2_model.output.model_dir,
        bed=rules.create_xist_bed.output.bed,
    output:
        vcf=f"{RESULTS_BASE}/{{version}}/{{sample}}.vcf.gz",
    params:
        outdir=lambda wc: f"{WORK_BASE}/{wc.sample}/{wc.version}",
        out_vcf=lambda wc: f"{WORK_BASE}/{wc.sample}/{wc.version}/{wc.sample}.clair2.vcf",
        ctg=XIST["chrom"],
        ctg_start=XIST["start"],
        ctg_end=XIST["end"],
        threshold=config["clair2"].get("threshold", 0.125),
        min_coverage=config["clair2"].get("min_coverage", 4),
        extra_args=config["clair2"].get("extra_args", ""),
    threads: THREADS["clair2"]
    container:
        config["containers"]["clair2"]
    log:
        f"{LOG_BASE}/{{sample}}.{{version}}.call.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.vcf})
          rm -rf {params.outdir}
          mkdir -p {params.outdir}

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

          # Some Clair2 containers include pypy3 without intervaltree.
          # Try exposing Python3 site-packages to pypy first (no downloads).
          # Fall back to python for pypy-subtasks only if still unavailable.
          if ! "$pypy_bin" -c 'import intervaltree' >/dev/null 2>&1; then
            echo "[INFO] pypy missing intervaltree; trying PYTHONPATH bridge." >&2
            pypy_py_path=""
            for site_dir in /usr/local/lib/python*/site-packages /opt/conda/lib/python*/site-packages /usr/lib/python*/site-packages; do
              if [[ -d "$site_dir" ]] && PYTHONPATH="$site_dir${{PYTHONPATH:+:$PYTHONPATH}}" "$pypy_bin" -c 'import intervaltree' >/dev/null 2>&1; then
                pypy_py_path="$site_dir"
                break
              fi
            done

            if [[ -n "$pypy_py_path" ]]; then
              export PYTHONPATH="$pypy_py_path${{PYTHONPATH:+:$PYTHONPATH}}"
              echo "[INFO] using pypy via PYTHONPATH bridge: $pypy_py_path" >&2
            else
              pypy_bin="$python_bin"
              echo "[WARNING] pypy intervaltree unavailable; falling back to python for --pypy." >&2
            fi
          fi

          "$python_bin" "$clair_bin" callVarBam \
            --chkpnt_fn {input.model}/model \
            --ref_fn {input.ref} \
            --bed_fn {input.bed} \
            --bam_fn {input.bam} \
            --call_fn {params.out_vcf} \
            --ctgName {params.ctg} \
            --ctgStart {params.ctg_start} \
            --ctgEnd {params.ctg_end} \
            --sampleName {wildcards.sample} \
            --threads {threads} \
            --threshold={params.threshold} \
            --minCoverage={params.min_coverage} \
            --samtools "$samtools_bin" \
            --pypy "$pypy_bin" \
            {params.extra_args}

          if [[ ! -s {params.out_vcf} ]]; then
            echo "Expected Clair2 output VCF missing: {params.out_vcf}" >&2
            find {params.outdir} -maxdepth 3 -type f >&2 || true
            exit 1
          fi

          if command -v bgzip >/dev/null 2>&1; then
            bgzip -f -c {params.out_vcf} > {output.vcf}
          else
            gzip -f -c {params.out_vcf} > {output.vcf}
          fi
        ) > {log} 2>&1
        """


rule compare_xist_to_previous:
    input:
        candidate=rules.call_clair2_xist.output.vcf,
        truth=truth_xist_path,
        script="workflow/scripts/benchmark_xist_match.py",
    output:
        summary=f"{RESULTS_BASE}/{{version}}/{{sample}}.comparison.tsv",
        threshold_scan=f"{RESULTS_BASE}/{{version}}/{{sample}}.threshold_scan.tsv",
        missing_truth=f"{RESULTS_BASE}/{{version}}/{{sample}}.missing_truth.tsv",
    log:
        f"{LOG_BASE}/{{sample}}.{{version}}.compare.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python {input.script} \
            --sample {wildcards.sample} \
            --version {wildcards.version} \
            --candidate-vcf {input.candidate} \
            --truth-vcf {input.truth} \
            --out-summary {output.summary} \
            --out-threshold-scan {output.threshold_scan} \
            --out-missing-truth {output.missing_truth}
        ) > {log} 2>&1
        """


rule rank_versions:
    input:
        summaries=lambda wc: expand(
            f"{RESULTS_BASE}/{{version}}/{{sample}}.comparison.tsv",
            sample=wc.sample,
            version=VERSION_NAMES,
        ),
        script="workflow/scripts/rank_xist_benchmark.py",
    output:
        ranking=f"{RESULTS_BASE}/{{sample}}.ranking.tsv",
    log:
        f"{LOG_BASE}/{{sample}}.ranking.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python {input.script} \
            --summaries {input.summaries} \
            --output {output.ranking}
        ) > {log} 2>&1
        """
