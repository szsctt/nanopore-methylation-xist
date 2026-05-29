import csv

configfile: "config/clair_benchmark.yaml"

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
BENCHMARK_NAME = config.get("benchmark_name", "clair_xist")
RESULTS_BASE = f"results/{{sample}}/benchmarks/{BENCHMARK_NAME}"
WORK_BASE = f"work/benchmark_{BENCHMARK_NAME}"
LOG_BASE = f"logs/benchmark_{BENCHMARK_NAME}"
MODEL_NAME = config["clair3_model"]["name"]
MODEL_ARCHIVE = f"resources/models/{BENCHMARK_NAME}/{MODEL_NAME}.tar.gz"
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
        bed=f"{WORK_BASE}/regions/xist.bed",
    params:
        chrom=XIST["chrom"],
        start=XIST["start"],
        end=XIST["end"],
    shell:
        """
        mkdir -p $(dirname {output.bed})
        # BED is 0-based half-open; config stores 1-based inclusive coordinates.
        printf "%s\t%d\t%d\n" {params.chrom} $(({params.start}-1)) {params.end} > {output.bed}
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


rule download_clair3_model:
    output:
        archive=MODEL_ARCHIVE,
    params:
        url=config["clair3_model"]["source_url"],
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
            echo "Missing clair3_model.source_url in config/clair_benchmark.yaml" >&2
            exit 1
          fi

          if [[ "$model_url" == *.tar.gz ]]; then
            if command -v curl >/dev/null 2>&1; then
              curl -L --fail --retry 5 --retry-delay 10 -o {output.archive} "$model_url"
            else
              wget -O {output.archive} "$model_url"
            fi
          else
            # Clair3 v2 PyTorch endpoints are model directories with .pt files.
            tmpdir=$(mktemp -d)
            trap 'rm -rf "$tmpdir"' EXIT

            if [[ "$model_url" != */ ]]; then
              model_url="${{model_url}}/"
            fi

            for model_file in pileup.pt full_alignment.pt; do
              if command -v curl >/dev/null 2>&1; then
                curl -L --fail --retry 5 --retry-delay 10 -o "$tmpdir/$model_file" "${{model_url}}${{model_file}}"
              else
                wget -O "$tmpdir/$model_file" "${{model_url}}${{model_file}}"
              fi
            done

            tar -czf {output.archive} -C "$tmpdir" .
          fi
        ) > {log} 2>&1
        """


rule prepare_clair3_model:
    input:
        archive=rules.download_clair3_model.output.archive,
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

          tar -xzf {input.archive} -C "$tmpdir"
          model_candidate_file=$(find "$tmpdir" -type f \\( -name '*.pt' -o -name '*.index' -o -name 'pileup*' -o -name 'full_alignment*' \\) | head -n 1)
          if [[ -n "$model_candidate_file" ]]; then
            model_candidate=$(dirname "$model_candidate_file")
          else
            model_candidate=$(find "$tmpdir" -mindepth 1 -maxdepth 1 -type d | head -n 1)
          fi

          if [[ -z "$model_candidate" ]]; then
            echo "Could not locate model contents in archive: {input.archive}" >&2
            exit 1
          fi

          rm -rf {output.model_dir}
          mkdir -p {output.model_dir}
          cp -a "$model_candidate"/. {output.model_dir}/
        ) > {log} 2>&1
        """


rule call_clair3_xist:
    input:
        bam=lambda wc: f"results/{wc.sample}/alignment/{wc.sample}.sorted.bam",
        bai=lambda wc: f"results/{wc.sample}/alignment/{wc.sample}.sorted.bam.bai",
        ref=config["reference_fasta"],
        fai=rules.index_reference.output.fai,
        model=rules.prepare_clair3_model.output.model_dir,
        bed=rules.create_xist_bed.output.bed,
    output:
        vcf=f"{RESULTS_BASE}/{{version}}/{{sample}}.vcf.gz",
    params:
        outdir=lambda wc: f"{WORK_BASE}/{wc.sample}/{wc.version}",
        platform=config["clair3"]["platform"],
        ctg=XIST["chrom"],
        container_image=lambda wc: f"docker://docker.io/hkubal/clair3:{wc.version}",
        bind_path=config.get("bind_path", ""),
        extra_args=config["clair3"].get("extra_args", ""),
    threads: THREADS["clair3"]
    log:
        f"{LOG_BASE}/{{sample}}.{{version}}.call.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.vcf})
          rm -rf {params.outdir}
          mkdir -p {params.outdir}

          runtime=""
          if command -v apptainer >/dev/null 2>&1; then
            runtime="apptainer"
          elif command -v singularity >/dev/null 2>&1; then
            runtime="singularity"
          else
            echo "Neither apptainer nor singularity is available in PATH" >&2
            exit 1
          fi

          bind_opt=""
          if [[ -n "{params.bind_path}" ]]; then
            bind_opt="-B {params.bind_path}"
          fi

          $runtime exec --cleanenv $bind_opt {params.container_image} /bin/bash -lc "
            runner=''
            runner_mode='exec'
            if command -v run_clair3.sh >/dev/null 2>&1; then
              runner=\$(command -v run_clair3.sh)
            elif command -v run_clair3 >/dev/null 2>&1; then
              runner=\$(command -v run_clair3)
            else
              for p in /opt/bin/run_clair3.sh /opt/bin/scripts/run_clair3.sh /usr/local/bin/run_clair3.sh /opt/conda/bin/run_clair3.sh; do
                if [[ -x \"\$p\" ]]; then
                  runner=\"\$p\"
                  break
                fi
              done
              if [[ -z \"\$runner\" ]]; then
                for p in /opt/bin/run_clair3.py /opt/bin/scripts/run_clair3.py; do
                  if [[ -f \"\$p\" ]]; then
                    runner=\"\$p\"
                    runner_mode='python'
                    break
                  fi
                done
              fi
            fi

            if [[ -z \"\$runner\" ]]; then
              echo 'Could not locate Clair3 launcher in container.' >&2
              echo 'PATH:' \"\$PATH\" >&2
              ls -la /opt/bin >&2 || true
              exit 1
            fi

            pypy_bin='pypy3'
            if ! command -v \"\$pypy_bin\" >/dev/null 2>&1; then
              pypy_bin='python3'
            fi

            pick_bin() {{
              local default_name=\"\$1\"
              shift
              for candidate in \"\$@\"; do
                if [[ -x \"\$candidate\" ]]; then
                  echo \"\$candidate\"
                  return 0
                fi
              done
              if command -v \"\$default_name\" >/dev/null 2>&1; then
                command -v \"\$default_name\"
                return 0
              fi
              return 1
            }}

            samtools_bin=\$(pick_bin samtools /opt/conda/envs/clair3/bin/samtools /opt/conda/envs/clair3_v2/bin/samtools /opt/conda/envs/clair3_v1/bin/samtools /opt/conda/bin/samtools /usr/local/bin/samtools /usr/bin/samtools /bin/samtools /opt/bin/samtools || true)
            parallel_bin=\$(pick_bin parallel /opt/conda/envs/clair3/bin/parallel /opt/conda/envs/clair3_v2/bin/parallel /opt/conda/envs/clair3_v1/bin/parallel /opt/conda/bin/parallel /usr/local/bin/parallel /usr/bin/parallel /bin/parallel /opt/bin/parallel || true)
            whatshap_bin=\$(pick_bin whatshap /opt/conda/envs/clair3/bin/whatshap /opt/conda/envs/clair3_v2/bin/whatshap /opt/conda/envs/clair3_v1/bin/whatshap /opt/conda/bin/whatshap /usr/local/bin/whatshap /usr/bin/whatshap /bin/whatshap /opt/bin/whatshap || true)
            python_bin=\$(pick_bin python3 /opt/conda/envs/clair3/bin/python /opt/conda/envs/clair3/bin/python3 /opt/conda/envs/clair3_v2/bin/python /opt/conda/envs/clair3_v2/bin/python3 /opt/conda/envs/clair3_v1/bin/python /opt/conda/envs/clair3_v1/bin/python3 /opt/conda/bin/python /opt/conda/bin/python3 /usr/local/bin/python3 /usr/bin/python3 /bin/python3 /opt/bin/python3 || true)
            if [[ "\$pypy_bin" == 'python3' ]]; then
              for candidate in /opt/conda/envs/clair3/bin/pypy3 /opt/conda/envs/clair3_v2/bin/pypy3 /opt/conda/envs/clair3_v1/bin/pypy3 /opt/conda/bin/pypy3 /usr/local/bin/pypy3 /usr/bin/pypy3 /bin/pypy3 /opt/bin/pypy3; do
                if [[ -x "\$candidate" ]]; then
                  pypy_bin="\$candidate"
                  break
                fi
              done
            fi

            # Old Clair3 releases require TensorFlow; ensure chosen python has it.
            if [[ -n "\$python_bin" ]] && ! "\$python_bin" -c 'import tensorflow' >/dev/null 2>&1; then
              for candidate in /opt/conda/envs/clair3/bin/python /opt/conda/envs/clair3/bin/python3 /opt/conda/envs/clair3_v1/bin/python /opt/conda/envs/clair3_v1/bin/python3 /opt/conda/envs/clair3_v2/bin/python /opt/conda/envs/clair3_v2/bin/python3; do
                if [[ -x "\$candidate" ]] && "\$candidate" -c 'import tensorflow' >/dev/null 2>&1; then
                  python_bin="\$candidate"
                  break
                fi
              done
            fi

            if [[ -z \"\$python_bin\" || -z \"\$samtools_bin\" || -z \"\$parallel_bin\" || -z \"\$whatshap_bin\" ]]; then
              echo 'Missing required binaries inside container.' >&2
              echo \"python_bin=\$python_bin\" >&2
              echo \"samtools_bin=\$samtools_bin\" >&2
              echo \"parallel_bin=\$parallel_bin\" >&2
              echo \"whatshap_bin=\$whatshap_bin\" >&2
              echo 'PATH:' \"\$PATH\" >&2
              exit 1
            fi
            if ! "\$python_bin" -c 'import tensorflow' >/dev/null 2>&1; then
              echo 'Selected python lacks tensorflow, Clair3 v0/v1 cannot run.' >&2
              echo \"python_bin=\$python_bin\" >&2
              exit 1
            fi

            if [[ \"\$runner_mode\" == 'python' ]]; then
              \"\$python_bin\" \"\$runner\" \
                --bam_fn={input.bam} \
                --ref_fn={input.ref} \
                --output={params.outdir} \
                --threads={threads} \
                --platform={params.platform} \
                --model_path={input.model} \
                --ctg_name={params.ctg} \
                --bed_fn={input.bed} \
                --python=\"\$python_bin\" \
                --samtools=\"\$samtools_bin\" \
                --parallel=\"\$parallel_bin\" \
                --whatshap=\"\$whatshap_bin\" \
                --pypy=\"\$pypy_bin\" {params.extra_args}
            else
              \"\$runner\" \
                --bam_fn={input.bam} \
                --ref_fn={input.ref} \
                --output={params.outdir} \
                --threads={threads} \
                --platform={params.platform} \
                --model_path={input.model} \
                --ctg_name={params.ctg} \
                --bed_fn={input.bed} \
                --python=\"\$python_bin\" \
                --samtools=\"\$samtools_bin\" \
                --parallel=\"\$parallel_bin\" \
                --whatshap=\"\$whatshap_bin\" \
                --pypy=\"\$pypy_bin\" {params.extra_args}
            fi
          "

          result_vcf=""
          if [[ -s {params.outdir}/merge_output.vcf.gz ]]; then
            result_vcf="{params.outdir}/merge_output.vcf.gz"
          elif [[ -s {params.outdir}/pileup.vcf.gz ]]; then
            result_vcf="{params.outdir}/pileup.vcf.gz"
          elif [[ -s {params.outdir}/pileup.vcf ]]; then
            if command -v bgzip >/dev/null 2>&1; then
              bgzip -f -c {params.outdir}/pileup.vcf > {params.outdir}/merge_output.vcf.gz
            else
              gzip -f -c {params.outdir}/pileup.vcf > {params.outdir}/merge_output.vcf.gz
            fi
            result_vcf="{params.outdir}/merge_output.vcf.gz"
          fi

          if [[ -z "$result_vcf" ]]; then
            echo "No final Clair3 VCF found under {params.outdir}" >&2
            find {params.outdir} -maxdepth 4 -type f >&2 || true
            exit 1
          fi

          cp "$result_vcf" {output.vcf}
        ) > {log} 2>&1
        """


rule compare_xist_to_previous:
    input:
        candidate=rules.call_clair3_xist.output.vcf,
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
