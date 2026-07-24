import csv
from pathlib import Path


configfile: "config/clair2_model_minimap_raw.yaml"


sample_table = config.get("sample_table", "config/samples.tsv")
with open(sample_table, newline="") as fh:
    samples = [row["sample"] for row in csv.DictReader(fh, delimiter="\t") if row.get("sample")]

mm2_versions = config["minimap2_versions"]
models = config["models"]
model_labels = list(models.keys())

chrom = config.get("chrom", "chrX")
reference_fasta = config["reference_fasta"]
truth_vcf = config["truth_vcf"]
threads_clair = int(config.get("threads_clair", 12))
bind_path = config.get("bind_path", "")
clair_container = config["clair_container"]
clair_threshold = config.get("clair_threshold", 0.125)
clair_min_coverage = config.get("clair_min_coverage", 4)
clair_extra_args = config.get("clair_extra_args", "")
prepared_model_dir = config.get("prepared_model_dir", "resources/models/clair_prepared")
outdir_template = config.get("outdir", "results/{sample}/benchmarks/clair2_model_minimap_raw")
work_root = config.get("workdir", "work/benchmark_clair2_model_minimap_raw")
log_root = config.get("logdir", "logs/benchmark_clair2_model_minimap_raw")

prepared_pattern = prepared_model_dir + "/{model}"
vcf_pattern = outdir_template + "/{mm2ver}/{model}/{sample}.vcf.gz"
comparison_pattern = outdir_template + "/{mm2ver}/{model}/{sample}.comparison.tsv"
threshold_pattern = outdir_template + "/{mm2ver}/{model}/{sample}.threshold_scan.tsv"
missing_pattern = outdir_template + "/{mm2ver}/{model}/{sample}.missing_truth.tsv"
ranking_pattern = outdir_template + "/{sample}.ranking.tsv"

prepare_log_pattern = log_root + "/prepare.{model}.log"
clair_log_pattern = log_root + "/{sample}.{mm2ver}.{model}.clair.log"
compare_log_pattern = log_root + "/{sample}.{mm2ver}.{model}.compare.log"
ranking_log_pattern = log_root + "/{sample}.ranking.log"


for label, archive in models.items():
    if not Path(archive).exists():
        raise ValueError(f"Missing Clair2 model archive for {label}: {archive}")

clair_image = str(clair_container).lower()
if "clair3" in clair_image or "hkubal/clair3" in clair_image:
    raise ValueError(
        f"clair_container is set to a Clair3 image ({clair_container}). "
        "This benchmark workflow is Clair2-only."
    )
if "biocontainers/clair:" not in clair_image:
    raise ValueError(
        f"clair_container must point to a Clair2 biocontainers image, got: {clair_container}"
    )


rule all:
    input:
        expand(ranking_pattern, sample=samples),


rule prepare_clair2_model_archive:
    input:
        archive=lambda wc: models[wc.model],
    output:
        model_dir=directory(prepared_pattern),
    log:
        prepare_log_pattern
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
            echo "Prepared Clair2 model is missing model.index: {output.model_dir}" >&2
            exit 1
          fi
        ) > {log} 2>&1
        """


rule call_clair2_model_minimap_raw:
    input:
        bam=lambda wc: config["bam_template"].format(sample=wc.sample, mm2ver=wc.mm2ver),
        bai=lambda wc: config["bam_template"].format(sample=wc.sample, mm2ver=wc.mm2ver) + ".bai",
        ref=reference_fasta,
        model=rules.prepare_clair2_model_archive.output.model_dir,
        runner="workflow/scripts/run_clair2_container.sh",
    output:
        vcf=vcf_pattern,
    params:
        outdir=lambda wc: f"{work_root}/{wc.sample}/{wc.mm2ver}/{wc.model}/clair2",
        clair_container=clair_container,
        bind_path=bind_path,
        chrom=chrom,
        threshold=clair_threshold,
        min_coverage=clair_min_coverage,
        extra_args=clair_extra_args,
    threads: threads_clair
    log:
        clair_log_pattern
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bash {input.runner} \
            --container-image {params.clair_container} \
            --bam-fn {input.bam} \
            --ref-fn {input.ref} \
            --model-dir {input.model} \
            --out-vcfgz {output.vcf} \
            --outdir {params.outdir} \
            --ctg-name {params.chrom} \
            --sample-name {wildcards.sample} \
            --threads {threads} \
            --threshold {params.threshold} \
            --min-coverage {params.min_coverage} \
            --extra-args '{params.extra_args}' \
            --bind-path {params.bind_path}
        ) > {log} 2>&1
        """


rule compare_clair2_model_minimap_raw:
    input:
        candidate=rules.call_clair2_model_minimap_raw.output.vcf,
        truth=truth_vcf,
        script="workflow/scripts/compare_vcf_sets.py",
    output:
        comparison=comparison_pattern,
        threshold=threshold_pattern,
        missing=missing_pattern,
    params:
        chrom=chrom,
        run_label=lambda wc: f"{wc.mm2ver}__{wc.model}",
    log:
        compare_log_pattern
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python3 {input.script} \
            --candidate-vcf {input.candidate} \
            --truth-vcf {input.truth} \
            --chrom {params.chrom} \
            --sample {wildcards.sample} \
            --version {params.run_label} \
            --comparison-tsv {output.comparison} \
            --threshold-tsv {output.threshold} \
            --missing-tsv {output.missing}
        ) > {log} 2>&1
        """


rule rank_clair2_model_minimap_raw:
    input:
        comparisons=lambda wc: expand(
            comparison_pattern,
            sample=wc.sample,
            mm2ver=mm2_versions,
            model=model_labels,
        ),
        script="workflow/scripts/rank_clair_benchmark.py",
    output:
        ranking=ranking_pattern,
    log:
        ranking_log_pattern
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python3 {input.script} --sort-by any --inputs {input.comparisons} --output {output.ranking}
        ) > {log} 2>&1
        """
