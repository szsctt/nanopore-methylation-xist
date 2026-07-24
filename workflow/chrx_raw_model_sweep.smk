import csv
from pathlib import Path


configfile: "config/chrx_raw_model_sweep.yaml"


sample_table = config.get("sample_table", "config/samples.tsv")
with open(sample_table, newline="") as fh:
    samples = [row["sample"] for row in csv.DictReader(fh, delimiter="\t") if row.get("sample")]

models = config["models"]
model_labels = list(models.keys())
chrom = config.get("chrom", "chrX")
truth_vcf = config["truth_vcf"]
reference_fasta = config["reference_fasta"]
threads = int(config.get("threads", 12))
bind_path = config.get("bind_path", "")
outdir_template = config.get("outdir", "results/{sample}/benchmarks/chrx_raw_model_sweep")
work_root = config.get("workdir", "work/benchmark_chrx_raw_model_sweep")
log_root = config.get("logdir", "logs/benchmark_chrx_raw_model_sweep")
vcf_pattern = outdir_template + "/{model}/{sample}.vcf.gz"
comparison_pattern = outdir_template + "/{model}/{sample}.comparison.tsv"
threshold_pattern = outdir_template + "/{model}/{sample}.threshold_scan.tsv"
missing_pattern = outdir_template + "/{model}/{sample}.missing_truth.tsv"
ranking_pattern = outdir_template + "/{sample}.ranking.tsv"
call_log_pattern = log_root + "/{sample}.{model}.call.log"
compare_log_pattern = log_root + "/{sample}.{model}.compare.log"
ranking_log_pattern = log_root + "/{sample}.ranking.log"


def outdir_for_sample(sample):
    return outdir_template.format(sample=sample)


for label, spec in models.items():
    caller = spec.get("caller")
    if caller not in {"clair2", "clair3"}:
        raise ValueError(f"Model '{label}' has unsupported caller: {caller}")

    model_path = Path(spec["model_path"])
    if caller == "clair2" and not (model_path / "model.index").exists():
        raise ValueError(f"Clair2 model '{label}' is missing model.index under {model_path}")
    if caller == "clair3" and not ((model_path / "pileup.index").exists() or (model_path / "pileup.pt").exists()):
        raise ValueError(f"Clair3 model '{label}' is missing pileup model files under {model_path}")


rule all:
    input:
        expand(ranking_pattern, sample=samples),


rule call_clair2_raw_model:
    input:
        bam=lambda wc: config["bam_template"].format(sample=wc.sample),
        bai=lambda wc: config["bam_template"].format(sample=wc.sample) + ".bai",
        ref=reference_fasta,
        fai=f"{reference_fasta}.fai",
        runner="workflow/scripts/run_clair2_container.sh",
    output:
        vcf=vcf_pattern,
    params:
        caller=lambda wc: models[wc.model]["caller"],
        model_dir=lambda wc: models[wc.model]["model_path"],
        container_image=lambda wc: models[wc.model]["container"],
        outdir=lambda wc: f"{work_root}/{wc.sample}/{wc.model}",
        chrom=chrom,
        threshold=config.get("clair2_threshold", 0.125),
        min_coverage=config.get("clair2_min_coverage", 4),
        extra_args=config.get("clair2_extra_args", ""),
        bind_path=bind_path,
    threads: threads
    log:
        call_log_pattern
    wildcard_constraints:
        model="|".join([label for label, spec in models.items() if spec["caller"] == "clair2"])
    shell:
        """
        mkdir -p $(dirname {log})
        (
          if [[ "{params.caller}" != "clair2" ]]; then
            echo "Wrong caller for Clair2 rule: {params.caller}" >&2
            exit 2
          fi

          mkdir -p $(dirname {output.vcf})
          rm -rf {params.outdir}
          mkdir -p {params.outdir}

          bash {input.runner} \
            --container-image {params.container_image} \
            --bam-fn {input.bam} \
            --ref-fn {input.ref} \
            --model-dir {params.model_dir} \
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


rule call_clair3_raw_model:
    input:
        bam=lambda wc: config["bam_template"].format(sample=wc.sample),
        bai=lambda wc: config["bam_template"].format(sample=wc.sample) + ".bai",
        ref=reference_fasta,
        fai=f"{reference_fasta}.fai",
        runner="workflow/scripts/run_clair3_container.sh",
    output:
        vcf=vcf_pattern,
    params:
        caller=lambda wc: models[wc.model]["caller"],
        model_path=lambda wc: models[wc.model]["model_path"],
        container_image=lambda wc: models[wc.model]["container"],
        outdir=lambda wc: f"{work_root}/{wc.sample}/{wc.model}",
        chrom=chrom,
        platform=config.get("clair3_platform", "ont"),
        bind_path=bind_path,
    threads: threads
    log:
        call_log_pattern
    wildcard_constraints:
        model="|".join([label for label, spec in models.items() if spec["caller"] == "clair3"])
    shell:
        """
        mkdir -p $(dirname {log})
        (
          if [[ "{params.caller}" != "clair3" ]]; then
            echo "Wrong caller for Clair3 rule: {params.caller}" >&2
            exit 2
          fi

          mkdir -p {params.outdir}

          if [[ -s {params.outdir}/merge_output.vcf.gz ]]; then
            echo "Reusing existing Clair3 output: {params.outdir}/merge_output.vcf.gz"
          else
            rm -rf {params.outdir}
            mkdir -p {params.outdir}

            bash {input.runner} \
              --container-image {params.container_image} \
              --bam-fn {input.bam} \
              --ref-fn {input.ref} \
              --model-path {params.model_path} \
              --outdir {params.outdir} \
              --ctg-name {params.chrom} \
              --threads {threads} \
              --platform {params.platform} \
              --bind-path {params.bind_path}
          fi

          if [[ ! -s {params.outdir}/merge_output.vcf.gz ]]; then
            echo "Expected Clair3 output missing: {params.outdir}/merge_output.vcf.gz" >&2
            find {params.outdir} -maxdepth 3 -type f >&2 || true
            exit 1
          fi
          mkdir -p $(dirname {output.vcf})
          cp -f {params.outdir}/merge_output.vcf.gz {output.vcf}
        ) > {log} 2>&1
        """


rule compare_raw_model:
    input:
        candidate=vcf_pattern,
        truth=truth_vcf,
        script="workflow/scripts/compare_vcf_sets.py",
    output:
        comparison=comparison_pattern,
        threshold=threshold_pattern,
        missing=missing_pattern,
    params:
        chrom=chrom,
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
            --version {wildcards.model} \
            --comparison-tsv {output.comparison} \
            --threshold-tsv {output.threshold} \
            --missing-tsv {output.missing}
        ) > {log} 2>&1
        """


rule rank_raw_models:
    input:
        comparisons=lambda wc: expand(
            comparison_pattern,
            sample=wc.sample,
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
