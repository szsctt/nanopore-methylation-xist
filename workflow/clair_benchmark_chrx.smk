import csv

configfile: "config/clair_benchmark_chrx.yaml"

sample_table = config.get("sample_table", "config/samples.tsv")
with open(sample_table, newline="") as fh:
    samples = [row["sample"] for row in csv.DictReader(fh, delimiter="\t") if row.get("sample")]

versions = config["versions"]
container_by_version = config["container_by_version"]
chrom = config.get("chrom", "chrX")
model_path = config["model_path"]
truth_vcf = config["truth_vcf"]
threads = int(config.get("threads", 12))
threshold = config.get("threshold", 0.125)
min_coverage = config.get("min_coverage", 4)
extra_args = config.get("extra_args", "")
bind_path = config.get("bind_path", "")

if not versions:
    raise ValueError("config:versions is empty")

for v in versions:
    if v not in container_by_version:
        raise ValueError(f"Missing config:container_by_version entry for version '{v}'")
    image = str(container_by_version[v]).lower()
    if "clair3" in image or "hkubal/clair3" in image:
        raise ValueError(
            f"Version '{v}' is configured with a Clair3 image ({container_by_version[v]}). "
            "This benchmark workflow is Clair2-only."
        )
    if "biocontainers/clair:" not in image:
        raise ValueError(
            f"Version '{v}' does not point to a Clair2 biocontainers image: {container_by_version[v]}"
        )


rule all:
    input:
        expand("results/{sample}/benchmarks/clair_chrx/{sample}.ranking.tsv", sample=samples),


rule call_clair2_chrx:
    input:
        bam=lambda wc: f"results/{wc.sample}/alignment/{wc.sample}.sorted.bam",
        bai=lambda wc: f"results/{wc.sample}/alignment/{wc.sample}.sorted.bam.bai",
        ref=config["reference_fasta"],
        fai=lambda wc: f"{config['reference_fasta']}.fai",
        model=model_path,
        runner="workflow/scripts/run_clair2_container.sh",
    output:
        vcf="results/{sample}/benchmarks/clair_chrx/{version}/{sample}.vcf.gz",
    params:
        outdir=lambda wc: f"work/benchmark_clair_chrx/{wc.sample}/{wc.version}",
        chrom=chrom,
        threshold=threshold,
        min_coverage=min_coverage,
        extra_args=extra_args,
        bind_path=bind_path,
        container_image=lambda wc: container_by_version[wc.version],
    threads: threads
    log:
        "logs/benchmark_clair_chrx/{sample}.{version}.call.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.vcf})
          rm -rf {params.outdir}
          mkdir -p {params.outdir}

          bash {input.runner} \
            --container-image {params.container_image} \
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


rule compare_clair2_chrx:
    input:
        candidate=rules.call_clair2_chrx.output.vcf,
        truth=truth_vcf,
        script="workflow/scripts/compare_vcf_sets.py",
    output:
        comparison="results/{sample}/benchmarks/clair_chrx/{version}/{sample}.comparison.tsv",
        threshold="results/{sample}/benchmarks/clair_chrx/{version}/{sample}.threshold_scan.tsv",
        missing="results/{sample}/benchmarks/clair_chrx/{version}/{sample}.missing_truth.tsv",
    params:
        chrom=chrom,
    log:
        "logs/benchmark_clair_chrx/{sample}.{version}.compare.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python3 {input.script} \
            --candidate-vcf {input.candidate} \
            --truth-vcf {input.truth} \
            --chrom {params.chrom} \
            --sample {wildcards.sample} \
            --version {wildcards.version} \
            --comparison-tsv {output.comparison} \
            --threshold-tsv {output.threshold} \
            --missing-tsv {output.missing}
        ) > {log} 2>&1
        """


rule rank_clair2_chrx:
    input:
        comparisons=lambda wc: expand(
            "results/{sample}/benchmarks/clair_chrx/{version}/{sample}.comparison.tsv",
            sample=wc.sample,
            version=versions,
        ),
        script="workflow/scripts/rank_clair_benchmark.py",
    output:
        ranking="results/{sample}/benchmarks/clair_chrx/{sample}.ranking.tsv",
    log:
        "logs/benchmark_clair_chrx/{sample}.ranking.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python3 {input.script} --inputs {input.comparisons} --output {output.ranking}
        ) > {log} 2>&1
        """
