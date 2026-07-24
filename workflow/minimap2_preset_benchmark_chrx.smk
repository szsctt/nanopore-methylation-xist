import csv

configfile: "config/minimap2_preset_benchmark_chrx.yaml"

sample_table = config.get("sample_table", "config/samples.tsv")
sample_to_fastq = {}
with open(sample_table, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("sample"):
            sample_to_fastq[row["sample"]] = row["fastq"]

samples = sorted(sample_to_fastq)
mm2_versions = config["minimap2_versions"]
mm2_presets = config["minimap2_presets"]
mm2_container = config["minimap2_container_by_version"]

for v in mm2_versions:
    if v not in mm2_container:
        raise ValueError(f"Missing minimap2_container_by_version entry for {v}")

chrom = config.get("chrom", "chrX")
threads_align = int(config.get("threads_align", 12))
threads_sort = int(config.get("threads_sort", 8))
threads_clair = int(config.get("threads_clair", 12))
bind_path = config.get("bind_path", "")
minimap2_extra_args = config.get("minimap2_extra_args", "")

reference_fasta = config["reference_fasta"]
truth_vcf = config["truth_vcf"]
samtools_container = config["samtools_container"]
clair_container = config["clair_container"]
clair_model_dir = config["clair_model_dir"]
clair_threshold = config.get("clair_threshold", 0.125)
clair_min_coverage = config.get("clair_min_coverage", 4)
clair_extra_args = config.get("clair_extra_args", "")

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


def combos_for_sample(sample):
    return [(sample, v, p) for v in mm2_versions for p in mm2_presets]


rule all:
    input:
        expand("results/{sample}/benchmarks/minimap2_chrx/{sample}.ranking.tsv", sample=samples),


rule align_minimap2_benchmark:
    input:
        fastq=lambda wc: sample_to_fastq[wc.sample],
        ref=reference_fasta,
        runner="workflow/scripts/run_minimap2_container.sh",
    output:
        sam=temp("work/benchmark_minimap2_chrx/{sample}/{mm2ver}/{preset}/{sample}.sam"),
    params:
        mm2_container=lambda wc: mm2_container[wc.mm2ver],
        bind_path=bind_path,
        preset=lambda wc: wc.preset,
        extra_args=minimap2_extra_args,
    threads: threads_align
    log:
        "logs/benchmark_minimap2_chrx/{sample}.{mm2ver}.{preset}.align.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bash {input.runner} \
            --container-image {params.mm2_container} \
            --ref-fn {input.ref} \
            --fastq-fn {input.fastq} \
            --out-sam {output.sam} \
            --preset {params.preset} \
            --threads {threads} \
            --extra-args '{params.extra_args}' \
            --bind-path {params.bind_path}
        ) > {log} 2>&1
        """


rule sort_index_benchmark_bam:
    input:
        sam=rules.align_minimap2_benchmark.output.sam,
        runner="workflow/scripts/run_samtools_sort_index_container.sh",
    output:
        bam="results/{sample}/benchmarks/minimap2_chrx/{mm2ver}/{preset}/{sample}.sorted.bam",
        bai="results/{sample}/benchmarks/minimap2_chrx/{mm2ver}/{preset}/{sample}.sorted.bam.bai",
    params:
        samtools_container=samtools_container,
        bind_path=bind_path,
    threads: threads_sort
    log:
        "logs/benchmark_minimap2_chrx/{sample}.{mm2ver}.{preset}.sort.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bash {input.runner} \
            --container-image {params.samtools_container} \
            --in-sam {input.sam} \
            --out-bam {output.bam} \
            --out-bai {output.bai} \
            --threads {threads} \
            --bind-path {params.bind_path}
        ) > {log} 2>&1
        """


rule call_clair2_benchmark:
    input:
        bam=rules.sort_index_benchmark_bam.output.bam,
        bai=rules.sort_index_benchmark_bam.output.bai,
        ref=reference_fasta,
        model=clair_model_dir,
        runner="workflow/scripts/run_clair2_container.sh",
    output:
        vcf="results/{sample}/benchmarks/minimap2_chrx/{mm2ver}/{preset}/{sample}.vcf.gz",
    params:
        outdir=lambda wc: f"work/benchmark_minimap2_chrx/{wc.sample}/{wc.mm2ver}/{wc.preset}/clair2",
        clair_container=clair_container,
        bind_path=bind_path,
        chrom=chrom,
        threshold=clair_threshold,
        min_coverage=clair_min_coverage,
        extra_args=clair_extra_args,
    threads: threads_clair
    log:
        "logs/benchmark_minimap2_chrx/{sample}.{mm2ver}.{preset}.clair.log",
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


rule compare_benchmark_vcf:
    input:
        candidate=rules.call_clair2_benchmark.output.vcf,
        truth=truth_vcf,
        script="workflow/scripts/compare_vcf_sets.py",
    output:
        comparison="results/{sample}/benchmarks/minimap2_chrx/{mm2ver}/{preset}/{sample}.comparison.tsv",
        threshold="results/{sample}/benchmarks/minimap2_chrx/{mm2ver}/{preset}/{sample}.threshold_scan.tsv",
        missing="results/{sample}/benchmarks/minimap2_chrx/{mm2ver}/{preset}/{sample}.missing_truth.tsv",
    params:
        chrom=chrom,
        run_label=lambda wc: f"{wc.mm2ver}__{wc.preset}",
    log:
        "logs/benchmark_minimap2_chrx/{sample}.{mm2ver}.{preset}.compare.log",
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


rule rank_minimap2_benchmark:
    input:
        comparisons=lambda wc: expand(
            "results/{sample}/benchmarks/minimap2_chrx/{mm2ver}/{preset}/{sample}.comparison.tsv",
            sample=wc.sample,
            mm2ver=mm2_versions,
            preset=mm2_presets,
        ),
        script="workflow/scripts/rank_clair_benchmark.py",
    output:
        ranking="results/{sample}/benchmarks/minimap2_chrx/{sample}.ranking.tsv",
    log:
        "logs/benchmark_minimap2_chrx/{sample}.ranking.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python3 {input.script} --inputs {input.comparisons} --output {output.ranking}
        ) > {log} 2>&1
        """
