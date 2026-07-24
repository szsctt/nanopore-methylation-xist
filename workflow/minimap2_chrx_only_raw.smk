import csv


configfile: "config/minimap2_chrx_only_raw.yaml"


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

source_reference_fasta = config["source_reference_fasta"]
reference_fasta = config["reference_fasta"]
truth_vcf = config["truth_vcf"]
samtools_container = config["samtools_container"]
clair_container = config["clair_container"]
clair_model_dir = config["clair_model_dir"]
clair_threshold = config.get("clair_threshold", 0.125)
clair_min_coverage = config.get("clair_min_coverage", 4)
clair_extra_args = config.get("clair_extra_args", "")
outdir_template = config.get("outdir", "results/{sample}/benchmarks/minimap2_chrx_only_raw")
work_root = config.get("workdir", "work/benchmark_minimap2_chrx_only_raw")
log_root = config.get("logdir", "logs/benchmark_minimap2_chrx_only_raw")

sam_pattern = work_root + "/{sample}/{mm2ver}/{preset}/{sample}.sam"
bam_pattern = outdir_template + "/{mm2ver}/{preset}/{sample}.sorted.bam"
bai_pattern = outdir_template + "/{mm2ver}/{preset}/{sample}.sorted.bam.bai"
vcf_pattern = outdir_template + "/{mm2ver}/{preset}/{sample}.vcf.gz"
comparison_pattern = outdir_template + "/{mm2ver}/{preset}/{sample}.comparison.tsv"
threshold_pattern = outdir_template + "/{mm2ver}/{preset}/{sample}.threshold_scan.tsv"
missing_pattern = outdir_template + "/{mm2ver}/{preset}/{sample}.missing_truth.tsv"
ranking_pattern = outdir_template + "/{sample}.ranking.tsv"

extract_log = log_root + "/extract_reference.log"
align_log_pattern = log_root + "/{sample}.{mm2ver}.{preset}.align.log"
sort_log_pattern = log_root + "/{sample}.{mm2ver}.{preset}.sort.log"
clair_log_pattern = log_root + "/{sample}.{mm2ver}.{preset}.clair.log"
compare_log_pattern = log_root + "/{sample}.{mm2ver}.{preset}.compare.log"
ranking_log_pattern = log_root + "/{sample}.ranking.log"


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


rule extract_chrx_reference:
    input:
        source=source_reference_fasta,
        script="workflow/scripts/extract_fasta_contig.py",
    output:
        fasta=reference_fasta,
        fai=reference_fasta + ".fai",
    params:
        chrom=chrom,
    log:
        extract_log
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python3 {input.script} \
            --input-fasta {input.source} \
            --contig {params.chrom} \
            --output-fasta {output.fasta} \
            --output-fai {output.fai}
        ) > {log} 2>&1
        """


rule align_minimap2_chrx_only:
    input:
        fastq=lambda wc: sample_to_fastq[wc.sample],
        ref=rules.extract_chrx_reference.output.fasta,
        fai=rules.extract_chrx_reference.output.fai,
        runner="workflow/scripts/run_minimap2_container.sh",
    output:
        sam=temp(sam_pattern),
    params:
        mm2_container=lambda wc: mm2_container[wc.mm2ver],
        bind_path=bind_path,
        preset=lambda wc: wc.preset,
        extra_args=minimap2_extra_args,
    threads: threads_align
    log:
        align_log_pattern
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


rule sort_index_chrx_only_bam:
    input:
        sam=rules.align_minimap2_chrx_only.output.sam,
        runner="workflow/scripts/run_samtools_sort_index_container.sh",
    output:
        bam=bam_pattern,
        bai=bai_pattern,
    params:
        samtools_container=samtools_container,
        bind_path=bind_path,
    threads: threads_sort
    log:
        sort_log_pattern
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


rule call_clair2_chrx_only:
    input:
        bam=rules.sort_index_chrx_only_bam.output.bam,
        bai=rules.sort_index_chrx_only_bam.output.bai,
        ref=rules.extract_chrx_reference.output.fasta,
        fai=rules.extract_chrx_reference.output.fai,
        model=clair_model_dir,
        runner="workflow/scripts/run_clair2_container.sh",
    output:
        vcf=vcf_pattern,
    params:
        outdir=lambda wc: f"{work_root}/{wc.sample}/{wc.mm2ver}/{wc.preset}/clair2",
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


rule compare_chrx_only:
    input:
        candidate=rules.call_clair2_chrx_only.output.vcf,
        truth=truth_vcf,
        script="workflow/scripts/compare_vcf_sets.py",
    output:
        comparison=comparison_pattern,
        threshold=threshold_pattern,
        missing=missing_pattern,
    params:
        chrom=chrom,
        run_label=lambda wc: f"chrx_only__{wc.mm2ver}__{wc.preset}",
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


rule rank_chrx_only:
    input:
        comparisons=lambda wc: expand(
            comparison_pattern,
            sample=wc.sample,
            mm2ver=mm2_versions,
            preset=mm2_presets,
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
