import csv

configfile: "config/config.yaml"


shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")


def load_samples(tsv_path):
    samples = {}
    with open(tsv_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "fastq", "fast5_dir"}
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
SAMPLE_NAMES = sorted(SAMPLES.keys())

LOCI = config["loci"]
LOCUS_NAMES = sorted(LOCI.keys())
HAPS = ["HP1", "HP2"]

THREADS = config["threads"]
CONTAINERS = config["containers"]
REFERENCE_FASTA = config["reference_fasta"]
REFERENCE_URL = config["reference_url"]
CLAIR_MODEL_NAME = config["clair_model"]["name"]
CLAIR_MODEL_ARCHIVE = f"resources/models/clair/{CLAIR_MODEL_NAME}.tar"
CLAIR_MODEL_DIR = f"resources/models/clair/{CLAIR_MODEL_NAME}"
ANNOTATION_GTF = config["annotation"]["gtf"]


wildcard_constraints:
    sample="|".join(SAMPLE_NAMES),
    locus="|".join(LOCUS_NAMES),
    hap="HP1|HP2",
    window="phase_region|dma_region|plot_region",


def sample_field(wildcards, field, default=""):
    value = SAMPLES[wildcards.sample].get(field, "")
    return value if value else default


def window(locus, key):
    """Return (chrom, start, end) for one of a locus' three windows."""
    entry = LOCI[locus]
    span = entry[key]
    return entry["chrom"], int(span["start"]), int(span["end"])


def region_str(locus, key):
    chrom, start, end = window(locus, key)
    return f"{chrom}:{start}-{end}"


def locus_targets(pattern, **kwargs):
    return expand(pattern, sample=SAMPLE_NAMES, locus=LOCUS_NAMES, **kwargs)


rule all:
    input:
        expand("results/{sample}/alignment/{sample}.sorted.bam", sample=SAMPLE_NAMES),
        expand("results/{sample}/alignment/{sample}.sorted.bam.bai", sample=SAMPLE_NAMES),
        locus_targets("results/{sample}/{locus}/variants/{sample}.{locus}.vcf.gz"),
        locus_targets("results/{sample}/{locus}/variants/{sample}.{locus}.vcf.gz.tbi"),
        locus_targets("results/{sample}/{locus}/variants/{sample}.{locus}.phased.vcf.gz"),
        locus_targets("results/{sample}/{locus}/variants/{sample}.{locus}.phased.vcf.gz.tbi"),
        locus_targets("results/{sample}/{locus}/variants/{sample}.{locus}.filt.vcf.gz"),
        locus_targets("results/{sample}/{locus}/variants/{sample}.{locus}.filt.vcf.gz.tbi"),
        locus_targets("results/{sample}/{locus}/qc/{sample}.{locus}.vcf.summary.tsv"),
        locus_targets("results/{sample}/{locus}/qc/{sample}.{locus}.phased.vcf.summary.tsv"),
        locus_targets("results/{sample}/{locus}/qc/{sample}.{locus}.filt.vcf.summary.tsv"),
        locus_targets("results/{sample}/{locus}/nanopolish/{sample}.{locus}.methylation_calls.tsv"),
        locus_targets(
            "results/{sample}/{locus}/nanomethphase/{sample}_{locus}_NanoMethPhase_{hap}_MethylFrequency.bedgraph.gz",
            hap=HAPS,
        ),
        locus_targets("results/{sample}/{locus}/dma/{sample}_{locus}_DMLtest.txt.gz"),
        locus_targets("results/{sample}/{locus}/dma/{sample}_{locus}_callDMR.txt.gz"),
        locus_targets("results/{sample}/{locus}/dma/{sample}_{locus}_DMR_summary.tsv"),
        locus_targets("results/{sample}/{locus}/dma/{sample}_{locus}_DMR_DML_sites.tsv"),
        expand("results/{sample}/figures/{sample}.figure1.pdf", sample=SAMPLE_NAMES),
        expand("results/{sample}/figures/{sample}.figure1.png", sample=SAMPLE_NAMES),


# ---------------------------------------------------------------------------
# Sample-wide steps (locus independent)
# ---------------------------------------------------------------------------


rule stage_fastq:
    input:
        fastq=lambda wc: sample_field(wc, "fastq"),
    output:
        fastq="work/{sample}/{sample}.fastq.gz",
    log:
        "logs/{sample}/stage_fastq.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.fastq})
          ln -sf $(realpath {input.fastq}) {output.fastq}
        ) > {log} 2>&1
        """


rule download_reference:
    output:
        fasta=REFERENCE_FASTA,
    params:
        url=REFERENCE_URL,
    log:
        "logs/reference/download_reference.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.fasta})
          tmp=$(mktemp)
          trap 'rm -f "$tmp"' EXIT

          if command -v curl >/dev/null 2>&1; then
            curl -L --fail --retry 5 --retry-delay 10 -o "$tmp" {params.url}
          else
            wget -O "$tmp" {params.url}
          fi

          if [[ "{params.url}" == *.gz ]]; then
            gunzip -c "$tmp" > {output.fasta}
          else
            cp "$tmp" {output.fasta}
          fi
        ) > {log} 2>&1
        """


rule index_reference:
    input:
        fasta=rules.download_reference.output.fasta,
    output:
        fai=REFERENCE_FASTA + ".fai",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/reference/index_reference.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          samtools faidx {input.fasta}
        ) > {log} 2>&1
        """


rule download_annotation:
    output:
        gtf=ANNOTATION_GTF,
    params:
        url=config["annotation"]["gtf_url"],
    log:
        "logs/reference/download_annotation.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.gtf})
          if command -v curl >/dev/null 2>&1; then
            curl -L --fail --retry 5 --retry-delay 10 -o {output.gtf} {params.url}
          else
            wget -O {output.gtf} {params.url}
          fi
        ) > {log} 2>&1
        """


rule align_minimap2:
    input:
        fastq=rules.stage_fastq.output.fastq,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
    output:
        sam=temp("work/{sample}/alignment/{sample}.sam"),
    params:
        preset=config["minimap2"]["preset"],
    threads: THREADS["align"]
    container:
        CONTAINERS["minimap2"]
    log:
        "logs/{sample}/align_minimap2.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.sam})
          minimap2 -ax {params.preset} -t {threads} {input.ref} {input.fastq} > {output.sam}
        ) > {log} 2>&1
        """


rule sort_bam:
    input:
        sam=rules.align_minimap2.output.sam,
    output:
        bam="results/{sample}/alignment/{sample}.sorted.bam",
    threads: THREADS["samtools_sort"]
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/sort_bam.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.bam})
          samtools sort -@ {threads} -o {output.bam} {input.sam}
        ) > {log} 2>&1
        """


rule index_bam:
    input:
        bam=rules.sort_bam.output.bam,
    output:
        bai="results/{sample}/alignment/{sample}.sorted.bam.bai",
    threads: 4
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/index_bam.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          samtools index -@ {threads} {input.bam} {output.bai}
        ) > {log} 2>&1
        """


rule nanopolish_index:
    input:
        fastq=rules.stage_fastq.output.fastq,
    output:
        done="work/{sample}/nanopolish/{sample}.index.done",
    params:
        fast5_dir=lambda wc: sample_field(wc, "fast5_dir"),
        sequencing_summary=lambda wc: sample_field(wc, "sequencing_summary"),
    container:
        CONTAINERS["nanopolish"]
    log:
        "logs/{sample}/nanopolish_index.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.done})
          SEQ_SUMMARY_OPT=""
          if [[ -n "{params.sequencing_summary}" && -f "{params.sequencing_summary}" ]]; then
            SEQ_SUMMARY_OPT="-s {params.sequencing_summary}"
          fi

          nanopolish index -d {params.fast5_dir} $SEQ_SUMMARY_OPT {input.fastq}
          touch {output.done}
        ) > {log} 2>&1
        """


rule download_clair_model:
    output:
        archive=CLAIR_MODEL_ARCHIVE,
    params:
        url=config["clair_model"]["source_url"],
    log:
        "logs/reference/download_clair_model.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.archive})
          model_url="{params.url}"
          if [[ -z "$model_url" ]]; then
            echo "Missing clair_model.source_url in config/config.yaml" >&2
            exit 1
          fi

          if [[ -f "$model_url" ]]; then
            cp "$model_url" {output.archive}
          elif [[ "$model_url" == file://* ]]; then
            cp "${{model_url#file://}}" {output.archive}
          elif command -v curl >/dev/null 2>&1; then
            curl -L --fail --retry 5 --retry-delay 10 -o {output.archive} "$model_url"
          else
            wget -O {output.archive} "$model_url"
          fi
        ) > {log} 2>&1
        """


rule prepare_clair_model:
    input:
        archive=rules.download_clair_model.output.archive,
    output:
        model_dir=directory(CLAIR_MODEL_DIR),
    container:
        CONTAINERS["clair"]
    log:
        "logs/reference/prepare_clair_model.log",
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
            echo "Prepared Clair model is missing model.index at {output.model_dir}" >&2
            exit 1
          fi
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Per-locus steps
# ---------------------------------------------------------------------------


rule subset_locus_bam:
    input:
        bam=rules.sort_bam.output.bam,
        bai=rules.index_bam.output.bai,
    output:
        bam="results/{sample}/{locus}/alignment/{sample}.{locus}.bam",
        bai="results/{sample}/{locus}/alignment/{sample}.{locus}.bam.bai",
    params:
        region=lambda wc: region_str(wc.locus, "phase_region"),
    threads: 4
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/subset_locus_bam.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.bam})
          samtools view -@ {threads} -b {input.bam} {params.region} -o {output.bam}
          samtools index -@ {threads} {output.bam} {output.bai}
        ) > {log} 2>&1
        """


rule call_clair:
    input:
        bam=rules.sort_bam.output.bam,
        bai=rules.index_bam.output.bai,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
        model=rules.prepare_clair_model.output.model_dir,
    output:
        vcf="results/{sample}/{locus}/variants/{sample}.{locus}.vcf.gz",
    params:
        outdir=lambda wc: f"work/{wc.sample}/{wc.locus}/clair",
        raw_vcf=lambda wc: f"work/{wc.sample}/{wc.locus}/clair/{wc.sample}.vcf",
        ctg=lambda wc: window(wc.locus, "phase_region")[0],
        # Restrict calling to the locus phase_region
        ctg_range_args=lambda wc: "--ctgStart {} --ctgEnd {}".format(
            *window(wc.locus, "phase_region")[1:]
        ),
        threshold=config["clair"].get("threshold", 0.125),
        min_coverage=config["clair"].get("min_coverage", 4),
        extra_args=config["clair"].get("extra_args", ""),
    threads: THREADS["clair"]
    container:
        CONTAINERS["clair"]
    log:
        "logs/{sample}/{locus}/call_clair.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.vcf})
          rm -rf {params.outdir}
          mkdir -p {params.outdir}

          # Clair's pypy subtasks import intervaltree, which the image installs for
          # python3 only, so put python3's site-packages on PYTHONPATH. Resolved at
          # run time rather than hard-coded, so the rule survives an image rebuild.
          export PYTHONPATH="$(python3 -c 'import site; print(site.getsitepackages()[0])')${{PYTHONPATH:+:$PYTHONPATH}}"
          if ! pypy3 -c 'import intervaltree' >/dev/null 2>&1; then
            echo "ERROR: pypy3 cannot import intervaltree; PYTHONPATH=$PYTHONPATH" >&2
            echo "Expected a Clair 2.1.1 image; see containers.clair in config." >&2
            exit 1
          fi

          python3 "$(command -v clair.py)" callVarBam \
            --chkpnt_fn {input.model}/model \
            --ref_fn {input.ref} \
            --bam_fn {input.bam} \
            --call_fn {params.raw_vcf} \
            --ctgName {params.ctg} {params.ctg_range_args} \
            --sampleName {wildcards.sample} \
            --threads {threads} \
            --threshold={params.threshold} \
            --minCoverage={params.min_coverage} \
            --samtools "$(command -v samtools)" \
            --pypy "$(command -v pypy3)" {params.extra_args}

          if [[ ! -s {params.raw_vcf} ]]; then
            echo "Expected Clair output VCF missing: {params.raw_vcf}" >&2
            find {params.outdir} -maxdepth 3 -type f >&2 || true
            exit 1
          fi

          bgzip -f -c {params.raw_vcf} > {output.vcf}
        ) > {log} 2>&1
        """


rule index_called_vcf:
    input:
        vcf=rules.call_clair.output.vcf,
    output:
        tbi="results/{sample}/{locus}/variants/{sample}.{locus}.vcf.gz.tbi",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/index_called_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          tabix -f -p vcf {input.vcf}
        ) > {log} 2>&1
        """


rule filter_snps_for_phasing:
    input:
        vcf=rules.call_clair.output.vcf,
        tbi=rules.index_called_vcf.output.tbi,
    output:
        vcf=temp("results/{sample}/{locus}/variants/{sample}.{locus}.snps.filtered.vcf"),
    params:
        min_q=config["snp_filter"]["min_q"],
        biallelic_snp_only=config["snp_filter"].get("biallelic_snp_only", False),
        biallelic_snp_args=(
            "--remove-indels --min-alleles 2 --max-alleles 2"
            if config["snp_filter"].get("biallelic_snp_only", False)
            else ""
        ),
    container:
        CONTAINERS["vcftools"]
    log:
        "logs/{sample}/{locus}/filter_snps_for_phasing.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          tmp_vcf={output.vcf}.tmp
          vcftools \
            --gzvcf {input.vcf} \
            --minQ {params.min_q} \
            {params.biallelic_snp_args} \
            --recode \
            --recode-INFO-all \
            --stdout \
          > "$tmp_vcf"

          site_count=$(grep -vc '^#' "$tmp_vcf" || true)
          if [[ "${{site_count:-0}}" -eq 0 ]]; then
            echo "No variants passed filter_snps_for_phasing." >&2
            echo "Filter: biallelic_snp_only={params.biallelic_snp_only}, minQ={params.min_q}." >&2
            echo "Lower config:snp_filter:min_q for Clair QUAL scale." >&2
            rm -f "$tmp_vcf"
            exit 1
          fi

          mv "$tmp_vcf" {output.vcf}
        ) > {log} 2>&1
        """


rule phase_variants_whatshap:
    input:
        vcf=rules.filter_snps_for_phasing.output.vcf,
        bam=rules.sort_bam.output.bam,
        bai=rules.index_bam.output.bai,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
    output:
        vcf=temp("results/{sample}/{locus}/variants/{sample}.{locus}.phased.vcf"),
    threads: THREADS["whatshap"]
    container:
        CONTAINERS["whatshap"]
    log:
        "logs/{sample}/{locus}/phase_variants_whatshap.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          whatshap phase \
            --ignore-read-groups \
            --reference {input.ref} \
            --output {output.vcf} \
            {input.vcf} {input.bam}
        ) > {log} 2>&1
        """


rule bgzip_phase_vcf:
    input:
        vcf=rules.phase_variants_whatshap.output.vcf,
    output:
        vcf="results/{sample}/{locus}/variants/{sample}.{locus}.phased.vcf.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/bgzip_phase_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bgzip -c {input.vcf} > {output.vcf}
        ) > {log} 2>&1
        """


rule index_phase_vcf:
    input:
        vcf=rules.bgzip_phase_vcf.output.vcf,
    output:
        tbi="results/{sample}/{locus}/variants/{sample}.{locus}.phased.vcf.gz.tbi",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/index_phase_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          tabix -f -p vcf {input.vcf}
        ) > {log} 2>&1
        """


rule locus_filter_vcf:
    input:
        vcf=rules.bgzip_phase_vcf.output.vcf,
        tbi=rules.index_phase_vcf.output.tbi,
    output:
        vcf=temp("results/{sample}/{locus}/variants/{sample}.{locus}.filt.vcf"),
    params:
        chrom=lambda wc: window(wc.locus, "dma_region")[0],
        start=lambda wc: window(wc.locus, "dma_region")[1],
        end=lambda wc: window(wc.locus, "dma_region")[2],
        min_q=config["snp_filter"]["min_q"],
    container:
        CONTAINERS["vcftools"]
    log:
        "logs/{sample}/{locus}/locus_filter_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          vcftools \
            --gzvcf {input.vcf} \
            --chr {params.chrom} \
            --from-bp {params.start} \
            --to-bp {params.end} \
            --minQ {params.min_q} \
            --recode \
            --recode-INFO-all \
            --stdout \
          > {output.vcf}
        ) > {log} 2>&1
        """


rule bgzip_locus_vcf:
    input:
        vcf=rules.locus_filter_vcf.output.vcf,
    output:
        vcf="results/{sample}/{locus}/variants/{sample}.{locus}.filt.vcf.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/bgzip_locus_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bgzip -c {input.vcf} > {output.vcf}
        ) > {log} 2>&1
        """


rule index_locus_vcf:
    input:
        vcf=rules.bgzip_locus_vcf.output.vcf,
    output:
        tbi="results/{sample}/{locus}/variants/{sample}.{locus}.filt.vcf.gz.tbi",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/index_locus_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          tabix -f -p vcf {input.vcf}
        ) > {log} 2>&1
        """


rule summarize_called_vcf:
    input:
        vcf=rules.call_clair.output.vcf,
        tbi=rules.index_called_vcf.output.tbi,
        script="workflow/scripts/summarize_vcf_metrics.py",
    output:
        summary="results/{sample}/{locus}/qc/{sample}.{locus}.vcf.summary.tsv",
    log:
        "logs/{sample}/{locus}/summarize_called_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.summary})
          python {input.script} --vcf {input.vcf} --out-summary {output.summary}
        ) > {log} 2>&1
        """


rule summarize_phased_vcf:
    input:
        vcf=rules.bgzip_phase_vcf.output.vcf,
        tbi=rules.index_phase_vcf.output.tbi,
        script="workflow/scripts/summarize_vcf_metrics.py",
    output:
        summary="results/{sample}/{locus}/qc/{sample}.{locus}.phased.vcf.summary.tsv",
    log:
        "logs/{sample}/{locus}/summarize_phased_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.summary})
          python {input.script} --vcf {input.vcf} --out-summary {output.summary}
        ) > {log} 2>&1
        """


rule summarize_locus_filtered_vcf:
    input:
        vcf=rules.bgzip_locus_vcf.output.vcf,
        tbi=rules.index_locus_vcf.output.tbi,
        script="workflow/scripts/summarize_vcf_metrics.py",
    output:
        summary="results/{sample}/{locus}/qc/{sample}.{locus}.filt.vcf.summary.tsv",
    log:
        "logs/{sample}/{locus}/summarize_locus_filtered_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.summary})
          python {input.script} --vcf {input.vcf} --out-summary {output.summary}
        ) > {log} 2>&1
        """


rule nanopolish_call_methylation:
    input:
        fastq=rules.stage_fastq.output.fastq,
        bam=rules.subset_locus_bam.output.bam,
        bai=rules.subset_locus_bam.output.bai,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
        indexed=rules.nanopolish_index.output.done,
    output:
        tsv="results/{sample}/{locus}/nanopolish/{sample}.{locus}.methylation_calls.tsv",
    threads: THREADS["nanopolish"]
    container:
        CONTAINERS["nanopolish"]
    log:
        "logs/{sample}/{locus}/nanopolish_call_methylation.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.tsv})
          nanopolish call-methylation \
            -t {threads} \
            -q cpg \
            -r {input.fastq} \
            -b {input.bam} \
            -g {input.ref} \
          > {output.tsv}
        ) > {log} 2>&1
        """


rule methyl_call_processor:
    input:
        calls=rules.nanopolish_call_methylation.output.tsv,
    output:
        bed="results/{sample}/{locus}/nanomethphase/{sample}.{locus}.MethylationCall.bed.gz",
        tbi="results/{sample}/{locus}/nanomethphase/{sample}.{locus}.MethylationCall.bed.gz.tbi",
    threads: THREADS["methyl_call_processor"]
    container:
        CONTAINERS["nanomethphase"]
    log:
        "logs/{sample}/{locus}/methyl_call_processor.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.bed})
          python3 -m nanomethphase methyl_call_processor \
            -mc {input.calls} \
            -t {threads} \
          | sort -k1,1 -k2,2n -k3,3n \
          | bgzip -c > {output.bed}

          tabix -f -p bed {output.bed}
        ) > {log} 2>&1
        """


rule nanomethphase_phase:
    input:
        bam=rules.subset_locus_bam.output.bam,
        bai=rules.subset_locus_bam.output.bai,
        vcf=rules.bgzip_phase_vcf.output.vcf,
        vcf_tbi=rules.index_phase_vcf.output.tbi,
        meth_bed=rules.methyl_call_processor.output.bed,
        meth_bed_tbi=rules.methyl_call_processor.output.tbi,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
    output:
        hp1="results/{sample}/{locus}/nanomethphase/{sample}_{locus}_NanoMethPhase_HP1_MethylFrequency.tsv",
        hp2="results/{sample}/{locus}/nanomethphase/{sample}_{locus}_NanoMethPhase_HP2_MethylFrequency.tsv",
    params:
        prefix=lambda wc: f"results/{wc.sample}/{wc.locus}/nanomethphase/{wc.sample}_{wc.locus}",
        min_snv=config["nanomethphase"].get("min_snv", 2),
    threads: THREADS["nanomethphase_phase"]
    container:
        CONTAINERS["nanomethphase"]
    log:
        "logs/{sample}/{locus}/nanomethphase_phase.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.hp1})

          # --min_SNV is the v1.0 spelling of the read-assignment option, default 2.
          # v1.1.0+ renamed it --min_variant and defaulted it to 1
          python3 -m nanomethphase phase \
            -b {input.bam} \
            -v {input.vcf} \
            -mc {input.meth_bed} \
            -r {input.ref} \
            -o {params.prefix} \
            -of methylcall \
            -ms {params.min_snv} \
            -t {threads} \
            --overwrite

          test -s {output.hp1}
          test -s {output.hp2}
        ) > {log} 2>&1
        """


rule methylfreq_to_bedgraph:
    input:
        tsv="results/{sample}/{locus}/nanomethphase/{sample}_{locus}_NanoMethPhase_{hap}_MethylFrequency.tsv",
    output:
        bedgraph="results/{sample}/{locus}/nanomethphase/{sample}_{locus}_NanoMethPhase_{hap}_MethylFrequency.bedgraph.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/methylfreq_to_bedgraph_{hap}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $1,$2,$3,$7}}' {input.tsv} | bgzip -c > {output.bedgraph}
        ) > {log} 2>&1
        """


rule subset_dma_methylfreq:
    input:
        tsv="results/{sample}/{locus}/nanomethphase/{sample}_{locus}_NanoMethPhase_{hap}_MethylFrequency.tsv",
    output:
        tsv="results/{sample}/{locus}/dma/{sample}_{locus}_{hap}_MethylFrequency.tsv",
    params:
        chrom=lambda wc: window(wc.locus, "dma_region")[0],
        start=lambda wc: window(wc.locus, "dma_region")[1],
        end=lambda wc: window(wc.locus, "dma_region")[2],
    log:
        "logs/{sample}/{locus}/subset_dma_methylfreq_{hap}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.tsv})
          awk -v chrom={params.chrom} -v start={params.start} -v end={params.end} 'BEGIN{{FS=OFS="\t"}} NR==1 {{print; next}} $1==chrom && $2>=start && $3<=end {{print}}' {input.tsv} > {output.tsv}
        ) > {log} 2>&1
        """


rule nanomethphase_dma:
    input:
        hp1="results/{sample}/{locus}/dma/{sample}_{locus}_HP1_MethylFrequency.tsv",
        hp2="results/{sample}/{locus}/dma/{sample}_{locus}_HP2_MethylFrequency.tsv",
    output:
        dml="results/{sample}/{locus}/dma/{sample}_{locus}_DMLtest.txt",
        dmr="results/{sample}/{locus}/dma/{sample}_{locus}_callDMR.txt",
        # dss inputs 
        ready_case="results/{sample}/{locus}/dma/{sample}_{locus}_ReadyForDSS_case1.tsv",
        ready_control="results/{sample}/{locus}/dma/{sample}_{locus}_ReadyForDSS_control1.tsv",
    params:
        outdir=lambda wc: f"results/{wc.sample}/{wc.locus}/dma",
        prefix=lambda wc: f"{wc.sample}_{wc.locus}",
        columns=config["dma"]["columns"],
        rscript=config["dma"].get("rscript", "Rscript"),
    threads: THREADS["dma"]
    container:
        CONTAINERS["nanomethphase"]
    log:
        "logs/{sample}/{locus}/nanomethphase_dma.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p {params.outdir}
          python3 -m nanomethphase dma \
            -c {params.columns} \
            -ca {input.hp1} \
            -co {input.hp2} \
            -o {params.outdir} \
            -op {params.prefix} \
            --Rscript {params.rscript} \
            -ow

          test -s {output.dml}
          test -s {output.dmr}
          test -s {output.ready_case}
          test -s {output.ready_control}
        ) > {log} 2>&1
        """


rule bgzip_dmltest:
    input:
        dml=rules.nanomethphase_dma.output.dml,
    output:
        dml="results/{sample}/{locus}/dma/{sample}_{locus}_DMLtest.txt.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/bgzip_dmltest.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bgzip -c {input.dml} > {output.dml}
        ) > {log} 2>&1
        """


rule bgzip_calldmr:
    input:
        dmr=rules.nanomethphase_dma.output.dmr,
    output:
        dmr="results/{sample}/{locus}/dma/{sample}_{locus}_callDMR.txt.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/{locus}/bgzip_calldmr.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bgzip -c {input.dmr} > {output.dmr}
        ) > {log} 2>&1
        """


rule supplementary_dml_table:
    # report raw values as well as values from DSS 
    input:
        dmr=rules.nanomethphase_dma.output.dmr,
        dml=rules.nanomethphase_dma.output.dml,
        case=rules.nanomethphase_dma.output.ready_case,
        control=rules.nanomethphase_dma.output.ready_control,
        script="workflow/scripts/make_dml_supplementary_table.py",
    output:
        dmr="results/{sample}/{locus}/dma/{sample}_{locus}_DMR_summary.tsv",
        dml="results/{sample}/{locus}/dma/{sample}_{locus}_DMR_DML_sites.tsv",
    params:
        p_threshold=config["dma"].get("pval_cutoff", 0.001),
        fdr_threshold=config["dma"].get("fdr_report_threshold", 0.05),
    log:
        "logs/{sample}/{locus}/supplementary_dml_table.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python {input.script} \
            --locus {wildcards.locus} \
            --dmr {input.dmr} \
            --dmltest {input.dml} \
            --case {input.case} \
            --control {input.control} \
            --p-threshold {params.p_threshold} \
            --fdr-threshold {params.fdr_threshold} \
            --out-dmr {output.dmr} \
            --out-dml {output.dml}
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Figure 1
# ---------------------------------------------------------------------------


rule figure1:
    input:
        methylfreq=lambda wc: expand(
            "results/{sample}/{locus}/nanomethphase/{sample}_{locus}_NanoMethPhase_{hap}_MethylFrequency.tsv",
            sample=wc.sample,
            locus=LOCUS_NAMES,
            hap=HAPS,
        ),
        dmr=lambda wc: expand(
            "results/{sample}/{locus}/dma/{sample}_{locus}_callDMR.txt",
            sample=wc.sample,
            locus=LOCUS_NAMES,
        ),
        gtf=ANNOTATION_GTF,
        config="config/config.yaml",
        script="workflow/scripts/plot_figure1.py",
    output:
        pdf="results/{sample}/figures/{sample}.figure1.pdf",
        png="results/{sample}/figures/{sample}.figure1.png",
        panel_data="results/{sample}/figures/{sample}.figure1.source_data.tsv",
    conda:
        "workflow/envs/figures.yaml"
    log:
        "logs/{sample}/figure1.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.pdf})
          python {input.script} \
            --config {input.config} \
            --sample {wildcards.sample} \
            --gtf {input.gtf} \
            --out-pdf {output.pdf} \
            --out-png {output.png} \
            --out-source-data {output.panel_data}
        ) > {log} 2>&1
        """
