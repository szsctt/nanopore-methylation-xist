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

XIST = config["xist_region"]
THREADS = config["threads"]
CONTAINERS = config["containers"]
REFERENCE_FASTA = config["reference_fasta"]
REFERENCE_URL = config["reference_url"]
CLAIR3_MODEL_NAME = config["clair3_model"]["name"]
CLAIR3_MODEL_ARCHIVE = f"resources/models/clair3/{CLAIR3_MODEL_NAME}.tar.gz"
CLAIR3_MODEL_DIR = f"resources/models/clair3/{CLAIR3_MODEL_NAME}"


def sample_field(wildcards, field, default=""):
    value = SAMPLES[wildcards.sample].get(field, "")
    return value if value else default


rule all:
    input:
        expand("results/{sample}/alignment/{sample}.sorted.bam", sample=SAMPLE_NAMES),
        expand("results/{sample}/alignment/{sample}.sorted.bam.bai", sample=SAMPLE_NAMES),
        expand("results/{sample}/variants/{sample}.vcf.gz", sample=SAMPLE_NAMES),
        expand("results/{sample}/variants/{sample}.vcf.gz.tbi", sample=SAMPLE_NAMES),
        expand("results/{sample}/variants/{sample}.phased.vcf.gz", sample=SAMPLE_NAMES),
        expand("results/{sample}/variants/{sample}.phased.vcf.gz.tbi", sample=SAMPLE_NAMES),
        expand("results/{sample}/variants/{sample}.XIST.filt.vcf.gz", sample=SAMPLE_NAMES),
        expand("results/{sample}/variants/{sample}.XIST.filt.vcf.gz.tbi", sample=SAMPLE_NAMES),
        expand("results/{sample}/qc/{sample}.vcf.summary.tsv", sample=SAMPLE_NAMES),
        expand("results/{sample}/qc/{sample}.phased.vcf.summary.tsv", sample=SAMPLE_NAMES),
        expand("results/{sample}/qc/{sample}.XIST.filt.vcf.summary.tsv", sample=SAMPLE_NAMES),
        expand("results/{sample}/nanopolish/{sample}.methylation_calls.tsv", sample=SAMPLE_NAMES),
        expand("results/{sample}/nanomethphase/{sample}.gtf.gz", sample=SAMPLE_NAMES),
        expand(
            "results/{sample}/nanomethphase/{sample}_NanoMethPhase_{hap}_MethylFrequency.bedgraph.gz",
            sample=SAMPLE_NAMES,
            hap=["HP1", "HP2"],
        ),
        expand("results/{sample}/dma/{sample}_XIST_DMLtest.txt.gz", sample=SAMPLE_NAMES),
        expand("results/{sample}/dma/{sample}_XIST_callDMR.txt.gz", sample=SAMPLE_NAMES),


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


rule download_clair3_model:
    output:
        archive=CLAIR3_MODEL_ARCHIVE,
    params:
        url=config["clair3_model"]["source_url"],
    container:
        CONTAINERS["clair3"]
    log:
        "logs/reference/download_clair3_model.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.archive})
          if command -v curl >/dev/null 2>&1; then
            curl -L --fail --retry 5 --retry-delay 10 -o {output.archive} {params.url}
          else
            wget -O {output.archive} {params.url}
          fi
        ) > {log} 2>&1
        """


rule prepare_clair3_model:
    input:
        archive=rules.download_clair3_model.output.archive,
    output:
        model_dir=directory(CLAIR3_MODEL_DIR),
    container:
        CONTAINERS["clair3"]
    log:
        "logs/reference/prepare_clair3_model.log",
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
            echo "Could not locate Clair3 model contents in archive: {input.archive}" >&2
            exit 1
          fi

          rm -rf {output.model_dir}
          mkdir -p {output.model_dir}
          cp -a "$model_candidate"/. {output.model_dir}/
        ) > {log} 2>&1
        """


rule call_clair3:
    input:
        bam=rules.sort_bam.output.bam,
        bai=rules.index_bam.output.bai,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
        model=rules.prepare_clair3_model.output.model_dir,
    output:
        vcf="results/{sample}/variants/{sample}.vcf.gz",
    params:
        outdir=lambda wc: f"work/{wc.sample}/clair3",
        platform=config["clair3"]["platform"],
    threads: THREADS["clair3"]
    container:
        CONTAINERS["clair3"]
    log:
        "logs/{sample}/call_clair3.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.vcf})
          rm -rf {params.outdir}
          mkdir -p {params.outdir}

          run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.ref} \
            --output={params.outdir} \
            --threads={threads} \
            --platform={params.platform} \
            --model_path={input.model}

          cp {params.outdir}/merge_output.vcf.gz {output.vcf}
        ) > {log} 2>&1
        """


rule index_called_vcf:
    input:
        vcf=rules.call_clair3.output.vcf,
    output:
        tbi="results/{sample}/variants/{sample}.vcf.gz.tbi",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/index_called_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          tabix -f -p vcf {input.vcf}
        ) > {log} 2>&1
        """


rule filter_snps_for_phasing:
    input:
        vcf=rules.call_clair3.output.vcf,
        tbi=rules.index_called_vcf.output.tbi,
    output:
        vcf=temp("results/{sample}/variants/{sample}.snps.filtered.vcf"),
    params:
        min_q=config["snp_filter"]["min_q"],
    container:
        CONTAINERS["vcftools"]
    log:
        "logs/{sample}/filter_snps_for_phasing.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          vcftools \
            --gzvcf {input.vcf} \
            --minQ {params.min_q} \
            --recode \
            --recode-INFO-all \
            --stdout \
          > {output.vcf}
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
        vcf=temp("results/{sample}/variants/{sample}.phased.vcf"),
    threads: THREADS["whatshap"]
    container:
        CONTAINERS["whatshap"]
    log:
        "logs/{sample}/phase_variants_whatshap.log",
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
        vcf="results/{sample}/variants/{sample}.phased.vcf.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/bgzip_phase_vcf.log",
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
        tbi="results/{sample}/variants/{sample}.phased.vcf.gz.tbi",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/index_phase_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          tabix -f -p vcf {input.vcf}
        ) > {log} 2>&1
        """


rule xist_filter_vcf:
    input:
        vcf=rules.bgzip_phase_vcf.output.vcf,
        tbi=rules.index_phase_vcf.output.tbi,
    output:
        vcf=temp("results/{sample}/variants/{sample}.XIST.filt.vcf"),
    params:
        chrom=XIST["chrom"],
        start=XIST["start"],
        end=XIST["end"],
        min_q=config["snp_filter"]["min_q"],
    container:
        CONTAINERS["vcftools"]
    log:
        "logs/{sample}/xist_filter_vcf.log",
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


rule bgzip_xist_vcf:
    input:
        vcf=rules.xist_filter_vcf.output.vcf,
    output:
        vcf="results/{sample}/variants/{sample}.XIST.filt.vcf.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/bgzip_xist_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bgzip -c {input.vcf} > {output.vcf}
        ) > {log} 2>&1
        """


rule index_xist_vcf:
    input:
        vcf=rules.bgzip_xist_vcf.output.vcf,
    output:
        tbi="results/{sample}/variants/{sample}.XIST.filt.vcf.gz.tbi",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/index_xist_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          tabix -f -p vcf {input.vcf}
        ) > {log} 2>&1
        """


rule summarize_called_vcf:
    input:
        vcf=rules.call_clair3.output.vcf,
        tbi=rules.index_called_vcf.output.tbi,
        script="workflow/scripts/summarize_vcf_metrics.py",
    output:
        summary="results/{sample}/qc/{sample}.vcf.summary.tsv",
    log:
        "logs/{sample}/summarize_called_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python {input.script} --vcf {input.vcf} --out-summary {output.summary}
        ) > {log} 2>&1
        """


rule summarize_phased_vcf:
    input:
        vcf=rules.bgzip_phase_vcf.output.vcf,
        tbi=rules.index_phase_vcf.output.tbi,
        script="workflow/scripts/summarize_vcf_metrics.py",
    output:
        summary="results/{sample}/qc/{sample}.phased.vcf.summary.tsv",
    log:
        "logs/{sample}/summarize_phased_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python {input.script} --vcf {input.vcf} --out-summary {output.summary}
        ) > {log} 2>&1
        """


rule summarize_xist_filtered_vcf:
    input:
        vcf=rules.bgzip_xist_vcf.output.vcf,
        tbi=rules.index_xist_vcf.output.tbi,
        script="workflow/scripts/summarize_vcf_metrics.py",
    output:
        summary="results/{sample}/qc/{sample}.XIST.filt.vcf.summary.tsv",
    log:
        "logs/{sample}/summarize_xist_filtered_vcf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python {input.script} --vcf {input.vcf} --out-summary {output.summary}
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


rule nanopolish_call_methylation:
    input:
        fastq=rules.stage_fastq.output.fastq,
        bam=rules.sort_bam.output.bam,
        bai=rules.index_bam.output.bai,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
        indexed=rules.nanopolish_index.output.done,
    output:
        tsv="results/{sample}/nanopolish/{sample}.methylation_calls.tsv",
    threads: THREADS["nanopolish"]
    container:
        CONTAINERS["nanopolish"]
    log:
        "logs/{sample}/nanopolish_call_methylation.log",
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
        bed="results/{sample}/nanomethphase/{sample}.MethylationCall.bed.gz",
        tbi="results/{sample}/nanomethphase/{sample}.MethylationCall.bed.gz.tbi",
    threads: THREADS["methyl_call_processor"]
    container:
        CONTAINERS["nanomethphase"]
    log:
        "logs/{sample}/methyl_call_processor.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.bed})
          python nanomethphase.py methyl_call_processor \
            -mc {input.calls} \
            -t {threads} \
          | sort -k1,1 -k2,2n -k3,3n \
          | bgzip -c > {output.bed}

          tabix -f -p bed {output.bed}
        ) > {log} 2>&1
        """


rule nanomethphase_phase:
    input:
        bam=rules.sort_bam.output.bam,
        bai=rules.index_bam.output.bai,
        vcf=rules.bgzip_phase_vcf.output.vcf,
        vcf_tbi=rules.index_phase_vcf.output.tbi,
        meth_bed=rules.methyl_call_processor.output.bed,
        meth_bed_tbi=rules.methyl_call_processor.output.tbi,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
    output:
        hp1="results/{sample}/nanomethphase/{sample}_NanoMethPhase_HP1_MethylFrequency.tsv",
        hp2="results/{sample}/nanomethphase/{sample}_NanoMethPhase_HP2_MethylFrequency.tsv",
        gtf="results/{sample}/nanomethphase/{sample}.gtf",
    params:
        prefix=lambda wc: f"results/{wc.sample}/nanomethphase/{wc.sample}",
    threads: THREADS["nanomethphase_phase"]
    container:
        CONTAINERS["nanomethphase"]
    log:
        "logs/{sample}/nanomethphase_phase.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p $(dirname {output.hp1})
          python nanomethphase.py phase \
            -b {input.bam} \
            -v {input.vcf} \
            -mc {input.meth_bed} \
            -r {input.ref} \
            -o {params.prefix} \
            -of methylcall \
            -t {threads}

          test -s {output.hp1}
          test -s {output.hp2}
          test -s {output.gtf}
        ) > {log} 2>&1
        """


rule bgzip_gtf:
    input:
        gtf=rules.nanomethphase_phase.output.gtf,
    output:
        gtf="results/{sample}/nanomethphase/{sample}.gtf.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/bgzip_gtf.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bgzip -c {input.gtf} > {output.gtf}
        ) > {log} 2>&1
        """


rule methylfreq_to_bedgraph:
    input:
        tsv="results/{sample}/nanomethphase/{sample}_NanoMethPhase_{hap}_MethylFrequency.tsv",
    output:
        bedgraph="results/{sample}/nanomethphase/{sample}_NanoMethPhase_{hap}_MethylFrequency.bedgraph.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/methylfreq_to_bedgraph_{hap}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $1,$2,$3,$7}}' {input.tsv} | bgzip -c > {output.bedgraph}
        ) > {log} 2>&1
        """


rule subset_xist_methylfreq:
    input:
        tsv="results/{sample}/nanomethphase/{sample}_NanoMethPhase_{hap}_MethylFrequency.tsv",
    output:
        tsv="results/{sample}/dma/{sample}_XIST_{hap}_MethylFrequency.tsv",
    params:
        chrom=XIST["chrom"],
        start=XIST["start"],
        end=XIST["end"],
    log:
        "logs/{sample}/subset_xist_methylfreq_{hap}.log",
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
        hp1="results/{sample}/dma/{sample}_XIST_HP1_MethylFrequency.tsv",
        hp2="results/{sample}/dma/{sample}_XIST_HP2_MethylFrequency.tsv",
    output:
        dml="results/{sample}/dma/{sample}_XIST_DMLtest.txt",
        dmr="results/{sample}/dma/{sample}_XIST_callDMR.txt",
    params:
        outdir=lambda wc: f"results/{wc.sample}/dma",
        prefix=lambda wc: f"{wc.sample}_XIST",
        columns=config["dma"]["columns"],
    threads: THREADS["dma"]
    container:
        CONTAINERS["nanomethphase"]
    log:
        "logs/{sample}/nanomethphase_dma.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          mkdir -p {params.outdir}
          python nanomethphase.py dma \
            -c {params.columns} \
            -ca {input.hp1} \
            -co {input.hp2} \
            -o {params.outdir} \
            -op {params.prefix}

          test -s {output.dml}
          test -s {output.dmr}
        ) > {log} 2>&1
        """


rule bgzip_dmltest:
    input:
        dml=rules.nanomethphase_dma.output.dml,
    output:
        dml="results/{sample}/dma/{sample}_XIST_DMLtest.txt.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/bgzip_dmltest.log",
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
        dmr="results/{sample}/dma/{sample}_XIST_callDMR.txt.gz",
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/bgzip_calldmr.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          bgzip -c {input.dmr} > {output.dmr}
        ) > {log} 2>&1
        """
