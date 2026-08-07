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
ANALYSIS = config.get("analysis_region", {"chrom": XIST["chrom"], "start": "", "end": ""})
THREADS = config["threads"]
CONTAINERS = config["containers"]
REFERENCE_FASTA = config["reference_fasta"]
REFERENCE_URL = config["reference_url"]
CLAIR_MODEL_NAME = config["clair_model"]["name"]
CLAIR_MODEL_ARCHIVE = f"resources/models/clair/{CLAIR_MODEL_NAME}.tar"
CLAIR_MODEL_DIR = f"resources/models/clair/{CLAIR_MODEL_NAME}"


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


rule subset_region_bam:
    # Genome-wide alignment is retained (reads must map against the whole
    # reference); everything downstream operates on this analysis-region subset.
    input:
        bam=rules.sort_bam.output.bam,
        bai=rules.index_bam.output.bai,
    output:
        bam="results/{sample}/alignment/{sample}.region.bam",
        bai="results/{sample}/alignment/{sample}.region.bam.bai",
    params:
        region=(
            f"{ANALYSIS['chrom']}:{ANALYSIS['start']}-{ANALYSIS['end']}"
            if ANALYSIS.get("start", "") != "" and ANALYSIS.get("end", "") != ""
            else ANALYSIS["chrom"]
        ),
    threads: 4
    container:
        CONTAINERS["samtools"]
    log:
        "logs/{sample}/subset_region_bam.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          samtools view -@ {threads} -b {input.bam} {params.region} -o {output.bam}
          samtools index -@ {threads} {output.bam} {output.bai}
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


rule call_clair:
    input:
        bam=rules.sort_bam.output.bam,
        bai=rules.index_bam.output.bai,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
        model=rules.prepare_clair_model.output.model_dir,
    output:
        vcf="results/{sample}/variants/{sample}.vcf.gz",
    params:
        outdir=lambda wc: f"work/{wc.sample}/clair",
        raw_vcf=lambda wc: f"work/{wc.sample}/clair/{wc.sample}.vcf",
        ctg=ANALYSIS["chrom"],
        # Restrict calling to the analysis region
        ctg_range_args=(
            f"--ctgStart {ANALYSIS['start']} --ctgEnd {ANALYSIS['end']}"
            if ANALYSIS.get("start", "") != "" and ANALYSIS.get("end", "") != ""
            else ""
        ),
        threshold=config["clair"].get("threshold", 0.125),
        min_coverage=config["clair"].get("min_coverage", 4),
        extra_args=config["clair"].get("extra_args", ""),
    threads: THREADS["clair"]
    container:
        CONTAINERS["clair"]
    log:
        "logs/{sample}/call_clair.log",
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
            echo "Missing required binary in Clair container." >&2
            echo "clair_bin=$clair_bin" >&2
            echo "python_bin=$python_bin" >&2
            echo "samtools_bin=$samtools_bin" >&2
            echo "pypy_bin=$pypy_bin" >&2
            exit 1
          fi

          # Some Clair2 images include pypy without intervaltree.
          # Try exposing Python3 site-packages to pypy first (no downloads).
          # Fall back to python for pypy-subtasks only if still unavailable.
          if ! "$pypy_bin" -c 'import intervaltree' >/dev/null 2>&1; then
            echo "[INFO] pypy missing intervaltree; trying PYTHONPATH bridge." >&2
            pypy_py_path=""
            for site_dir in /opt/conda/envs/*/lib/python*/site-packages /usr/local/lib/python*/site-packages /opt/conda/lib/python*/site-packages /usr/lib/python*/site-packages; do
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
          else
            echo "[INFO] pypy intervaltree import OK with default environment." >&2
          fi
          echo "[INFO] --pypy binary set to: $pypy_bin" >&2

          "$python_bin" "$clair_bin" callVarBam \
            --chkpnt_fn {input.model}/model \
            --ref_fn {input.ref} \
            --bam_fn {input.bam} \
            --call_fn {params.raw_vcf} \
            --ctgName {params.ctg} {params.ctg_range_args} \
            --sampleName {wildcards.sample} \
            --threads {threads} \
            --threshold={params.threshold} \
            --minCoverage={params.min_coverage} \
            --samtools "$samtools_bin" \
            --pypy "$pypy_bin" {params.extra_args}

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
        vcf=rules.call_clair.output.vcf,
        tbi=rules.index_called_vcf.output.tbi,
    output:
        vcf=temp("results/{sample}/variants/{sample}.snps.filtered.vcf"),
    params:
        min_q=config["snp_filter"]["min_q"],
        biallelic_snp_only=config["snp_filter"].get("biallelic_snp_only", False),
        biallelic_snp_args=(
            "--remove-indels --min-alleles 2 --max-alleles 2"
            if config["snp_filter"].get("biallelic_snp_only", False)
            else ""
        ),
        # No explicit region filter needed: call_clair already restricts calling
        # to the analysis_region, so the phasing input is region-scoped, which is
        # what keeps distant extra calls out of whatshap's read selection.
    container:
        CONTAINERS["vcftools"]
    log:
        "logs/{sample}/filter_snps_for_phasing.log",
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
        vcf=rules.call_clair.output.vcf,
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
        bam=rules.subset_region_bam.output.bam,
        bai=rules.subset_region_bam.output.bai,
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
        bam=rules.subset_region_bam.output.bam,
        bai=rules.subset_region_bam.output.bai,
        vcf=rules.bgzip_phase_vcf.output.vcf,
        vcf_tbi=rules.index_phase_vcf.output.tbi,
        meth_bed=rules.methyl_call_processor.output.bed,
        meth_bed_tbi=rules.methyl_call_processor.output.tbi,
        ref=rules.download_reference.output.fasta,
        fai=rules.index_reference.output.fai,
    output:
        hp1="results/{sample}/nanomethphase/{sample}_NanoMethPhase_HP1_MethylFrequency.tsv",
        hp2="results/{sample}/nanomethphase/{sample}_NanoMethPhase_HP2_MethylFrequency.tsv",
    params:
        prefix=lambda wc: f"results/{wc.sample}/nanomethphase/{wc.sample}",
        min_snv=config["nanomethphase"].get("min_snv", 2),
        expected_version=config["nanomethphase"].get("expected_version", "1.0"),
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

          # Guard the pinned build. v1.0 has no --version, so verify the v1.0
          # read-assignment option instead: v1.0 exposes --min_SNV, which v1.1.0+
          # renamed to --min_variant. See config nanomethphase / reproduction_progress.md.
          help_out=$(python3 -m nanomethphase phase --help 2>&1 || true)
          if [[ "$help_out" != *"--min_SNV"* ]]; then
            echo "ERROR: nanomethphase build lacks --min_SNV; expected v{params.expected_version}." >&2
            echo "v1.1.0+ renamed it --min_variant. Pin containers:nanomethphase to a v1.0 build." >&2
            exit 1
          fi
          echo "[INFO] nanomethphase v{params.expected_version} confirmed (--min_SNV present)" >&2

          # v1.0 uses --min_SNV (default 2), the historical read-assignment.
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
        dss_version=config["dma"].get("dss_version", "2.58.0"),
        rscript=config["dma"].get("rscript", "Rscript"),
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
          {params.rscript} -e "v <- as.character(packageVersion('DSS')); cat('DSS version:', v, '\\n'); if (!identical(v, '{params.dss_version}')) stop(paste('Expected DSS', '{params.dss_version}', 'got', v))"
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
