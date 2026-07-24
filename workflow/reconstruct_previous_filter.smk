configfile: "config/reconstruct_previous_filter.yaml"

sample = config.get("sample", "OTC3")
outdir = config["outdir"]


rule all:
    input:
        f"{outdir}/{sample}.filter_grid.tsv",
        f"{outdir}/{sample}.filter_top.tsv",
        f"{outdir}/{sample}.best_filter.yaml",
        f"{outdir}/{sample}.filtered.vcf.gz",
        f"{outdir}/{sample}.only_raw.tsv",
        f"{outdir}/{sample}.only_target.tsv",
        f"{outdir}/{sample}.boundary_records.tsv",


rule reconstruct_previous_filter:
    input:
        raw=config["raw_vcf"],
        target=config["target_vcf"],
        script="workflow/scripts/reconstruct_vcf_filter.py",
    output:
        grid=f"{outdir}/{sample}.filter_grid.tsv",
        top=f"{outdir}/{sample}.filter_top.tsv",
        best=f"{outdir}/{sample}.best_filter.yaml",
        filtered=f"{outdir}/{sample}.filtered.vcf.gz",
        only_raw=f"{outdir}/{sample}.only_raw.tsv",
        only_target=f"{outdir}/{sample}.only_target.tsv",
        boundary=f"{outdir}/{sample}.boundary_records.tsv",
    params:
        chrom=config.get("chrom", "chrX"),
        top_n=int(config.get("top_n", 100)),
    log:
        f"logs/reconstruct_previous_filter/{sample}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        (
          python3 {input.script} \
            --raw-vcf {input.raw} \
            --target-vcf {input.target} \
            --chrom {params.chrom} \
            --grid-tsv {output.grid} \
            --top-tsv {output.top} \
            --best-yaml {output.best} \
            --filtered-vcf-gz {output.filtered} \
            --only-raw-tsv {output.only_raw} \
            --only-target-tsv {output.only_target} \
            --boundary-tsv {output.boundary} \
            --top-n {params.top_n}
        ) > {log} 2>&1
        """
