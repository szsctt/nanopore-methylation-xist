#!/usr/bin/env python3
"""Plot Figure 1B and 1C: allele-specific 5mC methylation at XIST and OTC.

One panel per locus. Each panel stacks the two NanoMethPhase haplotype
MethylFrequency tracks over a GENCODE gene model, with the DMR shaded.

NanoMethPhase labels haplotypes arbitrarily within each WhatsHap phase block, so
a locus whose block is labelled opposite to the published run can be drawn in
swapped track order via `display_swap_haplotypes` in the locus config. 
"""

import argparse
import csv
import gzip
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import yaml
from matplotlib.patches import Rectangle

HAP_COLOUR = {"HP1": "#1f5fa9", "HP2": "#c1392b"}
DMR_FACE = "#d9d9d9"
DMR_EDGE = "#2b2b2b"
PANEL_LABELS = ["B", "C", "D", "E", "F"]


def read_methylfreq(path, chrom, start, end):
    """Return sorted [(pos, methyl_freq, n_calls)] within the window."""
    points = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["chromosome"] != chrom:
                continue
            pos = int(row["start"])
            if not (start <= pos <= end):
                continue
            points.append((pos, float(row["MethylFreq"]), int(row["NumOfAllCalls"])))
    points.sort()
    return points


def read_dmrs(path, chrom, start, end):
    dmrs = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["chr"] != chrom:
                continue
            dmr_start, dmr_end = int(float(row["start"])), int(float(row["end"]))
            if dmr_end < start or dmr_start > end:
                continue
            dmrs.append(
                (
                    dmr_start,
                    dmr_end,
                    int(float(row["nCG"])),
                    float(row["meanMethy1"]),
                    float(row["meanMethy2"]),
                )
            )
    return sorted(dmrs)


def read_gene_model(gtf_path, gene_name, chrom, start, end):
    """Return (gene_span, strand, exons) for the longest transcript of a gene."""
    gene_re = re.compile(r'gene_name "([^"]+)"')
    tx_re = re.compile(r'transcript_id "([^"]+)"')
    transcripts = {}
    gene_span = None
    strand = "+"

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] != chrom:
                continue
            f_start, f_end = int(fields[3]), int(fields[4])
            if f_end < start or f_start > end:
                continue
            match = gene_re.search(fields[8])
            if not match or match.group(1) != gene_name:
                continue
            strand = fields[6]
            if fields[2] == "gene":
                gene_span = (f_start, f_end)
            elif fields[2] == "exon":
                tx = tx_re.search(fields[8])
                if tx:
                    transcripts.setdefault(tx.group(1), []).append((f_start, f_end))

    if not transcripts:
        return gene_span, strand, []

    longest = max(transcripts.values(), key=lambda exons: sum(e - s for s, e in exons))
    return gene_span, strand, sorted(longest)


def shade_dmrs(ax, dmrs):
    """Shade each DMR, then re-draw its outline above the data.

     The fill is beneath the data and the outline on top.
    """
    for dmr_start, dmr_end, *_ in dmrs:
        width = dmr_end - dmr_start
        ax.add_patch(
            Rectangle(
                (dmr_start, 0),
                width,
                1,
                transform=ax.get_xaxis_transform(),
                facecolor=DMR_FACE,
                edgecolor="none",
                zorder=0,
            )
        )
        ax.add_patch(
            Rectangle(
                (dmr_start, 0),
                width,
                1,
                transform=ax.get_xaxis_transform(),
                facecolor="none",
                edgecolor=DMR_EDGE,
                linewidth=1.0,
                zorder=6,
            )
        )


def draw_methyl_track(ax, points, hap, chrom, start, end, dmrs):
    shade_dmrs(ax, dmrs)
    if points:
        ax.vlines(
            [p[0] for p in points],
            0,
            [p[1] for p in points],
            color=HAP_COLOUR[hap],
            linewidth=0.55,
            zorder=2,
        )
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1.02)
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel(f"{hap}\n5mC freq.", fontsize=7)
    ax.tick_params(axis="y", labelsize=6)
    ax.tick_params(axis="x", length=0, labelbottom=False)
    ax.spines[["top", "right"]].set_visible(False)


def draw_gene_track(ax, locus, gene_model, chrom, start, end, dmrs):
    shade_dmrs(ax, dmrs)
    gene_span, strand, exons = gene_model
    min_exon_width = (end - start) / 500

    if gene_span:
        ax.hlines(0.55, gene_span[0], gene_span[1], color="black", linewidth=0.7, zorder=2)
        for exon_start, exon_end in exons:
            ax.add_patch(
                Rectangle(
                    (exon_start, 0.3),
                    max(exon_end - exon_start, min_exon_width),
                    0.5,
                    facecolor="black",
                    edgecolor="none",
                    zorder=3,
                )
            )
        promoter = gene_span[0] if strand == "+" else gene_span[1]
        ax.annotate(
            "",
            xy=(promoter, 0.9),
            xytext=(promoter, 1.9),
            arrowprops=dict(arrowstyle="-|>", color="black", linewidth=0.9),
            annotation_clip=False,
        )
        label_x = 0.012 if strand == "+" else 0.988
        ax.text(
            label_x,
            0.06,
            f"{locus} ({strand})",
            transform=ax.transAxes,
            ha="left" if strand == "+" else "right",
            va="bottom",
            fontsize=6.5,
            style="italic",
        )

    ax.set_xlim(start, end)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=6)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.set_xlabel(f"{chrom}:{start:,}-{end:,}", fontsize=6.5, labelpad=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out-pdf", required=True)
    parser.add_argument("--out-png", required=True)
    parser.add_argument("--out-source-data", required=True)
    parser.add_argument("--results-root", default="results")
    args = parser.parse_args()

    with open(args.config) as handle:
        config = yaml.safe_load(handle)

    loci = config["loci"]
    order = [name for name in ("XIST", "OTC") if name in loci]
    order += [name for name in sorted(loci) if name not in order]

    fig = plt.figure(figsize=(7.2, 3.1 * len(order)))
    outer = fig.add_gridspec(len(order), 1, hspace=0.42, top=0.93, bottom=0.07)
    source_rows = []

    for index, locus in enumerate(order):
        spec = loci[locus]
        chrom = spec["chrom"]
        start, end = spec["plot_region"]["start"], spec["plot_region"]["end"]
        swap = bool(spec.get("display_swap_haplotypes", False))

        hap_points = {}
        for hap in ("HP1", "HP2"):
            path = (
                f"{args.results_root}/{args.sample}/{locus}/nanomethphase/"
                f"{args.sample}_{locus}_NanoMethPhase_{hap}_MethylFrequency.tsv"
            )
            hap_points[hap] = read_methylfreq(path, chrom, start, end)
            for pos, freq, calls in hap_points[hap]:
                source_rows.append((locus, hap, chrom, pos, pos + 1, freq, calls))

        dmrs = read_dmrs(
            f"{args.results_root}/{args.sample}/{locus}/dma/"
            f"{args.sample}_{locus}_callDMR.txt",
            chrom,
            start,
            end,
        )
        gene_model = read_gene_model(args.gtf, locus, chrom, start, end)

        inner = outer[index].subgridspec(3, 1, height_ratios=[1, 1, 0.34], hspace=0.12)
        ax_top = fig.add_subplot(inner[0])
        ax_bottom = fig.add_subplot(inner[1])
        ax_gene = fig.add_subplot(inner[2])

        row_order = ("HP2", "HP1") if swap else ("HP1", "HP2")
        for ax, hap in zip((ax_top, ax_bottom), row_order):
            draw_methyl_track(ax, hap_points[hap], hap, chrom, start, end, dmrs)
        draw_gene_track(ax_gene, locus, gene_model, chrom, start, end, dmrs)

        heading = f"{PANEL_LABELS[index]}   {locus}"
        if dmrs:
            heading += "    DMR " + "; ".join(
                f"{chrom}:{s:,}-{e:,}  {n} CpG  "
                f"mean 5mC {m1:.2f} vs {m2:.2f}"
                for s, e, n, m1, m2 in dmrs
            )
        ax_top.set_title(heading, fontsize=8, loc="left", pad=4)

    fig.suptitle(
        f"{args.sample}: allele-specific 5mC methylation at XIST and OTC "
        "(reproduction of Figure 1B, 1C)",
        fontsize=9,
    )
    footnote = (
        "Grey band: DSS-called DMR. Arrow: transcription start. "
        "Haplotype labels are arbitrary within each WhatsHap phase block."
    )
    fig.text(0.5, 0.012, footnote, ha="center", fontsize=6, color="#444444")

    fig.savefig(args.out_pdf, bbox_inches="tight")
    fig.savefig(args.out_png, dpi=300, bbox_inches="tight")

    with open(args.out_source_data, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["locus", "haplotype", "chrom", "start", "end", "methyl_freq", "n_calls"]
        )
        writer.writerows(source_rows)


if __name__ == "__main__":
    main()
