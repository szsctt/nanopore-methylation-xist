#!/usr/bin/env python3
"""Build a supplementary DMR/DML table carrying raw and model-based values.

DSS reports `mu1`/`mu2` as smoothed, dispersion-weighted, shrunk estimates when
`DMLtest(..., smoothing = TRUE)` is used. The raw values are the counts DSS was given: 
coverage and methylated-call count per CpG per haplotype. This script reports both 
for each CpG inside a called DMR.

Counts come from the ReadyForDSS tables, in which NanoMethPhase has pooled the two
strands of each CpG (sum of the `+` and `-` rows of the MethylFrequency file at the 
shared coordinate). Coverage and methylated calls are therefore both-strand totals, 
one row per CpG, and `nCG` counts CpGs rather than rows.

Region methylation is reported two ways: `raw_pooled`, total methylated calls over
total calls across the DMR, which weights each CpG by its coverage; and DSS
`meanMethy`, the mean of the smoothed `mu` values that the published statistics
use.
"""

import argparse
import csv
import sys


def read_dss_input(path):
    """ReadyForDSS: chr, pos, coverage (N), methylated calls (X), no header."""
    values = {}
    with open(path, newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0] == "chr":
                continue
            values[(row[0], int(row[1]))] = (int(row[2]), int(row[3]))
    return values


def read_table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value, places=6):
    if value is None:
        return "NA"
    return f"{value:.{places}f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locus", required=True)
    parser.add_argument("--dmr", required=True, help="DSS callDMR.txt")
    parser.add_argument("--dmltest", required=True, help="DSS DMLtest.txt")
    parser.add_argument("--case", required=True, help="ReadyForDSS case file (HP1)")
    parser.add_argument("--control", required=True, help="ReadyForDSS control file (HP2)")
    parser.add_argument("--p-threshold", type=float, default=0.001)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--out-dml", required=True)
    parser.add_argument("--out-dmr", required=True)
    args = parser.parse_args()

    hp1 = read_dss_input(args.case)
    hp2 = read_dss_input(args.control)
    by_pos = {(r["chr"], int(r["pos"])): r for r in read_table(args.dmltest)}
    dmrs = read_table(args.dmr)

    dml_rows = []
    dmr_rows = []

    for dmr in dmrs:
        chrom = dmr["chr"]
        start, end = int(float(dmr["start"])), int(float(dmr["end"]))
        inside = sorted(
            (key for key in by_pos if key[0] == chrom and start <= key[1] <= end),
            key=lambda k: k[1],
        )

        raw1_num = raw1_den = raw2_num = raw2_den = 0
        mu1_sum = mu2_sum = 0.0
        n_sig_p = n_sig_fdr = 0

        for key in inside:
            row = by_pos[key]
            cov1, meth1 = hp1.get(key, (0, 0))
            cov2, meth2 = hp2.get(key, (0, 0))
            pval, fdr = float(row["pval"]), float(row["fdr"])
            mu1, mu2 = float(row["mu1"]), float(row["mu2"])

            raw1_num += meth1
            raw1_den += cov1
            raw2_num += meth2
            raw2_den += cov2
            mu1_sum += mu1
            mu2_sum += mu2
            n_sig_p += pval < args.p_threshold
            n_sig_fdr += fdr < args.fdr_threshold

            dml_rows.append(
                [
                    args.locus,
                    f"{chrom}:{start}-{end}",
                    chrom,
                    key[1],
                    cov1,
                    meth1,
                    fmt(meth1 / cov1) if cov1 else "NA",
                    cov2,
                    meth2,
                    fmt(meth2 / cov2) if cov2 else "NA",
                    fmt(mu1),
                    fmt(mu2),
                    fmt(float(row["diff"])),
                    fmt(float(row["diff.se"])),
                    fmt(float(row["stat"])),
                    fmt(float(row["phi1"])),
                    fmt(float(row["phi2"])),
                    f"{pval:.3e}",
                    f"{fdr:.3e}",
                    "yes" if pval < args.p_threshold else "no",
                    "yes" if fdr < args.fdr_threshold else "no",
                ]
            )

        n = len(inside)
        # callDMR sets meanMethy to the unweighted mean of mu over the region.
        # Recomputing it here guards against reading the wrong DMLtest rows.
        for mu_sum, reported in ((mu1_sum, "meanMethy1"), (mu2_sum, "meanMethy2")):
            if n and abs(mu_sum / n - float(dmr[reported])) > 1e-6:
                sys.exit(
                    f"{args.locus} {chrom}:{start}-{end}: recomputed {reported} "
                    f"{mu_sum / n:.6f} does not match callDMR {dmr[reported]}"
                )

        dmr_rows.append(
            [
                args.locus,
                chrom,
                start,
                end,
                dmr["length"],
                dmr["nCG"],
                n,
                n_sig_p,
                n_sig_fdr,
                fmt(raw1_num / raw1_den) if raw1_den else "NA",
                fmt(raw2_num / raw2_den) if raw2_den else "NA",
                fmt(float(dmr["meanMethy1"])),
                fmt(float(dmr["meanMethy2"])),
                fmt(float(dmr["diff.Methy"])),
                fmt(float(dmr["areaStat"]), 4),
            ]
        )

    with open(args.out_dml, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "locus", "dmr", "chr", "pos",
                "HP1_coverage", "HP1_methylated_calls", "HP1_raw_methylation",
                "HP2_coverage", "HP2_methylated_calls", "HP2_raw_methylation",
                "mu1_smoothed", "mu2_smoothed", "diff_smoothed", "diff.se",
                "stat", "phi1", "phi2", "pval", "fdr",
                f"p_lt_{args.p_threshold}", f"fdr_lt_{args.fdr_threshold}",
            ]
        )
        writer.writerows(dml_rows)

    with open(args.out_dmr, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "locus", "chr", "start", "end", "length", "nCG",
                "cpg_in_dmr", f"n_p_lt_{args.p_threshold}", f"n_fdr_lt_{args.fdr_threshold}",
                "HP1_raw_pooled", "HP2_raw_pooled",
                "meanMethy1_DSS", "meanMethy2_DSS", "diff.Methy_DSS", "areaStat_DSS",
            ]
        )
        writer.writerows(dmr_rows)


if __name__ == "__main__":
    main()
