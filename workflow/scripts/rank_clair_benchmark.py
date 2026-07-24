#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def to_float(v: str) -> float:
    try:
        return float(v)
    except Exception:
        return 0.0


def main() -> None:
    p = argparse.ArgumentParser(description='Rank Clair benchmark comparison rows.')
    p.add_argument('--inputs', nargs='+', required=True)
    p.add_argument('--output', required=True)
    p.add_argument(
        '--sort-by',
        choices=['snp', 'any'],
        default='snp',
        help='Rank by SNP metrics or all-variant metrics. Default: snp',
    )
    args = p.parse_args()

    rows = []
    for path in args.inputs:
        with open(path, 'r', newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                rows.append(row)

    if args.sort_by == 'any':
        rows.sort(
            key=lambda r: (
                to_float(r.get('f1_any', '0')),
                to_float(r.get('recall_any', '0')),
                to_float(r.get('precision_any', '0')),
                -to_float(r.get('candidate_any', '0')),
            ),
            reverse=True,
        )
    else:
        rows.sort(
            key=lambda r: (
                to_float(r.get('f1_snp', '0')),
                to_float(r.get('recall_snp', '0')),
                to_float(r.get('precision_snp', '0')),
                -to_float(r.get('candidate_snp', '0')),
            ),
            reverse=True,
        )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        'rank', 'sample', 'version', 'chrom',
        'candidate_any', 'truth_any', 'shared_any', 'precision_any', 'recall_any', 'f1_any',
        'candidate_snp', 'truth_snp', 'shared_snp', 'precision_snp', 'recall_snp', 'f1_snp',
        'best_min_q_snp', 'best_candidate_snp', 'best_shared_snp', 'best_precision_snp', 'best_recall_snp', 'best_f1_snp',
        'max_min_q_keep_all_shared_snp',
    ]

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for i, row in enumerate(rows, start=1):
            out = dict(row)
            out['rank'] = i
            writer.writerow(out)


if __name__ == '__main__':
    main()
