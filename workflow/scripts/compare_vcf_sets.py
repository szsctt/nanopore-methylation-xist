#!/usr/bin/env python3
import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

VariantKey = Tuple[str, int, str, str]


@dataclass
class Variant:
    key: VariantKey
    qual: float
    is_snp: bool


def open_text(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


def parse_qual(raw: str) -> float:
    if raw == '.' or raw == '':
        return 0.0
    try:
        return float(raw)
    except ValueError:
        return 0.0


def is_snp(ref: str, alt: str) -> bool:
    if ',' in alt:
        return False
    return len(ref) == 1 and len(alt) == 1


def load_variants(path: str, chrom: str) -> List[Variant]:
    variants: List[Variant] = []
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            c = parts[0]
            if c != chrom:
                continue
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            qual = parse_qual(parts[5])
            key = (c, pos, ref, alt)
            variants.append(Variant(key=key, qual=qual, is_snp=is_snp(ref, alt)))
    return variants


def load_truth_keys(path: str, chrom: str, snp_only: bool) -> Set[VariantKey]:
    keys: Set[VariantKey] = set()
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            c = parts[0]
            if c != chrom:
                continue
            ref = parts[3]
            alt = parts[4]
            if snp_only and not is_snp(ref, alt):
                continue
            key = (c, int(parts[1]), ref, alt)
            keys.add(key)
    return keys


def f1(precision: float, recall: float) -> float:
    if precision + recall == 0:
        return 0.0
    return 2.0 * precision * recall / (precision + recall)


def metric_row(candidate: Set[VariantKey], truth: Set[VariantKey]) -> Dict[str, float]:
    shared = len(candidate & truth)
    precision = shared / len(candidate) if candidate else 0.0
    recall = shared / len(truth) if truth else 0.0
    return {
        'candidate': len(candidate),
        'truth': len(truth),
        'shared': shared,
        'precision': precision,
        'recall': recall,
        'f1': f1(precision, recall),
    }


def write_threshold_scan(path: str, variants: List[Variant], truth_keys: Set[VariantKey], snp_only: bool) -> None:
    filtered = [v for v in variants if (v.is_snp if snp_only else True)]
    total = len(filtered)
    by_qual = sorted(filtered, key=lambda v: v.qual, reverse=True)
    thresholds: List[float] = []
    seen = set()
    for v in by_qual:
        if v.qual not in seen:
            thresholds.append(v.qual)
            seen.add(v.qual)

    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow([
            'min_q', 'candidate', 'truth', 'shared', 'precision', 'recall', 'f1', 'fraction_candidates_kept'
        ])
        for q in thresholds:
            kept = {v.key for v in filtered if v.qual >= q}
            m = metric_row(kept, truth_keys)
            frac = (len(kept) / total) if total else 0.0
            writer.writerow([
                f'{q:g}',
                m['candidate'],
                m['truth'],
                m['shared'],
                f"{m['precision']:.6g}",
                f"{m['recall']:.6g}",
                f"{m['f1']:.6g}",
                f"{frac:.6g}",
            ])


def write_missing_truth(path: str, candidate: Set[VariantKey], truth: Set[VariantKey]) -> None:
    missing = sorted(truth - candidate)
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['chrom', 'pos', 'ref', 'alt'])
        for row in missing:
            writer.writerow(row)


def main() -> None:
    p = argparse.ArgumentParser(description='Compare candidate and truth VCF sets in one chromosome.')
    p.add_argument('--candidate-vcf', required=True)
    p.add_argument('--truth-vcf', required=True)
    p.add_argument('--chrom', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--version', required=True)
    p.add_argument('--comparison-tsv', required=True)
    p.add_argument('--threshold-tsv', required=True)
    p.add_argument('--missing-tsv', required=True)
    args = p.parse_args()

    variants = load_variants(args.candidate_vcf, args.chrom)
    truth_any = load_truth_keys(args.truth_vcf, args.chrom, snp_only=False)
    truth_snp = load_truth_keys(args.truth_vcf, args.chrom, snp_only=True)

    cand_any = {v.key for v in variants}
    cand_snp = {v.key for v in variants if v.is_snp}

    m_any = metric_row(cand_any, truth_any)
    m_snp = metric_row(cand_snp, truth_snp)

    # Best SNP threshold (for legacy matching where truth is SNP-only).
    unique_q = sorted({v.qual for v in variants if v.is_snp}, reverse=True)
    best = {
        'min_q': None,
        'candidate': 0,
        'shared': 0,
        'precision': 0.0,
        'recall': 0.0,
        'f1': 0.0,
    }
    best_keep_all = {'min_q': None}

    for q in unique_q:
        kept = {v.key for v in variants if v.is_snp and v.qual >= q}
        m = metric_row(kept, truth_snp)
        if (m['f1'] > best['f1']) or (m['f1'] == best['f1'] and m['precision'] > best['precision']):
            best = {
                'min_q': q,
                'candidate': m['candidate'],
                'shared': m['shared'],
                'precision': m['precision'],
                'recall': m['recall'],
                'f1': m['f1'],
            }
        if m['shared'] == m_snp['shared']:
            best_keep_all['min_q'] = q

    write_threshold_scan(args.threshold_tsv, variants, truth_snp, snp_only=True)
    write_missing_truth(args.missing_tsv, cand_snp, truth_snp)

    Path(args.comparison_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.comparison_tsv, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow([
            'sample', 'version', 'chrom',
            'candidate_any', 'truth_any', 'shared_any', 'precision_any', 'recall_any', 'f1_any',
            'candidate_snp', 'truth_snp', 'shared_snp', 'precision_snp', 'recall_snp', 'f1_snp',
            'best_min_q_snp', 'best_candidate_snp', 'best_shared_snp', 'best_precision_snp', 'best_recall_snp', 'best_f1_snp',
            'max_min_q_keep_all_shared_snp'
        ])
        writer.writerow([
            args.sample, args.version, args.chrom,
            m_any['candidate'], m_any['truth'], m_any['shared'], f"{m_any['precision']:.6g}", f"{m_any['recall']:.6g}", f"{m_any['f1']:.6g}",
            m_snp['candidate'], m_snp['truth'], m_snp['shared'], f"{m_snp['precision']:.6g}", f"{m_snp['recall']:.6g}", f"{m_snp['f1']:.6g}",
            '' if best['min_q'] is None else f"{best['min_q']:g}",
            best['candidate'], best['shared'], f"{best['precision']:.6g}", f"{best['recall']:.6g}", f"{best['f1']:.6g}",
            '' if best_keep_all['min_q'] is None else f"{best_keep_all['min_q']:g}",
        ])


if __name__ == '__main__':
    main()
