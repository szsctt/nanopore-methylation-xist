#!/usr/bin/env python3
"""Compare candidate XIST VCF variants to a truth XIST VCF and scan minQ thresholds."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple


@dataclass
class Metrics:
    candidate: int
    truth: int
    shared: int

    @property
    def precision(self) -> float:
        if self.candidate == 0:
            return 0.0
        return self.shared / self.candidate

    @property
    def recall(self) -> float:
        if self.truth == 0:
            return 0.0
        return self.shared / self.truth

    @property
    def f1(self) -> float:
        p = self.precision
        r = self.recall
        if p + r == 0:
            return 0.0
        return 2 * p * r / (p + r)


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_truth_keys(vcf_path: str) -> Set[str]:
    keys: Set[str] = set()
    with open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            key = "\t".join((fields[0], fields[1], fields[3], fields[4]))
            keys.add(key)
    return keys


def parse_float(raw: str) -> float:
    try:
        return float(raw)
    except (TypeError, ValueError):
        return 0.0


def read_candidate_keys_and_qual(vcf_path: str) -> Dict[str, float]:
    qual_by_key: Dict[str, float] = {}
    with open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            key = "\t".join((fields[0], fields[1], fields[3], fields[4]))
            qual = parse_float(fields[5])
            previous = qual_by_key.get(key)
            if previous is None or qual > previous:
                qual_by_key[key] = qual
    return qual_by_key


def format_float(value: Optional[float]) -> str:
    if value is None:
        return "."
    return f"{value:.6g}"


def threshold_scan(qual_by_key: Dict[str, float], truth_keys: Set[str]) -> List[Dict[str, float]]:
    grouped: Dict[float, List[str]] = defaultdict(list)
    for key, qual in qual_by_key.items():
        grouped[qual].append(key)

    rows: List[Dict[str, float]] = []
    selected = 0
    shared = 0
    candidate_keys = len(qual_by_key)
    truth_n = len(truth_keys)

    for min_q in sorted(grouped.keys(), reverse=True):
        for key in grouped[min_q]:
            selected += 1
            if key in truth_keys:
                shared += 1

        m = Metrics(candidate=selected, truth=truth_n, shared=shared)
        rows.append(
            {
                "min_q": min_q,
                "candidate": selected,
                "truth": truth_n,
                "shared": shared,
                "precision": m.precision,
                "recall": m.recall,
                "f1": m.f1,
                "fraction_candidates_kept": selected / candidate_keys if candidate_keys else 0.0,
            }
        )
    return rows


def choose_best_row(rows: List[Dict[str, float]]) -> Optional[Dict[str, float]]:
    if not rows:
        return None
    # Max F1, then recall, then precision, then higher min_q.
    return max(rows, key=lambda r: (r["f1"], r["recall"], r["precision"], r["min_q"]))


def write_threshold_scan(path: str, rows: List[Dict[str, float]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "min_q",
                "candidate",
                "truth",
                "shared",
                "precision",
                "recall",
                "f1",
                "fraction_candidates_kept",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    format_float(row["min_q"]),
                    int(row["candidate"]),
                    int(row["truth"]),
                    int(row["shared"]),
                    format_float(row["precision"]),
                    format_float(row["recall"]),
                    format_float(row["f1"]),
                    format_float(row["fraction_candidates_kept"]),
                ]
            )


def write_missing_truth(path: str, truth_keys: Set[str], candidate_keys: Set[str]) -> None:
    missing = sorted(truth_keys - candidate_keys)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT"])
        for key in missing:
            writer.writerow(key.split("\t"))


def write_summary(
    path: str,
    sample: str,
    version: str,
    any_m: Metrics,
    best_row: Optional[Dict[str, float]],
    max_min_q_keep_all_shared: Optional[float],
) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "version",
                "candidate_xist",
                "truth_xist",
                "shared_xist_any",
                "precision_any",
                "recall_any",
                "f1_any",
                "max_recall_possible",
                "best_min_q",
                "best_candidate_xist",
                "best_shared_xist",
                "best_precision",
                "best_recall",
                "best_f1",
                "max_min_q_keep_all_shared",
            ]
        )
        writer.writerow(
            [
                sample,
                version,
                any_m.candidate,
                any_m.truth,
                any_m.shared,
                format_float(any_m.precision),
                format_float(any_m.recall),
                format_float(any_m.f1),
                format_float(any_m.recall),
                "." if best_row is None else format_float(best_row["min_q"]),
                0 if best_row is None else int(best_row["candidate"]),
                0 if best_row is None else int(best_row["shared"]),
                "." if best_row is None else format_float(best_row["precision"]),
                "." if best_row is None else format_float(best_row["recall"]),
                "." if best_row is None else format_float(best_row["f1"]),
                format_float(max_min_q_keep_all_shared),
            ]
        )


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--candidate-vcf", required=True)
    parser.add_argument("--truth-vcf", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-threshold-scan", required=True)
    parser.add_argument("--out-missing-truth", required=True)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()

    truth_keys = read_truth_keys(args.truth_vcf)
    qual_by_key = read_candidate_keys_and_qual(args.candidate_vcf)
    candidate_keys = set(qual_by_key.keys())
    shared_keys = candidate_keys & truth_keys

    any_m = Metrics(candidate=len(candidate_keys), truth=len(truth_keys), shared=len(shared_keys))
    rows = threshold_scan(qual_by_key, truth_keys)
    best = choose_best_row(rows)

    if shared_keys:
        max_min_q_keep_all_shared = min(qual_by_key[key] for key in shared_keys)
    else:
        max_min_q_keep_all_shared = None

    write_threshold_scan(args.out_threshold_scan, rows)
    write_missing_truth(args.out_missing_truth, truth_keys, candidate_keys)
    write_summary(
        path=args.out_summary,
        sample=args.sample,
        version=args.version,
        any_m=any_m,
        best_row=best,
        max_min_q_keep_all_shared=max_min_q_keep_all_shared,
    )


if __name__ == "__main__":
    main()
