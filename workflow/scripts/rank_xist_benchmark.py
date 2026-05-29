#!/usr/bin/env python3
"""Rank Clair benchmark versions by match quality to truth XIST variants."""

from __future__ import annotations

import argparse
import csv
from typing import Iterable, List, Optional


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summaries", nargs="+", required=True, help="Input summary TSVs")
    parser.add_argument("--output", required=True, help="Output ranking TSV")
    return parser.parse_args(argv)


def load_rows(paths: List[str]) -> List[dict]:
    rows: List[dict] = []
    for path in paths:
        with open(path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                rows.append(row)
    return rows


def main() -> None:
    args = parse_args()
    rows = load_rows(args.summaries)

    rows.sort(
        key=lambda r: (
            -to_float(r.get("best_f1", "0")),
            -to_float(r.get("f1_any", "0")),
            -to_float(r.get("best_recall", "0")),
            -to_float(r.get("recall_any", "0")),
            r.get("version", ""),
        )
    )

    fieldnames = ["rank"] + (list(rows[0].keys()) if rows else [])
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for idx, row in enumerate(rows, start=1):
            out_row = {"rank": idx}
            out_row.update(row)
            writer.writerow(out_row)


if __name__ == "__main__":
    main()
