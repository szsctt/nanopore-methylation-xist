#!/usr/bin/env python3
"""Summarize QUAL, DP, and AF fields from a VCF/VCF.GZ file."""

from __future__ import annotations

import argparse
import gzip
import os
from dataclasses import dataclass
from typing import Dict, Iterable, Optional


@dataclass
class MetricStats:
    total_records: int = 0
    present: int = 0
    missing: int = 0
    value_sum: float = 0.0
    value_min: Optional[float] = None
    value_max: Optional[float] = None

    def add(self, value: Optional[float]) -> None:
        self.total_records += 1
        if value is None:
            self.missing += 1
            return

        self.present += 1
        self.value_sum += value
        if self.value_min is None or value < self.value_min:
            self.value_min = value
        if self.value_max is None or value > self.value_max:
            self.value_max = value

    @property
    def mean(self) -> Optional[float]:
        if self.present == 0:
            return None
        return self.value_sum / self.present


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _parse_float(value: str) -> Optional[float]:
    if value in {"", ".", "NA", "nan", "NaN"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _extract_from_info(info: str, key: str) -> Optional[float]:
    prefix = f"{key}="
    for token in info.split(";"):
        if token.startswith(prefix):
            raw = token[len(prefix) :].split(",")[0]
            return _parse_float(raw)
    return None


def _extract_sample_field(fmt: str, sample: str, field: str) -> Optional[str]:
    fmt_keys = fmt.split(":")
    sample_values = sample.split(":")
    if field not in fmt_keys:
        return None

    idx = fmt_keys.index(field)
    if idx >= len(sample_values):
        return None
    return sample_values[idx]


def summarize_vcf(path: str) -> Dict[str, MetricStats]:
    stats = {"QUAL": MetricStats(), "DP": MetricStats(), "AF": MetricStats()}

    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            qual = _parse_float(fields[5])
            info = fields[7]
            fmt = fields[8] if len(fields) > 8 else ""
            sample = fields[9] if len(fields) > 9 else ""

            dp_val = None
            af_val = None

            if fmt and sample:
                dp_raw = _extract_sample_field(fmt, sample, "DP")
                af_raw = _extract_sample_field(fmt, sample, "AF")
                if dp_raw is not None:
                    dp_val = _parse_float(dp_raw)
                if af_raw is not None:
                    af_val = _parse_float(af_raw.split(",")[0])

            if dp_val is None:
                dp_val = _extract_from_info(info, "DP")
            if af_val is None:
                af_val = _extract_from_info(info, "AF")

            stats["QUAL"].add(qual)
            stats["DP"].add(dp_val)
            stats["AF"].add(af_val)

    return stats


def write_summary(path: str, input_vcf: str, stats: Dict[str, MetricStats]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w", encoding="utf-8") as out:
        out.write(
            "vcf\tmetric\ttotal_records\tpresent\tmissing\tmin\tmax\tmean\n"
        )
        for metric in ("QUAL", "DP", "AF"):
            m = stats[metric]
            min_v = "." if m.value_min is None else f"{m.value_min:.6g}"
            max_v = "." if m.value_max is None else f"{m.value_max:.6g}"
            mean_v = "." if m.mean is None else f"{m.mean:.6g}"
            out.write(
                f"{input_vcf}\t{metric}\t{m.total_records}\t{m.present}\t{m.missing}\t{min_v}\t{max_v}\t{mean_v}\n"
            )


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ path")
    parser.add_argument("--out-summary", required=True, help="Output TSV summary")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    stats = summarize_vcf(args.vcf)
    write_summary(args.out_summary, args.vcf, stats)


if __name__ == "__main__":
    main()
