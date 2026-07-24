#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


VariantKey = Tuple[str, int, str, str]


@dataclass(frozen=True)
class Record:
    index: int
    key: VariantKey
    qual: float
    filter_value: str
    gt: str
    gq: Optional[float]
    dp: Optional[float]
    af: Optional[float]
    is_snp: bool
    is_indel: bool
    line: str


@dataclass(frozen=True)
class FilterSpec:
    variant_mode: str
    filter_mode: str
    gt_mode: str
    min_qual: Optional[float]
    min_gq: Optional[float]
    min_dp: Optional[float]
    max_dp: Optional[float]
    min_af: Optional[float]
    max_af: Optional[float]


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_float(raw: str) -> Optional[float]:
    if raw in {"", "."}:
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def parse_qual(raw: str) -> float:
    value = parse_float(raw)
    return 0.0 if value is None else value


def parse_sample(format_raw: str, sample_raw: str) -> Dict[str, str]:
    keys = format_raw.split(":") if format_raw else []
    vals = sample_raw.split(":") if sample_raw else []
    return {k: vals[i] if i < len(vals) else "" for i, k in enumerate(keys)}


def is_snp(ref: str, alt: str) -> bool:
    if "," in alt:
        return False
    return len(ref) == 1 and len(alt) == 1


def is_indel(ref: str, alt: str) -> bool:
    for allele in alt.split(","):
        if allele.startswith("<") or len(allele) != len(ref):
            return True
    return False


def parse_record(line: str, index: int, chrom: str) -> Optional[Record]:
    if not line or line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 10 or parts[0] != chrom:
        return None

    sample = parse_sample(parts[8], parts[9])
    key = (parts[0], int(parts[1]), parts[3], parts[4])
    return Record(
        index=index,
        key=key,
        qual=parse_qual(parts[5]),
        filter_value=parts[6],
        gt=sample.get("GT", ""),
        gq=parse_float(sample.get("GQ", "")),
        dp=parse_float(sample.get("DP", "")),
        af=parse_float(sample.get("AF", "")),
        is_snp=is_snp(parts[3], parts[4]),
        is_indel=is_indel(parts[3], parts[4]),
        line=line.rstrip("\n"),
    )


def load_raw_vcf(path: str, chrom: str) -> Tuple[List[str], List[Record]]:
    header: List[str] = []
    records: List[Record] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
                continue
            rec = parse_record(line, len(records), chrom)
            if rec is not None:
                records.append(rec)
    return header, records


def load_target(path: str, chrom: str) -> Dict[VariantKey, str]:
    target: Dict[VariantKey, str] = {}
    with open_text(path) as fh:
        for line in fh:
            rec = parse_record(line, len(target), chrom)
            if rec is not None:
                target[rec.key] = rec.line
    return target


def bit_for(index: int) -> int:
    return 1 << index


def build_condition_masks(records: List[Record], thresholds: Dict[str, List[Optional[float]]]) -> Dict[Tuple[str, object], int]:
    masks: Dict[Tuple[str, object], int] = {}
    all_mask = (1 << len(records)) - 1

    masks[("variant_mode", "ANY")] = all_mask
    masks[("variant_mode", "SNP_ONLY")] = 0
    masks[("variant_mode", "EXCLUDE_INDEL")] = 0

    masks[("filter_mode", "ANY")] = all_mask
    masks[("filter_mode", "DOT")] = 0
    masks[("filter_mode", "PASS")] = 0

    for gt_mode in ["ANY", "CALLED", "HET", "HET_UNPHASED", "HET_PHASED", "HOM_ALT"]:
        masks[("gt_mode", gt_mode)] = 0 if gt_mode != "ANY" else all_mask

    for field in ["min_qual", "min_gq", "min_dp", "max_dp", "min_af", "max_af"]:
        for threshold in thresholds[field]:
            masks[(field, threshold)] = all_mask if threshold is None else 0

    for rec in records:
        b = bit_for(rec.index)
        if rec.is_snp:
            masks[("variant_mode", "SNP_ONLY")] |= b
        if not rec.is_indel:
            masks[("variant_mode", "EXCLUDE_INDEL")] |= b

        if rec.filter_value == ".":
            masks[("filter_mode", "DOT")] |= b
        if rec.filter_value == "PASS":
            masks[("filter_mode", "PASS")] |= b

        gt = rec.gt.replace("|", "/")
        alleles = gt.split("/") if gt else []
        called = bool(alleles) and all(a not in {"", "."} for a in alleles)
        het = called and len(alleles) == 2 and alleles[0] != alleles[1]
        hom_alt = called and len(alleles) == 2 and alleles[0] == alleles[1] and alleles[0] != "0"
        if called:
            masks[("gt_mode", "CALLED")] |= b
        if het:
            masks[("gt_mode", "HET")] |= b
        if rec.gt in {"0/1", "1/0"}:
            masks[("gt_mode", "HET_UNPHASED")] |= b
        if rec.gt in {"0|1", "1|0"}:
            masks[("gt_mode", "HET_PHASED")] |= b
        if hom_alt:
            masks[("gt_mode", "HOM_ALT")] |= b

        values = {
            "min_qual": rec.qual,
            "min_gq": rec.gq,
            "min_dp": rec.dp,
            "max_dp": rec.dp,
            "min_af": rec.af,
            "max_af": rec.af,
        }
        for field, value in values.items():
            for threshold in thresholds[field]:
                if threshold is None or value is None:
                    continue
                if field.startswith("min_") and value >= threshold:
                    masks[(field, threshold)] |= b
                elif field.startswith("max_") and value <= threshold:
                    masks[(field, threshold)] |= b

    return masks


def numeric_thresholds(values: Iterable[Optional[float]], defaults: Iterable[float], around: Iterable[float]) -> List[Optional[float]]:
    present = sorted({v for v in values if v is not None and not math.isnan(v)})
    chosen = {float(v) for v in defaults}
    for value in around:
        for delta in [-2, -1, 0, 1, 2]:
            candidate = value + delta
            if candidate >= 0:
                chosen.add(float(candidate))
    if present:
        chosen.update({present[0], present[-1]})
    return [None] + sorted(chosen)


def build_thresholds(records: List[Record], target_keys: set) -> Dict[str, List[Optional[float]]]:
    target_records = [r for r in records if r.key in target_keys]
    target_min_qual = min((r.qual for r in target_records), default=0)
    target_min_gq = min((r.gq for r in target_records if r.gq is not None), default=target_min_qual)
    target_min_dp = min((r.dp for r in target_records if r.dp is not None), default=0)

    return {
        "min_qual": numeric_thresholds(
            (r.qual for r in records),
            [0, 1, 2, 10, 20, 50, 100, 120, 150, 500, 700, 740, 750, 751, 752, 800, 850, 886, 900, 1000],
            [target_min_qual],
        ),
        "min_gq": numeric_thresholds(
            (r.gq for r in records),
            [0, 1, 2, 10, 20, 50, 100, 120, 150, 500, 700, 740, 750, 751, 752, 800, 850, 900, 1000],
            [target_min_gq],
        ),
        "min_dp": numeric_thresholds(
            (r.dp for r in records),
            [0, 3, 4, 5, 8, 10, 12, 15, 20, 25, 30],
            [target_min_dp],
        ),
        "max_dp": [None, 20, 25, 30, 40, 50, 75, 100, 250, 500, 1000],
        "min_af": [None, 0.0, 0.01, 0.05, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30, 0.40, 0.50],
        "max_af": [None, 0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0],
    }


def generate_specs(thresholds: Dict[str, List[Optional[float]]]) -> Iterable[FilterSpec]:
    variant_modes = ["ANY", "SNP_ONLY", "EXCLUDE_INDEL"]
    filter_modes = ["ANY", "DOT"]
    gt_modes = ["ANY", "HET_UNPHASED"]
    qual_thresholds = [None, 0.0, 500.0, 700.0, 750.0, 751.0, 752.0, 800.0]

    for variant_mode in variant_modes:
        for filter_mode in filter_modes:
            for gt_mode in gt_modes:
                for min_qual in qual_thresholds:
                    yield FilterSpec(variant_mode, filter_mode, gt_mode, min_qual, None, None, None, None, None)

    # Focus deeper FORMAT sweeps around the plausible Clair filtering boundary.
    for variant_mode in ["SNP_ONLY"]:
        for gt_mode in ["ANY", "HET_UNPHASED"]:
            for min_qual in [750.0, 751.0, 752.0]:
                for min_gq in [750.0, 751.0, 752.0]:
                    yield FilterSpec(variant_mode, "ANY", gt_mode, min_qual, min_gq, None, None, None, None)
                for min_dp in [3.0, 4.0, 5.0]:
                    yield FilterSpec(variant_mode, "ANY", gt_mode, min_qual, None, min_dp, None, None, None)
                for max_dp in [20.0, 25.0, 30.0]:
                    yield FilterSpec(variant_mode, "ANY", gt_mode, min_qual, None, None, max_dp, None, None)
                for min_af in [0.10, 0.20, 0.25]:
                    yield FilterSpec(variant_mode, "ANY", gt_mode, min_qual, None, None, None, min_af, None)
                for max_af in [0.80, 0.90, 1.0]:
                    yield FilterSpec(variant_mode, "ANY", gt_mode, min_qual, None, None, None, None, max_af)


def spec_complexity(spec: FilterSpec) -> int:
    score = 0
    score += spec.variant_mode != "ANY"
    score += spec.filter_mode != "ANY"
    score += spec.gt_mode != "ANY"
    score += spec.min_qual is not None
    score += spec.min_gq is not None
    score += spec.min_dp is not None
    score += spec.max_dp is not None
    score += spec.min_af is not None
    score += spec.max_af is not None
    return int(score)


def spec_mask(spec: FilterSpec, masks: Dict[Tuple[str, object], int]) -> int:
    mask = masks[("variant_mode", spec.variant_mode)]
    mask &= masks[("filter_mode", spec.filter_mode)]
    mask &= masks[("gt_mode", spec.gt_mode)]
    mask &= masks[("min_qual", spec.min_qual)]
    mask &= masks[("min_gq", spec.min_gq)]
    mask &= masks[("min_dp", spec.min_dp)]
    mask &= masks[("max_dp", spec.max_dp)]
    mask &= masks[("min_af", spec.min_af)]
    mask &= masks[("max_af", spec.max_af)]
    return mask


def fmt_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def write_grid(path: str, rows: List[Dict[str, object]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "rank",
        "exact",
        "variant_mode",
        "filter_mode",
        "gt_mode",
        "min_qual",
        "min_gq",
        "min_dp",
        "max_dp",
        "min_af",
        "max_af",
        "candidate",
        "target",
        "shared",
        "only_raw",
        "only_target",
        "precision",
        "recall",
        "f1",
        "complexity",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for i, row in enumerate(rows, start=1):
            out = {k: fmt_value(row.get(k)) for k in fields}
            out["rank"] = i
            writer.writerow(out)


def write_best_yaml(path: str, row: Dict[str, object]) -> None:
    lines = [
        f"exact: {str(row['exact']).lower()}",
        "filter:",
        f"  variant_mode: {row['variant_mode']}",
        f"  filter_mode: {row['filter_mode']}",
        f"  gt_mode: {row['gt_mode']}",
        f"  min_qual: {fmt_value(row['min_qual']) or 'null'}",
        f"  min_gq: {fmt_value(row['min_gq']) or 'null'}",
        f"  min_dp: {fmt_value(row['min_dp']) or 'null'}",
        f"  max_dp: {fmt_value(row['max_dp']) or 'null'}",
        f"  min_af: {fmt_value(row['min_af']) or 'null'}",
        f"  max_af: {fmt_value(row['max_af']) or 'null'}",
        "metrics:",
        f"  candidate: {row['candidate']}",
        f"  target: {row['target']}",
        f"  shared: {row['shared']}",
        f"  only_raw: {row['only_raw']}",
        f"  only_target: {row['only_target']}",
        f"  precision: {row['precision']}",
        f"  recall: {row['recall']}",
        f"  f1: {row['f1']}",
    ]
    Path(path).write_text("\n".join(lines) + "\n")


def iter_mask_indices(mask: int) -> Iterable[int]:
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def write_filtered_vcf(path: str, header: List[str], records: List[Record], keep_keys: set) -> None:
    with gzip.open(path, "wt") as out:
        for line in header:
            out.write(line + "\n")
        for rec in records:
            if rec.key in keep_keys:
                out.write(rec.line + "\n")


def write_diff(path: str, keys: Iterable[VariantKey], record_lines: Dict[VariantKey, str]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["chrom", "pos", "ref", "alt", "record"])
        for key in sorted(keys):
            writer.writerow([key[0], key[1], key[2], key[3], record_lines.get(key, "")])


def write_boundary(path: str, records: List[Record], target_keys: set, best: FilterSpec) -> None:
    thresholds = {
        "qual": best.min_qual,
        "gq": best.min_gq,
        "dp_min": best.min_dp,
        "dp_max": best.max_dp,
        "af_min": best.min_af,
        "af_max": best.max_af,
    }
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["class", "chrom", "pos", "ref", "alt", "QUAL", "GT", "GQ", "DP", "AF", "record"])
        for rec in records:
            near = False
            for field, threshold in thresholds.items():
                if threshold is None:
                    continue
                value = {
                    "qual": rec.qual,
                    "gq": rec.gq,
                    "dp_min": rec.dp,
                    "dp_max": rec.dp,
                    "af_min": rec.af,
                    "af_max": rec.af,
                }[field]
                if value is not None and abs(value - threshold) <= (0.02 if field.startswith("af") else 2):
                    near = True
            if near:
                row_class = "target" if rec.key in target_keys else "raw_only"
                writer.writerow([row_class, rec.key[0], rec.key[1], rec.key[2], rec.key[3], rec.qual, rec.gt, rec.gq, rec.dp, rec.af, rec.line])


def row_sort_key(row: Dict[str, object]):
    variant_preference = {"SNP_ONLY": 0, "EXCLUDE_INDEL": 1, "ANY": 2}
    return (
        0 if row["exact"] else 1,
        row["only_raw"] + row["only_target"],
        row["only_target"],
        row["only_raw"],
        row["complexity"],
        variant_preference.get(row["variant_mode"], 9),
        -row["f1"],
        row["candidate"],
    )


def gt_matches(gt: str, mode: str) -> bool:
    if mode == "ANY":
        return True
    normalized = gt.replace("|", "/")
    alleles = normalized.split("/") if normalized else []
    called = bool(alleles) and all(a not in {"", "."} for a in alleles)
    het = called and len(alleles) == 2 and alleles[0] != alleles[1]
    if mode == "CALLED":
        return called
    if mode == "HET":
        return het
    if mode == "HET_UNPHASED":
        return gt in {"0/1", "1/0"}
    if mode == "HET_PHASED":
        return gt in {"0|1", "1|0"}
    if mode == "HOM_ALT":
        return called and len(alleles) == 2 and alleles[0] == alleles[1] and alleles[0] != "0"
    raise ValueError(f"Unknown gt_mode: {mode}")


def record_passes(rec: Record, spec: FilterSpec) -> bool:
    if spec.variant_mode == "SNP_ONLY" and not rec.is_snp:
        return False
    if spec.variant_mode == "EXCLUDE_INDEL" and rec.is_indel:
        return False
    if spec.filter_mode == "DOT" and rec.filter_value != ".":
        return False
    if spec.filter_mode == "PASS" and rec.filter_value != "PASS":
        return False
    if not gt_matches(rec.gt, spec.gt_mode):
        return False
    if spec.min_qual is not None and rec.qual < spec.min_qual:
        return False
    if spec.min_gq is not None and (rec.gq is None or rec.gq < spec.min_gq):
        return False
    if spec.min_dp is not None and (rec.dp is None or rec.dp < spec.min_dp):
        return False
    if spec.max_dp is not None and (rec.dp is None or rec.dp > spec.max_dp):
        return False
    if spec.min_af is not None and (rec.af is None or rec.af < spec.min_af):
        return False
    if spec.max_af is not None and (rec.af is None or rec.af > spec.max_af):
        return False
    return True


def main() -> None:
    p = argparse.ArgumentParser(description="Reconstruct a historical VCF filter from raw and published VCFs.")
    p.add_argument("--raw-vcf", required=True)
    p.add_argument("--target-vcf", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--grid-tsv", required=True)
    p.add_argument("--top-tsv", required=True)
    p.add_argument("--best-yaml", required=True)
    p.add_argument("--filtered-vcf-gz", required=True)
    p.add_argument("--only-raw-tsv", required=True)
    p.add_argument("--only-target-tsv", required=True)
    p.add_argument("--boundary-tsv", required=True)
    p.add_argument("--top-n", type=int, default=100)
    args = p.parse_args()

    header, records = load_raw_vcf(args.raw_vcf, args.chrom)
    target = load_target(args.target_vcf, args.chrom)
    target_keys = set(target)
    key_to_record = {r.key: r for r in records}

    thresholds = build_thresholds(records, target_keys)

    rows: List[Dict[str, object]] = []
    seen_specs = set()
    for spec in generate_specs(thresholds):
        if spec in seen_specs:
            continue
        seen_specs.add(spec)
        candidate_count = 0
        shared = 0
        for rec in records:
            if record_passes(rec, spec):
                candidate_count += 1
                if rec.key in target_keys:
                    shared += 1
        only_raw = candidate_count - shared
        only_target = len(target_keys) - shared
        precision = shared / candidate_count if candidate_count else 0.0
        recall = shared / len(target_keys) if target_keys else 0.0
        f1 = 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
        rows.append({
            "exact": only_raw == 0 and only_target == 0,
            "variant_mode": spec.variant_mode,
            "filter_mode": spec.filter_mode,
            "gt_mode": spec.gt_mode,
            "min_qual": spec.min_qual,
            "min_gq": spec.min_gq,
            "min_dp": spec.min_dp,
            "max_dp": spec.max_dp,
            "min_af": spec.min_af,
            "max_af": spec.max_af,
            "candidate": candidate_count,
            "target": len(target_keys),
            "shared": shared,
            "only_raw": only_raw,
            "only_target": only_target,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "complexity": spec_complexity(spec),
            "_spec": spec,
        })

    rows.sort(key=row_sort_key)
    write_grid(args.grid_tsv, rows)
    write_grid(args.top_tsv, rows[: args.top_n])

    best = rows[0]
    best_spec = best["_spec"]
    best_keys = {rec.key for rec in records if record_passes(rec, best_spec)}
    write_best_yaml(args.best_yaml, best)
    write_filtered_vcf(args.filtered_vcf_gz, header, records, best_keys)

    raw_only = best_keys - target_keys
    target_only = target_keys - best_keys
    write_diff(args.only_raw_tsv, raw_only, {k: key_to_record[k].line for k in raw_only if k in key_to_record})
    write_diff(args.only_target_tsv, target_only, target)
    write_boundary(args.boundary_tsv, records, target_keys, best_spec)


if __name__ == "__main__":
    main()
