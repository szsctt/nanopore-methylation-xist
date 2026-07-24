#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def extract_contig(source: Path, contig: str) -> str:
    in_target = False
    chunks = []

    with open_text(source) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if in_target and name != contig:
                    break
                in_target = name == contig
                continue

            if in_target:
                chunks.append(line.strip())

    sequence = "".join(chunks)
    if not sequence:
        raise SystemExit(f"Contig not found or empty in FASTA: {contig}")
    return sequence


def write_fasta_and_fai(sequence: str, contig: str, fasta: Path, fai: Path, line_bases: int) -> None:
    fasta.parent.mkdir(parents=True, exist_ok=True)
    fai.parent.mkdir(parents=True, exist_ok=True)

    header = f">{contig}\n"
    with fasta.open("w") as handle:
        handle.write(header)
        for i in range(0, len(sequence), line_bases):
            handle.write(sequence[i : i + line_bases] + "\n")

    offset = len(header.encode("utf-8"))
    line_width = line_bases + 1
    with fai.open("w") as handle:
        handle.write(f"{contig}\t{len(sequence)}\t{offset}\t{line_bases}\t{line_width}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract one contig from a FASTA and write a matching .fai.")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--contig", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-fai", required=True)
    parser.add_argument("--line-bases", type=int, default=60)
    args = parser.parse_args()

    sequence = extract_contig(Path(args.input_fasta), args.contig)
    write_fasta_and_fai(
        sequence,
        args.contig,
        Path(args.output_fasta),
        Path(args.output_fai),
        args.line_bases,
    )


if __name__ == "__main__":
    main()
