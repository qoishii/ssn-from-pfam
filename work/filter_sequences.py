#!/usr/bin/env python3
"""Filter UniProt FASTA sequences and matching metadata for SSN construction."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class FastaRecord:
    accession: str
    header: str
    sequence: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter FASTA and metadata before sequence similarity network construction."
    )
    parser.add_argument("--fasta", type=Path, required=True, help="Input FASTA file.")
    parser.add_argument("--metadata", type=Path, required=True, help="Input UniProt metadata TSV.")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--min-length", type=int, default=80, help="Minimum sequence length.")
    parser.add_argument("--max-length", type=int, default=2000, help="Maximum sequence length.")
    parser.add_argument(
        "--max-ambiguous-fraction",
        type=float,
        default=0.05,
        help="Maximum fraction of non-standard amino acid characters. Default: 0.05.",
    )
    parser.add_argument(
        "--keep-fragments",
        action="store_true",
        help="Keep entries whose FASTA header or protein name marks them as fragments.",
    )
    return parser.parse_args()


def accession_from_header(header: str) -> str:
    first = header.split()[0].lstrip(">")
    parts = first.split("|")
    if len(parts) >= 2:
        return parts[1]
    return first


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                sequence = "".join(chunks).upper()
                records.append(FastaRecord(accession_from_header(header), header, sequence))
            header = line
            chunks = []
        else:
            chunks.append(line)

    if header is not None:
        sequence = "".join(chunks).upper()
        records.append(FastaRecord(accession_from_header(header), header, sequence))

    return records


def read_metadata(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Entry" not in reader.fieldnames:
            raise ValueError("Metadata TSV must include a UniProt 'Entry' column.")
        rows = {row["Entry"]: row for row in reader}
        return reader.fieldnames, rows


def is_fragment(record: FastaRecord, metadata_row: dict[str, str] | None) -> bool:
    haystack = record.header.lower()
    if metadata_row:
        haystack += " " + metadata_row.get("Protein names", "").lower()
    return "fragment" in haystack


def ambiguous_fraction(sequence: str) -> float:
    if not sequence:
        return 1.0
    ambiguous = sum(1 for residue in sequence if residue not in STANDARD_AA)
    return ambiguous / len(sequence)


def wrap_sequence(sequence: str, width: int = 60) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def write_fasta(records: list[FastaRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f"{record.header}\n")
            handle.write(f"{wrap_sequence(record.sequence)}\n")


def write_metadata(
    records: list[FastaRecord],
    metadata_fields: list[str],
    metadata_by_accession: dict[str, dict[str, str]],
    path: Path,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_fields, delimiter="\t")
        writer.writeheader()
        for record in records:
            row = metadata_by_accession.get(record.accession)
            if row:
                writer.writerow(row)


def main() -> int:
    args = parse_args()
    if args.min_length < 1:
        raise ValueError("--min-length must be at least 1.")
    if args.max_length < args.min_length:
        raise ValueError("--max-length must be greater than or equal to --min-length.")
    if not 0 <= args.max_ambiguous_fraction <= 1:
        raise ValueError("--max-ambiguous-fraction must be between 0 and 1.")

    records = read_fasta(args.fasta)
    metadata_fields, metadata_by_accession = read_metadata(args.metadata)

    kept: list[FastaRecord] = []
    seen_sequences: dict[str, str] = {}
    reasons: dict[str, int] = {
        "too_short": 0,
        "too_long": 0,
        "fragment": 0,
        "ambiguous": 0,
        "duplicate_sequence": 0,
        "missing_metadata": 0,
    }

    for record in records:
        metadata_row = metadata_by_accession.get(record.accession)
        length = len(record.sequence)
        if metadata_row is None:
            reasons["missing_metadata"] += 1
            continue
        if length < args.min_length:
            reasons["too_short"] += 1
            continue
        if length > args.max_length:
            reasons["too_long"] += 1
            continue
        if not args.keep_fragments and is_fragment(record, metadata_row):
            reasons["fragment"] += 1
            continue
        if ambiguous_fraction(record.sequence) > args.max_ambiguous_fraction:
            reasons["ambiguous"] += 1
            continue
        if record.sequence in seen_sequences:
            reasons["duplicate_sequence"] += 1
            continue
        seen_sequences[record.sequence] = record.accession
        kept.append(record)

    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_out = args.outdir / "sequences.filtered.fasta"
    metadata_out = args.outdir / "metadata.filtered.tsv"
    summary_out = args.outdir / "filter_summary.txt"

    write_fasta(kept, fasta_out)
    write_metadata(kept, metadata_fields, metadata_by_accession, metadata_out)

    summary_lines = [
        "Sequence filtering summary",
        "==========================",
        f"Input FASTA: {args.fasta}",
        f"Input metadata: {args.metadata}",
        f"Input records: {len(records)}",
        f"Kept records: {len(kept)}",
        f"Minimum length: {args.min_length}",
        f"Maximum length: {args.max_length}",
        f"Maximum ambiguous fraction: {args.max_ambiguous_fraction}",
        f"Fragments kept: {args.keep_fragments}",
        "Removed records:",
    ]
    summary_lines.extend(f"  {reason}: {count}" for reason, count in reasons.items())
    summary_lines.extend(
        [
            f"Filtered FASTA: {fasta_out}",
            f"Filtered metadata: {metadata_out}",
            "",
        ]
    )
    summary = "\n".join(summary_lines)
    summary_out.write_text(summary, encoding="utf-8")
    print(summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
