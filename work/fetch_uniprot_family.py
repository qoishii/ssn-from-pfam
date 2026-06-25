#!/usr/bin/env python3
"""Fetch protein sequences and metadata from UniProt for an enzyme family query."""

from __future__ import annotations

import argparse
import csv
import sys
import textwrap
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_FIELDS = [
    "accession",
    "id",
    "protein_name",
    "organism_name",
    "length",
    "reviewed",
    "ec",
    "gene_names",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download FASTA sequences and metadata from UniProt.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples:
              python work/fetch_uniprot_family.py \\
                --query 'protein_name:laccase AND reviewed:true' \\
                --max-records 50 \\
                --outdir outputs/uniprot_laccase

              python work/fetch_uniprot_family.py \\
                --query 'xref:pfam-PF00128 AND reviewed:true' \\
                --max-records 200 \\
                --outdir outputs/pfam_pf00128
            """
        ),
    )
    parser.add_argument(
        "--query",
        required=True,
        help="UniProt query, such as 'protein_name:laccase AND reviewed:true'.",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=200,
        help="Maximum records to retrieve. Default: 200.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Directory where sequences.fasta, metadata.tsv, and summary.txt are written.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="HTTP timeout in seconds. Default: 60.",
    )
    return parser.parse_args()


def fetch_uniprot(query: str, size: int, output_format: str, timeout: int) -> str:
    params = {
        "query": query,
        "format": output_format,
        "size": size,
    }
    if output_format == "tsv":
        params["fields"] = ",".join(DEFAULT_FIELDS)

    url = f"{UNIPROT_SEARCH_URL}?{urlencode(params)}"
    request = Request(url, headers={"User-Agent": "codex-ssn-workflow/0.1"})

    try:
        with urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"UniProt request failed with HTTP {exc.code}: {detail}") from exc
    except URLError as exc:
        raise RuntimeError(f"Could not reach UniProt: {exc.reason}") from exc


def count_fasta_records(fasta_text: str) -> int:
    return sum(1 for line in fasta_text.splitlines() if line.startswith(">"))


def count_tsv_records(tsv_text: str) -> int:
    rows = list(csv.reader(tsv_text.splitlines(), delimiter="\t"))
    if not rows:
        return 0
    return max(0, len(rows) - 1)


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def validate_inputs(query: str, max_records: int) -> None:
    if not query.strip():
        raise ValueError("--query cannot be empty.")
    if max_records < 1:
        raise ValueError("--max-records must be at least 1.")
    if max_records > 500:
        raise ValueError(
            "--max-records is capped at 500 for this first workflow step. "
            "Raise this later once batching is implemented."
        )


def nonempty_lines(text: str) -> Iterable[str]:
    for line in text.splitlines():
        if line.strip():
            yield line


def main() -> int:
    args = parse_args()
    validate_inputs(args.query, args.max_records)

    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_path = args.outdir / "sequences.fasta"
    metadata_path = args.outdir / "metadata.tsv"
    summary_path = args.outdir / "summary.txt"

    metadata = fetch_uniprot(args.query, args.max_records, "tsv", args.timeout)
    fasta = fetch_uniprot(args.query, args.max_records, "fasta", args.timeout)

    metadata_count = count_tsv_records(metadata)
    fasta_count = count_fasta_records(fasta)

    if metadata_count == 0 or fasta_count == 0:
        raise RuntimeError(
            "UniProt returned no records. Try a broader query or remove reviewed:true."
        )

    if metadata_count != fasta_count:
        raise RuntimeError(
            f"Metadata/FASTA record mismatch: metadata={metadata_count}, fasta={fasta_count}."
        )

    write_text(metadata_path, metadata)
    write_text(fasta_path, fasta)

    summary = "\n".join(
        [
            "UniProt retrieval summary",
            "=========================",
            f"Query: {args.query}",
            f"Requested max records: {args.max_records}",
            f"Retrieved records: {fasta_count}",
            f"Metadata fields: {', '.join(DEFAULT_FIELDS)}",
            f"FASTA file: {fasta_path}",
            f"Metadata file: {metadata_path}",
            "",
        ]
    )
    write_text(summary_path, summary)

    print(summary)
    print("First metadata lines:")
    print("\n".join(list(nonempty_lines(metadata))[:3]))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
