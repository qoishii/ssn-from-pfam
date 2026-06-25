#!/usr/bin/env python3
"""Compute all-vs-all protein similarity for SSN construction."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


BLAST_COLUMNS = [
    "query",
    "target",
    "pident",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "score",
    "bitscore",
    "qlen",
    "slen",
]

BLAST_OUTFMT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "score",
    "bitscore",
    "qlen",
    "slen",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run all-vs-all protein similarity using BLASTP or MMseqs2."
    )
    parser.add_argument("--fasta", type=Path, required=True, help="Filtered FASTA input.")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory.")
    parser.add_argument(
        "--method",
        choices=["auto", "mmseqs", "blastp"],
        default="auto",
        help="Similarity engine. Default: auto.",
    )
    parser.add_argument(
        "--min-seq-id",
        type=float,
        default=0.3,
        help="Minimum sequence identity for MMseqs2 prefiltering. Default: 0.3.",
    )
    parser.add_argument(
        "--evalue",
        default="1e-5",
        help="Maximum e-value threshold. Default: 1e-5.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=2,
        help="Worker threads for external aligner. Default: 2.",
    )
    return parser.parse_args()


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing input FASTA: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Input FASTA is empty: {path}")


def choose_method(requested: str) -> str:
    if requested == "mmseqs" and not shutil.which("mmseqs"):
        raise RuntimeError("MMseqs2 was requested, but 'mmseqs' is not installed.")
    if requested == "blastp" and not shutil.which("blastp"):
        raise RuntimeError("BLASTP was requested, but 'blastp' is not installed.")
    if requested != "auto":
        return requested

    if shutil.which("blastp") and shutil.which("makeblastdb"):
        return "blastp"
    if shutil.which("mmseqs"):
        return "mmseqs"

    raise RuntimeError(
        "No supported protein similarity engine is installed. Install MMseqs2 "
        "or NCBI BLAST+. Recommended for EFI-like SSNs: brew install blast"
    )


def run(command: list[str]) -> None:
    print("+ " + " ".join(command))
    subprocess.run(command, check=True)


def normalize_blast_like_table(input_path: Path, output_path: Path) -> int:
    kept = 0
    with input_path.open("r", encoding="utf-8", newline="") as src, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        writer = csv.writer(dst, delimiter="\t")
        writer.writerow(BLAST_COLUMNS + ["qcov", "scov", "coverage"])
        for raw in src:
            if not raw.strip():
                continue
            row = raw.rstrip("\n").split("\t")
            if len(row) != len(BLAST_COLUMNS):
                raise ValueError(
                    f"Expected {len(BLAST_COLUMNS)} columns in {input_path}, got {len(row)}."
                )
            query, target = row[0], row[1]
            if query == target:
                continue
            aln_len = float(row[3])
            qlen = float(row[12])
            slen = float(row[13])
            qcov = aln_len / qlen if qlen else 0.0
            scov = aln_len / slen if slen else 0.0
            coverage = min(qcov, scov)
            writer.writerow(row + [f"{qcov:.6f}", f"{scov:.6f}", f"{coverage:.6f}"])
            kept += 1
    return kept


def run_mmseqs(args: argparse.Namespace, outdir: Path) -> tuple[Path, int]:
    workdir = outdir / "mmseqs_work"
    tmpdir = workdir / "tmp"
    workdir.mkdir(parents=True, exist_ok=True)
    tmpdir.mkdir(parents=True, exist_ok=True)

    seqdb = workdir / "seqdb"
    result = workdir / "result"
    raw_tsv = outdir / "similarity.raw.tsv"
    normalized_tsv = outdir / "similarity.tsv"

    run(["mmseqs", "createdb", str(args.fasta), str(seqdb)])
    run(
        [
            "mmseqs",
            "search",
            str(seqdb),
            str(seqdb),
            str(result),
            str(tmpdir),
            "--min-seq-id",
            str(args.min_seq_id),
            "-e",
            str(args.evalue),
            "--threads",
            str(args.threads),
        ]
    )
    run(
        [
            "mmseqs",
            "convertalis",
            str(seqdb),
            str(seqdb),
            str(result),
            str(raw_tsv),
            "--format-output",
            "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,raw,bits,qlen,tlen",
        ]
    )
    hit_count = normalize_blast_like_table(raw_tsv, normalized_tsv)
    return normalized_tsv, hit_count


def run_blastp(args: argparse.Namespace, outdir: Path) -> tuple[Path, int]:
    db_prefix = outdir / "blastdb" / "seqdb"
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    raw_tsv = outdir / "similarity.raw.tsv"
    normalized_tsv = outdir / "similarity.tsv"

    run(["makeblastdb", "-in", str(args.fasta), "-dbtype", "prot", "-out", str(db_prefix)])
    outfmt = "6 " + " ".join(BLAST_OUTFMT_FIELDS)
    run(
        [
            "blastp",
            "-query",
            str(args.fasta),
            "-db",
            str(db_prefix),
            "-out",
            str(raw_tsv),
            "-outfmt",
            outfmt,
            "-evalue",
            str(args.evalue),
            "-num_threads",
            str(args.threads),
        ]
    )
    hit_count = normalize_blast_like_table(raw_tsv, normalized_tsv)
    return normalized_tsv, hit_count


def main() -> int:
    args = parse_args()
    require_file(args.fasta)
    args.outdir.mkdir(parents=True, exist_ok=True)

    method = choose_method(args.method)
    if method == "mmseqs":
        similarity_path, hit_count = run_mmseqs(args, args.outdir)
    else:
        similarity_path, hit_count = run_blastp(args, args.outdir)

    summary_path = args.outdir / "similarity_summary.txt"
    summary = "\n".join(
        [
            "Similarity computation summary",
            "==============================",
            f"Input FASTA: {args.fasta}",
            f"Method: {method}",
            f"E-value threshold: {args.evalue}",
            f"MMseqs2 min sequence identity, if used: {args.min_seq_id}",
            f"Non-self hits: {hit_count}",
            f"Similarity table: {similarity_path}",
            "",
        ]
    )
    summary_path.write_text(summary, encoding="utf-8")
    print(summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
