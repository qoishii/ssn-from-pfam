#!/usr/bin/env python3
"""Run the full SSN-from-PFam workflow."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the complete SSN workflow.")
    parser.add_argument("--query", required=True, help="UniProt query, such as xref:interpro-IPR024001.")
    parser.add_argument("--name", required=True, help="Short output name, such as IPR024001.")
    parser.add_argument("--max-records", type=int, default=500, help="Maximum UniProt records to download.")
    parser.add_argument("--min-length", type=int, default=80, help="Minimum protein length.")
    parser.add_argument("--max-length", type=int, default=2000, help="Maximum protein length.")
    parser.add_argument("--max-ambiguous-fraction", type=float, default=0.05)
    parser.add_argument("--min-identity", type=float, default=50.0, help="Minimum percent identity for SSN edges.")
    parser.add_argument("--min-coverage", type=float, default=0.70, help="Minimum alignment coverage for SSN edges.")
    parser.add_argument("--evalue", default="1e-5", help="BLASTP e-value threshold.")
    parser.add_argument("--threads", type=int, default=2, help="BLASTP worker threads.")
    return parser.parse_args()


def run(command: list[str], env: dict[str, str] | None = None) -> None:
    print("+ " + " ".join(command))
    subprocess.run(command, check=True, cwd=ROOT, env=env)


def ensure_blast_available() -> None:
    blastp = subprocess.run(
        ["blastp", "-version"],
        check=False,
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    makeblastdb = subprocess.run(
        ["makeblastdb", "-version"],
        check=False,
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if blastp.returncode == 0 and makeblastdb.returncode == 0:
        return
    raise RuntimeError(
        "BLAST+ is required and must be available on PATH. Install it with "
        "'brew install blast'."
    )


def main() -> int:
    args = parse_args()
    ensure_blast_available()

    raw_dir = Path("outputs") / f"uniprot_{args.name}_all"
    filtered_dir = Path("outputs") / f"uniprot_{args.name}_all_filtered"
    similarity_dir = Path("outputs") / f"uniprot_{args.name}_similarity"
    network_dir = Path("outputs") / f"uniprot_{args.name}_network"

    run(
        [
            sys.executable,
            "work/fetch_uniprot_family.py",
            "--query",
            args.query,
            "--max-records",
            str(args.max_records),
            "--outdir",
            str(raw_dir),
        ]
    )
    run(
        [
            sys.executable,
            "work/filter_sequences.py",
            "--fasta",
            str(raw_dir / "sequences.fasta"),
            "--metadata",
            str(raw_dir / "metadata.tsv"),
            "--outdir",
            str(filtered_dir),
            "--min-length",
            str(args.min_length),
            "--max-length",
            str(args.max_length),
            "--max-ambiguous-fraction",
            str(args.max_ambiguous_fraction),
        ]
    )
    run(
        [
            sys.executable,
            "work/compute_similarity.py",
            "--fasta",
            str(filtered_dir / "sequences.filtered.fasta"),
            "--outdir",
            str(similarity_dir),
            "--method",
            "blastp",
            "--evalue",
            str(args.evalue),
            "--threads",
            str(args.threads),
        ]
    )
    run(
        [
            sys.executable,
            "work/build_network_plot.py",
            "--similarity",
            str(similarity_dir / "similarity.tsv"),
            "--metadata",
            str(filtered_dir / "metadata.filtered.tsv"),
            "--outdir",
            str(network_dir),
            "--min-identity",
            str(args.min_identity),
            "--min-coverage",
            str(args.min_coverage),
            "--max-evalue",
            str(args.evalue),
        ]
    )

    print(f"Workflow complete. See {network_dir / 'network.svg'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
