"""
Download sequences from NCBI GenBank Nucleotide using a search query.
Saves results in GenBank format with optional checkpoint/resume.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.gb_utils import fetch_seq_from_Nucleotide


def main():
    parser = argparse.ArgumentParser(
        description="Download sequences from GenBank Nucleotide by search query."
    )
    parser.add_argument(
        "query",
        help="NCBI search query (e.g., '\"txid39733\"[Organism]')"
    )
    parser.add_argument(
        "outfile",
        help="Output file name (GenBank format)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Number of IDs per efetch request (max ~200, default 100)"
    )
    parser.add_argument(
        "--checkpoint",
        help="Checkpoint file to resume interrupted download"
    )
    args = parser.parse_args()

    fetch_seq_from_Nucleotide(
        query=args.query,
        outfile=args.outfile,
        batch_size=args.batch_size,
        checkpoint_file=args.checkpoint
    )


if __name__ == "__main__":
    main()