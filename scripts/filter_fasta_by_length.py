#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def filter_fasta_by_length(
    input_fasta: str,
    output_fasta: str,
    min_len: int = 0,
    max_len: int | None = None
):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_len = len(record.seq)
            if seq_len >= min_len and (max_len is None or seq_len <= max_len):
                SeqIO.write(record, out, "fasta")

def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences by length"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file"
    )
    parser.add_argument(
        "-m", "--min-len",
        type=int,
        required=True,
        help="Minimum sequence length"
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=None,
        help="Maximum sequence length (optional)"
    )

    args = parser.parse_args()

    filter_fasta_by_length(
        input_fasta=args.input,
        output_fasta=args.output,
        min_len=args.min_len,
        max_len=args.max_len
    )


if __name__ == "__main__":
    main()
