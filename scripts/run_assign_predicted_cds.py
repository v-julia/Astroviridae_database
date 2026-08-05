#!/usr/bin/env python3
"""
Assign ORF types to Prodigal‑predicted ORFs using Pfam domain hits.
Calls assign_predicted_orfs from src.hmmer_utils.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.hmmer_utils import assign_predicted_orfs


def main():
    parser = argparse.ArgumentParser(
        description="Assign ORF types to Prodigal‑predicted ORFs."
    )
    parser.add_argument(
        "protein_fasta",
        help="Protein FASTA file from Prodigal (e.g., *_no_cds_sequences_proteins.faa)"
    )
    parser.add_argument(
        "domtbl",
        help="hmmscan --domtblout file for predicted ORFs"
    )
    parser.add_argument(
        "domain_map",
        help="CSV mapping domain names to ORF types"
    )
    parser.add_argument(
        "output_tsv",
        help="Output TSV with assigned ORF and coordinates"
    )
    args = parser.parse_args()

    assign_predicted_orfs(
        args.protein_fasta,
        args.domtbl,
        args.domain_map,
        args.output_tsv
    )


if __name__ == "__main__":
    main()