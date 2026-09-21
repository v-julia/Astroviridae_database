"""
Append nucleotide sequences of fusion-CDS accessions to the no_cds FASTA.

Reads the original GenBank file, extracts sequences for accessions listed
in the fusion table, and writes them into *_no_cds_sequences.fasta so
Prodigal can predict their CDS.
"""
import argparse
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent))


def read_fusion_accessions(fusion_tsv):
    accs = set()
    with open(fusion_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        if "Accession" not in header:
            return accs
        idx = header.index("Accession")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[idx]:
                accs.add(parts[idx])
    return accs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank",     help="Original GenBank file")
    parser.add_argument("fusion_tsv",  help="Tsv table with accessions which annotated CDS are in facted fused ORFs;*_fusion_orfs.tsv")
    parser.add_argument("no_cds_fasta",help="File with nt genome sequences of entries with no annotaded CDS; *_no_cds_sequences.fasta (appended in-place)")
    args = parser.parse_args()

    fusion_accs = read_fusion_accessions(args.fusion_tsv)
    if not fusion_accs:
        print("No fusion accessions — nothing appended.")
        return

    # Read existing accessions to avoid duplicates
    existing = set()
    if Path(args.no_cds_fasta).exists():
        for rec in SeqIO.parse(args.no_cds_fasta, "fasta"):
            existing.add(rec.id.split(".")[0])

    to_add = [rec for rec in SeqIO.parse(args.genbank, "genbank")
              if rec.id.split(".")[0] in fusion_accs and rec.id not in existing]

    if not to_add:
        print("All fusion accessions already present — nothing appended.")
        return

    mode = "a" if Path(args.no_cds_fasta).exists() else "w"
    with open(args.no_cds_fasta, mode) as out:
        SeqIO.write(to_add, out, "fasta")

    print(f"Appended {len(to_add)} fusion sequences to {args.no_cds_fasta}")


if __name__ == "__main__":
    main()