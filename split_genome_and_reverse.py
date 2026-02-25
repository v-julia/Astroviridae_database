import argparse
import os
import copy
import pandas as pd
from Bio import SeqIO


def split_genome_and_reverse(input_file, coord_file):

    coord_df = pd.read_csv(coord_file, sep=",", index_col=0)
    orfs = ['1A', '1B', '2']
    dict_orfs = {orf: [] for orf in orfs}
    seqs = SeqIO.parse(input_file, format='fasta')
    for seq in seqs:
        acc = seq.id.split("_")[0]
        if acc in ['NC', 'AC']:
            acc = '_'.join(seq.id.split("_")[:2])
        for orf in orfs:
            if acc in coord_df.index and coord_df.loc[acc, orf] != 'NA-NA':
                cd = coord_df.loc[acc, orf]
                s, e = map(int, cd.split('-'))
                seq_orf = copy.deepcopy(seq)
                seq_orf.seq = seq.seq[s:e]
                strand_col = f"{orf}-strand"
                if strand_col in coord_df.columns:
                    strand_value = coord_df.loc[acc, strand_col]
                    if strand_value == -1:
                        seq_orf.seq = seq_orf.seq.reverse_complement()
                        seq_orf.id += "_reverse"
                        seq_orf.description = seq_orf.id

                dict_orfs[orf].append(seq_orf)

    out_prefix = os.path.splitext(input_file)[0]

    for orf, records in dict_orfs.items():
        outfile = f"{out_prefix}_{orf}.fasta"
        SeqIO.write(records, outfile, "fasta")
        print(f"Written: {outfile} ({len(records)} sequences)")


def main():
    parser = argparse.ArgumentParser(
        description="Split genome FASTA into ORFs using coordinates and strand info"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file with full genomes"
    )
    parser.add_argument(
        "-c", "--coords",
        required=True,
        help="CSV file with ORF coordinates and strand info"
    )

    args = parser.parse_args()

    split_genome_and_reverse(args.input, args.coords)


if __name__ == "__main__":
    main()
