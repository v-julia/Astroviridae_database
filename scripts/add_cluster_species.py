#!/usr/bin/env python3
"""
Propagate species and virus names from cluster members that have WG_2025 annotations
to all sequences in the same cluster. Prioritises non‑NA annotations.
"""
import argparse
import pandas as pd


def clean_value(val):
    if pd.isna(val) or str(val).strip() == '':
        return 'NA'
    return str(val).strip()

def parse_clustering(uc_file):
    """
    Parse mmseq2 .uc file to build:
      - seq_to_rep: mapping sequence -> cluster representative
      - clusters: mapping rep -> list of all sequences in that cluster
    """
    df = pd.read_csv(uc_file, sep='\t', header=None)
    df['seq'] = df[8].str.split('/').str[0]
    df['rep'] = df[9].str.split('/').str[0]

    seq_to_rep = dict(zip(df['seq'], df['rep']))

    clusters = {}
    for _, row in df.iterrows():
        rep = row['rep']
        seq = row['seq']
        clusters.setdefault(rep, set()).add(seq)
        clusters[rep].add(rep)          # ensure representative is included

    return seq_to_rep, {rep: list(seqs) for rep, seqs in clusters.items()}


def main():
    parser = argparse.ArgumentParser(
        description="Add cluster-based species and virus names to final table."
    )
    parser.add_argument('--annotated', required=True,
                        help='Final annotated TSV (metadata + coords + host taxonomy)')
    parser.add_argument('--uc', required=True,
                        help='mmseq2 .uc clustering file')
    parser.add_argument('--wg_metadata', required=True,
                        help='WG_2025 metadata TSV with columns: GenBank.Accession, species_binomial, virus_name')
    parser.add_argument('--output', required=True,
                        help='Output TSV with added Species_ICTV and virus_name columns')
    args = parser.parse_args()

    # 1. Parse clustering
    seq_to_rep, clusters = parse_clustering(args.uc)

    # 2. Read WG metadata and build annotation dict
    wg = pd.read_csv(args.wg_metadata, sep='\t')
    wg_dict = {}
    for _, row in wg.iterrows():
        acc = str(row['GenBank.Accession']).strip()
        wg_dict[acc] = {
                'species': clean_value(row.get('species_binomial', 'NA')),
                'virus_name': clean_value(row.get('virus name(10.1093/ve/veaf006)', 'NA'))
                }
    # 3. For each cluster, find the best annotation (prefer non‑NA)
    cluster_annotation = {}
    for rep, seqs in clusters.items():
        # First, try to find a member with both species and virus != 'NA'
        best_species = None
        best_virus = None
        for seq in seqs:
            if seq in wg_dict:
                sp = wg_dict[seq]['species']
                vn = wg_dict[seq]['virus_name']
                if sp != 'NA' and vn != 'NA':
                    best_species = sp
                    best_virus = vn
                    break
        # If not found, fallback to any annotated member (including 'NA')
        if best_species is None:
            for seq in seqs:
                if seq in wg_dict:
                    best_species = wg_dict[seq]['species']
                    best_virus = wg_dict[seq]['virus_name']
                    break
        if best_species is None:
            cluster_annotation[rep] = {'species': 'NA', 'virus_name': 'NA'}
        else:
            cluster_annotation[rep] = {'species': best_species, 'virus_name': best_virus}

    # 4. Build annotation for each sequence:
    #    - If sequence itself has an annotation in wg_dict, use that (overrides cluster)
    #    - Else if sequence is in clustering, use cluster annotation
    #    - Else, NA
    seq_annotation = {}
    all_seqs = set(seq_to_rep.keys()) | set(wg_dict.keys())
    for seq in all_seqs:
        
        if seq in wg_dict:
            seq_annotation[seq] = wg_dict[seq]
        elif seq in seq_to_rep:
            rep = seq_to_rep[seq]
            seq_annotation[seq] = cluster_annotation[rep]
        else:
            seq_annotation[seq] = {'species': 'NA', 'virus_name': 'NA'}
        if seq == "KJ495986":
            print(seq_annotation[seq])
    # 5. Read final annotated TSV
    final_df = pd.read_csv(args.annotated, sep='\t')

    # 6. Add new columns
    final_df['Species_ICTV'] = final_df['Accession'].map(
        lambda acc: seq_annotation.get(acc, {}).get('species', 'NA')
    )
    final_df['virus_name'] = final_df['Accession'].map(
        lambda acc: seq_annotation.get(acc, {}).get('virus_name', 'NA')
    )

    # 7. Write output
    final_df.to_csv(args.output, sep='\t', index=False)
    print(f"Updated annotation written to {args.output}")


if __name__ == '__main__':
    main()