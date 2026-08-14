"""
Merge metadata, final coordinates, and host taxonomy into one annotated table.
"""
import argparse
import sys
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

def main():
    parser = argparse.ArgumentParser(
        description="Merge metadata, final coordinates, and host taxonomy."
    )
    parser.add_argument("metadata_tsv", help="Metadata TSV (from fetch_metadata_from_gb)")
    parser.add_argument("coord_csv", help="Final coordinate CSV (e.g., *_orf-coords_full.csv)")
    parser.add_argument("host_taxonomy_tsv", help="Host taxonomy mapping TSV (from run_get_host_taxonomy)")
    parser.add_argument("output_tsv", help="Output merged TSV")
    args = parser.parse_args()

    # Read inputs
    meta = pd.read_csv(args.metadata_tsv, sep='\t')
    coords = pd.read_csv(args.coord_csv, sep='\t')
    tax = pd.read_csv(args.host_taxonomy_tsv, sep='\t')

    # Merge metadata with coordinates on Accession
    merged = meta.merge(coords, left_on="Accession", right_on="Accession", how="left")

    # Merge with taxonomy on Host (original, not cleaned)
    merged = merged.merge(tax, left_on="Host", right_on="host_original", how="left")
    # Drop the redundant 'host_original' column from tax
    merged.drop(columns=["host_original"], inplace=True, errors='ignore')

    taxonomy_cols = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    for col in taxonomy_cols:
        if col in merged.columns:
            merged.rename(columns={col: f"Host_{col}"}, inplace=True)

    # Reorder columns: Accession, metadata, then coordinates, then taxonomy
    # (optional)
    merged.to_csv(args.output_tsv, sep='\t', index=False)
    print(f"Merged file written to {args.output_tsv}")

if __name__ == "__main__":
    main()