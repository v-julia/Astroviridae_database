#!/usr/bin/env python3
"""
Extract ORF sequences from GenBank (or FASTA) using coordinates from a TSV file.
Generate sequence names based on metadata (host, country, date, etc.) with optional regex mappings.
"""
import argparse
import os
import re
import calendar
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.mapping_utils import load_mapping,clean_field,apply_mapping,standardize_date

def build_sequence_name(row, columns, host_compiled, country_compiled):
    """
    Build a sequence name by concatenating values from selected columns.
    Applies special handling to:
      - 'Host' -> regex mapping (if provided)
      - 'Geo location' -> regex mapping (if provided)
      - 'Collection date' -> standardize_date()
      - other columns -> clean_field()
    """
    parts = []
    for col in columns:
        val = row[col]
        if col == 'Host' and host_compiled:
            val = apply_mapping(val, host_compiled)
        elif col == 'Geo location' and country_compiled:
            val = apply_mapping(val, country_compiled)
        elif col == 'Collection date':
            val = standardize_date(val)
        else:
            val = clean_field(val)
        parts.append(val)
    return '/'.join(parts)


# ----------------------------------------------------------------------
# Main extraction function
# ----------------------------------------------------------------------
def extract_orfs(input_file, coord_tsv, meta_tsv, output_dir,
                 columns, input_format='gb', basename=None,
                 host_map=None, country_map=None, translate=False):
    """
    Extract ORF sequences using coordinates from a TSV file.
    Names are built from metadata using the provided columns.

    Parameters:
    -----------
    input_file : str
        Path to input file (GenBank or FASTA)
    coord_tsv : str
        TSV file with ORF coordinates; index: Accession, columns: 1A, 1B, 2, and strand columns
    meta_tsv : str
        Metadata TSV file (must contain all columns used in --columns)
    output_dir : str
        Directory where FASTA files will be written
    columns : str
        Comma-separated list of columns to use for building sequence names
    input_format : str
        'gb' for GenBank, 'fasta' for FASTA (default: 'gb')
    basename : str, optional
        Base name for output files (e.g., 'Astroviridae_30062026')
        If not provided, uses the stem of the input file.
    host_map : str (optional)
        TSV mapping for host names (regex_pattern, short_code)
    country_map : str (optional)
        TSV mapping for country names (regex_pattern, short_code)
    translate : bool
        If True, translate nucleotide sequences to proteins (output .faa instead of .fna)
    """
    # 1. Determine basename
    if basename is None:
        basename = os.path.splitext(os.path.basename(input_file))[0]

    # 2. Load regex mappings (if provided)
    host_compiled = load_mapping(host_map) if host_map else []
    country_compiled = load_mapping(country_map) if country_map else []

    # 3. Read coordinate file
    coords = pd.read_csv(coord_tsv, index_col=0, sep='\t')   # index is Accession
    orf_list = ['1A', '1B', '2']
    for orf in orf_list:
        if orf not in coords.columns:
            raise ValueError(f"Column '{orf}' not found in coordinate file")
        # Ensure strand columns exist; if not, assume all are +1
        if f'{orf}-strand' not in coords.columns:
            coords[f'{orf}-strand'] = 1

    # 4. Read metadata
    meta = pd.read_csv(meta_tsv, sep='\t')
    meta_dict = {row['Accession']: row for _, row in meta.iterrows()}

    # 5. Parse input file (GB or FASTA)
    if input_format.lower() not in ['gb', 'fasta']:
        raise ValueError("format must be 'gb' or 'fasta'")
    records = SeqIO.parse(input_file, input_format.lower())

    out_records = {orf: [] for orf in orf_list}
    skipped_no_coord = 0
    skipped_no_meta = 0

    for rec in records:
        acc = rec.id.split('.')[0]   # remove version suffix (e.g., ".1")
        if acc not in coords.index:
            skipped_no_coord += 1
            continue
        if acc not in meta_dict:
            skipped_no_meta += 1
            continue

        meta_row = meta_dict[acc]

        for orf in orf_list:
            coord_str = coords.loc[acc, orf]
            if coord_str == 'NA-NA' or pd.isna(coord_str):
                continue
            start, end = map(int, coord_str.split('-'))
            strand = coords.loc[acc, f'{orf}-strand']

            # Extract sequence
            if strand == 1:
                seq = rec.seq[start:end]
            else:
                seq = rec.seq[start:end].reverse_complement()

            # Optional translation
            if translate:
                seq = seq.translate(to_stop=True)

            # Build sequence name (without ORF in the name)
            name = build_sequence_name(meta_row, columns.split(','), host_compiled, country_compiled)
            # Replace '/' with '_' for safety, remove extra spaces
            name = name.replace('/', '_').strip()

            # Create SeqRecord
            new_seq = SeqRecord(seq, id=name, description='')
            out_records[orf].append(new_seq)

    # 6. Write output files
    os.makedirs(output_dir, exist_ok=True)
    suffix = 'faa' if translate else 'fna'
    for orf, seqs in out_records.items():
        if seqs:
            out_file = os.path.join(output_dir, f'{basename}_{orf}.{suffix}')
            SeqIO.write(seqs, out_file, 'fasta')
            print(f"Written {len(seqs)} sequences to {out_file}")
        else:
            print(f"No sequences for {orf}")

    if skipped_no_coord:
        print(f"Skipped {skipped_no_coord} records: accession not in coordinate file")
    if skipped_no_meta:
        print(f"Skipped {skipped_no_meta} records: accession not in metadata file")

    print("Done.")


# ----------------------------------------------------------------------
# Command-line entry point
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Extract ORFs from GenBank/FASTA using coordinates from TSV, with metadata-based naming."
    )
    parser.add_argument('--input', required=True,
                        help='Input file (GenBank or FASTA)')
    parser.add_argument('--format', default='gb', choices=['gb', 'fasta'],
                        help='Input format: gb (GenBank) or fasta (default: gb)')
    parser.add_argument('--coords', required=True,
                        help='TSV file with ORF coordinates (index: Accession, columns: 1A,1B,2,...-strand)')
    parser.add_argument('--metadata', required=True,
                        help='TSV metadata file (must contain all columns used in --columns)')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for FASTA files')
    parser.add_argument('--columns', required=True,
                        help='Comma-separated columns for sequence name (e.g., "Accession,Host_class,Host,Collection date,Country")')
    parser.add_argument('--basename', help='Base name for output files (default: stem of input file)')
    parser.add_argument('--host_map', help='TSV mapping for host names (regex, short_code) – optional')
    parser.add_argument('--country_map', help='TSV mapping for country names (regex, short_code) – optional')
    parser.add_argument('--translate', action='store_true',
                        help='Translate nucleotide sequences to proteins (output .faa instead of .fna)')
    args = parser.parse_args()

    extract_orfs(
        input_file=args.input,
        coord_tsv=args.coords,
        meta_tsv=args.metadata,
        output_dir=args.output_dir,
        columns=args.columns,
        input_format=args.format,
        basename=args.basename,
        host_map=args.host_map,
        country_map=args.country_map,
        translate=args.translate
    )


if __name__ == '__main__':
    main()