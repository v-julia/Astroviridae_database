#!/usr/bin/env python3
"""
Extract ORF sequences from GenBank (or FASTA) using coordinates from a CSV file.
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


# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def load_mapping(mapping_file):
    """
    Load a TSV mapping file with columns: regex_pattern, short_code.
    Returns a list of (compiled_regex, short_code) tuples.
    """
    if not mapping_file or not os.path.exists(mapping_file):
        return []
    df = pd.read_csv(mapping_file, sep='\t', header=None, names=['pattern', 'code'])
    compiled = []
    for _, row in df.iterrows():
        try:
            pat = re.compile(str(row['pattern']), re.IGNORECASE)
            compiled.append((pat, row['code']))
        except re.error as e:
            print(f"Warning: invalid regex '{row['pattern']}': {e}")
    return compiled


def apply_mapping(value, compiled_patterns, default='NA'):
    """
    Apply regex patterns to a value and return the first matching short code.
    If no pattern matches, return a cleaned version of the original value.
    """
    if pd.isna(value):
        return default
    for pat, code in compiled_patterns:
        if pat.search(str(value)):
            return code
    # Fallback: clean the original value
    return str(value).replace(' ', '-').replace('/', '-')


def clean_field(value):
    """Replace problematic characters for safe use in filenames."""
    if pd.isna(value):
        return 'NA'
    return str(value).replace('/', '-').replace(' ', '-').replace(',', '-').replace('(', '-').replace(')', '-')


def standardize_date(date_str):
    """
    Convert various date formats to ISO 8601-like string:
      - YYYY-MM-DD if day available
      - YYYY-MM if only month available
      - YYYY if only year available
    Handles ranges (e.g., 08-Jul-2014/11-Oct-2016) by joining with '_'.
    """
    if pd.isna(date_str):
        return 'NA'
    s = str(date_str).strip()

    # Handle date range: split on '/' and standardize each part, then join with '_'
    if '/' in s and not re.match(r'^\d{4}$', s) and not re.match(r'^\d{4}-\d{2}$', s):
        parts = s.split('/')
        std_parts = [standardize_date(p) for p in parts]
        return '_'.join(std_parts)

    # DD-MMM-YYYY (e.g., 13-Jun-2024)
    match = re.match(r'^(\d{1,2})-([A-Za-z]{3})-(\d{4})$', s)
    if match:
        day, mon, year = match.groups()
        month_num = list(calendar.month_abbr).index(mon.capitalize())
        return f"{year}-{month_num:02d}-{int(day):02d}"

    # DD.MM.YYYY (e.g., 02.10.2019)
    match = re.match(r'^(\d{1,2})\.(\d{1,2})\.(\d{4})$', s)
    if match:
        day, mon, year = match.groups()
        return f"{year}-{int(mon):02d}-{int(day):02d}"

    # YYYY-MM (e.g., 2019-10)
    match = re.match(r'^(\d{4})-(\d{1,2})$', s)
    if match:
        year, mon = match.groups()
        return f"{year}-{int(mon):02d}"

    # MM-YYYY or MM.YYYY (e.g., 11-2025 or 11.2025)
    match = re.match(r'^(\d{1,2})[-.](\d{4})$', s)
    if match:
        mon, year = match.groups()
        return f"{year}-{int(mon):02d}"

    # MMM-YYYY (e.g., Nov-2025)
    match = re.match(r'^([A-Za-z]{3})-(\d{4})$', s)
    if match:
        mon, year = match.groups()
        month_num = list(calendar.month_abbr).index(mon.capitalize())
        return f"{year}-{month_num:02d}"

    # Just YYYY
    match = re.match(r'^(\d{4})$', s)
    if match:
        return s

    # If nothing matches, clean and return as is
    return clean_field(s)


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
def extract_orfs(input_file, coord_csv, meta_tsv, output_dir,
                 columns, input_format='gb', basename=None,
                 host_map=None, country_map=None, translate=False):
    """
    Extract ORF sequences using coordinates from a CSV file.
    Names are built from metadata using the provided columns.

    Parameters:
    -----------
    input_file : str
        Path to input file (GenBank or FASTA)
    coord_csv : str
        CSV file with ORF coordinates; index: Accession, columns: 1A, 1B, 2, and strand columns
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
    coords = pd.read_csv(coord_csv, index_col=0)   # index is Accession
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
            if coord_str == 'NA-NA':
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
        description="Extract ORFs from GenBank/FASTA using coordinates from CSV, with metadata-based naming."
    )
    parser.add_argument('--input', required=True,
                        help='Input file (GenBank or FASTA)')
    parser.add_argument('--format', default='gb', choices=['gb', 'fasta'],
                        help='Input format: gb (GenBank) or fasta (default: gb)')
    parser.add_argument('--coords', required=True,
                        help='CSV file with ORF coordinates (index: Accession, columns: 1A,1B,2,...-strand)')
    parser.add_argument('--metadata', required=True,
                        help='TSV metadata file (must contain all columns used in --columns)')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for FASTA files')
    parser.add_argument('--columns', required=True,
                        help='Comma-separated columns for sequence name (e.g., Host_class,Host,Accession,Collection_date,Country)')
    parser.add_argument('--basename', help='Base name for output files (default: stem of input file)')
    parser.add_argument('--host_map', help='TSV mapping for host names (regex, short_code) – optional')
    parser.add_argument('--country_map', help='TSV mapping for country names (regex, short_code) – optional')
    parser.add_argument('--translate', action='store_true',
                        help='Translate nucleotide sequences to proteins (output .faa instead of .fna)')
    args = parser.parse_args()

    extract_orfs(
        input_file=args.input,
        coord_csv=args.coords,
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