""""
Utilities to merge assigned ORFs (from problematic or predicted candidates)
into the master coordinate CSV.
"""


import csv
from collections import defaultdict

def update_coords_from_assignments(coord_csv, assignments_tsv, output_csv):
    """
    Read the coordinate CSV and a TSV of assigned ORFs.
    For each accession+ORF combination in the assignments file,
    if the ORF is currently 'NA-NA' in the coordinate CSV, fill it with
    the coordinates and strand from the best assignment (longest ORF).
    If the accession is not present in the coordinate CSV, add a new row.
    All other columns (other ORFs) remain unchanged.

    Expected columns in assignments_tsv:
        accession, orf, start, end, strand
    (orf must be '1A', '1B', or '2')
    """
    # --- Read assignments and keep the longest per (accession, orf) ---
    best = {}
    with open(assignments_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Ensure required columns exist
        required = {'record_name', 'assigned_orf', 'start', 'end', 'codon_start', 'strand'}
        if not required.issubset(reader.fieldnames):
            raise ValueError(f"Assignments TSV must have columns: {required}")
        for row in reader:
            orf = row['assigned_orf']
            if orf not in ['1A', '1B', '2']:
                continue
            acc = row['record_name']
            start = int(row['start']) + int(row['codon_start'])-1
            end = int(row['end'])
            strand = int(row['strand'])
            length = end - start
            key = (acc, orf)
            if key not in best or length > best[key]['length']:
                best[key] = {
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'length': length
                }

    # --- Read existing coordinate CSV ---
    with open(coord_csv) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    # Build a set of accessions present in the CSV
    accessions_in_csv = {row['Accession'] for row in rows}
    modified_count = 0  # will count rows that were updated

    # --- Update rows for existing accessions ---
    for row in rows:
        acc = row['Accession']
        row_modified = False
        for orf in ['1A', '1B', '2']:
            if row[orf] == 'NA-NA':
                key = (acc, orf)
                if key in best:
                    c = best[key]
                    row[orf] = f"{c['start']}-{c['end']}"
                    row[orf + '-strand'] = str(c['strand'])
                    row_modified = True
        if row_modified:
            modified_count += 1
    # --- Add rows for new accessions that appear in assignments but not in CSV ---
    # First, collect all accessions from assignments
    assigned_accessions = {acc for (acc, _) in best.keys()}
    new_accessions = assigned_accessions - accessions_in_csv
    print("New accessions:")
    print(new_accessions)
    for acc in new_accessions:
        new_row = {'Accession': acc}
        # Initialise all ORFs to NA-NA
        for orf in ['1A', '1B', '2']:
            new_row[orf] = 'NA-NA'
            new_row[orf + '-strand'] = 'NA'
        # Fill in the assigned ones
        for orf in ['1A', '1B', '2']:
            key = (acc, orf)
            if key in best:
                c = best[key]
                new_row[orf] = f"{c['start']}-{c['end']}"
                new_row[orf + '-strand'] = str(c['strand'])
        rows.append(new_row)

    # --- Write updated CSV ---
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Updated coordinate CSV written to: {output_csv}")
    print(f"  - Modified rows for existing accessions: {modified_count}")
    print(f"  - Added rows for new accessions: {len(new_accessions)}")