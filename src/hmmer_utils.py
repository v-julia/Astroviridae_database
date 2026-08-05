"""
HMMER results processing: parse domain tables, assign ORF types based on domain hits,
check mismatches between annotation and domain assignment, and handle
problematic / predicted CDS candidates.
"""

import csv
import re

def load_domain_map(csv_file):
    """
    Load domain-to-ORF mapping from a CSV file.
    Expected columns: domain_name, orf_type
    Returns a dict {domain_name: orf_type}.
    """
    domain_map = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            domain = row['domain_name'].strip()
            orf = row['orf_type'].strip()
            if domain and orf:
                domain_map[domain] = orf
    return domain_map

def parse_hmmscan_domtbl(domtbl_file, evalue_threshold=0.01):
    """
    Parse hmmscan --domtblout output.
    Returns a dict: {query_id: [list_of_domain_hits]}
    Each hit is a dict with keys:
        'domain_name', 'domain_accession', 'domain_description',
        'query_start', 'query_end', 'score', 'evalue', 'bitscore'
    """
    domain_hits = {}
    with open(domtbl_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 23:
                continue
            # Fields (0‑based indices):
            # 0: target name (domain)
            # 1: target accession
            # 2: tlen
            # 3: query name (your sequence ID)
            # 4: query accession (often '-')
            # 5: qlen
            # 6: evalue (full sequence)
            # 7: score (full sequence)
            # 8: bias
            # 9: domain number
            # 10: domain count in this hit
            # 11: domain evalue
            # 12: domain score
            # 13: domain bias
            # 14: hmm coord from
            # 15: hmm coord to
            # 16: ali coord from (query start)
            # 17: ali coord to (query end)
            # 18: env coord from
            # 19: env coord to
            # 20: acc (optional)
            # 21: description of target (the rest)
            query = fields[3]
            domain_name = fields[0]
            domain_acc = fields[1]
            domain_evalue = float(fields[11])
            if domain_evalue > evalue_threshold:
                continue
            query_start = int(fields[16])
            query_end = int(fields[17])
            bitscore = float(fields[13])
            # Full sequence description (may contain spaces)
            desc = ' '.join(fields[22:]) if len(fields) > 22 else ''
            
            hit = {
                'domain_name': domain_name,
                'domain_accession': domain_acc,
                'domain_description': desc,
                'query_start': query_start,
                'query_end': query_end,
                'score': bitscore,
                'evalue': domain_evalue
            }
            domain_hits.setdefault(query, []).append(hit)
    return domain_hits
def assign_orf_from_domains(domain_hits, domain_map, min_score=0, min_coverage=0.5):
    """
    Given a list of domain hits for a single sequence,
    return the most likely ORF type based on the best hit (lowest E-value).
    If no hit matches the map, return None.
    """
    if not domain_hits:
        return None
    # Sort by E-value (ascending)
    sorted_hits = sorted(domain_hits, key=lambda x: x['evalue'])
    for hit in sorted_hits:
        dom_name = hit['domain_name']
        if dom_name in domain_map:
            return domain_map[dom_name]
    return None

def check_annotated_orfs(coord_csv, domain_table_file, domain_map_file):
    # Parse domain hits
    domain_hits = parse_hmmscan_domtbl(domain_table_file)
    domain_map = load_domain_map(domain_map_file)
    mismatches = []
    with open(coord_csv) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    for row in rows:
        acc = row['Accession']
        for orf in ['1A', '1B', '2']:
            if row[orf] == 'NA-NA':
                continue
            start, end = row[orf].split('-')
            strand = row[orf + '-strand']
            # Build the exact header used in FASTA
            header = f"{acc}|{orf}|{start}-{end}|strand={strand}"  # no len, but you can include if needed
            # Find domain hits for this header
            hits = domain_hits.get(header, [])
            assigned = assign_orf_from_domains(hits, domain_map)
            if assigned and assigned != orf:
                mismatches.append((acc, orf, assigned, start, end, strand))
    return mismatches


def assign_problematic_orfs(candidate_tsv, domain_table_file, domain_map_file, output_tsv):
    """
    For each problematic candidate, use the domain hits to assign an ORF type.
    Writes a new TSV with an additional 'assigned_orf' column.
    """
    domain_map = load_domain_map(domain_map_file)
    domain_hits = parse_hmmscan_domtbl(domain_table_file)
    with open(candidate_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    not_assigned_rows = []
    assigned_rows = []
    for row in rows:
        # Reconstruct the FASTA header from the candidate info
        cid = row['candidate_id']
        header = f"{cid}|{row['accession']}|{row['gene']}|start={row['start']}|end={row['end']}|codon_start={row['codon_start']}|strand={row['strand']}"
        hits = domain_hits.get(header, [])
        orf_type = assign_orf_from_domains(hits, domain_map)
        row['assigned_orf'] = orf_type if orf_type else 'unknown'
        assigned_rows.append(row)
        if not orf_type:
            not_assigned_rows.append(row)
        
    print(f"{len(not_assigned_rows)} CDS were not assigned")
    
    # Write updated TSV
    with open(output_tsv, 'w', newline='') as f:
        fieldnames = reader.fieldnames + ['assigned_orf']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(assigned_rows)

def parse_prodigal_fasta_headers(fasta_file):
    """
    Parse Prodigal protein FASTA headers to extract coordinates.
    Returns a dict: {query_id: {'accession': acc, 'start': start, 'end': end, 'strand': strand}}
    start: 0-based, end: exclusive, strand: 1 or -1.
    
    Expected header format:
        >PQ110289.1_1 # 2 # 886 # 1 # ID=1_1;partial=10;...
    """
    coords = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                # Split by '#' to extract fields
                parts = [p.strip() for p in header.split('#')]
                if len(parts) >= 4:
                    query_id = parts[0]          # e.g., "PQ110289.1_1"
                    # Extract accession (remove the trailing gene ID)

                    id_parts = query_id.split('_')
                    if len(id_parts) >= 2:
                        gene_num = id_parts.pop()  # the gene number (e.g., "1")
                        accession = '_'.join(id_parts)
                    else:
                        accession = query_id

                    start = int(parts[1])        # 1-based start
                    end = int(parts[2])          # 1-based end (inclusive)
                    strand = int(parts[3])       # 1 or -1
                    # Convert to 0-based start and exclusive end (to match your CSV)
                    start0 = start - 1
                    end_excl = end + 1
                    coords[query_id] = {
                        'accession': accession,
                        'record_name': accession.split('.')[0],
                        'start': start0,
                        'end': end_excl,
                        'strand': strand,
                    }
    return coords


def assign_predicted_orfs(protein_fasta, domain_table_file, domain_map_file, output_tsv):
    """
    Parse Prodigal protein FASTA headers, combine with hmmscan domain hits,
    and assign ORF types to each predicted CDS.
    Writes a TSV with columns: accession, assigned_orf, start, end, strand, codon_start (1).
    """
    # Load domain map
    domain_map = load_domain_map(domain_map_file)
    # Parse domain hits
    domain_hits = parse_hmmscan_domtbl(domain_table_file)
    # Parse Prodigal FASTA headers
    coord_map = parse_prodigal_fasta_headers(protein_fasta)

    assigned_rows = []
    not_assigned = 0
    for query_id, coords in coord_map.items():
        hits = domain_hits.get(query_id, [])
        orf_type = assign_orf_from_domains(hits, domain_map)
        if not orf_type:
            orf_type = 'unknown'
            not_assigned += 1
        row = {
            'record_name': coords['record_name'],
            'accession': coords['accession'],
            'assigned_orf': orf_type,
            'start': coords['start'],
            'end': coords['end'],
            'strand': coords['strand'],
            'codon_start': 1,          # Prodigal uses 1-based with codon_start=1 for its predictions
        }
        assigned_rows.append(row)

    print(f"Total predicted CDS: {len(coord_map)}")
    print(f"Assigned: {len(assigned_rows) - not_assigned}")
    print(f"Unassigned: {not_assigned}")

    # Write TSV
    with open(output_tsv, 'w', newline='') as f:
        fieldnames = ['record_name','accession', 'assigned_orf', 'start', 'end', 'strand', 'codon_start']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(assigned_rows)