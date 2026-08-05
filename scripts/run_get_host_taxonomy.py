"""
Extract host taxonomy (lineage ranks) from NCBI for unique hosts in metadata.
Uses a manually curated mapping file to handle common names and ambiguous cases.
"""
import argparse
import sys
import time
import re
import pandas as pd
from pathlib import Path
from Bio import Entrez

# ---------------------------
# Manual mapping for common names 
# ---------------------------

def clean_host_name(host):
    """
    Clean a host name to make it NCBI‑queryable.
    Returns (cleaned_name, has_valid_name).
    """
    if pd.isna(host) or not str(host).strip():
        return None, False

    original = str(host).strip()

    # 1. Remove text after semicolon or comma (often sex/age/lab info)
    #    Take the first part before ; or ,
    cleaned = re.split(r'[;,]\s*', original)[0].strip()

    # 2. Remove parenthesised content
    cleaned = re.sub(r'\([^)]*\)', '', cleaned).strip()

    # 3. Remove trailing ' sp.' or ' sp' (with a space)
    cleaned = re.sub(r'\s+sp\.?\s*$', '', cleaned).strip()

    # 4. Remove age/sex descriptors and everything after them
    #    We remove from the first occurrence of common descriptors
    cleaned = re.sub(
        r'\s+(sex|age|year|yr|y|month|mo|old|child|boy|girl|male|female|fetus|infant|baby|newborn|adult|breed|genotype|line|strain)\s*.*$',
        '',
        cleaned,
        flags=re.IGNORECASE
    ).strip()

    # 5. Remove descriptors like "breed", "genotype", "line", "strain"
    cleaned = re.sub(r'(?i)\b(breed|genotype|line|strain)\s*.*$', '', cleaned).strip()

    # 6. Remove trailing lab codes (e.g., "FMNH 228875")
    cleaned = re.sub(r'\s+\b[A-Z]*\d+\b$', '', cleaned).strip()

    # 7. Remove trailing numbers (e.g., "Homo sapiens 4" -> "Homo sapiens")
    cleaned = re.sub(r'\s+\d+$', '', cleaned).strip()

    # 8. If empty or too short
    if len(cleaned) < 2:
        print(cleaned)
        return None, False

    return cleaned, True

def get_lineage_with_ranks(organism_name):
    """Query NCBI taxonomy and return a dict of ranks -> scientific names."""
    try:
        stream = Entrez.esearch(db="taxonomy", term=organism_name, retmode="xml")
        record = Entrez.read(stream)
        stream.close()
        idlist = record.get("IdList", [])
        if not idlist:
            # Try with a broader search: use the name as a title?
            # We'll just return empty.
            return {}
        taxid = idlist[0]
        stream = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(stream)
        stream.close()
        lineage_list = records[0].get("LineageEx", [])
        lineage_dict = {
            item["Rank"]: item["ScientificName"]
            for item in lineage_list
            if item["Rank"] != "no rank"
        }
        lineage_dict[records[0]["Rank"]] = records[0]["ScientificName"]
        return lineage_dict
    except Exception as e:
        print(f"Error querying {organism_name}: {e}", file=sys.stderr)
        return {}


def load_mapping(mapping_file):
    """Load a CSV with columns 'raw' and 'cleaned'."""
    if not mapping_file or not Path(mapping_file).exists():
        return {}
    df = pd.read_csv(mapping_file)
    # Convert to dict, lowercasing raw for case‑insensitive matching
    mapping = {}
    for _, row in df.iterrows():
        raw = str(row['raw']).strip().lower()
        cleaned = str(row['cleaned']).strip()
        mapping[raw] = cleaned
    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Get NCBI taxonomy ranks for each unique host in metadata."
    )
    parser.add_argument("metadata_tsv", help="Input metadata TSV (must have 'Host' column)")
    parser.add_argument("output_tsv", help="Output TSV with columns: host_original, host_clean, class, family, ...")
    parser.add_argument("--email", default="A.N.Other@example.com", help="NCBI email")
    parser.add_argument("--mapping", help="Optional CSV mapping file (raw, cleaned)")
    args = parser.parse_args()

    Entrez.email = args.email

    # Load mapping
    mapping = load_mapping(args.mapping)

    # Read metadata
    df = pd.read_csv(args.metadata_tsv, sep='\t')
    raw_hosts = df["Host"].dropna().unique().tolist()
    # Create a mapping from raw host to cleaned name
    host_to_clean = {}
    for raw in raw_hosts:
        # 1. Check raw string in mapping (case‑insensitive)
        raw_lower = raw.lower()
        if raw_lower in mapping:
            cleaned = mapping[raw_lower]
            host_to_clean[raw] = cleaned
        else:
            # 2. Clean the raw host
            cleaned, ok = clean_host_name(raw)
            if not ok:
                host_to_clean[raw] = None
                print(f"Could not clean {raw}")
                continue
            # 3. Check cleaned name in mapping (lowercase)
            if cleaned.lower() in mapping:
                cleaned = mapping[cleaned.lower()]
        host_to_clean[raw] = cleaned
    # Collect unique cleaned names
    unique_cleaned = set()
    for raw, cleaned in host_to_clean.items():
        if cleaned is not None:
            unique_cleaned.add(cleaned)
    # Query each unique cleaned name once
    lineage_cache = {}
    total = len(unique_cleaned)
    print(f"Querying {total} unique host names from NCBI...", file=sys.stderr)
    for i, cleaned_name in enumerate(sorted(unique_cleaned), 1):
        print(f"{i}/{total}: {cleaned_name}", file=sys.stderr)
        lineage = get_lineage_with_ranks(cleaned_name)
        lineage_cache[cleaned_name] = lineage
        time.sleep(1)

    # Build output rows
    out_rows = []
    for raw_host in raw_hosts:
        cleaned = host_to_clean.get(raw_host)
        if cleaned is None:
            row = {"host_original": raw_host, "host_clean": None}
            for rank in ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']:
                row[rank] = 'NA'
            out_rows.append(row)
            continue
        lineage = lineage_cache.get(cleaned, {})
        row = {"host_original": raw_host, "host_clean": cleaned}
        for rank in ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']:
            row[rank] = lineage.get(rank, 'NA')
        out_rows.append(row)

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.output_tsv, sep='\t', index=False)

    # Print a summary
    total_orig = len(raw_hosts)
    succeeded = sum(1 for r in out_rows if r['host_clean'] is not None and r['species'] != 'NA')
    print(f"Processed {total_orig} unique hosts.", file=sys.stderr)
    print(f"Successfully mapped/queried: {succeeded}", file=sys.stderr)
    print(f"Failed (no taxonomy): {total_orig - succeeded}", file=sys.stderr)

if __name__ == "__main__":
    main()