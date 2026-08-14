"""
Check for discrepancies between annotated ORF types and Pfam domain hits.
Calls check_annotated_orfs from src.hmmer_utils.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.hmmer_utils import check_annotated_orfs


def main():
    parser = argparse.ArgumentParser(
        description="Check mismatches between annotated ORFs and Pfam domain assignments."
    )
    parser.add_argument(
        "coord_tsv",
        help="Coordinate TSV file (e.g., *_orf-coords.tsv)"
    )
    parser.add_argument(
        "domtbl",
        help="hmmscan --domtblout file for annotated ORFs"
    )
    parser.add_argument(
        "domain_map",
        help="CSV mapping domain names to ORF types (domain_name, orf_type)"
    )
    parser.add_argument(
        "out_report",
        help="Output TSV report of mismatches (or 'No mismatches found.')"
    )
    args = parser.parse_args()

    mismatches = check_annotated_orfs(args.coord_tsv, args.domtbl, args.domain_map)

    with open(args.out_report, 'w') as f:
        if mismatches:
            f.write("accession\tannotated_orf\tdomain_orf\tstart\tend\tstrand\n")
            for acc, orf, assigned, start, end, strand in mismatches:
                f.write(f"{acc}\t{orf}\t{assigned}\t{start}\t{end}\t{strand}\n")
        else:
            f.write("No mismatches found.\n")


if __name__ == "__main__":
    main()