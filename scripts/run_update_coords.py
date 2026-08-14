"""
Merge assigned ORFs (from problematic or predicted candidates) into the master
coordinate CSV.
Calls update_coords_from_assignments from src.update_utils.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.update_utils import update_coords_from_assignments


def main():
    parser = argparse.ArgumentParser(
        description="Update coordinate TSV with assigned ORFs from a TSV file."
    )
    parser.add_argument(
        "coord_tsv",
        help="Input coordinate TSV (e.g., *_orf-coords.tsv)"
    )
    parser.add_argument(
        "assignments_tsv",
        help="TSV file with assignments (columns: record_name, assigned_orf, start, end, strand, codon_start)"
    )
    parser.add_argument(
        "output_csv",
        help="Output CSV with updated coordinates"
    )
    args = parser.parse_args()

    update_coords_from_assignments(
        args.coord_tsv,
        args.assignments_tsv,
        args.output_csv
    )


if __name__ == "__main__":
    main()