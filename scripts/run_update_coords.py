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
        description="Update coordinate CSV with assigned ORFs from a TSV file."
    )
    parser.add_argument(
        "coord_csv",
        help="Input coordinate CSV (e.g., *_orf-coords.csv)"
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
        args.coord_csv,
        args.assignments_tsv,
        args.output_csv
    )


if __name__ == "__main__":
    main()