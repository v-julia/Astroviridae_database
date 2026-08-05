"""
Wrapper script for ORF coordinate extraction.
Calls orf_coord_updated from the src.orf_core module.
"""
import argparse
import sys
from pathlib import Path

# Allow importing from the parent src directory
sys.path.insert(0, str(Path(__file__).parent.parent))
from src.orf_core import orf_coord_updated


def main():
    parser = argparse.ArgumentParser(
        description="Extract ORF coordinates (1A, 1B, 2) from a GenBank file."
    )
    parser.add_argument(
        "input_gb",
        help="Input GenBank file (e.g., data/raw/Astroviridae_30062026.gb)"
    )
    parser.add_argument(
        "orf_map",
        help="CSV file mapping ORF names to codes (e.g., ORF1a -> 1A)"
    )
    parser.add_argument(
        "output_dir",
        help="Directory where all output files will be written (coords, logs, etc.)"
    )
    args = parser.parse_args()

    orf_coord_updated(args.input_gb, args.orf_map, args.output_dir)


if __name__ == "__main__":
    main()