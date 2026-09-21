"""
Extract metadata (accession, taxonomy, references, features, etc.)
from a GenBank file and write a TSV file.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.gb_utils import fetch_metadata_from_gb


def main():
    parser = argparse.ArgumentParser(
        description="Extract metadata from a GenBank file and save as TSV."
    )
    parser.add_argument(
        "input_gb",
        help="Input GenBank file"
    )
    parser.add_argument(
        "output_dir",
        help="Output directory for the TSV file (file name will be based on input)"
    )
    parser.add_argument(
        "country_map",
        help="TSV file with correspondense of countries and their 3 letter codes"
    )
    args = parser.parse_args()

    fetch_metadata_from_gb(args.input_gb, args.output_dir, args.country_map)


if __name__ == "__main__":
    main()