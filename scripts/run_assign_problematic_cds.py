"""
Assign ORF types to problematic candidates using Pfam domain hits.
Calls assign_problematic_orfs from src.hmmer_utils.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.hmmer_utils import assign_problematic_orfs


def main():
    parser = argparse.ArgumentParser(
        description="Assign ORF types to problematic candidates (conflicts / no target annotation)."
    )
    parser.add_argument(
        "candidate_tsv",
        help="TSV file with problematic candidates (from *_problematic_candidates.tsv)"
    )
    parser.add_argument(
        "domtbl",
        help="hmmscan --domtblout file for problematic candidates"
    )
    parser.add_argument(
        "domain_map",
        help="CSV mapping domain names to ORF types"
    )
    parser.add_argument(
        "output_tsv",
        help="Output TSV with an additional 'assigned_orf' column"
    )
    args = parser.parse_args()

    assign_problematic_orfs(
        args.candidate_tsv,
        args.domtbl,
        args.domain_map,
        args.output_tsv
    )


if __name__ == "__main__":
    main()