"""
Correct ORF annotations in *_orf-coords.tsv using mismatches report.

For each (accession, annotated_orf) mismatch:
  * if domain_orf is a normal ORF (1A, 1B, 2, ...) -> move coords into correct column
  * if domain_orf is a fusion (contains "_")       -> clear column + write to fusion table

Fusion accessions are also exported as a list so their nucleotide sequences
can be appended to the *_no_cds_sequences.fasta for Prodigal prediction.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

#ORF_NAMES = ["1A", "1B", "2"]      # extend if your project has more
COORD_NA        = "NA-NA"
STRAND_NA = "NA"


def load_mismatches(path):
    """Read mismatches report. Returns empty DataFrame if 'No mismatches found.'"""
    with open(path) as f:
        first = f.readline()
    if first.startswith("No mismatches found"):
        return pd.DataFrame(columns=["accession", "annotated_orf", "domain_orf",
                                     "start", "end", "strand"])
    return pd.read_csv(path, sep="\t")


def build_corrections(mismatch_df):
    """
    Return  {accession: [ {wrong_orf, correct_orf, coords, strand}, ... ]}
    """
    per_acc = {}
    for _, r in mismatch_df.iterrows():
        per_acc.setdefault(r["accession"], []).append({
            "wrong_orf":   r["annotated_orf"],
            "correct_orf": r["domain_orf"],
            "coords":      f"{r['start']}-{r['end']}",
            "strand":      r["strand"],
        })
    return per_acc

def apply_corrections(coord_df, corrections):
    """
    Returns (corrected_coord_df, fusion_df).

    corrections: {accession: [fix, ...]} where fix =
        {"wrong_orf", "correct_orf", "coords", "strand"}
    """
    coord_df = coord_df.set_index("Accession", drop=False)
    n_corrected = 0
    fusion_rows = []

    for acc, fixes in corrections.items():
        if acc not in coord_df.index:
            continue
        
        sources, targets = set(), set()

        for fix in fixes:
            wrong   = fix["wrong_orf"]
            correct = fix["correct_orf"]

            # Fusion -> record, clear source
            if "_" in str(correct):
                fusion_rows.append({
                    "Accession":   acc,
                    "wrong_orf":   wrong,
                    "correct_orf": correct,
                    "coordinates": fix["coords"],
                    "strand":      fix["strand"],
                })
                coord_df.loc[acc, wrong]             = COORD_NA
                coord_df.loc[acc, f"{wrong}-strand"] = STRAND_NA
                sources.add(wrong)
                continue

            # Normal ORF
            coord_df.loc[acc, correct]             = fix["coords"]
            coord_df.loc[acc, f"{correct}-strand"] = fix["strand"]
            sources.add(wrong)
            targets.add(correct)
            n_corrected += 1

        for src in sources - targets:
            coord_df.loc[acc, src]             = COORD_NA
            coord_df.loc[acc, f"{src}-strand"] = STRAND_NA

    coord_df = coord_df.reset_index(drop=True)
    fusion_df = pd.DataFrame(
        fusion_rows,
        columns=["Accession", "wrong_orf", "correct_orf",
                 "coordinates", "strand"],
    )
    return coord_df, fusion_df, n_corrected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coord_tsv",   help="*_orf-coords.tsv (original)")
    parser.add_argument("mismatch_tsv",help="*_mismatches_report.tsv")
    parser.add_argument("out_corrected",help="*_orf-coords_corrected.tsv")
    parser.add_argument("out_fusion",  help="*_fusion_orfs.tsv")
    args = parser.parse_args()

    coord_df    = pd.read_csv(args.coord_tsv, sep="\t")
    mismatch_df = load_mismatches(args.mismatch_tsv)
    corrections = build_corrections(mismatch_df)

    
    coord_df_corrected, fusion_df, n_corrected = apply_corrections(coord_df, corrections)
    
    coord_df_corrected.to_csv(args.out_corrected, sep="\t", index=False,
                                                            na_rep="NA")

    fusion_df.to_csv(args.out_fusion, sep="\t", index=False, na_rep="NA")

    print(f"Corrected accessions : {n_corrected}")
    print(f"Entries with fused ORFs: {len(fusion_df)}")


if __name__ == "__main__":
    main()
