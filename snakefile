import os
import yaml
from pathlib import Path

with open("config.yaml") as f:
    config = yaml.safe_load(f)

GB = config["paths"]["genbank"]
ORF_MAP = config["paths"]["orf_map"]
DOMAIN_MAP = config["paths"]["domain_map"]
PFAM = config["paths"]["pfam_db"]
ROOT = Path(config["paths"]["output_root"])

BASE = Path(GB).stem

COORDS = ROOT / "coords"
UPDATED = ROOT / "updated_coords"
PRODIGAL = ROOT / "prodigal"
HMMER_ANNOT = ROOT / "hmmer" / f"{BASE}annotated_orfs"
HMMER_PROB = ROOT / "hmmer" / f"{BASE}problematic_orfs"
HMMER_PRED = ROOT / "hmmer" / f"{BASE}predicted_orfs"


FINAL_COORD_CSV = UPDATED / f"{BASE}_orf-coords_full.csv"


METADATA_TSV = COORDS / f"{BASE}_metadata.tsv"


# Paths for host taxonomy
HOST_TAXONOMY_DIR = ROOT / "host_taxonomy"
HOST_TAXONOMY_TSV = HOST_TAXONOMY_DIR / f"{BASE}_host_taxonomy.tsv"
HOST_MAP = Path("data/annotations/host_map.csv")


ANNOTATED_METADATA = ROOT / "final" / f"{BASE}_annotated.tsv"


rule all:
    input:
        UPDATED / f"{BASE}_orf-coords_full.csv"
        FINAL_ANNOTATED_TSV

# ----------------------------------------------------------------------
rule fetch_genbank:
    output:
        gb = GB,
    params:
        query = config["paths"]["search_query"],
        checkpoint = "data/raw/download_checkpoint.txt",
    shell:
        """
        mkdir -p $(dirname {params.checkpoint})
        python scripts/run_fetch_genbank.py "{params.query}" {output.gb} --checkpoint {params.checkpoint}
        """

# ----------------------------------------------------------------------
rule extract_metadata:
    input:
        gb = GB,
    output:
        tsv = METADATA_TSV,
    params:
        coords_dir = COORDS,
    shell:
        "python scripts/run_fetch_metadata.py {input.gb} {params.coords_dir}"

# ----------------------------------------------------------------------
rule orf_extraction:
    input:
        gb = GB,
        orf_map = ORF_MAP,
    output:
        coord_csv = COORDS / f"{BASE}_orf-coords.csv",
        orfs_faa = COORDS / f"{BASE}_orfs.faa",
        prob_faa = COORDS / f"{BASE}_problematic_candidates.faa",
        prob_tsv = COORDS / f"{BASE}_problematic_candidates.tsv",
        no_cds_fasta = COORDS / f"{BASE}_no_cds_sequences.fasta",
        log_dir = directory(COORDS / "logs"),
    params:
        coords_dir = COORDS,
    shell:
        "python scripts/run_orf_extraction.py {input.gb} {input.orf_map} {params.coords_dir}"

# ----------------------------------------------------------------------
rule hmmscan_annotated:
    input:
        fasta = COORDS / f"{BASE}_orfs.faa",
        db = PFAM,
    output:
        domtbl = HMMER_ANNOT / "annotated_domains.tbl",
        log = HMMER_ANNOT / "hmmscan.log",
    shell:
        "{config[tools][hmmscan]} --domtblout {output.domtbl} --noali {input.db} {input.fasta} > {output.log} 2>&1"

# ----------------------------------------------------------------------
rule check_mismatches:
    input:
        coord_csv = COORDS / f"{BASE}_orf-coords.csv",
        domtbl = HMMER_ANNOT / "annotated_domains.tbl",
        domain_map = DOMAIN_MAP,
    output:
        report = HMMER_ANNOT / "mismatches_report.tsv",
    shell:
        "python scripts/run_check_mismatches.py {input.coord_csv} {input.domtbl} {input.domain_map} {output.report}"

# ----------------------------------------------------------------------
rule hmmscan_problematic:
    input:
        fasta = COORDS / f"{BASE}_problematic_candidates.faa",
        db = PFAM,
    output:
        domtbl = HMMER_PROB / "problematic_domains.tbl",
        log = HMMER_PROB / "hmmscan.log",
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "Input FASTA empty. Creating empty domtbl." > {output.log}
            touch {output.domtbl}
        else
            {config[tools][hmmscan]} --domtblout {output.domtbl} --noali {input.db} {input.fasta} > {output.log} 2>&1
        fi
        """

# ----------------------------------------------------------------------
rule assign_problematic:
    input:
        cand_tsv = COORDS / f"{BASE}_problematic_candidates.tsv",
        domtbl = HMMER_PROB / "problematic_domains.tbl",
        domain_map = DOMAIN_MAP,
    output:
        assigned_tsv = HMMER_PROB / "problematic_assigned.tsv",
    shell:
        "python scripts/run_assign_problematic_cds.py {input.cand_tsv} {input.domtbl} {input.domain_map} {output.assigned_tsv}"

# ----------------------------------------------------------------------
rule update_with_problematic:
    input:
        coord_csv = COORDS / f"{BASE}_orf-coords.csv",
        assigned_tsv = HMMER_PROB / "problematic_assigned.tsv",
    output:
        updated_csv = UPDATED / f"{BASE}_orf_coords_with_problematic.csv",
    shell:
        "python scripts/run_update_coords.py {input.coord_csv} {input.assigned_tsv} {output.updated_csv}"

# ----------------------------------------------------------------------
rule prodigal:
    input:
        fasta = COORDS / f"{BASE}_no_cds_sequences.fasta",
    output:
        proteins = PRODIGAL / f"{BASE}_no_cds_sequences_proteins.faa",
        gff = PRODIGAL / f"{BASE}_no_cds_predgenes.gff",
        gbk = PRODIGAL / f"{BASE}_no_cds_predgenes.gbk",
        log = PRODIGAL / "prodigal.log",
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "Input FASTA is empty. Skipping Prodigal." > {output.log}
            touch {output.proteins} {output.gff} {output.gbk}
        else
            {config[tools][prodigal]} -i {input.fasta} -a {output.proteins} -d {output.fna} -o {output.gff} -p meta > {output.log} 2>&1
        fi
        """

# ----------------------------------------------------------------------
rule hmmscan_predicted:
    input:
        fasta = PRODIGAL / f"{BASE}_no_cds_sequences_proteins.faa",
        db = PFAM,
    output:
        domtbl = HMMER_PRED / "predicted_domains.tbl",
        log = HMMER_PRED / "hmmscan.log",
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "Input FASTA empty. Creating empty domtbl." > {output.log}
            touch {output.domtbl}
        else
            {config[tools][hmmscan]} --domtblout {output.domtbl} --noali {input.db} {input.fasta} > {output.log} 2>&1
        fi
        """

# ----------------------------------------------------------------------
rule assign_predicted:
    input:
        proteins = PRODIGAL / f"{BASE}_no_cds_sequences_proteins.faa",
        domtbl = HMMER_PRED / "predicted_domains.tbl",
        domain_map = DOMAIN_MAP,
    output:
        assigned_tsv = HMMER_PRED / "predicted_assigned.tsv",
    shell:
        "python scripts/run_assign_predicted_cds.py {input.proteins} {input.domtbl} {input.domain_map} {output.assigned_tsv}"

# ----------------------------------------------------------------------
rule update_full:
    input:
        coord_csv = UPDATED / f"{BASE}_orf_coords_with_problematic.csv",
        assigned_tsv = HMMER_PRED / "predicted_assigned.tsv",
    output:
        final_csv = UPDATED / f"{BASE}_orf-coords_full.csv",
    shell:
        "python scripts/run_update_coords.py {input.coord_csv} {input.assigned_tsv} {output.final_csv}"


rule get_host_taxonomy:
    input:
        metadata = METADATA_TSV,
        mapping = HOST_MAP,
    output:
        tax_tsv = HOST_TAXONOMY_TSV,
    params:
        email = config.get("email", "A.N.Other@example.com"),
    shell:
        "python scripts/run_get_host_taxonomy.py {input.metadata} {output.tax_tsv} --email {params.email} --mapping {input.mapping}"


rule merge_final:
    input:
        metadata = METADATA_TSV,
        coords = FINAL_COORD_CSV,
        taxonomy = HOST_TAXONOMY_TSV,   # now generated by the pipeline
    output:
        annotated = ANNOTATED_METADATA,
    shell:
        "python scripts/run_merge_metadata_coords.py {input.metadata} {input.coords} {input.taxonomy} {output.annotated}"