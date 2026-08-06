# Astroviridae Database

This repository contains a **Snakemake pipeline** for annotating astrovirus CDSs from GenBank entries, along with a curated dataset of astrovirus and astro‑like sequences deposited in GenBank.

The pipeline automates:

- ORF coordinate extraction (ORF1a, ORF1b, ORF2)
- Frame validation and correction
- Verification of CDS ORF assignments using Pfam‑A domain searches
- Retrieval of host lineage information from NCBI Taxonomy
- Propagation of ICTV species designations

---

## Dataset

`data/results/final/Astroviridae_30062026_annotated_MAstV_species.tsv` is a curated, tab‑separated TSV file containing metadata extracted from GenBank entries, coordinates of astrovirus ORFs, host lineage information, and ICTV species assignments for classified mamastroviruses.

The file can be imported into Excel, R, Python (pandas), or any other spreadsheet software.

- **16,448 sequences** downloaded on **30 June 2026**
- Includes astroviruses and astro‑like viruses (bastroviruses, picorna‑like sequences, and other entries occasionally assigned to Astroviridae by GenBank staff)

### Dataset columns

| Column | Description |
|--------|-------------|
| `Accession` | GenBank accession (without version suffix) |
| `Version` | Full versioned accession (e.g., `NC_123456.1`) |
| `GenBank title` | Definition line from GenBank |
| `Release date` | Date the record was released |
| `Organism name` | Organism scientific name |
| `Virus Lineage` | Full taxonomic lineage |
| `Species` | Virus species (from NCBI taxonomy) |
| `Family` | Virus family |
| `Length` | Genome length (bp) |
| `Reference_title`, `Submitters`, `PUBMED_ID` | Reference information |
| `Strain`, `Isolate` | Strain/isolate names |
| `Geo location` | Geographic location |
| `Collection date` | Collection date |
| `Tissue/Specimen/Source` | Isolation source |
| `Host` | Original host name |
| `Environmental` | `yes` if environmental sample |
| `CDS count` | Number of annotated CDSs |
| `Artificial` | `yes` if sequence is artificial (patent, nonfunctional, etc.) |
| `UNVERIFIED` | `yes` if sequence is marked UNVERIFIED |
| `1A`, `1B`, `2` | ORF coordinates (`NA‑NA` if absent) |
| `1A‑strand`, `1B‑strand`, `2‑strand` | Strand orientation (`1` or `-1`) |
| `Host_species`, `Host_genus`, `Host_family`, `Host_order`, `Host_class`, `Host_phylum`, `Host_kingdom` | Full host lineage (from NCBI Taxonomy) |
| `Species_ICTV` | ICTV species designation (for classified mamastroviruses) |
| `virus_name` | Virus name (from 10.1093/ve/veaf006) |

All ORF coordinates are **0‑based and exclusive**, consistent with Biopython's `feature.location` representation.

---

## Pipeline Overview

The pipeline is implemented as a **Snakemake workflow** and consists of the following steps:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `run_fetch_genbank.py` | Download GenBank entries from NCBI using a user‑defined query |
| 2 | `run_fetch_metadata.py` | Extract all metadata (taxonomy, references, features) from the GenBank file |
| 3 | `run_orf_extraction.py` | Extract ORF1a, ORF1b, and ORF2 coordinates; handle conflicting/problematic annotations |
| 4 | `orf_core.py` (internal) | Extracts ORFs coordinates, validates and corrects reading frames; identifies CDSs with ambiguous annotations and genomes lacking annotated CDSs |
| 5 | `run_hmmscan.py` | Run HMMER against Pfam‑A on annotated ORFs, CDSs with ambigous annotation and predicted ORFs |
| 6 | `run_prodigal.py` | Predict ORFs with Prodigal‑gv for sequences lacking CDS annotation |
| 7 | `run_assign_problematic.py` / `run_assign_predicted.py` | Assign ORF types (ORF1A/ORF1B/ORF2) using domain‑to‑ORF mapping |
| 8 | `run_update_coords.py` | Merge original coordinates, problematic assignments, and predicted ORFs into one table |
| 9 | `run_get_host_taxonomy.py` | Query NCBI Taxonomy for each host and retrieve full lineage |
| 10 | `run_merge_all.py` | Combine metadata, coordinates, and host taxonomy into a single TSV |
| 11 | `add_cluster_species.py` | Propagate ICTV species and virus names from reference clusters using mmseq2 |
| 12 | `run_extract_orfs_named.py` | Extract individual ORF FASTA files with metadata‑derived headers |

---

## Auxiliary Scripts

The following scripts are **not** required for the main pipeline but provide additional functionality:

| Script | Description |
|--------|-------------|
| `run_extract_orfs_named.py` | Extract individual ORF FASTA files with metadata‑derived headers (host class, host, accession, collection date, country). Useful for downstream phylogenetic analyses. |


## Host lineages

Host lineage is retrieved from the NCBI Taxonomy database. The dataset includes sequences from the following host classes:

- **Mammalia** – mammals
- **Aves** – birds
- **Lepidosauria** – lizards and snakes
- **Amphibia** – amphibians
- **Actinopteri**, **Chondrichthyes**, **Hyperoartia** – fish
- **Arachnida**, **Insecta**, **Malacostraca** – arthropods
- **Clitellata** – annelid worms
- **Bivalvia** – mollusks
- **Magnoliopsida** – flowering plants

---

## Species Assignment

The `Species_ICTV` column contains ICTV‑recognised species designations propagated from a reference dataset of complete genomes of the genus *Mamastrovirus* using mmseq2 clustering (17% p-distance in nucleotide ORF1b). 
Species are assigned to all sequences that cluster with a reference classified sequence. This column is being expanded as the *Astroviridae* Study Group progresses.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/v-julia/Astroviridae_database.git
cd Astroviridae_database
```

### 2. Create a Conda environment

```bash
conda install -c bioconda snakemake biopython pandas pyyaml hmmer prodigal-gv
```

### 3. Prepare the Pfam database

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

Update `pfam_db` in `config.yaml` to point to `Pfam-A.hmm`.

---

## Configuration (`config.yaml`)

```yaml
paths:
  genbank: "data/raw/Astroviridae_30062026.gb"
  search_query: '"txid39733"[Organism]'
  orf_map: "data/annotations/ORF_names.csv"
  domain_map: "data/annotations/domain_map.csv"
  pfam_db: "/path/to/Pfam-A.hmm"
  output_root: "data/results"
tools:
  prodigal: "prodigal-gv"
  hmmscan: "hmmscan"
email: "your.email@example.com"
```

---

## Running the Pipeline

### Dry run

```bash
snakemake -n
```

### Full run with 8 cores

```bash
snakemake --cores 8 --use-conda
```

---

## Input Files (User‑Provided)

| File | Location | Description |
|------|----------|-------------|
| `ORF_names.csv` | `data/annotations/` | ORF name mapping (e.g., `ORF1a,1A`) |
| `domain_map.csv` | `data/annotations/` | Pfam domain → ORF type mapping |
| `host_mapping.csv` | `data/annotations/` | Raw host → scientific name mapping |
| `host_map.tsv` | `data/annotations/` | Regex patterns for host naming |
| `country_map.tsv` | `data/annotations/` | Regex patterns for country names |
| `exceptions.txt` | `data/annotations/` | Accessions to skip |
| `bastro_manual_ORF1b.txt` | `data/annotations/` | Manual overrides for bastroviruses |
| `MAstV_WG_2025_metadata.tsv` | `data/annotations/` | Reference metadata for ICTV species |
| `ORF1B_upd2_17nt_mmseq2.uc` | `data/results/clustering/` | mmseq2 clustering output |

---

## Output Files

| File | Location | Description |
|------|----------|-------------|
| `{base}_orf‑coords.csv` | `data/results/coords/` | Coordinate table (original ORF1a/1b/2 coordinates) |
| `{base}_orfs.faa` | `data/results/coords/` | Protein FASTA of extracted ORFs |
| `{base}_problematic_candidates.faa` | `data/results/coords/` | Protein FASTA of problematic CDSs (conflicting or missing annotations) |
| `{base}_problematic_candidates.tsv` | `data/results/coords/` | Mapping file for problematic CDSs (coordinates, gene, product, etc.) |
| `{base}_no_cds_sequences.fasta` | `data/results/coords/` | Sequences lacking any CDS annotation (input for Prodigal) |
| `{base}_metadata.tsv` | `data/results/coords/` | Metadata extracted from GenBank (taxonomy, references, features, etc.) |
| `logs/` (directory) | `data/results/coords/logs/` | Log files from ORF extraction (lists of skipped entries, no‑annotation summaries, etc.) |
| `annotated_domains.tbl` | `data/results/hmmer/annotated_orfs/` | HMMER domain table for annotated ORFs (ORF1a/1b/2) |
| `annotated_hmmscan.log` | `data/results/hmmer/annotated_orfs/` | HMMER execution log for annotated ORFs |
| `mismatches_report.tsv` | `data/results/hmmer/annotated_orfs/` | Report of annotation/domain mismatches (if any) |
| `problematic_domains.tbl` | `data/results/hmmer/problematic_orfs/` | HMMER domain table for CDS with ambiguous annotation |
| `problematic_hmmscan.log` | `data/results/hmmer/problematic_orfs/` | HMMER execution log for CDS with ambiguous annotation |
| `problematic_assigned.tsv` | `data/results/hmmer/problematic_orfs/` | Assigned ORF types for CDS with ambiguous annotation (based on Pfam domains) |
| `predicted_domains.tbl` | `data/results/hmmer/predicted_orfs/` | HMMER domain table for Prodigal‑predicted ORFs |
| `predicted_hmmscan.log` | `data/results/hmmer/predicted_orfs/` | HMMER execution log for predicted ORFs |
| `predicted_assigned.tsv` | `data/results/hmmer/predicted_orfs/` | Assigned ORF types for predicted ORFs (based on Pfam domains) |
| `{base}_no_cds_sequences_proteins.faa` | `data/results/prodigal/` | Protein FASTA from Prodigal (predicted ORFs) |
| `{base}_no_cds_predgenes.gff` | `data/results/prodigal/` | Gene predictions in GFF format |
| `{base}_no_cds_predgenes.fna` | `data/results/prodigal/` | Gene predictions in FASTA format |
| `prodigal.log` | `data/results/prodigal/` | Prodigal execution log |
| `{base}_orf‑coords_with_problematic.csv` | `data/results/updated_coords/` | Coordinates after merging with updated coordinates for problematic CDS |
| `{base}_orf‑coords_full.csv` | `data/results/updated_coords/` | Final coordinate table (including predicted ORFs) |
| `{base}_host_taxonomy.tsv` | `data/results/host_taxonomy/` | Host taxonomy mapping (species, genus, family, order, class, phylum, kingdom) |
| `{base}_annotated.tsv` | `data/results/final/` | All metadata + coordinates + host taxonomy |
| `{base}_annotated_MAstV_species.tsv` | `data/results/final/` | **Main result** – same as `annotated.tsv` plus `Species_ICTV` and `virus_name` columns |
| `{base}_{orf}.fna` (or `.faa`) | `data/results/orfs_named/` | FASTA files for separate ORFs. Sequences names are defined by user (generated by auxiliary script) |

---

## Known Limitations

- Host taxonomy retrieval may fail for poorly formatted host strings. The pipeline includes a mapping file (`host_mapping.csv`) to handle common cases.
- The Pfam database is large (~2.5 GB). Ensure sufficient disk space.


---

## Citation

If you use this pipeline or dataset, please cite:

- Snakemake: Köster & Rahmann, 2018
- Biopython: Cock et al., 2009
- HMMER: Eddy, 2011
- Prodigal‑gv: Hyatt et al., 2010

For the curated dataset, please cite the original GenBank entries and acknowledge the *Astroviridae* Study Group.

---

## Contact

For questions or issues, please open an [issue](https://github.com/v-julia/Astroviridae_database/issues) or contact the maintainer at `vjulia94@gmail.com`.


