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
| `{run_id}/coords/{base}_orf‑coords.csv` | `data/results/{run_id}/coords/` | Coordinate table (original ORF1a/1b/2 coordinates) |
| `{run_id}/coords/{base}_orfs.faa` | `data/results/{run_id}/coords/` | Protein FASTA of extracted ORFs |
| `{run_id}/coords/{base}_problematic_candidates.faa/.tsv` | `data/results/{run_id}/coords/` | FASTA/mapping for problematic CDSs |
| `{run_id}/coords/{base}_no_cds_sequences.fasta` | `data/results/{run_id}/coords/` | Sequences lacking CDS (input for Prodigal) |
| `{run_id}/coords/{base}_metadata.tsv` | `data/results/{run_id}/coords/` | Metadata extracted from GenBank |
| `{run_id}/coords/logs/` | `data/results/{run_id}/coords/logs/` | ORF extraction logs |
| `{run_id}/hmmer/annotated_orfs/annotated_domains.tbl` | HMMER domain table for annotated ORFs |
| `{run_id}/hmmer/annotated_orfs/mismatches_report.tsv` | Mismatch report |
| `{run_id}/hmmer/problematic_orfs/problematic_domains.tbl` | HMMER domain table for CDS with ambiguous annotation |
| `{run_id}/hmmer/problematic_orfs/problematic_assigned.tsv` | Assigned ORF types (CDS with ambiguous annotation) |
| `{run_id}/hmmer/predicted_orfs/predicted_domains.tbl` | HMMER domain table for predicted ORFs |
| `{run_id}/hmmer/predicted_orfs/predicted_assigned.tsv` | Assigned ORF types (predicted) |
| `{run_id}/prodigal/` | Prodigal outputs (protein FASTA, GFF, GBK, log) |
| `{run_id}/updated_coords/{base}_orf‑coords_full.csv` | Final coordinate table (including predicted ORFs) |
| `{run_id}/host_taxonomy/{base}_host_taxonomy.tsv` | Host taxonomy mapping |
| `{run_id}/final/{base}_annotated.tsv` | All metadata + coordinates + host taxonomy |
| `{run_id}/final/{base}_annotated_MAstV_species.tsv` | **Main result** – same plus `Species_ICTV` and `virus_name` |
| `{base}_{orf}.fna` (or `.faa`) | `data/results/orfs_named/` | Named ORF FASTA files (auxiliary script) |
| `latest -> {run_id}` | `data/results/latest` | Symlink pointing to the most recent run |

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


