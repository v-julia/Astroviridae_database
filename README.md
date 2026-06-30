# Astroviridae Database

This repository contains scripts for annotating astrovirus CDSs from GenBank entries, along with a curated dataset of astrovirus and astro‑like sequences deposited in GenBank.

## Dataset
`data/Astroviridae_15102025_curated.csv` is a curated, semicolon-separated CSV file containing metadata extracted from GenBank entries, coordinates of astrovirus ORFs, and host class information. For some sequences, the species designation according to the ICTV Astroviridae Study Group is also indicated.

The dataset includes 16,012 sequences downloaded on 15 October 2025. The file can be imported into Excel or any other spreadsheet software.

This repository contains scripts for annotation of astrovirus CDSs from genbank entries and a curated dataset of astroviruses and astro-like viruses deposited in GenBank.

### Columns extracted from GenBank entries

| Column | Description |
|--------|-------------|
| `ID` | GenBank accession number |
| `Version` | GenBank entry version |
| `GenBank title` | Entry title |
| `Release date` | Date of release |
| `Organism name` | Organism name as in GenBank |
| `Virus Lineage` | Taxonomic lineage from the GenBank taxonomy field |
| `Species` | Virus species from the taxonomy field (may contain errors) |
| `Family` | Virus family |
| `Length` | Sequence length |
| `Submitters` | Submitter information |
| `Strain` | Strain name (if available) |
| `Isolate` | Isolate name (if available) |
| `Geo location` | Geographic location |
| `Collection date` | Collection date |
| `Tissue/Specimen/Source` | Value of the `isolation_source` qualifier |
| `Host` | Host organism |
| `Environmental` | `yes` if the entry has an `environmental_sample` qualifier |
| `CDS count` | Number of annotated CDSs in the sequence |
| `Artificial` | `yes` if the description contains “modified microbial nucleic acid”, “patent”, or “oligonucleotide” |
| `UNVERIFIED` | `yes` if the entry title contains “UNVERIFIED” |

### ORF coordinates

For each sequence, the table provides the coordinates of the three main astrovirus ORFs when present:

- `1A` – coordinates of ORF1a (`NA-NA` if absent)  
- `1A-strand` – strand orientation  
- `1B` – coordinates of ORF1b  
- `1B-strand` – strand orientation  
- `2` – coordinates of ORF2  
- `2-strand` – strand orientation  

### Host classification

The column `Host class` contains the class-rank taxonomic lineage of the host, retrieved from the NCBI Taxonomy database. The dataset includes sequences from the following host classes:

- Mammals (*Mammalia*)
- Birds (*Aves*)
- Lizards and snakes (*Lepidosauria*)
- Amphibians (*Amphibia*)
- Fish (*Actinopteri*, *Chondrichthyes*, *Hyperoartia*)
- Arthropods (*Arachnida*, *Insecta*, *Malacostraca*)
- Annelid worms (*Clitellata*)
- Mollusks (*Bivalvia*)
- Flowering plants (*Magnoliopsida*)


### Species assignment

The column `Species_ICTV` contains the ICTV species designation for a subset of classified mamastroviruses. This column will be expanded as the study group progresses.

### Known limitations

ORF coordinate annotations may contain inconsistencies and errors. The coordinates were semi-automatically inspected and curated for entries with ORF1b sequences (>1300 nt). Coordinates for other ORFs may still contain inaccuracies. The annotation pipeline is currently being refined to enable rapid curation of newly deposited sequences.
