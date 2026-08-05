"""
Functions for downloading GenBank records from NCBI, parsing GenBank files,
and extracting metadata (taxonomy, references, features, etc.).
"""



import os
import time
#import sys
from Bio import SeqIO, Entrez
import pandas as pd


'''
Downloads entries from GenBank Nucleotide by given query in batches with checkpoint/resume capability
Saves entries to file in GenBank format


Input:
        query (str): search query for GenBank Nucleotide.
        outfile (str): output filename (GenBank format).
        batch_size (int): number of IDs per efetch request (max ~200).
        checkpoint_file (str, optional): file to save progress (index of last downloaded ID).
Output:
    file with sequences in genbank format
'''
def fetch_seq_from_Nucleotide(query, outfile, batch_size=100, checkpoint_file=None):

    Entrez.email = "A.N.Other@example.com"

    # ---------------------------
    # 1. Get list of IDs from query
    # ---------------------------
    print("Query to GenBank Nucleotide database:\"{}\"".format(query))

    # list with ids obtained by query

    handle = Entrez.esearch(db="nucleotide", term=query, idtype="acc", RetMax=1000000)
    record = Entrez.read(handle)
    handle.close()
    id_list = record['IdList']
    total = len(id_list)

    print("Number of records found: {}".format(record["Count"]))
    print("{} records will be downloaded".format(total))

    # ---------------------------
    # 2. Determine starting point (resume)
    # ---------------------------

    start_index = 0
    if checkpoint_file and os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as cf:
            try:
                start_index = int(cf.read().strip())
                print(f"Resuming from index {start_index} (already downloaded {start_index} records)")
            except:
                print("Invalid checkpoint, starting from beginning")


    # ---------------------------
    # 3. Open output file in append mode if resuming
    # ---------------------------
    mode = 'a' if start_index > 0 else 'w'
    with open(os.path.abspath(outfile), mode) as file_out:
        # ---------------------------
        # 4. Process in batches
        # ---------------------------
        for i in range(start_index, total, batch_size):
            batch_ids = id_list[i:i+batch_size]
            ids_str = ",".join(batch_ids)
            print(f"Downloading batch {i//batch_size + 1}: IDs {i+1} to {min(i+batch_size, total)}")

            # Retry logic for HTTP errors
            retries = 5
            for attempt in range(retries):
                try:
                    handle = Entrez.efetch(db="nucleotide", id=ids_str,
                                           rettype="genbank", retmode="text")
                    for line in handle:
                        file_out.write(line)
                    handle.close()
                    break  # success, exit retry loop
                except Exception as e:
                    print(f"Error: {e}. Attempt {attempt+1}/{retries}")
                    if attempt < retries - 1:
                        sleep_time = 2 ** attempt  # exponential backoff: 1,2,4,8 sec
                        print(f"Waiting {sleep_time} seconds before retry...")
                        time.sleep(sleep_time)
                    else:
                        # If all retries fail, skip this batch and continue
                        print(f"Failed to download batch starting at index {i}. Skipping.")
                        # Optionally log the failed IDs for later manual recovery
                        with open("failed_batches.log", "a") as f:
                            f.write(f"{i}: {ids_str}\n")
                        # Continue to next batch

            # Update checkpoint after each successful batch
            if checkpoint_file:
                with open(checkpoint_file, 'w') as cf:
                    cf.write(str(i + len(batch_ids)))

            # Respect NCBI rate limit: at most 3 requests per second
            time.sleep(0.34)

    print(f"Finished. Output saved to {os.path.abspath(outfile)}")
    return(os.path.abspath(outfile))

def fetch_metadata_from_gb(input_file, output_dir):
    '''
    For each entry in GenBank file retrieves the following data:

    # Identifiers
    Accession
    GenBank title (DEFINITION)

    # Reference
    Reference_title
    Submitters
    PUBMED_ID
    
    # Classification
    Organism name
    Species
    Isolate
    Strain
    Family
    Lineage

    #Sequence quality
    Length
    Strand

    Num of annotated CDS

    ORF1a
    ORF1b
    ORF2
    coords
    
    #Sources
    Geo location
    Tissue/Specimen/Source
    Collection date
    Release date
    
    '''

    '''
     dictionary
     entries_data[GenBank Accession] = {}
     entries_data[GenBank Accession].keys() = [
                                                'Version',
                                                'GenBank title'
                                                'Organism name',
                                                'Species',
                                                'Family',
                                                'Virus Lineage',
                                                'Length',
                                                'Isolate',
                                                'Strain',
                                                'Strand',
                                                'Geo location'
                                                'Country',
                                                'Tissue/Specimen/Source',
                                                'Host',
                                                'Collection date',
                                                'Release date',
                                                'Reference_title',
                                                'Submitters',
                                                'PUBMED_ID',
                                                'Sequencing Technology'
                                                'CDS number',
                                                'Artificial',
                                                'UNVERIFIED'
                                                ]
    '''
    entries_data = {}

    source_qualifiers = ["strain", "isolate", "geo_loc_name", "country", "collection_date", "isolation_source", "host"]
    source_qualifiers_columns = ["Strain", "Isolate", "Geo location", "Country", "Collection date", "Tissue/Specimen/Source", "Host"]

    entries = SeqIO.parse(os.path.abspath(input_file), "genbank")

    for entry in entries:
        entry_name = entry.id.split(".")[0]
        entries_data[entry_name] = {}
        entries_data[entry_name]['LOCUS'] = entry.name
        entries_data[entry_name]['Version'] = entry.id
        entries_data[entry_name]['GenBank title'] = entry.description

        entries_data[entry_name]['Release date'] = entry.annotations['date']
        entries_data[entry_name]['Organism name'] = entry.annotations['organism']
        
        entries_data[entry_name]['Virus Lineage'] = ';'.join(entry.annotations['taxonomy'])
        num_taxa = len(entry.annotations['taxonomy'])
        if num_taxa == 9:
            entries_data[entry_name]['Species'] = entry.annotations['taxonomy'][-1]
        else:
             entries_data[entry_name]['Species'] = 'NA'
        try:
            entries_data[entry_name]['Family'] = entry.annotations['taxonomy'][6]
        except:
            entries_data[entry_name]['Family'] = ""
        
        entries_data[entry_name]['Length'] = len(entry)


        entries_data[entry_name]['Reference_title'] = ""
        entries_data[entry_name]['Submitters'] = ""
        entries_data[entry_name]['PUBMED_ID'] = ""
        
        if len(entry.annotations['references']) == 1:
            entries_data[entry_name]['Reference_title'] = entry.annotations['references'][0].title
            entries_data[entry_name]['Submitters'] = entry.annotations['references'][0].authors
            entries_data[entry_name]['PUBMED_ID'] = entry.annotations['references'][0].pubmed_id
        else:
            for reference in entry.annotations['references'][:-1]:
                entries_data[entry_name]['Reference_title'] += (reference.title + ';')
                entries_data[entry_name]['Submitters'] += (reference.authors  + ';')
                entries_data[entry_name]['PUBMED_ID'] += (reference.pubmed_id  + ';')

            entries_data[entry_name]['Reference_title'] = entries_data[entry_name]['Reference_title'].strip(";")
            entries_data[entry_name]['Submitters'] = entries_data[entry_name]['Submitters'].strip(";")
            entries_data[entry_name]['PUBMED_ID'] = entries_data[entry_name]['PUBMED_ID'].strip(";")
        #for reference in entry.annotations['references']:
        #    if reference.title == 'Direct Submission':
        #        entries_data[entry_name]['Submitters'] = reference.authors

        count_cds = 0
        for feature in entry.features:
            if feature.type == 'source':
                for qualif, colname in zip(source_qualifiers, source_qualifiers_columns):
                    if qualif in feature.qualifiers:
                        entries_data[entry_name][colname] = feature.qualifiers[qualif][0]
                    else:
                        entries_data[entry_name][colname] = 'NA'
                if 'environmental_sample' in feature.qualifiers:
                    entries_data[entry_name]['Environmental'] = 'Yes'
                else:
                    entries_data[entry_name]['Environmental'] = 'No'
                    
            if feature.type == 'CDS':
                
                count_cds +=1
                
        entries_data[entry_name]['CDS count'] = count_cds

        artificial_condition0 = ('patent' in entry.annotations['references'][-1].journal.lower())
        artificial_condition1 = ("nonfunctional" in entry.description.lower())

        if artificial_condition0 or artificial_condition1:

            entries_data[entry_name]['Artificial'] = "Yes"
        else:
            entries_data[entry_name]['Artificial'] = "No"


        if 'UNVERIFIED' in entry.description:
            entries_data[entry_name]['UNVERIFIED'] = "Yes"
        else:
            entries_data[entry_name]['UNVERIFIED'] = "No"
    meta_dataframe = pd.DataFrame.from_dict(entries_data, orient='index')


    base_name = os.path.splitext(os.path.basename(input_file))[0]
    if output_dir is None:
        # Original behaviour: save next to input file
        out_dir = os.path.dirname(os.path.abspath(input_file))
        out_file_name_temp = os.path.join(out_dir, base_name)
    else:
        # save into specified output directory, creating it if needed
        out_dir = os.path.abspath(output_dir)
        os.makedirs(out_dir, exist_ok=True)
        out_file_name_temp = os.path.join(out_dir, base_name)
    
    meta_dataframe.to_csv(out_file_name_temp + '_metadata.tsv', sep='\t', index_label="Accession")
    #country_map - file with abbreviations of countries
    return meta_dataframe
