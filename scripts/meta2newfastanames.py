import argparse
import os
import re
from Bio import SeqIO
import pandas as pd


# SETTING SHORT NAMES FOR COLUMNS
                                                                      
                                                                   

def create_mapping_dict(mapping_file_path):
    """
    Create a dictionary mapping regex patterns to short codes.
    
    Parameters:
    -----------
    mapping_file_path : str
        Path to CSV file with regex patterns and short codes
        
    Returns:
    --------
    dict : {regex_pattern: short_code}
    list : [(compiled_regex_pattern, short_code), ...]
    """
    # Read the mapping file
    mapping_df = pd.read_csv(mapping_file_path, sep='\t', header=None, 
                           names=['regex_pattern', 'short_code'])
    
    # Create two versions:
    # 1. Simple dict for reference
    pattern_to_code = dict(zip(mapping_df['regex_pattern'], mapping_df['short_code']))
    #print(mapping_df)
    # 2. List of compiled regex patterns with their codes
    compiled_patterns = []
    for _, row in mapping_df.iterrows():
        try:
            # Compile the regex pattern
            pattern = re.compile(str(row['regex_pattern']), re.IGNORECASE)
            compiled_patterns.append((pattern, row['short_code']))
        except re.error as e:
            print(f"Warning: Invalid regex pattern '{row['regex_pattern']}': {e}")

    return pattern_to_code, compiled_patterns


def map_value_to_short_code(value, compiled_patterns, default="NA"):
    """
    Map a value to its short code using regex patterns.
    
    Parameters:
    -----------
    value : str
        The value to map
    compiled_patterns : list
        List of (compiled_regex_pattern, short_code) tuples
    default : str
        Default value if no pattern matches
        
    Returns:
    --------
    str : The short code
    """
    if pd.isna(value):
        return default
    

    
    # Try each pattern in order
    found = 0
    for pattern, short_code in compiled_patterns:
        if pattern.search(value):
            #print(pattern, value,short_code )
            #print(pattern.search(value))
            found = 1
            return short_code
            
    if found == 0:
        print("{} was not found in mapping file".format(value))
    
    return value.replace(" ", "-").replace(",", "-")


# Map values in column
                                                              

def process_column(meta_df, mapping_file_path, column="Host"):
    """
    Process the column to replace long names with short codes.
    
    Parameters:
    -----------
    meta_df : pandas.DataFrame
        Metadata dataframe
    mapping_file_path : str
        Path to regex mapping TSV file
    host_column : str
        Name of the host column
        
    Returns:
    --------
    pandas.DataFrame : Updated dataframe
    dict : Accession to string mapping
    """
    # Create regex mapping
    pattern_to_code, compiled_patterns = create_mapping_dict(mapping_file_path)
    
    # Create a copy to avoid modifying original
    meta_processed = meta_df.copy()
    
    # Apply the mapping to Host column
    meta_processed[column + "_original"] = meta_processed[column]  # Keep original
    meta_processed[column] = meta_processed[column].apply(
        lambda x: map_value_to_short_code(x, compiled_patterns)
    )
    
    return meta_processed

                  

# EXTACTING YEAR FROM COLLECTION DATE

def extract_year_from_date(date_string):
    """
    Extract year (4 digits) from a date string.
    Returns the year as string, or 'NA' if not found.
    """
    if pd.isna(date_string):
        return "NA"
    
    flexible_pattern = r'\b\d{4}\b'
    matches = re.findall(flexible_pattern, date_string)
    
    found = 0
    if matches:
        # Filter to reasonable years (e.g., 1900-2099)
        for match in matches:
            year_int = int(match)
            if 1900 <= year_int <= 2099:
                found = 1
                return match
        # If no reasonable year found, return the first 4-digit number
        return matches[0]
    if found == 0:
        print("Year was not found")
        print(date_string)
    return date_string


def clean_value(value):
    """Clean a single value: handle NaN and replace slashes"""
    if pd.isna(value):
        return "NA"
    return str(value).replace('/', '-').replace(',', '-').replace(' ', '-').replace(')','-').replace('(','-').strip()
                      
                           
                           
                           
                           
                           
                


def create_string(row, columns):
    """Create string from selected columns"""
    parts = [str(row[col]) for col in columns]
    return "|".join(parts)


                    
                                       


def meta2fastaname(fasta_path, meta_path, countrymap_path, meta_columns, out_path):
    """Changes sequences names in fasta file"""
    meta = pd.read_csv(meta_path,sep="\t")

                                               
    meta = process_column(meta, countrymap_path, column="Geo location")

    meta['Year'] = meta['Collection date'].apply(extract_year_from_date)
    
    
    columns_to_include = meta_columns.split(",")#["Accession", "Geo location", "Isolate", "Host", "class", "Year"]

    for col in columns_to_include:
        meta[col] = meta[col].apply(clean_value)

    # Create the strings as a new column
    meta['info_string'] = meta.apply(
        lambda row: create_string(row, columns_to_include),
        axis=1
    )
    

    new_seq_names = dict(zip(meta["Accession"], meta["info_string"]))
    
    
    seqs = list(SeqIO.parse(fasta_path, 'fasta'))
    print('Total number of records: {}'.format(len(seqs)))
    
    new_seqs = []
    k=0

    for seq in seqs:
        ac = seq.name.split('/')[0]
        seq.name = new_seq_names[ac]
        seq.id = new_seq_names[ac]
        #print(seq.name)
        seq.description = ''

        new_seqs.append(seq)
        k+=1
    print('Total number of processed records: {}'.format(k))
    
    base_name = os.path.basename(fasta_path)
    #SeqIO.write(new_seqs, os.path.splitext(fasta_path)[0] + '_ren.fasta', 'fasta')
    with open(os.path.join(out_path, base_name), 'w', encoding='utf-8') as output_handle:
        SeqIO.write(new_seqs, output_handle, 'fasta')

                                  


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-in_fasta", "--input_fasta", type=str,
                        help="Input alignment with sequences in fasta format.", required=True)
    parser.add_argument("-in_meta", "--input_meta", type=str,
                        help="Input metafile", required=True)
    parser.add_argument("-country_map", "--country_map", type=str,
                        help="TSV file with regular expressions for countries and their ISO codes", required=True)
    parser.add_argument("-columns", "--columns", type=str,
                        help="Columns to be included in new sequence names, should be separated by comma. Example: 'Accession,Geo location,Isolate,Host,class,Year'", required=True)
    parser.add_argument("-out", "--out_path", type=str,
                        help="Path to output directory", required=True)
    args = parser.parse_args()
    meta2fastaname(args.input_fasta, args.input_meta, args.country_map, args.columns, args.out_path)

                   
                      
                     
                      
                         
                    
     