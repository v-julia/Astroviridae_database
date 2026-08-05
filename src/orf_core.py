"""
Core ORF extraction: identifies coordinates of ORF1a, ORF1b, and ORF2 from
GenBank records, validates translations, corrects frames, and exports
protein FASTA files and problematic candidate lists.
"""

import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq

from .annotation_utils import read_exception_file,read_csv,map_feature,map_keywords


#===== Validation of frames for CDSs in records that belong to target CDS ========
def validate_and_correct_record(record, coord_dict, transl_table=1, max_stops=0):
    """
    For each ORF in coord_dict, translate the nucleotide sequence using the stored
    start/end/strand. If internal stop codons are found, try shifting the reading
    frame by 0,1,2 bases. If a clean frame is found, update the coordinates;
    otherwise remove the ORF (so it will be reported as NA).
    """
    for orf in ['1A', '1B', '2']:
        if orf not in coord_dict or orf+'-strand' not in coord_dict:
            continue
        start, end = coord_dict[orf]
        strand = coord_dict[orf+'-strand']
        if start is None or end is None:
            continue

        # Extract the coding sequence (reverse‑complement for minus strand)
        if strand == 1:
            seq = record.seq[start:end]
        else:
            seq = record.seq[start:end].reverse_complement()

        # Translate (include stop codon to check internal stops)
        trans = seq.translate(table=transl_table, to_stop=False)
        internal_stops = [i for i, aa in enumerate(trans) if aa == '*' and i != len(trans)-1]

        if len(internal_stops) == 0:
            # Already clean
            continue

        # Try frame shifts (skip 0,1,2 bases at the beginning)
        best_offset = None
        best_trans = None
        best_stop_count = None
        for offset in range(3):
            if offset >= len(seq):
                continue
            subseq = seq[offset:]
            subtrans = subseq.translate(table=transl_table, to_stop=False)
            stops = [i for i, aa in enumerate(subtrans) if aa == '*' and i != len(subtrans)-1]
            stop_count = len(stops)
            if stop_count <= max_stops:
                if (best_offset is None or
                    stop_count < best_stop_count or
                    (stop_count == best_stop_count and len(subtrans) > len(best_trans))):
                    best_offset = offset
                    best_trans = subtrans
                    best_stop_count = stop_count

        if best_offset is not None and best_stop_count == 0:
            # Correction found – adjust genomic coordinate
            if strand == 1:
                # Shift start forward
                coord_dict[orf] = [start + best_offset, end]
            else:
                # Shift end backward (5' end is at end coordinate for minus strand)
                coord_dict[orf] = [start, end - best_offset]
            # strand remains unchanged
            print(f"Corrected {orf} for {record.id}: shifted by {best_offset} bases")
        else:
            # Cannot obtain a clean ORF – remove it
            print(f"Removing {orf} for {record.id} – cannot resolve internal stop codons")
            del coord_dict[orf]
            del coord_dict[orf+'-strand']

#===== Export target CDSs from validated records ========
def export_orfs_to_fasta(records, coord_dict, fasta_filename, transl_table=1):
    """
    Extracts nucleotide sequences according to coordinates in coord_dict,
    translates them to protein, and writes a FASTA file.
    
    Parameters:
        records       – list of SeqRecord objects already loaded from the GenBank file
        coord_dict    – dictionary {accession: {orf: [start, end], orf+'-strand': strand}}
        fasta_filename – output FASTA file name (protein sequences)
        transl_table  – genetic code number (default 1)
    """
    with open(fasta_filename, 'w') as out_fasta:
        for rec in records:
            acc = rec.name
            if acc not in coord_dict:
                continue
            orf_data = coord_dict[acc]
            for orf in ['1A', '1B', '2']:
                if orf not in orf_data or orf+'-strand' not in orf_data:
                    continue
                start, end = orf_data[orf]
                strand = orf_data[orf+'-strand']
                if start is None or end is None:
                    continue
                # Extract sequence
                if strand == 1:
                    seq = rec.seq[start:end]
                else:
                    seq = rec.seq[start:end].reverse_complement()
                # Translate (stop codon removed)
                protein = str(seq.translate(table=transl_table, to_stop=True))
                # FASTA header
                header = (f">{acc}|{orf}|{start}-{end}|strand={strand}|len={len(protein)}")
                out_fasta.write(header + "\n")
                # Write protein in lines of 60 characters
                for i in range(0, len(protein), 60):
                    out_fasta.write(protein[i:i+60] + "\n")

#======= Analyses CDS feature which annotation does not correspond to target CDS =======
def analyze_cds(record, feature, transl_table=1):
    """
    Extracts translation and coordinates from a CDS feature.
    Returns a dict with details, including internal stop count.
    """
    
    codon_start = int(feature.qualifiers.get("codon_start", [1])[0])
    offset = codon_start - 1


    # Get parts – if not a compound location, treat as single part
    if len(feature.location.parts) > 1:
        parts = feature.location.parts
    else:
        parts = [feature.location]

    
    candidates = []
    for idx, part in enumerate(parts):
        raw_seq = part.extract(record.seq)
        # codon_start offset only applies to the first part
        if idx == 0:
            cds_seq = raw_seq[offset:]
        else:
            cds_seq = raw_seq

        # Translate, include stop codon
        translation = cds_seq.translate(table=transl_table, to_stop=False)
        annotated_translation = feature.qualifiers.get("translation", [""])[0]
        
        # Count internal stops (* not at last position)
        stops = [i for i, aa in enumerate(translation) if aa == "*" and i != len(translation)-1]
        internal_stop_count = len(stops)
        
        if translation and translation[-1] == "*":
            translation = translation[:-1]

        if len(parts) == 1 and translation != annotated_translation:
            print(f"For {record.id} the annotated translation does not correspond to actual translation of coding sequence")
            print("Actual translation: " + translation)
            print("Annotated translation: " + annotated_translation)


        start = int(part.start)
        end = int(part.end)
        
        cand = {
            "record_name": record.name,
            "accession": record.id,
            "location": str(part),                # individual part location
            "codon_start": str(codon_start if idx == 0 else 1),
            "start": start,
            "end": end,
            "strand": part.strand,
            "gene": feature.qualifiers.get("gene", ["unknown"])[0],
            "product": feature.qualifiers.get("product", ["unknown"])[0],
            "transl_table": transl_table,
            "raw_nucleotide": str(raw_seq),
            "cds_sequence": str(cds_seq),
            "translation": str(translation),
            "translation_annot": str(annotated_translation),
            "internal_stop_count": internal_stop_count,
            "is_clean": internal_stop_count == 0,
            "part_index": idx,   # optional, for debugging
        }
        candidates.append(cand)
    return candidates

#======= Validation and correction of frame in CDS feature which annotation does not correspond to target CDS =======
def correct_cds_frame(cds_info, max_stops=2):
    """
    Attempt to correct the translation frame of a CDS.
    Tries all 3 reading frames on the raw nucleotide sequence.
    Returns the best translation and the offset to apply.

    Parameters:
    cds_info - dictionary with information of CDS feature, output of analyze_cds funciton
    max_stops - maximal number of stop codons
    """
    raw = cds_info["raw_nucleotide"]  # use the untrimmed raw sequence
    transl_table = cds_info["transl_table"]
    
    best = None
    best_offset = None
    
    for frame in range(3):  # 0, 1, 2
        seq = Seq(raw[frame:])
        
        trans = seq.translate(table=transl_table, to_stop=False)
        # count internal stops
        internal_stops = [i for i, aa in enumerate(trans) if aa == "*" and i != len(trans)-1]
        stop_count = len(internal_stops)
        
        # We prefer fewer stops, and if equal, longer length
        if stop_count <= max_stops:
            if best is None or (stop_count < best["stop_count"]) or \
               (stop_count == best["stop_count"] and len(trans) > len(best["translation"])):
                best = {
                    "translation": str(trans),
                    "stop_count": stop_count,
                    "frame_offset": frame,  # number of bases to skip from the raw start
                    "nucleotide_used": str(seq),
                }
                best_offset = frame
    return best, best_offset


def is_duplicate_candidate(candidate, coord_dict):
    """
    Return True if candidate's coordinates and strand match an ORF already
    stored in coord_dict for the same accession.
    """
    acc = candidate['accession']
    if acc not in coord_dict:
        return False
    orf_data = coord_dict[acc]
    for orf in ['1A', '1B', '2']:
        if orf in orf_data and orf+'-strand' in orf_data:
            start, end = orf_data[orf]
            if candidate['start'] == start and candidate['end'] == end:
                if candidate['strand'] == orf_data[orf+'-strand']:
                    return True
    return False
    

# ======= write CDS with no target or conflicting annotation
def write_fasta_and_mapping(candidates, fasta_file, mapping_file):
    """
    Writes two files:
      - fasta_file: amino acid sequences with header containing accession and some metadata.
      - mapping_file: TSV (or CSV) with detailed coordinates/strand per candidate.
    """
    with open(fasta_file, "w") as f_fasta, open(mapping_file, "w", newline="") as f_map:
        writer = csv.writer(f_map, delimiter="\t")
        writer.writerow(["candidate_id",
            "record_name", "accession", "start", "end", "codon_start","strand",
            "gene", "product", "is_corrected", "correction_offset", "original_translation"
        ])
        for idx, c in enumerate(candidates, 1):
            # FASTA header: include ID and some quick info
            cid = f"cand_{idx:04d}"
            header = f">{cid}|{c['accession']}|{c['gene']}|start={c['start']}|end={c['end']}|codon_start={c['codon_start']}|strand={c['strand']}"
            f_fasta.write(header + "\n")
            trans = c['translation']
            for i in range(0, len(trans), 60):
                f_fasta.write(trans[i:i+60] + "\n")
            
            # Write mapping line
            writer.writerow([
                cid,
                c['record_name'],
                c['accession'],
                c['start'],
                c['end'],
                c['codon_start'],
                c['strand'],
                c['gene'],
                c['product'],
                c.get('is_corrected', False),
                c.get('correction_frame_offset', ''),
                c.get('original_translation', '')
            ])

# ======= Get coordinates of CDS that belong to ORF1a, ORF1b, ORF2, write translated sequences to fasta. Gether CDS with conflicts or no target annotation =====
def orf_coord_updated(input_file, orf_map, output_dir):
    '''
    Retrieves the coordinates of ORFs
    
    Input:
        input_file - file with nucleotide sequences in genbank-format
        orf_map - csv file annotation of orfs and their codes
    Output:
        coord_file - file with coordinates
    '''

    # ====== LOAD ANNOTATIONS
    # dictionary with possiple names of astrovirus ORFs indicated in note, gene and product qualifiers
    # orf_dict["ORF1a"] = "1A"
    orf_dict = read_csv(orf_map)

    # ORFs which coordinates this script will extract
    orf_types_final = ['1A', '1B', '2']
    # ORFs that can be met
    orf_types = ['1A', '1B', '1AB', '2']
 
    # Keywords in description that can indicate which ORF is covered  by the sequence
    description_keywords = {"ORF2-like": "2", 
                            "ORF2": "2",
                            "precapsid":"2", 
                            "CP sequence": "2", 
                            "capsid": "2",
                            "ORF1a-like": "1A", 
                            "nsp1a-like": "1A", 
                            "ORF1A": "1A",
                            "1b-like": "1B", 
                            "RNA polymerase": "1B", 
                            "RdRp": "1B",
                            "ORF1b": "1B",
                            "polymerase":"1B"
                            }

    # entries to escape
    exceptions = read_exception_file("data/annotations/exceptions.txt")

    # Files with GBAC of entries that contain signle or no annotated CDS
    # with vague product annotation such as 'hypotetical protein'.
    # These entries were manually checked
    #exceptions_ORF1a =  read_exception_file("data/annotations/ORF1a_manual.txt")
    #exceptions_ORF1b =  read_exception_file("data/annotations/ORF1b_manual.txt")
    #exceptions_ORF2 =  read_exception_file("data/annotations/ORF2_manual.txt")

    # we will extract ORF1b of these  bastroviruses 
    bastro_ORF1b = read_exception_file("data/annotations/bastro_manual_ORF1b.txt")
    
    ## ENTRIES with known errors in annotation
    #conflict_annot_ORF2 = ["MN725025", "KY047739", "S68561", "KP404151","KP404152"] # the last CDS in conflict ORF1b ORF2 is ORF2
    # PV999257 is actually ORF1a
    #conflict_annot_1a = ["PV999257"]
    #AB518702 AB518703 is actually RdRp
    #conflict_annot_1b = []


    # ===== OUTPUT DIRECTORIES ====
    
    # basename for output files
    base_name = os.path.splitext(os.path.basename(input_file))[0]

    # set output directory
    if output_dir is None:
        # save next to input file
        out_dir = os.path.dirname(os.path.abspath(input_file))
        out_file_name_temp = os.path.join(out_dir, base_name)
    else:
        # save into specified output directory, creating it if needed
        out_dir = os.path.abspath(output_dir)
        os.makedirs(out_dir, exist_ok=True)
        out_file_name_temp = os.path.join(out_dir, base_name)
        
    out_file_name = out_file_name_temp  + '_orf-coords.csv'
    print(out_file_name)

    # Create a subfolder for logs (inside the same output directory)
    log_dir = os.path.join(out_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    print(f"Output files will be saved to: {out_dir}")

    # ====== PARSE GENBANK
    
    # list with dictionaries cds_info for CDS features that have conflicting annotation or no annotation of target CDS
    problematic_candidates = []

    # entries with no target CDS
    entries_with_notarget_cds = set()
    # dictionary with ORF coordinates
    dict_coord = {}


    with open(input_file) as handle:
        records = list(SeqIO.parse(handle, 'gb'))

        # entries with no annotations in product, gene, note qualifiers
        no_annot = 0
        
        
        # CDS with no annotations of target ORFs (ORF1a, ORF1b, ORF2)
        notarget_annot_CDS = 0
        # CDS with no annotations in product, gene, note qualifiers
        no_annot_CDS = 0

        # Entries with CDS that has annotation "nonstructural protein"
        entries_onlynonstrprotein = []
        entries_onlynonstrprotein_short = []

        # Entries with conflicting annotations
        entries_conflicting_annotations = []

        # Entries with no annotated CDSs at all. entry accession and description are written to this list
        entries_noannot = []
        # accessions of entries with no annotated CDS
        no_cds_accessions = []

        # Entries with no annotated target ORFs
        entries_notarget_annot = []
        
        for rec in records:
            rec_name = rec.id.split(".")[0]
            if rec_name in exceptions:
                continue
            
            # Omit unverified sequences
            if "UNVERIFIED" in rec.description:
                continue
            'patent' in rec.annotations['references'][-1].journal.lower()

            artificial_condition0 = ('patent' in rec.annotations['references'][-1].journal.lower())
            artificial_condition1 = ("nonfunctional" in rec.description.lower())
            
             # Omit artificial sequences
            if artificial_condition0 or artificial_condition1:
                continue

            #CDS in records that contain no target CDS annotation

            # Dictionary with coordinates for ORFs in record
            dict_coord[rec_name] = {}
            # iterate over record features

            # number of CDS
            cds_count = 0
            # number of annotated CDS with target ORFs
            cds_annot_count = 0

            # Set coordinates for nonstructural protein which type cannot be inferred from  annotations
            NSP_start, NSP_end, NSP_strand = None, None, None
            for feature in rec.features:

                # check primers to identify the amplified region
                # !!! REVISE removing Mon269_Mon270
                if feature.type == 'source':
                    Mon340_Mon348 = False
                    Mon269_Mon270 = False
                    if feature.qualifiers.get("PCR_primers") != None:
                        #print(feature.qualifiers["PCR_primers"])
                        if "mon340" in feature.qualifiers["PCR_primers"][0].lower() or "mon348" in feature.qualifiers["PCR_primers"][0].lower():
                            Mon340_Mon348 = True
                            Mon269_Mon270 = False
                        elif "mon269" in feature.qualifiers["PCR_primers"][0].lower() or "mon270" in feature.qualifiers["PCR_primers"][0].lower():
                            Mon340_Mon348 = False
                            Mon269_Mon270 = True
        
                if feature.type == 'CDS':
                    cds_count += 1
                    
                    # Check codon start
                    if 'codon_start' in feature.qualifiers.keys():
                        cod_start = int(feature.qualifiers['codon_start'][0]) - 1
                        if cod_start<0:
                            print('Codon start < 0, something is wrong')
                    else:
                        cod_start = 0

                    # check all annotations in gene, product, note qualifiers
                    CDS_raw_annotations = []
                    for el in ('product', 'gene', 'note'):
                        CDS_raw_annotations.extend(feature.qualifiers.get(el, []))
                    
                    if len(CDS_raw_annotations) != 0:
                        # now we will leave annotations of target ORFs
                        CDS_annotations = list(map(lambda x: map_feature(x, orf_dict), CDS_raw_annotations))
                        CDS_annotations = list(filter(lambda x: x in orf_types, CDS_annotations))

                        # IF TARGET ORF ANNOTATIONS WERE FOUND
                        if len(CDS_annotations) != 0:
                                                    
                            #print("Annotations were successfully found")
                            #print(CDS_annotations)

                            # check for consistency of CDS annotations
                            uniq_annotations = list(set(CDS_annotations))
                            if len(uniq_annotations) == 1:
                                CDS_product = uniq_annotations[0]
                                # CDS is a target ORF, the coordinates are not joined
                                if CDS_product in orf_types_final:
                                    dict_coord[rec_name][CDS_product] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                    dict_coord[rec_name][CDS_product + '-strand'] = feature.location.strand
                                    cds_annot_count +=1
                                elif CDS_product == "1AB":

                                    # check whether the coordinates are joined
                                    if len(feature.location.parts) == 2:
                                        dict_coord[rec_name]['1A'] = [int(feature.location.parts[0].start) + cod_start, int(feature.location.parts[0].end)]
                                        dict_coord[rec_name]['1A-strand'] = feature.location.parts[0].strand
                                        dict_coord[rec_name]['1B'] = [int(feature.location.parts[1].start), int(feature.location.parts[1].end)]
                                        dict_coord[rec_name]['1B-strand'] = feature.location.parts[1].strand
                                        cds_annot_count +=1
                                    elif len(feature.location.parts) == 1:
                                        dict_coord[rec_name]['1B'] = [int(feature.location.parts[0].start) + cod_start, int(feature.location.parts[0].end)]
                                        dict_coord[rec_name]['1B-strand'] = feature.location.parts[0].strand
                                        cds_annot_count +=1
                                    else:
                                        print("strange locations in {}".format(rec_name))
                                        print(feature.location)
                            else:
                                print("Conflicting CDS annotations for {}".format(rec_name))
                                #print(CDS_annotations)
                                #print(uniq_annotations)

                                present_1A = ('1A' in uniq_annotations)
                                present_1B = ('1B' in uniq_annotations)
                                present_1A_1B = present_1A and present_1B
                                present_1AB = ('1AB' in uniq_annotations)
                                present_1A_1AB = present_1A and present_1AB
                                present_1B_1AB = present_1B and present_1AB
                                present_2 = ('2' in uniq_annotations)
                                present_1AB_2 = (present_1A and present_2) or (present_1B and present_2) or (present_1AB and present_2) 
                                
                                
                                if present_1AB_2:
                                    #if rec.name in conflict_annot_ORF2:
                                    #    dict_coord[rec.name]['2'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                    #    dict_coord[rec.name]['2-strand'] = feature.location.strand
                                    #    cds_annot_count +=1
                                    #    print("Resolved for ORF2")
                                    #else:
                                    if rec_name in bastro_ORF1b:
                                        dict_coord[rec_name]['1B'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                        dict_coord[rec_name]['1B-strand'] = feature.location.strand
                                        cds_annot_count +=1
                                    else:
                                        print("Could not resolve the conflict")
                                        entries_conflicting_annotations.append(rec.id + ': ORF1ab + ORF2' + ':' + rec.description + '\n')

                                        cds_infos = analyze_cds(rec, feature)
                                        for cds_info in cds_infos:
                                            if not cds_info["is_clean"]:
                                                corrected, offset = correct_cds_frame(cds_info, max_stops=1)
                                                if corrected is not None:
                                                    cds_info["is_corrected"] = True
                                                    cds_info["corrected_translation"] = corrected["translation"]
                                                    cds_info["original_translation"] = cds_info["translation"]
                                                    cds_info["translation"] = corrected["translation"]
                                                    cds_info["correction_frame_offset"] = offset
                                                else:
                                                    cds_info["is_corrected"] = False
                                                    cds_info["correction_failed"] = True
                                            problematic_candidates.append(cds_info)

                                # If ORF1a and ORF1b are present, we choose ORF1b
                                #elif present_1A_1B:
                                #    if rec.name in conflict_annot_1a:
                                #        dict_coord[rec.name]['1A'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                #        dict_coord[rec.name]['1A-strand'] = feature.location.strand
                                #        cds_annot_count +=1
                                #        print("Resolved for ORF1a")
                                #    else:
                                #        dict_coord[rec.name]['1B'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                #        dict_coord[rec.name]['1B-strand'] = feature.location.strand
                                #        cds_annot_count +=1
                                #        #entries_conflicting_annotations.append(rec.name + ': ORF1a + ORF1b' + ':' + rec.description + '\n')
                                #        print("Resolved for ORF1b")
                                # ORF1ab with joined coordinates
                                elif len(feature.location.parts) == 2:
                                    dict_coord[rec_name]['1A'] = [int(feature.location.parts[0].start) + cod_start, int(feature.location.parts[0].end)]
                                    dict_coord[rec_name]['1A-strand'] = feature.location.parts[0].strand
                                    dict_coord[rec_name]['1B'] = [int(feature.location.parts[1].start), int(feature.location.parts[1].end)]
                                    dict_coord[rec_name]['1B-strand'] = feature.location.parts[1].strand
                                    cds_annot_count +=1
                                    print("Resolved for ORF1ab")
                                elif present_1A_1AB and ('1A' not in dict_coord[rec_name].keys()):
                                    dict_coord[rec_name]['1A'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                    dict_coord[rec_name]['1A-strand'] = feature.location.strand
                                    cds_annot_count +=1
                                    print("Resolved for ORF1a")
                                elif present_1B_1AB and ('1B' not in dict_coord[rec_name].keys()):
                                    dict_coord[rec_name]['1B'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                    dict_coord[rec_name]['1B-strand'] = feature.location.strand
                                    cds_annot_count +=1
                                    print("Resolved for ORF1b")
                                elif '1B' in dict_coord[rec_name].keys() or '1A' in dict_coord[rec_name].keys():
                                    print("The coordinates've been extracted from another CDS")
                            if rec_name == "MF973517":
                                print("MF973517")
                                print(cds_annot_count)
                                print(dict_coord[rec_name])
                        # NO TARGET ORF ANNOTATIONS WERE FOUND
                        else:            
                            #print("Could no find annotation of target ORF for CDS in {}".format(rec.name))
                            #print(CDS_raw_annotations)
                            #print(rec.description)
                            notarget_annot_CDS +=1

                            # Gather CDS info
                            cds_infos = analyze_cds(rec, feature)
                            for cds_info in cds_infos:
                                if not cds_info["is_clean"]:
                                    corrected, offset = correct_cds_frame(cds_info, max_stops=1)
                                    if corrected is not None:
                                        cds_info["is_corrected"] = True
                                        cds_info["corrected_translation"] = corrected["translation"]
                                        cds_info["original_translation"] = cds_info["translation"]
                                        cds_info["translation"] = corrected["translation"]
                                        cds_info["correction_frame_offset"] = offset
                                    else:
                                        cds_info["is_corrected"] = False
                                        cds_info["correction_failed"] = True
                                problematic_candidates.append(cds_info)
                                entries_with_notarget_cds.add(";".join([rec.id, rec.description]))
                            #if rec.name in exceptions_ORF1a:
                            #    dict_coord[rec.name]["1A"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                            #    dict_coord[rec.name]['1A-strand'] = feature.location.strand
                            #    cds_annot_count +=1
                            #    cds_count +=1
                            #if rec.name in exceptions_ORF1b:
                            #    dict_coord[rec.name]["1B"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                            #    dict_coord[rec.name]['1B-strand'] = feature.location.strand
                            #    cds_annot_count +=1
                            #    cds_count +=1
                            #if rec.name in exceptions_ORF2:
                            #    dict_coord[rec.name]["2"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                            #    dict_coord[rec.name]['2-strand'] = feature.location.strand
                            #    cds_annot_count +=1
                            #    cds_count +=1
                            #
                            if "nonstructural protein" in CDS_raw_annotations:
                                #print("CDS of {} has the only annotation \'nonstructural protein\' ".format(rec_name))
                                if len(rec) > 300:
                                    if "orf1a" in rec.description.lower():
                                        dict_coord[rec_name]["1A"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                        dict_coord[rec_name]['1A-strand'] = feature.location.strand
                                        cds_annot_count +=1
                                        cds_count +=1
                                    # sequences covers ORF1b and ORF2
                                    elif  len(rec)< 700 and ("capsid" in rec.description.lower()):
                                        dict_coord[rec_name]["1B"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                        dict_coord[rec_name]['1B-strand'] = feature.location.strand
                                        cds_annot_count +=1
                                        cds_count +=1
                                    else:
                                        # ORF1b with joined coordinates
                                        #print(feature.location.parts)
                                        if len(feature.location.parts) == 2:
                                            dict_coord[rec_name]['1A'] = [int(feature.location.parts[0].start) + cod_start, int(feature.location.parts[0].end)]
                                            dict_coord[rec_name]['1A-strand'] = feature.location.parts[0].strand
                                            dict_coord[rec_name]['1B'] = [int(feature.location.parts[1].start), int(feature.location.parts[1].end)]
                                            dict_coord[rec_name]['1B-strand'] = feature.location.parts[1].strand
                                            cds_annot_count +=1
                                        #else:
                                        #    NSP_start = int(feature.location.start) + cod_start
                                        #    NSP_end = int(feature.location.parts[0].end)
                                        #    NSP_strand = feature.location.strand
                                        #    #entries_onlynonstrprotein.append(rec_name + ':' + str(len(rec)) + ":"  +rec.description + '\n')
                                        #    # Тут записать координаты рамки и в конце проверить, нет ли других аннотаций

                                else:
                                    if Mon340_Mon348 == True:
                                        dict_coord[rec_name]["1A"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                        dict_coord[rec_name]['1A-strand'] = feature.location.strand
                                        cds_annot_count +=1
                                        cds_count +=1

                    else:
                        print("No annotations of CDS in {}".format(rec_name))
                        #print(rec.description)
                        #print(feature)
                        no_annot_CDS +=1


            if cds_annot_count == 0 or cds_count == 0:
                
                # Check exception files. I comment this part because I hope that future pipeline will determine the ORF in such exceptions automatically
                '''
                if rec.name in exceptions_ORF1a:
                    dict_coord[rec.name]["1A"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                    dict_coord[rec.name]['1A-strand'] = feature.location.strand
                    cds_annot_count +=1
                    cds_count +=1
                    print("{} was found in ORF1a exceptions".format(rec.name))
                    continue
                if rec.name in exceptions_ORF1b:
                    dict_coord[rec.name]["1B"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                    dict_coord[rec.name]['1B-strand'] = feature.location.strand
                    cds_annot_count +=1
                    cds_count +=1
                    print("{} was found in ORF1b exceptions".format(rec.name))
                    continue
                if rec.name in exceptions_ORF2:
                    dict_coord[rec.name]["2"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                    dict_coord[rec.name]['2-strand'] = feature.location.strand
                    cds_annot_count +=1
                    cds_count +=1
                    print("{} was found in ORF2 exceptions".format(rec.name))
                    continue
                '''
                # Check the keywords in DEFINITION
                if len(rec) < 2700:
                    ORFs = map_keywords(rec.description,description_keywords)
                    if len(ORFs) == 1:
                        ORF_code = ORFs[0]
                        for feature in rec.features:
                            if feature.type == 'misc_feature':
                                dict_coord[rec_name][ORF_code] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                dict_coord[rec_name][ORF_code + '-strand'] = feature.location.strand
                                cds_annot_count +=1
                                cds_count +=1
                    elif len(ORFs) == 0:
                        misc_f_count = 0
                        for feature in rec.features:
                            if feature.type == 'misc_feature':
                                note = feature.qualifiers.get("note")[0]
                                ORFs = map_keywords(note,description_keywords)
                                if len(ORFs) !=0:
                                    ORF_code = ORFs[0]
                                    dict_coord[rec_name][ORF_code] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                    dict_coord[rec_name][ORF_code + '-strand'] = feature.location.strand
                                    misc_f_count +=1
                                    cds_annot_count +=1
                                    cds_count +=1
                        if misc_f_count == 0:
                            print("No keywords were found in description for {}".format(rec_name))
                    else:
                        print("Conflicting keywords were found in description for {}".format(rec_name))
                        print(ORFs)
                else:
                    # Check misc features for annotations
                    for feature in rec.features:
                        if feature.type == "misc_feature":
                            miscfeature_raw_annotations = []
                            for el in ('product', 'gene', 'note'):
                                miscfeature_raw_annotations.extend(feature.qualifiers.get(el, []))
                            #print(miscfeature_raw_annotations)
                            
                            miscfeature_annotations = []
                            for f in miscfeature_raw_annotations:
                                ORFs = map_keywords(f,description_keywords)
                                miscfeature_annotations.extend(ORFs)
                            miscfeature_annotations = list(set(miscfeature_annotations))
                            if len(miscfeature_annotations) != 1:
                                '''
                                if rec.name in ["PV793828","OR871061","OR871062"]:
                                    dict_coord[rec.name]["1A"] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                    dict_coord[rec.name]['1A-strand'] = feature.location.strand
                                    cds_annot_count +=1
                                    cds_count +=1
                                else:
                                '''
                                print("Conflicting annotations in misc_feature for {}".format(rec_name))
                                print(miscfeature_raw_annotations)
                                print(miscfeature_annotations)
                            else:
                                ORF_code = miscfeature_annotations[0]
                                dict_coord[rec_name][ORF_code] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                dict_coord[rec_name][ORF_code + '-strand'] = feature.location.strand
                                cds_annot_count +=1
                                cds_count +=1
            
            
            if cds_count == 0:
                no_annot +=1
                entries_noannot.append(rec.id + ':' + str(len(rec)) + ':' + rec.description + '\n')
                no_cds_accessions.append(rec.id)
           
            validate_and_correct_record(rec, dict_coord[rec_name])



        # Write extracted coordinates to file
        with open(out_file_name, 'w') as out_file:
            line = "Accession"
            for ORF in orf_types_final:
                line = line + ',' + ORF + ',' + ORF + '-strand'
            out_file.write(line + '\n')
            for rec_id, values in dict_coord.items():
                s = rec_id
                for orf in orf_types_final:
                    if orf in values:
                        s += f",{values[orf][0]}-{values[orf][1]},{values[orf+'-strand']}"
                    else:
                        s += ",NA-NA,NA"
                s += "\n"
                out_file.write(s)
        out_file.close()

        # Write all ORFs with determined coordinates to file
        fasta_out = out_file_name_temp + '_orfs.faa'
        export_orfs_to_fasta(records, dict_coord, fasta_out)

        if problematic_candidates:
        # Remove duplicates that are already in coord_dict
            unique_candidates = [c for c in problematic_candidates if not is_duplicate_candidate(c, dict_coord)]
            if unique_candidates:
                fasta_cand = out_file_name_temp + '_problematic_candidates.faa'
                mapping_cand = out_file_name_temp + '_problematic_candidates.tsv'
                write_fasta_and_mapping(unique_candidates, fasta_cand, mapping_cand)
            else:
                print("All problematic candidates are already represented in main ORFs.")

        # Write sequences with no annotated CDS to file in fasta format
        fasta_no_cds = out_file_name_temp + '_no_cds_sequences.fasta'
        with open(fasta_no_cds, 'w') as f:
            for rec in records:
                if rec.id in no_cds_accessions:
                    f.write(f">{rec.id} {rec.description}\n")
                    # Wrap sequence in lines of 60 characters
                    seq = str(rec.seq)
                    for i in range(0, len(seq), 60):
                        f.write(seq[i:i+60] + '\n')

        
        print("Entries with no CDS annotation: {}".format(no_annot))
        print("CDS with no target ORFs: {}".format(notarget_annot_CDS))
        print("CDS with no annotations: {}".format(no_annot_CDS))

        # Write information about entries with no annoyated CDS
        with open(os.path.join(log_dir, base_name + '_noannotCDS.txt'), 'w') as file:
            file.writelines(entries_noannot)

        with open(os.path.join(log_dir, base_name + '_entries_with_notarget_cds.txt'), 'w') as f:
            for acc in sorted(entries_with_notarget_cds):
                f.write(acc + '\n')

        with open(os.path.join(log_dir, base_name + '_conflictingannot.txt'), 'w') as file:
            file.writelines(entries_conflicting_annotations)

        #with open(out_file_name_temp + '_onlyNSP_long.txt', 'w') as file:
        #    file.writelines(entries_onlynonstrprotein)
        #    
        #with open(out_file_name_temp + '_onlyNSP_short.txt', 'w') as file:
        #    file.writelines(entries_onlynonstrprotein_short)

        
        return 0