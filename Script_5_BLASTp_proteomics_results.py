"""
GF GCMS Data Analysis Script 5: BLASTp to align proteomics ProteinNames to Mycocosm proteinIDs and Annotations
Lazarina Butkovich 1/7/25

- Perform BLASTp on proteomics results vs Mycocosm aa fasta file to align ProteinName (from proteomics data) with Mycocosm proteinID (from fasta file)
- Align functional annotations to results

*** Note, run in 2 parts. First part prepares files for commandline BLASTp. The second part aligns results with annotations.

"""

import pandas as pd
from os.path import join as pjoin

"""
Functions
"""
def parse_fasta_file(f):
    """Return a dict of {id:gene_seq} pairs based on the sequences in the input FASTA file
    input_file -- a file handle for an input fasta file
    """
    parsed_seqs = {}
    curr_seq_id = None
    curr_seq = []

    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if curr_seq_id is not None:
                parsed_seqs[curr_seq_id] = ''.join(curr_seq)

            curr_seq_id = line[1:]
            curr_seq = []
            continue

        curr_seq.append(line)

    #Add the final sequence to the dict
    parsed_seqs[curr_seq_id] = ''.join(curr_seq)
    return parsed_seqs

def clean_ID_column_FASTA(dict_fasta):
    """ 
    Input a dictionary from a parsed FASTA file and output the dictionary with the keys formatted a specific way. In this case, the dict_fasta keys are formatted as 'jgi|Neosp1|$|...', where $ is the only value we want to keep. Take the value between the 2nd and 3rd '|' and use that as the new key.
    """
    new_dict = {}
    for k,v in dict_fasta.items():
        # Split the key by '|'
        split_key = k.split('|')
        # Get the value between the 2nd and 3rd '|'
        new_key = split_key[2]
        # Add the new key and value to the new dictionary
        new_dict[new_key] = v
    return new_dict

def check_ID_column(pd_df, id_name="proteinID"):
    """
    The function check_ID_column checks if the ID column id_name is present in the pandas dataframe, pd_df.

    inputs:
    pd_df: pandas dataframe
    id_name: name of id column, default set to 'proteinID' (string)

    output:
    boolean: True if id_name column in pd_df, exit if not
    """
    if id_name in pd_df.columns:
        return
    else:
        print("Error: no expected " + id_name + " column in pd_df")
        exit()

def make_proteinID_annot_dict(pd_annot, annot_col, id_name="proteinID"):
    """
    The function make_proteinID_annot_dict makes a dictionary with proteinIDs as keys and annotations as values from the annot_col column of the pandas dataframe pd_annot. The function removes duplicate annotations for a given ID. Multiple annotations for a single ID are all kept and separated by a comma in a single value string. Each proteinID key will therefore only have 1 value string. The purpose of this formatting is to prepare the data for addition to the summary dataframe, with the single string as a value in a column.

    inputs:
    pd_annot: pandas dataframe with annotations
    - must have id_name (string) column
    annot_col: specific column name (string) with annotations
    id_name: name of id column, default set to 'proteinID' (string)

    output:
    proteinID_annot_dict: dictionary with IDs as keys and annotations listed (sep=', ') as the values. Remove duplicate annotations for a given ID.
    """
    proteinID_annot_dict = {}
    # if no proteinID (or id_name) column in pd_annot, print error message
    check_ID_column(pd_annot, id_name)
    for i in range(len(pd_annot)):
        proteinID = pd_annot[id_name][i]
        annot = pd_annot[annot_col][i]
        if proteinID in proteinID_annot_dict:
            # if annot is different from proteinID_annot_dict[proteinID] ...
            # (I don't want duplicates of same annotation values)
            if annot not in proteinID_annot_dict[proteinID]:
                # ... then append annot to proteinID_annot_dict[proteinID]
                proteinID_annot_dict[proteinID].append(annot)
        else:
            proteinID_annot_dict[proteinID] = [annot]
    for key in proteinID_annot_dict:
        if len(proteinID_annot_dict[key]) > 1:
            vals = list(map(str, proteinID_annot_dict[key]))
            proteinID_annot_dict[key] = ",".join(vals)
        else:
            proteinID_annot_dict[key] = proteinID_annot_dict[key][0]
    return proteinID_annot_dict

def add_to_df(df, X_annot_pd, annot_cols_list, shared_col='proteinID'):
    """
    The function add_to_df adds annotation columns from the pandas dataframe X_annot_pd to the pandas dataframe df. The function returns df with the new columns. For example, df is the DGE_summary dataframe and X_annot_pd could be the KOG annotations dataframe.

    inputs: 
    df: pandas dataframe with proteinIDs
    X_annot_pd: pandas dataframe with proteinIDs and annotations
    annot_cols_list: list of column names in X_annot_pd to add to DGE_summary (list of strings)
    shared_col: column name in df and X_annot_pd that is shared (string) default set to proteinID

    output: df with new columns (pandas dataframe)
    """
    X_dicts_list = []
    for col in annot_cols_list:
        X_dicts_list.append(make_proteinID_annot_dict(X_annot_pd, col, shared_col))

    # Add annotations to df
    annot_col_num = 0
    for annot_dict in X_dicts_list:
        df[annot_cols_list[annot_col_num]] = df[shared_col].map(annot_dict)
        annot_col_num += 1

    return df

"""
Part 1: Prepare files for BLASTp
"""

"""
***Values to change for Part 1***
"""
INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'

# # Comment out the one not being used
# For CC:
STUDY_NAME = 'CC'
PROTEOMICS_DATA_FILENAME = 'EMSL_50386_OMall_RFS_Cc_FirstHitsResults.xlsx'
MYCOCOSM_FASTA_FILENAME = 'Caecom1_GeneCatalog_proteins_20171213.aa.fasta'

# # For AR:
# STUDY_NAME = 'AR'
# PROTEOMICS_DATA_FILENAME = 'EMSL_50386_OMall_RFS_Ar_S4_FirstHitsResults.xlsx'
# MYCOCOSM_FASTA_FILENAME = 'Anasp1_FilteredModels3_deflines.aa.fasta'


# For CC and AR:
PROTEOMICS_DATA_SHEETNAME = 'Protein Crosstab wSeqs'
PROTEOMICS_ID_COLUMN_NAME = 'ProteinName'
PROTEOMICS_SEQ_COLUMN_NAME = 'Sequence'


"""
Import data files
"""
# Parse the fasta file into a pd dataframe
fasta_data = open(pjoin(INPUT_FOLDER, MYCOCOSM_FASTA_FILENAME), 'r')
parsed_all_seqs = parse_fasta_file(fasta_data)
fasta_data.close()

# Import proteomics data into a pd dataframe. Use the given sheetname. The column PROTEOMICS_ID_COLUMN_NAME contains the protein names.
proteomics_data = pd.read_excel(pjoin(INPUT_FOLDER, PROTEOMICS_DATA_FILENAME), sheet_name=PROTEOMICS_DATA_SHEETNAME)


"""
Write the db file for BLASTp
"""
# Cleanup the id column (extract proteinID value), if CC STUDY_NAME
if STUDY_NAME == 'CC':
    parsed_all_seqs_cleaned = clean_ID_column_FASTA(parsed_all_seqs)
else :
    parsed_all_seqs_cleaned = parsed_all_seqs

# Write the cleaned fasta data to a new .fasta file in the temp folder
# note the .fasta file type
with open(pjoin(TEMP_FOLDER, MYCOCOSM_FASTA_FILENAME.removesuffix('.aa.fasta') + '_db_BLASTp.fasta'), 'w') as f:
    for k,v in parsed_all_seqs_cleaned.items():
        f.write('>' + k + '\n' + v + '\n')

"""
Write the query file for BLASTp
"""
# Write and simplify the proteomics data to .txt file in fasta format 
with open(pjoin(TEMP_FOLDER, PROTEOMICS_DATA_FILENAME.removesuffix('.xlsx') + '_query_BLASTp.txt'), 'w') as f:
    # Iterate through rows
    # Value in the PROTEOMICS_ID_COLUMN_NAME column is the protein name to go after each '>'
    # The value in the PROTEOMICS_SEQ_COLUMN_NAME column is the sequence to go after each '>'
    for index, row in proteomics_data.iterrows():
        f.write('>' + row[PROTEOMICS_ID_COLUMN_NAME] + '\n' + row[PROTEOMICS_SEQ_COLUMN_NAME] + '\n')





"""
***
Manual: run commandline BLASTp. Copy-paste the output into the input folder and run the below script to add the correct headers and align annotations.
***
"""




"""
***Values to change for Part 2***
"""
# # Comment out the one not being used
# For CC:
BLASTP_RESULTS = 'output_CC_proteomics_BLASTp_updated.txt'
MYCOCOSM_KOG_ANNOTATIONS_FILENAME = 'CC_KOG_annotations.xlsx'
MYCOCOSM_KOG_PROTEINID_COLUMN_NAME = 'proteinId'

# # For AR:
# BLASTP_RESULTS = 'output_AR_proteomics_BLASTp_updated.txt'
# MYCOCOSM_KOG_ANNOTATIONS_FILENAME = 'AR_KOG_annotations.xlsx'
# MYCOCOSM_KOG_PROTEINID_COLUMN_NAME = 'proteinId'


BLASTP_HEADER_NAMES = ['query_seq_id', 'subject_seq_id', 'percent_identity', 'alignment_length', 'Query Coverage Per Subject (for all HSPs)', 'Query Coverage Per HSP', 'mismatches', 'gap_openings', 'query_start', 'query_end', 'subject_start', 'subject_end', 'e_value', 'bit_score']
PIDENT_CUTOFF = 90
QUERY_COVERAGE_CUTOFF = 75


"""
Import BLASTp output data file
"""
# Adjust the sep delimiter as necessary, no index column
blastp_results = pd.read_csv(pjoin(INPUT_FOLDER, BLASTP_RESULTS), sep=',', header=None, names=BLASTP_HEADER_NAMES, index_col=False)

# Sort blastp_results by subject_seq_id and sub-sort by query_seq_id using merge sort
blastp_results = blastp_results.sort_values(by=['subject_seq_id', 'query_seq_id'], kind='mergesort')

# Reset index
blastp_results = blastp_results.reset_index(drop=True)


"""
Filter BLASTp Results
"""
# First filter by percent identity and query coverage thresholds
filtered_df = blastp_results[
    (blastp_results['percent_identity'] >= PIDENT_CUTOFF) & 
    (blastp_results['Query Coverage Per HSP'] >= QUERY_COVERAGE_CUTOFF)
].copy()

# Sort by quality metrics in priority order
filtered_df = filtered_df.sort_values(
    by=['query_seq_id', 'alignment_length', 'percent_identity', 'bit_score'],
    ascending=[True, False, False, False]
)

# Keep only the best match for each query sequence
blastp_results_filtered = filtered_df.drop_duplicates(subset=['query_seq_id'], keep='first')

# Verify no duplicates remain
if blastp_results_filtered['query_seq_id'].duplicated().any():
    print("Warning: Still found duplicate query sequences after filtering")
    # Additional cleanup if needed
    blastp_results_filtered = blastp_results_filtered.drop_duplicates(subset=['query_seq_id'], keep='first')

# Sort final results
blastp_results_filtered = blastp_results_filtered.sort_values(
    by=['query_seq_id', 'subject_seq_id'], 
    kind='mergesort'
).reset_index(drop=True)

# Make a new dataframe with just the needed columns
blastp_results_filtered = blastp_results_filtered[
    ['query_seq_id', 'subject_seq_id', 'percent_identity', 'Query Coverage Per HSP']
].copy()

"""
Align Annotations: Import Annotation File
"""
# Based on proteinID (subject_seq_id), align the annotations from the KOG annotations file to the blastp_results_filtered dataframe
kog_annotations = pd.read_excel(pjoin(INPUT_FOLDER, MYCOCOSM_KOG_ANNOTATIONS_FILENAME))


"""
Pandas Dataframes: Rename ID Columns to proteinID or ProteinName, as appropriate
"""
# Change the ID column name in kog_annotations to 'proteinID'
kog_annotations = kog_annotations.rename(columns={MYCOCOSM_KOG_PROTEINID_COLUMN_NAME: 'proteinID'})

# Change the subject_seq_id column name in blastp_results_filtered to 'proteinID'
blastp_results_filtered = blastp_results_filtered.rename(columns={'subject_seq_id': 'proteinID'})

# Change the query_seq_id column name in blastp_results_filtered to 'ProteinName'
blastp_results_filtered = blastp_results_filtered.rename(columns={'query_seq_id': 'ProteinName'})

"""
Use proteinID ID values to align additional columns to blastp_results_filtered: KOG Annotations
"""
# Append all columns from kog_annotations to blastp_results_filtered
annot_cols_list = kog_annotations.columns.tolist()
# remove proteinID column from annot_cols_list
annot_cols_list.remove('proteinID')
blastp_results_filtered = add_to_df(blastp_results_filtered, kog_annotations, annot_cols_list, 'proteinID')

"""
Use ProteinName ID values to align additional columns to blastp_results_filtered: proteomics data
"""
# Append all unused data columns from proteomics_data to blastp_results_filtered
proteomics_cols_list = proteomics_data.columns.tolist()
# remove ProteinName column from proteomics_cols_list
proteomics_cols_list.remove('ProteinName')
blastp_results_filtered = add_to_df(blastp_results_filtered, proteomics_data, proteomics_cols_list, 'ProteinName')

"""
Export the final dataframe to an excel file
"""
blastp_results_filtered.to_excel(pjoin(OUTPUT_FOLDER, STUDY_NAME + '_BLASTp_proteomics_results_final.xlsx'), index=False)

# print how many entries there are in the original fasta file
print("Number of entries in the original fasta file: ", len(parsed_all_seqs))