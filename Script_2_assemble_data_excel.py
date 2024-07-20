"""
GF GCMS Data Analysis Script 1
Lazarina Butkovich 7/19/24

This script compiles outputs from different GC-MS data analysis tools to create a combined data excel to reference.

Tool outputs:
MS-DIAL
GNPS

***Prior to using output data tables:
MS-DIAL: re-format so that the table does not have the top rows that are inconsistent with the rest of the format. Rename the average and standard deviation columns for samples to prevent them from having the same name (ie: rewrite as 'AR_avg' and'AR_std" instead of 'AR' and 'AR').

"""

import pandas as pd
import numpy as np
from os.path import join as pjoin


"""""""""""""""""""""""""""""""""""""""""""""
Functions
"""""""""""""""""""""""""""""""""""""""""""""
def check_key_col(df1, df2, col_name1, col_name2, key_type):
    """
    Check that the key columns are in the dataframes. If not, print an error message and exit the script. Change the key column type to match the other dataframe if necessary.
    input:
    df1: DataFrame 1
    df2: DataFrame 2
    col_name1: Name of the key column in df1
    col_name2: Name of the key column in df2
    key_type: 'int' or 'str'

    output: None
    """
    if col_name1 not in df1.columns:
        print("Error: no " + col_name1 + " column in main dataframe")
        exit()

    if col_name2 not in df2.columns:
        print("Error: no " + col_name2 + " column in dataframe to add cols from")
        exit()
    
    if key_type == 'int':
        df1[col_name1] = df1[col_name1].astype(int)
        df2[col_name2] = df2[col_name2].astype(int)
    elif key_type == 'str':
        df1[col_name1] = df1[col_name1].astype(str)
        df2[col_name2] = df2[col_name2].astype(str)
    else:
        print("Error: key_type must be 'int' or 'str'")
    return

def check_cols_to_add(df_main, df_add_cols_to_keep, key_col_main, key_col_add):
    """
    Checks that the key columns are not in df_add_cols_to_keep. If they are, remove them. Print a warning that previous values in df_main will be replaced with values from df_add for any matching cols in df_add_cols_to_keep.

    Inputs
    df_main: DataFrame to add columns to
    df_add_cols_to_keep: List of column names to add to df_main
    key_col_main: Name of the key column in df_main
    key_col_add: Name of the key column in df_add

    Outputs
    return: None
    """
    # If key column is in df_add_cols_to_keep, remove it
    if key_col_add in df_add_cols_to_keep:
        df_add_cols_to_keep.remove(key_col_add)
    if key_col_main in df_add_cols_to_keep:
        df_add_cols_to_keep.remove(key_col_main)
    # Print warning statement if any column names from df_add_cols_to_keep are already in df_main. This function will replace the previous values in df_main with the values from df_add.
    for col in df_add_cols_to_keep:
        if col in df_main.columns:
            print("Warning: " + col + " column already in df_main. This function is replacing the previous values in df_main with the values from df_add")
    return

def add_columns_to_df(df_main, df_add, df_add_cols_to_keep, key_col_main, key_col_add):
    """
    Function used by combine_dfs to add columns from df_add to df_main.

    Inputs
    df_main: DataFrame to add columns to
    df_add: DataFrame to add columns from
    df_add_cols_to_keep: List of column names to add to df_main
    key_col_main: Name of the key column in df_main
    key_col_add: Name of the key column in df_add

    Outputs
    return: None    
    """
    # Prepare df_main for merging by adding new blank columns for the columns to add from df_add
    for col in df_add_cols_to_keep:
        df_main[col] = np.nan

    # Iterate through each row of key values in df_main
    for index, row in df_main.iterrows():
        key_val = row[key_col_main]
        # Check that the key value is in df_add
        if key_val not in df_add[key_col_add].values:
            # skip
            continue
        else:
            # Get the row from df_add that matches the key value
            add_row = df_add.loc[df_add[key_col_add] == key_val]
            # Add the columns from df_add to df_main
            for col in df_add_cols_to_keep:
                df_main.at[index, col] = add_row[col].values[0]
    return

def combine_dfs(df_main, df_add, df_add_cols_to_keep, key_col_main, key_col_add, key_type = 'int'):
    """
    Main function to add columns from df_add to df_main.
    
    Use key_col_main for df_main and key_col_add for df_add to match rows. Need to iterate through each row of key values in df_main since the key values may not be in the same order in df_main and df_add, and not all key values from one table will be in another table.

    Inputs
    df_main: DataFrame to add columns to
    df_add: DataFrame to add columns from
    df_add_cols_to_keep: List of column names to add to df_main
    key_col_main: Name of the key column in df_main
    key_col_add: Name of the key column in df_add
    key_type: 'int' or 'str', desired data type for the key columns
    
    Outputs
    return: None
    """
    # Change the key column types between the dataframes
    check_key_col(df_main, df_add, key_col_main, key_col_add, key_type)

    # Check columns to keep and key columns
    check_cols_to_add(df_main, df_add_cols_to_keep, key_col_main, key_col_add)    

    # Add columns to df_main
    add_columns_to_df(df_main, df_add, df_add_cols_to_keep, key_col_main, key_col_add)

    return

def write_table_to_excel(writer, df, sheet_name):
    """
    Write a dataframe to an excel sheet. The column width will be set to the size of the header text.

    Inputs
    writer: ExcelWriter object
    df: DataFrame to write to the sheet
    sheet_name: Name of the sheet to write to

    Outputs
    return: None
    """
    df.to_excel(writer, sheet_name = sheet_name, index = False)
    # Format the excel sheets so that the column width matches the size of the header text
    workbook = writer.book
    worksheet = writer.sheets[sheet_name]
    for idx, col in enumerate(df):  # loop through all columns
        series = df[col]
        # Set max_len to the length of only the header text
        max_len = len(str(series.name)) + 1
        worksheet.set_column(idx, idx, max_len)  # set column width
    return
   
def format_column(worksheet, df):
    """
    Format excel sheet column width to match the size of the header text.

    Inputs
    worksheet: ExcelWriter worksheet object
    df: DataFrame to format

    Outputs
    return: None
    """
    for idx, col in enumerate(df):  # loop through all columns
        series = df[col]
        # Set max_len to the length of only the header text
        max_len = len(str(series.name)) + 1
        worksheet.set_column(idx, idx, max_len)  # set column width
        # Make the top header "sticky" so that it is always visible
        worksheet.freeze_panes(1, 0)
    return


"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'

# Key column to keep consistent across datasets, unless otherwise specified
KEY_COL = 'shared name'

# MS-DIAL output to GNPS. Note that we want to include a key column in the columns to keep, becasue we will use MSDIAL_OUTPUT table as the base table for the summary table. For all additional tables, we do not want to keep adding the key columns. The average peak intensities were normalized by TIC sum. 'Alignment ID' will get converted to a 'shared name' column by adding 1 to the values.
MSDIAL_OUTPUT_FILENAME = 'MSDIAL_norm_TIC_output.xlsx'
MSDIAL_OUTPUT_COLS_TO_KEEP = ['shared name', 'Average Rt(min)', 'Metabolite name', 'SMILES', 'BLANK_avg', 'FAMES_avg', 'AR_avg', 'CC_avg', 'MC_avg', 'RF_avg']

# GNPS outputs for all library hits including singletons; singletons without library hits are excluded by GNPS
GNPS_ALL_LIB_MATCHES_FILENAME = 'GNPS_all_lib_matches.xlsx'
GNPS_ALL_LIB_MATCHES_COLS_TO_KEEP = ['Compound_Name', 'MQScore', 'Precursor_MZ', 'molecular_formula', 'npclassifier_superclass', 'npclassifier_class', 'npclassifier_pathway', 'Smiles', 'Compound_Source', 'Data_Collector', 'Instrument', 'INCHI']
KEY_COL_GNPS_LIB_MATCHES = 'Scan_num'

# Compound matches from FienLib + NIST 13 (from PNNL)
CMPD_IDS_PNNL_FILENAME = 'Compound_ids_PNNL.xlsm'
CMPD_IDS_PNNL_COLS_TO_KEEP = ['cmpd_id_nist', 'Metabolite', 'Kegg ID', 'Metabolite Class','Confidence']
# Note, at the moment, no key column to match to other data...

# Cell pellet weight data for direct comparison of CC and AR relative peak intensities
CELL_PELLET_WEIGHTS_FILENAME = 'GF_cell_pellet_weights.xlsx'
# Column names: 'Sample', 'Sample Mass mg'

OUTPUT_FILENAME = 'GF_GCMS_summary_table_temp.xlsx'

FINAL_COLS_ORDER = ['shared name','Average Rt(min)', 'Precursor_MZ', 'Compound_Name','MQScore', 'Smiles', 'INCHI', 'Metabolite name', 'SMILES','molecular_formula', 'npclassifier_superclass', 'npclassifier_class', 'npclassifier_pathway', 'Compound_Source', 'Data_Collector', 'Instrument', 'BLANK_avg', 'BLANK_avg_log10', 'FAMES_avg', 'FAMES_avg_log10', 'AR_avg', 'AR_avg_log10', 'CC_avg', 'CC_avg_log10', 'MC_avg', 'MC_avg_log10', 'RF_avg', 'RF_avg_log10']

# FINAL_COLS_ORDER_SIMPLE = ['shared name','Average Rt(min)', 'Precursor_MZ', 'Compound_Name','MQScore', 'Smiles', 'Metabolite name', 'SMILES','molecular_formula', 'npclassifier_superclass', 'npclassifier_class', 'npclassifier_pathway', 'BLANK_avg', 'BLANK_avg_log10', 'FAMES_avg', 'FAMES_avg_log10', 'AR_avg', 'AR_avg_log10', 'CC_avg', 'CC_avg_log10', 'MC_avg', 'MC_avg_log10', 'RF_avg', 'RF_avg_log10']

FINAL_COLS_ORDER_SIMPLE = ['shared name','Average Rt(min)', 'Precursor_MZ', 'Compound_Name','MQScore', 'Smiles','molecular_formula', 'npclassifier_superclass', 'npclassifier_class', 'npclassifier_pathway', 'BLANK_avg', 'BLANK_avg_log10', 'FAMES_avg', 'FAMES_avg_log10', 'AR_avg', 'AR_avg_log10', 'CC_avg', 'CC_avg_log10', 'MC_avg', 'MC_avg_log10', 'RF_avg', 'RF_avg_log10']


COLS_NAME_CONVERTER = {'Average Rt(min)':'RT', 'Precursor_MZ':'EI spectra quant mass', 'Compound_Name':'Compounds_Name_GNPS','MQScore':'MQScore_GNPS', 'Smiles':'SMILES_GNPS', 'Metabolite name': 'Metabolite name MSDIAL', 'SMILES':'SMILES MSDIAL', 'INCHI':'INCHI_GNPS', 'molecular_formula':'molecular_formula_GNPS', 'npclassifier_superclass':'npclassifier_superclass_GNPS', 'npclassifier_class':'npclassifier_class_GNPS', 'npclassifier_pathway':'npclassifier_pathway_GNPS','Compound_Source':'Compound Source GNPS', 'Data_Collector':'Data Collector GNPS', 'Instrument':'Instrument_GNPS'}


"""""""""""""""""""""""""""""""""""""""""""""
Main
"""""""""""""""""""""""""""""""""""""""""""""
"""
Import data tables
"""
# MS-DIAL output to GNPS
msdial_output = pd.read_excel(pjoin(INPUT_FOLDER, MSDIAL_OUTPUT_FILENAME))

# GNPS outputs for all library hits including singletons; singletons without library hits are excluded by GNPS
gnps_all_lib_hits = pd.read_excel(pjoin(INPUT_FOLDER, GNPS_ALL_LIB_MATCHES_FILENAME))

# Compound matches from FienLib + NIST 13 (from PNNL)
cmpd_ids_pnnl = pd.read_excel(pjoin(INPUT_FOLDER, CMPD_IDS_PNNL_FILENAME))

# Cell pellet weight data for direct comparison of CC and AR relative peak intensities
cell_pellet_weights = pd.read_excel(pjoin(INPUT_FOLDER, CELL_PELLET_WEIGHTS_FILENAME))

"""
For script 2 use: Add shared name key column to MS-DIAL output table and export to TEMP folder
"""
# Create the shared name column and data values, where the values are 1 plus the 'Alignment ID' values
msdial_output[KEY_COL] = msdial_output['Alignment ID'] + 1

# Move KEY_COL to the first column
cols = msdial_output.columns.tolist()
cols = cols[-1:] + cols[:-1]
msdial_output = msdial_output[cols]

# For values in the 'Metabolite name MSDIAL column', if the value is "Unknown", change to a blank value
msdial_output.loc[msdial_output['Metabolite name'] == 'Unknown', 'Metabolite name'] = ''

# Set aside table to export to TEMP folder
msdial_output_export = msdial_output.copy()

# Use COLS_NAME_CONVERTER to rename the columns
msdial_output_export.rename(columns=COLS_NAME_CONVERTER, inplace=True)

# Export to TEMP folder
msdial_output.to_excel(pjoin(TEMP_FOLDER, 'MSDIAL_output_updated.xlsx'), index = False)


"""
Initialize Summary Data Table
"""
# Use MSDIAL output as base table because it has 1 row per feature (720 total). Remove index. Keep only indicated columns.
summary_table = msdial_output[MSDIAL_OUTPUT_COLS_TO_KEEP].copy()


"""
Filter GNPS All Library Hits Table for Best Matches and Add to Summary Data Table
"""
# Add GNPS all library hits data
# Before combining this table, we need to filter gnps_all_lib_hits for the best compound matches for any given feature. We will use the MQScore to determine the best match.

# Filter gnp_all_lib_hits for the best compound match for each feature
gnps_all_lib_hits_best_match = gnps_all_lib_hits.loc[gnps_all_lib_hits.groupby(KEY_COL_GNPS_LIB_MATCHES)['MQScore'].idxmax()]

# Combine the filtered table with the summary table
combine_dfs(summary_table, gnps_all_lib_hits_best_match, GNPS_ALL_LIB_MATCHES_COLS_TO_KEEP, KEY_COL, KEY_COL_GNPS_LIB_MATCHES, 'str')


"""
Add log10 average peak intensity columns
"""
# Create log10 average peak intensity columns. For original values of 0, set the log10 value to nan, to avoid divide by 0 warning
for col in ['BLANK_avg', 'FAMES_avg', 'AR_avg', 'CC_avg', 'MC_avg', 'RF_avg']:
    for index, row in summary_table.iterrows():
        if row[col] == 0:
            summary_table.at[index, col + '_log10'] = np.nan
        else:
            summary_table.at[index, col + '_log10'] = np.log10(row[col])
            # Round values to 2 decimal places
            summary_table.at[index, col + '_log10'] = round(summary_table.at[index, col + '_log10'], 2)


"""
Format Summary Data Table
"""
# Convert the key column to int, since shared name is a number
summary_table[KEY_COL] = summary_table[KEY_COL].astype(int)

# Filter the summary table to only include the columns in FINAL_COLS_ORDER
summary_table = summary_table[FINAL_COLS_ORDER]

# Create a simple copy of the summary table with the columns in FINAL_COLS_ORDER_SIMPLE
summary_table_simple = summary_table[FINAL_COLS_ORDER_SIMPLE].copy()

# Convert the column names to the new names using COLS_NAME_CONVERTER
# in summary_table 
summary_table.rename(columns=COLS_NAME_CONVERTER, inplace=True)
# in summary_table_simple
summary_table_simple.rename(columns=COLS_NAME_CONVERTER, inplace=True)


"""
Export to Excel
"""
# Write results to excel
writer = pd.ExcelWriter(pjoin(TEMP_FOLDER, OUTPUT_FILENAME), engine='xlsxwriter')

# write summary_table
write_table_to_excel(writer, summary_table, 'Summary_Table')
workbook = writer.book
worksheet = writer.sheets['Summary_Table']
format_column(worksheet, summary_table)

# write summary_table_simple
write_table_to_excel(writer, summary_table_simple, 'Summary_Table_Simple')
workbook = writer.book
worksheet = writer.sheets['Summary_Table_Simple']
format_column(worksheet, summary_table_simple)

writer.close()






