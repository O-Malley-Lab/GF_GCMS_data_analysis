"""
GF GCMS Data Analysis Script 1
Lazarina Butkovich 7/19/24

This script compiles outputs from different GC-MS data analysis tools to create a combined data excel to reference.

Tool outputs:
MS-DIAL
GNPS

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
TEMP_OVERALL_FOLDER = r'temp'

# Key column to keep consistent across datasets, unless otherwise specified
KEY_COL = 'shared name'

# MS-DIAL output to GNPS. Note that we want to include the key column in the columns to keep, becasue we will use MSDIAL_OUTPUT table as the base table for the summary table. For all additional tables, we do not want to keep adding the key columns.
MSDIAL_OUTPUT_FILENAME = 'MSDIAL_output.xlsx'
MSDIAL_OUTPUT_COLS_TO_KEEP = [KEY_COL, 'Average Rt(min)', 'Metabolite name MSDIAL', 'SMILES MSDIAL']
MSDIAL_OUTPUT_COLS_TO_KEEP_WEIGHT_NORM = ['Average Rt(min)', 'Metabolite name MSDIAL', 'SMILES MSDIAL', 'Weight Normalized Peak Area', 'OMALL_RFS_AR_S4_1_M','OMALL_RFS_AR_S4_2_M','OMALL_RFS_AR_S4_3_M','OMALL_RFS_AR_S4_4_M', 'OMALL_RFS_CC1_M', 'OMALL_RFS_CC2_M', 'OMALL_RFS_CC3_M', 'OMALL_RFS_CC4_M']

# GNPS output for non-singletons
GNPS_NODE_TABLE_FILENAME = 'GNPS_node_table.xlsx'
GNPS_NODE_TABLE_COLS_TO_KEEP = ['Compound_Name', 'MQScore', 'precursor mass', 'RTMean','GNPSGROUP:CC', 'GNPSGROUP:AR', 'GNPSGROUP:MC', 'GNPSGROUP:RF', 'GNPSGROUP:BLANK', 'GNPSGROUP:FAMES', 'Smiles', 'Compound_Source', 'Data_Collector', 'Instrument', 'INCHI', 'GNPSLibraryURL']
# GNPS outputs for all library hits including singletons; singletons without library hits are excluded by GNPS
GNPS_ALL_LIB_MATCHES_FILENAME = 'GNPS_all_lib_matches.xlsx'
GNPS_ALL_LIB_MATCHES_COLS_TO_KEEP = ['molecular_formula', 'npclassifier_superclass', 'npclassifier_class', 'npclassifier_pathway', 'Compound_Name', 'MQScore', 'precursor mass', 'RTMean','GNPSGROUP:CC', 'GNPSGROUP:AR', 'GNPSGROUP:MC', 'GNPSGROUP:RF', 'GNPSGROUP:BLANK', 'GNPSGROUP:FAMES', 'Smiles', 'Compound_Source', 'Data_Collector', 'Instrument', 'INCHI', 'GNPSLibraryURL']
KEY_COL_GNPS_LIB_MATCHES = 'Scan_num'

# Compound matches from FienLib + NIST 13 (from PNNL)
CMPD_IDS_PNNL_FILENAME = 'Compound_ids_PNNL.xlsm'
CMPD_IDS_PNNL_COLS_TO_KEEP = ['cmpd_id_nist', 'Metabolite', 'Kegg ID', 'Metabolite Class','Confidence']
# Note, at the moment, no key column to match to other data...

# Cell pellet weight data for direct comparison of CC and AR relative peak intensities
CELL_PELLET_WEIGHTS_FILENAME = 'GF_cell_pellet_weights.xlsx'
# Column names: 'Sample', 'Sample Mass mg'

"""
Import data tables
"""
# GNPS output for non-singletons
gnps_node_table = pd.read_excel(pjoin(INPUT_FOLDER, GNPS_NODE_TABLE_FILENAME))
# GNPS outputs for all library hits including singletons; singletons without library hits are excluded by GNPS
gnps_all_lib_hits = pd.read_excel(pjoin(INPUT_FOLDER, GNPS_ALL_LIB_MATCHES_FILENAME))

# Compound matches from FienLib + NIST 13 (from PNNL)
cmpd_ids_pnnl = pd.read_excel(pjoin(INPUT_FOLDER, CMPD_IDS_PNNL_FILENAME))

# Cell pellet weight data for direct comparison of CC and AR relative peak intensities
cell_pellet_weights = pd.read_excel(pjoin(INPUT_FOLDER, CELL_PELLET_WEIGHTS_FILENAME))

# MS-DIAL output to GNPS
msdial_output = pd.read_excel(pjoin(INPUT_FOLDER, MSDIAL_OUTPUT_FILENAME))


"""
Compile Summary Data Table
"""
# Use MSDIAL output as base table because it has 1 row per feature (720 total). Remove index. Keep only indicated columns.
summary_table = msdial_output[MSDIAL_OUTPUT_COLS_TO_KEEP].copy()

# Add GNPS node table data
combine_dfs(summary_table, gnps_node_table, GNPS_NODE_TABLE_COLS_TO_KEEP, KEY_COL, KEY_COL, 'str')




