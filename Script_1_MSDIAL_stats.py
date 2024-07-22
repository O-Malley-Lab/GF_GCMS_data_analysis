"""
GF GCMS Data Analysis Script 1
Lazarina Butkovich 7/19/24

This script uses MS-DIAL outputs for GC-MS analysis to generate statistics to describe significant differences of metabolite amounts between samples. Two normalization methods will be used:
(1) Normalize CC and AR peak area values based on cell pellet weight data. This script takes the peak area values from 'MSDIAL_area_output.xlsx' and the cell pellet data from 'GF_cell_pellet_weights.xlsx' to do this normalization.
(2) Normalize all peak area values based on TIC sum (already done by MS-DIAL; these are the data in 'MSDIAL_norm_TIC_output.xlsx')

As output, this script exports a clean excel file with relevant statistics:
(1) Using cell pellet weight normalization, generate t test p-values to describe which metabolites are significantly present in CC and not AR (p_val_CC_vs_AR)
(2) Using TIC sum normalization, generate t test p-values to describe which metabolites are significantly present in CC and not MC (p_val_CC_vs_MC), and which metabolites are significantly present in AR and not MC (p_val_AR_vs_MC).


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

# MSDIAL output file with TIC normalized data
FILENAME_MSDIAL_OUTPUT_NORM_TIC = 'MSDIAL_norm_TIC_output.xlsx'

# MSDIAL output file with peak area data
FILENAME_MSDIAL_OUTPUT_AREA = 'MSDIAL_area_output.xlsx'

# Cell pellet weight data for direct comparison of CC and AR relative peak intensities
FILENAME_CELL_PELLET_WEIGHTS = 'GF_cell_pellet_weights.xlsx'
# Column names: 'Sample', 'Sample Mass mg'

# Dictionary of tuples to describe pre- and post- strings in sample names (the middle part is 1 to n, where n=3 or 4 for biological (for AR and CC) or technical (for BLANK, MC, RF) replicates)
SAMPLE_NAME_PRE_POST_STRS_DICT = {'CC':('OMALL_RFS_CC','_M'),'AR':('OMALL_RFS_AR_S4_','_M'),'MC':('OMALL_RFS_MC','_M'),'RF':('OMALL_RFS_RF','_M'), 'FAMES':('GCMS_FAMES_0','_GCMS01_20201209'), 'BLANK':('GCMS_BLANK_0','_GCMS01_20201209')}
REPLICATE_NUMS = {'CC':3, 'AR':3, 'MC':4, 'RF':4,'FAMES':1,'BLANK':3}

OUTPUT_FILENAME = 'MSDIAL_stats.xlsx'


"""""""""""""""""""""""""""""""""""""""""""""
Main
"""""""""""""""""""""""""""""""""""""""""""""
"""
Import data tables
"""
# Import data tables
df_msdial_norm_tic = pd.read_excel(pjoin(INPUT_FOLDER, FILENAME_MSDIAL_OUTPUT_NORM_TIC))
df_msdial_area = pd.read_excel(pjoin(INPUT_FOLDER, FILENAME_MSDIAL_OUTPUT_AREA))
df_cell_pellet_weights = pd.read_excel(pjoin(INPUT_FOLDER, FILENAME_CELL_PELLET_WEIGHTS))


"""
Generate Sample Name Columns Based on Templates
"""
# See SAMPLE_NAME_PRE_POST_STRS_DICT and REPLICATE_NUMS for how to assemble sample names. 
# Example: OMALL_RFS_CC1_M, OMALL_RFS_CC2_M, OMALL_RFS_CC3_M, OMALL_RFS_CC4_M

sample_cols = []
sample_groups_dict = {}
for key in SAMPLE_NAME_PRE_POST_STRS_DICT:
    pre_str = SAMPLE_NAME_PRE_POST_STRS_DICT[key][0]
    post_str = SAMPLE_NAME_PRE_POST_STRS_DICT[key][1]
    replicate_total = REPLICATE_NUMS[key]
    for i in range(1, replicate_total+1):
        col_name = pre_str + str(i) + post_str
        sample_cols.append(col_name)
        # Add dictinoary entry for key as key and list of sample names as values
        if key in sample_groups_dict:
            sample_groups_dict[key].append(col_name)
        else:
            sample_groups_dict[key] = [col_name]


"""
Normalize Peak Area Data Based on Cell Pellet Data
"""

"""
 To-do:

1. import cell pellet data
2. import peak area data
3. normalize peak area data based on cell pellet data
4. generate statistics for CC vs AR,
5. import TIC normalized data
6. generate statistics CC vs MC, AR vs MC
7. Organize statistics in a clean excel file
8. Export excel file and include an excel tab with filtered data
- filters: (p_val_sig = 0.05)

    a) p_val_CC_vs_MC < 0.05 --> metabolites significantly present in CC and not MC
    b) p_val_AR_vs_MC < 0.05 --> metabolites significantly present in AR and not MC
    c) p_val_CC_vs_AR < 0.05, CC_avg_cp_norm > AR_avg_cp_norm --> metabolites significantly more present in CC than AR
    d) p_val_CC_vs_AR < 0.05, AR_avg_cp_norm < AR_avg_cp_norm --> metabolites significantly more present in AR than CC
    
"""


