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
import matplotlib.pyplot as plt
# from bioinfokit import visuz


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

def color_excel_column(wksheet, df, col_name, min_val = 0.5, min_color = "#FFFFFF", max_val = 1, max_color = "#008000"):
    # Find the column index
    col_name_index = df.columns.get_loc(col_name)
    wksheet.conditional_format(1, col_name_index, len(df), col_name_index, {'type': '3_color_scale', 'min_type': 'num', 'min_value': min_val, 'min_color': min_color, 'max_type': 'num', 'max_value': max_val, 'max_color': max_color})
    return


"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'

# Key column to keep consistent across datasets, unless otherwise specified
KEY_COL = 'shared name'

# MS-DIAL output to GNPS. 'Alignment ID' was converted to  'shared name' column by adding 1 to the values in Script 1.
FILENAME_MSDIAL_OUTPUT = 'MSDIAL_stats.xlsx'

# GNPS outputs for all library hits including singletons; singletons without library hits are excluded by GNPS
FILENAME_GNPS_ALL_LIB_MATCHES = 'GNPS_all_lib_matches.xlsx'
COLS_TO_KEEP_GNPS_ALL_LIB_MATCHES = ['Compound_Name_GNPS', 'MQScore_GNPS', 'EI_spectra_quant_mass', 'molecular_formula_GNPS', 'npclassifier_superclass_GNPS', 'npclassifier_class_GNPS', 'npclassifier_pathway_GNPS', 'SMILES_GNPS', 'Compound_Source_GNPS', 'Data_Collector_GNPS', 'Instrument_GNPS', 'INCHI_GNPS']
KEY_COL_GNPS_LIB_MATCHES = 'Scan_num'

# Compound matches from FienLib + NIST 13 (from PNNL)
FILENAME_CMPD_IDS_PNNL = 'Compound_ids_PNNL.xlsm'
COLS_TO_KEEP_CMPD_IDS_PNNL = ['cmpd_id_nist', 'Metabolite', 'Kegg ID', 'Metabolite Class','Confidence']
# Note, at the moment, no key column to match to other data...

FILENAME_OUTPUT = 'GF_GCMS_stats_summary_table.xlsx'

FINAL_COLS_ORDER_SIMPLE = ['shared name', 'Alignment_ID_MSDIAL', 'RT', 'EI_spectra_quant_mass', 'Compound_Name_GNPS','MQScore_GNPS', 'SMILES_GNPS','molecular_formula_GNPS', 'npclassifier_superclass_GNPS', 'npclassifier_class_GNPS', 'npclassifier_pathway_GNPS', 'Metabolite_name_MSDIAL', 'SMILES_MSDIAL', 'Total_spectrum_similarity_MSDIAL',
'p_val_CC_vs_AR_cell_norm', 'p_val_CC_vs_MC','p_val_AR_vs_MC', 'p_val_CC_vs_BLANK', 'p_val_AR_vs_BLANK',
'CC_cell_norm_avg', 'CC_TIC_norm_avg', 'CC_TIC_norm_std', 'CC_avg_log10',
'AR_cell_norm_avg', 'AR_TIC_norm_avg', 'AR_TIC_norm_std', 'AR_avg_log10',
'MC_TIC_norm_avg', 'MC_TIC_norm_std', 'MC_avg_log10',
'RF_TIC_norm_avg', 'RF_TIC_norm_std', 'RF_avg_log10',
'FAMES_TIC_norm_avg', 'FAMES_TIC_norm_std', 'FAMES_avg_log10',
'BLANK_TIC_norm_avg', 'BLANK_TIC_norm_std', 'BLANK_avg_log10']

COLS_NAME_CONVERTER = {'Alignment ID': 'Alignment_ID_MSDIAL','Average Rt(min)':'RT', 'Precursor_MZ':'EI_spectra_quant_mass', 'Quant mass': 'Quant_mass', 'Compound_Name':'Compound_Name_GNPS','MQScore':'MQScore_GNPS', 'Smiles':'SMILES_GNPS', 'INCHI':'INCHI_GNPS', 'Metabolite name': 'Metabolite_name_MSDIAL', 'SMILES':'SMILES_MSDIAL', 'INCHI':'INCHI_GNPS', 'molecular_formula':'molecular_formula_GNPS', 'npclassifier_superclass':'npclassifier_superclass_GNPS', 'npclassifier_class':'npclassifier_class_GNPS', 'npclassifier_pathway':'npclassifier_pathway_GNPS','Compound_Source':'Compound_Source_GNPS', 'Data_Collector':'Data_Collector_GNPS', 'Instrument':'Instrument_GNPS', 'Total spectrum similarity': 'Total_spectrum_similarity_MSDIAL'}


"""""""""""""""""""""""""""""""""""""""""""""
Main
"""""""""""""""""""""""""""""""""""""""""""""
"""
Import data tables
"""
# MS-DIAL output to GNPS, read "Summary" tab
msdial_output = pd.read_excel(pjoin(TEMP_FOLDER, FILENAME_MSDIAL_OUTPUT),sheet_name='Summary')

# GNPS outputs for all library hits including singletons; singletons without library hits are excluded by GNPS
gnps_all_lib_hits = pd.read_excel(pjoin(INPUT_FOLDER, FILENAME_GNPS_ALL_LIB_MATCHES))
# Convert column names using COLS_NAME_CONVERTER
gnps_all_lib_hits.rename(columns=COLS_NAME_CONVERTER, inplace=True)

# Compound matches from FienLib + NIST 13 (from PNNL)
cmpd_ids_pnnl = pd.read_excel(pjoin(INPUT_FOLDER, FILENAME_CMPD_IDS_PNNL))


"""
Initialize Summary Data Table
"""
# Use MSDIAL output as base table because it has 1 row per feature (720 total). Remove index. Keep only indicated columns.
summary_table = msdial_output.copy()


"""
Filter GNPS All Library Hits Table for Best Matches and Add to Summary Data Table
"""
# Add GNPS all library hits data
# Before combining this table, we need to filter gnps_all_lib_hits for the best compound matches for any given feature. We will use the MQScore to determine the best match.

# Filter gnp_all_lib_hits for the best compound match for each feature
gnps_all_lib_hits_best_match = gnps_all_lib_hits.loc[gnps_all_lib_hits.groupby(KEY_COL_GNPS_LIB_MATCHES)['MQScore_GNPS'].idxmax()]

# Combine the filtered table with the summary table
combine_dfs(summary_table, gnps_all_lib_hits_best_match, COLS_TO_KEEP_GNPS_ALL_LIB_MATCHES, KEY_COL, KEY_COL_GNPS_LIB_MATCHES)


"""
Add log10 average peak intensity columns
"""
# Create log10 average peak intensity columns. For original values of 0, set the log10 value to nan, to avoid divide by 0 warning. Note, these are the non-normalized averages
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
# Create a simple copy of the summary table with the columns in FINAL_COLS_ORDER_SIMPLE
summary_table_simple = summary_table[FINAL_COLS_ORDER_SIMPLE].copy()


"""
Export to Excel
"""
# Write results to excel
writer = pd.ExcelWriter(pjoin(OUTPUT_FOLDER, FILENAME_OUTPUT), engine='xlsxwriter')

# write summary_table_simple
write_table_to_excel(writer, summary_table_simple, 'Summary Table Simple')
workbook = writer.book
worksheet = writer.sheets['Summary Table Simple']
format_column(worksheet, summary_table_simple)

# write summary table
write_table_to_excel(writer, summary_table, 'Summary Table')
workbook = writer.book
worksheet = writer.sheets['Summary Table']
format_column(worksheet, summary_table)

# Write filtered tables
# Write a simple filtered table with metabolite significantly present in CC and not MC, sorted by ascending p_val_CC_vs_MC:
# a) p_val_CC_vs_MC < P_VAL_SIG --> metabolites significantly present in CC and not MC
write_table_to_excel(writer, summary_table_simple.loc[(summary_table_simple['p_val_CC_vs_MC'] < 0.05) & ((summary_table_simple['CC_TIC_norm_avg'] > summary_table_simple['MC_TIC_norm_avg']))].sort_values(by='p_val_CC_vs_MC'), 'filter CC vs MC')

# Write a simple filtered table with metabolite significantly present in AR and not MC, sorted by ascending p_val_AR_vs_MC:
# b) p_val_AR_vs_MC < P_VAL_SIG --> metabolites significantly present in AR and not MC
write_table_to_excel(writer, summary_table_simple.loc[(summary_table_simple['p_val_AR_vs_MC'] < 0.05) & ((summary_table_simple['AR_TIC_norm_avg'] > summary_table_simple['MC_TIC_norm_avg']))].sort_values(by='p_val_AR_vs_MC'), 'filter AR vs MC')

# Write a simple filtered table with metabolites significantly more present in CC than AR, sorted by ascending p_val_CC_vs_AR:
# c) p_val_CC_vs_AR < 0.05, CC_cell_norm_avg > AR_cell_norm_avg --> metabolites significantly more present in CC than AR
write_table_to_excel(writer, summary_table_simple.loc[(summary_table_simple['p_val_CC_vs_AR_cell_norm'] < 0.05) & (summary_table_simple['CC_cell_norm_avg'] > summary_table_simple['AR_cell_norm_avg'])].sort_values(by='p_val_CC_vs_AR_cell_norm'), 'filter CC vs AR')

# Write a simple filtered table with metabolites significantly more present in AR than CC, sorted by ascending p_val_CC_vs_AR:
# d) p_val_CC_vs_AR < 0.05, AR_cell_norm_avg > CC_cell_norm_avg --> metabolites significantly more present in AR than CC
write_table_to_excel(writer, summary_table_simple.loc[(summary_table_simple['p_val_CC_vs_AR_cell_norm'] < 0.05) & (summary_table_simple['AR_cell_norm_avg'] > summary_table_simple['CC_cell_norm_avg'])].sort_values(by='p_val_CC_vs_AR_cell_norm'), 'filter AR vs CC')

# Write a simple filtered table for metabolites detected in FAMES sample. Determine this based on that top 100 features detected in FAMES (sort descending FAMES_TIC_norm_avg) (or, could filter for rows exceeding the theshold value of 0.0001 (semi-arbitrary, looked at histogram)).
summary_table_simple.sort_values(by='FAMES_TIC_norm_avg', ascending=False).head(100).to_excel(writer, sheet_name = 'filter FAMES', index = False)

# For each sheet in worksheet, color the MQScore_GNPS and Total_spectrum_similarity_MSDIAL columns. The color gradient will be from white (low) to green (high). 
for sheet in writer.sheets:
    worksheet = writer.sheets[sheet]
    # Add conditional formatting to the MQScore_GNPS column
    color_excel_column(worksheet, summary_table_simple, 'MQScore_GNPS')
    # Add conditional formatting to the Total_spectrum_similarity_MSDIAL column
    color_excel_column(worksheet, summary_table_simple, 'Total_spectrum_similarity_MSDIAL', min_val = 50, max_val = 100)

writer.close()


"""
Generate Histograms to Show Peak Intensity Distributions (use _avg_log10 values)
"""
# For each sample type, generate histograms of the log10 average peak intensities
for sample_type in ['CC', 'AR', 'MC', 'RF', 'FAMES', 'BLANK']:
    # Create a histogram of the log10 average peak intensities
    plt.hist(summary_table_simple[sample_type + '_avg_log10'], bins=20)
    plt.title(sample_type)
    plt.xlabel('log10 average peak intensity')
    plt.ylabel('Frequency')
    plt.savefig(pjoin(OUTPUT_FOLDER, 'histogram_' + sample_type + '_log10_avg_intensity.png'))
    plt.close()


"""
Generate Volcano Plots?
"""
