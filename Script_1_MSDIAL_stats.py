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
from scipy.stats import ttest_ind

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

def msdial_table_cleanup(df, key_col, cols_name_converter, original_key_col = 'Alignment ID', metab_col = 'Metabolite name', unknown_val = 'Unknown'):
    """
    Function specific to cleaning up data tables from MS-DIAL

    Inputs
    df: DataFrame to clean up
    key_col: Name of the shared name column
    cols_name_converter: Dictionary to rename the columns
    original_key_col: Name of the original key column
    metab_col: Name of the metabolite name column
    unknown_val: Value to replace in the metabolite name column

    Outputs
    return: DataFrame
    """
    # Create the shared name column and data values, where the values are 1 plus the 'Alignment ID' values
    df[key_col] = df[original_key_col] + 1

    # Move KEY_COL to the first column
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]

    # For values in the 'Metabolite name MSDIAL column', if the value is "Unknown", change to a blank value
    df.loc[df[metab_col] == unknown_val, metab_col] = ''

    # Use COLS_NAME_CONVERTER to rename the columns
    df.rename(columns=cols_name_converter, inplace=True)
    return df

def generate_pval_col(df_data, sample_groups_dict, grp1_name, grp2_name, suffix=''):
    """
    Generate a p-value column for the comparison of two sample groups. The p-value is generated using a t-test.

    Inputs
    df_data: DataFrame to add the p-value column to
    sample_groups_dict: Dictionary of sample groups
    grp1_name: Name of the first sample group
    grp2_name: Name of the second sample group

    Outputs
    return    
    """
    p_val_list = []
    for index, row in df_data.iterrows():
        # if either group has NaN values, append NaN to the p_val_list
        # if row[sample_groups_dict[grp1_name]].isnull().values.any() or row[sample_groups_dict[grp2_name]].isnull().values.any():
        #     p_val_list.append(np.nan)
        #     continue
        val = ttest_ind(row[sample_groups_dict[grp1_name]], row[sample_groups_dict[grp2_name]])[1]
        if np.isnan(val):
            p_val_list.append(np.nan)
            continue
        p_val_list.append(val)
    df_data['p_val_' + grp1_name + '_vs_' + grp2_name + suffix] = p_val_list

    return

def generate_log2_fc_col(df_data, grp1_name, grp2_name, data_col_prefix='_TIC_norm_avg', suffix=''):
    # if there is a divide by zero error, set the log2 fold change to NaN
    log2_fc_list = []
    for index, row in df_data.iterrows():
        # check if divide by zero
        if row[grp2_name + data_col_prefix] == 0:
            log2_fc_list.append(np.inf)
            continue
        # if log2(0), set to -inf
        if row[grp1_name + data_col_prefix] == 0:
            log2_fc_list.append(-np.inf)
            continue
        # check if nan values
        if np.isnan(row[grp1_name + data_col_prefix]) or np.isnan(row[grp2_name + data_col_prefix]):
            log2_fc_list.append(np.nan)
            continue      
        val = np.log2(row[grp1_name + data_col_prefix] / row[grp2_name + data_col_prefix])
        log2_fc_list.append(val)
    df_data['log2_FC_' + grp1_name + '_vs_' + grp2_name + suffix] = log2_fc_list
    return


"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'

KEY_COL = 'shared name'

# MSDIAL output file with TIC normalized data
FILENAME_MSDIAL_OUTPUT_NORM_TIC = 'MSDIAL_norm_TIC_output.xlsx'

# MSDIAL output file with peak area data
FILENAME_MSDIAL_OUTPUT_AREA = 'MSDIAL_area_output.xlsx'

# Cell pellet weight data for direct comparison of CC and AR relative peak intensities
FILENAME_CELL_PELLET_WEIGHTS = 'GF_cell_pellet_weights.xlsx'
# Column names: 'Sample', 'Sample Mass mg'

# Dictionary of tuples to describe pre- and post- strings in sample names (the middle part is 1 to n, where n=3 or 4 for biological (for AR and CC) or technical (for BLANK, MC, RF) replicates)
SAMPLE_NAME_PRE_POST_STRS_DICT = {'CC':('OMALL_RFS_CC','_M'),'AR':('OMALL_RFS_AR_S4_','_M'),'MC':('OMALL_RFS_MC','_M'),'RF':('OMALL_RFS_RF','_M'), 'FAMES':('GCMS_FAMES_0','_GCMS01_20201209'), 'BLANK':('GCMS_BLANK_0','_GCMS01_20201209')}
REPLICATE_NUMS = {'CC':4, 'AR':4, 'MC':4, 'RF':4,'FAMES':1,'BLANK':3}

OUTPUT_FILENAME = 'MSDIAL_stats.xlsx'

COLS_NAME_CONVERTER = {'Alignment ID': 'Alignment_ID_MSDIAL','Average Rt(min)':'RT', 'Precursor_MZ':'EI_spectra_quant_mass', 'Quant mass': 'Quant_mass', 'Compound_Name':'Compound_Name_GNPS','MQScore':'MQScore_GNPS', 'Smiles':'SMILES_GNPS', 'INCHI':'INCHI_GNPS', 'Metabolite name': 'Metabolite_name_MSDIAL', 'SMILES':'SMILES_MSDIAL', 'INCHI':'INCHI_GNPS', 'molecular_formula':'molecular_formula_GNPS', 'npclassifier_superclass':'npclassifier_superclass_GNPS', 'npclassifier_class':'npclassifier_class_GNPS', 'npclassifier_pathway':'npclassifier_pathway_GNPS','Compound_Source':'Compound_Source_GNPS', 'Data_Collector':'Data_Collector_GNPS', 'Instrument':'Instrument_GNPS', 'Total spectrum similarity': 'Total_spectrum_similarity_MSDIAL'}

COLS_TO_KEEP_SUMMARY_OUTPUT = ['shared name', 'Alignment_ID_MSDIAL', 'RT', 'Quant_mass', 'Metabolite_name_MSDIAL', 'Total_spectrum_similarity_MSDIAL',  'SMILES_MSDIAL', 
'p_val_CC_vs_AR_cell_norm', 'log2_FC_CC_vs_AR_cell_norm',
'p_val_CC_vs_AR', 'log2_FC_CC_vs_AR',
'p_val_CC_vs_MC', 'log2_FC_CC_vs_MC',
'p_val_AR_vs_MC', 'log2_FC_AR_vs_MC',
'p_val_CC_vs_BLANK', 'log2_FC_CC_vs_BLANK',
'p_val_AR_vs_BLANK', 'log2_FC_AR_vs_BLANK',
'p_val_FAMES_vs_BLANK', 'log2_FC_FAMES_vs_BLANK',
'CC_TIC_norm_avg', 'CC_TIC_norm_std', 'AR_TIC_norm_avg', 'AR_TIC_norm_std',
'MC_TIC_norm_avg', 'MC_TIC_norm_std', 'RF_TIC_norm_avg', 'RF_TIC_norm_std',
'CC_cell_norm_avg', 'CC_TIC_norm_avg', 'CC_TIC_norm_std',
'AR_cell_norm_avg', 'AR_TIC_norm_avg',  'AR_TIC_norm_std',
'MC_TIC_norm_avg', 'MC_TIC_norm_std', 
'RF_TIC_norm_avg', 'RF_TIC_norm_std', 'FAMES_TIC_norm_avg',
'FAMES_TIC_norm_std', 'BLANK_TIC_norm_avg', 'BLANK_TIC_norm_std'] 

P_VAL_SIG = 0.05


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

# Ensure data tables are ordered by 'Alignment ID' column
df_msdial_norm_tic.sort_values('Alignment ID', inplace=True)
df_msdial_area.sort_values('Alignment ID', inplace=True)


"""
Clean up MS-DIAL data tables
"""
# Clean up steps specific to the MS-DIAL data tables to make easier for downstream use (ie: add shared name column)
df_msdial_area = msdial_table_cleanup(df_msdial_area, KEY_COL, COLS_NAME_CONVERTER)
df_msdial_norm_tic = msdial_table_cleanup(df_msdial_norm_tic, KEY_COL, COLS_NAME_CONVERTER)


"""
Generate Sample Name Columns Based on Templates
"""
# See SAMPLE_NAME_PRE_POST_STRS_DICT and REPLICATE_NUMS for how to assemble sample names. 
# Example: OMALL_RFS_CC1_M, OMALL_RFS_CC2_M, OMALL_RFS_CC3_M, OMALL_RFS_CC4_M

sample_cols_all = []
sample_groups_dict = {}
for key in SAMPLE_NAME_PRE_POST_STRS_DICT:
    pre_str = SAMPLE_NAME_PRE_POST_STRS_DICT[key][0]
    post_str = SAMPLE_NAME_PRE_POST_STRS_DICT[key][1]
    replicate_total = REPLICATE_NUMS[key]
    for i in range(1, replicate_total+1):
        col_name = pre_str + str(i) + post_str
        sample_cols_all.append(col_name)
        # Add dictinoary entry for key as key and list of sample names as values
        if key in sample_groups_dict:
            sample_groups_dict[key].append(col_name)
        else:
            sample_groups_dict[key] = [col_name]


"""
Normalize Peak Area Data Based on Cell Pellet Data
"""
# Keep only the CC and AR columns from the peak area data, as well as the key column
# In sample_groups_dict, keep only column names that are values of the keys 'CC' and 'AR'
cols_to_keep_cell_norm = []
sample_cols_exp = []
# Add the key column to the list of columns to keep
cols_to_keep_cell_norm.append(KEY_COL)

for key in sample_groups_dict:
    if key == 'CC' or key == 'AR':
        cols_to_keep_cell_norm += sample_groups_dict[key]
        sample_cols_exp += sample_groups_dict[key]

# Keep only the columns to keep
df_msdial_area_cell_norm = df_msdial_area.copy()
df_msdial_area_cell_norm = df_msdial_area_cell_norm[cols_to_keep_cell_norm]

# Match the sample columns to the rows in the cell pellet data. Divide the data values in the sample columns by the 'Sample Mass mg' value in the cell pellet data.
for col in sample_cols_exp:
    mass = df_cell_pellet_weights.loc[df_cell_pellet_weights['Sample'] == col, 'Sample Mass mg'].values[0]
    df_msdial_area_cell_norm[col] = df_msdial_area_cell_norm[col]/mass


"""
Generate Statistics for CC and AR Peak Area Data Normalized by Cell Pellet Weight
"""
# Create CC_cell_norm_avg and AR_cell_norm_avg columns
df_msdial_area_cell_norm['CC_cell_norm_avg'] = df_msdial_area_cell_norm[sample_groups_dict['CC']].mean(axis=1)
df_msdial_area_cell_norm['AR_cell_norm_avg'] = df_msdial_area_cell_norm[sample_groups_dict['AR']].mean(axis=1)

# Create CC_cell_norm_std and AR_cell_norm_std columns
df_msdial_area_cell_norm['CC_cell_norm_std'] = df_msdial_area_cell_norm[sample_groups_dict['CC']].std(axis=1)
df_msdial_area_cell_norm['AR_cell_norm_std'] = df_msdial_area_cell_norm[sample_groups_dict['AR']].std(axis=1)

# # Generate t test p-values for CC vs AR; note that the p-values are generated using the cell normalized data
# p_val_CC_vs_AR = []
# for index, row in df_msdial_area_cell_norm.iterrows():
#     p_val = ttest_ind(row[sample_groups_dict['CC']], row[sample_groups_dict['AR']])[1]
#     p_val_CC_vs_AR.append(p_val)
# df_msdial_area_cell_norm['p_val_CC_vs_AR_cell_norm'] = p_val_CC_vs_AR

generate_pval_col(df_msdial_area_cell_norm, sample_groups_dict, 'CC', 'AR', suffix = '_cell_norm')
generate_log2_fc_col(df_msdial_area_cell_norm, 'CC', 'AR', data_col_prefix = '_cell_norm_avg', suffix = '_cell_norm')


"""
Generate Statistics for TIC Normalized Data (CC vs MC and AR vs MC)
"""
# Make copy df of df_msdial_norm_tic with only the key column and sample columns
df_msdial_norm_tic_stats = df_msdial_norm_tic.copy()
# Keep only the key column and sample columns
df_msdial_norm_tic_stats = df_msdial_norm_tic_stats[[KEY_COL] + sample_cols_all]

# Add average and standard deviation columns for all sample groups (label like CC_TIC_norm_avg and CC_TIC_norm_std)
for key in sample_groups_dict:
    df_msdial_norm_tic_stats[key + '_TIC_norm_avg'] = df_msdial_norm_tic_stats[sample_groups_dict[key]].mean(axis=1)
    df_msdial_norm_tic_stats[key + '_TIC_norm_std'] = df_msdial_norm_tic_stats[sample_groups_dict[key]].std(axis=1)

# Generate p-values and log2 fold-change values for CC vs AR (TIC), CC vs MC, AR vs MC, CC vs AR, CC vs BLANK, AR vs BLANK, FAMES vs BLANK
generate_pval_col(df_msdial_norm_tic_stats, sample_groups_dict, 'CC', 'AR')
generate_log2_fc_col(df_msdial_norm_tic_stats, 'CC', 'AR')

generate_pval_col(df_msdial_norm_tic_stats, sample_groups_dict, 'CC', 'MC')
generate_log2_fc_col(df_msdial_norm_tic_stats, 'CC', 'MC')

generate_pval_col(df_msdial_norm_tic_stats, sample_groups_dict, 'AR', 'MC')
generate_log2_fc_col(df_msdial_norm_tic_stats, 'AR', 'MC')

generate_pval_col(df_msdial_norm_tic_stats, sample_groups_dict, 'CC', 'BLANK')
generate_log2_fc_col(df_msdial_norm_tic_stats, 'CC', 'BLANK')

generate_pval_col(df_msdial_norm_tic_stats, sample_groups_dict, 'AR', 'BLANK')
generate_log2_fc_col(df_msdial_norm_tic_stats, 'AR', 'BLANK')

generate_pval_col(df_msdial_norm_tic_stats, sample_groups_dict, 'FAMES', 'BLANK')
generate_log2_fc_col(df_msdial_norm_tic_stats, 'FAMES', 'BLANK')


"""
Assemble Summary Excel with Relevant Statistics
"""
df_msdial_summary = df_msdial_area.copy()

# Use combine_dfs to add columns from df_msdial_norm_tic
cols_to_add_tic = ['shared name', 
    'p_val_CC_vs_AR', 'log2_FC_CC_vs_AR',
    'p_val_CC_vs_MC', 'log2_FC_CC_vs_MC',
    'p_val_AR_vs_MC', 'log2_FC_AR_vs_MC',
    'p_val_CC_vs_BLANK', 'log2_FC_CC_vs_BLANK',
    'p_val_AR_vs_BLANK', 'log2_FC_AR_vs_BLANK',
    'p_val_FAMES_vs_BLANK', 'log2_FC_FAMES_vs_BLANK',
    'CC_TIC_norm_avg', 'CC_TIC_norm_std',
    'AR_TIC_norm_avg', 'AR_TIC_norm_std',
    'MC_TIC_norm_avg', 'MC_TIC_norm_std',
    'RF_TIC_norm_avg', 'RF_TIC_norm_std',
    'FAMES_TIC_norm_avg', 'FAMES_TIC_norm_std',
    'BLANK_TIC_norm_avg', 'BLANK_TIC_norm_std']

combine_dfs(df_msdial_summary, df_msdial_norm_tic_stats, cols_to_add_tic, KEY_COL, KEY_COL)

# Add columns from df_msdial_area_cell_norm
cols_to_add_cell_norm = ['shared name',
    'p_val_CC_vs_AR_cell_norm', 'log2_FC_CC_vs_AR_cell_norm',
    'CC_cell_norm_avg', 'CC_cell_norm_std',
    'AR_cell_norm_avg', 'AR_cell_norm_std']

combine_dfs(df_msdial_summary, df_msdial_area_cell_norm, cols_to_add_cell_norm, KEY_COL, KEY_COL)


"""
Export excel files
"""
# Have all excels as different sheets of one excel output
writer = pd.ExcelWriter(pjoin(TEMP_FOLDER, OUTPUT_FILENAME), engine='xlsxwriter')

# Write the summary table to the excel file
write_table_to_excel(writer, df_msdial_summary, 'Summary')

# # Write simpler summary table with most relevant columns
df_msdial_summary_output = df_msdial_summary[COLS_TO_KEEP_SUMMARY_OUTPUT]
# to-do: unsure why, but write_table_to_excel function does not work for df_msdial_summary_output
# write_table_to_excel(writer, df_msdial_summary_output, 'Summary Simple')
df_msdial_summary_output.to_excel(writer, sheet_name = 'Summary Simple', index = False)

# Write the cell pellet normalized stats
write_table_to_excel(writer, df_msdial_area_cell_norm, 'Cell Norm Stats')

# Write the TIC normalized stats
write_table_to_excel(writer, df_msdial_norm_tic_stats, 'TIC Norm Stats')

writer.close()