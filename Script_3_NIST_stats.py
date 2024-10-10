"""
GF GCMS Data Analysis Script 3
Lazarina Butkovich 10/9/24

This script takes in AMDIS NIST processed GCMS data. Additionally, the script analyzes compound matches to (1) the PNNL compound library, with standards run at PNNL on same instrumentation, and (2) NIST search to NIST23 library.

For any given NIST report input for a sample, this script consolidates rows for the same PNNL compound sample name. NIST23 library matches are not included in subsequent analysis but are exported for manual inspection and comparison to potential features of interest.

Then, the script creates and overall dataframe, combining sample abundance data for each compound from each NIST report. Samples are combined based on 'Compound Name' match. A column is added to list the RT and RIs of sample rows combined. The overall dataframe has an 'Area_...' column for each sample (ie: 'Area_AR_1'). The script exports the overall dataframe as an excel file.

"""

import pandas as pd

import numpy as np
import os
from os.path import join as pjoin
import matplotlib.pyplot as plt
# from bioinfokit import visuz
import re
pd.options.mode.copy_on_write = True


"""""""""""""""""""""""""""""""""""""""""""""
Functions
"""""""""""""""""""""""""""""""""""""""""""""
def import_nist_reports_to_dict(report_names_list):
    """
    Import NIST reports as pandas dataframes.

    Input: 
    report_names_list: list of report names (strings) (no .txt)

    Output: 
    dictionary with report names as keys and pandas dataframes as values
    """
    nist_report_dict = {}
    for nist_report_name in report_names_list:
        nist_report_path = pjoin(INPUT_FOLDER, nist_report_name + '.txt')
        nist_report_dict[nist_report_name] = pd.read_csv(nist_report_path, sep='\t', header=0)
    return nist_report_dict

def generate_nist23_dfs_to_dict(report_names_list, nist_report_dict_raw):
    """
    Filter NIST reports to separate NIST23 library matches from PNNL library matches.

    Inputs: 
    report_names_list: list of report names (strings) (no .txt)
    nist_report_dict_raw: dictionary with report names as keys and pandas dataframes as values

    Output: dictionary with report names as keys and pandas dataframes as values
    """
    nist_report_NIST23_dict = {}
    for nist_report_name in report_names_list:
    # Set aside rows with NIST search results
        nist_report_NIST23 = nist_report_dict_raw[nist_report_name].copy()
        nist_report_NIST23 = nist_report_NIST23[nist_report_NIST23['Compound Name'].str.startswith('>')]
        nist_report_NIST23_dict[nist_report_name] = nist_report_NIST23
    return nist_report_NIST23_dict

def consolidate_rows_in_report(nist_report, area_col_name):
    """
    For each unique 'Compound Name', keep only the row with the highest 'Net' value. If that is equal, keep the row with the highest area_col_name value. If those match, keep the first row.

    Inputs:
    nist_report: pandas dataframe of a NIST report
    area_col_name: name of the area column for the report

    Outputs: cleaned up pandas dataframe    
    """
    # Sort the dataframe by 'Net' and area_col_name in descending order
    nist_report = nist_report.sort_values(by=['Net', area_col_name], ascending=[False, False])

    # Drop duplicates based on 'Compound Name' and keep the first row
    nist_report = nist_report.drop_duplicates(subset='Compound Name', keep='first')
    
    return nist_report

def cleanup_PNNL_report(nist_report, nist_report_name, vals_to_remove, cols_to_keep):
    """
    Pass in a dataframe of a PNNL report and clean up based on specific criteria. To add more filters, adjust this function.

    Inputs:
    nist_report: pandas dataframe of a PNNL report
    vals_to_remove: list of strings to remove from 'Compound Name' values
    nist_report_name: name of the report

    Output: cleaned up pandas dataframe
    """
    # vals_to_remove are variations of the '?' in 'Compound Name' that indicate confidence levels. However, we do not need this info in the 'Compound Name' column.
    nist_report = nist_report.sort_values(by=['RT', 'Net'], ascending=[True, False]).drop_duplicates(subset='RT')
    for val in vals_to_remove:
        nist_report['Compound Name'] = nist_report['Compound Name'].str.replace(val, '')
    
    # Use regular expressions to remove patterns '[...] ' or ' [...]'
    nist_report['Compound Name'] = nist_report['Compound Name'].str.replace(r'\[.*?\] ', '', regex=True)
    nist_report['Compound Name'] = nist_report['Compound Name'].str.replace(r' \[.*?\]', '', regex=True)

    # Remove any trailing ' #' from the 'Compound Name' values
    nist_report['Compound Name'] = nist_report['Compound Name'].str.replace(r' \d+$', '', regex=True)

    # keep only the following columns: Name, RT, RI, RI-RI(lib), Net, Weighted, Reverse, (m/z), Area
    nist_report = nist_report[cols_to_keep]

    # For each report, rename the 'Area' file to be in the form 'Area_AR_1', for example, where AR_1 comes from the report name (value before the 2nd '_' in the name).
    short = nist_report_name.split('_')[0] + '_' + nist_report_name.split('_')[1]
    area_col_name = 'Area_' + short
    nist_report.rename(columns={'Area': area_col_name}, inplace=True)

    # For each report, make all letters lower case in the 'Compound Name' column
    nist_report['Compound Name'] = nist_report['Compound Name'].str.lower()

    # Consolidate rows in the report based on 'Compound Name'. Keep only the row with the highest 'Net' value. If that is equal, keep the row with the highest area_col_name value. If those match, keep the first row.
    nist_report = consolidate_rows_in_report(nist_report, area_col_name)

    return nist_report

def generate_PNNL_dfs_to_dict(nist_report_dict_raw, report_names_list, vals_to_remove, cols_to_keep):
    """
    Clean up each NIST report separately for PNNL in-house library matches. This function uses cleanup_PNNL_report.
    Inputs: 
    nist_report_dict_raw: dictionary with report names as keys and pandas dataframes as values
    report_names_list: list of report names (strings) (no .txt)
    vals_to_remove: list of strings to remove from 'Compound Name' values
    cols_to_keep: list of columns to keep in the cleaned up report

    Output: dictionary with report names as keys and pandas dataframes as values   
    """
    nist_report_NIST23_dict = {}
    for nist_report_name in report_names_list:

        # nist_report_cleaned contains only PNNL library matches
        nist_report = nist_report_dict_raw[nist_report_name].copy()
        nist_report = nist_report[~nist_report['Compound Name'].str.startswith('>')]

        # For values in the 'Compound Name' column, replace any of the strings from vals_to_remove with ''
        nist_report = cleanup_PNNL_report(nist_report, nist_report_name, vals_to_remove, cols_to_keep)
        
        # Add the cleaned up PNNL report to the dictionary
        nist_report_NIST23_dict[nist_report_name] = nist_report

    return nist_report_NIST23_dict

def generate_list_all_compound_names(nist_report_dict, report_names_list):
    """
    Generate a list of all unique compound names detected in the NIST reports.

    Inputs:
    nist_report_dict: dictionary with report names as keys and pandas dataframes as values
    report_names_list: list of report names (strings) (no .txt)

    Output: list of all unique compound names detected in the NIST reports
    """
    all_compound_names_detected =[]
    for nist_report_name in report_names_list:
        all_compound_names_detected.extend(nist_report_dict[nist_report_name]['Compound Name'].tolist())

    all_compound_names_detected = list(set(all_compound_names_detected))
    return all_compound_names_detected

def add_data_to_compound_row(new_row, compound_name, nist_report_PNNL_dict, report_names_list, cols_overall_df_list_type, cols_corresponding_overall_df_list_type):
    """
    For a given compound name, add data to the new row in the overall_df.
    
    Inputs:
    new_row: pandas dataframe row
    compound_name: string
    nist_report_PNNL_dict: dictionary with report names as keys and pandas dataframes as values
    report_names_list: list of report names (strings) (no .txt)
    cols_overall_df_list_type: list of column names in the overall_df that are lists
    cols_corresponding_overall_df_list_type: list of column names in the PNNL report that correspond to the overall_df list columns

    Output: None
    """
    # For each report, add the data to the new row
    for nist_report_name in report_names_list:
        # If the compound is in the report, add the data to the new row
        if compound_name in nist_report_PNNL_dict[nist_report_name]['Compound Name'].tolist():
            # Get the row from the report
            row = nist_report_PNNL_dict[nist_report_name][nist_report_PNNL_dict[nist_report_name]['Compound Name'] == compound_name]
            # Add the data to the new row. For overall_df columns 'RT_list', 'RI_list', 'RI-RI(lib)_list', 'Net_list', 'Weighted_list', 'Reverse_list' add the nist_report values from columns 'RT', 'RI', 'RI-RI(lib)', 'Net', 'Weighted', 'Reverse', respectively, to the existing list in the new row.
            for col in cols_overall_df_list_type:
                # Add value. Use COLS_CORRESPONDING_TO_OVERALL_DF_LIST_TYPE to get the corresponding column name in the PNNL report.
                new_row[col] = new_row[col].apply(lambda x: x + row[cols_corresponding_overall_df_list_type[cols_overall_df_list_type.index(col)]].tolist())
            
            # Add the 'Area' value to the new row
            short = nist_report_name.split('_')[0] + '_' + nist_report_name.split('_')[1]
            new_row['Area_' + short] = row['Area_' + short].tolist()

        else:
            # If the compound is not in the report, add a NaN value to the new row for that 'Area_...' column
            short = nist_report_name.split('_')[0] + '_' + nist_report_name.split('_')[1]
            new_row['Area_' + short] = [np.nan]
    return

def generate_compound_row(compound_name, nist_report_PNNL_dict, overall_df, report_names_list, cols_overall_df, cols_overall_df_list_type, cols_corresponding_overall_df_list_type):
    """
    Generate a row for the overall_df for a given compound name.

    Inputs:
    compound_name: string
    nist_report_PNNL_dict: dictionary with report names as keys and pandas dataframes as values
    overall_df: pandas dataframe
    report_names_list: list of report names (strings) (no .txt)
    cols_overall_df: list of column names in the overall_df
    cols_overall_df_list_type: list of column names in the overall_df that are lists
    cols_corresponding_overall_df_list_type: list of column names in the PNNL report that correspond to the overall_df list columns

    Output: None
    """
    # Create a new row for the overall_df
    new_row = pd.DataFrame(columns=cols_overall_df)

    # Add the 'Compound Name' value to the new row
    new_row['Compound Name'] = [compound_name]

    # Initialize lists for each list type column in the new row
    for col in cols_overall_df_list_type:
        new_row[col] = [[]]

    # Add the data to the new row
    add_data_to_compound_row(new_row, compound_name, nist_report_PNNL_dict, report_names_list, cols_overall_df_list_type, cols_corresponding_overall_df_list_type)
    
    # Add the new row to the overall_df
    overall_df = pd.concat([overall_df, new_row], ignore_index=True)
    return overall_df

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

def sci_notation_excel_values(df, writer, workbook, sheet_name):
    """
    Put the values in the 'Area_...' columns in scientific notation.

    Inputs
    df: DataFrame to format
    writer: ExcelWriter object
    workbook: workbook object
    sheet_name: string

    Outputs
    return: None
    """
    for col_num, col_name in enumerate(df.columns):
        if 'Area' in col_name:
            # Set the column format to scientific notation
            writer.sheets[sheet_name].set_column(col_num, col_num, None, workbook.add_format({'num_format': '0.00E+00'}))
    return

def conditional_formatting_area_data_excel(df, writer, sheet_name):
    """
    Apply conditional formatting colors to data values in excel, for columns with 'Area_...".

    Inputs
    df: DataFrame to format
    writer: ExcelWriter object
    sheet_name: string

    Outputs
    return: None
    """
    # Get the max value in the 'Area_...' columns
    max_area = df[[col for col in df.columns if 'Area' in col]].max().max()
    # Get the min value in the 'Area_...' columns
    min_area = df[[col for col in df.columns if 'Area' in col]].min().min()
    # Get the midpoint value
    midpoint = (max_area + min_area) / 2
    # Set the color scale for the conditional formatting from white to green
    color_scale = [{'type': '3_color_scale', 'min_type': 'num', 'min_value': min_area, 'min_color': '#FFFFFF', 'mid_type': 'num', 'mid_value': midpoint, 'mid_color': '#FFEB84', 'max_type': 'num', 'max_value': max_area, 'max_color': '#63BE7B'}]
    # Apply the conditional formatting to the 'Area_...' columns
    for col_num, col_name in enumerate(df.columns):
        if 'Area' in col_name:
            writer.sheets[sheet_name].conditional_format(1, col_num, len(df), col_num, color_scale[0])
    
    return

def conditional_formatting_rt_std_excel(df, writer, sheet_name):
    """
    Apply conditional formatting to the 'RT_std' column in excel. Higher values are redder, and values closer to 0 are white.

    Inputs
    df: DataFrame to format
    writer: ExcelWriter object
    sheet_name: string

    Outputs
    return: None
    """
    # Get the max value in the 'RT_std' column
    max_rt_std = df['RT_std'].max()
    # Get the min value in the 'RT_std' column. Min color is white, not green
    min_rt_std = df['RT_std'].min()
    # Get the midpoint value
    midpoint = (max_rt_std + min_rt_std) / 2
    # Set the color scale for the conditional formatting
    color_scale = [{'type': '3_color_scale', 'min_type': 'num', 'min_value': min_rt_std, 'min_color': '#FFFFFF', 'mid_type': 'num', 'mid_value': midpoint, 'mid_color': '#FFEB84', 'max_type': 'num', 'max_value': max_rt_std, 'max_color': '#F8696B'}]
    # Apply the conditional formatting to the 'RT_std' column
    col_num = df.columns.get_loc('RT_std')
    writer.sheets[sheet_name].conditional_format(1, col_num, len(df), col_num, color_scale[0])
    
    return
"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'

# Input files: from AMDIS + NIST search generated reports for: AR reps 1-4, CC reps 1-4, MC reps 1-4, BLANK reps 1-3, and FAMES
NIST_REPORT_GROUP_NAMES = ['AR_1_NIST', 'AR_2_NIST', 'AR_3_NIST', 'AR_4_NIST', 'CC_1_NIST', 'CC_2_NIST', 'CC_3_NIST', 'CC_4_NIST', 'MC_1_NIST', 'MC_2_NIST', 'MC_3_NIST', 'MC_4_NIST', 'BLANK_1_NIST', 'BLANK_2_NIST', 'BLANK_3_NIST', 'FAMES_1_NIST']

# For comparing RI values between features, use a +/- 5 RI cutoff
RI_CUTOFF = 5

# For cleaning up the PNNL report 'Compound Name' column, remove the following values
VALS_TO_REMOVE_FROM_NAME_COL = ['? ', '?? ', '??? ', '?']

# For each PNNL report, keep only the following columns: 'Compound Name', 'RT', 'RI', 'RI-RI(lib)', 'Net', 'Weighted', 'Reverse', '(m/z)', 'Area'
COLS_TO_KEEP_PNNL_REPORT = ['Compound Name', 'RT', 'RI', 'RI-RI(lib)', 'Net', 'Weighted', 'Reverse', '(m/z)', 'Area']

COLS_TO_KEEP_OVERALL_DF = ['Compound Name', 'RT_list', 'RI_list', 'RI-RI(lib)_list', 'Net_list', 'Weighted_list', 'Reverse_list', '(m/z)_list', 'Area_AR_1', 'Area_AR_2', 'Area_AR_3', 'Area_AR_4', 'Area_CC_1', 'Area_CC_2', 'Area_CC_3', 'Area_CC_4', 'Area_MC_1', 'Area_MC_2', 'Area_MC_3', 'Area_MC_4', 'Area_BLANK_1', 'Area_BLANK_2', 'Area_BLANK_3', 'Area_FAMES_1']

COLS_TO_KEEP_OVERALL_DF_LIST_TYPE = ['RT_list', 'RI_list', 'RI-RI(lib)_list', 'Net_list', 'Weighted_list', 'Reverse_list', '(m/z)_list']

COLS_CORRESPONDING_TO_OVERALL_DF_LIST_TYPE = ['RT', 'RI', 'RI-RI(lib)', 'Net', 'Weighted', 'Reverse', '(m/z)']

FATTY_ACIDS_LIST = [
    'methyl caprylate' , 'methyl caprate', 'methyl laurate', 'methyl myristate', 'methyl palmitate', 'methyl stearate', 'methyl eicosanoate', 'methyl docosanoate', 'methyl linocerate', 'methyl hexacosanoate', 'methyl octacosanoate', 'capric acid', 'lauric acid', 'myristic acid', 'palmitic acid', 'palmitoleic acid', 'stearic acid', 'oleic acid', 'linoleic acid', 'arachidic acid', 'docosahexaenoic acid'
]


"""""""""""""""""""""""""""""""""""""""""""""
Main
"""""""""""""""""""""""""""""""""""""""""""""
"""
Import NIST reports as pandas dataframes
"""
# For each NIST report, add the tab-delimited text file data as a pandas dataframe to a dictionary, with the key as the report name
nist_report_dict_raw = import_nist_reports_to_dict(NIST_REPORT_GROUP_NAMES)

"""
Rename 'Name' columns to be 'Compound Name' for each report
"""
# For each report, rename the 'Name' column to be 'Compound Name'
for nist_report_name in NIST_REPORT_GROUP_NAMES:
    nist_report_dict_raw[nist_report_name].rename(columns={'Name': 'Compound Name'}, inplace=True)

"""
Clean up each NIST report separately for PNNL in-house library matches
"""
# For each report, create a new dataframe with the cleaned up data.
# Create a separate storage dictionaries for NIST23 library matches and PNNL library matches
nist_report_NIST23_dict = generate_nist23_dfs_to_dict(NIST_REPORT_GROUP_NAMES, nist_report_dict_raw)
nist_report_PNNL_dict = generate_PNNL_dfs_to_dict(nist_report_dict_raw, NIST_REPORT_GROUP_NAMES, VALS_TO_REMOVE_FROM_NAME_COL, COLS_TO_KEEP_PNNL_REPORT)


"""
Export NIST reports per sample, for matches to NIST23 library and PNNL in-house library
"""
# For each report, export as an excel file
for nist_report_name in NIST_REPORT_GROUP_NAMES:
    nist_report_NIST23_dict[nist_report_name].to_excel(pjoin(TEMP_FOLDER, nist_report_name + '23_matches.xlsx'), index=False)
    nist_report_PNNL_dict[nist_report_name].to_excel(pjoin(TEMP_FOLDER, nist_report_name + '_PNNL_matches.xlsx'), index=False)


"""
Consolidate Reports to Overall Dataframe
"""
# Consolidate the cleaned up PNNL reports into one overall dataframe (overall_df), matching features based on 'Compound Name'. 

# Create overall_df
overall_df = pd.DataFrame(columns=COLS_TO_KEEP_OVERALL_DF)

# Generate a list of all unique 'Compound Name' values from all reports.
all_PNNL_compound_names_detected = generate_list_all_compound_names(nist_report_PNNL_dict, NIST_REPORT_GROUP_NAMES)

# For each unique 'Compound Name', create a row in the overall_df
for compound_name in all_PNNL_compound_names_detected:
    # Generate a row for the compound in the overall_df
    overall_df = generate_compound_row(compound_name, nist_report_PNNL_dict, overall_df, NIST_REPORT_GROUP_NAMES, COLS_TO_KEEP_OVERALL_DF, COLS_TO_KEEP_OVERALL_DF_LIST_TYPE, COLS_CORRESPONDING_TO_OVERALL_DF_LIST_TYPE)

# Add RT_avg and RT_std columns after RT_list column
overall_df.insert(2, 'RT_avg', overall_df['RT_list'].apply(lambda x: np.mean(x)))
overall_df.insert(3, 'RT_std', overall_df['RT_list'].apply(lambda x: np.std(x)))


"""
Perform Fatty Acid Profiling for Each Sample
"""
fatty_acids_df = overall_df.copy()

# Filter fatty_acids_df for 'Compound Name' values in FATTY_ACIDS_LIST
fatty_acids_df = fatty_acids_df[fatty_acids_df['Compound Name'].isin(FATTY_ACIDS_LIST)]


"""
Export Dataframes in Excel File
"""
writer = pd.ExcelWriter(pjoin(OUTPUT_FOLDER, 'NIST_PNNL_lib_matches_summary.xlsx'), engine='xlsxwriter')

# Write each dataframe to a different sheet
overall_df.to_excel(writer, sheet_name='Summary', index=False)
fatty_acids_df.to_excel(writer, sheet_name='Fatty Acids', index=False)

 # Format the excel sheets so that the column width matches the size of the header text
workbook = writer.book
# For each table and corresponding excel tab, format width
format_column(writer.sheets['Summary'], overall_df)
format_column(writer.sheets['Fatty Acids'], fatty_acids_df)

# In each sheet, format the 'Area_...' columns to be in scientific notation
sci_notation_excel_values(overall_df, writer, workbook, 'Summary')
sci_notation_excel_values(fatty_acids_df, writer, workbook, 'Fatty Acids')

# Apply conditional formatting to 'Area_...' columns
conditional_formatting_area_data_excel(overall_df, writer, 'Summary')
conditional_formatting_area_data_excel(fatty_acids_df, writer, 'Fatty Acids')

# Apply conditional formatting to 'RT_std' column
conditional_formatting_rt_std_excel(overall_df, writer, 'Summary')
conditional_formatting_rt_std_excel(fatty_acids_df, writer, 'Fatty Acids')
    
writer.close()