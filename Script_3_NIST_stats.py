"""
GF GCMS Data Analysis Script 3
Lazarina Butkovich 10/9/24

This script takes in AMDIS NIST processed GCMS data. Additionally, the script analyzes compound matches to (1) the PNNL compound library, with standards run at PNNL on same instrumentation, and (2) NIST search to NIST23 library.

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

"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'

# Input files: from AMDIS + NIST search generated reports for: AR reps 1-4, CC reps 1-4, MC reps 1-4, and FAMES
NIST_REPORT_GROUP_NAMES = ['AR_1_NIST', 'AR_2_NIST', 'AR_3_NIST', 'AR_4_NIST', 'CC_1_NIST', 'CC_2_NIST', 'CC_3_NIST', 'CC_4_NIST', 'MC_1_NIST', 'MC_2_NIST', 'MC_3_NIST', 'MC_4_NIST', 'FAMES_1_NIST']

# For comparing RI values between features, use a +/- 5 RI cutoff
RI_CUTOFF = 5

"""""""""""""""""""""""""""""""""""""""""""""
Main
"""""""""""""""""""""""""""""""""""""""""""""
"""
Import NIST reports as pandas dataframes
"""
# For each NIST report, add the tab-delimited text file data as a pandas dataframe to a dictionary, with the key as the report name
nist_report_dict = {}
for nist_report in NIST_REPORT_GROUP_NAMES:
    nist_report_path = pjoin(INPUT_FOLDER, nist_report + '.txt')
    nist_report_dict[nist_report] = pd.read_csv(nist_report_path, sep='\t', header=0)


"""
Clean up each NIST report separately for PNNL in-house library matches
"""
# Header:
# 'FileName', 'CAS', 'Name', 'RT', 'RI', 'Width', 'Purity', 'Model', 'Min. Abund.', 'Amount', 'Scan', 'Peak Tailing', 'S/N (total)', 'Base Peak', 'Max. Amount', 'Area', 'Intgr.Signal', 'Max. Area', 'Extra Width', 'Models', 'Frac. Good', 'Expec. RT', 'RI-RI(lib)', 'Net', 'Weighted', 'Simple', 'Reverse', 'Corrections', '(m/z)', 'S/N (m/z)', 'Area % (m/z)', 'Conc.', 'RT-RT(lib)'

# For each report, we want to consolidate rows within a report, such that:
# 1) For PNNL library matches, we will remove all rows from the NIST search (rows with a value in the 'Name' column that starts with '>')
# 2) For all rows with any given 'RT' value, we will only keep the row with the highest Net value 
# 3) Next, for each row, we will look at the value in the 'Name' column and clean it up if needed. If any of the following are present, remove: '? ', '?? ', '??? ', '[#] ' and ' [#]' where # is any set of characters, usually numbers
# 4) Some 'Name' values end with a number, in the form of ' #'. This number is not necessary so we will remove the ' #' if it appears as the end of the 'Name' value

# For each report, we will create a new dataframe with the cleaned up data
nist_report_cleaned_dict = {}
vals_to_remove = ['? ', '?? ', '??? ', '?']
for nist_report in NIST_REPORT_GROUP_NAMES:
    nist_report_cleaned = nist_report_dict[nist_report].copy()
    nist_report_cleaned = nist_report_cleaned[~nist_report_cleaned['Name'].str.startswith('>')]
    nist_report_cleaned = nist_report_cleaned.sort_values(by=['RT', 'Net'], ascending=[True, False]).drop_duplicates(subset='RT')
    
    # For values in the 'Name' column, replace any of the strings from vals_to_remove with ''
    for val in vals_to_remove:
        nist_report_cleaned['Name'] = nist_report_cleaned['Name'].str.replace(val, '')
    
    # Use regular expressions to remove patterns '[...] ' or ' [...]'
    nist_report_cleaned['Name'] = nist_report_cleaned['Name'].str.replace(r'\[.*?\] ', '', regex=True)
    nist_report_cleaned['Name'] = nist_report_cleaned['Name'].str.replace(r' \[.*?\]', '', regex=True)

    # Remove any trailing ' #' from the 'Name' values
    nist_report_cleaned['Name'] = nist_report_cleaned['Name'].str.replace(r' \d+$', '', regex=True)
    
    nist_report_cleaned_dict[nist_report] = nist_report_cleaned

# For each report in nist_report_cleaned_dict, keep only the following columns: Name, RT, RI, RI-RI(lib), Net, Weighted, Reverse, (m/z)
for nist_report in NIST_REPORT_GROUP_NAMES:
    nist_report_cleaned_dict[nist_report] = nist_report_cleaned_dict[nist_report][['Name', 'RT', 'RI', 'RI-RI(lib)', 'Net', 'Weighted', 'Reverse', '(m/z)', 'Amount']]

# The values in the 'Amount' column are strings with a '%' at the end. We will remove the '%' and convert the values to floats. The value should still represent the actual percent value (ie: divide by 100)
for nist_report in NIST_REPORT_GROUP_NAMES:
    nist_report_cleaned_dict[nist_report].loc[:, 'Amount'] = nist_report_cleaned_dict[nist_report]['Amount'].str.replace('%', '').astype(float) / 100

# For each report, rename the 'Amount' file to be in the form 'Amount_AR_1', for example, where AR_1 comes from the report name (value before the 2nd '_' in the name).
for nist_report in NIST_REPORT_GROUP_NAMES:
    report_name = nist_report.split('_')[0] + '_' + nist_report.split('_')[1]
    nist_report_cleaned_dict[nist_report].rename(columns={'Amount': 'Amount_' + report_name}, inplace=True)


# """
# Export cleaned up NIST reports
# """
# # For each cleaned up report, export as an excel file
# for nist_report in NIST_REPORT_GROUP_NAMES:
#     nist_report_cleaned_path = pjoin(TEMP_FOLDER, nist_report + '_cleaned.xlsx')
#     nist_report_cleaned_dict[nist_report].to_excel(nist_report_cleaned_path, index=False)

"""
Perform Fatty Acid Profiling for Each Sample
"""




# """
# Consolidate Reports to Overall Dataframe
# """
# Next, we want to consolidate the cleaned up NIST reports into one dataframe. We will match features based on values in 'RI' (+/- RI_CUTOFF) and 'Name'. The overall dataframe (overall_df) will start as a copy of the first cleaned up report. For each subsequent report, we will add a corresponding "Amount_..." column and iterate through the rows in overall_df to find any matches. If there is a match, we will add the value from the 'Amount_...' column to the matching row. If there is not a match, we will add a new row to overall_df with the values from the report and 0 value for previous 'Amount_...' columns. We will iterate over all reports in NIST_REPORT_GROUP_NAMES to create the overall_df.

# # Start with the first cleaned up report
# overall_df = nist_report_cleaned_dict[NIST_REPORT_GROUP_NAMES[0]].copy()

# # For each subsequent report, add the 'Amount_...' column to the overall_df
# cols_to_add = ['Amount_' + nist_report.split('_')[0] + '_' + nist_report.split('_')[1] for nist_report in NIST_REPORT_GROUP_NAMES]

# i=0
# for col in cols_to_add:
#     if i == 0:
#         # skip for the first report
#         i+=1
#         continue
#     # Add the 'Amount_...' column to the overall_df. the values will be nan for now
#     overall_df[col] = np.nan


# # For each subsequent report, iterate through the rows in overall_df to find matches

# # Export the overall_df as an excel file
# overall_df_path = pjoin(TEMP_FOLDER, 'overall_PNNL_lib_matches.xlsx')
# overall_df.to_excel(overall_df_path, index=False)
