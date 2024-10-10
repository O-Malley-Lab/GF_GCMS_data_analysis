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

"""
Notes for Fatty Acids to Profile
"""
# FAMES related
# 'methyl caprylate' - C8 in FAMES
# 'methyl caprate' - C10 in FAMES
# 'methyl laurate' - C12 in FAMES
# 'methyl myristate' - C14 in FAMES
# 'methyl palmitate' - C16 in FAMES
# 'methyl stearate' - C18 in FAMES
# 'methyl eicosanoate' - C20 in FAMES
# 'methyl docosanoate' - C22 in FAMES
# 'methyl linocerate' - C24 in FAMES
# 'methyl hexacosanoate' - C26 in FAMES
# 'methyl octacosanoate' - C28 in FAMES

# PNNL library, FAs from Kar et al. GCMS profiling for anaerobic gut fungi
# 'capric acid' - C10:0, Decanoic
# 'lauric acid' -  C12:0, Lauric
# cannot find in PNNL lib - C13:0, Tridecanoic
# 'myristic acid' - C14:0, Myristic
# cannot find in PNNL lib - C15:0, Pentadecanoic
# 'palmitic acid' - C16:0, Palmitic
# 'palmitoleic acid' - C16:1, Palmitoleic
# 'stearic acid' - C18:0, Stearic
# 'oleic acid' - C18:1n9c, Oleic
# 'linoleic acid' - C18:2n6c, Linoleic
# 'arachidic acid' - C20:0, Arachidic
# cannot find in PNNL lib - C20:1, Eicosenoic
# cannot find in PNNL lib - C20:4n6, Arachidonic
# cannot find in PNNL lib - C21:0, Heneicosanoic
# cannot find in PNNL lib - C22:1n9, Erucic
# 'docosahexaenoic acid' - C22:6n3, Cis-4,7,10,13,16,19- Docosahexaenoic
# cannot find in PNNL lib - C24:1, Nervonic

# Optional to add: NIST library, FAs from Kar et al. GCMS profiling for anaerobic gut fungi
# '>Pentadecanoic acid, TMS derivative' - C15:0, Pentadecanoic
# '>Heneicosanoic acid, TMS derivative' - C21:0, Heneicosanoic
# ...

# Fatty Acids of Interest, as they appear in PNNL in-house library
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
# 5) Make all letters in the 'Name' column lower case

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
    nist_report_cleaned_dict[nist_report] = nist_report_cleaned_dict[nist_report][['Name', 'RT', 'RI', 'RI-RI(lib)', 'Net', 'Weighted', 'Reverse', '(m/z)', 'Area']]

# For each report, rename the 'Area' file to be in the form 'Area_AR_1', for example, where AR_1 comes from the report name (value before the 2nd '_' in the name).
for nist_report in NIST_REPORT_GROUP_NAMES:
    report_name = nist_report.split('_')[0] + '_' + nist_report.split('_')[1]
    nist_report_cleaned_dict[nist_report].rename(columns={'Area': 'Area_' + report_name}, inplace=True)

# For each report, make all letters lower case in the 'Name' column
for nist_report in NIST_REPORT_GROUP_NAMES:
    nist_report_cleaned_dict[nist_report]['Name'] = nist_report_cleaned_dict[nist_report]['Name'].str.lower()

"""
Export cleaned up NIST reports
"""
# For each cleaned up report, export as an excel file
for nist_report in NIST_REPORT_GROUP_NAMES:
    nist_report_cleaned_path = pjoin(TEMP_FOLDER, nist_report + '_cleaned.xlsx')
    nist_report_cleaned_dict[nist_report].to_excel(nist_report_cleaned_path, index=False)

"""
Perform Fatty Acid Profiling for Each Sample
"""



