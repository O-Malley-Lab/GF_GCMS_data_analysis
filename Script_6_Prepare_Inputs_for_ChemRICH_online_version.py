import os
import pandas as pd
"""
***Run Script in Separate Parts ***
"""
"""
Part 1*************************************************************
"""

"""
Functions
"""
def check_key_column(pd_df, key_col_name):
    """
    input:
    pd_df: pandas dataframe

    output:
    boolean: True if proteinID column in pd_df, False if not
    """
    if key_col_name in pd_df.columns:
        return
    else:
        print(f"Error: no {key_col_name} column in pandas dataframe")
        exit()


def make_key_val_dict(pd_df, key_col_name, val_col_name):
    """
    inputs:
    pd_df: pandas dataframe with annotations
    - must have key_col_name column
    key_col_name: specific column name (string) with keys
    val_col_name: specific column name (string) with values to match in its target row based on key match

    output:
    key_val_dict: dictionary with key-value matches
    """
    key_val_dict = {}
    # if no key_col_name column in pd_df, print error message
    check_key_column(pd_df, key_col_name)
    for i in range(len(pd_df)):
        key = pd_df[key_col_name][i]
        val = pd_df[val_col_name][i]
        if key in key_val_dict:
            # if val is different from key_val_dict[key] ...
            # I don't want duplicates of same values
            if val not in key_val_dict[key]:
                # append val to key_val_dict[key]
                key_val_dict[key].append(val)
        else:
            key_val_dict[key] = [val]
    for key_entry in key_val_dict:
        if len(key_val_dict[key_entry]) > 1:
            vals = list(map(str, key_val_dict[key_entry]))
            key_val_dict[key_entry] = ",".join(vals)
        else:
            key_val_dict[key_entry] = key_val_dict[key_entry][0]
    return key_val_dict

def add_to_df(pd_df, key_val_pd, val_cols_list, key_col_name):
    """
    inputs: 
    pd_df: pandas dataframe with proteinIDs
    key_val_pd: pandas dataframe with proteinIDs and annotations
    val_cols_list: list of column names in key_val_pd to add to DGE_summary (list of strings)
    key_col_name: column name shared between dataframes for mapping

    output: pd_df with new columns (pandas dataframe)
    """
    X_dicts_list = []
    for col in val_cols_list:
        X_dicts_list.append(make_key_val_dict(key_val_pd, key_col_name, col))

    # Add annotations to pd_df
    annot_col_num = 0
    for annot_dict in X_dicts_list:
        pd_df[val_cols_list[annot_col_num]] = pd_df[key_col_name].map(annot_dict)
        annot_col_num += 1

    return pd_df

"""
Values
"""
INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'

# Input excel filename in output folder
GF_GNPS_DATA_FILENAME = "GF_GCMS_stats_summary_table.xlsx"
GF_GNPS_DATA_FILE_TAB_NAME = "Summary Table"

# GF_GNPS_DATA_FILENAME excel has 4 columns of interest. The other columns, InchiKeys and Pubchem ID, need to be generated manually with the Pubchem identifier tool: https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi
COLS_OF_INTEREST = ["Compound_Name_GNPS", "SMILES_GNPS", "p_val_CC_vs_AR", "log2_FC_CC_vs_AR", "MQScore_GNPS"]

NEW_COL_NAMES_DICT = {"Compound_Name_GNPS": "Compound Name", "SMILES_GNPS": "SMILES", "p_val_CC_vs_AR": "pvalue", "log2_FC_CC_vs_AR": "foldchange"}

FINAL_COL_ORDER = ["Compound Name", "InChiKeys", "Pubchem ID", "SMILES", "pvalue", "foldchange"]

MQ_SCORE_COL_NAME = "MQScore_GNPS"
MQ_SCORE_CUTOFF = 0.7

SMILES_INCHIKEY_MATCHES_FILENAME = "SMILES_INCHIKEY_matches.txt"
SMILES_PUBCHEMID_MATCHES_FILENAME = "SMILES_PubchemID_matches.txt"

OUTPUT_FILENAME = "GF_GCMS_my_batch_ChemRICH_input.xlsx"

"""
Import data into pandas dataframe
"""
gf_gnps_data = pd.read_excel(os.path.join(OUTPUT_FOLDER, GF_GNPS_DATA_FILENAME))


"""
Apply Confidence Filter
"""
# Filter for rows with MQ score greater than or equal to MQ_SCORE_CUTOFF
gf_gnps_data = gf_gnps_data[gf_gnps_data[MQ_SCORE_COL_NAME] >= MQ_SCORE_CUTOFF]

# For duplicate rows (based on "Compound_Name_GNPS"), keep only rows with this highest MQ score. If there are multiple rows with the same MQ score, keep the first row.
gf_gnps_data = gf_gnps_data.sort_values(by=[MQ_SCORE_COL_NAME], ascending=False)
gf_gnps_data = gf_gnps_data.drop_duplicates(subset=["Compound_Name_GNPS"], keep="first")


"""
Filter for Columns of Interest and No Missing Values
"""
gf_gnps_data = gf_gnps_data[COLS_OF_INTEREST]

# Remove rows with missing values in any of the columns of interest
gf_gnps_data = gf_gnps_data.dropna()


"""
Export SMILES .txt file
"""
# Export SMILES list. Each SMILES entry is separated by a newline character.
with open(os.path.join(OUTPUT_FOLDER, "smiles_ChemRICH.txt"), "w") as f:
    for smiles in gf_gnps_data["SMILES_GNPS"]:
        f.write(smiles + "\n")


"""
***Manual Step: Use "smiles_ChemRICH.txt" to generate InChiKeys and Pubchem ID using the Pubchem identifier tool: https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi***
Save the output as SMILES_INCHIKEY_matches.txt in the output folder.

Repeat to generate Pubchem IDs.
Save the output as SMILES_PubchemID_matches.txt in the output folder.
"""


"""
Part 2*************************************************************
"""
"""
Add InChiKeys and Pubchem IDs to gf_gnps_data
"""
# Use add_to_df to add InChiKeys and Pubchem IDs to gf_gnps_data, using SMILES as the key_col_name column.
# Import SMILES_INCHIKEY_matches.txt and SMILES_PubchemID_matches.txt with column names
smiles_inchikey_matches = pd.read_csv(os.path.join(OUTPUT_FOLDER, SMILES_INCHIKEY_MATCHES_FILENAME), sep="\t", names=["SMILES_GNPS", "InChiKeys"])
smiles_pubchemid_matches = pd.read_csv(os.path.join(OUTPUT_FOLDER, SMILES_PUBCHEMID_MATCHES_FILENAME), sep="\t", names=["SMILES_GNPS", "Pubchem ID"])

# Add InChiKeys and Pubchem IDs to gf_gnps_data
gf_gnps_data = add_to_df(gf_gnps_data, smiles_inchikey_matches, ["InChiKeys"], "SMILES_GNPS")
gf_gnps_data = add_to_df(gf_gnps_data, smiles_pubchemid_matches, ["Pubchem ID"], "SMILES_GNPS")


"""
Remove Rows with Missing Values and Duplicate Values for ID columns
"""
# If any rows are missing any values, remove them
missing_values = gf_gnps_data.isnull().sum(axis=1)
missing_values_indices = missing_values[missing_values > 0].index
if len(missing_values_indices) > 0:
    gf_gnps_data = gf_gnps_data.dropna()

# For rows with an identical value in ID columns (Compound Name, SMILES, InChiKeys, or Pubchem ID) with 1 or more other row, keep only the row with higher MQ score.
gf_gnps_data = gf_gnps_data.sort_values(by=[MQ_SCORE_COL_NAME], ascending=False)
gf_gnps_data = gf_gnps_data.drop_duplicates(subset=["Compound_Name_GNPS"], keep="first")
gf_gnps_data = gf_gnps_data.drop_duplicates(subset=["SMILES_GNPS"], keep="first")
gf_gnps_data = gf_gnps_data.drop_duplicates(subset=["InChiKeys"], keep="first")
gf_gnps_data = gf_gnps_data.drop_duplicates(subset=["Pubchem ID"], keep="first")


"""
Format Datafame for Export and ChemRICH Input
"""
# Rename columns
gf_gnps_data = gf_gnps_data.rename(columns=NEW_COL_NAMES_DICT)

# Reorder columns
gf_gnps_data = gf_gnps_data[FINAL_COL_ORDER]

# Order rows based on Pubchem ID
gf_gnps_data = gf_gnps_data.sort_values(by=["Pubchem ID"])

# For each column, print any duplicated values
for col in gf_gnps_data.columns:
    duplicated_values = gf_gnps_data[gf_gnps_data.duplicated(subset=[col], keep=False)]
    if len(duplicated_values) > 0:
        print(f"Duplicated values in column {col}:")
        print(duplicated_values)

# Note for my data, pyrophosphate and Pyrophosphate are duplicates. Remove the row with Pyrophosphate in the compound name column (lower MQ score)
gf_gnps_data = gf_gnps_data.drop(gf_gnps_data[gf_gnps_data["Compound Name"] == "Pyrophosphate"].index)

# In foldchange column, change -inf to -10. change inf to 10
gf_gnps_data = gf_gnps_data.replace([float("-inf")], -10)
gf_gnps_data = gf_gnps_data.replace([float("inf")], 10)

# Export to excel file with adjusted column widths
with pd.ExcelWriter(os.path.join(OUTPUT_FOLDER, OUTPUT_FILENAME), engine='xlsxwriter') as writer:
    gf_gnps_data.to_excel(writer, index=False)
    worksheet = writer.sheets['Sheet1']
    for idx, col in enumerate(gf_gnps_data.columns):
        worksheet.set_column(idx, idx, len(col) + 2)  # +2 for padding