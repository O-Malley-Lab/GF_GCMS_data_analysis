from os.path import join as pjoin
import pandas as pd
import pubchempy as pcp

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

# Input MeSh Class Prediction Filename (created partway through code with Script 7)
MESH_CMPD_CLASSES_FILENAME = "MeSh_Prediction_Results_edited.xlsx"
MESH_PREDICTIONS_COL_NAME = "MeSH_Class"

# GF_GNPS_DATA_FILENAME excel has 4 columns of interest. The other columns, InchiKeys and Pubchem ID, need to be generated manually with the Pubchem identifier tool: https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi
COLS_OF_INTEREST = ["Compound_Name_GNPS", "SMILES_GNPS", "p_val_CC_vs_AR", "log2_FC_CC_vs_AR", "MQScore_GNPS"]

COL_NAME_CONVERT_DICT = {"Compound_Name_GNPS": "compound_name", "SMILES_GNPS": "smiles", "p_val_CC_vs_AR": "pvalue", "log2_FC_CC_vs_AR": "effect_size", "MeSH_Class": "set"}

MESH_COL_ORDER = ["compound_name", "pubchem_id", "smiles"]

FINAL_COL_ORDER = ["compound_name", "smiles", "pvalue", "effect_size", "set"]

MQ_SCORE_COL_NAME = "MQScore_GNPS"
MQ_SCORE_CUTOFF = 0.7

SMILES_PUBCHEMID_MATCHES_FILENAME = "SMILES_PubchemID_matches.txt"

OUTPUT_FILENAME_MESH_INPUT = "GF_GCMS_MESH_prediction_input.xlsx"
OUTPUT_FILENAME_CHEMRICH_INPUT = "GF_GCMS_my_batch_ChemRICH_R_input.xlsx"


"""
Import data into pandas dataframe
"""
gf_gnps_data = pd.read_excel(pjoin(OUTPUT_FOLDER, GF_GNPS_DATA_FILENAME))

"""
Filter for Columns of Interest and No Missing Values
"""
gf_gnps_data = gf_gnps_data[COLS_OF_INTEREST]

# Remove rows with missing values in any of the columns of interest
gf_gnps_data = gf_gnps_data.dropna()

"""
Convert Column Names to Standard Names
"""
gf_gnps_data = gf_gnps_data.rename(columns=COL_NAME_CONVERT_DICT)

"""
Apply Confidence Filter
"""
# Filter for rows with MQ score greater than or equal to MQ_SCORE_CUTOFF
gf_gnps_data = gf_gnps_data[gf_gnps_data[MQ_SCORE_COL_NAME] >= MQ_SCORE_CUTOFF]

# For duplicate rows (based on "compound_name"), keep only rows with this highest MQ score. If there are multiple rows with the same MQ score, keep the first row.
gf_gnps_data = gf_gnps_data.sort_values(by=[MQ_SCORE_COL_NAME], ascending=False)
gf_gnps_data = gf_gnps_data.drop_duplicates(subset=["compound_name"], keep="first")


"""
Use pubchempy to add Pubchem IDs to gf_gnps_data
"""
smiles_list = gf_gnps_data["smiles"].tolist()
smiles_to_cid_dict = {}

for smiles in smiles_list:
    try:
        compound = pcp.get_compounds(smiles, 'smiles')[0]
        smiles_to_cid_dict[smiles] = compound.cid
    except:
        smiles_to_cid_dict[smiles] = "NA"

gf_gnps_data["pubchem_id"] = gf_gnps_data["smiles"].map(smiles_to_cid_dict)


"""
Create mesh_input_df
"""
# Create mesh_input_df with columns in MESH_COL_ORDER from gf_gnps_data
mesh_input_df = gf_gnps_data[MESH_COL_ORDER]

# Export mesh_input_df to excel in output_folder
mesh_input_df.to_excel(pjoin(OUTPUT_FOLDER, OUTPUT_FILENAME_MESH_INPUT), index=False)


"""
Pause: Run MESH generation in Script 7 (in R)*************************************************************
"""

"""
Part 2*************************************************************
"""
# Import MeSh Prediction Results
mesh_cmpd_classes = pd.read_excel(pjoin(OUTPUT_FOLDER, MESH_CMPD_CLASSES_FILENAME))

# Use add_to_df function to add MeSh Prediction Results to gf_gnps_data, using "compound_name" as the key column
gf_gnps_data = add_to_df(gf_gnps_data, mesh_cmpd_classes, [MESH_PREDICTIONS_COL_NAME], "compound_name")

# Use COL_NAME_CONVERT_DICT to convert me names to standard names ("set")
gf_gnps_data = gf_gnps_data.rename(columns=COL_NAME_CONVERT_DICT)

# Reorder columns in gf_gnps_data to FINAL_COL_ORDER
gf_gnps_data = gf_gnps_data[FINAL_COL_ORDER]

# Print number of rows prior to filtering for MeSh Predictions
print(f"Number of rows in gf_gnps_data before filtering for MeSh Predictions: {len(gf_gnps_data)}")

# Remove rows without MeSh Predictions
gf_gnps_data = gf_gnps_data.dropna(subset=["set"])

# Print how many rows remain
print(f"Number of rows in gf_gnps_data after filtering for MeSh Predictions: {len(gf_gnps_data)}")

# Export gf_gnps_data to excel in output_folder
gf_gnps_data.to_excel(pjoin(OUTPUT_FOLDER, OUTPUT_FILENAME_CHEMRICH_INPUT), index=False)