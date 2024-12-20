"""
GF GCMS Data Analysis Script 6
Lazarina Butkovich 12/20/24

This script generates plots to visualize differences between AR and CC metabolite profiles.
Using normalized peak intensities from MSDIAL_stats.xlsx.

Outputs:
- FA (similar to PCA) plot showing AR vs CC sample grouping
- Loadings plot showing metabolites that contribute most to AR vs CC differences
"""

import pandas as pd
import numpy as np
import os
from os.path import join as pjoin
import matplotlib.pyplot as plt
np.float_ = np.float64
np.int_ = np.int64
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import FactorAnalysis  # Change from ProbabilisticPCA

"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
OUTPUT_FOLDER = r'output'

# Input files
FILENAME_DATA = 'GF_GCMS_stats_summary_table.xlsx'

# AR sample columns
AR_COLS = ['OMALL_RFS_AR_S4_1_M', 'OMALL_RFS_AR_S4_2_M', 
           'OMALL_RFS_AR_S4_3_M', 'OMALL_RFS_AR_S4_4_M']

# CC sample columns 
CC_COLS = ['OMALL_RFS_CC1_M', 'OMALL_RFS_CC2_M', 
           'OMALL_RFS_CC3_M', 'OMALL_RFS_CC4_M']

# Number of top loading features to label
N_TOP_FEATURES = 10

"""""""""""""""""""""""""""""""""""""""""""""
Main
"""""""""""""""""""""""""""""""""""""""""""""
"""
Parse Data
"""
# Import data
data_df = pd.read_excel(pjoin(OUTPUT_FOLDER, FILENAME_DATA), sheet_name='Summary Table')


# Extract AR and CC data
ar_data = data_df[AR_COLS]
cc_data = data_df[CC_COLS]
gnps_id_data = data_df['Compound_Name_GNPS']

# Combine data
combined_data = pd.concat([ar_data, cc_data], axis=1)

# Normalize the data
scaler = StandardScaler()
normalized_data = scaler.fit_transform(combined_data.T).T

"""
Perform Factor Analysis
"""
# Perform Factor Analysis instead of pPCA
fa = FactorAnalysis(n_components=2, random_state=42)
fa_results = fa.fit_transform(normalized_data.T)

# Create a DataFrame for FA results
fa_df = pd.DataFrame(fa_results, columns=['Factor1', 'Factor2'])
fa_df['Sample'] = ['AR'] * len(AR_COLS) + ['CC'] * len(CC_COLS)

# Get loadings
loadings = fa.components_.T

# Create loadings DataFrame
loadings_df = pd.DataFrame(loadings, columns=['Loadings1', 'Loadings2'])
loadings_df['Compound_Name_GNPS'] = gnps_id_data

# Calculate variance explained (approximate for FA)
var_explained = np.var(fa_results, axis=0) / np.sum(np.var(normalized_data.T, axis=0)) * 100

"""
Plot Factor Analysis
"""
# Plot FA results
plt.figure(figsize=(10, 7))
for sample, color in zip(['AR', 'CC'], ['darkblue', 'lightgreen']):
    subset = fa_df[fa_df['Sample'] == sample]
    plt.scatter(subset['Factor1'], subset['Factor2'], label=sample, color=color, alpha=0.7)
plt.xlabel(f'Loadings Factor1 ({var_explained[0]:.1f}% variance explained)')
plt.ylabel(f'Loadings Factor2 ({var_explained[1]:.1f}% variance explained)')
plt.title('FA Plot: AR vs CC')
plt.legend()
plt.savefig(pjoin(OUTPUT_FOLDER, 'FA_plot_batch_3.png'))
plt.show()

"""
Plot Loadings
"""
# Plot loadings
plt.figure(figsize=(12, 8))
plt.scatter(loadings_df['Loadings1'], loadings_df['Loadings2'], alpha=0.5)
plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
plt.axvline(x=0, color='k', linestyle='--', alpha=0.3)

# Label top contributing features
top_pos = loadings_df.nlargest(5, 'Loadings1')
top_neg = loadings_df.nsmallest(5, 'Loadings1')

for idx, row in pd.concat([top_pos, top_neg]).iterrows():
    plt.annotate(row['Compound_Name_GNPS'], 
                (row['Loadings1'], row['Loadings2']),
                xytext=(5, 5), textcoords='offset points')

plt.xlabel(f'Loadings Factor1 ({var_explained[0]:.1f}% variance explained)')
plt.ylabel(f'Loadings Factor2 ({var_explained[1]:.1f}% variance explained)')
plt.title('FA Loadings Plot')
plt.tight_layout()
plt.savefig(pjoin(OUTPUT_FOLDER, 'FA_loadings_plot_batch_3.png'))
plt.show()

# Export loadings_df to Excel
loadings_df.to_excel(pjoin(OUTPUT_FOLDER, 'Loadings_df_batch_3.xlsx'), index=False)

"""
Generate Bar Charts of Top Contributing Metabolites
"""
# Calculate separate contributions for CC and AR
loadings_feature = pd.DataFrame(
    fa.components_.T,
    columns=['Factor1', 'Factor2'],
    index=gnps_id_data
)

# Separate positive (CC) and negative (AR) loadings on Factor1
cc_features = loadings_feature[loadings_feature['Factor1'] > 0]
ar_features = loadings_feature[loadings_feature['Factor1'] < 0]

# Calculate contributions
cc_contributions = (cc_features ** 2).sum(axis=1) * 100
ar_contributions = (ar_features ** 2).sum(axis=1) * 100

# Get top 20 for each
top_cc = cc_contributions.nlargest(30)
top_ar = ar_contributions.nlargest(30)

# Create subplot figure
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot CC contributions
ax1.bar(range(len(top_cc)), top_cc.values, color='lightgreen')
ax1.set_xticks(range(len(top_cc)))
ax1.set_xticklabels(top_cc.index, rotation=45, ha='right')
ax1.set_ylabel('% Contribution')
ax1.set_title('Top 30 Contributing Metabolites in CC Samples')

# Plot AR contributions
ax2.bar(range(len(top_ar)), top_ar.values, color='darkblue')
ax2.set_xticks(range(len(top_ar)))
ax2.set_xticklabels(top_ar.index, rotation=45, ha='right')
ax2.set_ylabel('% Contribution')
ax2.set_title('Top 30 Contributing Metabolites in AR Samples')

# Adjust layout and save
plt.tight_layout()
plt.savefig(pjoin(OUTPUT_FOLDER, 'top_30_metabolites_by_group.png'))
plt.show()

# Export contributions to Excel
contributions_df = pd.DataFrame({
    'CC_Metabolite': top_cc.index,
    'CC_Contribution': top_cc.values,
    'AR_Metabolite': top_ar.index, 
    'AR_Contribution': top_ar.values
})
contributions_df.to_excel(pjoin(OUTPUT_FOLDER, 'FA_metabolite_contributions_by_group.xlsx'), index=False)