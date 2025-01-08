"""
GF GCMS Data Analysis Script 4a: pPCA and Heatmap Plots for Batch 3
Lazarina Butkovich 12/20/24
"""

import os

import pandas as pd
import numpy as np
from os.path import join as pjoin
import matplotlib.pyplot as plt
# Use probabilistic PCA
from ppca import PPCA
from scipy import stats
import seaborn as sns
from statsmodels.stats.multitest import fdrcorrection
from scipy.cluster import hierarchy
# Import pca module from sklearn
from sklearn.decomposition import PCA

"""
Functions
"""
def normalize_by_tic(data, sample_groups):
    """
    Normalize each data column by the total ion chromatogram, or the sum of all the values in the column.
    """
    # Normalize each sample_groups column by the sum of the column
    for sample_type, cols in sample_groups.items():
        data[cols] = data[cols].div(data[cols].sum(), axis=1)

    return data

def create_ppca_plot(data, sample_groups, colors, title="pPCA Analysis"):
    # Prepare data
    all_cols = []
    sample_labels = []
    for sample_type, cols in sample_groups.items():
        all_cols.extend(cols)
        sample_labels.extend([sample_type] * len(cols))

    # Convert to numpy array for proper indexing 
    sample_labels = np.array(sample_labels)
    
    # Extract and transpose data for PPCA (samples as rows)
    X = data[all_cols].T.values
    
    # Fit PPCA; d=2 for reduction to 2 dimensions
    ppca = PPCA()
    ppca.fit(data=X, d=2)
    transformed = ppca.transform()
        
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Plot each sample group with full names in legend
    for sample_type in np.unique(sample_labels):
        mask = (sample_labels == sample_type)
        plt.scatter(transformed[mask, 0], 
                   transformed[mask, 1],
                   c=colors[sample_type],
                   label=FULL_NAMES[sample_type], # Use full name from FULL_NAMES dict
                   alpha=0.7)

    # Update axis labels with variance explained
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    
    plt.title(title)
    plt.grid(True, alpha=0.3)
    
    # Add legend with italic font for species names
    legend = plt.legend(prop={'style': 'italic'})
    
    legend_labels = []
    handles = []
    for handle, label in zip(*plt.gca().get_legend_handles_labels()):
        # Get the original sample type from the full name
        sample_type = next(k for k, v in FULL_NAMES.items() if v == label)
        if ITALICIZE_NAMES[sample_type]:
            label = r'$\mathit{' + label + '}$'  # Use \mathit for proper LaTeX italics
        legend_labels.append(label)
        handles.append(handle)
        
    plt.legend(handles, legend_labels)
    
    plt.tight_layout()
    
    return plt.gcf()

def analyze_metabolites(data, sample_groups):
    """Perform ANOVA and calculate FDR-corrected q-values for each metabolite"""
    results = []
    
    # Get all samples except blanks
    groups_no_blank = {k:v for k,v in sample_groups.items() if k != 'BLANK'}
    
    # Store p-values for FDR correction
    p_values = []
    metabolites = []
    group_means = []
    
    for metabolite in data[CMPD_COL_NAME]:
        # Prepare groups for ANOVA
        groups = []
        for group in groups_no_blank.values():
            groups.append(data.loc[data[CMPD_COL_NAME] == metabolite, group].values[0])
        
        # Perform one-way ANOVA
        f_stat, p_val = stats.f_oneway(*groups)
        
        # Calculate mean values for each group
        means = {}
        for group_name, group_cols in groups_no_blank.items():
            means[group_name] = data.loc[data[CMPD_COL_NAME] == metabolite, group_cols].values[0].mean()
        
        metabolites.append(metabolite)
        p_values.append(p_val)
        group_means.append(means)
    
    # Calculate FDR-corrected q-values
    rejected, q_values = fdrcorrection(p_values, alpha=0.05, method='indep')
    
    # Combine results
    for metabolite, p_val, q_val, means in zip(metabolites, p_values, q_values, group_means):
        results.append({
            'Metabolite': metabolite,
            'p_value': p_val,
            'q_value': q_val,
            **means
        })
    
    return pd.DataFrame(results)

def create_metabolite_heatmap(data, sample_groups, cmpd_col='Metabolite'):
    # Get mean values for each sample group
    group_means = pd.DataFrame()
    for group, cols in sample_groups.items():
        group_means[group] = data[cols].mean(axis=1)
    
    # Set metabolite names as index
    group_means.index = data[cmpd_col]
    
    # Calculate z-scores
    z_scores = pd.DataFrame(
        stats.zscore(group_means, axis=1),
        columns=group_means.columns,
        index=group_means.index
    )

    # Create xticklabels with italicized names where necessary
    xticklabels = []
    for col in z_scores.columns:
        if ITALICIZE_NAMES[col]:
            xticklabels.append(r'$\mathit{' + FULL_NAMES[col] + '}$')
        else:
            xticklabels.append(FULL_NAMES[col])

    # Create clustermap
    g = sns.clustermap(
        z_scores,
        cmap='RdBu_r',
        center=0,
        robust=True,
        cbar_kws={'label': 'Z-score'},
        xticklabels=xticklabels,
        yticklabels=True,
        row_cluster=True,
        col_cluster=False,
        dendrogram_ratio=0.2,
        figsize=(8, len(z_scores)*0.3 + 2)
    )
    
    # Adjust styling
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    
    return g.figure

def create_pca_plot(data, sample_groups, colors, title="PCA Analysis"):
    # Prepare data
    all_cols = []
    sample_labels = []
    for sample_type, cols in sample_groups.items():
        all_cols.extend(cols)
        sample_labels.extend([sample_type] * len(cols))

    # Convert to numpy array
    sample_labels = np.array(sample_labels)
    
    # Extract and transpose data (samples as rows)
    X = data[all_cols].T.values
    
    # Fit PCA
    pca = PCA(n_components=2)
    transformed = pca.fit_transform(X)

    # Create plot    
    plt.figure(figsize=(10, 8))
    
    # Plot each sample group
    for sample_type in np.unique(sample_labels):
        mask = (sample_labels == sample_type)
        plt.scatter(transformed[mask, 0], 
                   transformed[mask, 1],
                   c=colors[sample_type],
                   label=FULL_NAMES[sample_type],
                   alpha=0.7)

    # Add variance explained to labels
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    
    plt.title(title)
    plt.grid(True, alpha=0.3)
    
    legend = plt.legend(prop={'style': 'italic'})
    
    legend_labels = []
    handles = []
    for handle, label in zip(*plt.gca().get_legend_handles_labels()):
        # Get the original sample type from the full name
        sample_type = next(k for k, v in FULL_NAMES.items() if v == label)
        if ITALICIZE_NAMES[sample_type]:
            label = r'$\mathit{' + label + '}$'  # Use \mathit for proper LaTeX italics
        legend_labels.append(label)
        handles.append(handle)
        
    plt.legend(handles, legend_labels)
    
    plt.tight_layout()
    
    return plt.gcf()

"""
Values
"""
SAMPLE_NAME_PRE_POST_STRS_DICT = {'CC':('OMALL_RFS_CC','_M'),'AR':('OMALL_RFS_AR_S4_','_M'),'MC':('OMALL_RFS_MC','_M'),'RF':('OMALL_RFS_RF','_M'), 'FAMES':('GCMS_FAMES_0','_GCMS01_20201209'), 'BLANK':('GCMS_BLANK_0','_GCMS01_20201209')}
REPLICATE_NUMS = {'CC':4, 'AR':4, 'MC':4, 'RF':4,'FAMES':1,'BLANK':3}

CMPD_COL_NAME = 'Metabolite'

DATA_FILENAME = 'Compound_ids_PNNL.xlsm'

INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'

# Sample colors
COLORS = {'CC':'lightgreen', 'AR':'darkblue', 'MC':'tan', 'RF':'dimgrey', 'FAMES':'pink', 'BLANK':'olive', 'G1': 'lightblue', 'S3': 'peru', 'PF': 'khaki'}

# AR and CC full names should be italicized in figure legends and labels
FULL_NAMES = {
    'AR': 'A. robustus',
    'CC': 'C. churrovis',
    'MC': 'Medium C',
    'RF': 'Rumen Fluid',
    'BLANK': 'Blank'
}

ITALICIZE_NAMES = {
    'AR': True,  # A. robustus should be italicized
    'CC': True,  # C. churrovis should be italicized
    'MC': False, # Medium C should not be italicized
    'RF': False, # Rumen Fluid should not be italicized
    'BLANK': False
}

"""
Import data
"""
data = pd.read_excel(pjoin(INPUT_FOLDER, DATA_FILENAME))

"""
Generate Sample Names
"""
# Generate column names, using the sample name prefixes and replicate numbers
ar_col_names = [f'{SAMPLE_NAME_PRE_POST_STRS_DICT["AR"][0]}{i}{SAMPLE_NAME_PRE_POST_STRS_DICT["AR"][1]}' for i in range(1, REPLICATE_NUMS['AR']+1)]
cc_col_names = [f'{SAMPLE_NAME_PRE_POST_STRS_DICT["CC"][0]}{i}{SAMPLE_NAME_PRE_POST_STRS_DICT["CC"][1]}' for i in range(1, REPLICATE_NUMS['CC']+1)]
mc_col_names = [f'{SAMPLE_NAME_PRE_POST_STRS_DICT["MC"][0]}{i}{SAMPLE_NAME_PRE_POST_STRS_DICT["MC"][1]}' for i in range(1, REPLICATE_NUMS['MC']+1)]
rf_col_names = [f'{SAMPLE_NAME_PRE_POST_STRS_DICT["RF"][0]}{i}{SAMPLE_NAME_PRE_POST_STRS_DICT["RF"][1]}' for i in range(1, REPLICATE_NUMS['RF']+1)]

# Define sample groups and colors
sample_groups = {
    'AR': ar_col_names,
    'CC': cc_col_names, 
    'MC': mc_col_names,
    'RF': rf_col_names,
}

"""
Normalize data by TIC
"""
# Normalize data by TIC
data = normalize_by_tic(data, sample_groups)

"""
Generate pPCA plot
"""
# Generate pPCA plot
fig = create_ppca_plot(data, sample_groups, COLORS)
plt.show()
# Save plot
fig.savefig(pjoin(OUTPUT_FOLDER, 'ppca_plot_batch_3.png'), dpi=600, bbox_inches='tight')

# """
# pPCA for AR and CC samples only --> causes LinAlgError
# """
# # Create a pPCA plot with only 'AR' and 'CC' samples
# # Create sample groups with only AR and CC
# ar_cc_groups = {k: v for k, v in sample_groups.items() if k in ['AR', 'CC']}

# # Generate pPCA plot for AR and CC only
# ar_cc_fig = create_ppca_plot(data, ar_cc_groups, COLORS, title="pPCA Analysis: A. robustus vs C. churrovis")
# plt.show()

# # Save plot
# ar_cc_fig.savefig(pjoin(OUTPUT_FOLDER, 'ppca_plot_ar_cc_batch_3.png'), dpi=600, bbox_inches='tight')

"""
PCA for AR and CC samples only
"""
# Create a PCA plot with only 'AR' and 'CC' samples
# Create sample groups with only AR and CC
ar_cc_groups = {k: v for k, v in sample_groups.items() if k in ['AR', 'CC']}
# Generate PCA plot for AR and CC only
pca_ar_cc_fig = create_pca_plot(data, ar_cc_groups, COLORS, title="PCA Analysis: A. robustus vs C. churrovis")
plt.show()

# Save plot
pca_ar_cc_fig.savefig(pjoin(OUTPUT_FOLDER, 'pca_plot_ar_cc_batch_3.png'), dpi=600, bbox_inches='tight')

"""
PCA for all samples
"""
# Generate PCA plot for all samples
pca_all_fig = create_pca_plot(data, sample_groups, COLORS)
plt.show()

# Save plot
pca_all_fig.savefig(pjoin(OUTPUT_FOLDER, 'pca_plot_all_batch_3.png'), dpi=600, bbox_inches='tight')


"""
Generate ANOVA analysis and heatmap
"""
# Create data_knowns, the data df filtered to remove rows with 'Unknown...' in the CMPD_COL_NAME column
data_knowns = data[~data[CMPD_COL_NAME].str.contains('Unknown')]

# Check that no rows have the same metabolite, if they do, end the code and print warning message
if data_knowns[CMPD_COL_NAME].duplicated().any():
    print("Warning: Duplicate metabolites detected in data_knowns. Please check the data. Duplicated metabolites:")
    print(data_knowns[data_knowns[CMPD_COL_NAME].duplicated(keep=False)][CMPD_COL_NAME])
    exit()

# Perform ANOVA analysis on all sample groups except BLANK (note, in this script BLANK is not included in the analysis)
sample_groups_no_blank = {k: v for k, v in sample_groups.items() if k != 'BLANK'}

anova_results = analyze_metabolites(data_knowns, sample_groups_no_blank)

# Save ANOVA results
anova_results.to_csv(pjoin(OUTPUT_FOLDER, 'anova_results_batch_3.csv'), index=False)


# Create and save heatmap
heatmap_fig = create_metabolite_heatmap(data_knowns, sample_groups_no_blank)
plt.show()
heatmap_fig.savefig(pjoin(OUTPUT_FOLDER, 'metabolite_heatmap_batch_3.png'), 
                    dpi=600, bbox_inches='tight')


# Create heatmap with low confidence metabolites removed (column 'Confidence' value is 'low')
data_high_confidence = data_knowns[data_knowns['Confidence'] != 'low']

heatmap_high_confidence_fig = create_metabolite_heatmap(data_high_confidence, sample_groups_no_blank)
plt.show()
heatmap_high_confidence_fig.savefig(pjoin(OUTPUT_FOLDER, 'metabolite_heatmap_high_confidence_batch_3.png'), 
                                    dpi=600, bbox_inches='tight')


# Create heatmap with only high confidence metabolites and no RF samples
data_high_confidence_no_rf = data_high_confidence.drop(columns=rf_col_names)
sample_groups_no_rf = {k: v for k, v in sample_groups_no_blank.items() if k != 'RF'}

# After removing RF samples, check if there are any metabolites (rows) with all zero or missing values in any of the remaining sample group columns
# If there are, remove them
data_high_confidence_no_rf = data_high_confidence_no_rf.loc[~(data_high_confidence_no_rf[sample_groups_no_rf['AR'] + sample_groups_no_rf['CC'] + sample_groups_no_rf['MC']].sum(axis=1) == 0)]

# Create and save heatmap
heatmap_high_confidence_no_rf_fig = create_metabolite_heatmap(data_high_confidence_no_rf, sample_groups_no_rf)
plt.show()
heatmap_high_confidence_no_rf_fig.savefig(pjoin(OUTPUT_FOLDER, 'metabolite_heatmap_high_confidence_no_rf_batch_3.png'), 
                                        dpi=600, bbox_inches='tight')

# Lastly, create heatmap with both low and high confidence metabolites but without RF samples
data_no_rf = data_knowns.drop(columns=rf_col_names)
sample_groups_no_rf = {k: v for k, v in sample_groups_no_blank.items() if k != 'RF'}

# After removing RF samples, check if there are any metabolites (rows) with all zero or missing values in any of the remaining sample group columns
# If there are, remove them
data_no_rf = data_no_rf.loc[~(data_no_rf[sample_groups_no_rf['AR'] + sample_groups_no_rf['CC'] + sample_groups_no_rf['MC']].sum(axis=1) == 0)]

# Create and save heatmap
heatmap_no_rf_fig = create_metabolite_heatmap(data_no_rf, sample_groups_no_rf)
plt.show()
heatmap_no_rf_fig.savefig(pjoin(OUTPUT_FOLDER, 'metabolite_heatmap_no_rf_batch_3.png'), 
                        dpi=600, bbox_inches='tight')

"""
Heatmaps for Metabolite Class Subsets
"""
# Iterate through the metabolite classes and create heatmaps for each. Include metabolites of low or high confidence. Remove RF samples.
# For the metabolite class named 'Amino acid/peptide', rename it so that the / does not cause issues with saving the file
data_knowns.loc[data_knowns['Metabolite Class'] == 'Amino acid/peptide', 'Metabolite Class'] = 'Amino acid-peptide'
metabolite_classes = data_knowns['Metabolite Class'].unique()

# If a metabolite class has 2 or less rows, remove it from consideration for heatmap generation
metabolite_classes = [m for m in metabolite_classes if data_knowns[data_knowns['Metabolite Class'] == m].shape[0] > 2]

# Create a folder in output folder called 'Metabolite Class Heatmaps', if it does not already exist
metabolite_class_heatmap_folder = pjoin(OUTPUT_FOLDER, 'Metabolite Class Heatmaps')
os.makedirs(metabolite_class_heatmap_folder, exist_ok=True)

# Remove RF samples
data_knowns_no_rf = data_knowns.drop(columns=rf_col_names)
sample_groups_no_rf = {k: v for k, v in sample_groups_no_blank.items() if k != 'RF'}

# After removing RF samples, check if there are any metabolites (rows) with all zero or missing values in any of the remaining sample group columns
# If there are, remove them
data_knowns_no_rf = data_knowns_no_rf.loc[~(data_knowns_no_rf[sample_groups_no_rf['AR'] + sample_groups_no_rf['CC'] + sample_groups_no_rf['MC']].sum(axis=1) == 0)]

for metabolite_class in metabolite_classes:
    # Filter data to only include metabolites in the current class
    data_class = data_knowns_no_rf[data_knowns_no_rf['Metabolite Class'] == metabolite_class]

    # Create and save heatmap
    heatmap_class_fig = create_metabolite_heatmap(data_class, sample_groups_no_rf)
    # plt.show()
    heatmap_class_fig.savefig(pjoin(metabolite_class_heatmap_folder, f'metabolite_heatmap_{metabolite_class}_batch_3.png'), 
                            dpi=600, bbox_inches='tight')
    