"""
Gut Fungal GC-MS Profiling, Script 4: PCA and Heatmap Plots
Lazarina Butkovich 2024

This script performs various analyses and creates visualizations for GC-MS data:
1. Principal Component Analysis (PCA) and probabilistic PCA (pPCA)
2. ANOVA analysis with FDR correction 
3. Heatmaps of metabolite abundances
4. Individual metabolite bar plots

Key plots generated:
- pPCA and PCA plots comparing all samples
- PCA comparing A. robustus vs C. churrovis
- Hierarchically clustered heatmaps of metabolite abundances
- Bar plots for metabolites of interest
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

    # Set font sizes
    plt.rcParams.update({
        'font.size': 14,
        'axes.linewidth': 2
    })

    # Create clustermap with increased linewidths
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
        figsize=(8, len(z_scores)*0.3 + 2),
        linewidths=0,  # Remove lines between cells
        tree_kws={'linewidths': 2.0}  # Increase dendrogram line width
    )
    
    # Adjust styling
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=12)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=12)
    
    # Increase tick width
    g.ax_heatmap.tick_params(width=2)
    
    # Rotate labels 180 degrees
    # g.ax_heatmap.yaxis.set_label_text('Metabolite', rotation=270)
    g.ax_cbar.set_ylabel('Z-score', rotation=270)
    
    # Adjust colorbar and font size
    g.ax_cbar.tick_params(labelsize=14, width=2)  # Increased from 10 to 14
    
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

def create_barplot_grid(data_knowns, barplot_grid_order, sample_groups, COLORS, FULL_NAMES, ITALICIZE_NAMES):
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    
    # Set consistent font sizes
    plt.rcParams.update({'font.size': 12})
    
    # Track missing metabolites
    missing_metabolites = []
    
    # Create subplots for each metabolite
    for i, metabolite in enumerate(barplot_grid_order):
        ax = plt.subplot(3, 4, i+1)
        
        # Filter data for current metabolite
        data_metabolite = data_knowns[data_knowns[CMPD_COL_NAME] == metabolite]
        
        if len(data_metabolite) == 0:
            missing_metabolites.append(metabolite)
            continue
            
        # Plot bars for each sample group
        bar_positions = np.arange(len(sample_groups))
        bar_width = 0.7
        for j, (sample_type, color) in enumerate(COLORS.items()):
            if sample_type in sample_groups:
                values = data_metabolite[sample_groups[sample_type]].values[0]
                mean = np.mean(values)
                std = np.std(values)
                ax.bar(j, mean, width=bar_width, color=color)
                ax.errorbar(j, mean, yerr=std, color='black', capsize=5, linewidth=1)
        
        # Customize plot style
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.tick_params(axis='y', direction='in')
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
        # Add title
        ax.set_title(metabolite, fontsize=10, pad=5)

    # Create legend in the last subplot
    ax = plt.subplot(3, 4, 12)
    ax.axis('off')
    
    # Add legend items
    legend_elements = []
    legend_labels = []
    for sample_type, color in COLORS.items():
        if sample_type in sample_groups:
            patch = plt.Rectangle((0,0), 1, 1, fc=color)
            legend_elements.append(patch)
            label = FULL_NAMES[sample_type]
            if ITALICIZE_NAMES[sample_type]:
                label = r'$\mathit{' + label + '}$'
            legend_labels.append(label)
    
    ax.legend(legend_elements, legend_labels, 
              loc='center', 
              frameon=False,
              fontsize=14,  # Increased from 12 to 14
              prop={'size': 14})  # Added explicit size setting

    # Adjust layout
    plt.tight_layout()
    
    return fig, missing_metabolites

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

# colors
COLORS = {'AR':'darkblue', 'CC':'lightgreen', 'MC':'lightgrey', 'RF':'dimgrey', 'FAMES':'pink', 'BLANK':'olive'}


# For metabolite bar charts
METABOLITES_OF_INTEREST_LIST = [
    'Myo-inositol',
    'Propylene glycol',
    'Lactic acid',
    '4-aminobutyric acid (GABA)',
    'Iminodiacetic acid',
    '4-hydroxy-3-methoxybenzoic acid (isovanillic acid)',
    '4-hydroxybenzoic acid (p-salicylic acid)',
    'Benzoic acid',
    '2-hydroxyglutaric acid',
    'D-malic acid',
    'Cellobiose',
    'Maltose',
    'Lactulose',
    'Glyceric acid',
    'D-fructose',
    'Sucrose',
    'Palatinose',
    'Melibiose',
    'Maltotriitol',
    'Maltotriose',
    'D-glucose-6-phosphate',
    'Trehalose',
    'D-xylitol',
    'D-glucose',
    'D-mannose',
    'D-sorbitol',
    'Fumaric acid',
    'Succinic acid',
    '2,3-dihydroxyisovaleric acid',
    'Arachidic acid',
    'Palmitoleic acid',
    'Behenic acid',
    'Lauric acid',
    'Oleic acid',
    'Stearic acid',
    '10-hydroxydecanoic acid',
    'Heptadecanoic acid',
    'Myristic acid',
    'Methyl oleate',
    'Palmitic acid',
    'Capric acid',
    '3-(4-hydroxyphenyl)propionic acid (phloretic acid)',
    'Xanthine',
    'Nicotinamide',
    '2,3-dihydroxypyridine',
    '2-hydroxypyridine',
    'Orotic acid',
    'p-cresol',
    'Threose',
    'DL-dihydrosphingosine'
    ]

FONT_SIZE = 20


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
Generate ANOVA analysis and heatmaps
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
# Iterate through the metabolite classes and create heatmaps for each. Include metabolites of only high confidence. Remove RF samples.
# For the metabolite class named 'Amino acid/peptide', rename it so that the / does not cause issues with saving the file
data_knowns.loc[data_knowns['Metabolite Class'] == 'Amino acid/peptide', 'Metabolite Class'] = 'Amino acid-peptide'
# Remove low confidence metabolites
data_knowns_high_conf = data_knowns[data_knowns['Confidence'] != 'low']
metabolite_classes = data_knowns_high_conf['Metabolite Class'].unique()

# If a metabolite class has 2 or less rows, remove it from consideration for heatmap generation
metabolite_classes = [m for m in metabolite_classes if data_knowns_high_conf[data_knowns_high_conf['Metabolite Class'] == m].shape[0] > 2]

# Create a folder in output folder called 'Metabolite Class Heatmaps', if it does not already exist
metabolite_class_heatmap_folder = pjoin(OUTPUT_FOLDER, 'Metabolite Class Heatmaps')
os.makedirs(metabolite_class_heatmap_folder, exist_ok=True)

# Remove RF samples
data_knowns_no_rf = data_knowns_high_conf.drop(columns=rf_col_names)
sample_groups_no_rf = {k: v for k, v in sample_groups_no_blank.items() if k != 'RF'}

# After removing RF samples and low confidence matches, check if there are any metabolites (rows) with all zero or missing values in any of the remaining sample group columns
# If there are, remove them
data_knowns_no_rf = data_knowns_no_rf.loc[~(data_knowns_no_rf[sample_groups_no_rf['AR'] + sample_groups_no_rf['CC'] + sample_groups_no_rf['MC']].sum(axis=1) == 0)]


for metabolite_class in metabolite_classes:
    # Filter data to only include metabolites in the current class
    data_class = data_knowns_no_rf[data_knowns_no_rf['Metabolite Class'] == metabolite_class]

    # Create and save heatmap
    heatmap_class_fig = create_metabolite_heatmap(data_class, sample_groups_no_rf)
    # plt.show()
    heatmap_class_fig.savefig(pjoin(metabolite_class_heatmap_folder, f'metabolite_heatmap_{metabolite_class}_high_conf_batch_3.png'), 
                            dpi=600, bbox_inches='tight')
    

"""
Create Barplots for Specific Metabolites
"""
# Create folder for bar charts
metabolite_barplot_folder = pjoin(OUTPUT_FOLDER, 'Metabolite Bar Charts')
os.makedirs(metabolite_barplot_folder, exist_ok=True)

# Track missing metabolites
missing_metabolites = []

# For each metabolite in METABOLITES_OF_INTEREST_LIST, create a barplot
# First create a separate legend figure
fig_legend = plt.figure(figsize=(3, 2))
ax_legend = fig_legend.add_subplot(111)

# Create dummy plots for legend
legend_labels = []
handles = []
for sample_type, color in COLORS.items():
    if sample_type in sample_groups:
        handle = ax_legend.plot([], [], 'o', color=color, markersize=10)[0]
        label = FULL_NAMES[sample_type]
        if ITALICIZE_NAMES[sample_type]:
            label = r'$\mathit{' + label + '}$'
        legend_labels.append(label)
        handles.append(handle)

ax_legend.axis('off')
legend = ax_legend.legend(handles, legend_labels, loc='center')
fig_legend.savefig(pjoin(metabolite_barplot_folder, 'legend.png'), 
                   dpi=600, bbox_inches='tight')
plt.close(fig_legend)

for metabolite in METABOLITES_OF_INTEREST_LIST:
    # Filter data to only include the current metabolite
    data_metabolite = data_knowns[data_knowns[CMPD_COL_NAME] == metabolite]
    
    # Skip if metabolite not found
    if len(data_metabolite) == 0:
        missing_metabolites.append(metabolite)
        continue
    
    # Print a note if the metabolite is low confidence
    if data_metabolite['Confidence'].values[0] == 'low':
        print(f"Metabolite '{metabolite}' is low confidence")

    # Create barplot with reduced width
    fig, ax = plt.subplots(figsize=(4, 6))  # Reduced width from 6 to 4
    
    # Plot metabolite composition for each sample group
    bar_positions = np.arange(len(sample_groups))
    bar_width = 0.7  # Increased bar width from 0.5 to 0.7
    for i, (sample_type, color) in enumerate(COLORS.items()):
        if sample_type in sample_groups:
            values = data_metabolite[sample_groups[sample_type]].values[0]
            mean = np.mean(values)
            std = np.std(values)
            ax.bar(i, mean, width=bar_width, color=color)
            ax.errorbar(i, mean, yerr=std, color='black', capsize=5, linewidth=2)
    
    # Remove border and keep only y-axis
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # Remove x-axis labels
    ax.set_xticks([])
    plt.rcParams.update({'font.size': FONT_SIZE})
    ax.tick_params(axis='y', labelsize=FONT_SIZE, length = 8, width = 2)
    ax.title.set_size(FONT_SIZE)
    
    # Set y-axis to have 5 ticks and face inward
    ax.tick_params(axis='y', direction='in')
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    
    ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=False))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    # # Add title
    # ax.set_title(metabolite)
    
    plt.tight_layout()
    
    # Save plot without legend
    fig.savefig(pjoin(metabolite_barplot_folder, f'metabolite_barplot_{metabolite}_batch_3.png'), 
                dpi=600, bbox_inches='tight')
    plt.close(fig)

# Print missing metabolites
if missing_metabolites:
    print("\nMetabolites not found in dataset:")
    for met in missing_metabolites:
        print(f"- {met}")


"""
Export Barplot Grid .png
"""
# create a grid of barplots in a 3 by 4 grid. The 12th plot will be the legend.
barplot_grid_order = ["D-malic acid", "Fumaric acid", "Glyceric acid", "Succinic acid", "D-glucose-6-phosphate", "D-mannose", "Lactulose", "DL-dihydrosphingosine", "3-(4-hydroxyphenyl)propionic acid (phloretic acid)", "4-hydroxy-3-methoxybenzoic acid (isovanillic acid)", "4-hydroxybenzoic acid (p-salicylic acid)"]

# Create and save the grid figure
grid_fig, missing = create_barplot_grid(data_knowns, barplot_grid_order, sample_groups, COLORS, FULL_NAMES, ITALICIZE_NAMES)
grid_fig.savefig(pjoin(OUTPUT_FOLDER, 'metabolite_barplots_grid_batch_3.png'), 
                 dpi=600, bbox_inches='tight')

# Print any missing metabolites
if missing:
    print("\nMetabolites not found in grid dataset:")
    for met in missing:
        print(f"- {met}")