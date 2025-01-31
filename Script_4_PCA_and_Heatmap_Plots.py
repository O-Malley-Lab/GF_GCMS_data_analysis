"""
Gut Fungal GC-MS Profiling, Script 4: PCA and Heatmap Plots
Lazarina Butkovich 2024

This script performs various analyses and creates visualizations for GC-MS data:
1. Principal Component Analysis (PCA) and probabilistic PCA (pPCA)
2. ANOVA analysis with FDR correction 
3. Heatmaps of metabolite abundances
4. Individual metabolite bar plots

- Compound matches were generated by NIST 14 analysis and in-house standards

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
from sklearn.cross_decomposition import PLSRegression
from adjustText import adjust_text

"""
Functions
"""
# def normalize_by_tic(data, sample_groups):
#     """
#     Normalize each data column by the total ion chromatogram, or the sum of all the values in the column.
#     """
#     # Normalize each sample_groups column by the sum of the column
#     for sample_type, cols in sample_groups.items():
#         data[cols] = data[cols].div(data[cols].sum(), axis=1)

#     return data

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
        
    # Create figure with larger size
    plt.figure(figsize=(12, 8))
    
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

    # Create figure with larger size
    plt.figure(figsize=(12, 8))
    
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
    fig = plt.figure(figsize=(9, 15))  # Swapped dimensions to match new grid layout
    
    # Set consistent font sizes
    plt.rcParams.update({'font.size': 12})
    
    # Track missing metabolites
    missing_metabolites = []
    
    # Create subplots for each metabolite
    for i, metabolite in enumerate(barplot_grid_order):
        ax = plt.subplot(5, 3, i+1)
        
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
        
        # Add title with simplified names
        title = metabolite
        if metabolite == "3-(4-hydroxyphenyl)propionic acid (phloretic acid)":
            title = "Phloretic acid"
        elif metabolite == "4-hydroxy-3-methoxybenzoic acid (isovanillic acid)":
            title = "Isovanillic acid"
        elif metabolite == "4-hydroxybenzoic acid (p-salicylic acid)":
            title = "p-salicylic acid"
            
        ax.set_title(title, fontsize=12, pad=5)

    # Create legend in the last subplot
    ax = plt.subplot(5, 3, 15)
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
              fontsize=14,
              prop={'size': 14})

    # Adjust layout
    plt.tight_layout()
    
    return fig, missing_metabolites

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

# Volcano plot cutoff for coefficient of variation (CV) when running t-test (avoid running t-test on metabolites with nearly identical values across all samples, as results are not meaningful/reliable)
VOLCANO_PLOT_CV_CUTOFF = 0.1

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
    'DL-dihydrosphingosine',
    'Loganin'
    ]

FONT_SIZE = 20
LINE_WIDTH = 2

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


# """
# Normalize data by TIC --> Do not normalize by TIC because this dataset was already previously normalized by global median for each sample
# """
# # Normalize data by TIC
# data = normalize_by_tic(data, sample_groups)


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


# """
# Partial Least Squares Discriminant Analysis (PLS-DA) and Loadings Plot
# """
# # Prepare data for PLS-DA
# ar_cc_groups = {k: v for k, v in sample_groups.items() if k in ['AR', 'CC']}

# # Get all sample columns and labels
# all_cols = []
# sample_labels = []
# for sample_type, cols in ar_cc_groups.items():
#     all_cols.extend(cols)
#     sample_labels.extend([sample_type] * len(cols))

# # Convert to numpy arrays
# X = data[all_cols].T.values
# y = np.array([1 if label == 'AR' else 0 for label in sample_labels])

# # Fit PLS-DA model
# plsda = PLSRegression(n_components=2)
# X_transformed = plsda.fit_transform(X, y)[0]

# # Create figure with two subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# # Plot scores (samples)
# for sample_type in ar_cc_groups:
#     mask = (np.array(sample_labels) == sample_type)
#     ax1.scatter(X_transformed[mask, 0], 
#                 X_transformed[mask, 1],
#                 c=COLORS[sample_type],
#                 label=FULL_NAMES[sample_type],
#                 alpha=0.7)

# ax1.set_xlabel('Component 1')
# ax1.set_ylabel('Component 2')
# ax1.set_title('PLS-DA Scores Plot: A. robustus vs C. churrovis')
# ax1.grid(True, alpha=0.3)

# # Add italicized legend to scores plot
# legend_labels = []
# handles = []
# for handle, label in zip(*ax1.get_legend_handles_labels()):
#     sample_type = next(k for k, v in FULL_NAMES.items() if v == label)
#     if ITALICIZE_NAMES[sample_type]:
#         label = r'$\mathit{' + label + '}$'
#     legend_labels.append(label)
#     handles.append(handle)
# ax1.legend(handles, legend_labels)

# # Plot loadings (features)
# loadings = plsda.x_loadings_
# metabolites = data[CMPD_COL_NAME]

# # Filter out unknowns from metabolites
# known_mask = ~data[CMPD_COL_NAME].str.contains('Unknown', na=False)
# known_loadings = loadings[known_mask]
# known_metabolites = metabolites[known_mask]

# # Calculate distance from origin for each loading
# distances = np.sqrt(known_loadings[:, 0]**2 + known_loadings[:, 1]**2)

# # Get indices of top 40 most important features
# top_n = 40
# top_indices = np.argsort(distances)[-top_n:]

# # Map back to original indices for plotting
# top_indices = np.where(known_mask)[0][top_indices]

# # Plot all loadings as small points
# ax2.scatter(loadings[:, 0], loadings[:, 1], c='gray', alpha=0.5, s=30)

# # For each top loading, determine if it contributes more to AR or CC
# # We'll use the sign of Component 1 since samples are separated along this axis
# for idx in top_indices:
#     color = COLORS['AR'] if loadings[idx, 0] > 0 else COLORS['CC']
#     ax2.scatter(loadings[idx, 0], loadings[idx, 1], c=color, s=50)
#     ax2.annotate(metabolites.iloc[idx], 
#                  (loadings[idx, 0], loadings[idx, 1]),
#                  xytext=(5, 5), textcoords='offset points',
#                  fontsize=8,
#                  bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))

# ax2.set_xlabel('Component 1')
# ax2.set_ylabel('Component 2')
# ax2.set_title('PLS-DA Loadings Plot')
# ax2.grid(True, alpha=0.3)

# # Add origin lines
# ax2.axhline(y=0, color='k', linestyle='-', alpha=0.3)
# ax2.axvline(x=0, color='k', linestyle='-', alpha=0.3)

# plt.tight_layout()
# plt.show()

# # Save plot
# fig.savefig(pjoin(OUTPUT_FOLDER, 'plsda_scores_and_loadings_ar_cc_batch_3.png'), 
#             dpi=600, bbox_inches='tight')


"""
Run T-test for AR vs CC and Generate Volcano Plot (-log10(raw p) vs log2(FC)) 
"""
# Run t-test for AR vs CC
ar_cc_groups = {k: v for k, v in sample_groups.items() if k in ['AR', 'CC']}
data_volcano = data.copy()
# Keep only the following column names: "Metabolite", "Kegg ID", "Metabolite Class", "Confidence"  and the columns in ar_cc_groups 
data_volcano = data_volcano[[CMPD_COL_NAME, 'Kegg ID', 'Metabolite Class', 'Confidence'] + ar_cc_groups['AR'] + ar_cc_groups['CC']]

# Get minimum non-zero data value in the data_volcano columns for ar_cc_groups. Use this divided by 5 to substitute all 0 values in the ar_cc_groups columns of data_volcano
data_volcano_find_min_will_delete = data_volcano.copy()
# Replace all 0 values with 0.1. We will ignore these. After substituting 0 values with the arbitrary value, find the non-zero minimum data value across the ar_cc_groups data columns.
data_volcano_find_min_will_delete[ar_cc_groups['AR'] + ar_cc_groups['CC']] = data_volcano_find_min_will_delete[ar_cc_groups['AR'] + ar_cc_groups['CC']].replace(0, 0.1)
# Find the minimum non-zero value across the ar_cc_groups data columns
min_non_zero_val = data_volcano_find_min_will_delete[ar_cc_groups['AR'] + ar_cc_groups['CC']].min().min()
# Delete data_volcano_find_min_df because we arbitrarily alterred the data and do not need the df anymore
del data_volcano_find_min_will_delete

# For data_volcano, replace all 0 values with min_non_zero_val/5
data_volcano[ar_cc_groups['AR'] + ar_cc_groups['CC']] = data_volcano[ar_cc_groups['AR'] + ar_cc_groups['CC']].replace(0, min_non_zero_val/5)

# Filter out data with nearly identical values across all samples (AR and CC). This is done by calculating the coefficient of variation (CV) for each metabolite across all samples. If the CV is less than VOLCANO_PLOT_CV_CUTOFF, the metabolite is considered to have nearly identical values and is filtered out of t-test analysis and volcano plot generation
# Calculate CV for each metabolite across all samples
# Add a cv column to data_volcano
data_volcano['cv'] = data_volcano[ar_cc_groups['AR'] + ar_cc_groups['CC']].std(axis=1) / data_volcano[ar_cc_groups['AR'] + ar_cc_groups['CC']].mean(axis=1)

# Export data_volcano to a csv file
data_volcano.to_csv(pjoin(OUTPUT_FOLDER, 'data_volcano_pre_cv_filter.csv'), index=False)

# Filter out rows with CV less VOLCANO_PLOT_CV_CUTOFF
data_volcano = data_volcano[data_volcano['cv'] >= VOLCANO_PLOT_CV_CUTOFF]

# Perform t-test for each metabolite
t_test_results = []
for metabolite in data_volcano[CMPD_COL_NAME]:
    # Get values for AR and CC
    ar_values = data_volcano.loc[data_volcano[CMPD_COL_NAME] == metabolite, ar_col_names].values[0]
    cc_values = data_volcano.loc[data_volcano[CMPD_COL_NAME] == metabolite, cc_col_names].values[0]

    # Perform t-test
    t_stat, p_val = stats.ttest_ind(ar_values, cc_values)
    
    # Calculate fold change (FC)
    fc = np.mean(ar_values) / np.mean(cc_values)
    
    # Store results
    t_test_results.append({
        'Metabolite': metabolite,
        't_stat': t_stat,
        'p_val': p_val,
        'FC': fc
    })

# Convert to DataFrame
t_test_df = pd.DataFrame(t_test_results)

# Calculate -log10(p-value) and log2(FC) column values for volcano plot
t_test_df['-log10(p_val)'] = -np.log10(t_test_df['p_val'])
t_test_df['log2(FC)'] = np.log2(t_test_df['FC'])

# Use add_to_df to add t_test_df columns to data_volcano (add columns t_stat, p_val, FC, -log10(p_val), log2(FC) to data_volcano). Use Metabolite as the key column
data_volcano = add_to_df(data_volcano, t_test_df, ['t_stat', 'p_val', 'FC', '-log10(p_val)', 'log2(FC)'], CMPD_COL_NAME)

# Create volcano plot
fig, ax = plt.subplots(figsize=(12, 8))

# Add dotted lines for cutoffs
p_cutoff = -np.log10(0.05)  # -log10(0.05) for p-value cutoff
fc_cutoff = 1  # log2(FC) cutoff of 1 (2-fold change)

# Add horizontal dotted line for p-value cutoff
ax.axhline(y=p_cutoff, color='black', linestyle='--', alpha=0.5, linewidth=LINE_WIDTH)
# Add vertical dotted lines for FC cutoffs
ax.axvline(x=fc_cutoff, color='black', linestyle='--', alpha=0.5, linewidth=LINE_WIDTH)
ax.axvline(x=-fc_cutoff, color='black', linestyle='--', alpha=0.5, linewidth=LINE_WIDTH)

# Plot all metabolites
ax.scatter(data_volcano['log2(FC)'], data_volcano['-log10(p_val)'], color='gray', alpha=0.5, s=FONT_SIZE, label='Not significant')
# Highlight significant metabolites (p < 0.05). For log2(FC) > 1, use the color for AR, for log2(FC) < -1, use the color for CC
sig_mask = data_volcano['p_val'] < 0.05
sig_ar = data_volcano[sig_mask & (data_volcano['log2(FC)'] > 1)]
sig_cc = data_volcano[sig_mask & (data_volcano['log2(FC)'] < -1)]

# Plot significant points separately for legend
ax.scatter(sig_ar['log2(FC)'], sig_ar['-log10(p_val)'], color=COLORS['AR'], alpha=0.7, s=30, 
          label=f'Enriched in {FULL_NAMES["AR"]}')
ax.scatter(sig_cc['log2(FC)'], sig_cc['-log10(p_val)'], color=COLORS['CC'], alpha=0.7, s=30,
          label=f'Enriched in {FULL_NAMES["CC"]}')

# Add data point labels for metabolites with extreme fold changes (|log2(FC)| > 3) that don't contain "Unknown" and with p_val < 0.05
extreme_metabolites = data_volcano[
    (~data_volcano[CMPD_COL_NAME].str.contains('Unknown', na=False)) & 
    (data_volcano['p_val'] < 0.05) &
    (data_volcano['log2(FC)'].abs() > 2)
]

# Use adjustText library to prevent label overlap
texts = []
for _, row in extreme_metabolites.iterrows():
    texts.append(ax.text(row['log2(FC)'], row['-log10(p_val)'], row[CMPD_COL_NAME],
                        fontsize=14, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')))

# Adjust text positions to prevent overlap.
adjust_text(texts, 
           arrowprops=dict(arrowstyle='-', color='black', lw=1, alpha=0.5),
           expand_points=(1.5, 1.5),
           force_points=(0.1, 0.1),
           avoid_self=True)

# Set linewidth for axes and ticks and increase tick length
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(LINE_WIDTH)
ax.tick_params(width=2, length=8, labelsize=FONT_SIZE)

# Add labels
ax.set_xlabel('log2(FC)', fontsize=FONT_SIZE)
ax.set_ylabel('-log10(p-value)', fontsize=FONT_SIZE)

# Add legend with italicized species names and increased font size
legend = ax.legend(prop={'size': 16}) 
for text in legend.get_texts():
    for species, full_name in FULL_NAMES.items():
        if full_name in text.get_text() and ITALICIZE_NAMES[species]:
            text.set_text(text.get_text().replace(full_name, r'$\mathit{' + full_name + '}$'))

# Save plot
fig.savefig(pjoin(OUTPUT_FOLDER, 'volcano_plot_NIST_matches_incl_unknowns_CV_filtered_ar_cc_batch_3.png'), dpi=600, bbox_inches='tight')

# Export data for volcano plot
data_volcano.to_csv(pjoin(OUTPUT_FOLDER, 'volcano_plot_data_NIST_matches_incl_unknowns_CV_filtered_ar_cc_batch_3.csv'), index=False)


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

# ANOVA analysis --- commented out:
# anova_results = analyze_metabolites(data_knowns, sample_groups_no_blank)

# # Save ANOVA results
# anova_results.to_csv(pjoin(OUTPUT_FOLDER, 'anova_results_batch_3.csv'), index=False)


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
# create a grid of barplots in a 3 by 5 grid. The 15th plot will be the legend.
barplot_grid_order = ["Maltotriitol", "Maltotriose", "Melibiose", "2,3-dihydroxyisovaleric acid", "D-malic acid", "Fumaric acid", "Glyceric acid", "Succinic acid", "D-mannose", "Lactulose", "DL-dihydrosphingosine", "3-(4-hydroxyphenyl)propionic acid (phloretic acid)", "4-hydroxy-3-methoxybenzoic acid (isovanillic acid)", "4-hydroxybenzoic acid (p-salicylic acid)"]

# Create and save the grid figure
grid_fig, missing = create_barplot_grid(data_knowns, barplot_grid_order, sample_groups, COLORS, FULL_NAMES, ITALICIZE_NAMES)
grid_fig.savefig(pjoin(OUTPUT_FOLDER, 'metabolite_barplots_grid_batch_3.png'), 
                 dpi=600, bbox_inches='tight')

# Print any missing metabolites
if missing:
    print("\nMetabolites not found in grid dataset:")
    for met in missing:
        print(f"- {met}")