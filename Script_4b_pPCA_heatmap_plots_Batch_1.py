"""
GF GCMS Data Analysis Script 4b: pPCA and Heatmap Plots for Batch 1
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

"""
Functions
"""
def create_combined_ppca_plot(data, sample_groups, colors, title="Combined pPCA Analysis"):
    # Prepare combined data
    all_cols = []
    sample_labels = []
    for sample_type, cols in sample_groups.items():
        all_cols.extend(cols)
        sample_labels.extend([sample_type] * len(cols))

    # Convert to numpy array for proper indexing 
    sample_labels = np.array(sample_labels)
    
    # Extract and transpose data for PPCA (samples as rows)
    X = data[all_cols].T.values
    
    # Fit PPCA
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

def create_metabolite_heatmap(anova_results, q_value_threshold=0.05):
    """Create heatmap of significant metabolites using FDR q-values with clustered y-axis"""
    # Filter significant metabolites
    significant = anova_results[anova_results['q_value'] < q_value_threshold].copy()
    significant['q_value'] = -np.log10(significant['q_value'])
    
    # Prepare data for heatmap
    heatmap_data = significant.set_index('Metabolite').drop(['p_value', 'q_value'], axis=1)
    
    # Rename columns with full species names
    heatmap_data = heatmap_data.rename(columns=FULL_NAMES)
    
    # Scale the data
    scaled_data = (heatmap_data - heatmap_data.mean()) / heatmap_data.std()
    
    # Perform hierarchical clustering on metabolites
    row_linkage = hierarchy.linkage(scaled_data, method='ward')
    
    # Create clustered heatmap
    g = sns.clustermap(scaled_data,
                   row_linkage=row_linkage,
                   col_cluster=False,
                   cmap='RdBu_r',
                   center=0,
                   xticklabels=True,
                   yticklabels=True,
                   cbar_kws={'label': 'Z-score', 'orientation': 'horizontal'},
                   cbar_pos=(0.15, 0.95, 0.4, 0.02),
                   figsize=(12, len(significant) * 0.3))
    
    # Set title
    g.ax_heatmap.set_title(f'Differentially Abundant Metabolites\n(FDR q < {q_value_threshold})', 
                          pad=20)
    
    # Italicize species names on x-axis
    plt.setp(g.ax_heatmap.get_xticklabels(), style='italic')
    
    # Rotate labels for better readability
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    
    return g.fig


"""
Values
"""
AR_COL_NAMES = ['Ar_01','Ar_02', 'Ar_03', 'Ar_04']
CC_COL_NAMES = ['Cc_01', 'Cc_02', 'Cc_03', 'Cc_04']
G1_COL_NAMES = ['Nc_G1_01', 'Nc_G1_02', 'Nc_G1_03', 'Nc_G1_04']
S3_COL_NAMES = ['NS3_01', 'NS3_02', 'NS3_03', 'NS3_04']
PF_COL_NAMES = ['Pf_01', 'Pf_02', 'Pf_03', 'Pf_04']
BLANK_COL_NAMES = ['Blank 01', 'Blank 02', 'Blank 03']
CMPD_COL_NAME = 'Metabolite'
SCORE_COL_NAME = 'Score'

DATA_FILENAME = 'Chaevien_batch_1.xlsx'

INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'

# Define sample groups and colors
SAMPLE_GROUPS = {
    'AR': AR_COL_NAMES,
    'CC': CC_COL_NAMES,
    'G1': G1_COL_NAMES,
    'S3': S3_COL_NAMES,
    'PF': PF_COL_NAMES,
    'BLANK': BLANK_COL_NAMES
}

# Sample colors
COLORS = {'CC':'lightgreen', 'AR':'darkblue', 'MC':'tan', 'RF':'brown', 'FAMES':'pink', 'BLANK':'olive', 'G1': 'lightblue', 'S3': 'peru', 'PF': 'khaki'}

FULL_NAMES = {
    'AR': 'A. robustus',
    'CC': 'C. churrovis',
    'G1': 'N. californiae',
    'S3': 'N. lanati',
    'PF': 'P. finnis',
    'BLANK': 'Blank'
}

"""
Import data
"""
data = pd.read_excel(pjoin(INPUT_FOLDER, DATA_FILENAME), sheet_name='Data')

"""
Generate pPCA plot
"""
# Generate combined pPCA plot
fig = create_combined_ppca_plot(data, SAMPLE_GROUPS, COLORS)
plt.show()
# Save plot
fig.savefig(pjoin(OUTPUT_FOLDER, 'ppca_plot_batch_1.png'), dpi=600, bbox_inches='tight')

"""
Generate ANOVA analysis and heatmap
"""
# Create data_knowns, the data df filtered to remove rows with 'Unknown...' in the CMPD_COL_NAME column
data_knowns = data[~data[CMPD_COL_NAME].str.contains('Unknown')]

# For rows with the same metabolite, remove rows with the lower Score values.
data_knowns = data_knowns.sort_values(by=[CMPD_COL_NAME, SCORE_COL_NAME], ascending=False).drop_duplicates(subset=CMPD_COL_NAME)

# Perform ANOVA analysis
anova_results = analyze_metabolites(data_knowns, SAMPLE_GROUPS)


heatmap_fig = create_metabolite_heatmap(anova_results)

# Show plot first
plt.show()

# Clean filename and save
output_file = os.path.join(OUTPUT_FOLDER, 'metabolite_heatmap_batch_1.png')
heatmap_fig.savefig(output_file, dpi=600, bbox_inches='tight')


# Save ANOVA results to Excel
anova_results.sort_values('p_value').to_excel(
    pjoin(OUTPUT_FOLDER, 'anova_results_batch_1.xlsx'),
    index=False
)
