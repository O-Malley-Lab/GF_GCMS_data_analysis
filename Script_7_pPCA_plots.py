"""
GF GCMS Data Analysis Script 7: Batch 1 pPCA Plots
Lazarina Butkovich 12/20/24
"""
import pandas as pd
import numpy as np
from os.path import join as pjoin
import matplotlib.pyplot as plt
# Use probabilistic PCA
from ppca import PPCA

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
    
    # Debug print
    print(f"Data shape: {X.shape}")
    print(f"Number of labels: {len(sample_labels)}")
    
    # Fit PPCA
    ppca = PPCA()
    ppca.fit(data=X, d=2)
    transformed = ppca.transform()
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Plot each sample group
    for sample_type in np.unique(sample_labels):
        mask = (sample_labels == sample_type)
        plt.scatter(transformed[mask, 0], 
                   transformed[mask, 1],
                   c=colors[sample_type],
                   label=sample_type,
                   alpha=0.7)

    # Update axis labels with variance explained
    plt.xlabel('Scores on PC1 ()')
    plt.ylabel('Scores on PC2 ()')
    
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    return plt.gcf()

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

COLORS = {
    'AR': 'darkblue',
    'CC': 'lightgreen', 
    'G1': 'lightblue',
    'S3': 'peru',
    'PF': 'khaki',
    'BLANK': 'olive'
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