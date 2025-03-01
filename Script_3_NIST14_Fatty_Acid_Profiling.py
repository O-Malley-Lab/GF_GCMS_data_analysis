"""
Gut Fungal GC-MS Profiling, Script 3: Fatty Acid Profiling from NIST14 Matches
Lazarina Butkovich 2024

This script takes manually curated fatty acid composition data from GC-MS analysis and does the following:
1. Imports fatty acid data from Excel sheet
2. Calculates relative fatty acid compositions for each sample
3. Calculates averages and standard deviations across replicates
4. Creates a scatter plot comparing fatty acid compositions between samples
5. Labels significantly different data points
6. Exports plot as high-resolution PNG
"""


import pandas as pd
from os.path import join as pjoin
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'

# Input data from NIST analysis
FA_DATA_FILENAME = 'Fatty_Acids_Chaevien_20241014.xlsx'
SAMPLE_GROUP_NAMES = ['AR', 'CC']
FA_DATA_SHEET_NAME = 'Batch 3'
REP_NUM = 4

# Sample group colors for plotting
COLORS = {'CC':'lightgreen', 'AR':'darkblue', 'MC':'lightgrey', 'RF':'dimgrey', 'FAMES':'pink', 'BLANK':'olive'}

# Label names
AR_LABEL = 'A. robustus'
CC_LABEL = 'C. churrovis'


"""""""""""""""""""""""""""""""""""""""""""""
Import data
"""""""""""""""""""""""""""""""""""""""""""""
# Import the data
fa_data = pd.read_excel(pjoin(INPUT_FOLDER, FA_DATA_FILENAME), sheet_name=FA_DATA_SHEET_NAME)
FONT_SIZE = 20
LINE_WIDTH = 2


"""""""""""""""""""""""""""""""""""""""""""""
Analysis for Fatty Acid Profiling
"""""""""""""""""""""""""""""""""""""""""""""
# For each sample replicate, sum the fatty acid values. 
# Sample name format: 'AR_1' etc.
# Store the sums in a dictionary with the sample name as the key.
fa_sums = {}
for sample_group in SAMPLE_GROUP_NAMES:
    for rep in range(1, REP_NUM+1):
        sample_name = sample_group + '_' + str(rep)
        fa_sums[sample_name] = fa_data[sample_name].sum()

# For each sample, generate a new column (ie: AR_1_composition) that is the percentage of each fatty acid in the total fatty acid sum for a sample_name
for sample_group in SAMPLE_GROUP_NAMES:
    for rep in range(1, REP_NUM+1):
        sample_name = sample_group + '_' + str(rep)
        fa_data[sample_name + '_composition'] = fa_data[sample_name] / fa_sums[sample_name]

# For each sample group, generate a column (ie: AR_avg_composition) that is the average of the fatty acid compositions for each sample in the group. Also generate a column (ie: AR_std_composition) that is the standard deviation of the fatty acid compositions for each sample in the group. Make the values percentages.
for sample_group in SAMPLE_GROUP_NAMES:
    sample_group_compositions = [sample_group + '_' + str(rep) + '_composition' for rep in range(1, REP_NUM+1)]
    fa_data[sample_group + '_avg_composition'] = fa_data[sample_group_compositions].mean(axis=1)*100
    fa_data[sample_group + '_std_composition'] = fa_data[sample_group_compositions].std(axis=1) * 100


"""""""""""""""""""""""""""""""""""""""""""""
Plot Fatty Acid Profiling
"""""""""""""""""""""""""""""""""""""""""""""
# Plot the average fatty acid compositions for each sample group . Make the plot a scatter plot with the x-axis as the fatty acid name and the y-axis as the average fatty acid composition. The error bars should be the standard deviation of the fatty acid compositions for each sample in the group. The x-axis can be ordered in descending average fatty acid composition for CC.
# Sort the fatty acids by the average composition in the 1st sample group (ie: AR) in descending order
fa_data = fa_data.sort_values(by=SAMPLE_GROUP_NAMES[0] + '_avg_composition', ascending=False)

# Plot the average fatty acid compositions for each sample group
# Create a dictionary of lists
dict_avg = {}
dict_std = {}
for sample_group in SAMPLE_GROUP_NAMES:
    dict_avg[sample_group] = fa_data[sample_group + '_avg_composition']
    dict_std[sample_group] = fa_data[sample_group + '_std_composition']

cmpd_list = fa_data['Compound Name']
label_dict = {'AR': AR_LABEL, 'CC': CC_LABEL}

# Set global font sizes and line widths
plt.rcParams.update({
    'font.size': FONT_SIZE,
    'axes.linewidth': LINE_WIDTH,
    'lines.linewidth': LINE_WIDTH,
    'axes.labelsize': FONT_SIZE,
    'xtick.major.width': LINE_WIDTH,
    'ytick.major.width': LINE_WIDTH,
    'xtick.labelsize': FONT_SIZE,
    'ytick.labelsize': FONT_SIZE
})

# Create figure with larger size
plt.figure(figsize=(12, 8))

# Plot the data points and error bars with increased line widths
for sample_group in SAMPLE_GROUP_NAMES:
    plt.scatter(cmpd_list, dict_avg[sample_group], 
               color=COLORS[sample_group], 
               label=label_dict[sample_group], 
               s=75)  # Increased marker size
    plt.errorbar(cmpd_list, dict_avg[sample_group], 
                yerr=dict_std[sample_group], 
                fmt='o', 
                color=COLORS[sample_group], 
                capsize=5, 
                markersize=5,
                capthick=3,
                elinewidth=3)

# For the data points that are significant (no overlap of the error bars), label the datapoints with the average value (rounded to 1 decimal place with a % sign). Place the labels such that they do not overlap anything else on the plot.
for i in range(len(cmpd_list)):
    if abs(dict_avg['AR'][i] - dict_avg['CC'][i]) > dict_std['AR'][i] + dict_std['CC'][i]:
        plt.text(cmpd_list[i], dict_avg['AR'][i], 
                '  ' + str(round(dict_avg['AR'][i], 1)) + '%', 
                ha='left', va='bottom', 
                fontsize=18)
        plt.text(cmpd_list[i], dict_avg['CC'][i], 
                '  ' + str(round(dict_avg['CC'][i], 1)) + '%', 
                ha='left', va='top', 
                fontsize=18)

# Update legend properties
font_properties = FontProperties()
font_properties.set_style('italic')
font_properties.set_size('large')
legend = plt.legend(prop=font_properties, markerscale=2)
for text in legend.get_texts():
    text.set_fontproperties(font_properties)

# Update axis labels with larger font size
plt.ylabel('% Fatty Acid Composition', fontsize=24)

# Update tick parameters
plt.xticks(rotation=-30, ha='left')
plt.tick_params(axis='both', which='major', labelsize=FONT_SIZE, width=3, length=15)

# Make the plot fit the figure
plt.tight_layout()

# Export as png (dpi=600) to output folder
plt.savefig(pjoin(OUTPUT_FOLDER, 'Fatty_Acid_Compositions_'+ FA_DATA_SHEET_NAME + '.png'), dpi=600)

plt.show()