"""
GF GCMS Data Analysis Script 3: Fatty Acid Profiling for Batches 1 and 3
Lazarina Butkovich 10/14/24

Take Chaevien's GCMS data and analyze the fatty acid profiling
"""

import pandas as pd
from os.path import join as pjoin
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

INPUT_FOLDER = r'input' 
TEMP_FOLDER = r'temp'
OUTPUT_FOLDER = r'output'


"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
FA_DATA_FILENAME = 'Fatty_Acids_Chaevien_20241014.xlsx'

# 'Batch 3' (my batch) or 'Batch 1'
FA_DATA_SHEET_NAME = 'Batch 1'

if FA_DATA_SHEET_NAME == 'Batch 3':
# Batch 3: AR, CC, MC, RF
    SAMPLE_GROUP_NAMES = ['AR', 'CC']
elif FA_DATA_SHEET_NAME == 'Batch 1':
    # Batch 1: AR, CC, G1, S3, PF
    SAMPLE_GROUP_NAMES = ['AR', 'CC', 'G1', 'S3', 'PF']
else: 
    raise ValueError('Indicate FA_DATA_SHEET_NAME. Must be "Batch 1" or "Batch 3"')

REP_NUM = 4

# colors
AR_COLOR = '#00008B'
CC_COLOR = '#00FF00'
G1_COLOR = 'lightblue'
S3_COLOR = 'indianred'
PF_COLOR = 'khaki'

# Label names
AR_LABEL = 'A. robustus'
CC_LABEL = 'C. churrovis'
G1_LABEL = 'N. californiae'
S3_LABEL = 'N. lanati'
PF_LABEL = 'P. finnis'


"""""""""""""""""""""""""""""""""""""""""""""
Import
"""""""""""""""""""""""""""""""""""""""""""""
# Import the data
fa_data = pd.read_excel(pjoin(INPUT_FOLDER, FA_DATA_FILENAME), sheet_name=FA_DATA_SHEET_NAME)


"""""""""""""""""""""""""""""""""""""""""""""
Analysis
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
color_dict = {'AR': AR_COLOR, 'CC': CC_COLOR, 'G1': G1_COLOR, 'S3': S3_COLOR, 'PF': PF_COLOR}
label_dict = {'AR': AR_LABEL, 'CC': CC_LABEL, 'G1': G1_LABEL, 'S3': S3_LABEL, 'PF': PF_LABEL}

for sample_group in SAMPLE_GROUP_NAMES:
        plt.scatter(cmpd_list, dict_avg[sample_group], color=color_dict[sample_group], label=label_dict[sample_group], s=5)
        plt.errorbar(cmpd_list, dict_avg[sample_group], yerr=dict_std[sample_group], fmt='o', color=color_dict[sample_group], capsize=5, markersize=5)

# For the data points that are significant (no overlap of the error bars), label the datapoints with the average value (rounded to 1 decimal place with a % sign). Place the labels such that they do not overlap anything else on the plot. Make the font size smaller.
# If batch 3:
if FA_DATA_SHEET_NAME == 'Batch 3':
    for i in range(len(cmpd_list)):
        if abs(dict_avg['AR'][i] - dict_avg['CC'][i]) > dict_std['AR'][i] + dict_std['CC'][i]:
            plt.text(cmpd_list[i], dict_avg['AR'][i], '  ' + str(round(dict_avg['AR'][i], 1)) + '%', ha='left', va='bottom', fontsize=8)
            plt.text(cmpd_list[i], dict_avg['CC'][i], '  ' + str(round(dict_avg['CC'][i], 1)) + '%', ha='left', va='top', fontsize=8)

# Legend labels: italicized LABEL_NAME values. Increase legend size and the dot size in the legend.
font_properties = FontProperties()
font_properties.set_style('italic')
font_properties.set_size('large')
legend = plt.legend(prop=font_properties, markerscale=2)
for text in legend.get_texts():
    text.set_fontproperties(font_properties)

# label y-axis '% Fatty Acid Composition'
plt.ylabel('% Fatty Acid Composition')

# Rotate the x-axis labels to be slanted slightly, read left to right
plt.xticks(rotation=-30, ha='left')

# Make the plot fit the figure
plt.tight_layout()

# Export as png (dpi=600) to output folder
plt.savefig(pjoin(OUTPUT_FOLDER, 'Fatty_Acid_Compositions_'+ FA_DATA_SHEET_NAME + '.png'), dpi=600)

plt.show()