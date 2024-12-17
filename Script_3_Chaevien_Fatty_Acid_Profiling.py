"""
GF GCMS Data Analysis Script 3, Option 3
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
Functions
"""""""""""""""""""""""""""""""""""""""""""""


"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
FA_DATA_FILENAME = 'Fatty_Acids_Chaevien_20241014.xlsx'

SAMPLE_GROUP_NAMES = ['AR', 'CC', 'MC', 'RF']

REP_NUM = 4


"""""""""""""""""""""""""""""""""""""""""""""
Import
"""""""""""""""""""""""""""""""""""""""""""""
# Import the data
fa_data = pd.read_excel(pjoin(INPUT_FOLDER, FA_DATA_FILENAME))


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


# Plot the average fatty acid compositions for sample groups CC and S4. Make the plot a scatter plot with the x-axis as the fatty acid name and the y-axis as the average fatty acid composition. The error bars should be the standard deviation of the fatty acid compositions for each sample in the group. The x-axis can be ordered in descending average fatty acid composition for CC.
# Sort the fatty acids by the average composition in AR
fa_data = fa_data.sort_values('CC_avg_composition', ascending=False)

# Plot the average fatty acid compositions for sample groups AR and CC
cc_avg_list= fa_data['CC_avg_composition']
cc_std_list = fa_data['CC_std_composition']
ar_avg_list = fa_data['AR_avg_composition']
ar_std_list = fa_data['AR_std_composition']
cmpd_list = fa_data['Compound Name']

# Color AR as blue, and CC as green
# Code for lime green: #00FF00
# Code for dark blue: #00008B
# Make the data points smaller
plt.scatter(cmpd_list, ar_avg_list, color='#00008B', label='A. robustus', s=5)
plt.errorbar(cmpd_list, ar_avg_list, yerr=ar_std_list, fmt='o', color='#00008B', capsize=5, markersize=5)
plt.scatter(cmpd_list, cc_avg_list, color='#00FF00', label='C. churrovis', s=5)
plt.errorbar(cmpd_list, cc_avg_list, yerr=cc_std_list, fmt='o', color='#00FF00', capsize=5, markersize=5)

# For the data points that are significant (no overlap of the error bars), label the datapoints with the average value (rounded to 1 decimal place with a % sign). Place the labels such that they do not overlap anything else on the plot. Make the font size smaller.
for i in range(len(cmpd_list)):
    if abs(ar_avg_list[i] - cc_avg_list[i]) > ar_std_list[i] + cc_std_list[i]:
        plt.text(cmpd_list[i], ar_avg_list[i], '  ' + str(round(ar_avg_list[i], 1)) + '%', ha='left', va='bottom', fontsize=8)
        plt.text(cmpd_list[i], cc_avg_list[i], '  ' + str(round(cc_avg_list[i], 1)) + '%', ha='left', va='top', fontsize=8)

# Legend labels: italicized "A. robustus" and "C. churrovis". Increase legend size and the dot size in the legend.
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
plt.savefig(pjoin(OUTPUT_FOLDER, 'Fatty_Acid_Compositions.png'), dpi=600)

plt.show()