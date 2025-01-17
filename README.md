# GC-MS Metabolic Profiling of Anaerobic Gut Fungi

## Background
This repository contains scripts for analyzing GC-MS metabolomics data from anaerobic gut fungi, specifically comparing metabolic profiles between A. robustus and C. churrovis species. The analysis pipeline processes MS-DIAL outputs, performs statistical comparisons, and generates visualizations.

## Associated Publication
[to be filled in]

## Setup

### Dependencies
- Python 3.x
- Core packages:
    - os
    - pandas
    - numpy 
    - matplotlib
    - seaborn
    - scipy
    - scikit-learn
- Statistics:
    - ppca
    - statsmodels
- File handling:
    - xlsxwriter
- Network visualization:
    - py4cytoscape (requires Cytoscape desktop application)

All dependencies can be installed via pip

### Additional Programs and Tools
- MS-DIAL
- GNPS
- Cytoscape
- BLAST
- Annotation Source (ie: JGI Mycocosm Portal)

### Required Input Files
- NIST analysis
- MS-DIAL outputs (peak areas and TIC-normalized data)
- GNPS library matches
- Cytoscape network files (.graphml)
- Proteomics data
- FASTA sequence files
- KOG annotation files

## Description of Scripts
| Script Name                                       | Description                               |
| ------------------------------------------------- | ----------------------------------------- |
| Script_1_MSDIAL_Statistics.py |Script 1 processes raw MS-DIAL outputs and TIC-normalized peak areas, calculates summary statistics across sample groups, generates statistical comparisons (t-tests, FDR-adjusted p-values), and creates a formatted output for input to MetaboAnalyst|
| Script_2_Assemble_Summary_Data_Excel.py | Script 2 combines data from MS-DIAL and GNPS analyses to create filtered tables of significant metabolites. The script also generates volcano plots, histograms, and a formatted metabolic network in Cytoscape. Note, prior to running this script, you must open the Cytoscape program|
| Script_3_NIST14_Fatty_Acid_Profiling.py | Script 3 takes manually curated fatty acid composition data from GC-MS analysis and calculates relative compositions between samples. The script generates scatter plots comparing fatty acid compositions and highlights significant differences |
| Script_4_PCA_and_Heatmap_Plots.py | Script 4 performs principal component analysis (PCA), probabilistic PCA (pPCA), and generates hierarchically clustered heatmaps to visualize metabolite abundances across sample groups. The script also creates individual metabolite bar plots |
| Script_5_BLASTp_and_Annotations_for_Proteomics_Results.py | Script 5 processes proteomics data by aligning protein sequences using BLASTp and combining results with KOG functional annotations. The script is run in two parts, with manual BLASTp execution required between parts |

## Summary of Output Files

### Script 1 (MSDIAL_Statistics.py)
Temp folder:
- `MSDIAL_stats.xlsx`: summary of statistics for MS-DIAL analysis
- `GF_GCMS_MetaboAnalyst_input.csv`: Formatted data for MetaboAnalyst analysis

### Script 2 (Assemble_Summary_Data_Excel.py)
Output folder:
- `GF_GCMS_stats_summary_table.xlsx`: Combined metabolite data with statistical analyses 
- `histogram_*_log10_avg_intensity.png`: Distribution plots showing log10 peak intensities across sample types
- `volcano_plot_*.png`: Volcano plots comparing different sample groups
- `GF_GCMS_cytoscape_batch_3.cys`: Cytoscape session file with metabolic network

### Script 3 (NIST14_Fatty_Acid_Profiling.py)
Output folder:
- `Fatty_Acid_Compositions_*.png`: boxplot depicting fatty acid profile for input gut fungi based on NIST compound identifications

### Script 4 (PCA_and_Heatmap_Plots.py)
Output folder:
- `pca_plot_*.png`, `ppca_plot_*.png`: PCA analysis visualizations
- `metabolite_heatmap_*.png`: Hierarchically clustered heatmaps 
- `anova_results_batch_3.csv`: Statistical comparisons between groups
- `/Metabolite Class Heatmaps/`: Folder containing class-specific heatmaps
- `/Metabolite Bar Charts/`: Folder containing individual metabolite plots

### Script 5 (BLASTp_and_Annotations_for_Proteomics_Results.py)
Temp folder:
- `*_query_BLASTp.txt`: Query sequences for BLASTp
- `*_db_BLASTp.fasta`: Database sequences for BLASTp
Output folder:
- `*_BLASTp_proteomics_results_final.xlsx`: Proteomics data with BLASTp alignments and JGI Mycocosm KOG annotations

## Support
For support with using these scripts, please contact lbutkovich@ucsb.edu.

## Authors and Acknowledgements
Primary author: Lazarina Butkovich (University of California, Santa Barbara)

Thank you to Fred Krauss for feedback and assistance in writing and formatting these scripts. 
