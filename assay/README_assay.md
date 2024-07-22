# README for bitter assay  

## Description  
Part 1 is a Jupyter notebook designed to create plots from raw assay data.  
Part 2 is an R script designed to create dose response curves  

## How to install and run  
Download script and open in a Jupyter environment or R, as appropriate.  

## How to use  
### Create heatmap for assay data  
Open Aki_analysis_v2.ipynb in a Jupyter environment  
  Ensure all necessary packages are available (see first cell) and up to date.  
  Update location of raw_file and converter_file.  The expected column names in raw_file are: well, compound, concentration, gene, group, AUC, file.  The only essential columns are gene & AUC.  
  We have included our version of converter_file (Converter_file.csv) which was used to convert between different gene nomenclatures used during early experiments.  If this is not necessary for your purposes, can comment this section out of code.  
  We have included two different functions for creating heatmaps.  The second (create_colored_max_p_values_byrow) is what was used in the paper.  In this function, data are sorted, raw p-values are calculated, and then p-values are adjusted with a BH correction within each gene.  Final color is based on fold change relative to buffer control with thresholds applied based on adjusted p-value.  The code used to call this function is shown in the subsequent cells.  We have also included an alternate function (create_colored_max_p_values_byall) which does the same thing except that p-value adjustment occurs with the entire dataset.  This function was created for purposes of comparison, but when the results were found to be similar, we decided to use create_colored_max_p_values_byrow going forward.    
  At the end, function plot_all_chems creates bar plots for all data for a given gene.  Several examples of the expected syntax are included.  This function was not used in the final manuscript.    
  
### Create dose response curves for selected compounds  
Open DR_curve.R in an R environment such as RStudio  
  Ensure all necessary packages are available and up to date.  
  Line 8: set working directory where plots will be made  
  Line 11: point towards raw data file  
  Lines 15-23: used to convert between different gene nomenclatures used during early experiments  

## Credits  
Underlying bitter assay developed by Y. Toda and Y. Ishimaru and implemented by A. Itoigawa.  
Code developed by K. Higgins and A. Itoigawa (authors) with programming advice from Allen B. Davis (acknowledged).  ChatGPT v3 and 4 used to draft early versions of this script.  
