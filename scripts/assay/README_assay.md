# README for bitter assay  

## Description  
Scripts necessary to analyze and create figures for assay data, as seen in Figure 6.

## How to install and run  
Download scripts and file tree structure as previously described.  Run scripts in python or R as appropriate.

## How to use
### 1. Create heatmap for assay data  
  Ensure data are in expected format, as shown with example data in results/assay/data_BitterAssay_ver240513.csv.  If necessary, modify the converter file located in the same folder.
  Ensure all necessary packages are available and up to date.  
  Additional functions are included at end of code (commented out) to perform slightly different analyses.  These do not appear in the manuscript.

### 2. Create dose response curves for selected compounds  
Open Lineplot_dose_response.R in an R environment such as RStudio  
  Ensure all necessary packages are available and up to date (tidyverse and scales).  
  The datasheet format is the same as for the datasheet for heatmap. The essential columns to crate dose-response curves are "gene", "coumpound", "concentration", and "AUC".

### 3.  barplots for agonist screenings from raw assay data (in Supplementary materials).
Open Barplot_agonist_screening.R in an R environment such as RStudio  
  Ensure all necessary packages are available and up to date (tidyverse).   
  The datasheet format is the same as for the datasheet for heatmap. The essential columns to crate dose-response curves are "gene", "coumpound", and "AUC".


## Credits  
Underlying bitter assay developed by Y. Toda and Y. Ishimaru and implemented by A. Itoigawa.  
Code developed by K. Higgins and A. Itoigawa (authors) with programming advice from Allen B. Davis (acknowledged).  ChatGPT v3 and 4 used to draft early versions of a python script.  
