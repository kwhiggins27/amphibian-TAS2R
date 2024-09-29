# README for coordinate analysis  

## Description  
Scripts analyzing TAS2R coordinates in several different ways, including: identification of clusters (several methods) and identifying surrounding repeat elements.  

### Identifying clusters  
As described in manuscript, clusters are defined either as all genes with no more than a given interval between genes (100kb, 200kb, 500kb, 1mb, 2mb, or 5mb) or by starting with a fixed interval and then subdividing based on medians.  From here, we performed additional calculations such as calculating the fraction of genes that had the most similar TAS2R sequence within the same cluster.  


## How to install and run  
Download all of these scripts and place in a single directory.  
These scripts assume that the gene identification pipeline has already been run for a list of accession numbers, generating a particular file tree structure containing output files of a certain format  

## How to use  
### Identifying clusters  
Run scripts in the order listed below.
Be sure to enable conda environment specified in scripts/phylogenetics/requirements.txt.

#### batch_find_clusters_window.sh runs find_clusters_window.py
  Line 13: update size of allowed window between genes (I used 100kb, 200kb, 500kb, 1mb, 2mb, or 5mb).  Note that this expects a number with no letter abbreviations (ex: 10000 rather than 10kb)  
  Line 21: update path to list of accession numbers to analyze  

#### batch_find_clusters_nearest_neighbor.sh runs cluster_pipeline_nearest_neighbor.sh which runs cluster_pipeline_nearest_neighbor.py   
  Line 13: update size of allowed window between genes (I used 100kb, 200kb, 500kb, 1mb, 2mb, or 5mb).  Note that this expects a number with no letter abbreviations (ex: 10000 rather than 10kb)  
  Line 22: update path to list of accession numbers to analyze

#### batch_find_clusters_median_method.sh runs find_clusters_median_method.py
  Line 13: update size of allowed window between genes (I used 100kb, 200kb, 500kb, 1mb, 2mb, or 5mb).  Note that this expects a number with no letter abbreviations (ex: 10000 rather than 10kb)  
  Line 22: update path to list of accession numbers to analyze

### Summarizing cluster coordinates and gene coordinates  
#### coordinates_of_singleton_cluster.py  
  Creates singletons.bed and clusters.bed in the original pipeline subdirectories  
  Line 14: update the name of your output file from your preferred cluster identification run (I used fixed window 100kb)  
  Line 65: update to reflect list of accessions to analyze  

#### expand_coordinates.py  
  Adds a window around each singleton or cluster.  Used for repeat element analysis (see below).  
  Line 6: update to reflect list of accessions to analyze   
  Lines 34, 43, 70, 79: update with desired window (here 100kb)  

#### chr_length_v2_simple.sh
Finds the length of each chromosome in each genome file.  This takes a while to run but it can process big genomes (just make sure you point to the original genome, not the one that has been split).  Output is a file called chromsizes.csv in each pipeline directory.    
  Line 28: point towards list of accessions to use  

#### check_where_genes_v4.py  
  Header changes as above  
  Line 16: set threshold below which a chromosome will be labeled as "short" and not included in the main analysis.  
  Line 17: currently set to 0.1.  This means that if the gene falls within the first 10% or last 10% of the length of the chromosome, it will be called a "near ends" gene   
  Line 23: list of accessions to include  

#### is_gene_near_end.R creates plots Fig. 2E  
  Line 7: set working directory where files will be saved  
  Line 10: point towards output of check_where_genes_v4.py  
  Line 23: point towards output of expand_coordinates.py  
  Line 50: set output  
  Line 229-238: create plot for amphibian genes (Fig. 2E left)  
  Lines 262-271: create plot for non-amphibian genes (Fig. 2E right)  

#### process_rerun_Sept.R
Creates plots Fig. 1A, 2B, 2C, and 2D  
  Line 9: set working directory while files will be saved  
  Line 19: set your own API key.  Only needed for initial analysis  
  Line 12: point towards spreadsheet summarizing details for each accession, critically number of TAS2Rs identified and desired plotting clade (we used mostly class level, with several subdivided into orders based on level of interest)
  Line 13: point towards output of cluster identification pipeline  
  Line 14: information about genome (optional, not needed to create 3 key figures)  
  Line 15: point towards spreadsheet showing position of each gene relative to the end of the chromosome  
  Line 16: point towards spreadsheet showing genome size  
  Line 22: filename of key intermediate file
  Lines 61-79: uncomment and run these lines of code to gather phylogenetic information about each species from NCBI.  This is slow and prone to crashing, so we recommend running it once, saving the results, and reloading  from the intermediate file in line 22 (at line 80) for future runs.
  Line 154-171: creates Fig. 1A panel  
  Lines 175-182: creates Fig. 2B  
  Lines 184-192: creates Fig. 2C  
  Lines 194-202: creates Fig. 2D  

#### GenePlot.pl  
  Written by Matthew Hill, shared with permission.  
  Dedicated readme in directory  
  Used to generate Fig. 2A  

### Repeat element analysis
Assumes you have already run coordinates_of_singleton_cluster.py and expand_coordinates.py, as above.

#### prep_nonamph_genomes.py
Generates a semi-random list of comparator species as described in paper  
  Line 10: update with location of spreadsheet containing each accession and associated class.    
  Lines 26 and 67: update with desired location of output failures.  These are not used again,
  Line 72: may wish to update number of threads with which RepeatModeler will run (I used 10)  

#### RepeatModeler.sh  
  Part 1 of repeat library generation
  Text generated by prep_nonamph_genomes.py and copied into file with appropriate header  
  Recommend running this in chunks in parallel, each containing perhaps 5-10 commands.
  Highly recommend using filenames matching common name for each genome, as per for_py2.py config file for each species.
#### RepeatModeler2.sh  
  Part 2 of repeat library generation
  Text generated by prep_nonamph_genomes.py and copied into file with appropriate header  
  Note that this analysis uses common names.  Check to make sure each unique.  
  Recommend running this in chunks in parallel, each containing perhaps 5 commands  
  Highly recommend using filenames matching common name for each genome, as per for_py2.py config file for each species.

#### RepeatMasker.sh
Identifies repeat elements in both regions surrounding TAS2Rs and in random regions of the genome.

#### tabulate_percentages.py
Makes a nice table summarizing the repeat analysis results for a given type (shown for "LTR elements").  I recommend running this for each type of repeat element of interest.  
  Line 6: update with output directory from RepeatMasker.sh  
  Lines 15, 17, 19: update with names of output files from RepeatMasker.sh  
  Line 35: specify type of repeat element to be studied  
  Line 51: output file name.  Strongly encourage using something with the repeat element type.  

## Credits  
Code developed by K. Higgins with experimental advice from R. Márquez, D.W. Bellott, A. Itoigawa, and J.-K. Weng (authors).  Script GenePlot.pl written by Matthew Hill (acknowledged) included with permission.  ChatGPT v3 and 4 used to draft early versions of these scripts.  
