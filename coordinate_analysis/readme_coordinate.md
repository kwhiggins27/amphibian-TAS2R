# README for coordinate analysis  

## Description  
Scripts analyzing TAS2R coordinates in several different ways, including: identification of clusters (several methods) and identifying surrounding repeat elements.  

### Identifying clusters  
As described in manuscript, clusters are defined either as all genes with no more than a given interval between genes (100kb, 200kb, 500kb, 1mb, 2mb, or 5mb) or by starting with a fixed interval and then subdividing based on medians.  From here, we performed additional calculations such as calculating the fraction of genes that had the most similar TAS2R sequence within the same cluster.  


## How to install and run  
Download all of these scripts and place in a single directory.  
These scripts assume that the gene identification pipeline has already been run for a list of accession numbers, generating a particular file tree structure containing output files of a certain format  

## How to use  
For all batched files:  
  Update email address (line 4) and desired notification settings (line 3)  
  Update location of output and error files (lines 10-11).  
  Don't change job name (line 2) without careful thought.  Many jobs reference each other using this identifier.  
### Identifying clusters  
batch_find_clusters_nearest_neighbor.sh runs cluster_pipeline_nearest_neighbor.sh which runs cluster_pipeline_nearest_neighbor.py for all accession numbers in accessions_to_keep.txt  
  batch_find_clusters_nearest_neighbor.sh:  
    Header changes as above  
    Line 15: update size of allowed window between genes (I used 100kb, 200kb, 500kb, 1mb, 2mb, or 5mb).  Note that this expects a number with no letter abbreviations (ex: 10000 rather than 10kb)  
    Line 26: update path to list of accession numbers to analyze  
    Both unix scripts need to be in the same directory.  If not, update line 21.  
    Line 24: update location of desired list of failures  
  cluster_pipeline_nearest_neighbor.sh  
    Header changes as above  
    Line 60: update desired location of summary output file 1 (which will be added to with every run)  
    Line 65: update location of python script.  
  cluster_pipeline_nearest_neighbor.py  
    Line 14: update location of pipeline subdirectories  
    Line 158: update desired location of summary output file 2 (similar to in previous script but with more information)  
batch_find_clusters_median_method.sh runs find_clusters_median_method.py for all accession numbers in accessions_to_keep.txt  
  batch_find_clusters_median_method.sh  
    Header changes as above  
    Line 23: update path to list of accession numbers to analyze  
    Both scripts need to be in the same directory.  If not, update line 18.  
    Line 21: update location of desired list of failures  
  find_clusters_median_method.py  
    Lines 12 and 15: update general location of pipeline output gtfs for each accession  
    Line 47: making a sorted gtf in same directory as pipeline output.  Make location match.  
    Line 289: set location of output file with full details on split clusters (not used often).  
    Line 292: set location of key output file.  

### Summarizing cluster coordinates and gene coordinates  
coordinates_of_singleton_cluster.py  
  Creates singletons.bed and clusters.bed in the original pipeline subdirectories  
  Line 10: update location of pipeline subdirectories  
  Line 14: update the name of your output file from your preferred cluster identification run (I used fixed window 100kb)  
  Line 15: if you changed the default name of your pipeline's gtf file, change here to match  
  Line 65: update to reflect list of accessions to analyze  

expand_coordinates.py  
  Adds a window around each singleton or cluster.  Used for repeat element analysis (see below).  
  Line 6: update to reflect list of accessions to analyze  
  Line 13 and 54: update to reflect location of output singletons and clusters files, respectively, from coordinates_of_singleton_cluster.py  
  Lines 48 and 84: update for location of desired output files.  
  Lines 34, 43, 70, 79: update with desired window (here 100kb)  

chr_length_v2_simple.sh finds the length of each chromosome in each genome file.  This takes a while to run but it can process big genomes (just make sure you point to the original genome, not the one that has been split).  Output is a file called chromsizes.csv in each pipeline directory.  
  Line 15: point towards generic location of gene identification pipeline directory  
  Line 21: point towards generic location of genome  
  Line 22: alternate used for big genomes  
  Line 36: point towards list of accessions to use  
  Line 37: alternate used for big genomes  

check_where_genes_v4.py  
  Header changes as above  
  Line 18: set threshold below which a chromosome will be labeled as "short" and not included in the main analysis.  
  Line 19: currently set to 0.1.  This means that if the gene falls within the first 10% or last 10% of the length of the chromosome, it will be called a "near ends" gene  
  Lines 21-22: set output file locations  
  Line 25: list of accessions to include  
  Line 30: point towards the generic folder output from your gene identification pipeline  

is_gene_near_end.R creates plots Fig. 2E  
  Line 7: set working directory where files will be saved  
  Line 10: point towards output of check_where_genes_v4.py  
  Line 23: point towards output of expand_coordinates.py  
  Line 50: set output  
  Line 229-238: create plot for amphibian genes (Fig. 2E left)  
  Lines 262-271: create plot for non-amphibian genes (Fig. 2E right)  

process_rerun.R creates plots Fig. 1A, 2B, 2C, and 2D  
  Line 7: set working directory while files will be saved  
  Line 9: set your own API key.  Only needed for initial analysis  
  Line 10: point towards spreadsheet summarizing details for each accession, critically number of TAS2Rs identified  
  Line 11: point towards output of cluster identification pipeline  
  Line 12: information about genome (optional, not needed to create 3 key figures)  
  Line 13: point towards spreadsheet showing position of each gene relative to the end of the chromosome  
  Line 14: point towards spreadsheet showing genome size  
  Line 100: creates key output file summarizing all data collected so far.  This will be used in future analyses.  
  Line 133-138: creates Fig. 1A panel  
  Lines 140-145: creates Fig. 2B  
  Lines 147-152: creates Fig. 2C  
  Lines 154-159: creates Fig. 2D  

GenePlot.pl  
  Written by Matthew Hill, shared with permission.  
  Dedicated readme in directory  
  Used to generate Fig. 2A  


### Repeat element analysis  
prep_nonamph_genomes.py generates a semi-random list of comparator species as described in paper  
  Line 10: update with location of spreadsheet containing each accession and associated class.  
  Line 32: update with location of pipeline subdirectories  
  Lines 26 and 67: update with desired location of output failures  
  Line 72: may wish to update number of threads with which RepeatModeler will run (I used 10)  

RepeatModeler.sh followed by RepeatModeler2.sh generates a repeat library for each genome  
RepeatModeler.sh  
  Header changes as above  
  Text generated by prep_nonamph_genomes.py and copied into file with appropriate header  
  Recommend running this in chunks in parallel, each containing perhaps 5-10 commands  
RepeatModeler2.sh  
  Header changes as above  
  Text generated by prep_nonamph_genomes.py and copied into file with appropriate header  
  Note that this analysis uses common names.  Check to make sure each unique.  
  Recommend running this in chunks in parallel, each containing perhaps 5 commands  

RepeatMasker.sh identifies repeat elements in both regions surrounding TAS2Rs and in random regions of the genome  
  Line 19: update with generic location for for_py2.py script (python config file from pipeline)  
  Lines 23 and 25: update with generic location of each genome  
  Lines 29-60: reformatting bed file.  May not be necessary for your analysis.  
  Line 67: update with location of output from RepeatModeler2.sh  
  Lines 69-75: update with coordinates of desired regions to search (i.e. singletons, clusters, random genomic regions).    
  Lines 89, 94, and 99: update with desired output directory  

tabulate_percentages.py makes a nice table summarizing the repeat analysis results for a given type (shown for "LTR elements").  I recommend running this for each type of repeat element of interest.  
  Line 6: update with output directory from RepeatMasker.sh  
  Lines 15, 17, 19: update with names of output files from RepeatMasker.sh  
  Line 35: specify type of repeat element to be studied  
  Line 51: output file name.  Strongly encourage using something with the repeat element type.  

## Credits  
Code developed by K. Higgins with experimental advice from R. Márquez, D.W. Bellott, A. Itoigawa, and J.-K. Weng (authors).  Script GenePlot.pl written by Matthew Hill (acand included with permission.  ChatGPT v3 and 4 used to draft early versions of these scripts.  
