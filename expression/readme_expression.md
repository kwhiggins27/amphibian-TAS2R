# README for expression

## Description
Scripts to accomplish three major steps:
1) Trimming raw RNA seq reads
2) Analyzing individual read files with STAR and featureCounts, and tabulating the results
3) Processing and visualizing the output

## How to install and run
Trimming (step 1) meant to be run on an LSF cluster.  It can be easily adapted to run on SLURM.
Bulk analysis (step 2) meant to be run on a slurm cluster.  The script exp_pipe_v1.sh runs other steps.  Be sure to save all scripts into the same directory.
Vizualization (step 3) is primarily an Jupyter notebook meant to be opened in a Jupyter environment.  PCA is also performed in R.

## How to use
For all batched files:
  Update email address (line 4) and desired notification settings (line 3)
  Update location of output and error files (lines 10-11).
  Don't change job name (line 2) without careful thought.  Many jobs reference each other using this identifier.

### Step 1: Trimming with trimmomatic
axolotl_trimmomatic_pt1.bat contains specifics about file names and submits batch jobs for all reads for a single species.  Despite the file name, this was used for all five species.
  Line 1: contains generic path to raw read files
  Lines 3 and 4: contains generic file suffix
  Lines 6-124: contains specific details for each species and a short name for calling them
  Lines 128-134: this is the actual command.  Loops through all short names (line 128) and submits one job of axolotl_trimmomatic_pt2.bat for each
  Lines 134-138: plug these into line 128 to run loop for different species
  To run this on slurm will need to change bsub command in line 130
axolotl_trimmomatic_pt2.bat runs trimmomatic for the forward and reverse reads of the sample specified in axolotl_trimmomatic_pt1.bat
  Header changes as above.  Note that this header is appropriate for LSF and will require substantial changes to adapt for SLURM.
  Line 13: contains generic path to raw read files
  Lines 14 and 15: contains generic file suffix
  Variable file name was passed as a variable ($1) by axolotl_trimmomatic_pt1.bat
  Line 20: update with location of trimmomatic .jar file
  Line 25: details about custom adapters used

### Step 2: Bulk analysis
exp_pipe_v1.sh is the core script that runs batch jobs for each species.  I included the code to run all five species, but I strongly recommend only submitting one at a time.  The code is written to allow restarting at intermediate steps by updating the go_to variable as noted.  
  Header changes as above.  Note that this header is appropriate for SLURM.
  Line 14-79: update all paths as appropriate.  Notes about key variables:
    species: the shortname you use for each species.  It is used in the tabulation section.
    genome: the location of the genome
    GTF: expected format is GTF output from gene discovery pipeline.
    Bitter_TS: also output from bitter pipeline.
    summary_table: what you want to call your output file
    jobname: what to name the STAR/FC run.  This should be a unique name for each species if you run multiple species in parallel.  The next step won't proceed until all of the jobnames for a species are complete.
    STAR_index: decide where to save your STAR index
    STAR_out: decide where to save your STAR output
    FC_out: decide where to save your featureCounts output
  Line 82: folder containing all of your trimmed reads
  Lines 84 and 85: generic suffix of trimmed files
  Lines 113, 115, 117, 119, 121: the unique part of the trimmed file names
  Line 144: update with your username on your cluster

star_index_pipeline.sh creates a STAR index
  Header changes as above.  Consider lowering memory request for smaller genomes.
  Variables optimized for large genomes

star_fc_pipeline.sh runs STAR and featureCounts on each forward/reverse pair.
  Header changes as above  
  Lines 13-20: describes expected variables to be received
  Lines 25 and 26: makes necessary directories if they don't already exist
  Line 28: ulimit appropriate for large genomes
  Lines 31-42: run STAR with settings appropriate for large genomes
  Lines 51-54: runs featureCounts
  Line 61: believed to reduce crashes, with poor support.  Consider removing.

tabulate_expression.py tabulates data for analysis in python.  Much of this is specific for our file names and data structure and may be of limited utility.
  Header changes as above
  Lines 38-53: assigns species to each unique file name
  Lines 105-125 assign names to each column of the output file.  This assumes three replicates (Counts_1.X, Counts_2.X, Counts_3.X) with seven unique tissues each (the number after the decimal).

### Step 3: visualization
Expression_analysis_colocalization.ipynb is a Jupyter notebook
  Ensure all packages are available and up to date
  This script contains a handful of direct file paths all of which contain "wengpj01".  I recommend searching for this string and updating all output files appropriately.  
  There are a lot of extra analyses here that were performed as sanity checks and never appeared in the manuscript.  They were included in case the reader has similar concerns.  Here's how to generate our key figures:

  Fig 4C: plot of fraction expressed by species.  See cell starting with "# Generates figure 4C".  Note that we're using function count_rows_above_threshold_min1 with a threshold of 0.01.  That means that we're counting all genes that are detected above 0.01 TPM in at least one sample.  We tried other thresholds and requiring a minimum of two independent detections as well.

  Fig 4D: plot of fraction expressed in each tissue.  See cell starting with "# Generates figure 4D".  If you follow this back several cells, you'll see that this is generated with calculate_percentage from combined_threshold_avg_1.  Note that now we're talking about having an average count above 0.001 TPM.  

  Fig 4E:plot of percent of receptors that are expressed in exactly one tissue. See cell starting with "# Generate figure 4E".  Calculated using functions "calculate_percentage" and "calculate_unique_expression_per_tissue".

  Fig 4F: scatter plot comparing the number of TAS2Rs with the genome with the fraction of TAS2Rs that are expressed in at least one extra-oral tissue (but not the tongue).  Trendline, correlation, and p-value shown (pearson correlation).  See cell starting with "# Generates figure 4F".  

  Fig 4G: heatmaps for all five species showing Spearman correlation between TAS2R repertoire between different pairs of tissues.   See cell starting with "# Generates figure 4G"



PCA.R is an R script designed to plot PCA components 1 and 2 for several key comparisons.  For Cane, Dart, and Clawed, all samples were fresh and passed RIN scores so the only relevant comparison is tissue.  For axolotl and bullfrog, a couple RIN scores were poor so this was included as an alternate comparator with shape.  For bullfrog, a couple frozen samples were also used so an additional plot is shown with variation by both tissue and fresh/frozen.
  Ensure all packages available
  Line 4: update working directory
  Lines 26-41: these are raw data inputed from our lab notebook

TreeTipLables_5amphibian.R is an R script designed to plot the phylogenetic tree of all 5 studied amphibians alongside a heatmap of their expression in different tissues as shown in the appendix.  This will result in two output files with lines in the same order, but they must be combined manually in Illustrator (or a similar program).  This was also subsetted for figure 5 as shown at the end.
  This script contains a bunch of absolute paths that will need to be updated.  I recommend searching for kwhiggins27 and wengpj01 and changing all of these at once.  
  Line 9: this should point to your tree
  Line 13: this should point to your annotation file (likely the output of Step 2)
  Lines 28-47: specify color in which different species will appear in your tree.
  Lines 131-137: this will put colored circles on the nodes based on the node.label, showing degree of confidence.
  Lines 173-216: this contains the code to make subtrees shown in Figure 5, as noted

## Credits
Code developed by K. Higgins with experimental advice from R. Márquez (authors).  ChatGPT v3 and 4 used to draft early versions of this script.
