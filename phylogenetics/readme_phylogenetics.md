# README for phylogenetics

## Description
A set of scripts to perform basic phylogenetic analyses, use BUSCO scores to identify syntanic regions, and identify copy number constrained orthologs.

## How to install and run
Download all scripts into a single directory.  Run in unix, python, or R as appropriate.


## How to use
For all batched files:
  Update email address (line 4) and desired notification settings (line 3)
  Update location of output and error files (lines 10-11).
  Don't change job name (line 2) without careful thought.  Many jobs reference each other using this identifier.

### Basic phylogenetics tools
iqtree.sh
  Header changes as above
  After -s: change path to fasta
  After -m: this was determined to be the optimal model using the MFP flag

mafft_skeleton.sh
  Header changes as above
  Several different variations have been included to show how mafft was used at different stages in this project
    Lines 22-25 used for creating skeleton tree using just genes in reference list (a few hundred)
    Lines 15-19 used for adding to this skeleton
    Lines 28-33 used for attempting a tree made with all genes from scratch.  This didn't work.

### Using BUSCO scores to identify syntanic regions
A novel method of using BUSCO genes as markers to define syntanic regions between rapidly evolving regions of DNA.
Run scripts in this order:
0) Run BUSCO on genomes.  Expected output has file structure busco/$accession/busco_output/run_vertebrata_odb10/full_table.tsv
1) Edit and run get_accepted_ranges.py
  Line 6: update to point to a list of accession numbers to include
  Lines 13 and 30: update to point to a bed file containing the coordinates of all clusters or singletons (respectively) with 1MB margins
  Line 77: update with desired location of output file
2) Edit and run getmnqk.py
  Header changes as above
  Line 14: link to output of step (2)
  Lines 35 and 6: update to point to generic busco output for any two accessions
  Lines 138-141: update with desired location of four output files
3) Edit and run hypergeometic_extraq.py
  Header changes as above
  Lines 21, 25, 29, 32: update with output from step (2)
  Line 38: update with output of step (1)
  Line 96: update with desired location of output file
4) Edit and run networks.py
  Header changes as above
  Line 18: update with output of step (2)
  Lines 93 and 112: update with desired location of output files
  Lines 26-46: attempt to run Benjamini-hotchberg correction, which never finished.  Went back to Bonferroni which appears in line 52.
5) Edit and run convert_names.py to display results with order (or some other property)
  Line 10: update to point to output of step (4)
  Line 39: update to point to a file which contains accessions, orders, etc.  This will be the basis of conversion.
  Line 64: update with desired location of output file
6) Edit and run cluster_age.R to calculate the minimum possible age of each cluster (Fig. 3C)
  Line 6: update to point to working directory containing output from step (5) and where the output from this analysis will be placed
  Line 9: update with path to output from step (5)
  Line 37: converter file containing accessions, latin names, orders, etc
  Line 91: a short list of species were excluded because they didn't have exact matches in our reference tree (see manuscript).  They were excluded here.
  Line 119: our phylogenetic tree file
  Lines 151, 164, 234, 243, 278, 281: location of output files (281 is the key figure we used in the paper)

### Identifying copy number constrained orthologs
single_copy_genes_local.R
  Line 9: set working directory
  Lines 20-24 and 89-101: subset of big all-TAS2R tree for different clades
  Line 40: the big all-TAS2R tree
  Line 42: a converter file that has accession number, order, and class
  Lines 340-363: set desired output parameters here.  How do you want to categorize each type of results?
  Lines 434-2406: manually inputed annotations of phylogenetic trees with first and last element.  I couldn't figure out a way to automate this in a robust manner.

  Much of the code after 2406 is unnecessary, or only necessary for specific functions.  Here are the broad brushstrokes:
    Lines 2409-2413: check to see if any sequences have been omitted
    Lines 2416-2507: create summary plot similar to Fig. 3B (multiple groups were combined/simplified for Fig. 3B)
    Lines 2514-2521: looking for genes that are duplicated more often than expected.  Didn't end up using this.
    Lines 2554-2605: calculate proportion of genes with different CNCO descriptors for different clades.
    Most of the rest of the code redoes calculations that are also present in other scripts.

### Pairwise orthologs between closely related species
Code to make Fig. 3A (some modification in Illustrator)

subset_tree.R takes the big all-TAS2R tree and subsets just those sequences from 3 accessions (2 for main comparison and 1 as outgroup)
  Line 5: update working directory
  Line 6: update path to big all-TAS2R tree
  Line 10: put 3 accession numbers here
  Lines 16, 19, 22: put unique accession number in each of these (corresponding to line 10)
  Lines 26-35: set thresholds for confidence value nodes
  Lines 42-45: creates a mini tree with nodes
  Lines 47-49: creates the same tree but without labeled nodes

Now there's a slow manual step where you open the trees and score each TAS2R as 1-to-1, 1-to-many, many-to-many, 1-to-0, or many-to-0.

plot_pairs.R creates nice color-coded plots for your comparisons
  Line 6: point to spreadsheet summarizing results of all trio comparisons
  Line 7: make sure column names match your spreadsheet
  Line 21: adjust colors if desired.
  Line 40: first species for comparison (not outgroup)
  Line 41: second species for comparison (not outgroup)
  Line 42: directory to save output graphic


## Credits
Code developed by K. Higgins (author) with programming advice from Allen B. Davis (acknowledged) and experimental advice from R. Márquez and D.W. Bellott (authors).  ChatGPT v3 and 4 used to draft early versions of this script.
