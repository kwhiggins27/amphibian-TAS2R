# README for TAS2R identification pipeline  

## Description  
The full code used for high throughput TAS2R gene identification in our manuscript.  
This program is only designed to run for intronless genes (like TAS2Rs)  
It is assumed that every genome has a unique accession number (ex: GCA_001660625.3).  
The pipeline is designed to run through all desired genomes first.  A fraction of them will crash for known reasons.  Genomes that crash because the genome was too large (at the initial blast stage) can be rerun with the "large genome" settings as described below.  The pipeline also crashes if the common name for the species (on NBCI) contains a special character, most commonly an ampersand (&) or apostrophe (').   

## How to install and run  
Download all scripts and place in a directory accessible by a slurm cluster.  
Ensure that you have the appropriate versions of all packages (listed at the beginning for py scripts), in addition to makeblastdb, tblastn, seqtk, cd-hit, tmbed.  Specifically, tmbed must be installed in the directory scripts/blast_pipeline/tmbed.


## How to use  

### Expected file structure:  
Ensure that directories: configs, genomes, logs, progress, references, scripts, and subdirs are downloaded and are in the configuration shown on the repository.  An example of the expected output (minus several large files) is shown for two genomes in subdirs.  Please note that this output also includes output from other scripts in this project.


### Iteratively run format_genome.sh to generate config files
  I ran this in many steps.  I created a script called format_genome_batch.sh which submits one job of format_genome.sh for each accession.  An example version of format_genome_batch.sh is included to illustrate syntax, but you'll need to generate your own (based on your unique paths).  
  Lines 30-53: update this to reflect the references you'll use for the reverse blast step.  For each reference, you'll need a genome with matched coordinates of known TAS2Rs (as a csv).  I've included information to download the genomes we used in references/reference_genomes and the corresponding coordinates in references/reference_coordinates.  It's very important that this contain ALL TAS2Rs in this genome with the correct coordinates because genes will be omitted from followup if they don't match something on this list.   
  Note that this will create config files with big_genome set to "no" by default.  This pipeline is designed to attempt to run on every genome and fail for unusually big genomes.  Plan to manually modify the config files for these and rerun  

### pipeline_v7.sh  
  Edit email in line 134 and 223.  

### blastF.sh, bitter_pt1.py, header_pull.txt, blastR.sh  
  No changes necessary.

### bitter_pt2.py  
  Line 442: creates arbitrary temp names for genes and transcripts.  Feel free to modify.   

### bitter_pt3_TMbed.py  
  Ensure tmbed file in the expected directory (in scripts/blast_pipeline/tmbed).  If necessary, modify line 253.  Note that the downloaded file tree also contains several nested directories called tmbed.  Errors occur if the wrong number of nested directories occur.

### Managing large genomes:  
  Of note, when this pipeline is first run, any large genomes will fail at the blastF.sh step (as of project development, tblastn can't handle large genomes)  
  For each of these, open the 3 config files and manually change the large genome flag to "yes"  

#### split_try3_pipe.sh, blast_big_pipe.sh, bitter_pt1_Big.py, header_big.txt
  No changes necessary.

### Edit batch submitter (or write your own) that will run pipeline_v7.sh for each accession number  
  ../progress/all_accessions_remaining.txt contains all accession numbers that need to run but have not been submitted yet.  If a job fails I want to rerun it, I manually add the accession number back to this list.

#### automated_add_new_jobs.sh  
  When the number of jobs running with the username kwh1 and the job name pipe (set in line 25) falls below 25 (set in line 17), submit 5 new jobs using add_new_jobs.sh.  Check every 30 sec (set in line 56).  Also updates progress bar using script progress.py

#### add_new_jobs.sh  
  Called by automated_add_new_jobs.sh, calls pipeline_v7.sh  
  Header changes as above  
  Script must be located in same directory as script for pipeline_v7.sh.  If not, update line 13  

#### progress.py   
  If different colors desired (ex: to make color-blind friendly) modify 58-61
  Output saved to ../progress

### comparison_with_prior_literature/gene_count_comparison.R creates Supplemetary Fig. 1
  This script is to compare the gene counts from the gene identidication pipeline to the previous studies.
  Necessary files (gene counts of this study and a previous study) are included in this folder (comparison_with_prior_literature).

## Credits  
Code developed by K. Higgins with programming advice from Allen B. Davis (acknowledged) and experimental advice from A. Itoigawa, J.-K. Weng, R. Márquez, and R. Anderson (authors) as well as Matthew Hill and George Bell (acknowledged).  ChatGPT v3 and 4 used to draft early versions of these scripts.  
