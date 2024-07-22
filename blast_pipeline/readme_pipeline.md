# README for TAS2R identification pipeline  

## Description  
The full code used for high throughput TAS2R gene identification in our manuscript.  
This program is only designed to run for intronless genes (like TAS2Rs)  
It is assumed that every genome has a unique accession number (ex: GCA_001660625.3).  Alternate code is provided for unusually large genomes.  

## How to install and run  
Download all scripts and place in a directory accessible by a slurm cluster.  
Ensure that you have the appropriate versions of all packages (listed at the beginning for py scripts), in addition to makeblastdb, tblastn, seqtk, cd-hit, tmbed.  


## How to use  
For all batched files:  
  Update email address (line 4) and desired notification settings (line 3)  
  Update location of output and error files (lines 10-11).  
  Don't change job name (line 2) without careful thought.  Many jobs reference each other using this identifier.  

Expected file structure:  
  You will need to create the following directories:  
    A folder containing all of the scripts and summary files (ex: /lab/wengpj01/vertebrate_pipeline/)  
    A folder for log and error files (ex: /lab/wengpj01/vertebrate_pipeline/logs/)  
    A folder for genomes (ex: /lab/wengpj01/genomes/ncbi-genomes-2023-06-07/)  
    A folder for subdirectories (ex: /lab/wengpj01/vertebrate_pipeline/subdirs/).  We will create one folder for each accession  
    A folder for config files (ex: /lab/wengpj01/vertebrate_pipeline/configs/)  


Create file tree and generate config files:  
  I ran this in many steps.  I created a script called format_genome_batch.sh which submits one job of format_genome.sh for each accession.  An example version of format_genome_batch.sh is included to illustrate syntax, but you'll need to generate your own (based on your unique paths).  
  Edit format_genome.sh  
    Header changes as above  
    Expected arguments shown in commented lines 13-18  
    Line 21: edit to reflect parent folder containing all target genomes  
    Lines 31-54:  
      Update this to reflect the references you'll use for the reverse blast step.  For each reference, you'll need a genome with matched coordinates of known TAS2Rs (as a csv).  As an example, we've included the gene list for humans (based off of hg38) in human_tas2r_reference.csv.  It's very important that this contain ALL TAS2Rs in this genome with the correct coordinates because genes will be omitted from followup if they don't match something on this list.  
    Line 66: update location where you want your config files to be sufficient_length_sorted  
    Line 69: update location where you want to create subdirectories for each accession.  Must end with subdirs/${shortname}  
    Lines 93 and 113: echo to match file location in 69  
    Lines 114: define where you want to create a summary file.  Recommend using something unique for each run.  This is not critical.  
  Note that this will create config files with big_genome set to "no" by default.  This pipeline is designed to attempt to run on every genome and fail for unusually big genomes.  Plan to manually modify the config files for these and rerun  
  Run format_genome.sh for each genome however desired (i.e. through batch submitter like format_genome_batch.sh)  
  Confirm that necessary folders and config files have been created.  

pipeline_v7.sh  
  Create text files for variables remaining, complete, attempted, and failed and update locations in lines 23-26  
  Edit email username in line 134 and 223.  
  Lines 99, 118, 119, 171, 172, 248, 264, 281: update location of scripts.  

  Edit scripts run by pipeline:  
    blastF.sh  
      Header changes as above  
      Line 27: update location of query fasta  
    bitter_pt1.py  
      Header changes as above  
      Script includes several different ways of expanding margins around putative TAS2R.    
      Line 137: update location of header_file (creates header for unix pull job to get region surrounding putative TAS2Rs)  
    header_pull.txt  
      Header changes as above  
    bitter_pt2.py  
      Header changes as above  
      Lines 316, 329: update location of list_complete.txt  
      Line 463: here I create arbitrary temp names for genes and transcripts.  In a moment of vanity, I used my own initials.  Feel free to choose your own nomenclature.  
    blastR.sh  
      Header changes as above  
    bitter_pt3_TMbed.py  
      Header changes as above  
      Lines 387, 388, 397: update location of subdirectories  
      Line 390: update location of tmbed  

Managing large genomes:  
  Of note, when this pipeline is first run, any large genomes will fail at the blastF.sh step (as of project development, tblastn can't handle large genomes)  
  For each of these, open the 3 config files and manually change the large genome flag to "yes"  
  Additional scripts to revise:  
    split_try3_pipe.sh  
      Header changes as above  
    blast_big_pipe.sh  
      Header changes as above  
      Line 35: update location of query fasta    
    bitter_pt1_Big.py  
      Header changes as above  
      Line 135: update location of header_big.txt  
    header_big.txt  
      Header changes as above  

Edit batch submitter (or write your own) that will run pipeline_v7.sh for each accession number  
  all_accessions_remaining.txt contains all accession numbers that need to run but have not been submitted yet.  If a job fails I want to rerun it, I manually add the accession number back to this list.  
  automated_add_new_jobs.sh  
    When the number of jobs running with the username kwh1 and the job name pipe (set in line 22) falls below 25 (set in line 14), submit 5 new jobs using add_new_jobs.sh.  Check every 30 sec (set in line 52).  Also updates progress bar using script progress.py  
    Header changes as above  
    Line 13: update location of remaining text file  
  add_new_jobs.sh  
    Called by automated_add_new_jobs.sh, calls pipeline_v7.sh  
    Header changes as above  
    Script must be located in same directory as script for pipeline_v7.sh.  If not, update line 13  
  progress.py  
    Lines 7-10 and 92: update location of files  
    If different colors desired (ex: to make color-blind friendly) modify 63-66  




## Credits  
Code developed by K. Higgins with programming advice from Allen B. Davis (acknowledged) and experimental advice from A. Itoigawa, J.-K. Weng, R. Márquez, and R. Anderson (authors) as well as Matthew Hill and George Bell (acknowledged).  ChatGPT v3 and 4 used to draft early versions of these scripts.  
