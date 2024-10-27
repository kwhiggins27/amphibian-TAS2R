# amphibian-TAS2R

## Description  
Repository accompanying the publication: Higgins, K. W., Itoigawa A., Toda Y., Bellot, D.B., Anderson, R., Ishimaru, Y., Márquez, R., Weng, J-K. Rapid Expansion and Specialization of the Bitter Taste Receptor Family in Amphibians. PLoS Genetics, in press.

## How to install and run  
This project is separated into six modules reflecting different sections of the paper.  It is expected that "pipeline" will be run earlier than the other sections (with the exception of "assay", which is independent) to generate the expected file structure used in the other directories.  For each folder, download the contents into a directory, and open the local README.md to see necessary changes and recommended order of run.  Most scripts are written for unix, python, or R, with a single script in perl.  Many scripts have headings that assume a SLURM job scheduler and will need to be revised for the user's system.

## How to use  
See specifics in README.md files in each subdirectory of scripts.

## Credits  
Majority of code developed by K. Higgins with additional scripts by A. Itoigawa and R. Márquez (authors) as noted.  One script each by Matthew Hill and Allen B. Davis (acknowledged), reproduced with permission.  Experimental design advice provided by all authors, in addition to Matthew Hill and George Bell (acknowledged).  Programming advice from Matthew Hill and Allen B. Davis (acknowledged).  ChatGPT v3 and 4 were used to draft early versions of many scripts, with substantial editing and revision.  
